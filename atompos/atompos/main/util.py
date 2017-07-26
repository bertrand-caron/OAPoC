import hashlib

import atompos
from django.core.cache import cache
import logging
import os
import re
from subprocess import Popen, PIPE
from tempfile import NamedTemporaryFile
from atb_api import API, HTTPError

# TODO: expand
SUPPORTED_FORMATS = [
  'smiles',
  'inchi',
  'pdb',
  'atb' # used for ATB IDs
]

CACHE_TIMEOUT = 60 * 60 * 24 * 365

OUTPUT_FORMAT = "mol2"
BABEL = "obabel"
BABEL_OPTS = "-o%s --gen2d" % OUTPUT_FORMAT
BABEL_OPTS_3D = "-o%s --gen3d" % OUTPUT_FORMAT
SUCCESS_MSG = "1 molecule converted\n"

try:
    ATB_API_TOKEN = os.environ['ATB_API_TOKEN']
except KeyError:
    raise Exception('Please export ATB_API_TOKEN in your environment.')
ATB_API = API(api_token=ATB_API_TOKEN)
ATB_ERROR_MSG = "Molecule unknown!"

logger = logging.getLogger('atompos')

class ValidationError(Exception):
  pass

class ConversionError(Exception):
  pass

class ATBLoadError(Exception):
  pass

class UnknownElementError(Exception):
  pass


class Molecule:
  def __init__(self, dataStr=None, molid=None):
    self.dataStr = dataStr
    self.molid = molid
    self.atoms = []
    self.bonds = []

  def get_atom(self, id):
    for atom in self.atoms:
      if atom.id == id:
        return atom

  def get_atom_bonds(self, atom, type=None):
    if type == None:
      bonds = self.bonds
    else:
      bonds = filter(lambda b: b.bondType == type, self.bonds)
    return filter(lambda b: b.a1 == atom or b.a2 == atom, bonds)

  def add_atom(self, id, element, elementID, x, y):
    self.atoms.append(Atom(self, id, element, elementID, x, y))

  def add_bond(self, id, a1, a2, bondType):
    self.bonds.append(Bond(self, id, a1, a2, bondType))

  def normalize_positions(self):
    if len(self.atoms) == 0:
      return

    # Shift the top left corner to 0,0
    xs = map((lambda a: a.x), self.atoms)
    ys = map((lambda a: a.y), self.atoms)
    mx = -min(xs)
    my = -min(ys)

    for a in self.atoms:
      a.x += mx
      a.y += my

    # Normalize the coordinates to values between 0 and 1
    xs = map((lambda a: a.x), self.atoms)
    ys = map((lambda a: a.y), self.atoms)
    mc = max(xs + ys)

    for a in self.atoms:
      a.x /= mc
      a.y /= mc

  @property
  def __json__(self):
    data = {
      "atoms": map(lambda a: a.__json__, self.atoms),
      "bonds": map(lambda b: b.__json__, self.bonds)
    }
    if self.dataStr:
      data["dataStr"] = self.dataStr
    if self.molid:
      data["molid"] = self.molid
    return data

class Atom:
  def __init__(self, molecule, id, element, elementID, x, y):
    self.molecule = molecule
    self.id = id
    self.element = element
    self.elementID = elementID
    self.x = x
    self.y = y

  def set_3d(self, x3d, y3d, z3d):
    self.x3d = x3d
    self.y3d = y3d
    self.z3d = z3d

  @property
  def iacm(self):
    bas = self.get_bonded_atoms()
    if self.element == 'C':
      bhs = filter(lambda a: a.element == 'H', bas)
      if len(bas) == 4 and len(bhs) == 0:
        return 13
      else:
        return 12
    elif self.element == 'H':
      if bas and bas[0].element == 'C':
        return 20
      else:
        return 21
    elif self.element == 'O':
      if len(filter(lambda a: a.element == 'C', bas)) == len(bas) and \
          len(bas) > 1:
        return 4
      elif len(bas) > 1:
        return 3
      elif bas and len(filter(lambda a: a.element == 'O' and \
          len(a.get_bonded_atoms()) == 1, bas[0].get_bonded_atoms())) > 1 and \
          bas != self.get_bonded_atoms(5):
        return 2
      else:
        return 1
    elif self.element == 'N':
      if len(bas) > 3:
        return 8
      elif len(bas) == 1:
        return 9
      elif len(self.get_bonded_atoms(5)) > 1:
        return 9
      elif len(filter(lambda a: a.element == 'H', bas)) < 2:
        return 6
      else:
        return 7
    elif self.element == 'S':
      if len(bas) > 2:
        return 42
      else:
        return 23
    elif self.element == 'P':
      return 30
    elif self.element == 'Si':
      return 31
    elif self.element == 'F':
      return 32
    elif self.element == 'Cl':
      return 33
    elif self.element == 'Br':
      return 34

    raise UnknownElementError("Encountered element of type %s" % self.element)

  def get_bonded_atoms(self, type=None):
    return map(
      lambda b: b.a2 if b.a1 == self else b.a1,
      self.molecule.get_atom_bonds(self, type)
    )

  @property
  def __json__(self):
    return {
      "id": self.id,
      "element": self.element,
      "elementID": self.elementID,
      "iacm": self.iacm,
      "x": self.x,
      "y": self.y,
      "x3d": self.x3d,
      "y3d": self.y3d,
      "z3d": self.z3d
    }

class Bond:
  def __init__(self, molecule, id, a1, a2, bondType):
    self.molecule = molecule
    self.id = id
    self.a1 = a1
    self.a2 = a2
    self.bondType = bondType

  @property
  def is_aromatic(self):
    return self.bondType == 5

  @property
  def __json__(self):
    return {
      "id": self.id,
      "a1": self.a1.id,
      "a2": self.a2.id,
      "bondType": self.bondType
    }


def generate_or_load_atb_pdb(molid):
  """Loads a PDB file from the ATB, tries to generate and load a new PDB file if that fails."""
  try:
    return load_atb_pdb(molid)
  except ATBLoadError:
    # Try generating the PDB file first
    logger.debug("Could not retrieve PDB for %s, trying to generate.." % molid)
    generate_atb_pdb(molid)
    return load_atb_pdb(molid)


def load_atb_pdb(molid, store_only=False):
  """Loads a PDB file from the ATB."""
  try:
    data = ATB_API.Molecules.download_file(
        molid=molid,
        atb_format='pdb_allatom_unoptimised',
    )
    first_line = data.split('\n')[0]
    if not re.search("HEADER", first_line):
      if re.search(ATB_ERROR_MSG, data):
        raise ATBLoadError("Molecule does not exist")
      else:
        raise ATBLoadError("Invalid PDB file returned")
  except HTTPError as e:
    if e.code == 404:
      msg = "Molecule does not exist"
    else:
      msg = "Server error: status %s (%s)" % (e.code, e)
    raise ATBLoadError("Could not retrieve PDB from ATB: %s" % msg)
  except Exception as e:
    raise ATBLoadError("Could not retrieve PDB from ATB: %s" % e)

  return get_data(data, "pdb")


def generate_atb_pdb(molid):
  """Generates a PDB file on the ATB."""
  try:
    ATB_API.Molecules.generate_mol_data(molid=molid)
  except HTTPError as e:
    raise ATBLoadError('Could not generate PDB on ATB: (Error was: "{0}")'.format(str(e.read())))
  except Exception as e:
    raise ATBLoadError('Could not generate PDB on ATB: (Error was: "{0}")'.format(str(e)))


def parse_mol2(mol2Str, dataStr=None):
  """Parses the babel output."""
  molecule = Molecule(dataStr)

  # Sections: 0 -> header, 1 -> atoms, 2 -> bonds, 3 -> footer
  section = 0
  for line in mol2Str.split('\n'):
    l = line.strip()
    if len(l) == 0:
      continue

    if re.search("ATOM", l):
      section = 1
      continue
    elif re.search("BOND", l):
      section = 2
      continue
    elif re.search("SUBSTRUCTURE", l):
      section = 3
      continue

    if section == 1:
      parts = re.split("\s+", l)
      # Strip off the atom index
      element = re.match("[A-Za-z]+", parts[5]).group(0)
      elementID = parts[1]
      molecule.add_atom(
        int(parts[0]),
        element,
        elementID,
        float(parts[2]),
        float(parts[3]),
      )
    elif section == 2:
      parts = re.split("\s+", l)

      t = parts[3]
      if t == "1":
        bondType = 1
      elif t == "2":
        bondType = 2
      elif t == "3":
        bondType = 3
      elif t == "am":
        bondType = 4
      elif t == "ar":
        bondType = 5
      elif t == "du":
        bondType = 6
      elif t == "un":
        bondType = 7
      else:
        bondType = 0

      molecule.add_bond(
        int(parts[0]),
        molecule.get_atom(int(parts[1])),
        molecule.get_atom(int(parts[2])),
        bondType
      )

  molecule.normalize_positions()
  return molecule


def parse_mol2_3d(mol2Str, molecule):
  """Adds the 3D coodinates from the babel output to the molecule."""

  # Sections: 0 -> header, 1 -> atoms, 2 -> bonds, 3 -> footer
  section = 0
  for line in mol2Str.split('\n'):
    l = line.strip()
    if len(l) == 0:
      continue

    if re.search("ATOM", l):
      section = 1
      continue
    elif re.search("BOND", l):
      section = 2
      continue
    elif re.search("SUBSTRUCTURE", l):
      section = 3
      continue

    if section == 1:
      parts = re.split("\s+", l)
      id = int(parts[0])
      for atom in molecule.atoms:
        if atom.id == id:
          atom.set_3d(
            float(parts[2]),
            float(parts[3]),
            float(parts[4]),
          )

  return molecule


def parse_pdb_3d(mol2Str, molecule):
  """Adds the 3D coodinates from the pdb data to the molecule."""

  # Sections: 0 -> header, 1 -> atoms, 2 -> bonds, 3 -> footer
  section = 0
  for line in mol2Str.split('\n'):
    l = line.strip()
    if len(l) == 0:
      continue

    if re.search("HETATM", l):
      parts = re.split("\s+", l)
      elementID = parts[2]
      for atom in molecule.atoms:
        if atom.elementID == elementID:
          atom.set_3d(
            float(parts[5]),
            float(parts[6]),
            float(parts[7]),
          )

  return molecule


def get_data(data, fmt=None):
  """Get the atom data from a string. Translates the string using openbabel."""
  data = data.strip()
  fmt = fmt.lower()

  p = Popen(
    "echo \"%s\" | %s -i%s %s -h" % (data, BABEL, fmt, BABEL_OPTS),
    shell=True,
    stdout=PIPE,
    stderr=PIPE
  )

  out, err = p.communicate()
  if len(err) > 0 and err != SUCCESS_MSG:
    raise ConversionError(err)
  molecule = parse_mol2(out.strip(), data)

  if fmt == "pdb":
    molecule = parse_pdb_3d(data, molecule)

  else:
    p = Popen(
      "echo \"%s\" | %s -i%s %s -h" % (data, BABEL, fmt, BABEL_OPTS_3D),
      shell=True,
      stdout=PIPE,
      stderr=PIPE
    )

    out, err = p.communicate()
    if len(err) > 0 and err != SUCCESS_MSG:
      raise ConversionError(err)
    molecule = parse_mol2_3d(out.strip(), molecule)

  return molecule.__json__


def infer_format(data):
  data = data.strip()
  first_line = data.split('\n')[0]
  if re.search("HEADER", first_line):
    fmt = "pdb"
  elif re.search("InChI", first_line):
    fmt = "inchi"
  elif data.isdigit():
    fmt = "atb"
  elif data.count('\n') == 0:
    fmt = "smiles"
  else:
    raise ConversionError("Could not identify data format")

  logger.debug("Inferred %s format" % fmt)
  return fmt


def get_atom_data_atb(args):
  """Returns and caches atom data. Loads the atom data from the ATB."""
  try:
    validate_args_atb(args)
  except ValidationError as e:
    return {'error': e.message}

  # This is safe now, as all have been validated
  molid = args.get("molid")

  # Check if data occurs in cache
  try:
    hash = ATB_API.Molecules.latest_topology_hash(molid=molid)
    if hash['status'] == 'success':
      try:
        cache_key = str(molid) + '_' + str(hash['latest_topology_hash'])
        return cache.get_or_set(cache_key, lambda: generate_or_load_atb_pdb(molid), CACHE_TIMEOUT)
      except (ATBLoadError, ConversionError, UnknownElementError) as e:
        return {'error': e.message}
    else:
      return {'error': 'Could not find molid %d' % molid}
  except HTTPError:
    return {'error': 'Could not find molid %d' % molid}


def validate_args_atb(args):
  molid = args.get("molid")

  if not molid:
    raise ValidationError("Missing molecule ID")

  try:
    int(molid)
  except ValueError:
    raise ValidationError("Molecule ID should be numeric")

  return True


def get_atom_data(args, cache_key=None):
  """Returns and caches atom data. Depending on the input format, loads the data from the ATB or translates it using openbabel."""
  try:
    validate_args(args)
  except ValidationError as e:
    return {'error': e.message}

  # This is safe now, as all have been validated
  fmt = args.get("fmt")
  data = args.get("data")

  try:
    if not fmt:
      fmt = infer_format(data)
    if fmt.lower() == "atb":
      return get_atom_data_atb({'molid': data})
    else:
      cache_key = cache_key or hashlib.md5(data).hexdigest()
      return cache.get_or_set(cache_key, lambda: get_data(data, fmt), CACHE_TIMEOUT)
  except (ConversionError, UnknownElementError) as e:
    return {'error': e.message}


def validate_args(args):
  fmt = args.get("fmt", None)
  data = args.get("data")

  if fmt and not fmt.lower() in SUPPORTED_FORMATS:
    raise ValidationError("Invalid data format")
  elif not data:
    raise ValidationError("Missing molecule data")

  return True

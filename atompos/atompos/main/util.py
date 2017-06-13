import atompos
from django.core.cache import cache
import json
import logging
import os
import re
import socket
from subprocess import Popen, PIPE
from tempfile import NamedTemporaryFile
from urllib2 import urlopen, HTTPError
from atb_api import API, HTTPError

# TODO: expand
SUPPORTED_FORMATS = [
  'smiles',
  'inchi',
  'pdb',
  'atb' # used for ATB IDs
]

OUTPUT_FORMAT = "mol2"
BABEL = "obabel"
BABEL_OPTS = "-o%s --gen2d" % OUTPUT_FORMAT
SUCCESS_MSG = "1 molecule converted\n"

ATB_DIR = os.path.normpath("%s/../atb_files/" % \
  os.path.dirname(atompos.__file__))
try:
    ATB_API_TOKEN = os.environ['ATB_API_TOKEN']
except KeyError:
    raise Exception('Please export ATB_API_TOKEN in your environment.')
ATB_API = API()
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
      "y": self.y
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


def load_atb_pdb(molid, store_only=False):
  file = "%s/%s.pdb" % (ATB_DIR, molid)
  if os.path.isfile(file):
    try:
      with open(file, 'r') as fp:
        data = fp.read()
      if store_only:
        logger.info("ATB PDB already stored for molid %s" % molid)
        return

      logger.debug("Loaded ATB PDB for molid %s from file" % molid)
      return get_atom_pos_atb(molid, data)
    except IOError:
      # Should not happen, but is possible when the file is deleted
      pass

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

  # Store the retrieved PDB file
  with open(file, 'w') as fp:
    fp.write(data)

  if store_only:
    logger.info("Stored ATB PDB for molid %s" % molid)
  else:
    return get_atom_pos_atb(molid, data)

def get_atom_pos_atb(molid, data):
  pos = get_atom_pos({"data": data, "fmt": "pdb"}, molid)

  if not "error" in pos:
    pos["dataStr"] = molid
    pos["molid"] = molid
    # Store the generated OAPoC Position Storage (.ops) file
    file = "%s/%s.ops" % (ATB_DIR, molid)
    with open(file, 'w') as fp:
      fp.write(json.dumps(pos))
    logger.debug("Stored OPS for ATB molid %s" % molid)

  return pos

def generate_atb_pdb(molid):
  try:
    ATB_API.Molecules.generate_mol_data(molid=molid)
  except HTTPError as e:
    raise ATBLoadError('Could not generate PDB on ATB: (Error was: "{0}")'.format(str(e.read())))
  except Exception as e:
    raise ATBLoadError('Could not generate PDB on ATB: (Error was: "{0}")'.format(str(e)))

def get_positions_atb(args):
  try:
    validate_args_atb(args)
  except ValidationError as e:
    return {'error': e.message}

  # This is safe now, as all have been validated
  molid = args.get("molid")

  # Check if data occurs in cache
  cachedPos = cache.get(molid)
  if cachedPos:
    return cachedPos

  # Check if data occurs in persistent storage (.ops file)
  file = "%s/%s.ops" % (ATB_DIR, molid)
  if os.path.isfile(file):
    try:
      with open(file, 'r') as fp:
        data = fp.read()
      logger.debug("Loaded positions from OPS file for molid %s" % molid)
      return json.loads(data)
    except (IOError, ValueError):
      # Should not happen, but is possible when the file is deleted
      pass

  try:
    pos = load_atb_pdb(molid)
  except ATBLoadError:
    # Try generating the PDB file first
    logger.debug("Could not retrieve PDB for %s, trying to generate.." % molid)
    try:
      generate_atb_pdb(molid)
      pos = load_atb_pdb(molid)
    except (ATBLoadError, ConversionError, UnknownElementError) as e:
      return {'error': e.message}
  except (ConversionError, UnknownElementError) as e:
    return {'error': e.message}

  return pos

def validate_args_atb(args):
  molid = args.get("molid")

  if not molid:
    raise ValidationError("Missing molecule ID")

  try:
    int(molid)
  except ValueError:
    raise ValidationError("Molecule ID should be numeric")

  return True


def parse_mol2(mol2Str, dataStr=None):
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
        float(parts[3])
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
  return molecule.__json__

def get_positions(data, fmt=None):
  data = data.strip()
  fmt = fmt.lower()
  if fmt == "atb":
    return get_positions_atb({"molid": data})

  elif fmt == "pdb":
    with NamedTemporaryFile() as fp:
      fp.write(data)
      fp.seek(0)

      p = Popen(
        "%s -ipdb %s %s" % (BABEL, fp.name, BABEL_OPTS),
        shell=True,
        stdout=PIPE,
        stderr=PIPE
      )

      out, err = p.communicate()
      if len(err) > 0 and err != SUCCESS_MSG:
        raise ConversionError(err)

    molecule = parse_mol2(out.strip(), data)

  else:
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

  return molecule

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

def get_atom_pos(args, cache_key=None):
  try:
    validate_args(args)
  except ValidationError as e:
    return {'error': e.message}

  # This is safe now, as all have been validated
  fmt = args.get("fmt")
  data = args.get("data")

  cache_key = cache_key or data
  if len(cache_key) < 1024:
    cachedPos = cache.get(cache_key)
    if cachedPos:
      return cachedPos

  try:
    if not fmt:
      fmt = infer_format(data)
    pos = get_positions(data, fmt)
  except (ConversionError, UnknownElementError) as e:
    return {'error': e.message}

  if len(cache_key) < 1024:
    # Cache for a year (infinitely enough..)
    cache.set(cache_key, pos, 60 * 60 * 24 * 365)

  return pos

def validate_args(args):
  fmt = args.get("fmt", None)
  data = args.get("data")

  if fmt and not fmt.lower() in SUPPORTED_FORMATS:
    raise ValidationError("Invalid data format")
  elif not data:
    raise ValidationError("Missing molecule data")

  return True

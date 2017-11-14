import hashlib
import traceback
from collections import Counter, defaultdict
from subprocess import PIPE

import itertools

import signal
from psutil import Popen
from rdkit import Chem
from rdkit.Chem import AllChem
from timeout_decorator import timeout, TimeoutError

import atompos
from django.core.cache import cache
import logging
import os
import re
from atb_api import API, HTTPError


SUPPORTED_FORMATS = [
  'smiles',
  'inchi',
  'pdb',
  'atb' # used for ATB IDs
]

CACHE_TIMEOUT = 60 * 60 * 24 * 365

TIMEOUT = 60
BABEL = "obabel"
SUCCESS_MSG = "1 molecule converted\n"

ELEMENTS_MAP = dict([(1, "O"), (2, "O"), (3, "O"), (4, "O"), (5, "O"), (6, "N"), (7, "N"), (8, "N"), (9, "N"), (10, "N"),
                     (11, "N"), (12, "C"), (13, "C"), (14, "C"), (15, "C"), (16, "C"), (17, "C"), (18, "C"), (19, "C"),
                     (20, "H"), (21, "H"), (23, "S"), (24, "Cu"), (25, "Cu"), (26, "Fe"), (27, "Zn"), (28, "Mg"), (29, "Ca"),
                     (30, "P"), (31, "Ar"), (32, "F"), (33, "Cl"), (34, "Br"), (35, "C"), (36, "O"), (37, "Na"), (38, "Cl"),
                     (39, "C"), (40, "Cl"), (41, "H"), (42, "S"), (43, "C"), (44, "O"), (45, "C"), (46, "Cl"), (47, "F"),
                     (48, "C"), (49, "C"), (50, "O"), (51, "C"), (52, "O"), (53, "N"), (54, "C"), (55, "I"), (56, "Cl"),
                     (57, "B"), (58, "Se"), (59, "H"), (60, "Cl"), (61, "Br")])

try:
    ATB_API_TOKEN = os.environ['ATB_API_TOKEN']
except KeyError:
    raise Exception('Please export ATB_API_TOKEN in your environment.')
#ATB_API = API(api_token=ATB_API_TOKEN, host='http://scmb-atb.biosci.uq.edu.au/atb-uqbcaron')
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

  def add_atom(self, id, iacm, element, elementID, x3d, y3d, z3d):
    self.atoms.append(Atom(self, id, iacm, element, elementID, x3d, y3d, z3d))

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
  def __init__(self, molecule, id, iacm, element, elementID, x3d, y3d, z3d):
    self.molecule = molecule
    self.id = id
    self.iacm = iacm
    self.element = element
    self.elementID = elementID
    self.x3d = x3d
    self.y3d = y3d
    self.z3d = z3d

  def set_2d(self, x, y):
    self.x = x
    self.y = y

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


def get_data_atb(molid):
  """Loads molecule data from the ATB."""
  try:
    lgf = ATB_API.Molecules.download_file(molid=molid, atb_format='lgf')
  except HTTPError as e:
    if e.code == 404:
      msg = "Molecule does not exist"
    else:
      msg = "Server error: status %s (%s)" % (e.code, e)
    raise ATBLoadError("Could not retrieve LGF from ATB: %s" % msg)
  except Exception as e:
    raise ATBLoadError("Could not retrieve LGF from ATB: %s" % e)

  try:
    molecule = Molecule(molid)
    parse_lgf(lgf, molecule)

    return molecule.__json__
  except:
    raise ConversionError("Invalid data format")


def parse_lgf(lgf, molecule):
  """Parses an lgf file to generate a molecule."""
  nodes = False
  edges = False
  node_columns = {}
  rdmol = Chem.RWMol(Chem.Mol())
  mol_2_rd = {}
  rd_2_mol = {}

  # parse lgf
  # add atoms to molecule
  # add atoms and bonds to rdkit molecule
  for line in lgf.splitlines():
    if '@nodes' in line:
      nodes = True
      edges = False
      continue
    if '@edges' in line:
      nodes = False
      edges = True
      continue
    if nodes:
      if len(node_columns) == 0:
        split = line.split()
        if not 'label' in split or not 'atomType' in split \
            or not 'coordX' in split or not 'coordY' in split or not 'coordZ' in split:
          continue
        node_columns['label'] = split.index('label')
        node_columns['label2'] = split.index('label2')
        node_columns['atomType'] = split.index('atomType')
        node_columns['coordX'] = split.index('coordX')
        node_columns['coordY'] = split.index('coordY')
        node_columns['coordZ'] = split.index('coordZ')

      else:
        split = line.split()
        id = int(split[node_columns['label']])
        iacm = int(split[node_columns['atomType']])
        element = ELEMENTS_MAP[iacm]

        molecule.add_atom(id,
                          iacm,
                          element,
                          split[node_columns['label2']],
                          float(split[node_columns['coordX']]),
                          float(split[node_columns['coordY']]),
                          float(split[node_columns['coordZ']]))

        atom = Chem.Atom(element)
        atom.SetNoImplicit(True)
        rd_id = rdmol.AddAtom(atom)
        mol_2_rd[id] = rd_id
        rd_2_mol[rd_id] = id

    if edges:
      if not 'label' in line:
        split = line.split()
        rdmol.AddBond(mol_2_rd[int(split[0])], mol_2_rd[int(split[1])], Chem.BondType.SINGLE)

  # compute 2D coordinates using rdkit
  AllChem.Compute2DCoords(rdmol)
  conf2d = rdmol.GetConformer()

  for atom in rdmol.GetAtoms():
    v = atom.GetIdx()
    pos2d = conf2d.GetAtomPosition(v)
    molecule.get_atom(rd_2_mol[v]).set_2d(pos2d.x, pos2d.y)

  molecule.normalize_positions()

  # translate rdkit molecule to pdb
  Chem.SanitizeMol(rdmol, sanitizeOps=Chem.SANITIZE_SETHYBRIDIZATION)

  pdb = Chem.MolToPDBBlock(rdmol)
  # use obabel to infer bond types
  parse_mol2(call_babel(pdb, "pdb", "mol2", add_h=False), molecule)


def parse_mol2(mol2Str, molecule):
  """Parses the babel output to get bonds and bond types."""

  bonds = False
  for line in mol2Str.split('\n'):
    l = line.strip()
    if len(l) == 0:
      continue

    if re.search("BOND", l):
      bonds = True
      continue

    if bonds:
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

      a1 = molecule.get_atom(int(parts[1]))
      a2 = molecule.get_atom(int(parts[2]))

      if a1 and a2:
        molecule.add_bond(
          int(parts[0]),
          molecule.get_atom(int(parts[1])),
          molecule.get_atom(int(parts[2])),
          bondType)


def rename_atoms(pdb):
  counts = defaultdict(int)
  def replace(line):
    if 'HETATM' in line:
      name = line[13:17].strip()
      line = line[:13] + '%-4s' % ('%s%d' % (name, counts[name])) + line[17:]
      counts[name] += 1
    return line

  return '\n'.join(itertools.imap(replace, pdb.split('\n')))


def call_babel(data, ifmt, ofmt, gen3d=False, add_h=True, timelimit=TIMEOUT):
  """Calls obabel."""

  if not gen3d:
    cmd = "echo \"%s\" | %s -i%s -o%s%s" % (data, BABEL, ifmt, ofmt, ' -h' if add_h else '')
  else:
    cmd = "echo \"%s\" | %s -i%s -o%s%s --gen3d" % (data, BABEL, ifmt, ofmt, ' -h' if add_h else '')

  p = Popen(
    cmd,
    shell=True,
    stdout=PIPE,
    stderr=PIPE,
    preexec_fn=os.setsid
  )

  @timeout(seconds=timelimit, use_signals=False)
  def comm(p):
    return p.communicate()

  try:
    out, err = comm(p)
  except TimeoutError as e:
    os.killpg(os.getpgid(p.pid), signal.SIGTERM)
    raise ConversionError(e)

  if len(err) > 0 and err != SUCCESS_MSG:
    raise ConversionError(err)
  return out.strip()


def check_atom_pos(rdmol):
  conf = rdmol.GetConformer()
  def get_pos(i):
    pos = conf.GetAtomPosition(i)
    return (pos.x, pos.y, pos.z)

  if any(itertools.imap(lambda c: c[1] > 1,
                        Counter([get_pos(i) for i in xrange(rdmol.GetNumAtoms())]).iteritems())):
    raise ValueError()

def get_data(data, fmt=None):
  """Get the atom data from a string. Translates the string using obabel,
  the ATB lgf generator and rdkit."""
  data = str(data.strip())
  fmt = fmt.lower()

  try:
    # read input data, translate to pdb with obabel
    if fmt == "inchi":
      pdb = call_babel(data, "inchi", "pdb")
    elif fmt == "smiles":
      pdb = call_babel(data, "smiles", "pdb")
    else:
      pdb = call_babel(data, "pdb", "pdb")

    # translate to rdkit molecule
    rdmol = Chem.MolFromPDBBlock(pdb, removeHs=False, sanitize=False)

  except:
    raise ConversionError("Invalid data format")

  # generate 3D coordinates with rdkit
  # translate back to pdb
  try:
    # raise error if two atoms have the same coordinates
    check_atom_pos(rdmol)
    pdb = rename_atoms(Chem.MolToPDBBlock(rdmol))
  except:

    @timeout(seconds=TIMEOUT, use_signals=False)
    def gen3drdkit(rdmol, useRandomCoords):
      AllChem.EmbedMolecule(rdmol, useRandomCoords=useRandomCoords)
      AllChem.MMFFSanitizeMolecule(rdmol)
      AllChem.MMFFOptimizeMolecule(rdmol)
      check_atom_pos(rdmol)
      return rename_atoms(Chem.MolToPDBBlock(rdmol))

    try:
      # generate 3D coordinates with RDKIT and eigenvalues start
      pdb = gen3drdkit(rdmol, False)
    except:
      try:
        # generate 3D coordinates with RDKIT and random start
        pdb = gen3drdkit(rdmol, True)
      except:
        # rdkit failed to generate 3D coordinates, try obabel
        try:
          pdb = call_babel(pdb, "pdb", "pdb", gen3d=True, add_h=False)
          rdmol = Chem.MolFromPDBBlock(pdb, removeHs=False, sanitize=False)
          check_atom_pos(rdmol)
          pdb = rename_atoms(pdb)
        except:
          raise ConversionError("Could not generate 3D coordinates")

  # generate lgf with IACM elements using the ATB
  try:
    ret = ATB_API.Molecules.lgf(pdb=pdb)
    if not ret or not isinstance(ret, dict):
      raise ATBLoadError('Internal Server Error')
    lgf = ret['lgf']
  except HTTPError as e:
    msg = "Server error: status %s (%s)" % (e.code, e)
    raise ATBLoadError("Could not retrieve LGF from ATB: %s" % msg)
  except Exception as e:
    raise ATBLoadError("Could not retrieve LGF from ATB: %s" % e)

  try:
    molecule = Molecule(data)
    parse_lgf(lgf, molecule)

    return molecule.__json__
  except:
    traceback.print_exc()
    raise ConversionError("Invalid data format")


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
    fmt = "pdb"

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
        return cache.get_or_set(cache_key, lambda: get_data_atb(molid), CACHE_TIMEOUT)
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
  except (ConversionError, UnknownElementError, ATBLoadError) as e:
    return {'error': e.message}


def validate_args(args):
  fmt = args.get("fmt", None)
  data = args.get("data")

  if fmt and not fmt.lower() in SUPPORTED_FORMATS:
    raise ValidationError("Invalid data format")
  elif not data:
    raise ValidationError("Missing molecule data")

  return True

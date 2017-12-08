import hashlib
import itertools
import json
import logging
import os
import re
import signal
import traceback
from collections import Counter, defaultdict
from subprocess import PIPE

import requests
from atb_api import API, HTTPError
from cStringIO import StringIO
from django.core.cache import cache
from psutil import Popen
from rdkit import Chem
from rdkit.Chem import AllChem
from timeout_decorator import timeout, TimeoutError

from atompos.main.pdb import PDB_Atom, str_for_pdb_atom, pdb_conect_line, is_pdb_atom_line, get_attribute_from_pdb_line, get_coords_from_pdbline, is_pdb_connect_line

SUPPORTED_FORMATS = [
  'smiles',
  'inchi',
  'pdb',
  'atb' # used for ATB IDs
]

CACHE_TIMEOUT = 60 * 60 * 24 * 365
OBABEL_TIMEOUT = 10
RDKIT_TIMEOUT = 120
REQUEST_TIMEOUT = 45
REQUEST_URL = 'http://fragments.atb.uq.edu.au/fdb/fragments/molecules/'

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

class LoadError(Exception):
  pass

class UnknownElementError(Exception):
  pass


class Molecule:
  def __init__(self, dataStr=None, molid=None):
    self.dataStr = dataStr
    self.molid = molid
    self.atoms = []
    self.bonds = []

  def get_atom_by_id(self, id):
    for atom in self.atoms:
      if atom.id == id:
        return atom

  def get_atom(self, elementID):
    for atom in self.atoms:
      if atom.elementID == elementID:
        return atom

  def get_bond(self, a1, a2):
    for bond in self.bonds:
      if (bond.a1.id == a1.id and bond.a2.id == a2.id) or (bond.a1.id == a2.id and bond.a2.id == a1.id):
        return bond

  def get_atom_bonds(self, atom, type=None):
    if type == None:
      bonds = self.bonds
    else:
      bonds = filter(lambda b: b.bondType == type, self.bonds)
    return filter(lambda b: b.a1 == atom or b.a2 == atom, bonds)

  def add_atom(self, id, element, elementID, x3d, y3d, z3d, iacm=None):
    self.atoms.append(Atom(self, id, element, elementID, x3d, y3d, z3d, iacm))

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

  def parse_lgf(self, lgf):
    """Parses an lgf file to generate a molecule."""
    nodes = False
    edges = False
    node_columns = {}

    # parse lgf
    # add atoms and bonds to molecule
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
          elementID = split[node_columns['label2']]
          x, y, z = float(split[node_columns['coordX']]), \
                    float(split[node_columns['coordY']]), \
                    float(split[node_columns['coordZ']])

          self.add_atom(id, element, elementID, x, y, z, iacm)

      if edges:
        if not 'label' in line:
          split = line.split()
          u, v = int(split[0]), int(split[1])
          self.add_bond(int(split[2]), self.get_atom_by_id(u), self.get_atom_by_id(v), 1)

  def parse_pdb(self, pdb):
    counts = defaultdict(int)
    k = 0
    for line in pdb.splitlines():
      if is_pdb_atom_line(line):
        id = int(get_attribute_from_pdb_line('atom_index', line))
        element = get_attribute_from_pdb_line('element', line)
        x, y, z = get_coords_from_pdbline(line)
        counts[element] += 1
        elementID = '%s%d' % (element, counts[element])
        self.add_atom(id, element, elementID, x, y, z)
      if is_pdb_connect_line(line):
        spl = line.split()
        a1 = self.get_atom_by_id(int(spl[1]))
        for id2 in map(int, spl[2:]):
          a2 = self.get_atom_by_id(id2)
          if not self.get_bond(a1, a2) and not self.get_bond(a2, a1):
            k += 1
            self.add_bond(k, a1, a2, 1)

  def compute_2d_coordinates(self):
    # build pdb
    pdb_buffer = StringIO()
    pdb_buffer.write('COMPND    UNNAMED   \n')

    # add atoms to pdb
    for atom in sorted(self.atoms, key=lambda a: a.id):
      pdb_buffer.write(str_for_pdb_atom(PDB_Atom(
          index=atom.id,
          name=atom.elementID,
          coordinates=(atom.x3d, atom.y3d, atom.z3d),
          element=atom.element,
          charge=None),
        residue_name='UNL'))
      pdb_buffer.write('\n')

    # add bonds to pdb
    for atom in sorted(self.atoms, key=lambda a: a.id):
      connected = sorted(map(lambda a: a.id, atom.get_bonded_atoms()))

      if len(connected) > 0:
        connected = [connected[i:i + 3] for i in xrange(0, len(connected), 3)]
        for vals in connected:
          vals.insert(0, atom.id)
          pdb_buffer.write(pdb_conect_line(vals))
          pdb_buffer.write('\n')

    master = 'MASTER        0    0    0    0    0    0    0    0    %-4d    0    %-4d    0\n' \
             % (len(self.atoms), len(self.bonds))
    pdb_buffer.write(master)
    pdb_buffer.write('END   \n')
    pdb = pdb_buffer.getvalue()
    pdb_buffer.close()

    try:
      # compute 2D coordinates and infer bond types using obabel
      self.parse_mol2(call_babel(pdb, "pdb", "mol2", gen2d=True, add_h=False))
    except:

      @timeout(seconds=RDKIT_TIMEOUT, use_signals=False)
      def comp_rdkit():
        # compute 2D coordinates using rdkit
        rdmol = Chem.RWMol(Chem.Mol())
        mol_2_rd = {}
        rd_2_mol = {}
        pos = {}

        for a in self.atoms:
          atom = Chem.Atom(a.element)
          atom.SetIntProp('id', a.id)
          rd_id = rdmol.AddAtom(atom)
          mol_2_rd[a.id] = rd_id
          rd_2_mol[rd_id] = a.elementID

        for e in self.bonds:
          rdmol.AddBond(mol_2_rd[e.a1.id],
                        mol_2_rd[e.a2.id], Chem.BondType.SINGLE)

        AllChem.Compute2DCoords(rdmol)
        conf2d = rdmol.GetConformer()

        for atom in rdmol.GetAtoms():
          v = atom.GetIdx()
          pos2d = conf2d.GetAtomPosition(v)
          pos[rd_2_mol[v]] = (pos2d.x, pos2d.y)

        return pos

      pos = comp_rdkit()
      for elementID, pos2d in pos.iteritems():
        self.get_atom(elementID).set_2d(pos2d[0], pos2d[1])

      # infer bond types using obabel
      self.parse_mol2(call_babel(pdb, "pdb", "mol2", add_h=False), set2d=False)

    self.normalize_positions()

  def parse_mol2(self, mol2Str, set2d=True):
    """Parses the babel output to get bonds and bond types."""

    id_2_el = dict()
    atoms = False
    bonds = False
    for line in mol2Str.split('\n'):
      l = line.strip()
      if len(l) == 0:
        continue

      if re.search("ATOM", l):
        atoms, bonds = True, False
        continue

      if re.search("BOND", l):
        atoms, bonds = False, True
        continue

      if atoms:
        parts = re.split("\s+", l)
        id_2_el[parts[0]] = parts[1]
        if set2d:
          self.get_atom(parts[1]).set_2d(float(parts[2]), float(parts[3]))

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

        a1 = self.get_atom(id_2_el[parts[1]])
        a2 = self.get_atom(id_2_el[parts[2]])

        if a1 and a2:
          self.get_bond(a1, a2).bondType = bondType

  def infer_iacm(self):
    for atom in self.atoms:
      atom.set_iacm()


class Atom:
  def __init__(self, molecule, id, element, elementID, x3d, y3d, z3d, iacm=None):
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

  def set_iacm(self):
    bas = self.get_bonded_atoms()
    if self.element == 'C':
      bhs = filter(lambda a: a.element == 'H', bas)
      if len(bas) == 4 and len(bhs) == 0:
        self.iacm = 13
        return
      else:
        self.iacm = 12
        return
    elif self.element == 'H':
      if bas and bas[0].element == 'C':
        self.iacm = 20
        return
      else:
        self.iacm = 21
        return
    elif self.element == 'O':
      if len(filter(lambda a: a.element == 'C', bas)) == len(bas) and \
              len(bas) > 1:
        self.iacm = 4
        return
      elif len(bas) > 1:
        self.iacm = 3
        return
      elif bas and len(filter(lambda a: a.element == 'O' and \
                                        len(a.get_bonded_atoms()) == 1, bas[0].get_bonded_atoms())) > 1 and \
              bas != self.get_bonded_atoms(5):
        self.iacm = 2
        return
      else:
        self.iacm = 1
        return
    elif self.element == 'N':
      if len(bas) > 3:
        self.iacm = 8
        return
      elif len(bas) == 1:
        self.iacm = 9
        return
      elif len(self.get_bonded_atoms(5)) > 1:
        self.iacm = 9
        return
      elif len(filter(lambda a: a.element == 'H', bas)) < 2:
        self.iacm = 6
        return
      else:
        self.iacm = 7
        return
    elif self.element == 'S':
      if len(bas) > 2:
        self.iacm = 42
        return
      else:
        self.iacm = 23
        return
    elif self.element == 'P':
      self.iacm = 30
      return
    elif self.element == 'Si':
      self.iacm = 31
      return
    elif self.element == 'F':
      self.iacm = 32
      return
    elif self.element == 'Cl':
      self.iacm = 33
      return
    elif self.element == 'Br':
      self.iacm = 34
      return

    raise UnknownElementError("Encountered element of type %s" % self.element)

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
    raise LoadError("Could not retrieve LGF from ATB: %s" % msg)
  except Exception as e:
    raise LoadError("Could not retrieve LGF from ATB: %s" % e)

  try:
    molecule = Molecule(molid=molid)
    molecule.parse_lgf(lgf)
    molecule.compute_2d_coordinates()

    return molecule.__json__
  except:
    raise ConversionError("Invalid data format")


def get_data_fdb(molid):
  """Loads molecule data from the FDB."""
  try:
    r = requests.post(REQUEST_URL + 'json/', timeout=REQUEST_TIMEOUT, data={'molid': molid})
    if r.status_code != requests.codes.ok:
      r.raise_for_status()
    ret = json.loads(r.content)
    if 'status' in ret and ret['status'] == 'success' and 'json' in ret:
      return json.loads(ret['json'])
    else:
      raise LoadError(molid)
  except HTTPError as e:
    if e.code == 404:
      msg = "Molecule does not exist"
    else:
      msg = "Server error: status %s (%s)" % (e.code, e)
    raise LoadError("Could not retrieve data from FDB: %s" % msg)
  except Exception as e:
    raise LoadError("Could not retrieve data from FDB: %s" % e)


def rename_atoms(pdb):
  counts = defaultdict(int)
  def replace(line):
    if 'HETATM' in line or 'ATOM' in line:
      name = line[12:16].strip()
      line = line[:12] + '%-4s' % ('%s%d' % (name, counts[name])) + line[16:]
      counts[name] += 1
    return line

  return '\n'.join(itertools.imap(replace, pdb.split('\n')))


def call_babel(data, ifmt, ofmt, gen2d=False, add_h=True, timelimit=OBABEL_TIMEOUT):
  """Calls obabel."""

  if gen2d:
    cmd = "echo \"%s\" | %s -i%s -o%s%s --gen2d" % (data, BABEL, ifmt, ofmt, ' -h' if add_h else '')
  else:
    cmd = "echo \"%s\" | %s -i%s -o%s%s" % (data, BABEL, ifmt, ofmt, ' -h' if add_h else '')

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
  except TimeoutError:
    os.killpg(os.getpgid(p.pid), signal.SIGTERM)
    raise

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
      pdb = call_babel(data, "pdb", "pdb", add_h=False)

    molecule = Molecule(dataStr=data)
    molecule.parse_pdb(pdb)
    molecule.compute_2d_coordinates()
    molecule.infer_iacm()

    return molecule.__json__
  except TimeoutError:
    raise ConversionError("Timeout")
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


def get_atom_data_fdb(args):
  """Returns and caches atom data. Loads the atom data from the ATB."""
  try:
    validate_args_atb(args)
  except ValidationError as e:
    return {'error': e.message}

  # Check if data occurs in cache
  try:
    # This is safe now, as all have been validated
    molid = int(args.get("molid"))

    r = requests.post(REQUEST_URL + 'hash/', timeout=REQUEST_TIMEOUT, data={'molid': molid})
    if r.status_code != requests.codes.ok:
      r.raise_for_status()
    hash = json.loads(r.content)
    if 'status' in hash and hash['status'] == 'success' and 'hash' in hash:
      try:
        cache_key = str(molid) + '_' + str(hash['hash'])
        return cache.get_or_set(cache_key, lambda: get_data_fdb(molid), CACHE_TIMEOUT)
      except (LoadError, ConversionError, UnknownElementError) as e:
        return {'error': e.message}
    else:
      return {'error': 'Could not find molid %d' % molid}
  except:
    return {'error': 'Could not find molid %d' % molid}


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
      except (LoadError, ConversionError, UnknownElementError) as e:
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
  except (ConversionError, UnknownElementError, LoadError) as e:
    return {'error': e.message}


def validate_args(args):
  fmt = args.get("fmt", None)
  data = args.get("data")

  if fmt and not fmt.lower() in SUPPORTED_FORMATS:
    raise ValidationError("Invalid data format")
  elif not data:
    raise ValidationError("Missing molecule data")

  return True

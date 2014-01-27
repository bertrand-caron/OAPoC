from django.core.cache import cache
import logging
import os
import re
import socket
from subprocess import Popen, PIPE
from tempfile import NamedTemporaryFile
from urllib2 import urlopen

# TODO: expand
SUPPORTED_FORMATS = [
  'smiles',
  'inchi',
  'pdb'
]

OUTPUT_FORMAT = "mol2"
BABEL = "obabel"
BABEL_OPTS = "-o%s --gen2d -h" % OUTPUT_FORMAT
SUCCESS_MSG = "1 molecule converted\n"

ATB_PDB_URL = "http://compbio.biosci.uq.edu.au/atb/download.py?outputType=" + \
  "v2Top&file=pdb_allatom_unoptimised&molid="

logger = logging.getLogger('atompos')

class ValidationError(Exception):
  pass

class ConversionError(Exception):
  pass

class ATBLoadError(Exception):
  pass

class Atom(object):
  id = 0
  element = ""
  elementID = 0
  x = 0.
  y = 0.
  
  def __init__(self, id, element, elementID, x, y):
    self.id = id
    self.element = element
    self.elementID = elementID
    self.x = x
    self.y = y

  def __repr__(self):
    return str(self.__dict__)

class Bond(object):
  id = 0
  a1 = 0
  a2 = 0
  bondType = 0
  
  def __init__(self, id, a1, a2, bondType):
    self.id = id
    self.a1 = a1
    self.a2 = a2
    self.bondType = bondType

  def __repr__(self):
    return str(self.__dict__)


def load_atb_pdb(molid):
  url = "%s%s" % (ATB_PDB_URL, molid)
  try:
    data = urlopen(url).read()
    first_line = data.split('\n')[0]
    if not re.search("HEADER", first_line):
      raise ATBLoadError("Invalid PDB file returned")
  except Exception as e:
    raise ATBLoadError("Could not retrieve PDB from ATB: %s" % (e.message))

  return get_atom_pos({"data": data, "fmt": "pdb"}, molid)

def get_positions_atb(args):
  try:
    validate_args_atb(args)
  except ValidationError as e:
    return {'error': e.message}

  # This is safe now, as all have been validated
  molid = args.get("molid")

  cachedPos = cache.get(molid)
  if cachedPos:
    return cachedPos

  try:
    pos = load_atb_pdb(molid)
  except (ATBLoadError, ConversionError) as e:
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


def parse_atoms_bonds(mol2Str):
  atoms = []
  bonds = []
  
  # Sections: 0 -> header, 1 -> atoms, 2 -> bonds, 3 -> footer
  section = 0
  for line in mol2Str.split('\n'):
    l = line.strip()
    if len(l) == 0:
      continue

    logger.debug("Line: %s" % l)
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
      atoms.append(Atom(
        int(parts[0]),
        element,
        elementID,
        float(parts[2]),
        float(parts[3])
      ))
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

      bonds.append(Bond(int(parts[0]), int(parts[1]), int(parts[2]), bondType))
  
  return atoms, bonds

def normalize_positions(atoms):
  if len(atoms) == 0:
    return atoms

  # Shift the top left corner to 0,0
  xs = map((lambda a: a.x), atoms)
  ys = map((lambda a: a.y), atoms)
  mx = -min(xs)
  my = -min(ys)
  
  for a in atoms:
    a.x += mx
    a.y += my
  
  # Normalize the coordinates to values between 0 and 1
  xs = map((lambda a: a.x), atoms)
  ys = map((lambda a: a.y), atoms)
  mc = max(xs + ys)
  
  for a in atoms:
    a.x /= mc
    a.y /= mc

  return atoms

def get_positions(data, fmt=None):
  data = data.strip()
  if not fmt:
    first_line = data.split('\n')[0]
    if re.search("HEADER", first_line):
      fmt = "pdb"
    elif re.search("InChI", first_line):
      fmt = "inchi"
    elif data.count('\n') == 0:
      fmt = "smiles"
    else:
      raise ValidationError("Could not identify data format")

    logger.debug("Assumed %s format" % fmt)

  fmt = fmt.lower()
  if fmt == "pdb":
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

    atoms, bonds = parse_atoms_bonds(out.strip())

  else:
    p = Popen(
      "echo \"%s\" | %s -i%s %s" % (data, BABEL, fmt, BABEL_OPTS),
      shell=True,
      stdout=PIPE,
      stderr=PIPE
    )

    out, err = p.communicate()
    if len(err) > 0 and err != SUCCESS_MSG:
      raise ConversionError(err)

    atoms, bonds = parse_atoms_bonds(out)

  return {'dataStr': data, 'atoms': normalize_positions(atoms), 'bonds': bonds}

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
    pos = get_positions(data, fmt)
  except ConversionError as e:
    return {'error': e.message}

  if len(cache_key) < 1024:
    # Cache for a year (infinitely enough..)
    cache.set(cache_key, pos, 60 * 60 * 24 * 365)

  return pos

def validate_args(args):
  fmt = args.get("fmt")
  data = args.get("data")
  
  if fmt and not fmt.lower() in SUPPORTED_FORMATS:
    raise ValidationError("Invalid data format")
  elif not data:
    raise ValidationError("Missing molecule data")
  
  return True

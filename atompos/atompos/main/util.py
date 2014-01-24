from django.core.cache import cache
import logging
import os
import re
from subprocess import Popen, PIPE
from tempfile import NamedTemporaryFile

# TODO: expand
SUPPORTED_FORMATS = [
  'smiles',
  'mol'
]

LICENSE_ERROR = \
  "Mol2Export: License not found for charge plugin group. No charges will " + \
  "be written." + os.linesep

if os.name == 'nt':
  MOLCONVERT = "S:\\Programs\\MarvinBeans\\bin\\molconvert.bat"
else:
  MOLCONVERT = "/usr/local/lib/ChemAxon/MarvinBeans/bin/molconvert"

logger = logging.getLogger('atompos')

class ValidationError(Exception):
  pass

class ConversionError(Exception):
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


def parse_atoms_bonds(mol2Str):
  atoms = []
  bonds = []
  
  # Sections: 0 -> header, 1 -> atoms, 2 -> bonds, 3 -> footer
  section = 0
  for l in mol2Str.split('\n'):
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
      element = re.match("[A-Za-z]+", parts[1]).group(0)
      elementID = re.match(".*([0-9])", parts[1]).group(1)
      atoms.append(Atom(
        int(parts[0]),
        element,
        elementID,
        float(parts[2]),
        float(parts[3])
      ))
    elif section == 2:
      parts = re.split("\s+", l)
      try:
        bondType = int(parts[3])
      except ValueError:
        # Use type 4 for aromatic bonds
        bondType = 4
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

def fix_pdb(data):
  fixed = ""
  for line in data.split('\n'):
    parts = line.split()
    if len(parts) > 3 and parts[0] == "HETATM" and len(parts[2]) != 3:
      fix = "%s0%s" % (parts[2][0], parts[2][1])
      logger.debug("Fixing %s: %s" % (line, fix))
      fixed += "%s %4s %4s %4s %4s %11s %7s %7s %5s %5s %11s\n" % \
        (parts[0], parts[1], fix, parts[3], parts[4], parts[5], parts[6], \
         parts[7], parts[8], parts[9], parts[10])
    elif len(line) > 0:
      fixed += line + "\n"
  return fixed

def get_positions(fmt, data):
  if data.count('\n') > 1:
    logger.debug("Assuming PDB")

    with NamedTemporaryFile() as fp:
      fp.write(fix_pdb(data))
      fp.seek(0)

      p = Popen(
        "%s -2 mol2 %s" % (MOLCONVERT, fp.name),
        shell=True,
        stdout=PIPE,
        stderr=PIPE
      )

      out, err = p.communicate()
      if len(err) > 0 and err != LICENSE_ERROR:
        raise ConversionError(err)

    atoms, bonds = parse_atoms_bonds(out)

  else:
    logger.debug("%s mol2 -2 -s \"%s\"" % (MOLCONVERT, data))
    p = Popen(
      "%s mol2 -2 -s \"%s\"" % (MOLCONVERT, data),
      shell=True,
      stdout=PIPE,
      stderr=PIPE
    )

    out, err = p.communicate()
    if len(err) > 0 and err != LICENSE_ERROR:
      raise ConversionError(err)

    atoms, bonds = parse_atoms_bonds(out)

  return {'dataStr': data, 'atoms': normalize_positions(atoms), 'bonds': bonds}

def get_atom_pos(args):
  try:
    validate_args(args)
  except ValidationError as e:
    return {'error': e.message}

  # This is safe now, as all have been validated
  fmt = args.get("fmt").lower()
  data = args.get("data")

  if len(data) < 1024:
    cachedPos = cache.get(data)
    if cachedPos:
      return cachedPos

  try:
    pos = get_positions(fmt, data)
  except ConversionError as e:
    return {'error': e.message}

  if len(data) < 1024:
    # Cache for a year (infinitely enough..)
    cache.set(data, pos, 60 * 60 * 24 * 365)

  return pos

def validate_args(args):
  fmt = args.get("fmt")
  data = args.get("data")
  
  if not fmt or not fmt.lower() in SUPPORTED_FORMATS:
    raise ValidationError("Invalid data format")
  elif not data:
    raise ValidationError("Missing molecule data")
  
  return True

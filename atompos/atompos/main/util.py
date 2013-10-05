from subprocess import Popen, PIPE
import re

# TODO: expand
SUPPORTED_FORMATS = [
  'smiles',
  'mol'
]

LICENSE_ERROR = \
  "Mol2Export: License not found for charge plugin group. No charges will " + \
  "be written.\n"

MOLCONVERT = "/usr/local/lib/ChemAxon/MarvinBeans/bin/molconvert"

class ValidationError(Exception):
  pass

class ConversionError(Exception):
  pass

class Atom(object):
  id = 0
  element = ""
  element_id = 0
  x = 0.
  y = 0.
  
  def __init__(self, id, element, element_id, x, y):
    self.id = id
    self.element = element
    self.element_id = element_id
    self.x = x
    self.y = y

  def __repr__(self):
    return str(self.__dict__)

class Bond(object):
  a1 = 0
  a2 = 0
  bond_type = 0
  
  def __init__(self, a1, a2, bond_type):
    self.a1 = a1
    self.a2 = a2
    self.bond_type = bond_type

  def __repr__(self):
    return str(self.__dict__)


def parse_atoms_bonds(mol2_str):
  atoms = []
  bonds = []
  
  # Sections: 0 -> header, 1 -> atoms, 2 -> bonds, 3 -> footer
  section = 0
  for l in mol2_str.split('\n'):
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
      print "p1", parts[1]
      element = re.match("[A-Za-z]+", parts[1]).group(0)
      element_id = re.match(".*([0-9])", parts[1]).group(1)
      atoms.append(Atom(
        int(parts[0]),
        element,
        element_id,
        float(parts[2]),
        float(parts[3])
      ))
    elif section == 2:
      parts = re.split("\s+", l)
      try:
        bond_type = int(parts[3])
      except ValueError:
        # Use type 4 for aromatic rings
        bond_type = 4
      # Given indices start at 1, correct this.
      bonds.append(Bond(int(parts[1]) - 1, int(parts[2]) - 1, bond_type))
  
  return atoms, bonds

def normalize_positions(atoms):
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

def get_positions(fmt, data):
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
  
  return {'atoms': normalize_positions(atoms), 'bonds': bonds}

def get_atom_pos(args):
  try:
    validate_args(args)
  except ValidationError as e:
    return {'error': e.message}

  # This is safe now, as all have been validated
  fmt = args.get("fmt").lower()
  data = args.get("data")

  try:
    pos = get_positions(fmt, data)
  except ConversionError as e:
    return {'error': e.message}

  return pos

def validate_args(args):
  fmt = args.get("fmt")
  data = args.get("data")
  
  if not fmt or not fmt.lower() in SUPPORTED_FORMATS:
    raise ValidationError("Invalid data format")
  elif not data:
    raise ValidationError("Missing molecule data")
  
  return True

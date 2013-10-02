import re
from subprocess import Popen, PIPE

# TODO: expand
SUPPORTED_FORMATS = [
  'smiles',
  'mol'
]

LICENSE_ERROR = "Mol2Export: License not found for charge plugin group. No charges will be written.\n"

class ValidationError(Exception):
  pass

class ConversionError(Exception):
  pass

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
      name = re.match("[A-Za-z]", parts[1]).group(0)
      atoms.append((name, float(parts[2]), float(parts[3])))
    elif section == 2:
      parts = re.split("\s+", l)
      # Given indices start at 1, correct this.
      bonds.append((int(parts[1]) - 1, int(parts[2]) - 1, parts[3]))
  
  return atoms, bonds

def normalize_positions(pos):
  xs = map((lambda l: l[1]), pos)
  ys = map((lambda l: l[2]), pos)
  mx = -min(xs)
  my = -min(ys)
  
  pos = map((lambda p: (p[0], p[1] + mx, p[2] + my)), pos)
  
  xs = map((lambda l: l[1]), pos)
  ys = map((lambda l: l[2]), pos)
  cs = xs + ys
  mc = max(cs)
  
  return map((lambda p: (p[0], p[1] / mc, p[2] / mc)), pos)

def get_positions(fmt, data):
  p = Popen(
    "molconvert mol2 -2 -s \"%s\"" % (data),
    shell=True,
    stdout=PIPE,
    stderr=PIPE
  )
  
  out, err = p.communicate()
  if len(err) > 0 and err != LICENSE_ERROR:
    raise ConversionError(err)

  atoms, bonds = parse_atoms_bonds(out)
  
  return {'pos': normalize_positions(atoms), 'bonds': bonds}

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


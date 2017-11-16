from itertools import groupby
from operator import itemgetter

from typing import Optional, Tuple, NamedTuple

Coordinate = Tuple[float, float, float]

S = ('>', '')
CHARGE = ('<', '')
F = ('', '.3f')
I = ('>', '')

PDB_ATOM_INDEX_FIELD = (7, 11)
PDB_ATOM_NAME_FIELD = (13, 16)
PDB_COORD_FIELDS = ((31,38), (39, 46), (47, 54)) # Taken from ftp://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_Letter.pdf, page 194
PDB_ELEMENT_FIELD = (77, 78)
PDB_CHARGE = (79, 80)

PDB_INT = int
PDB_STR = lambda x: str(x).strip()
def PDB_MAYBE_INT_WITH_SIGN(x):
    if x.strip() == '':
        return None
    else:
        charge_int, charge_sign = x
        if charge_sign == '-':
            sign = -1
        elif charge_sign == '+':
            sign = +1
        else:
            raise Exception('Unexpected charge_sign: "{0}" (x="{1}")'.format(charge_sign, x))
        return sign * int(charge_int)

PDB_FIELD_AND_FORMATTER_FOR_ATTRIBUTE = {
    'atom_index': (PDB_ATOM_INDEX_FIELD, PDB_INT),
    'atom_name': (PDB_ATOM_NAME_FIELD, PDB_STR),
    'element': (PDB_ELEMENT_FIELD, PDB_STR),
    'charge': (PDB_CHARGE, PDB_MAYBE_INT_WITH_SIGN),
}

# Taken from: ftp://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_A4.pdf, p183
HETATM_SPECS = (
    (1, 6) + S,   # Record name
    PDB_ATOM_INDEX_FIELD + S,  # Atom serial number
    PDB_ATOM_NAME_FIELD + S, # Atom name
#    (17, 17) + S, # Alternate location indicator
    (18, 20) + S, # Residue name
    (22, 22) + S, # Chain identifier
    (23, 26) + S, # Residue sequence number
#    (27, 27) + S, # Code for insertion of residues
    PDB_COORD_FIELDS[0] + F, # Orthogonal coordinates for X
    PDB_COORD_FIELDS[1] + F, # Orthogonal coordinates for Y
    PDB_COORD_FIELDS[2] + F, # Orthogonal coordinates for Z
    (55, 60) + S, # Occupancy
    (61, 66) + S, # Temperature factor
    PDB_ELEMENT_FIELD + S, # Element symbol; right-justified
    PDB_CHARGE + CHARGE, # Charge on the atom
)

FIRST_COORDS_INDEX = HETATM_SPECS.index(PDB_COORD_FIELDS[0] + F)

PDB_CHARGE_START = HETATM_SPECS[-1][0] - 1

CONECT_SPECS = (
    (1, 6) + S, # Record name
    (7, 11) + I, # Atom serial number
    (12, 16) + I, # Serial number of bonded atom
    (17, 21) + I, # Serial number of bonded atom
    (22, 26) + I, # Serial number of bonded atom
    (27, 31) + I, # Serial number of bonded atom
)

def PDB_FORMAT_STR(pdb_specs):
    all_fields = []

    for i, fields in enumerate(pdb_specs):
        if i == 0:
            diff = 0
        else:
            diff = (fields[0] -1) - (pdb_specs[i - 1][1])
        if diff != 0:
            all_fields.append(' ' * diff)

        all_fields.append(
            '{' + '{i}:{left_or_right}{len}{formatter}'.format(
            i=i,
            len=(fields[1] - fields[0] + 1),
            left_or_right=fields[2],
            formatter=fields[3],
        ) + '}'
        )
    return ''.join(all_fields)

PDB_TEMPLATE = PDB_FORMAT_STR(HETATM_SPECS)

CONECT_TEMPLATE = PDB_FORMAT_STR(CONECT_SPECS)

NO_CHAIN = ' '

PDB_ATOM_RECORDS = ('ATOM  ', 'HETATM')

PDB_CONNECT_RECORDS = ('CONECT',)

PDB_Atom = NamedTuple(
    'PDB_Atom',
    [
        ('index', int),
        ('name', str),
        ('coordinates', Tuple[float, float, float]),
        ('element', str),
        ('charge', Optional[int]),
    ]
)

assert list(map(len, PDB_ATOM_RECORDS)) == list(map(len, [ PDB_ATOM_RECORDS[0] ]*len(PDB_ATOM_RECORDS))), '{0} != {1}'.format(list(map(len, PDB_ATOM_RECORDS)), list(map(len, [ PDB_ATOM_RECORDS[0] ]*len(PDB_ATOM_RECORDS))))

def is_pdb_atom_line(line):
    return line[0:len(PDB_ATOM_RECORDS[0])] in PDB_ATOM_RECORDS

def is_pdb_connect_line(line):
    return line[0:len(PDB_CONNECT_RECORDS[0])] in PDB_CONNECT_RECORDS

def pdb_atoms_in(pdb_str):
    return [
        PDB_Atom(
            index=get_attribute_from_pdb_line('atom_index', line),
            name=get_attribute_from_pdb_line('atom_name', line),
            coordinates=get_coords_from_pdbline(line),
            element=get_attribute_from_pdb_line('element', line),
            charge=get_attribute_from_pdb_line('charge', line),
        )
        for line in pdb_str.splitlines()
        if is_pdb_atom_line(line)
    ]

def str_for_pdb_atom(atom, residue_name = 'RES', chain_id = '', residue_number = 1):
    return PDB_TEMPLATE.format(
        PDB_ATOM_RECORDS[0],
        atom.index,
        atom.name,
        residue_name,
        chain_id,
        residue_number,
        atom.coordinates[0],
        atom.coordinates[1],
        atom.coordinates[2],
        '',
        '',
        atom.element,
        '' if atom.charge is None else (str(atom.charge) + ('+' if atom.charge >= 0 else '-')),
    )

def get_coords_from_pdbstr(pdb_str, filter_empty=False):
    optional_filtering = lambda list_of_lines: list_of_lines if not filter_empty else list(filter(bool, list_of_lines))

    return optional_filtering(
        [
            get_coords_from_pdbline(line)
            for line in pdb_str.splitlines()
        ]
    )

def get_coords_from_pdbline(line):
    if is_pdb_atom_line(line):
        return tuple(
            map(
                float,
                [
                    line[pdb_coord_field[0] - 1: pdb_coord_field[1]]
                    for pdb_coord_field in PDB_COORD_FIELDS
                ],
            ),
        )
    else:
        return None

def get_attribute_from_pdb_line(attribute, line):
    pdb_field, formatter = PDB_FIELD_AND_FORMATTER_FOR_ATTRIBUTE[attribute]
    return formatter(line[pdb_field[0] - 1: pdb_field[1]])

def get_elements_from_pdbstr(pdb_str):
    return [
        get_attribute_from_pdb_line('element', line)
        for line in pdb_str.splitlines()
        if is_pdb_atom_line(line)
    ]

def pdb_fields(line):
    return [
        line[field_start - 1: field_stop]
        for (field_start, field_stop, _, _) in HETATM_SPECS
    ]

def pdb_conect_line(fields):
    return CONECT_TEMPLATE.format(
        *list(PDB_CONNECT_RECORDS) + fields + [''] * (len(CONECT_SPECS) - len(fields) - 1)
    )

def substitute_coordinates_in(pdb_line, coords):
    assert len(coords) == 3, coords
    fields =  tuple(pdb_fields(pdb_line))
    return PDB_TEMPLATE.format(*fields[0:FIRST_COORDS_INDEX] + coords + fields[FIRST_COORDS_INDEX + 3:])

def pdb_atom_lines_number(pdb_str):
    return len([ 1 for line in pdb_str.splitlines() if is_pdb_atom_line(line) ])

def normalised_element(line):
    return line[76:78].replace(' ', '').upper()

def pdb_charge_str(line):
    assert line[79] in '+- ', 'Detected formal atomic charge larger than 9. Line was: "{0}"'.format(line)
    charge_chars = line[78:80]
    if charge_chars.strip() == '':
        return '0+'
    else:
        return charge_chars

def pdb_formula(pdb_str):
    return [
        (element, len(list(occurences)))
        for (element, occurences) in
        groupby(
            sorted([normalised_element(line) for line in pdb_str.splitlines() if is_pdb_atom_line(line)]),
        )
    ]

on_sign = itemgetter(1)

def pdb_total_charge(pdb_str):
    return sum([
        int(sign + '1') * sum([int(x[0]) for x in occurences])
        for (sign, occurences) in
        groupby(
            sorted(
                [pdb_charge_str(line) for line in pdb_str.splitlines() if is_pdb_atom_line(line)],
                key=on_sign,
            ),
            key=on_sign,
        )
    ])

def pdb_formula_string(pdb_str, add_charge):
    if add_charge:
        total_charge = pdb_total_charge(pdb_str)
        if total_charge == 0:
            add_charge = False

    return ''.join(
        [
            '{element}{occurence}'.format(element=element.title(), occurence=(occurence if occurence != 1 else ''))
            for (element, occurence) in pdb_formula(pdb_str)
        ]
        +  [ '{sign}{total_charge}'.format(total_charge=abs(total_charge), sign=('+' if total_charge > 0 else '-')) if add_charge else '']
    )

def remove_charge_from(line):
    if is_pdb_atom_line(line):
        return line[0:PDB_CHARGE_START]
    else:
        return line

def remove_pdb_charges(pdb_str):
    return '\n'.join([remove_charge_from(line)for line in pdb_str.splitlines()])

def replace_ATOM_by_HETATM(x):
    return x.replace(*PDB_ATOM_RECORDS)

if __name__ == '__main__':
    print(PDB_TEMPLATE)
    print(CONECT_TEMPLATE)
    print(pdb_conect_line([1, 2, 3]))

    pdb_str = '''
HETATM   36  H15 UNK  1120      -3.379  -2.942  -1.191  1.00  0.00           B1-
HETATM   36  H15 UNK     0      -3.379  -2.942  -1.191  1.00  0.00           A1-
HETATM   36  H15 UNK     0      -3.379  -2.942  -1.191  1.00  0.00          CL1-
HETATM   36  H15 UNK     0      -3.379  -2.942  -1.191  1.00  0.00           H9-
HETATM   36  H15 UNK     0      -3.379  -2.942  -1.191  1.00  0.00          CL3+
HETATM   36  H15 UNK     0      -3.379  -2.942  -1.191  1.00  0.00          CL5+
HETATM   18  C17 RTL A 201      14.349  -6.384 -14.255  1.00  0.00           C  
'''
    print(pdb_formula_string(pdb_str))
    print(pdb_total_charge(pdb_str))

    for line in pdb_str.splitlines():
        if is_pdb_atom_line(line):
            print(line)
            print(substitute_coordinates_in(line, (0., 0., 0.)))
            print()

    print(remove_pdb_charges(pdb_str))

    print(pdb_atoms_in(pdb_str))
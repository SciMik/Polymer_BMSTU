import sys
import random as rand
import numpy as np

data = open(sys.argv[1], "r").readlines() # GO file with Bonds
double_data = open("double_C_bonds.dat", "r").readlines()

bond_C = 1.73
bond_O = 1.4
bond_H = 1

double_c_dict = dict()
for line in double_data:
    atom1 = int(line.split()[0])
    atom2 = int(line.split()[1])
    double_c_dict.update({atom1:atom2})

atoms_num = int(data[2].split()[0])
bonds_num = int(data[4].split()[0])

atom_pos, bond_pos = 0, 0
for i, line in enumerate(data):
    try:
        if line.split()[0] == "Atoms":
            atom_pos = i
        elif line.split()[0] == "Bonds":
            bond_pos = i
    except:
        nothing = 1

atoms_type_dict = dict()
atoms_coord_dict = dict()
atoms_mol_dict = dict()
for line in data[atom_pos+2:atom_pos+2+atoms_num]:
    atom = int(line.split()[0])
    mol = int(line.split()[1])
    type = int(line.split()[2])
    x, y, z = float(line.split()[3]), float(line.split()[4]), float(line.split()[5])
    atoms_type_dict.update({atom:type})
    atoms_coord_dict.update({atom:[x, y, z]})
    atoms_mol_dict.update({atom:mol})

bonds_dict = dict()
#print(bond_pos)
for line in data[bond_pos+2:bond_pos+2+bonds_num]:
    atom1 = int(line.split()[2])
    atom2 = int(line.split()[3])
#    if atom1 == 546 or atom2 == 546:
#        print("GOOD ATOM")
    try:
        bonds_dict[atom1].append(atom2)
    except:
        bonds_dict.update({atom1:[atom2]})
    try:
        bonds_dict[atom2].append(atom1)
    except:
        bonds_dict.update({atom2:[atom1]})
#print(bonds_num)
#print(len(bonds_dict))
#print(bonds_dict[546])

###### 1 - c, 2 - co, 3 - coh, 4 - ch, 5 - ch2, 7 - ho, 8 - hc, 10 - h2c, 5 - o, 6 - oh

# 1 - c, 2 - co, 3 - coh, 4 - ch, 5 - ch2, 6 - ho, 7 - hc, 8 - h2c, 9 - o, 10 - oh another prog version

# 1 - c, 2 - h, 3 - co, 4 - oh, 5 - ho, 6 - o

type_change_dict = dict()
type_change_dict.update({1:1, 2:3, 3:3, 4:1, 5:1, 6:5, 7:2, 8:2, 9:6, 10:4, 11:11})

h_max = 1
h_num = 0


new_atoms_type = dict()
used_atoms = []
H_atoms = dict()
H_num = atoms_num + 1
H_bonds_dict = dict()

cooh_num = 1
cooh_atoms_dict = dict()
cooh_bonds_dict = dict()
cooh_types_dict = dict()
cooh_types_charges = dict()
cooh_types_charges.update({8:0.635, 9:-0.44, 10:-0.53, 5:0.45})

for atom in bonds_dict:
    atom_type = atoms_type_dict[atom]
    if atom_type == 1 and atom not in used_atoms:
        check_o = 0
        for bonded_atom in bonds_dict[atom]:
            if atoms_type_dict[bonded_atom] == 3:
                check_h = 0
                check_o = 1
                for atom_tmp in bonds_dict[bonded_atom]:
                    if atoms_type_dict[atom_tmp] == 2:
                        new_atoms_type.update({bonded_atom:10})
                        new_atoms_type.update({atom_tmp:6})
                        new_atoms_type.update({atom:3})
                        check_h = 1
                if check_h == 0:
                    new_atoms_type.update({bonded_atom: 9})
                    new_atoms_type.update({atom: 2})
        if check_o == 1:
            used_atoms.append(atom)
            atom2 = double_c_dict[atom]
            used_atoms.append(atom2)
            check_o_2 = 0
            for bonded_atom in bonds_dict[atom2]:
                if atoms_type_dict[bonded_atom] == 3:
                    check_h = 0
                    check_o_2 = 1
                    for atom_tmp in bonds_dict[bonded_atom]:
                        if atoms_type_dict[atom_tmp] == 2:
                            new_atoms_type.update({bonded_atom: 10})
                            new_atoms_type.update({atom_tmp: 6})
                            new_atoms_type.update({atom2: 3})
                            check_h = 1
                    if check_h == 0:
                        new_atoms_type.update({bonded_atom: 9})
                        new_atoms_type.update({atom2: 2})
            if check_o_2 == 0:
                bonds_atom_num = len(bonds_dict[atom2])
                if bonds_atom_num == 3:
                    number = rand.randint(1, 2)
                    if number == 1:
                        H_x = atoms_coord_dict[atom2][0]
                        H_y = atoms_coord_dict[atom2][1]
                        H_z = atoms_coord_dict[atom2][2] + 0.9
                    else:
                        H_x = atoms_coord_dict[atom2][0]
                        H_y = atoms_coord_dict[atom2][1]
                        H_z = atoms_coord_dict[atom2][2] - 0.9
                    H_atoms.update({H_num:[H_x, H_y, H_z]})
                    H_bonds_dict.update({H_num:atom2})
                    new_atoms_type.update({atom2:4})
                    new_atoms_type.update({H_num:7})
                    H_num += 1
                elif bonds_atom_num == 2 and h_num <= h_max:  # add H atom
                    h_num += 1
                    number = rand.randint(1, 2)
                    if number == 1:
                        H_x = atoms_coord_dict[atom2][0]
                        H_y = atoms_coord_dict[atom2][1]
                        H_z = atoms_coord_dict[atom2][2] + 0.9
                    else:
                        H_x = atoms_coord_dict[atom2][0]
                        H_y = atoms_coord_dict[atom2][1]
                        H_z = atoms_coord_dict[atom2][2] - 0.9
                    H_atoms.update({H_num: [H_x, H_y, H_z]})
                    H_bonds_dict.update({H_num: atom2})
                    new_atoms_type.update({atom2: 4})
                    new_atoms_type.update({H_num: 7})
                    H_num += 1
                elif bonds_atom_num == 2 and h_num > h_max:
                    vec1 = [atoms_coord_dict[atom2][0] - atoms_coord_dict[bonds_dict[atom2][0]][0],
                            atoms_coord_dict[atom2][1] - atoms_coord_dict[bonds_dict[atom2][0]][1]]
                    vec2 = [atoms_coord_dict[atom2][0] - atoms_coord_dict[bonds_dict[atom2][1]][0],
                            atoms_coord_dict[atom2][1] - atoms_coord_dict[bonds_dict[atom2][1]][1]]
                    vec_sum = np.array([vec1[0] + vec2[0], vec1[1] + vec2[1]])
                    vec_norm = vec_sum / np.linalg.norm(vec_sum)
                    C_x = atoms_coord_dict[atom2][0] + vec_norm[0] * bond_C
                    C_y = atoms_coord_dict[atom2][1] + vec_norm[1] * bond_C
                    C_z = atoms_coord_dict[atom2][2]
                    cooh_atoms_dict.update({cooh_num: [C_x, C_y, C_z]})
                    cooh_types_dict.update({cooh_num: 8})
                    cooh_bonds_dict.update({cooh_num: [atom2]})
                    new_atoms_type.update({atom2: 11})
                    number_cooh = rand.randint(1, 2)
                    if number_cooh == 1:
                        O_z = C_z - bond_O
                        cooh_atoms_dict.update({cooh_num + 1: [C_x, C_y, O_z]})
                        cooh_types_dict.update({cooh_num + 1: 9})
                        cooh_bonds_dict.update({cooh_num + 1: [cooh_num, 123]})
                        OH_z = C_z + bond_O
                        cooh_atoms_dict.update({cooh_num + 2: [C_x, C_y, OH_z]})
                        cooh_types_dict.update({cooh_num + 2: 10})
                        cooh_bonds_dict.update({cooh_num + 2: [cooh_num, 123]})
                        H_z = C_z + bond_O + bond_H
                        cooh_atoms_dict.update({cooh_num + 3: [C_x, C_y, H_z]})
                        cooh_types_dict.update({cooh_num + 3: 5})
                        cooh_bonds_dict.update({cooh_num + 3: [cooh_num + 2, 123]})
                    else:
                        O_z = C_z + bond_O
                        cooh_atoms_dict.update({cooh_num + 1: [C_x, C_y, O_z]})
                        cooh_types_dict.update({cooh_num + 1: 9})
                        cooh_bonds_dict.update({cooh_num + 1: [cooh_num, 123]})
                        OH_z = C_z - bond_O
                        cooh_atoms_dict.update({cooh_num + 2: [C_x, C_y, OH_z]})
                        cooh_types_dict.update({cooh_num + 2: 10})
                        cooh_bonds_dict.update({cooh_num + 2: [cooh_num, 123]})
                        H_z = C_z - bond_O - bond_H
                        cooh_atoms_dict.update({cooh_num + 3: [C_x, C_y, H_z]})
                        cooh_types_dict.update({cooh_num + 3: 5})
                        cooh_bonds_dict.update({cooh_num + 3: [cooh_num + 2, 123]})
                    cooh_num += 4
                    #H_x = atoms_coord_dict[atom2][0]
                    #H_y = atoms_coord_dict[atom2][1]
                    #H_z = atoms_coord_dict[atom2][2] + 0.9
                    #H_atoms.update({H_num: [H_x, H_y, H_z]})
                    #H_bonds_dict.update({H_num: atom2})
                    #new_atoms_type.update({atom2: 4})
                    #new_atoms_type.update({H_num: 7})
                    #H_num += 1
                    #H_x = atoms_coord_dict[atom2][0]
                    #H_y = atoms_coord_dict[atom2][1]
                    #H_z = atoms_coord_dict[atom2][2] - 0.9
                    #H_atoms.update({H_num: [H_x, H_y, H_z]})
                    #H_bonds_dict.update({H_num: atom2})
                    #new_atoms_type.update({atom2: 4})
                    #new_atoms_type.update({H_num: 7})
                    #H_num += 1



for atom in double_c_dict:
    #if atom == 86:
    #    print(bonds_dict[atom])
    if atom not in used_atoms:
        bonds_atom_num = len(bonds_dict[atom])
        if bonds_atom_num == 3:
            new_atoms_type.update({atom:1})
        elif bonds_atom_num == 2:
            h_cooh_val = rand.randint(1,2)
            if h_cooh_val == 3 and h_num <= h_max:  # add H atom
                h_num += 1
                number = rand.randint(1, 2)
                if number == 1:
                    H_x = atoms_coord_dict[atom][0]
                    H_y = atoms_coord_dict[atom][1]
                    H_z = atoms_coord_dict[atom][2] + 0.9
                else:
                    H_x = atoms_coord_dict[atom][0]
                    H_y = atoms_coord_dict[atom][1]
                    H_z = atoms_coord_dict[atom][2] - 0.9
                H_atoms.update({H_num: [H_x, H_y, H_z]})
                H_bonds_dict.update({H_num: atom})
                new_atoms_type.update({atom: 4})
                new_atoms_type.update({H_num: 7})
                H_num += 1
            else:
                vec1 = [atoms_coord_dict[atom][0] - atoms_coord_dict[bonds_dict[atom][0]][0], atoms_coord_dict[atom][1] - atoms_coord_dict[bonds_dict[atom][0]][1]]
                vec2 = [atoms_coord_dict[atom][0] - atoms_coord_dict[bonds_dict[atom][1]][0], atoms_coord_dict[atom][1] - atoms_coord_dict[bonds_dict[atom][1]][1]]
                vec_sum = np.array([vec1[0] + vec2[0], vec1[1] + vec2[1]])
                vec_norm = vec_sum / np.linalg.norm(vec_sum)
                C_x = atoms_coord_dict[atom][0] + vec_norm[0]*bond_C
                C_y = atoms_coord_dict[atom][1] + vec_norm[1]*bond_C
                C_z = atoms_coord_dict[atom][2]
                cooh_atoms_dict.update({cooh_num : [C_x, C_y, C_z]})
                cooh_types_dict.update({cooh_num : 8})
                cooh_bonds_dict.update({cooh_num : [atom]})
                new_atoms_type.update({atom: 11})
                number_cooh = rand.randint(1, 2)
                if number_cooh == 1:
                    O_z = C_z - bond_O
                    cooh_atoms_dict.update({cooh_num + 1: [C_x, C_y, O_z]})
                    cooh_types_dict.update({cooh_num + 1: 9})
                    cooh_bonds_dict.update({cooh_num + 1: [cooh_num, 123]})
                    OH_z = C_z + bond_O
                    cooh_atoms_dict.update({cooh_num + 2: [C_x, C_y, OH_z]})
                    cooh_types_dict.update({cooh_num + 2: 10})
                    cooh_bonds_dict.update({cooh_num + 2: [cooh_num, 123]})
                    H_z = C_z + bond_O + bond_H
                    cooh_atoms_dict.update({cooh_num + 3: [C_x, C_y, H_z]})
                    cooh_types_dict.update({cooh_num + 3: 5})
                    cooh_bonds_dict.update({cooh_num + 3: [cooh_num + 2, 123]})
                else:
                    O_z = C_z + bond_O
                    cooh_atoms_dict.update({cooh_num + 1: [C_x, C_y, O_z]})
                    cooh_types_dict.update({cooh_num + 1: 9})
                    cooh_bonds_dict.update({cooh_num + 1: [cooh_num, 123]})
                    OH_z = C_z - bond_O
                    cooh_atoms_dict.update({cooh_num + 2: [C_x, C_y, OH_z]})
                    cooh_types_dict.update({cooh_num + 2: 10})
                    cooh_bonds_dict.update({cooh_num + 2: [cooh_num, 123]})
                    H_z = C_z - bond_O - bond_H
                    cooh_atoms_dict.update({cooh_num + 3: [C_x, C_y, H_z]})
                    cooh_types_dict.update({cooh_num + 3: 5})
                    cooh_bonds_dict.update({cooh_num + 3: [cooh_num + 2, 123]})
                cooh_num += 4

        atom2 = double_c_dict[atom]
        bonds_atom_num = len(bonds_dict[atom2])
        if bonds_atom_num == 3:
            new_atoms_type.update({atom2: 1})
        elif bonds_atom_num == 2:
            h_cooh_val = rand.randint(1, 2)
            if h_cooh_val == 3 and h_num <= h_max:  # add H atom
                h_num += 1
                number = rand.randint(1, 2)
                if number == 1:
                    H_x = atoms_coord_dict[atom2][0]
                    H_y = atoms_coord_dict[atom2][1]
                    H_z = atoms_coord_dict[atom2][2] + 0.9
                else:
                    H_x = atoms_coord_dict[atom2][0]
                    H_y = atoms_coord_dict[atom2][1]
                    H_z = atoms_coord_dict[atom2][2] - 0.9
                H_atoms.update({H_num: [H_x, H_y, H_z]})
                H_bonds_dict.update({H_num: atom2})
                new_atoms_type.update({atom2: 4})
                new_atoms_type.update({H_num: 7})
                H_num += 1
            else:
                vec1 = [atoms_coord_dict[atom2][0] - atoms_coord_dict[bonds_dict[atom2][0]][0],
                        atoms_coord_dict[atom2][1] - atoms_coord_dict[bonds_dict[atom2][0]][1]]
                vec2 = [atoms_coord_dict[atom2][0] - atoms_coord_dict[bonds_dict[atom2][1]][0],
                        atoms_coord_dict[atom2][1] - atoms_coord_dict[bonds_dict[atom2][1]][1]]
                vec_sum = np.array([vec1[0] + vec2[0], vec1[1] + vec2[1]])
                vec_norm = vec_sum / np.linalg.norm(vec_sum)
                C_x = atoms_coord_dict[atom2][0] + vec_norm[0] * bond_C
                C_y = atoms_coord_dict[atom2][1] + vec_norm[1] * bond_C
                C_z = atoms_coord_dict[atom2][2]
                cooh_atoms_dict.update({cooh_num: [C_x, C_y, C_z]})
                cooh_types_dict.update({cooh_num: 8})
                cooh_bonds_dict.update({cooh_num: [atom2]})
                new_atoms_type.update({atom2: 11})
                number_cooh = rand.randint(1, 2)
                if number_cooh == 1:
                    O_z = C_z - bond_O
                    cooh_atoms_dict.update({cooh_num + 1: [C_x, C_y, O_z]})
                    cooh_types_dict.update({cooh_num + 1: 9})
                    cooh_bonds_dict.update({cooh_num + 1: [cooh_num, 123]})
                    OH_z = C_z + bond_O
                    cooh_atoms_dict.update({cooh_num + 2: [C_x, C_y, OH_z]})
                    cooh_types_dict.update({cooh_num + 2: 10})
                    cooh_bonds_dict.update({cooh_num + 2: [cooh_num, 123]})
                    H_z = C_z + bond_O + bond_H
                    cooh_atoms_dict.update({cooh_num + 3: [C_x, C_y, H_z]})
                    cooh_types_dict.update({cooh_num + 3: 5})
                    cooh_bonds_dict.update({cooh_num + 3: [cooh_num + 2, 123]})
                else:
                    O_z = C_z + bond_O
                    cooh_atoms_dict.update({cooh_num + 1: [C_x, C_y, O_z]})
                    cooh_types_dict.update({cooh_num + 1: 9})
                    cooh_bonds_dict.update({cooh_num + 1: [cooh_num, 123]})
                    OH_z = C_z - bond_O
                    cooh_atoms_dict.update({cooh_num + 2: [C_x, C_y, OH_z]})
                    cooh_types_dict.update({cooh_num + 2: 10})
                    cooh_bonds_dict.update({cooh_num + 2: [cooh_num, 123]})
                    H_z = C_z - bond_O - bond_H
                    cooh_atoms_dict.update({cooh_num + 3: [C_x, C_y, H_z]})
                    cooh_types_dict.update({cooh_num + 3: 5})
                    cooh_bonds_dict.update({cooh_num + 3: [cooh_num + 2, 123]})
                cooh_num += 4
        used_atoms.append(atom)
        used_atoms.append(atom2)

outputfile = open(sys.argv[2], "w")

outputfile.write('''LAMMPS data file via write_data, version 7 Jan 2022, timestep = 1

%d atoms
11 atom types
%d bonds
14 bond types
0 angles
36 angle types
0 dihedrals
103 dihedral types
''' % (H_num - 1 + cooh_num - 1, bonds_num + len(H_bonds_dict) + len(cooh_bonds_dict)))

for line in data[10:17]:
    outputfile.write(line)
###### 1 - c, 2 - co, 3 - coh, 4 - ch, 5 - o, 6 - oh, 7 - ho, 8 - hc

# 1 - c, 2 - co, 3 - coh, 4 - ch, 5 - ch2, 6 - ho, 7 - hc, 8 - h2c, 9 - o, 10 - oh

# 1 - c, 2 - h, 3 - co, 4 - oh, 5 - ho, 6 - o
outputfile.write('''1 12.011000 # c
2 1.008000 # h
3 12.011000 # co
4 15.999000 # oh
5 1.008000 # ho
6 15.999000 # o
7 15.999000 # och2
8 12.011000 # co2
9 15.999000 # o2
10 15.999000 # oh2
11 12.011000 # c2

Pair Coeffs

1 0.070000 3.550000
2 0.030000 2.420000
3 0.066000 3.500000
4 0.170000 3.120000
5 0.000000 0.000000
6 0.140000 2.900000
7 0.170000 3.070000
8 0.105000 3.750000
9 0.210000 2.960000
10 0.170000 3.000000
11 0.070000 3.550000

Bond Coeffs

1   268.0   1.529 # 3 - 3
2   317.0   1.51 # 1 - 3
3   469.0   1.4 # 1 - 1
4   367.0   1.08 # 1 - 2
5   320.0   1.41 # 3 - 6
6   320.0   1.41 # 3 - 4
7   553.0   0.945 # 4 - 5
8   400.0   1.49 # 11 - 8
9   570.0   1.229 # 8 - 9
10  450.0   1.364 # 8 - 10
11  469.0   1.4 # 1 - 11
12  317.0   1.51 # 3 - 11
13  469.0   1.4 # 11 - 11
14  553.0   0.945 # 10 - 5

Angle Coeffs

1       58.350000       112.700000 # 3 - 3 - 3
2       50.000000       109.500000 # 3 - 3 - 4
3       63.000000       114.000000 # 3 - 3 - 1
4       50.000000       109.500000 # 3 - 3 - 6
5       50.000000       109.500000 # 1 - 3 - 6
6       70.000000       130.000000 # 3 - 1 - 3
7       70.000000       120.000000 # 1 - 1 - 3
8       50.000000       109.500000 # 1 - 3 - 4
9       63.000000       120.000000 # 1 - 1 - 1
10      40.000000       109.500000 # 1 - 3 - 1
11      63.000000       120.000000 # 11 - 11 - 1
12      85.000000       120.000000 # 11 - 11 - 8
13      85.000000       120.000000 # 1 - 11 - 8
14      63.000000       120.000000 # 1 - 1 - 11
15      35.000000       120.000000 # 1 - 1 - 2
16      63.000000       120.000000 # 1 - 11 - 1
17      35.000000       117.000000 # 3 - 1 - 2
18      63.000000       114.000000 # 3 - 3 - 11
19      40.000000       109.500000 # 1 - 3 - 11
20      50.000000       109.500000 # 11 - 3 - 6
21      70.000000       119.700000 # 3 - 11 - 8
22      70.000000       120.000000 # 3 - 11 - 1
23      80.000000       120.400000 # 11 - 8 - 9
24      70.000000       120.000000 # 11 - 8 - 10
25      80.000000       121.000000 # 9 - 8 - 10
26      35.000000       120.000000 # 2 - 1 - 11
27      70.000000       120.000000 # 11 - 1 - 3
28      50.000000       109.500000 # 11 - 3 - 4
29      63.000000       120.000000 # 11 - 1 - 11
30      35.000000       113.000000 # 8 - 10 - 5
31      60.000000       109.500000 # 3 - 6 - 3
32      55.000000       108.500000 # 3 - 4 - 5
33      63.000000       120.000000 # 3 - 11 - 11
34      63.000000       120.000000 # 11 - 11 - 11
35      63.000000       120.000000 # 3 - 11 - 3
36      63.000000       120.000000 # 11 - 3 - 11


Dihedral Coeffs

1       1.300000        -0.050000       0.200000        0.000000 # 3 3 3 3
2       1.300000        -0.050000       0.200000        0.000000 # 3 3 3 6
3       -1.552000       0.000000        0.000000        0.000000 # 4 3 3 3
4       4.319000        0.000000        0.000000        0.000000 # 4 3 3 6
5       1.711000        -0.500000       0.663000        0.000000 # 3 3 3 1
6       -4.344000       -1.714000       0.000000        0.000000 # 6 3 3 1
7       2.817000        -0.169000       0.543000        0.000000 # 3 3 1 3
8       1.711000        -0.500000       0.663000        0.000000 # 6 3 1 3
9       -4.344000       -1.714000       0.000000        0.000000 # 1 3 3 1
10      0.000000        0.000000        0.366000        0.000000 # 1 3 3 4
11      0.000000        0.000000        0.000000        0.000000 # 3 3 1 1
12      0.000000        0.000000        0.000000        0.000000 # 4 3 1 1
13      1.711000        -0.500000       0.663000        0.000000 # 4 3 1 3
14      -0.550000       0.000000        0.000000        0.000000 # 6 3 3 6
15      9.508000        0.000000        0.000000        0.000000 # 4 3 3 4
16      0.000000        7.250000        0.000000        0.000000 # 3 1 1 1
17      0.000000        7.250000        0.000000        0.000000 # 3 1 1 3
18      0.000000        7.250000        0.000000        0.000000 # 1 1 1 1
19      0.000000        0.000000        0.000000        0.000000 # 1 1 3 1
20      0.000000        0.000000        0.000000        0.000000 # 6 3 1 1
21      0.000000        -8.000000       0.000000        0.000000 # 3 1 3 1
22      0.000000        7.250000        0.000000        0.000000 # 1 1 1 11
23      0.000000        7.250000        0.000000        0.000000 # 1 1 11 11
24      0.000000        7.250000        0.000000        0.000000 # 1 1 11 8
25      0.000000        7.250000        0.000000        0.000000 # 1 1 1 2
26      0.000000        7.250000        0.000000        0.000000 # 11 1 1 2
27      0.000000        7.250000        0.000000        0.000000 # 1 1 11 1
28      -4.344000       -1.714000       0.000000        0.000000 # 6 3 3 11
29      1.711000        -0.500000       0.663000        0.000000 # 3 3 3 11
30      0.000000        0.000000        0.468000        0.000000 # 3 3 1 2
31      0.000000        0.000000        0.000000        0.000000 # 11 3 1 1
32      0.000000        -8.000000       0.000000        0.000000 # 11 3 1 2
33      0.000000        0.000000        0.468000        0.000000 # 6 3 1 2
34      0.500000        0.000000        0.000000        0.000000 # 3 3 11 8
35      0.000000        0.000000        0.000000        0.000000 # 3 3 11 1
36      0.500000        0.000000        0.000000        0.000000 # 1 3 11 8
37      0.000000        0.000000        0.000000        0.000000 # 1 3 11 1
38      0.500000        0.000000        0.000000        0.000000 # 6 3 11 8
39      0.000000        0.000000        0.000000        0.000000 # 6 3 11 1
40      0.900000        0.230000        -0.505000       0.000000 # 3 11 8 9
41      0.900000        0.230000        -0.505000       0.000000 # 3 11 8 10
42      0.000000        2.100000        0.000000        0.000000 # 1 11 8 9
43      0.000000        2.100000        0.000000        0.000000 # 1 11 8 10
44      0.000000        2.100000        0.000000        0.000000 # 11 11 8 9
45      0.000000        2.100000        0.000000        0.000000 # 11 11 8 10
46      0.000000        7.250000        0.000000        0.000000 # 3 11 1 1
47      0.000000        7.250000        0.000000        0.000000 # 3 11 1 2
48      0.000000        7.250000        0.000000        0.000000 # 8 11 1 2
49      0.000000        7.250000        0.000000        0.000000 # 3 1 1 2
50      0.000000        7.250000        0.000000        0.000000 # 3 1 1 11
51      0.000000        -8.000000       0.000000        0.000000 # 1 3 1 2
52      0.000000        0.000000        0.468000        0.000000 # 4 3 1 2
53      0.000000        0.000000        0.000000        0.000000 # 6 3 1 11
54      0.000000        0.000000        0.000000        0.000000 # 3 3 1 11
55      0.000000        7.250000        0.000000        0.000000 # 3 1 11 3
56      0.000000        0.000000        0.366000        0.000000 # 4 3 3 11
57      0.000000        0.000000        0.000000        0.000000 # 4 3 11 1
58      0.000000        0.000000        0.000000        0.000000 # 4 3 1 11
59      0.000000        7.250000        0.000000        0.000000 # 3 1 11 1
60      0.000000        7.250000        0.000000        0.000000 # 3 1 11 8
61      0.000000        0.000000        0.000000        0.000000 # 1 3 1 11
62      0.000000        7.250000        0.000000        0.000000 # 2 1 1 2
63      -4.344000       -1.714000       0.000000        0.000000 # 1 3 3 11
64      0.000000        7.250000        0.000000        0.000000 # 2 1 11 1
65      0.000000        7.250000        0.000000        0.000000 # 11 1 11 1
66      0.000000        7.250000        0.000000        0.000000 # 11 1 11 8
67      0.000000        7.250000        0.000000        0.000000 # 11 1 11 11
68      0.000000        7.250000        0.000000        0.000000 # 1 11 11 1
69      0.000000        7.250000        0.000000        0.000000 # 1 11 11 8
70      4.000000        5.000000        0.000000        0.000000 # 11 8 10 5
71      0.000000        5.000000        0.000000        0.000000 # 9 8 10 5
72      0.650000        -0.250000       0.670000        0.000000 # 3 3 6 3
73      -0.521000       -2.018000       1.996000        0.000000 # 1 3 6 3
74      -0.356000       -0.174000       0.492000        0.000000 # 3 3 4 5
75      -0.900000       0.000000        0.000000        0.000000 # 1 3 4 5
76      -0.521000       -2.018000       1.996000        0.000000 # 11 3 6 3
77      -0.900000       0.000000        0.000000        0.000000 # 11 3 4 5
78      0.000000        0.000000        0.000000        0.000000 # 11 3 1 3
79      0.000000        0.000000        0.366000        0.000000 # 4 3 11 8
80      0.000000        7.250000        0.000000        0.000000 # 11 1 1 11
81      0.000000        7.250000        0.000000        0.000000 # 11 11 1 3
82      0.000000        7.250000        0.000000        0.000000 # 8 11 11 8
83      0.000000        7.250000        0.000000        0.000000 # 11 1 3 11
84      0.000000        7.250000        0.000000        0.000000 # 3 11 1 11
85      0.000000        7.250000        0.000000        0.000000 # 11 11 1 2
86      0.000000        7.250000        0.000000        0.000000 # 11 11 3 3
87      0.000000        0.000000        0.366000        0.000000 # 4 3 11 11
88      0.000000        7.250000        0.000000        0.000000 # 3 11 11 3
89      0.000000        7.250000        0.000000        0.000000 # 3 11 11 8
90      0.000000        7.250000        0.000000        0.000000 # 1 3 11 11
91      0.000000        0.000000        0.366000        0.000000 # 6 3 11 11
92      0.000000        7.250000        0.000000        0.000000 # 3 11 3 3
93      0.000000        0.000000        0.366000        0.000000 # 4 3 11 3
94      0.000000        7.250000        0.000000        0.000000 # 3 11 3 1
95      0.000000        7.250000        0.000000        0.000000 # 11 3 11 8
96      0.000000        7.250000        0.000000        0.000000 # 11 3 11 3
97      0.000000        7.250000        0.000000        0.000000 # 11 3 11 1
98      0.000000        7.250000        0.000000        0.000000 # 1 11 11 11
99      0.000000        7.250000        0.000000        0.000000 # 11 11 11 8
100      0.000000        7.250000        0.000000        0.000000 # 3 11 11 1
101      0.000000        7.250000        0.000000        0.000000 # 3 11 11 11
102      0.000000        7.250000        0.000000        0.000000 # 11 3 3 11
103      -4.344000       -1.714000       0.000000        0.000000 # 3 11 3 6

Atoms # full

''')

charges_dict = dict()
charges = [0.000, 0.200, 0.265, -0.115, -0.120, 0.418, 0.115, 0.060, -0.400, -0.683, -0.115, -0.115]
#charges = [0.000, 0.115, 0.2, -0.683, 0.418, -0.4]
for i in range(11):
    charges_dict.update({i+1:charges[i]})
#print(charges_dict[2])

atomsIDs = [i+1 for i in range(atoms_num)]
for atom in atomsIDs:
    #print(new_atoms_type[atom])
    #print(charges_dict[new_atoms_type[atom]])
    #print(atom)
    outputfile.write("%d 1 %d %f %f %f %f 0 0 0\n" % (atom, type_change_dict[new_atoms_type[atom]], charges_dict[new_atoms_type[atom]], atoms_coord_dict[atom][0], atoms_coord_dict[atom][1], atoms_coord_dict[atom][2]))
for atom in H_atoms:
    outputfile.write("%d 1 %d %f %f %f %f 0 0 0\n" % (atom, type_change_dict[new_atoms_type[atom]], charges_dict[new_atoms_type[atom]], H_atoms[atom][0], H_atoms[atom][1], H_atoms[atom][2]))
for atom in cooh_atoms_dict:
    outputfile.write("%d 1 %d %f %f %f %f 0 0 0\n" % (H_num - 1 + atom, cooh_types_dict[atom], cooh_types_charges[cooh_types_dict[atom]], cooh_atoms_dict[atom][0], cooh_atoms_dict[atom][1], cooh_atoms_dict[atom][2]))
# atom types 1 - c, 2 - h, 3 - co, 4 - oh, 5 - ho, 6 - o
# bond types 1 CO-CO 2 C-CO 3 C-C 4 C-H 5 CO-O 6 CO-OH 7 OH-H
# 1 - 1, 2 - 2, 3 - 3, 4 - 7, 5 - 14, 6 - 15, 7 - 16


bond_types_dict = dict()
bond_types_dict.update({"3, 3" : 1, "1, 3" : 2, "3, 1" : 2, "1, 1" : 3, "1, 2" : 4, "2, 1" : 4, "3, 6" : 5, "6, 3" : 5, "3, 4" : 6, "4, 3" : 6, "4, 5" : 7, "5, 4" : 7})
bond_types_dict.update({"11, 8" : 8, "8, 11" : 8, "8, 9" : 9, "9, 8" : 9, "8, 10" : 10, "10, 8" : 10, "1, 11" : 11, "11, 1" : 11, "3, 11" : 12, "11, 3" : 12, "11, 11" : 13, "10, 5" : 14, "5, 10" : 14})
outputfile.write('''
Bonds

''')
for line in data[bond_pos+2:bond_pos+2+bonds_num]:
    bond_id = int(line.split()[0])
    atom1 = int(line.split()[2])
    atom2 = int(line.split()[3])
    atom_type1 = type_change_dict[new_atoms_type[atom1]]
    atom_type2 = type_change_dict[new_atoms_type[atom2]]
    bond_tmp = str(atom_type1) + ", " + str(atom_type2)
    outputfile.write("%d %d %d %d\n" % (bond_id, bond_types_dict[bond_tmp], atom1, atom2))

b_count = bonds_num+1
for atom in H_bonds_dict:
    outputfile.write("%d 4 %d %d\n" % (b_count, H_bonds_dict[atom], atom))
    b_count += 1
for atom in cooh_bonds_dict:
    if len(cooh_bonds_dict[atom]) > 1:
        atom_type1 = cooh_types_dict[atom]
        atom_type2 = cooh_types_dict[cooh_bonds_dict[atom][0]]
        bond_tmp = str(atom_type1) + ", " + str(atom_type2)
        outputfile.write("%d %d %d %d\n" % (b_count, bond_types_dict[bond_tmp], cooh_bonds_dict[atom][0] + H_num - 1, atom + H_num - 1))
    else:
        atom_type1 = cooh_types_dict[atom]
        atom_type2 = new_atoms_type[cooh_bonds_dict[atom][0]]
        bond_tmp = str(atom_type1) + ", " + str(atom_type2)
        outputfile.write("%d %d %d %d\n" % (b_count, bond_types_dict[bond_tmp], cooh_bonds_dict[atom][0], atom + H_num - 1))
    b_count += 1
outputfile.close()









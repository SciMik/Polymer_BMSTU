import sys

def create_bonds(bonds_d, atoms_a, new_bonds_d, IDs, coords):
    for atom in IDs:
        if atom not in atoms_a:
            atoms_a.append(atom)
            xyz = coords[atom]
            #if atom == 193:
            #    print(atom, xyz)
            for atom_tmp in bonds_d[atom]:
                if atom_tmp not in atoms_a:
                    xyz_tmp = coords[atom_tmp]
                    if xyz[1] == xyz_tmp[1]:
                        atoms_a.append(atom_tmp)
                        new_bonds_d.update({atom:atom_tmp})
                        new_bonds_d.update({atom_tmp:atom})

data = open(sys.argv[1], "r").readlines()

type_dict = dict()

atoms_num = int(data[2].split()[0])
atoms_IDs = [i+1 for i in range(atoms_num)]
bonds_num = int(data[4].split()[0])
box = [[float(data[7].split()[0]), float(data[7].split()[1])],
           [float(data[8].split()[0]), float(data[8].split()[1])],
           [float(data[9].split()[0]), float(data[9].split()[1])]]
lbox = [box[0][1] - box[0][0], box[1][1] - box[1][0], box[2][1] - box[2][0]]


bonds_pos = 0
atom_pos = 0
for i, line in enumerate(data):
    try:
        if line.split()[0] == "Bonds":
            bonds_pos = i
            break
        elif line.split()[0] == "Atoms":
            atom_pos = i
    except:
        nothong = 1

atom_dict = dict()
atom_coords_dict = dict()
for line in data[atom_pos+2:atom_pos+2+atoms_num]:
    type_dict.update({int(line.split()[0]):int(line.split()[1])})
    atom_coords_dict.update({int(line.split()[0]):[float(line.split()[2]), float(line.split()[3]), float(line.split()[4])]})

#print(atom_coords_dict[193])
#print(atom_coords_dict[194])

spec_bonds_dict = dict() #Bonds to add
bonds_dict = dict()
for line in data[bonds_pos+2:bonds_pos+2+bonds_num]:
    atom1 = int(line.split()[2])
    atom2 = int(line.split()[3])
    try:
        bonds_dict[atom1].append(atom2)
    except:
        bonds_dict.update({atom1:[atom2]})
    try:
        bonds_dict[atom2].append(atom1)
    except:
        bonds_dict.update({atom2:[atom1]})
    if type_dict[atom1] == 1 and type_dict[atom2] == 3:
        spec_bonds_dict.update({atom1:atom2})
    elif type_dict[atom1] == 3 and type_dict[atom2] == 1:
        spec_bonds_dict.update({atom1:atom2})

#print(bonds_dict[81])
#print(len(spec_bonds_dict))
#spec_data = open("SpecBonds.dat", "r").readlines()
spec_data = []
for line in spec_data:
    if len(spec_data) != 0:
        atom1 = int(line.split()[0])
        atom2 = int(line.split()[1])
        spec_bonds_dict.update({atom1:atom2})

new_bonds_dict = dict()
atoms_arr = []

for atom in spec_bonds_dict:
    atoms_arr.append(atom)
    atoms_arr.append(spec_bonds_dict[atom])
    new_bonds_dict.update({atom:spec_bonds_dict[atom]})
    new_bonds_dict.update({spec_bonds_dict[atom]:atom})
create_bonds(bonds_dict, atoms_arr, new_bonds_dict, atoms_IDs, atom_coords_dict)
#string1 = ""
#for atom in atoms_arr:
#    string1 += "ParticleIdentifier == " + str(atom) + " || "
#print(string1)

outfile = open(sys.argv[2], "w")
for line in data[:bonds_pos+2]:
    outfile.write(line)
counter = 1
used_atoms = []
for atom in new_bonds_dict:
    if atom not in used_atoms:
        used_atoms.append(atom)
        used_atoms.append(new_bonds_dict[atom])
        outfile.write("%d 1 %d %d\n" % (counter, atom, new_bonds_dict[atom]))
        counter += 1
outfile.close()

outfile = open("double_C_bonds.dat", "w")
for atom in new_bonds_dict:
    outfile.write("%d %d\n" % (atom, new_bonds_dict[atom]))
outfile.close()
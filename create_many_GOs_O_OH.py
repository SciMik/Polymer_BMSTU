import sys
import random as rand
import os

def create_go(data, outfile_name, file_id):

    z_alone_O = 1.2
    bond_O = 1.4
    bond_H = 1
    OH_groups_num = int(sys.argv[5])
    O_groups_num = int(sys.argv[6])

    atoms_num = int(data[2].split()[0])
    bonds_num = int(data[4].split()[0])

    box = [[float(data[7].split()[0]), float(data[7].split()[1])],
           [float(data[8].split()[0]), float(data[8].split()[1])],
           [float(data[9].split()[0]), float(data[9].split()[1])]]
    lbox = [box[0][1] - box[0][0], box[1][1] - box[1][0], box[2][1] - box[2][0]]
    #print(atoms_num, bonds_num)
    #print(box, lbox)

    atoms_coord = dict()

    atoms_gr, bonds_gr = 0, 0

    for i, line in enumerate(data):
        try:
            if line.split()[0] == 'Atoms':
                #print(line)
                atoms_gr = i
            elif line.split()[0] == 'Bonds':
                #print(line)
                bonds_gr = i
        except:
            nothing = 1

    for line in data[atoms_gr + 2:atoms_gr + 2 + atoms_num]:
        atom_id = int(line.split()[0])
        x = float(line.split()[2])
        y = float(line.split()[3])
        z = float(line.split()[4])
        #x = float(line.split()[3])
        #y = float(line.split()[4])
        #z = float(line.split()[5])
        atoms_coord.update({atom_id: [x, y, z]})

    bonds_dict = dict()
    # print(data[bonds_gr+2],data[bonds_gr+1+bonds_num])
    bonds_C_num = 0
    for line in data[bonds_gr + 2:bonds_gr + 2 + bonds_num]:
        atom_IDs = [int(line.split()[2]), int(line.split()[3])]
        bonds_C_num += 1
        try:
            bonds_dict[atom_IDs[0]].append(atom_IDs[1])
        except:
            bonds_dict.update({atom_IDs[0]: [atom_IDs[1]]})
        try:
            bonds_dict[atom_IDs[1]].append(atom_IDs[0])
        except:
            bonds_dict.update({atom_IDs[1]: [atom_IDs[0]]})

    # print("Length of bonds_dict is ", len(bonds_dict))
    # print(bonds_dict)

    OH_array = []
    added = 0
    while added != OH_groups_num:
        number = rand.randint(1, atoms_num)
#        if number not in OH_array and atoms_coord[number][0] < 12.0:
        if number not in OH_array:
            OH_array.append(number)
            added += 1

    O_array = []
    O_array_all = []

    added = 0
    while added != O_groups_num:
        number = rand.randint(1, atoms_num)
#        if number not in OH_array and number not in O_array_all and atoms_coord[number][0] < 12.0:
        if number not in OH_array and number not in O_array_all:
            checker = 0
            id_ban = []
            for id in bonds_dict[number]:
                if id in OH_array or id in O_array_all:
                    id_ban.append(id)
            if len(id_ban) < len(bonds_dict[number]):
                var = rand.randint(0, len(bonds_dict[number])-1)
                #print(bonds_dict[number], var)
                number2 = bonds_dict[number][var]
                while number2 in id_ban:
                    var = rand.randint(0, len(bonds_dict[number])-1)
                    number2 = bonds_dict[number][var]
                O_array.append([number, number2])
                O_array_all.append(number)
                O_array_all.append(number2)
                added += 1

    # diff_data = open(sys.argv[3], "r").readlines()

    # print(OH_array)
    # print(O_array)
    # for elem in OH_array:
    # for elem2 in O_array:
    #    if elem2[0] in OH_array or elem2[1] in OH_array:
    #        print("Error carbon ID is ", elem2)

    O_coords = dict()
    H_coords = dict()
    alone_O_coords = dict()
    counter = atoms_num
    mol_num = 1
    # print("Hydroxide\n")
    bonds_OH = [] # [[bondtype, atom_1, atom_2]]
    for num in OH_array:
        number = rand.randint(1, 2)
        # print(num, atoms_coord[num])
        # print(counter, num, atoms_coord[num])
        if number == 1:
            mol_num += 1
            z_O = atoms_coord[num][2] + bond_O
            z_H = atoms_coord[num][2] + bond_O + bond_H
            counter += 1
            O_coords.update({counter: [atoms_coord[num][0], atoms_coord[num][1], z_O, mol_num]})
            H_coords.update({counter + OH_groups_num: [atoms_coord[num][0], atoms_coord[num][1], z_H, mol_num]})
            bonds = [[2,num,counter],[3,counter,counter + OH_groups_num]]
            for bond_OH in bonds:
                bonds_OH.append(bond_OH)
        elif number == 2:
            mol_num += 1
            z_O = atoms_coord[num][2] - bond_O
            z_H = atoms_coord[num][2] - bond_O - bond_H
            counter += 1
            O_coords.update({counter: [atoms_coord[num][0], atoms_coord[num][1], z_O, mol_num]})
            H_coords.update({counter + OH_groups_num: [atoms_coord[num][0], atoms_coord[num][1], z_H, mol_num]})
            bonds = [[2, num, counter], [3, counter, counter + OH_groups_num]]
            for bond_OH in bonds:
                bonds_OH.append(bond_OH)
    # print("Oxygen\n")
    bonds_O = []
    for elem in O_array:
        number = rand.randint(1, 2)
        num = elem[0]
       # print(num, atoms_coord[num])
        num2 = elem[1]
        # print(num2, atoms_coord[num2])
        dist_x = abs(atoms_coord[num][0] - atoms_coord[num2][0])
        if dist_x < 4:
            x_O = (atoms_coord[num][0] + atoms_coord[num2][0]) / 2.0
        else:
            x_O = (atoms_coord[num][0] + atoms_coord[num2][0] + lbox[0]) / 2.0
            if x_O > box[0][1]:
                x_O -= box[0][1]
            elif x_O < box[0][0]:
                x_O += box[0][1]
        dist_y = abs(atoms_coord[num][1] - atoms_coord[num2][1])
        if dist_y < 4:
            y_O = ((atoms_coord[num][1] + atoms_coord[num2][1]) / 2.0)
        else:
            y_O = (atoms_coord[num][1] + atoms_coord[num2][1] + lbox[1]) / 2.0
            if y_O > box[1][1]:
                y_O -= box[1][1]
            elif y_O < box[1][0]:
                y_O += box[1][1]
       # print(num, atoms_coord[num])
        # print(num2, atoms_coord[num2])
        if number == 1:
            mol_num += 1
            z_O = atoms_coord[num][2] + z_alone_O
            counter += 1
            # print(counter + OH_groups_num, [x_O, y_O, z_O], num, atoms_coord[num], num2, atoms_coord[num2])
            # print(counter + OH_groups_num, [x_O, y_O, z_O, mol_num])
            alone_O_coords.update({counter + OH_groups_num: [x_O, y_O, z_O, mol_num]})

        if number == 2:
            mol_num += 1
            z_O = atoms_coord[num][2] - z_alone_O
            counter += 1
            # print(counter + OH_groups_num, [x_O, y_O, z_O], num, atoms_coord[num], num2, atoms_coord[num2])
            # print(counter + OH_groups_num, [x_O, y_O, z_O, mol_num])
            alone_O_coords.update({counter + OH_groups_num: [x_O, y_O, z_O, mol_num]})
        bonds = [[4,num,counter + OH_groups_num], [4,num2,counter + OH_groups_num]]
        for bond_O in bonds:
            bonds_O.append(bond_O)

    #print(counter, OH_groups_num, O_groups_num)
    counter += OH_groups_num
    #counter += O_groups_num
    # counter += O_groups_num
    bonds_num = bonds_C_num + len(bonds_OH) + len(bonds_O)
    output_file = open(outfile_name + "." + str(file_id) + ".data", "w")
    output_file.write(data[0])
    output_file.write(
        "\n%d atoms\n3 atom types\n%d bonds\n4 bond types\n0 angles\n3 angle types\n0 dihedrals\n3 dihedral types\n\n" % (
            counter, bonds_num))
    for i in range(7, 13):
        output_file.write(data[i])
    output_file.write('''1 12.0107
2 1.00794
3 15.9994

Pair Coeffs # lj/cut

1 0.07 3.55
2 0.01 1.2
3 0.21 2.96

Bond Coeffs # harmonic

1 310 1.526
2 232 1.42
3 350 0.946
4 232 1.42

Angle Coeffs # harmonic

1 40 114
2 35 109
3 35 109

Dihedral Coeffs # opls

1 1.74 -0.157 0.279 0
2 1.711 -0.5 0.663 0
3 -0.356 -0.174 0.492 0
    
Atoms # molecular
    
''')
    for atom1 in atoms_coord:
        output_file.write(
            "%d 1 1 %f %f %f 0 0 0\n" % (atom1, atoms_coord[atom1][0], atoms_coord[atom1][1], atoms_coord[atom1][2]))
    for atom2 in O_coords:
        output_file.write("%d %d 3 %f %f %f 0 0 0\n" % (
        atom2, O_coords[atom2][3], O_coords[atom2][0], O_coords[atom2][1], O_coords[atom2][2]))
    for atom3 in H_coords:
        output_file.write("%d %d 2 %f %f %f 0 0 0\n" % (
        atom3, H_coords[atom3][3], H_coords[atom3][0], H_coords[atom3][1], H_coords[atom3][2]))
    for atom4 in alone_O_coords:
        # print(atom4, alone_O_coords[atom4])
        output_file.write("%d %d 3 %f %f %f 0 0 0\n" % (
        atom4, alone_O_coords[atom4][3], alone_O_coords[atom4][0], alone_O_coords[atom4][1], alone_O_coords[atom4][2]))
    output_file.write("\nBonds\n\n")
    bond_counter = 1
    for line in data[bonds_gr + 2:bonds_gr + 2 + bonds_num]:
        output_file.write(line)
        bond_counter += 1
    for elem in bonds_OH:
        output_file.write("%d %d %d %d\n" % (bond_counter, elem[0], elem[1], elem[2]))
        bond_counter += 1
    for elem in bonds_O:
        output_file.write("%d %d %d %d\n" % (bond_counter, elem[0], elem[1], elem[2]))
        bond_counter += 1


    main_arrays = [OH_array, O_array_all]
    #main_arrays = [OH_array]
    return main_arrays

path = "./"

data = open(sys.argv[1], "r").readlines()
outfile_name = sys.argv[2]
max_ID = 0
for name in os.listdir(path):
    #print(name)
    if name.split(".")[0] == outfile_name and name.split(".")[1].isnumeric():
            if int(name.split(".")[1]) > max_ID:
                max_ID = int(name.split(".")[1])
    #except:
    #    print("No names ", outfile_name, " found.")
arrays = dict()
id = 0
#arrays.update({0:[]})
try:
    arr_data = open(sys.argv[3], "r").readlines()
    for line in arr_data:
        tmp_arr1 = line.split(',')[0]
        tmp_arr2 = line.split(',')[1]
        arrays.update({id:[]})
        tmp_arr_big = []
        for elem in tmp_arr1.split():
            tmp_arr_big.append(int(elem))
        arrays[id].append(tmp_arr_big)
        tmp_arr_big = []
        for elem in tmp_arr2.split():
            tmp_arr_big.append(int(elem))
        arrays[id].append(tmp_arr_big)
        id += 1
except:
    print("No data")


arr_file = open(sys.argv[3], "a")
file_num = int(sys.argv[4])
i = 0

while i < file_num:
    file_id = i+1+max_ID
    #print(file_id)
    tmp_arrays = create_go(data, outfile_name, file_id)
    flag = 1
    if id == 0:
        arrays.update({id:tmp_arrays})
        id += 1
        i += 1
        for elem in tmp_arrays[0]:
            arr_file.write(str(elem) + " ")
        arr_file.write(",")
        for elem in tmp_arrays[1]:
            arr_file.write(str(elem) + " ")
        arr_file.write("\n")
        flag = 0
    else:
        for elem in arrays:
            checker = 0
            for num in arrays[elem][0]:
                if num not in tmp_arrays[0]:
                    checker = 1
                    break
            if checker == 0:
                for num in arrays[elem][1]:
                    if num not in tmp_arrays[1]:
                        checker = 1
                        break
            if checker == 0:
                flag = 0

    if flag:
        arrays.update({id: tmp_arrays})
        for elem in tmp_arrays[0]:
            arr_file.write(str(elem) + " ")
        arr_file.write(",")
        for elem in tmp_arrays[1]:
            arr_file.write(str(elem) + " ")
        arr_file.write("\n")
        i += 1
arr_file.close()


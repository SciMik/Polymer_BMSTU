import sys
import random as rand

data = open(sys.argv[1], "r").readlines()

atoms_num = int(data[2].split()[0])
bonds_num = int(data[4].split()[0])

bonds_list = [i for i in range(1,bonds_num+1)]

atoms_bond = dict()
atoms_bonded = dict()

atoms_type = dict()

bonds_pos = 0
atoms_pos = 0
for i, line in enumerate(data):
    try:
        if line.split()[0] == "Atoms":
            atoms_pos = i
        if line.split()[0] == "Bonds":
            bonds_pos = i
    except:
        nothing = 1

for line in data[atoms_pos+2:atoms_pos+2+atoms_num]:
    atomID = int(line.split()[0])
    atomType = int(line.split()[2])
    atoms_type.update({atomID:atomType})

for line in data[bonds_pos+2:bonds_pos+2+bonds_num]: ####BONDS
    var1 = int(line.split()[0])
    var2 = int(line.split()[1])
    atoms_bond.update({var1:[set(), var2]})
    atoms_bond[int(line.split()[0])][0].add(int(line.split()[2]))
    atoms_bond[int(line.split()[0])][0].add(int(line.split()[3]))
    try:
        atoms_bonded[int(line.split()[2])].append(int(line.split()[3]))
    except:
        atoms_bonded.update({int(line.split()[2]) : [int(line.split()[3])]})
    try:
        atoms_bonded[int(line.split()[3])].append(int(line.split()[2]))
    except:
        atoms_bonded.update({int(line.split()[3]) : [int(line.split()[2])]})

angles = []
angles_num = 0
angle_type_dict = dict()
#angle_type_dict.update({"3, 3, 3" : 1, "3, 3, 4" : 2,"3, 3, 1" : 3, "3, 3, 6" : 4, "1, 3, 6" : 5, "3, 1, 3" : 6, "1, 1, 3" : 7, "1, 3, 4" : 8})
#angle_type_dict.update({"1, 1, 1" : 9, "1, 3, 1" : 10, "1, 1, 2" : 11, "3, 1, 2" : 12, "3, 6, 3" : 13, "3, 4, 5" : 14, "2, 1, 2" : 15, "6, 3, 4" : 16})
angle_type_dict.update({"3, 3, 3" : 1, "3, 3, 4" : 2, "3, 3, 1" : 3, "3, 3, 6" : 4, "1, 3, 6" : 5, "3, 1, 3" : 6, "1, 1, 3" : 7, "1, 3, 4" : 8, "1, 1, 1" : 9, "1, 3, 1" : 10, "11, 11, 1" : 11, "11, 11, 8" : 12, "1, 11, 8" : 13, "1, 1, 11" : 14, "1, 1, 2" : 15, "1, 11, 1" : 16, "3, 1, 2" : 17, "3, 3, 11" : 18, "1, 3, 11" : 19, "11, 3, 6" : 20, "3, 11, 8" : 21, "3, 11, 1" : 22, "11, 8, 9" : 23, "11, 8, 10" : 24, "9, 8, 10" : 25, "2, 1, 11" : 26, "11, 1, 3" : 27, "11, 3, 4" : 28, "11, 1, 11" : 29, "8, 10, 5" : 30, "3, 6, 3" : 31, "3, 4, 5" : 32, "3, 11, 11" : 33, "11, 11, 11" : 34, "3, 11, 3" : 35, "11, 3, 11" : 36})
#angle_type_list = []
#for elem in angle_type_dict:
#    types_tmp = [int(elem.split(",")[0]), int(elem.split(",")[1]), int(elem.split(",")[2])]
#    angle_type_list.append(types_tmp)
for i in range(len(bonds_list)):
    counter = 0 ### max check bonds
    for j in range(i+1,len(bonds_list)):
        if atoms_bond[bonds_list[i]][0].intersection(atoms_bond[bonds_list[j]][0]):
            middle_atom = int(list(atoms_bond[bonds_list[i]][0].intersection(atoms_bond[bonds_list[j]][0]))[0])
            angles_num += 1
            counter += 1
            angle_set = list()
            for elem in atoms_bond[bonds_list[i]][0]:
                if elem != middle_atom:
                    angle_set.append(elem)
            angle_set.append(middle_atom)
            for elem in atoms_bond[bonds_list[j]][0]:
                if elem != middle_atom:
                    angle_set.append(elem)
                #angle_set = list(angle_set)
            types = [atoms_type[angle_set[0]], atoms_type[angle_set[1]], atoms_type[angle_set[2]]]
            types_string = str(types[0]) + ", " + str(types[1]) + ", " + str(types[2])
            types_string_reverse = str(types[2]) + ", " + str(types[1]) + ", " + str(types[0])
            error_check = True
            for elem in angle_type_dict:
                if types_string == elem:
                    angleType = angle_type_dict[types_string]
                    angles.append([angles_num, angleType, angle_set[0], angle_set[1], angle_set[2]])
                    error_check = False
                elif types_string_reverse == elem:
                    angleType = angle_type_dict[types_string_reverse]
                    angles.append([angles_num, angleType, angle_set[2], angle_set[1], angle_set[0]])
                    error_check = False
            if error_check:
                print(angle_set)
                print(types_string)
        if counter == 4:
            break

angles_list = [i for i in range(1,angles_num+1)]

angles_dict = dict()
angles_dict_list = dict()

for angle in angles:
    angles_dict.update({angle[0]:[set(), angle[1]]})
    angles_dict_list.update({angle[0]: [angle[2], angle[3], angle[4]]})
    angles_dict[angle[0]][0].add(angle[2])
    angles_dict[angle[0]][0].add(angle[3])
    angles_dict[angle[0]][0].add(angle[4])

dihedrals_type_dict = dict()
#dihedrals_type_dict.update({"3, 3, 3, 3" : 1, "3, 3, 3, 6" : 2, "4, 3, 3, 3" : 3, "4, 3, 3, 6" : 4, "3, 3, 3, 1" : 5, "6, 3, 3, 1" : 6})
#dihedrals_type_dict.update({"3, 3, 1, 3" : 7, "6, 3, 1, 3" : 8, "1, 3, 3, 1" : 9, "1, 3, 3, 4" : 10, "3, 3, 1, 1" : 11, "4, 3, 1, 1" : 12})
#dihedrals_type_dict.update({ "4, 3, 1, 3" : 13, "6, 3, 3, 6" : 14, "4, 3, 3, 4" : 15, "3, 1, 1, 1" : 16, "3, 1, 1, 3" : 17, "1, 1, 1, 1" : 18})
#dihedrals_type_dict.update({"1, 1, 3, 1" : 19, "6, 3, 1, 1" : 20, "3, 1, 3, 1" : 21, "1, 1, 1, 2" : 22, "3, 3, 1, 2" : 23, "6, 3, 1, 2" : 24})
#dihedrals_type_dict.update({"3, 1, 1, 2" : 25, "1, 3, 1, 2" : 26, "4, 3, 1, 2" : 27, "2, 1, 1, 2" : 28, "3, 3, 6, 3" : 29, "1, 3, 6, 3" : 30, "3, 3, 4, 5" : 31, "1, 3, 4, 5" : 32})
#dihedrals_type_dict.update({"6, 3, 4, 5" : 33})

dihedrals_type_dict.update({"3, 3, 3, 3" : 1, "3, 3, 3, 6" : 2, "4, 3, 3, 3" : 3, "4, 3, 3, 6" : 4, "3, 3, 3, 1" : 5, "6, 3, 3, 1" : 6, "3, 3, 1, 3" : 7, "6, 3, 1, 3" : 8, "1, 3, 3, 1" : 9, "1, 3, 3, 4" : 10, "3, 3, 1, 1" : 11, "4, 3, 1, 1" : 12, "4, 3, 1, 3" : 13, "6, 3, 3, 6" : 14, "4, 3, 3, 4" : 15, "3, 1, 1, 1" : 16, "3, 1, 1, 3" : 17, "1, 1, 1, 1" : 18, "1, 1, 3, 1" : 19, "6, 3, 1, 1" : 20, "3, 1, 3, 1" : 21, "1, 1, 1, 11" : 22, "1, 1, 11, 11" : 23, "1, 1, 11, 8" : 24, "1, 1, 1, 2" : 25, "11, 1, 1, 2" : 26, "1, 1, 11, 1" : 27, "6, 3, 3, 11" : 28, "3, 3, 3, 11" : 29, "3, 3, 1, 2" : 30, "11, 3, 1, 1" : 31, "11, 3, 1, 2" : 32, "6, 3, 1, 2" : 33, "3, 3, 11, 8" : 34, "3, 3, 11, 1" : 35, "1, 3, 11, 8" : 36, "1, 3, 11, 1" : 37, "6, 3, 11, 8" : 38, "6, 3, 11, 1" : 39, "3, 11, 8, 9" : 40, "3, 11, 8, 10" : 41, "1, 11, 8, 9" : 42, "1, 11, 8, 10" : 43, "11, 11, 8, 9" : 44, "11, 11, 8, 10" : 45, "3, 11, 1, 1" : 46, "3, 11, 1, 2" : 47, "8, 11, 1, 2" : 48, "3, 1, 1, 2" : 49, "3, 1, 1, 11" : 50, "1, 3, 1, 2" : 51, "4, 3, 1, 2" : 52, "6, 3, 1, 11" : 53, "3, 3, 1, 11" : 54, "3, 1, 11, 3" : 55, "4, 3, 3, 11" : 56, "4, 3, 11, 1" : 57, "4, 3, 1, 11" : 58, "3, 1, 11, 1" : 59, "3, 1, 11, 8" : 60, "1, 3, 1, 11" : 61, "2, 1, 1, 2" : 62, "1, 3, 3, 11" : 63, "2, 1, 11, 1" : 64, "11, 1, 11, 1" : 65, "11, 1, 11, 8" : 66, "11, 1, 11, 11" : 67, "1, 11, 11, 1" : 68, "1, 11, 11, 8" : 69, "11, 8, 10, 5" : 70, "9, 8, 10, 5" : 71, "3, 3, 6, 3" : 72, "1, 3, 6, 3" : 73, "3, 3, 4, 5" : 74, "1, 3, 4, 5" : 75, "11, 3, 6, 3" : 76, "11, 3, 4, 5" : 77, "11, 3, 1, 3" : 78, "4, 3, 11, 8" : 79, "11, 1, 1, 11" : 80, "11, 11, 1, 3" : 81, "8, 11, 11, 8" : 82, "11, 1, 3, 11" : 83, "3, 11, 1, 11" : 84, "11, 11, 1, 2" : 85, "11, 11, 3, 3" : 86, "4, 3, 11, 11" : 87, "3, 11, 11, 3" : 88, "3, 11, 11, 8" : 89, "1, 3, 11, 11" : 90, "6, 3, 11, 11" : 91, "3, 11, 3, 3" : 92, "4, 3, 11, 3" : 93, "3, 11, 3, 1" : 94, "11, 3, 11, 8" : 95, "11, 3, 11, 3" : 96, "11, 3, 11, 1" : 97, "1, 11, 11, 11" : 98, "11, 11, 11, 8" : 99, "3, 11, 11, 1" : 100, "3, 11, 11, 11" : 101, "11, 3, 3, 11" : 102, "3, 11, 3, 6" : 103})

dihedrals = []
dihedrals_num = 0
checker = 0
for i in range(len(angles_list)):
    counter = 0 ### max check angles
    for j in range(i+1,len(angles_list)):
        if angles_dict[angles_list[i]][0].intersection(angles_dict[angles_list[j]][0]) and len(angles_dict[angles_list[i]][0].intersection(angles_dict[angles_list[j]][0])) == 2:
            middle_elems = list(angles_dict[angles_list[i]][0].intersection(angles_dict[angles_list[j]][0]))
            tail_1 = 0
            tail_2 = 0
            for elem in angles_dict_list[angles_list[i]]:
                if elem not in middle_elems:
                    tail_1 = elem
            for elem in angles_dict_list[angles_list[j]]:
                if elem not in middle_elems:
                    tail_2 = elem
            if checker == 0:
                #print(atoms_bonded[751])
                #print(atoms_bonded[766])
                checker += 1
            flag = True
            for bond_atom in atoms_bonded[tail_1]:
                if bond_atom in atoms_bonded[tail_2]:
                    flag = False
            if flag:
                dihedrals_num += 1
                counter += 1
                dihedral_set = list()
                #if angles_dict_list[angles_list[i]][1] in middle_elems and angles_dict_list[angles_list[i]][2] in middle_elems:
                count_tmp = 0
                for iter, elem in enumerate(angles_dict_list[angles_list[i]]):
                    if elem not in middle_elems:
                        count_tmp = iter
                        #dihedral_set.append(elem)
                if count_tmp == 0:
                    dihedral_set=[angles_dict_list[angles_list[i]][0], angles_dict_list[angles_list[i]][1], angles_dict_list[angles_list[i]][2]]
                elif count_tmp == 2:
                    dihedral_set = [angles_dict_list[angles_list[i]][2], angles_dict_list[angles_list[i]][1],
                                    angles_dict_list[angles_list[i]][0]]
                for elem in angles_dict_list[angles_list[j]]:
                    if elem not in middle_elems:
                        dihedral_set.append(elem)
                #else:
                #    print("Error!")
                #    print(middle_elems)
                #    print(angles_dict_list[angles_list[i]])
                #    print(angles_dict_list[angles_list[j]])
                #dihedral_set = list(dihedral_set)
                types = [atoms_type[dihedral_set[0]], atoms_type[dihedral_set[1]], atoms_type[dihedral_set[2]], atoms_type[dihedral_set[3]]]
                types_string = str(types[0]) + ", " + str(types[1]) + ", " + str(types[2]) + ", " + str(types[3])
                types_string_reverse = str(types[3]) + ", " + str(types[2]) + ", " + str(types[1]) + ", " + str(types[0])
                error_check = True
                for elem in dihedrals_type_dict:
                    if types_string == elem:
                        dihedralType = dihedrals_type_dict[types_string]
                        dihedrals.append([dihedrals_num, dihedralType, dihedral_set[0], dihedral_set[1], dihedral_set[2], dihedral_set[3]])
                        error_check = False
                    elif types_string_reverse == elem:
                        dihedralType = dihedrals_type_dict[types_string_reverse]
                        dihedrals.append([dihedrals_num, dihedralType, dihedral_set[3], dihedral_set[2], dihedral_set[1], dihedral_set[0]])
                        error_check = False
                if error_check:
                    #print(types)
                    #print("Error")
                    print(dihedral_set)
                    print(types_string)
        if counter == 5:
            break



output_file = open(sys.argv[2], "w")
output_file.write(data[0])
output_file.write("\n%d atoms\n11 atom types\n%d bonds\n14 bond types\n%d angles\n36 angle types\n%d dihedrals\n103 dihedral types\n\n" % (atoms_num, bonds_num, angles_num, dihedrals_num))
#for i in range(11,21):
#    output_file.write(data[i])
for line in data[11:]:
    output_file.write(line)
output_file.write("\nAngles\n\n")
for angle in angles:
    output_file.write("%d %d %d %d %d\n" % (angle[0], angle[1], angle[2], angle[3], angle[4]))
output_file.write("\nDihedrals\n\n")
for dihedral in dihedrals:
    output_file.write("%d %d %d %d %d %d\n" % (dihedral[0], dihedral[1], dihedral[2], dihedral[3], dihedral[4], dihedral[5]))

#print(rand.randint(1,100))
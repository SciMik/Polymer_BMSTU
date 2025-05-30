import numpy as np
import sys

# Function to parse LAMMPS data file
def parse_lammps_data(file_name):
    with open(file_name, "r") as f:
        lines = f.readlines()
    print("Start parses")
    # Parse box dimensions
    for i, line in enumerate(lines):
        if "xlo xhi" in line:
            xlo, xhi = float(line.split()[0]), float(line.split()[1])
        elif "ylo yhi" in line:
            ylo, yhi = float(line.split()[0]), float(line.split()[1])
        elif "zlo zhi" in line:
            zlo, zhi = float(line.split()[0]), float(line.split()[1])
        elif "Masses" in line:
            masses_start = i + 2
        elif "Atoms" in line:
            atoms_start = i + 2
        elif "Velocities" in line:
            vels_start = i+2
        elif "Bonds" in line:
            bonds_start = i + 2
        elif "Angles" in line:
            angles_start = i + 2
        elif "Dihedrals" in line:
            dihedrals_start = i + 2
        elif "Pair Coeffs" in line:
            pair_coeffs_start = i + 2
        elif "Bond Coeffs" in line:
            bond_coeffs_start = i + 2
        elif "Angle Coeffs" in line:
            angle_coeffs_start = i + 2
        elif "Dihedral Coeffs" in line:
            dihedral_coeffs_start = i + 2

    # Box dimensions
    box = np.array([[xlo, xhi], [ylo, yhi], [zlo, zhi]])

    # Parse masses
    masses = {}
    atom_type_name = {}
    for line in lines[masses_start:]:
        if not line.strip():
            break
        atom_type, mass, symb, atom_n = line.split()
        atom_name = atom_n.split(",")[0]
        atom_type_name[int(atom_type)] = atom_name
        masses[int(atom_type)] = float(mass)
    print("Masses done")

    # Parse pair coefficients
    pair_coeffs = {}
    for line in lines[pair_coeffs_start:]:
        if not line.strip():
            break
        atom_type, epsilon, sigma, *_ = line.split()
        pair_coeffs[int(atom_type)] = {"epsilon": float(epsilon), "sigma": float(sigma)}
    print("Pair coeffs done")

    # Parse bond coefficients
    bond_coeffs = {}
    for line in lines[bond_coeffs_start:]:
        if not line.strip():
            break
        bond_type, k, r0, *_ = line.split()
        bond_coeffs[int(bond_type)] = {"k": float(k), "r0": float(r0)}
    print("Bond coeffs done")

    # Parse angle coefficients
    angle_coeffs = {}
    for line in lines[angle_coeffs_start:]:
        if not line.strip():
            break
        angle_type, k, theta0, *_ = line.split()
        angle_coeffs[int(angle_type)] = {"k": float(k), "theta0": float(theta0)}
    print("Angle coeffs done")

    # Parse dihedral coefficients
    dihedral_coeffs = {}
    for line in lines[dihedral_coeffs_start:]:
        if not line.strip():
            break
        fields = line.split()
        dihedral_type = int(fields[0])
        C_OPLS = [float(fields[1]), float(fields[2]), float(fields[3]), float(fields[4])]
        # print(C_Fourier)
        coeffs_RB = dict()
        for i in range(6):
            coeffs_RB.update({i: float(0.0)})
        coeffs_RB[0] = C_OPLS[1] + 1 / 2.0 * (C_OPLS[0] + C_OPLS[2])
        coeffs_RB[1] = 1 / 2.0 * (- C_OPLS[0] + 3 * C_OPLS[2])
        coeffs_RB[2] = - C_OPLS[1] + 4 * C_OPLS[3]
        coeffs_RB[3] = -2 * C_OPLS[2]
        coeffs_RB[4] = -4 * C_OPLS[3]
        coeffs_RB[5] = 0.0

        dihedral_coeffs[dihedral_type] = coeffs_RB

    print("Dihedral coeffs done")

    # Parse atoms
    atoms = []
    molecules = []
    molecule_id1 = []
    total_charge = 0.0
    for line in lines[atoms_start:]:
        if not line.strip():
            break
        fields = line.split()
        atom_id = int(fields[0])
        mol_id = int(fields[1])
        if mol_id == 1:
            molecule_id1.append(atom_id)
        molecules.append(mol_id)
        atom_type = int(fields[2])
        charge = float(fields[3])
        total_charge += float(charge)
        x, y, z = map(float, fields[4:7])
        nx, ny, nz = map(float, fields[7:10])
        atoms.append((atom_id, mol_id, atom_type, charge, x, y, z, nx, ny, nz))
        #atoms.append((atom_id, mol_id, atom_type, charge, x, y, z))
    #atoms = np.array(atoms)
    #atoms.sort(axis=0)
    atoms.sort(key=lambda x: x[0])
    #atoms = list(atoms)
    print("Total charge == %f" % (total_charge))
    print("Atoms done")

    # Parse bonds
    bonds = []
    bonds_dict = dict()
    for line in lines[bonds_start:]:
        if not line.strip():
            break
        bond_id, bond_type, atom1, atom2, *_ = line.split()
        bonds.append((int(bond_type), int(atom1), int(atom2)))
        try:
            bonds_dict[int(atom1)].append(int(atom2))
        except:
            bonds_dict.update({int(atom1):[int(atom2)]})
        try:
            bonds_dict[int(atom2)].append(int(atom1))
        except:
            bonds_dict.update({int(atom2):[int(atom1)]})
    bonds.sort(key=lambda x: x[1])
    print("Bonds done")
    try:
        vels = dict()
        for line in lines[vels_start:]:
            if not line.strip():
                break
            atom_id, vx, vy, vz, *_ = line.split()
            vels.update({int(atom_id):[float(vx), float(vy), float(vz)]})
        print("Vels done")
    except:
        vels=None

    # Parse angles
    angles = []
    for line in lines[angles_start:]:
        if not line.strip():
            break
        angle_id, angle_type, atom1, atom2, atom3, *_ = line.split()
        angles.append((int(angle_type), int(atom1), int(atom2), int(atom3)))
    angles.sort(key=lambda x: x[1])
    print("Angles done")

    # Parse dihedrals
    dihedrals = []
    for line in lines[dihedrals_start:]:
        if not line.strip():
            break
        dihedral_id, dihedral_type, atom1, atom2, atom3, atom4, *_ = line.split()
        dihedrals.append((int(dihedral_type), int(atom1), int(atom2), int(atom3), int(atom4)))
    dihedrals.sort(key=lambda x: x[1])
    print("Dihedrals done")

    return {
        "box": box,
        "masses": masses,
        "atom_type_name": atom_type_name,
        "pair_coeffs": pair_coeffs,
        "bond_coeffs": bond_coeffs,
        "angle_coeffs": angle_coeffs,
        "dihedral_coeffs": dihedral_coeffs,
        "atoms": atoms,
        "velocities" : vels,
        "bonds": bonds,
        "angles": angles,
        "dihedrals": dihedrals,
        "molnum": np.max(molecules),
        "bondsDict": bonds_dict,
        "mol1": molecule_id1,
    }

# Function to write GROMACS .gro file
def write_gro(data, NP_name, gro_file):
    with open(gro_file, "w") as f:
        f.write("Generated by LAMMPS to GROMACS converter\n")
        f.write(f"{int(len(data['atoms']))}\n")
        #f.write(f"{int(len(data['atoms']) / data['molnum'] * 2)}\n")
        el_counter = {"H1": [0 for i in range(data['molnum'])], "C": [0 for i in range(data['molnum'])], "Cl": [0 for i in range(data['molnum'])], "H2": [0 for i in range(data['molnum'])]}
        box = data["box"]
        Lx, Ly, Lz = box[0, 1] - box[0, 0], box[1, 1] - box[1, 0], box[2, 1] - box[2, 0]
        for atom in data["atoms"]:
            atom_id, mol_id, atom_type, charge, x, y, z, nx, ny, nz = atom
            #atom_id, mol_id, atom_type, charge, x, y, z = atom
            x = x + int(nx) * Lx
            y = y + int(ny) * Ly
            z = z + int(nz) * Lz
            if data["velocities"] != None:
                vx, vy, vz = data["velocities"][atom_id][0], data["velocities"][atom_id][1], data["velocities"][atom_id][2]
            element = data['atom_type_name'][atom_type]
            #f.write(f"{mol_id:5d}   {element:<5s}{atom_id:5d}{x/10:8.3f}{y/10:8.3f}{z/10:8.3f}{vx*100:8.3f}{vy*100:8.3f}{vz*100:8.3f}\n")
            #if mol_id == 1:
                #f.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n" % (mol_id, "MOL", hex(el_counter[element][mol_id - 1] + 1)[2:] + element, atom_id, x/10, y/10, z/10, vx*100, vy*100, vz*100))
                #f.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n" % (mol_id, "MOL", element, atom_id, x / 10, y / 10, z / 10, vx * 100, vy * 100, vz * 100))
            #    f.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" % (mol_id, "UNL", element, atom_id, x / 10, y / 10, z / 10))
            #    el_counter[element][mol_id - 1] += 1
            #if mol_id == 1:
            if data["velocities"] != None:
                f.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n" % (mol_id, NP_name, element, atom_id, x / 10, y / 10, z / 10, vx * 100, vy * 100, vz * 100))
            else:
                f.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" % (mol_id, NP_name, element, atom_id, x / 10, y / 10, z / 10))
                #f.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n" % (mol_id, "MOL", element, atom_id, x / 10, y / 10, z / 10, vx * 100, vy * 100, vz * 100))
                #f.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n" % (mol_id, "MOL", hex(el_counter[element][mol_id - 1] + 1)[2:] + element, atom_id, x/10, y/10, z/10, vx*100, vy*100, vz*100))
                #el_counter[element][mol_id - 1] += 1


        f.write(f"{(box[0, 1] - box[0, 0])/10:10.5f}{(box[1, 1] - box[1, 0])/10:10.5f}{(box[2, 1] - box[2, 0])/10:10.5f}\n")

# Function to write GROMACS .top file
def write_top(data, NP_name, top_file):
    with open(top_file, "w") as f:
        f.write("[ defaults ]\n")
        f.write("; nbfunc        comb-rule       gen-pairs       fudgeLJ     fudgeQQ\n")
        f.write("1               2               yes             0.5         0.8333\n\n")

        f.write("[ atomtypes ]\n")
        f.write("; name      mass      charge      ptype      sigma       epsilon\n")
        for atom_type, params in data["pair_coeffs"].items():
            element = data['atom_type_name'][atom_type]
            sigma = params["sigma"] * 0.1  # Convert Å to nm
            epsilon = params["epsilon"] * 4.184  # Convert kcal/mol to kJ/mol
            f.write(f"{element:<10s}    {data['masses'][atom_type]:10.6f} 0.0 A {sigma:10.6f} {epsilon:10.6f}\n")
            #f.write()

        f.write("\n[ moleculetype ]\n")
        f.write("; name            nrexcl\n")
        f.write("%s           3\n" % (NP_name))

        f.write("\n[ atoms ]\n")
        f.write("; id    type    resnr    residuename    atomname  cgnr  charge  mass  typeB    chargeB   massB\n")
        counter = 0
        #el_counter = {"H1": 0, "C": 0, "Cl": 0, "H2": 0}
        for atom in data["atoms"]:
            atom_id, mol_id, atom_type, charge, x, y, z, nx, ny, nz = atom
            #atom_id, mol_id, atom_type, charge, x, y, z = atom
            element = data['atom_type_name'][atom_type]
            if mol_id == 1:
                #f.write(f"{atom_id:6d}  {element:<10s}{mol_id:6d}   MOL    {hex(el_counter[element] + 1)[2:] + element:<6s}  1  {charge:10.20f}  {data['masses'][atom_type]:10.4f}\n")
                f.write(f"{atom_id:6d}  {element:<10s}{mol_id:6d}   {NP_name:<6s}    {element:<6s}  {atom_id:6d}  {charge:10.20f} {data['masses'][atom_type]:10.6f} ; qtot {charge:10.20f}\n")
            #f.write(f"{atom_id:6d}  {element:<10s}    1   MOL    {hex(el_counter[element] + 1)[2:] + element:<6s}  1  {charge:10.20f}\n")
                #el_counter[element] += 1
            counter += 1
        print(counter)

        f.write("\n[ bonds ]\n")
        f.write("; ai    aj    funct    c0    c1    c2   c3\n")
        for bond in data["bonds"]:
            bond_type, atom1, atom2 = bond
            if atom1 in data['mol1'] and atom2 in data['mol1']:
                params = data["bond_coeffs"][bond_type]
                r0 = params["r0"] * 0.1  # Convert Å to nm
                k = params["k"] * 418.4  # Convert kcal/mol/Å² to kJ/mol/nm²
                f.write(f"{atom1:6d}{atom2:6d}   1   {r0:10.6f}  {k:10.6f}\n")

        f.write("\n[ pairs ]\n")
        f.write("; ai    aj    funct   c0   c1   c2   c3\n")
        atoms_num = len(data['atoms'])
        for i in range(atoms_num):
            atom_id = i + 1
            if atom_id in data['mol1']:
                for atom_1 in data['bondsDict'][atom_id]:
                    for atom_2 in data['bondsDict'][atom_1]:
                        #if atom_2 > atom_id:
                        #    f.write(f"{atom_id:6d}  {atom_2:6d}   1\n")
                        if atom_2 != atom_id:
                            for atom_3 in data['bondsDict'][atom_2]:
                                if atom_3 != atom_1 and atom_3 > atom_id:
                                    f.write(f"{atom_id:6d}  {atom_3:6d}   1\n")

        f.write("\n[ angles ]\n")
        f.write("; ai    aj    ak    funct    c0    c1   c2   c3\n")
        for angle in data["angles"]:
            angle_type, atom1, atom2, atom3 = angle
            if atom1 in data['mol1'] and atom2 in data['mol1'] and atom3 in data['mol1']:
                params = data["angle_coeffs"][angle_type]
                theta0 = params["theta0"]
                k = params["k"] * 4.184  # Convert kcal/mol/radian² to kJ/mol/radian²
                f.write(f"{atom1:6d}{atom2:6d}{atom3:6d}  1  {theta0:10.6f}   {k:10.6f}\n")

        f.write("\n[ dihedrals ]\n")
        f.write("; ai    aj    ak    al    funct    c0    c1    c2    c3    c4    c5\n")
        for dihedral in data["dihedrals"]:
            dihedral_type, atom1, atom2, atom3, atom4 = dihedral
            if atom1 in data['mol1'] and atom2 in data['mol1'] and atom3 in data['mol1'] and atom4 in data['mol1']:
                coeffs = data["dihedral_coeffs"][dihedral_type]
                f.write(f"{atom1:6d}{atom2:6d}{atom3:6d}{atom4:6d}  3  ")
                # if len(coeffs) > 4:
                #    print("ERROR")
                #    break
                for coeff in coeffs:
                    f.write(f"{coeffs[coeff] * 4.184:10.6f}")
                f.write("\n")

        f.write("\n[ system ]\n")
        f.write("; Name\n")
        f.write("Nanoparticles\n")
        f.write("\n[ molecules ]\n")
        f.write("; Compound            nmols\n")
        #f.write("UNL           2\n\n")
        f.write("%s            %d\n\n" % (NP_name, data['molnum']))

# Main script
if __name__ == "__main__":
    lammps_file = sys.argv[1]
    NP_name = sys.argv[2]
    gro_file = sys.argv[3] + ".gro"
    top_file = sys.argv[3] + ".top"

    data = parse_lammps_data(lammps_file)
    write_gro(data, NP_name, gro_file)
    write_top(data, NP_name, top_file)
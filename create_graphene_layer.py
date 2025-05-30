import numpy as np
import sys

# INPUT PARAMETERS:

lenX = float(sys.argv[1]) #length of x in Angstroms
lenY = float(sys.argv[2]) #length of y in Angstroms

nx = int(lenX / 5.1)+1
ny = int(lenY / 2.94)+1

a  = 1.7   # interatomic distance
#nx = 2    # number of repetitions in the x direction
#ny = 2     # number of repetitions in the y direction
filename = 'graphene_'+str(nx)+"_"+str(ny)+'.dat'

filename = 'graphene_'+str(lenX)+"_A_"+str(lenY)+"_A_"+str(nx)+"_"+str(ny)+'.dat'

writeresults = True
#showresults = true


# Size of the unit cell
A = 3*a
B = np.sqrt(3)*a

#Coordinates of the 4 atoms in the unit cell
base = [[0.0 , 0.0 , 0.0 ],[a/2 , B/2 , 0.0],[A/2 , B/2 , 0.0],[2*a , 0.0 , 0.0]]
#bonds = [[0,1],[1,2],[2,3]]
# Total number of atoms
#N = len(list(base))*(nx+ny+1)

# Calculate the coordinates of the atoms in the layer
coords = []
bonds = []
id = 0

#atomsDelY = 0
#atomsDelX = 0

for i in range(nx):
    coords_layer = []
    for j in range(ny):
        step = [i * A, j * B, 0]
        #print(step)
        for iatom in range(len(base)):
            id += 1
            coord_tmp = []
            for k in range(3):
                coord_tmp.append(base[iatom][k] + step[k])
            #if coord_tmp[0] > lenX:
            #    atomsDelX += 1
            #if coord_tmp[1] > lenY:
            #    atomsDelY += 1
            #if coord_tmp[0] <= lenX and coord_tmp[1] <= lenY: # For exact size of graphene layer (no replication further)
            coords.append(coord_tmp)
            coords_layer.append(coord_tmp)
                #print(base[iatom][k])
    for p, atom_c in enumerate(coords_layer):
        step_fw = i*ny*4
        id = p+1+step_fw
        if id%4 == 1:
            bonds.append([id,id+1])
        elif id%4 == 2:
            bonds.append([id,id+1])
            if ny*4+step_fw - id > 4:
                bonds.append([id,id+3])
            else:
                bonds.append([id,ny*4+step_fw-p-2])
        elif id%4 == 3:
            bonds.append([id,id+1])
            if ny * 4 + step_fw - id > 4:
                bonds.append([id,id+5])
            else:
                bonds.append([id, ny * 4 + step_fw+2 - p])
        elif id%4 == 0 and i != nx-1:
            bonds.append([id,id+ny*4-3])
        elif id%4 == 0 and i == nx-1:
            bonds.append([id,id-(nx-1)*ny*4-3])

#print(coords)
#print(len(coords), N)

#dnx = 0
#if atomsDelX/ny%3 == 0:
#    dnx = 0.5*a
#elif atomsDelX/ny%3 == 2:
#    dnx = 1.5*a
#elif atomsDelX/ny%3 == 1:
#    dnx = 2*a

#dny = 0
#if atomsDelY/nx%4 == 2:
#    dny = B/2

if writeresults:
    wFile = open(filename,'w')
    wFile.write('''graphene a=%d\n\n%d atoms\n1 atom types\n%d bonds\n1 bond types\n\n0 %f xlo xhi\n0 %f ylo yhi\n-10 10 zlo zhi\n\nMasses\n\n1 12.0107\n\nAtoms\n\n''' %
                (a, len(coords), len(bonds), A*nx, B*ny))
    #wFile.write('''graphene a=%d\n\n%d atoms\n1 atom types\n\n0 %f xlo xhi\n0 %f ylo yhi\n-10 10 zlo zhi\n\nMasses\n\n1 12.0107\n\nAtoms\n\n''' %
                #(a, len(coords), A * nx, B * ny))
    for i in range(len(coords)):
        #if coords[i][0] <= lenX and coords[i][1] <= lenY:
        wFile.write('%d 1 %f %f %f\n' % (i+1, coords[i][0], coords[i][1], coords[i][2]))
        #print(i)
    wFile.write('''\nBonds\n\n''')
    for k in range(len(bonds)):
        wFile.write('%d 1 %d %d\n' % (k+1, bonds[k][0], bonds[k][1]))
    wFile.close()
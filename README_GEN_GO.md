Example of a script to generate Graphene Oxide surface.

python create_graphene_layer.py 50 50 ### create graphene surface data-file with bonds and size of Lx=50 A Ly=50 A ---- graphene_50_A_50_A*.dat
python create_double_C.py graphene_50_A_50_A* double_bonds.data #### create file double_C_bonds.dat with atoms having double bond. The double_bonds.data is file showing double bonds (need to rewrite number of bonds to visualize)
python create_many_GOs_O_OH.py graphene_50_A_50_A* GO_layer_easy_50_A_50_A arr_file 1 300 140 #### Creating series of GO_layers with unique random positions of -OH and -O groups
##arr_file - name of file containing unique positions, 1 - number of files, 300 - number of -OH groups, 140 - number of -O groups
python create_H_COOH_atoms_v1.py GO_layer_easy_50_A_50_A.1.data GO_layer_medium_50_A_50_A.1.data # creating -H groups on surface and -COOH groups on the edges (if no bonds there) to the "lonely" C atoms
##also setting the parameters of the OPLS-AA forcefield
python create_go_angles_dihedrals_new_types_COOH_v1.py GO_layer_medium_50_A_50_A.1.data GO_layer_FULL_50_A_50_A.1.data # create angles and dihedrals, make some cosmetic changes with types

##GO_layer_FULL_50_A_50_A.1.data is the final LAMMPS data-file to use

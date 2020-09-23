# get_xtb_info README
Extract useful results from xTB calculations

This python script finds results from xTB output files (geometry optimizations and frequency calculations).
The results are printed out in a line in csv format.

Syntaxis:

python get_xtb_info.py [path_to_xTB_output_file] 

optional arguments: 
* [-X number_of_atom_X_in_input_geometry] 
* [-H number_of_atom_H_in_input_geometry] 
* [-xyz path_to_xtbopt.xyz_file]
* [--header]

A valid path to an output file from xTB calculations is required.
* Extract additional info regarding a specific atom X by adding -X and then the number of the atom in the input geometry, so NOT its atomic number.
* Extract additional info regarding an atom bound to atom X by adding -H and then the number of atom H in the input geometry.
* A path to the xyz file written by the program with the optimized geometry can be supplied to increase precision of X-H distance (X and H need to be added).
* Adding the --header option prints the header for the csv line instead of the line, which shows what info is stored in what position.

The printed line starts with a comma, so you have to manually add the index number before each line when looping multiple files.

When only a path to an output file is given, the following properties are extracted:
*HOMO,LUMO,Total_energy,Gibbs_energy

When atom number X is added, this additional info is extracted regarding atom X:
*species(C or N for example),atom_nr,charge_in_geometry,coordination_nr,sum_of_atomnrs_of_coordinating_atoms

When atom H is specified as well, it also extracts:
*charge_of_atom_H_in_geometry,X-H_distance

The line is printed as follows (add --header option for latest structure):
*,HOMO,LUMO,Etotal,Egibbs,Xspecies,Xatomnr,Xcharge,Xcoordination,Xsumcoord,Hcharge,XH
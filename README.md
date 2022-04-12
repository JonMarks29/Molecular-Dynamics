# Molecular Dynamics

gro_handler.py: A class to read in Gro files and handle individual protein/lipid residues and atoms. Facilitates file reading/writing, finding and removing solvent clashing atoms, removing/adding new lipids at specific locations, or replacing lipids at specific locations.
Basic uage:

g = Gro(pathtofile, exclude=[])   # exclude specific atoms
g.residue_count() # count number of protein residues

g.create_resiude(coordinates, new_residue) # lipid residues
g.remove_residue(residue_number)
g.replace_residue(coordinates, original_residue, new_residue, number=1, flip=False, dm=0, solvent_clash=True):

g.tidy_numbers() # check that atom numbering is correct and clean after file editing
g.write(outputname) # write out modified system to new file 

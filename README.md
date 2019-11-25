# PDB-XYZ Conversion

This directory contains my Python scripts for converting from PDB to XYZ and 
from XYZ to PDB. 
These scripts are a workaround for TINKER's `pdbxyz` and `xyzpdb`, which can 
fail with nucleic acids and solvated systems with numerous WAT/HOH residues.

## Notes:
`pdbxyz-for-amber.py` works for AMBER parameter sets.

`xyzpdb-for-amber.py` works for both AMBER and AMOEBA parameter sets.
- Residue writing is limited to 9999 total residues, before starting back at 1.
This doesn't seem to be an issue for visualization with VMD, but may be 
important for other programs.

## How
The first few lines specify the `infile`, `outfile`, and `prmfile`. Change these to reflect your system.

These scripts work by taking the atom lines from an input parameter file and 
reading them into a pandas dataframe.
The residue name from the name string (ex: "Isoleucine CB") is converted to 
its 3-letter code (or some general 3-letter code), and that field is split into
multiple columns (ResName and AtomName). 

**Important: Name strings cannot have more than 1 space, as this will throw
an error in parsing the parameter file.**

For best results, give any non-standard residues unique atom names, and use
the wanted 3-letter code in their name string.

In `pdbxyz`, the residue names and atom names are matched from this table, and 
the atom.mass is set to the TINKER atom type.
In `xyzpdb`, the TINKER atom type is used to set the residue and atom names.

`xyzpdb` for AMOEBA has the added challenge of protein residue types defaulting 
to alanine for the backbone (excluding proline and glycine). 
To address this, the insertion\_code attribute is set to the residue number.
Later on, the insertion\_code is removed so that the final PDB is numbered 
correctly.

The reverse of this issue is that some non-standard residues fall back on known
atom types. In that case, then you'll likey want to add some `if` statements
after the comment that says `Catch non-standard problems here!` in `pdbxyz`.

## Incomplete Conversion
If `pdbxyz` cannot find an atom type, the residue/atom it had issues with will
be printed to the Terminal. Additionally, the final XYZ will have 
`ATOM TYPE NOT FOUND` following the connectivity.

If `xyzpdb` cannot find a residue name, the residue name will be assigned 
as `FIX`. A warning that this occurred will be printed to the terminal.

If you're using `xyzpdb` for a visualization program, check that each segment has the requisite `TER` cards.

## Future plans:
- Adapt pdbxyz-for-amber.py to run with AMOEBA paramter sets as well.
- Add argparse / something for command-line input.
- Investigate limits to residue numbering/ways to get around that.
- Fix TER cards after protein in xyzpdb.
- Maybe actually formulate this into a Python package with tests?

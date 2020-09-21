# PDB-XYZ Conversion

This directory contains my Python scripts for converting from PDB to XYZ and
from XYZ to PDB.
These scripts are a workaround for TINKER's `pdbxyz` and `xyzpdb`, which can
fail with nucleic acids and solvated systems with numerous WAT/HOH residues.

## Notes:
`pdbxyz-for-amber.py` works for AMBER parameter sets.

`pdbxyz-for-amoeba.py` works for AMOEBA bio, pro, and nuc parameter sets.
- This assumes a generic `HIS` is `HIE` (protonated at N-epsilon). 
If you only use `HIS` in your PDB, make sure you modify different protonation
states to `HIP` (both N-epsilon and N-delta) or `HID` (N-delta) naming.

`xyzpdb-for-amber.py` works for both AMBER and AMOEBA parameter sets.
- Residue writing is limited to 9999 total residues, before starting back at 1.
This doesn't seem to be an issue for visualization with VMD, but may be
important for other programs.

`pdbxyz4amber-pmd-params.py` is meant to work with prmtop-based parameter files
created using
[`generate_TINKER_parameters.py`](https://github.com/emleddin/research-scripts/blob/main/tinker-params/generate_TINKER_parameters.py).
This is because the `.prm` files the script generates are written in a
standardized format, where atom name strings are composed of only the residue's
3-letter code and atom name.
The bulk of the original Python scripts were processing the parameter
files distributed with the TINKER program.
`xyzpdb4amber-pmd-params.py` will do the conversion back to a PDB using the
same prmtop-based parameter files.

### A Note on AMBER Parameters
- The TINKER distribution of `amber96.prm`, `amber99.prm`, `amber99sb.prm` uses
`Glutamic Acid` to refer to the AMBER `GLU` NOT `GLH`. As such, `GLH` is not
included in the parameter sets, and will need to beadded as a "non-standard."
`amber94.prm` and `amber98.prm` both use `Glutamic Acid (COOH)` for `GLH`.

## How
The first few lines specify the `infile`, `outfile`, and `prmfile`.
Change these to reflect your system.

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

If you're using `xyzpdb` for a visualization program, check that each segment
has the requisite `TER` cards.

## Errors Encountered with Analyze
Sometimes when using AMBER parameter files, TINKER `analyze` will say there are
missing angles or torsions.
These parameters are likely to be defined in the
[relevant AMBER parm files](http://http://ambermd.org/AmberModels.php).
Add them to your TINKER parameter file and you should be golden.

## Future plans:
- Consolidate for-amber and for-amoeba programs.
- Add argparse / something for command-line input.
- Investigate limits to residue numbering/ways to get around that.
- Fix TER cards after protein in xyzpdb.
- Maybe actually formulate this into a Python package with tests?
- Fix DX5 HO5' atoms classed incorrectly with AMBER prm files (1244 instead of
1245).

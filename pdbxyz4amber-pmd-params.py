import parmed as pmd
import numpy as np
import pandas as pd
import re
import sys

infile="my_amber_sim.pdb"
outfile="my_amber_sim_convert.xyz"

param_file="TINKER_my_amber_sim.prm"
atom_lines="atom-lines.txt"

test_csv="test.csv"

def clean_atoms(temp):
    ## Deal with C and N Terminals
    ## If there's a terminal OXT for the protein, use as CTERM
    for residue in temp.residues:
        for atom in residue.atoms:
            if atom.name == 'OXT':
                residue.name = 'C'+residue.name
            ## Standardize ions
            ## If the +/- is last character, don't escape it.
            elif atom.name in ('K+', 'K', 'K\+1', 'k', 'k+'):
                atom.name = 'K'
                residue.name = 'K'
            elif atom.name in ('MG', 'MG2', 'MG2+', 'Mg2', 'MG\+2', 'Mg2+',
             'Mg\+2'):
                atom.name = 'MG'
                residue.name = 'MG'
            elif atom.name in ('ZN', 'ZN2', 'ZN2+', 'Zn2', 'ZN\+2', 'Zn2+',
             'Zn\+2'):
                atom.name = 'ZN'
                residue.name = 'ZN'
            elif atom.name in ('CL', 'CL-', 'CL1-', 'Cl', 'CL\-1', 'Cl1-',
             'Cl\-1'):
                atom.name = 'CL'
                residue.name = 'CL'
            else:
                continue
        ## If there's an H3 in a protein residue, use NTERM
        if residue.name in ('ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'CYX', 'GLN',
         'GLU', 'GLY', 'HID', 'HIE', 'HIP', 'HIS', 'ILE', 'LEU', 'LYS', 'MET',
         'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'):
            for atom in residue.atoms:
                if atom.name == 'H3':
                    residue.name = 'N'+residue.name


## Mercileslly taken from Mark's pdbtinker.py
def load_pdb(filename):
    """Loads in PDB using parmed and sets atom masses to zero. Atom masses are
    then used to store the Tinker atom types for XYZ conversion."""
    system = pmd.load_file(filename)
    for atom in system.atoms:
        atom.mass = 0
    clean_atoms(system)
    return system

## grep "^atom" param_file > atom_lines
def read_prm(param_file, atom_lines):
    """Based on using grep to locate atom lines."""
    parm_atom_lines = []
    pattern = '(?i)^atom ' # case insensitive, starts line
                           # keep the space to not match "atomic" in AMOEBA
    for line in open(param_file).readlines():
        if re.match(pattern, line):
            parm_atom_lines.append(line)
    ## Write a file of all the atom lines for use later
    with open(atom_lines, 'w') as filehandle:
        filehandle.writelines("%s" % line for line in parm_atom_lines)
    filehandle.close()

def fix_params(atom_lines,test_csv):
    """Read the atom_lines file into a pandas object. Rewrite the story
    (quoted) section into PDB residue names and atom names through as series
    of string replacements. This will work well for AMBER or AMOEBA sets."""
    lines = pd.read_csv(atom_lines, sep='[\s]{2,}', header=None,
     names=["what","T_type","T_atom_class","A_atom_type","A_names","element",
     "mass","connectivity"], engine='python')
    ##
    ## Remove the double quotes in the A_names column
    lines.A_names = lines.A_names.replace({'"':''}, regex=True)
    ##
    ## Determine if AMOEBA or AMBER based on aspartic acid
    ## While they vary by the space, this *should* check user-adjusted lines
    ## Labeled Aspartate in AMOEBA sets, but aspartic acid in AMBER
    ## AMBER version needs to be ASP, AMOEBA is ASH
    trial = lines[lines['A_names'].str.contains(r'(?i)Aspartate', regex=True)]
    if trial.A_names.empty == False:
        print("Processing as AMOEBA parameters.")
        AMOEBA = True
        # lines.A_names = lines.A_names.str.replace(r'(?i)XXX', 'YYY',
        #   regex=True)
    else:
        print("Processing as AMBER parameters.")
        AMOEBA = False
    ##
    ## Replace names of special residues (use the ?i regex to ignore case)
    #
    # lines.A_names = lines.A_names.str.replace(r'(?i)XXX', 'YYY',
    #   regex=True)
    ##
    ## The distributied TINKER params list incorrectly uses Glutamic Acid to
    ## mean GLU NOT glutamate.
    ## If you need the actual glutamic acid, you need to add params for GLH
    ##
    ## Now standardize the ions for string search
    lines.A_names = lines.A_names.str.replace(r'(?i)Li\+', 'LI', regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Na\+', 'NA', regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)K\+', 'K', regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Rb\+', 'RB', regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Cs\+', 'CS', regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Be\+', 'BE', regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Mg\+2', 'MG', regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Ca\+2', 'CA', regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Zn\+2', 'ZN', regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)F\-', 'F', regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Cl\-', 'CL', regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Br\-', 'BR', regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)I\-', 'I', regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Ba\+2', 'BA', regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Sr\+2', 'SR', regex=True)
    ##
    ## User-defined
    ## The search string is what is listed in the string of the prm file's
    ## atom line for a given atom type and the replace string is the ResName
    ## in the PDB
    ##
    ## From LEAPRC files
    lines.A_names = lines.A_names.str.replace(r'(?i)TP3', 'WAT',
      regex=True)
    ##
    ## Print the new lines for testing
    #lines.to_csv(test_csv, index=False, encoding='utf8')
    ## Now split the column into two columns
    convert = lines['A_names'].str.split(" ", n=1, expand=True)
    ## Add those new columns back into the atom lines
    lines['ResName'] = convert[0]
    lines['AtomName'] = convert[1]
    ## Print the new lines for testing
    lines.to_csv(test_csv, index=False, encoding='utf8')
    return lines, AMOEBA

def convert_names(system, lines, AMOEBA):
    """For every atom, find lines matching the residue name in the
    pandas_object. From those lines, check for lines that match the atom name.
    If a match isn't found, check through the known naming problems. Update the
    atom mass with the TINKER type if a match is found; if no match is found,
    keep the atom mass as zero.
    """
    print(
    '''If I didn't find residues, they'll be listed here:
    Residue Name | Atom Name | Search ResName | Search Atom Name
    ''')
    for residue in system.residues:
        for atom in residue.atoms:
            test_name = atom.name
            res_test = lines[lines.ResName == residue.name]
            atom_test = res_test[res_test.AtomName == test_name]
            ## Address problem residues!
            # if atom_test.empty == True:
            #     if residue.name in ('AAA'):
            #         test_RN = 'XXX'
            #         test_name = 'J'
            #         res_test = lines[lines.ResName == test_RN]
            #         atom_test = res_test[res_test.AtomName == test_name]
            #     elif residue.name in ('BBB'):
            #         test_RN = 'YYY'
            #         test_name = 'K'
            #         res_test = lines[lines.ResName == test_RN]
            #         atom_test = res_test[res_test.AtomName == test_name]
            #     ### And so on and so forth
            if atom_test.empty == True:
                ## Prints out what you need to fix :)
                try:
                    print(residue.name, atom.name, test_RN, test_name)
                    atom.mass = int(0)
                except UnboundLocalError:
                    print(residue.name, atom.name, residue.name, test_name)
                atom.mass = int(0)
            else:
                try:
                    atom.mass += int(atom_test['T_type'])
                except TypeError:
                    print("""
                    Oof, please check the parameter file.
                       I have too many matches...
                    """)
                    print(residue.name, atom.name, "ResID:", residue.number)
                # atom.mass += atom_test['T_type']
    return system

## Taken from Mark's PDBTinker
def write_xyz(system, outfile):
    ## Add XYZ extension
    if outfile.split(".")[-1] != "xyz":
        outfile = outfile+".xyz"
    ## Write to outfile
    f = open(outfile,"w+")
    f.write(str(len(system.atoms))+"\n")
    for atom in system.atoms:
        bondstring = ""
        bondlist = []
        for i in atom.bonds:
            if i.atom1.idx == atom.idx:
                bondlist.append(i.atom2.idx+1)
                # bondstring = bondstring + str(i.atom2.idx+1) + "\t"
            elif i.atom2.idx == atom.idx:
                bondlist.append(i.atom1.idx+1)
                # bondstring = bondstring + str(i.atom1.idx+1) + "\t"
        if atom.name == "SS" and atom.residue.name == "CYX" and len(bondlist) \
         < 2:
            for atom2 in system.atoms:
                if atom2.residue.name == "CYX" and atom2.name == "SS":
                    ia_dist = interatomic_distance(atom,atom2)
                    if ia_dist < 2 and ia_dist > 0:
                        bondlist.append(atom2.idx)
        bondlist.sort()
        bondstring = ''.join(str('{:>8}'.format(x)) for x in bondlist)
        index = atom.idx+1
        name = atom.name
        x = atom.xx
        y = atom.xy
        z = atom.xz
        atomtype = atom.mass
        linestring = "{:>6}  {:<4} {:>12.6f} {:>11.6f} {:>11.6f} {:>5}{}".\
         format(str(index),str(name),x,y,z,str(atomtype),bondstring)
        if atomtype == 0:
            linestring = linestring + " ATOM TYPE NOT FOUND"
        linestring = linestring + "\n"
        f.write(linestring)
    f.close()

## Begin function calls
read_prm(param_file, atom_lines)

try:
    lines, AMOEBA = fix_params(atom_lines,test_csv)
except pd.errors.ParserError:
    print("""
          ,-~~-.___.
     / |  '     \              Something's fishy with your parameter file.
    (  )         0             In order for me to do my job, I need the string
     \_/-, ,----'              portion (what's in quotes) of the lines starting
        ====           //      with "atom" to have no more than 1 space between
       /  \-'~;    /~~~(O)     words.
      /  __/~|   /       |
    =(  _____| (_________|     Once you fix that, it *should* be smooth sailing.
    No promises, though. I am a Python script after all.
    """)
    sys.exit()

system = load_pdb(infile)

system = convert_names(system, lines, AMOEBA)

write_xyz(system, outfile)

# ## Write out the XYZ
# system.save(outfile, overwrite=True)

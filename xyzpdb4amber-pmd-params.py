import parmed as pmd
import numpy as np
import pandas as pd
import re
import sys

infile="my_amber_sim_convert.xyz"
outfile="my_amber_sim.pdb"

param_file="TINKER_my_amber_sim.prm"

## These are written out to help ID any problems
## atom_lines is what it gleans from the parameter file
## test_csv is how the program matches atom and residues
atom_lines="atom-lines.txt"
test_csv="test.csv"

## Note: This script assumes that your XYZ has each new AA residue starting with
## a nitrogen and each new nucleic acid residue starting with a phosphate or
## terminal hydrogen.
## Any untyped residues will have the residue name "FIX"

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

## Mercileslly taken from Mark's pdbtinker.py
def load_xyz(filename):
    """Loads in XYZ using parmed and sets atom masses to zero. Atom masses are
    then used to store the Tinker atom types for PDB conversion."""
    system = pmd.load_file(filename)
    for atom in system.atoms:
        atom.mass = 0
    return system

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

def clean_atom_names(system):
    """Clean up the names of some ions for easier conversion."""
    for atom in system.atoms:
        ## CA for alpha carbon and Calcium can be tricky bastards
        if atom.name in ('Ca', 'Ca+', 'CA2'):
            atom.name = 'CA2'
        elif atom.name in ('MG2', 'Mg2', 'MG'):
            atom.name = 'MG'
        elif atom.name in ('Cl', 'CL-', 'CL'):
            atom.name = 'CL'
        elif atom.name in ('K+', 'K'):
            atom.name = 'K'
        elif atom.name in ('ZN2', 'ZN2+', "Zn", "Zn2"):
            atom.name = 'ZN'
    return system

def convert_names(system, lines, AMOEBA):
    """For every atom, find lines matching the atom type in the
    pandas_object. From those lines, get the residue name.
    If a match isn't found, check through the known naming problems. Update the
    atom mass with the TINKER type if a match is found; if no match is found,
    keep the atom mass as zero.
    """
    ## build an empty structure
    new_pdb = pmd.structure.Structure()
    resnumber = 0
    for atom in system.atoms:
        for type in atom.type:
            test_name = int(atom.type)
            type_test = lines[lines.T_type == test_name]
            try:
                my_res = type_test['ResName'].values[0]
                #if AMOEBA == True:
                #    ## Address misguided sidechain atom names
                #    if atom.name in ('N'):
                #        atom.name = type_test['AtomName'].values[0]
                #else:
                if my_res not in ('CA', 'MG', 'K', 'ZN',
                 'CL', 'NTE', 'CTE'):
                    atom.name = type_test['AtomName'].values[0]
            except IndexError:
                print("""
                Oops, I couldn't figure out the residue name.
                Check your output for the FIX resiude(s).""")
                my_res = 'FIX'
        ## Deal with HEM Non-standard (Because it has N's that'll get flagged)
        if atom.name in ('FE') and my_res in ('HEM'):
            resnumber += 1
        ## Figure out the residue counting based on N and HO5' residues
        # if atom.name in ('N', 'HO5\'', 'P', 'MG', 'K', 'ZN', 'CL'):
        elif atom.name in ('N', 'HO5\'', 'P', 'MG', 'K', 'ZN', 'CL', 'CA2'):
            resnumber += 1
        elif my_res == 'WAT' and atom.name == 'O':
            resnumber += 1
        #############################################################
        ##  May need to add in a clause here for any non-standards ##
        #############################################################
        ## Give it a 3-letter
        if my_res in ('NALA', 'NARG', 'NASN', 'NASP', 'NCYS', 'NCYX', 'NGLN',
         'NGLU', 'NGLY', 'NHID', 'NHIE', 'NHIP', 'NHIS', 'NILE', 'NLEU', 'NLYS',
         'NMET', 'NPHE', 'NPRO', 'NSER', 'NTHR', 'NTRP', 'NTYR', 'NVAL', 'CALA',
         'CARG', 'CASN', 'CASP', 'CCYS', 'CCYX', 'CGLN', 'CGLU', 'CGLY', 'CHID',
         'CHIE', 'CHIP', 'CHIS', 'CILE', 'CLEU', 'CLYS', 'CMET', 'CPHE', 'CPRO',
         'CSER', 'CTHR', 'CTRP', 'CTYR', 'CVAL'):
         my_res = my_res[1:]
        # print(atom.name, my_res, resnumber)
        new_pdb.add_atom(atom=atom, resname=my_res, resnum=int(resnumber),
         inscode=str(resnumber))
    return system, new_pdb

def amoeba_name_fix(new_pdb, AMOEBA):
    """AMOEBA uses ALA as a base residue for everything except GLY and PRO.
    Save the residue numbers as a dictionary, then reassign the residue numbers
    to the residues."""
    if AMOEBA == True:
        my_res_dict = {}
        for residue in new_pdb.residues:
            ## Ensure certain mid-residue definitions aren't skipped
            if residue.insertion_code in my_res_dict:
                if my_res_dict[str(residue.insertion_code)] == 'CYX':
                    continue
                elif my_res_dict[str(residue.insertion_code)] in ('DC3', 'DG3',
                'DA3', 'DT3', 'RC3', 'RG3', 'RA3', 'RU3'):
                    residue.ter = True
                    continue
                elif residue.name == 'CTE':
                    residue.ter = True
                    continue
                ## Might need non-standard here
                elif residue.name == 'CTP':
                    continue
                else:
                    my_res_dict.update({residue.insertion_code : residue.name})
            else:
                my_res_dict.update({residue.insertion_code : residue.name})
        # print(my_res_dict)

        ## Forcefully reassign residue names
        for residue in new_pdb.residues:
            for atom in residue.atoms:
                residue.name = my_res_dict[str(residue.insertion_code)]

        amoeba_pdb = pmd.structure.Structure()
        resnumber = 0
        for residue in new_pdb.residues:
            for atom in residue.atoms:
                if residue.name in ('HEM') and atom.name in ('FE'):
                    resnumber += 1
                ## Figure out the residue counting based on N and HO5' residues
                # if atom.name in ('N', 'HO5\'', 'P', 'MG', 'K', 'ZN', 'CL'):
                elif atom.name in ('N', 'HO5\'', 'P', 'MG', 'K', 'ZN', 'CL',
                 'CA2'):
                    resnumber += 1
                elif residue.name == 'WAT' and atom.name == 'O':
                    resnumber += 1
                amoeba_pdb.add_atom(atom=atom, resname=residue.name,
                 resnum=(resnumber))
        new_pdb = amoeba_pdb
    else:
        ## Remove the placeholder insertion codes
        for residue in new_pdb.residues:
            residue.insertion_code = ''
    return new_pdb

## Mercileslly taken from Mark's pdbtinker.py
def write_pdb(new_pdb, filename):
    """Writes out the converted PDB."""
    ## Python may start at zero, but resnumbers don't!
    counter = 0
    for residue in new_pdb.residues:
        if residue.name in ('DX5', 'RX5'):
            residue.name = new_pdb.residues[counter+1].name + "5"
        elif residue.name in ('DX3', 'RX3') and residue.atoms[0].name in \
         ('P', 'OP1', 'OP2'):
            residue.name = new_pdb.residues[counter+1].name + "3"
        elif residue.name in ('DX3', 'RX3') and residue.atoms[0].name not in \
         ('P', 'OP1', 'OP2'):
            ## shouldn't have 3, because the previous ones are renamed
            residue.name = new_pdb.residues[counter-1].name
            ## Add a TER card
            residue.ter = True
        elif residue.name in ('DX', 'RX'):
            residue.name = new_pdb.residues[counter+1].name
        elif residue.name in ('DC', 'DG', 'DA', 'DT') and residue.atoms[0].name not in \
         ('P'):
            if new_pdb.residues[counter-1].name in ('DC5', 'DG5', 'DA5', 'DT5'):
                residue.name = new_pdb.residues[counter-1].name
            elif new_pdb.residues[counter+1].name and \
             new_pdb.residues[counter-1].name in ('DC3', 'DG3', 'DA3', 'DT3'):
                residue.name = new_pdb.residues[counter-1].name
        elif residue.name in ('RC', 'RG', 'RA', 'RU') and residue.atoms[0].name not in \
         ('P'):
            if new_pdb.residues[counter-1].name in ('RC5', 'RG5', 'RA5', 'RU5'):
                residue.name = new_pdb.residues[counter-1].name
            elif new_pdb.residues[counter+1].name and \
             new_pdb.residues[counter-1].name in ('RC3', 'RG3', 'RA3', 'RU3'):
                residue.name = new_pdb.residues[counter-1].name
        elif residue.name in ('MG', 'K', 'ZN', 'CL', 'CA',
         'DC3', 'DG3', 'DA3', 'DT3', 'RC3', 'RG3', 'RA3', 'RU3'):
            residue.ter = True

        if AMOEBA == False:
            for atom in residue.atoms:
                if atom.name == 'OXT':
                    residue.ter = True

        counter += 1
    ## Renumber the atoms so they're not all -1
    a_count = 1
    for atom in new_pdb.atoms:
        atom.number = a_count
        a_count += 1
    # new_pdb.write_pdb(filename)
    ## Use our explicitly labeled residue numbers; don't renumber everything
    new_pdb.write_pdb(filename, renumber=False)
    print("Please check the TER lines before using in a visualization program.")
    return new_pdb

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

system = load_xyz(infile)

system = clean_atom_names(system)

system, new_pdb = convert_names(system, lines, AMOEBA)

new_pdb = amoeba_name_fix(new_pdb, AMOEBA)

write_pdb(new_pdb, outfile)

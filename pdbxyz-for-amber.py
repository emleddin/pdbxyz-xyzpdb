import parmed as pmd
import numpy as np
import pandas as pd
import re
import sys

infile="crambin.pdb"
outfile="crambin_convert.xyz"

param_file="amoebapro13.prm"
atom_lines="atom-lines.txt"

test_csv="test.csv"

def clean_atoms(temp):
    """ILE, LEU, ASN, THR, GLN, and VAL all have screwy 4-letter hydrogens,
    as do ALL DNA residues. This fixes that."""
    for atom in temp.atoms:
        if atom.name == 'HD11':
            atom.name = '1HD1'
        if atom.name == 'HD12':
            atom.name = '2HD1'
        if atom.name == 'HD13':
            atom.name = '3HD1'
        if atom.name == 'HD21':
            atom.name = '1HD2'
        if atom.name == 'HD22':
            atom.name = '2HD2'
        if atom.name == 'HD23':
            atom.name = '3HD2'
        if atom.name == 'HG11':
            atom.name = '1HG1'
        if atom.name == 'HG12':
            atom.name = '2HG1'
        if atom.name == 'HG13':
            atom.name = '3HG1'
        if atom.name == 'HG21':
            atom.name = '1HG2'
        if atom.name == 'HG22':
            atom.name = '2HG2'
        if atom.name == 'HG23':
            atom.name = '3HG2'
        if atom.name == 'HE21':
            atom.name = '1HE2'
        if atom.name == 'HE22':
            atom.name = '2HE2'
        # if atom.name == 'H2\'1':
        #     atom.name = '1H2\''
        # if atom.name == 'H2\'2':
        #     atom.name = '2H2\''
    ## Terminal hydrogens on DNA are different
    ## in the stinking parameter file, so fix it.
    for atom in temp.atoms:
        if atom.name == 'H5\'':
            atom.name = 'H5\'1'
            # atom.name = '1H5\''
        if atom.name == 'H5\'\'':
            atom.name = 'H5\'2'
            # atom.name = '2H5\''
        if atom.name == 'H2\'':
            atom.name = 'H2\'1'
            # atom.name = '1H2\''
        if atom.name == 'H2\'\'':
            atom.name = 'H2\'2'
            # atom.name = '2H2\''
        if atom.name == 'HO3\'':
            atom.name = 'H3T'
        if atom.name == 'HO2\'':
            atom.name = 'HO\'2'
        # if atom.name == 'H5\'1':
        #     atom.name = '1H5\''
        # if atom.name == 'H5\'2':
        #     atom.name = '2H5\''
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

# class Typing:
#     def __init__(self, T_type, T_atom_class, A_atom_type, res_name, atom_name):
#         self.T_type = T_type
#         self.T_atom_class = T_atom_class
#         self.A_atom_type = A_atom_type
#         self.res_name = res_name
#         self.atom_name = atom_name
#
# a1 = Typing(1, 14, "N", "Glycine", "N")
# print(a1.T_type, a1.T_atom_class, a1.A_atom_type, a1.res_name, a1.atom_name)

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
        lines.A_names = lines.A_names.str.replace(r'(?i)Aspartic Acid', 'ASH',
          regex=True)
        lines.A_names = lines.A_names.str.replace(r'(?i)Aspartate', 'ASP',
         regex=True)
        lines.A_names = lines.A_names.str.replace(r'(?i)Glutamate', 'GLU',
         regex=True)
        lines.A_names = lines.A_names.str.replace(r'(?i)N\-MeAmide Cap', 'XXX',
         regex=True)
        lines.A_names = lines.A_names.str.replace(r'(?i)Amide Cap', 'XXX',
         regex=True)
        lines.A_names = lines.A_names.str.replace(r'(?i)N\-Terminal PRO',
         'NPRO', regex=True)
        lines.A_names = lines.A_names.str.replace(r'(?i)C\-Terminal COOH',
         'XXX', regex=True)
        lines.A_names = lines.A_names.str.replace(r'(?i)N\-Terminal',
         'NTE', regex=True)
        lines.A_names = lines.A_names.str.replace(r'(?i)C\-Terminal',
         'CTE', regex=True)
    else:
        print("Processing as AMBER parameters.")
        AMOEBA = False
    ##
    ## Replace names of special residues (use the ?i regex to ignore case)
    ## Histidine
    lines.A_names = lines.A_names.str.replace(r'(?i)Histidine \(HD\)', 'HID',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Histidine \(HE\)', 'HIE',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Histidine \(\+\)', 'HIP',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)HIS \(HD\)', 'HID',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)HIS \(HE\)', 'HIE',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)HIS \(\+\)', 'HIP',
      regex=True)
    ## Cysteine
    lines.A_names = lines.A_names.str.replace(r'(?i)Cysteine \(\-SH\)', 'CYS',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Cystine \(\-SS\-\)', 'CYX',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Cysteine Anion', 'CYM',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Cystine', 'CYX',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)CYS \(\-SH\)', 'CYS',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)CYS \(\-SS\-\)', 'CYX',
      regex=True)
    ## Lysine
    lines.A_names = lines.A_names.str.replace(r'(?i)Lysine \(NH2\)', 'LYN',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Lysine \(Neutral\)', 'LYN',
      regex=True)
    ## Tyrosine
    lines.A_names = lines.A_names.str.replace(r'(?i)Tyrosine Anion', 'TYD',
      regex=True)
    ## Asp
    lines.A_names = lines.A_names.str.replace(r'(?i)Aspartic Acid \(COOH\)',
      'ASH', regex=True)
    ## Caps and Termini
    lines.A_names = lines.A_names.str.replace(r'(?i)Acetyl Cap', 'ACE',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)N\-MeAmide Cap', 'NME',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)N\-MeAmide', 'NME',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)N\-Term ', 'N',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)C\-Term ', 'C',
      regex=True)
    ## Atypical residues
    lines.A_names = lines.A_names.str.replace(r'(?i)Ornithine', 'ORN',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)MethylAlanine', 'AIB',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Pyroglutamate', 'PCA',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Acetyl', 'ACE',
      regex=True)
    ##
    ## Replace standards with 3 letter codes
    lines.A_names = lines.A_names.str.replace(r'(?i)Aspartic Acid', 'ASP',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Glutamic Acid', 'GLH',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Phenylalanine', 'PHE',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Alanine', 'ALA',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Arginine', 'ARG',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Asparagine', 'ASN',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Cysteine', 'CYS',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Glutamine', 'GLN',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Glycine', 'GLY',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Histidine', 'HIS',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Isoleucine', 'ILE',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Leucine', 'LEU',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Lysine', 'LYS',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Methionine', 'MET',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Proline', 'PRO',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Serine', 'SER',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Threonine', 'THR',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Tryptophan', 'TRP',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Tyrosine', 'TYR',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Valine', 'VAL',
      regex=True)
    ##
    ## RNA
    lines.A_names = lines.A_names.str.replace(r'(?i)R-Adenosine', 'RA',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)R-Guanosine', 'RG',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)R-Cytosine', 'RC',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)R-Uracil', 'RU',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)R\-Phosphodiester', 'RX',
      regex=True)
    ## The O5' and HO5' have same name in prm file; address it
    lines.A_names = lines.A_names.str.replace(r'(?i)R\-5\'\-Hydroxyl O5\'',
      'RX5 HO5\'', regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)R\-5\'\-Hydroxyl', 'RX5',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)R\-5\'\-Phosphate', 'RX5',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)R\-3\'\-Phosphate', 'RX3',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)R\-3\'\-Hydroxyl O3\'',
      'RX3 HO3\'', regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)R\-3\'\-Hydroxyl', 'RX3',
      regex=True)
    ##
    ## DNA
    lines.A_names = lines.A_names.str.replace(r'(?i)D-Adenosine', 'DA',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)D-Guanosine', 'DG',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)D-Cytosine', 'DC',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)D-Thymine', 'DT',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)D-Phosphodiester', 'DX',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)D\-5\'\-Hydroxyl O5\'',
      'DX5 HO5\'', regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)D\-5\'\-Hydroxyl', 'DX5',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)D\-5\'\-Phosphate', 'DX5',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)D\-3\'\-Hydroxyl O3\'',
      'DX3 HO3\'', regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)D\-3\'\-Hydroxyl', 'DX3',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)D\-3\'\-Phosphate', 'DX3',
      regex=True)
    ##
    ## AMOEBA nucleics...
    lines.A_names = lines.A_names.replace({'\(CT\)':''}, regex=True)
    lines.A_names = lines.A_names.replace({'\(CU\)':''}, regex=True)
    lines.A_names = lines.A_names.replace({'\(AG\)':''}, regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Deoxyribose', 'DX',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Ribose', 'RX',
      regex=True)
    ##
    ## Do the same with water
    lines.A_names = lines.A_names.str.replace(r'(?i)TIP3P Oxygen', 'WAT O',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)TIP3P Hydrogen', 'WAT H',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)AMOEBA Water', 'WAT',
      regex=True)
    ##
    ## Remove the word "Ion"
    lines.A_names = lines.A_names.replace({' Ion':''}, regex=True)
    ## Do the same with ions
    lines.A_names = lines.A_names.str.replace(r'(?i)Lithium', 'LI',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Sodium', 'NA', regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Potassium', 'K',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Rubidium', 'RB', regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Cesium', 'CS', regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Beryllium', 'BE',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Magnesium', 'MG',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Calcium', 'CA',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Zinc', 'ZN', regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Fluoride', 'F', regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Chloride', 'CL',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Bromide', 'BR', regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Iodide', 'I', regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Barium', 'BA', regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Strontium', 'SR',
      regex=True)
    ## Now standardize them for string search
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
    # lines.A_names = lines.A_names.str.replace(r'(?i)DUP -Phosphate', 'DUP',
    #   regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)DUP-Uracil', 'DUP',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)DGP and DCP PO4', 'DXP',
      regex=True)
    # lines.A_names = lines.A_names.str.replace(r'(?i)DUP-Uracil', 'CTP',
    #   regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)DUP -Phosphate', 'CTP',
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
            ## Address DNA, RNA, and CTERM/NTERM
            if atom_test.empty == True:
                if residue.name in ('NALA', 'NARG', 'NASN', 'NASP', 'NCYS',
                 'NCYX', 'NGLN', 'NGLU', 'NGLY', 'NHID', 'NHIE', 'NHIP',
                 'NHIS', 'NILE', 'NLEU', 'NLYS', 'NMET', 'NPHE', 'NPRO',
                 'NSER', 'NTHR', 'NTRP', 'NTYR', 'NVAL'):
                    if test_name in ('H1', 'H2', 'H3'):
                        # It's HN for N-Terminal
                        test_name = 'HN'
                        res_test = lines[lines.ResName == residue.name]
                        atom_test = res_test[res_test.AtomName == test_name]
                    else:
                        test_RN = residue.name[1:]
                        test_name = atom.name
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                        ##### Fix the duplicate/triplicates here
                        if atom_test.empty == True:
                            if test_name == ('H'):
                                test_name = 'HN'
                                atom_test = res_test[res_test.AtomName == test_name]
                            elif test_name in ('HA2', 'HA3'):
                                test_name = 'HA'
                                atom_test = res_test[res_test.AtomName == test_name]
                            elif test_name in ('HB1', 'HB2', 'HB3'):
                                test_name = 'HB'
                                atom_test = res_test[res_test.AtomName == test_name]
                            elif test_name in ('HG1', 'HG2', 'HG3'):
                                test_name = 'HG'
                                atom_test = res_test[res_test.AtomName == test_name]
                            elif test_name in ('1HG1', '2HG1', '3HG1'):
                                test_name = 'HG1'
                                atom_test = res_test[res_test.AtomName == test_name]
                            elif test_name in ('1HG2', '2HG2', '3HG2'):
                                test_name = 'HG2'
                                atom_test = res_test[res_test.AtomName == test_name]
                            elif test_name in ('HD1', 'HD2', 'HD3'):
                                test_name = 'HD'
                                atom_test = res_test[res_test.AtomName == test_name]
                            elif test_name in ('1HD1', '2HD1', '3HD1'):
                                test_name = 'HD1'
                                atom_test = res_test[res_test.AtomName == test_name]
                                # HD condition
                                if atom_test.empty == True:
                                    test_name = 'HD'
                                    atom_test = res_test[res_test.AtomName == test_name]
                            elif test_name in ('1HD2', '2HD2', '3HD2'):
                                test_name = 'HD2'
                                atom_test = res_test[res_test.AtomName == test_name]
                            elif test_name in ('HE1', 'HE2', 'HE3'):
                                test_name = 'HE'
                                atom_test = res_test[res_test.AtomName == test_name]
                            elif test_name in ('1HE1', '2HE1', '3HE1'):
                                test_name = 'HE1'
                                atom_test = res_test[res_test.AtomName == test_name]
                            elif test_name in ('1HE2', '2HE2', '3HE2'):
                                test_name = 'HE2'
                                atom_test = res_test[res_test.AtomName == test_name]
                            elif test_name in ('HH11', 'HH12'):
                                test_name = 'HH1'
                                atom_test = res_test[res_test.AtomName == test_name]
                                if atom_test.empty == True:
                                    test_name = 'HH'
                                    atom_test = res_test[res_test.AtomName == test_name]
                            elif test_name in ('HH21', 'HH22'):
                                test_name = 'HH2'
                                atom_test = res_test[res_test.AtomName == test_name]
                                if atom_test.empty == True:
                                    test_name = 'HH'
                                    atom_test = res_test[res_test.AtomName == test_name]
                            elif test_name in ('HZ1', 'HZ2', 'HZ3'):
                                test_name = 'HZ'
                                atom_test = res_test[res_test.AtomName == test_name]
                            elif test_name in ('NH1', 'NH2'):
                                test_name = 'NH'
                                atom_test = res_test[res_test.AtomName == test_name]
                            elif test_name in ('OD1', 'OD2'):
                                test_name = 'OD'
                                atom_test = res_test[res_test.AtomName == test_name]
                            elif test_name in ('OE1', 'OE2'):
                                test_name = 'OE'
                                atom_test = res_test[res_test.AtomName == test_name]
                            # CD condition
                            elif test_name in ('CD1', 'CD2'):
                                test_name = 'CD'
                                atom_test = res_test[res_test.AtomName == test_name]
                            # CE vs CE1 condition
                            elif test_name in ('CE1', 'CE2'):
                                test_name = 'CE'
                                atom_test = res_test[res_test.AtomName == test_name]
                            elif test_name in ('H71', 'H72', 'H73'):
                                test_name = 'H7'
                                atom_test = res_test[res_test.AtomName == test_name]
                            else:
                                atom_test = res_test[res_test.AtomName == test_name]
                elif residue.name in ('CALA', 'CARG', 'CASN', 'CASP', 'CCYS',
                 'CCYX', 'CGLN', 'CGLU', 'CGLY', 'CHID', 'CHIE', 'CHIP',
                 'CHIS', 'CILE', 'CLEU', 'CLYS', 'CMET', 'CPHE', 'CPRO',
                 'CSER', 'CTHR', 'CTRP', 'CTYR', 'CVAL'):
                    test_RN = residue.name[1:]
                    test_name = atom.name
                    res_test = lines[lines.ResName == test_RN]
                    atom_test = res_test[res_test.AtomName == test_name]
                    ### Fix the duplicates/triplicates here
                    if atom_test.empty == True:
                        if AMOEBA == True:
                            test_RN = 'CTE'
                            res_test = lines[lines.ResName == test_RN]
                            if atom.name in ('C'):
                                atom_test = res_test[res_test.A_atom_type == atom.name]
                            elif atom.name in ('OXT'):
                                test_name = 'O'
                                atom_test = res_test[res_test.A_atom_type == test_name]
                            elif atom.name in ('N','CA','O'):
                                test_RN = 'ALA'
                                res_test = lines[lines.ResName == test_RN]
                                atom_test = res_test[res_test.AtomName == atom.name]
                        elif test_name == ('H'):
                            test_name = 'HN'
                            atom_test = res_test[res_test.AtomName == test_name]
                        elif test_name in ('HA2', 'HA3'):
                            test_name = 'HA'
                            atom_test = res_test[res_test.AtomName == test_name]
                        elif test_name in ('HB1', 'HB2', 'HB3'):
                            test_name = 'HB'
                            atom_test = res_test[res_test.AtomName == test_name]
                        elif test_name in ('HG1', 'HG2', 'HG3'):
                            test_name = 'HG'
                            atom_test = res_test[res_test.AtomName == test_name]
                        elif test_name in ('1HG1', '2HG1', '3HG1'):
                            test_name = 'HG1'
                            atom_test = res_test[res_test.AtomName == test_name]
                        elif test_name in ('1HG2', '2HG2', '3HG2'):
                            test_name = 'HG2'
                            atom_test = res_test[res_test.AtomName == test_name]
                        elif test_name in ('HD1', 'HD2', 'HD3'):
                            test_name = 'HD'
                            atom_test = res_test[res_test.AtomName == test_name]
                        elif test_name in ('1HD1', '2HD1', '3HD1'):
                            test_name = 'HD1'
                            atom_test = res_test[res_test.AtomName == test_name]
                            # HD condition
                            if atom_test.empty == True:
                                test_name = 'HD'
                                atom_test = res_test[res_test.AtomName == test_name]
                        elif test_name in ('1HD2', '2HD2', '3HD2'):
                            test_name = 'HD2'
                            atom_test = res_test[res_test.AtomName == test_name]
                        elif test_name in ('HE1', 'HE2', 'HE3'):
                            test_name = 'HE'
                            atom_test = res_test[res_test.AtomName == test_name]
                        elif test_name in ('1HE1', '2HE1', '3HE1'):
                            test_name = 'HE1'
                            atom_test = res_test[res_test.AtomName == test_name]
                        elif test_name in ('1HE2', '2HE2', '3HE2'):
                            test_name = 'HE2'
                            atom_test = res_test[res_test.AtomName == test_name]
                        elif test_name in ('HH11', 'HH12'):
                            test_name = 'HH1'
                            atom_test = res_test[res_test.AtomName == test_name]
                            if atom_test.empty == True:
                                test_name = 'HH'
                                atom_test = res_test[res_test.AtomName == test_name]
                        elif test_name in ('HH21', 'HH22'):
                            test_name = 'HH2'
                            atom_test = res_test[res_test.AtomName == test_name]
                            if atom_test.empty == True:
                                test_name = 'HH'
                                atom_test = res_test[res_test.AtomName == test_name]
                        elif test_name in ('HZ1', 'HZ2', 'HZ3'):
                            test_name = 'HZ'
                            atom_test = res_test[res_test.AtomName == test_name]
                        elif test_name in ('NH1', 'NH2'):
                            test_name = 'NH'
                            atom_test = res_test[res_test.AtomName == test_name]
                        elif test_name in ('OD1', 'OD2'):
                            test_name = 'OD'
                            atom_test = res_test[res_test.AtomName == test_name]
                        elif test_name in ('OE1', 'OE2'):
                            test_name = 'OE'
                            atom_test = res_test[res_test.AtomName == test_name]
                        # CD condition
                        elif test_name in ('CD1', 'CD2'):
                            test_name = 'CD'
                            atom_test = res_test[res_test.AtomName == test_name]
                        # CE vs CE1 condition
                        elif test_name in ('CE1', 'CE2'):
                            test_name = 'CE'
                            atom_test = res_test[res_test.AtomName == test_name]
                        elif test_name in ('H71', 'H72', 'H73'):
                            test_name = 'H7'
                            atom_test = res_test[res_test.AtomName == test_name]
                        else:
                            atom_test = res_test[res_test.AtomName == test_name]
                    #####
                elif residue.name in ('DC5', 'DG5', 'DA5', 'DT5'):
                    if atom.name in ('O5\'', 'HO5\'', 'H5T', 'P'):
                        test_RN = 'DX5'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == atom.name]
                    elif atom.name in ('OP1', 'OP2', 'O1P', 'O2P'):
                        test_name = 'OP'
                        atom_test = res_test[res_test.AtomName == test_name]
                    else:
                        test_RN = residue.name[:-1]
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == atom.name]
                        if atom.name in ('H71', 'H72', 'H73'):
                            test_name = 'H7'
                            atom_test = res_test[res_test.AtomName == test_name]
                elif residue.name in ('DC3', 'DG3', 'DA3', 'DT3'):
                    test_RN = 'DX3'
                    res_test = lines[lines.ResName == test_RN]
                    if atom.name in ('O3\'', 'HO3\'', 'H3T', 'P'):
                        atom_test = res_test[res_test.AtomName == atom.name]
                    elif atom.name in ('OP1', 'OP2', 'O1P', 'O2P'):
                        test_name = 'OP'
                        atom_test = res_test[res_test.AtomName == test_name]
                    else:
                        test_RN = residue.name[:-1]
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == atom.name]
                        if atom.name in ('H71', 'H72', 'H73'):
                            test_name = 'H7'
                            atom_test = res_test[res_test.AtomName == test_name]
                elif residue.name in ('DC', 'DG', 'DA', 'DT'):
                    if atom.name in ('H71', 'H72', 'H73'):
                        test_name = 'H7'
                        atom_test = res_test[res_test.AtomName == test_name]
                    else:
                        test_RN = 'DX'
                        res_test = lines[lines.ResName == test_RN]
                        if atom.name in ('P'):
                            atom_test = res_test[res_test.AtomName == atom.name]
                        elif atom.name in ('OP1', 'OP2', 'O1P', 'O2P'):
                            test_name = 'OP'
                            atom_test = res_test[res_test.AtomName == test_name]
                elif residue.name in ('C5', 'G5', 'A5', 'U5'):
                    test_RN = 'RX5'
                    res_test = lines[lines.ResName == test_RN]
                    if atom.name in ('O5\'', 'HO5\'', 'H5T', 'P'):
                        atom_test = res_test[res_test.AtomName == atom.name]
                    elif atom.name in ('OP1', 'OP2', 'O1P', 'O2P'):
                        test_name = 'OP'
                        atom_test = res_test[res_test.AtomName == test_name]
                    else:
                        test_RN = 'R' + residue.name[:-1]
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == atom.name]
                elif residue.name in ('C3', 'G3', 'A3', 'U3'):
                    test_RN = 'RX3'
                    res_test = lines[lines.ResName == test_RN]
                    if atom.name in ('O3\'', 'HO3\'', 'H3T', 'P'):
                        atom_test = res_test[res_test.AtomName == atom.name]
                    elif atom.name in ('OP1', 'OP2', 'O1P', 'O2P'):
                        test_name = 'OP'
                        atom_test = res_test[res_test.AtomName == test_name]
                    else:
                        test_RN = 'R' + residue.name[:-1]
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == atom.name]
                elif residue.name in ('C', 'G', 'A', 'U'):
                    test_RN = 'RX'
                    res_test = lines[lines.ResName == test_RN]
                    if atom.name in ('P'):
                        atom_test = res_test[res_test.AtomName == atom.name]
                    elif atom.name in ('OP1', 'OP2', 'O1P', 'O2P'):
                        test_name = 'OP'
                        atom_test = res_test[res_test.AtomName == test_name]
                    else:
                        test_RN = 'R' + residue.name
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                elif residue.name in ('RC5', 'RG5', 'RA5', 'RU5'):
                    test_RN = 'RX5'
                    res_test = lines[lines.ResName == test_RN]
                    if atom.name in ('O5\'', 'HO5\'', 'H5T', 'P'):
                        atom_test = res_test[res_test.AtomName == atom.name]
                    elif atom.name in ('OP1', 'OP2', 'O1P', 'O2P'):
                        test_name = 'OP'
                        atom_test = res_test[res_test.AtomName == test_name]
                    else:
                        test_RN = residue.name[:-1]
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == atom.name]
                elif residue.name in ('RC3', 'RG3', 'RA3', 'RU3'):
                    test_RN = 'RX3'
                    res_test = lines[lines.ResName == test_RN]
                    if atom.name in ('O3\'', 'HO3\'', 'H3T', 'P'):
                        atom_test = res_test[res_test.AtomName == atom.name]
                    elif atom.name in ('OP1', 'OP2', 'O1P', 'O2P'):
                        test_name = 'OP'
                        atom_test = res_test[res_test.AtomName == test_name]
                    else:
                        test_RN = residue.name[:-1]
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == atom.name]
                elif residue.name in ('RC', 'RG', 'RA', 'RU'):
                    test_RN = 'RX'
                    res_test = lines[lines.ResName == test_RN]
                    if atom.name in ('P'):
                        atom_test = res_test[res_test.AtomName == atom.name]
                    elif atom.name in ('OP1', 'OP2', 'O1P', 'O2P'):
                        test_name = 'OP'
                        atom_test = res_test[res_test.AtomName == test_name]
                elif residue.name in ('ACE'):
                    if atom.name in ('HH31', 'HH32', 'HH33'):
                        test_name = 'HA'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('CH3'):
                        test_name = 'CA'
                        atom_test = res_test[res_test.AtomName == test_name]
                elif residue.name in ('NME'):
                    if atom.name in ('HH31', 'HH32', 'HH33'):
                        test_name = 'HC'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('CH3'):
                        test_name = 'C'
                        atom_test = res_test[res_test.AtomName == test_name]
                #############################################
                ##    Deal with AMOEBA ALA Standard...     ##
                #############################################
                elif AMOEBA == True:
                    # IGNORE Proline, Glycine, Cysteine Anion (CYM), & N-TERM
                    if residue.name in ('ASN','ASH','CYS','CYX','GLH',\
                    'GLN','HID','HIE','HIP','HIS','LYS','MET',\
                    'SER','THR','TRP'):
                        if atom.name in ('N','CA','C','O'):
                            test_RN = 'ALA'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
                    elif residue.name in ('ARG'):
                        if atom.name in ('N','CA','C','O'):
                            test_RN = 'ALA'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
                        elif atom.name in ('NH1', 'NH2'):
                            test_name = 'NH'
                            atom_test = res_test[res_test.AtomName == test_name]
                        ## Might need condition for HH
                    elif residue.name in ('ASP'):
                        if atom.name in ('N','CA','C','O'):
                            test_RN = 'ALA'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
                        elif atom.name in ('OD1', 'OD2'):
                            test_name = 'OD'
                            atom_test = res_test[res_test.AtomName == test_name]
                    elif residue.name in ('GLU'):
                        if atom.name in ('N','CA','C','O'):
                            test_RN = 'ALA'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
                        elif atom.name in ('OE1', 'OE2'):
                            test_name = 'OE'
                            atom_test = res_test[res_test.AtomName == test_name]
                    elif residue.name in ('ILE'):
                        if atom.name in ('N','CA','C','O'):
                            test_RN = 'ALA'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
                        elif atom.name in ('CD1', 'CD2'):
                            test_name = 'CD'
                            atom_test = res_test[res_test.AtomName == test_name]
                        ## Might need condition for HD
                    elif residue.name in ('LEU'):
                        if atom.name in ('N','CA','C','O'):
                            test_RN = 'ALA'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
                        elif atom.name in ('CD1', 'CD2'):
                            test_name = 'CD'
                            atom_test = res_test[res_test.AtomName == test_name]
                        ## Might need condition for HD
                    elif residue.name in ('PHE'):
                        if atom.name in ('N','CA','C','O'):
                            test_RN = 'ALA'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
                        elif atom.name in ('CE1', 'CE2'):
                            test_name = 'CE'
                            atom_test = res_test[res_test.AtomName == test_name]
                        ## Might need condition for HE
                        elif atom.name in ('CD1', 'CD2'):
                            test_name = 'CD'
                            atom_test = res_test[res_test.AtomName == test_name]
                        ## Might need condition for HD
                    elif residue.name in ('TYR'):
                        if atom.name in ('N','CA','C','O'):
                            test_RN = 'ALA'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
                        elif atom.name in ('CD1', 'CD2'):
                            test_name = 'CD'
                            atom_test = res_test[res_test.AtomName == test_name]
                        ## Might need condition for HD
                        elif atom.name in ('CE1', 'CE2'):
                            test_name = 'CE'
                            atom_test = res_test[res_test.AtomName == test_name]
                        ## Might need condition for HE
                    elif residue.name in ('VAL'):
                        if atom.name in ('N','CA','C','O'):
                            test_RN = 'ALA'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
                        elif atom.name in ('CG1', 'CG2'):
                            test_name = 'CG'
                            atom_test = res_test[res_test.AtomName == test_name]
                        ## Might need condition for HG
                #############################################
                ##    Catch non-standard problems here!    ##
                #############################################
                ## Reinterpret CTP
                elif residue.name == 'CTP':
                    test_RN = 'DC'
                    res_test = lines[lines.ResName == test_RN]
                    atom_test = res_test[res_test.AtomName == atom.name]
                    if atom_test.empty == True:
                        test_RN = 'DUP'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == atom.name]
                ## For the deoxy residues
                ## Deoxy 5mC == 5CM
                elif residue.name == '5CM':
                    test_RN = 'DC'
                    res_test = lines[lines.ResName == test_RN]
                    atom_test = res_test[res_test.AtomName == atom.name]
                    if atom_test.empty == True:
                        test_RN = 'DX'
                        res_test = lines[lines.ResName == test_RN]
                        if atom.name in ('P'):
                            atom_test = res_test[res_test.AtomName == atom.name]
                        elif atom.name in ('OP1', 'OP2', 'O1P', 'O2P'):
                            test_name = 'OP'
                            atom_test = res_test[res_test.AtomName == test_name]
                ## Deoxy 5hmC == 5HC
                elif residue.name == '5HC':
                    test_RN = 'DC'
                    res_test = lines[lines.ResName == test_RN]
                    atom_test = res_test[res_test.AtomName == atom.name]
                    if atom_test.empty == True:
                        test_RN = 'DX'
                        res_test = lines[lines.ResName == test_RN]
                        if atom.name in ('P'):
                            atom_test = res_test[res_test.AtomName == atom.name]
                        elif atom.name in ('OP1', 'OP2', 'O1P', 'O2P'):
                            test_name = 'OP'
                            atom_test = res_test[res_test.AtomName == test_name]
                ## Ribo residues
                ## Ribo 5hmC == 5hC
                elif residue.name == '5hC':
                    test_RN = 'RC'
                    res_test = lines[lines.ResName == test_RN]
                    atom_test = res_test[res_test.AtomName == atom.name]
                    if atom_test.empty == True:
                        test_RN = 'RX'
                        res_test = lines[lines.ResName == test_RN]
                        if atom.name in ('P'):
                            atom_test = res_test[res_test.AtomName == atom.name]
                        elif atom.name in ('OP1', 'OP2', 'O1P', 'O2P'):
                            test_name = 'OP'
                            atom_test = res_test[res_test.AtomName == test_name]
                ## Ribo 5mC == 5mC
                elif residue.name == '5mC':
                    test_RN = 'RC'
                    res_test = lines[lines.ResName == test_RN]
                    atom_test = res_test[res_test.AtomName == atom.name]
                    if atom_test.empty == True:
                        test_RN = 'RX'
                        res_test = lines[lines.ResName == test_RN]
                        if atom.name in ('P'):
                            atom_test = res_test[res_test.AtomName == atom.name]
                        elif atom.name in ('OP1', 'OP2', 'O1P', 'O2P'):
                            test_name = 'OP'
                            atom_test = res_test[res_test.AtomName == test_name]
                ########################################################
                ##    Now we're just a typical protein residue lol    ##
                ########################################################
                else:
                    if test_name == ('H'):
                        test_name = 'HN'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif test_name in ('HA2', 'HA3'):
                        test_name = 'HA'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif test_name in ('H1', 'H2', 'H3'):
                        test_name = 'HA'
                        atom_test = res_test[res_test.AtomName == test_name]
                        if atom_test.empty == True:
                            test_name = 'H'
                            atom_test = res_test[res_test.AtomName == test_name]
                    elif test_name in ('HB1', 'HB2', 'HB3'):
                        test_name = 'HB'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif test_name in ('HG1', 'HG2', 'HG3'):
                        test_name = 'HG'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif test_name in ('1HG1', '2HG1', '3HG1'):
                        test_name = 'HG1'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif test_name in ('1HG2', '2HG2', '3HG2'):
                        test_name = 'HG2'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif test_name in ('HD1', 'HD2', 'HD3'):
                        test_name = 'HD'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif test_name in ('1HD1', '2HD1', '3HD1'):
                        test_name = 'HD1'
                        atom_test = res_test[res_test.AtomName == test_name]
                        # HD condition
                        if atom_test.empty == True:
                            test_name = 'HD'
                            atom_test = res_test[res_test.AtomName == test_name]
                    elif test_name in ('1HD2', '2HD2', '3HD2'):
                        test_name = 'HD2'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif test_name in ('HE1', 'HE2', 'HE3'):
                        test_name = 'HE'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif test_name in ('1HE1', '2HE1', '3HE1'):
                        test_name = 'HE1'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif test_name in ('1HE2', '2HE2', '3HE2'):
                        test_name = 'HE2'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif test_name in ('HH11', 'HH12'):
                        test_name = 'HH1'
                        atom_test = res_test[res_test.AtomName == test_name]
                        if atom_test.empty == True:
                            test_name = 'HH'
                            atom_test = res_test[res_test.AtomName == test_name]
                    elif test_name in ('HH21', 'HH22'):
                        test_name = 'HH2'
                        atom_test = res_test[res_test.AtomName == test_name]
                        if atom_test.empty == True:
                            test_name = 'HH'
                            atom_test = res_test[res_test.AtomName == test_name]
                    elif test_name in ('HZ1', 'HZ2', 'HZ3'):
                        test_name = 'HZ'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif test_name in ('NH1', 'NH2'):
                        test_name = 'NH'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif test_name in ('OD1', 'OD2'):
                        test_name = 'OD'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif test_name in ('OE1', 'OE2'):
                        test_name = 'OE'
                        atom_test = res_test[res_test.AtomName == test_name]
                    # CD condition
                    elif test_name in ('CD1', 'CD2'):
                        test_name = 'CD'
                        atom_test = res_test[res_test.AtomName == test_name]
                    # CE vs CE1 condition
                    elif test_name in ('CE1', 'CE2'):
                        test_name = 'CE'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif test_name in ('H71', 'H72', 'H73'):
                        test_name = 'H7'
                        atom_test = res_test[res_test.AtomName == test_name]
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

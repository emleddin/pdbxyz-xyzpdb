## This is the version meant for AMOEBA parameter sets!

## !! NOTE !!
## This version will only work if the nucleic acid strands are explicitly
## labeled for 5' and 3' ends (e.g., using A3, A5, RA3, RA5, DA3, DA5)

## In the future, I would love for the convert_names call to have
## subfunctions for N Term/C Term, DNA/RNA, NSA
## And even have an AMOBA call too (so it's just one function with if statements)
## Probably want to consider N/C term as common with the except of proline & N/H

import parmed as pmd
import numpy as np
import pandas as pd
import re
import sys

infile="PDB_test_NC.pdb"
outfile="PDB_test_NC.xyz"

# param_file="amoebapro13.prm"
# param_file="amoebabio09.prm"
param_file="amoebabio18.prm"
atom_lines="atom-lines.txt"

test_csv="test.csv"

def clean_atoms_AMOEBA(temp):
    """AMOEBA doesn't use 4-letter hydrogens because everything is awful."""
    for atom in temp.atoms:
        if atom.name in ('HD11', 'HD12', 'HD13'):
            atom.name = 'HD1'
        if atom.name in ('HD21', 'HD22', 'HD23'):
            atom.name = 'HD2'
        if atom.name in ('HG11', 'HG12', 'HG13'):
            atom.name = 'HG1'
        if atom.name in ('HG21', 'HG22', 'HG23'):
            atom.name = 'HG2'
        if atom.name in ('HE21', 'HE22'):
            atom.name = 'HE2'
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
        ## AMOEBA18 doesn't distinguish this
        if atom.name == 'H5\'\'':
            atom.name = 'H5\'2'
            # atom.name = '2H5\''
        if atom.name == 'H2\'':
            atom.name = 'H2\'1'
            # atom.name = '1H2\''
        ## AMOEBA19 doesn't distinguish this
        if atom.name == 'H2\'\'':
            atom.name = 'H2\'2'
            # atom.name = '2H2\''
        if atom.name == 'HO3\'':
            atom.name = 'H3T'
        if atom.name == 'HO5\'':
            atom.name = 'H5T'
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
            ## Standardize ions
            ## If the +/- is last character, don't escape it.
            if atom.name in ('K+', 'K', 'K\+1', 'k', 'k+'):
                atom.name = 'K'
                residue.name = 'K'
            elif atom.name in ('NA', 'NA+', 'NA\+1', 'NA1+', 'Na', 'Na+',
             'Na1+'):
                atom.name = 'NA'
                residue.name = 'NA'
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
        ## Assume HIE as default
        if residue.name in ('HIS'):
            print("Assuming HIE as the default protonation for HIS", residue.number)
            print(" If you get an error for HIE HD1, please verify if that")
            print(" residue should be either HID or HIP and modify your input.\n")
            residue.name = 'HIE'
        ## This just makes the deoxy/ribo distinction cleaner & addresses the
        ## "new"/AMBER naming scheme
        elif residue.name in ('A'):
            residue.name = 'RA'
        elif residue.name in ('A3'):
            residue.name = 'RA3'
        elif residue.name in ('A5'):
            residue.name = 'RA5'
        elif residue.name in ('C'):
            residue.name = 'RC'
        elif residue.name in ('C3'):
            residue.name = 'RC3'
        elif residue.name in ('C5'):
            residue.name = 'RC5'
        elif residue.name in ('G'):
            residue.name = 'RG'
        elif residue.name in ('G3'):
            residue.name = 'RG3'
        elif residue.name in ('G5'):
            residue.name = 'RG5'
        elif residue.name in ('U'):
            residue.name = 'RU'
        elif residue.name in ('U3'):
            residue.name = 'RU3'
        elif residue.name in ('U5'):
            residue.name = 'RU5'
        ## If there's an H3 in a protein residue, use NTERM
        if residue.name in ('ALA', 'ARG', 'ASH', 'ASN', 'ASP', 'CYM', 'CYS',
         'CYX', 'GLH', 'GLN', 'GLU', 'GLY', 'HID', 'HIE', 'HIP', 'ILE',
         'LEU', 'LYN', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYD',
         'TYR', 'VAL'):
            for atom in residue.atoms:
                if atom.name == 'H3':
                    residue.name = 'N'+residue.name
                ## And CTERM for OXT
                elif atom.name == 'OXT':
                    residue.name = 'C'+residue.name

## Mercileslly taken from Mark's pdbtinker.py
def load_pdb(filename, AMOEBA):
    """Loads in PDB using parmed and sets atom masses to zero. Atom masses are
    then used to store the Tinker atom types for XYZ conversion."""
    system = pmd.load_file(filename)
    for atom in system.atoms:
        atom.mass = 0
    if AMOEBA == True:
        clean_atoms_AMOEBA(system)
    else:
        # clean_atoms(system)
        print("Try another version of this script, please.\n\
        This one hasn't passed all the AMBER-based tests.")
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
        print("Processing as AMOEBA parameters.\n")
        AMOEBA = True
        lines.A_names = lines.A_names.str.replace(r'(?i)Aspartic Acid', 'ASH',
          regex=True)
        lines.A_names = lines.A_names.str.replace(r'(?i)Aspartate', 'ASP',
         regex=True)
        lines.A_names = lines.A_names.str.replace(r'(?i)Glutamate', 'GLU',
         regex=True)
        lines.A_names = lines.A_names.str.replace(r'(?i)Glutamic Acid', 'GLH',
          regex=True)
        lines.A_names = lines.A_names.str.replace(r'(?i)N\-MeAmide Cap', 'XXX',
         regex=True)
        lines.A_names = lines.A_names.str.replace(r'(?i)Amide Cap', 'XXX',
         regex=True)
        lines.A_names = lines.A_names.str.replace(r'(?i)N\-Terminal PRO',
         'NPRO', regex=True)
        lines.A_names = lines.A_names.str.replace(r'(?i)C\-Terminal COOH',
         'XXX', regex=True)
        lines.A_names = lines.A_names.str.replace(r'(?i)N\-Terminal NH3\+',
         'NTE N', regex=True)
        lines.A_names = lines.A_names.str.replace(r'(?i)N\-Terminal H3N\+',
         'NTE HN', regex=True)
        lines.A_names = lines.A_names.str.replace(r'(?i)C\-Terminal COO-',
         'CTEO', regex=True)
        lines.A_names = lines.A_names.str.replace(r'(?i)C\-Terminal',
         'CTE', regex=True)
        ### FIGURE OUT DNA/RNA
        ## Start with Ns that change
        lines.A_names = lines.A_names.str.replace(r'(?i)Adenine N9 RNA',
         'RA N9', regex=True)
        lines.A_names = lines.A_names.str.replace(r'(?i)Adenine N9 DNA',
         'DA N9', regex=True)
        lines.A_names = lines.A_names.str.replace(r'(?i)Cytosine N1 RNA',
         'RC N1', regex=True)
        lines.A_names = lines.A_names.str.replace(r'(?i)Cytosine N1 DNA',
         'DC N1', regex=True)
        lines.A_names = lines.A_names.str.replace(r'(?i)Guanine N9 RNA',
         'RG N9', regex=True)
        lines.A_names = lines.A_names.str.replace(r'(?i)Guanine N9 DNA',
         'DG N9', regex=True)
        ## General
        lines.A_names = lines.A_names.str.replace(r'(?i)Adenine',
         'DRA', regex=True)
        lines.A_names = lines.A_names.str.replace(r'(?i)Cytosine',
         'DRC', regex=True)
        lines.A_names = lines.A_names.str.replace(r'(?i)Guanine',
         'DRG', regex=True)
        lines.A_names = lines.A_names.str.replace(r'(?i)Thymine',
         'DT', regex=True)
        lines.A_names = lines.A_names.str.replace(r'(?i)Uracil',
         'RU', regex=True)
        ## Specific to amoebabio18
        lines.A_names = lines.A_names.str.replace(r'(?i)R\-5\'\-Hydroxyl O5\'T',
         'RX5 O5T', regex=True)
        lines.A_names = lines.A_names.str.replace(r'(?i)R\-3\'\-Hydroxyl O3\'T',
         'RX3 O3T', regex=True)
        lines.A_names = lines.A_names.str.replace(r'(?i)D\-5\'\-Hydroxyl O5\'T',
         'DX5 O5T', regex=True)
        lines.A_names = lines.A_names.str.replace(r'(?i)D\-3\'\-Hydroxyl O3\'T',
         'DX3 O3T', regex=True)
        lines.A_names = lines.A_names.str.replace(r'(?i)R-Phosphodiester',
         'RP', regex=True)
        lines.A_names = lines.A_names.str.replace(r'(?i)D-Phosphodiester',
         'DP', regex=True)
        lines.A_names = lines.A_names.str.replace(r'(?i)Phosphodiester',
         'DRP', regex=True)
        #
    else:
        print("Processing as AMBER parameters.")
        AMOEBA = False
        ## The default params list incorrectly uses Glutamic Acid to mean GLU
        ## NOT glutamate. If you need the actual glutamic acid, you need to add
        ## params for GLH
        lines.A_names = lines.A_names.str.replace(r'(?i)Glutamic Acid', 'GLU',
          regex=True)
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
    # lines.A_names = lines.A_names.str.replace(r'(?i)R\-Phosphodiester', 'RX',
    #   regex=True)
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
    # lines.A_names = lines.A_names.str.replace(r'(?i)D-Phosphodiester', 'DX',
    #   regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)D\-5\'\-Hydroxyl O5\'T',
      'DX5 O5\'', regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)D\-5\'\-Hydroxyl', 'DX5',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)D\-5\'\-Phosphate', 'DX5',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)D\-3\'\-Hydroxyl O3\'T',
      'DX3 O3\'', regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)D\-3\'\-Hydroxyl', 'DX3',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)D\-3\'\-Phosphate', 'DX3',
      regex=True)
    ## Change how HO'5 and HO'3 appear
    lines.A_names = lines.A_names.str.replace(r'(?i)HO\'5', 'HO5\'',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)HO\'3', 'HO3\'',
      regex=True)

    ##
    ## AMOEBA nucleics...
    ##
    ## Deal with Amoebabio09 pyrimadine and purines
    ## Based on https://stackoverflow.com/questions/45308569/python-regex-match-and-replace-beginning-and-end-of-string-but-keep-the-middle
    lines.A_names = lines.A_names.replace(to_replace='^(.*? )?(.*?)( \(CU\))$', value=r'RCU \2', regex=True)
    lines.A_names = lines.A_names.replace(to_replace='^R(.*? )?(.*?)( \(AG\))$', value=r'RAG \2', regex=True)
    lines.A_names = lines.A_names.replace(to_replace='^(.*? )?(.*?)( \(CT\))$', value=r'DCT \2', regex=True)
    lines.A_names = lines.A_names.replace(to_replace='^D(.*? )?(.*?)( \(AG\))$', value=r'DAG \2', regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Deoxyribose', 'DX',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)Ribose', 'RX',
      regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)HO\'5', 'H5T', regex=True)
    lines.A_names = lines.A_names.str.replace(r'(?i)HO\'3', 'H3T', regex=True)
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

def convert_names_AMOEBA(system, lines):
    """For every atom, find lines matching the residue name in the
    pandas_object. From those lines, check for lines that match the atom name.
    If a match isn't found, check through the known naming problems. Update the
    atom mass with the TINKER type if a match is found; if no match is found,
    keep the atom mass as zero. This is intended for AMOEBA parameter sets,
    because they default to ALA.
    """
    print(
    '''If I didn't find residues, they'll be listed here:
    Residue Name | Atom Name | Search ResName | Search Atom Name
    ''')
    for residue in system.residues:
        for atom in residue.atoms:
            ## Set RES_NAME so that you don't mess with weird copies and can
            ## Overwrite it for N/C TERM atoms in loop
            RES_NAME = residue.name
            test_name = atom.name
            #############################################
            ##    Deal with NTERM...     ##
            #############################################
            if residue.name in ('NALA', 'NARG', 'NASH', 'NASN', 'NASP', 'NCYM',
             'NCYS', 'NCYX', 'NGLH', 'NGLN', 'NGLU', 'NGLY', 'NHID', 'NHIE',
             'NHIP', 'NHIS', 'NILE', 'NLEU', 'NLYN', 'NLYS', 'NMET', 'NPHE',
             'NSER', 'NTHR', 'NTRP', 'NTYD', 'NTYR', 'NVAL'):
                # It's HN for N-Terminal
                if atom.name in ('H1', 'H2', 'H3'):
                    res_test = lines[lines.ResName == 'NTE']
                    atom_test = res_test[res_test.AtomName == 'HN']
                elif atom.name == 'N':
                    res_test = lines[lines.ResName == 'NTE']
                    atom_test = res_test[res_test.AtomName == 'N']
                else:
                    RES_NAME = residue.name[1:]
            elif residue.name == 'NPRO':
                if atom.name in ('N'):
                    res_test = lines[lines.ResName == 'NPRO']
                    atom_test = res_test[res_test.AtomName == 'NH2+']
                elif atom.name in ('H1', 'H2', 'H3'):
                    res_test = lines[lines.ResName == 'NPRO']
                    atom_test = res_test[res_test.AtomName == 'H2N+']
                elif atom.name in ('CA', 'C', 'O', 'HA', 'CD', 'HD'):
                    test_name = atom.name
                    res_test = lines[lines.ResName == residue.name]
                    atom_test = res_test[res_test.AtomName == test_name]
                else:
                    RES_NAME = residue.name[1:]
            elif residue.name in ('CALA', 'CARG', 'CASH', 'CASN', 'CASP',
             'CCYM', 'CCYS', 'CCYX', 'CGLH', 'CGLN', 'CGLU', 'CGLY', 'CHID',
             'CHIE', 'CHIP', 'CHIS', 'CILE', 'CLEU', 'CLYN', 'CLYS', 'CMET',
             'CPHE', 'CPRO', 'CSER', 'CTHR', 'CTRP', 'CTYD', 'CTYR', 'CVAL'):
                if atom.name in ('OXT','C'):
                    # It's HN for N-Terminal
                    test_RN = 'CTEO'
                    test_name = atom.name
                    res_test = lines[lines.ResName == test_RN]
                    ## Use the first letter for OXT vs C because of how param
                    ## file is written....
                    atom_test = res_test[res_test.A_atom_type == test_name[0]]
                else:
                    RES_NAME = residue.name[1:]
            if RES_NAME in ('ALA', 'ARG', 'ASH', 'ASN', 'ASP', 'CYS',\
             'CYM', 'CYX', 'GLH', 'GLN', 'GLU', 'GLY', 'HID', 'HIE', 'HIP',\
             'HIS', 'ILE', 'LEU', 'LYN', 'LYS', 'MET', 'PHE', 'PRO',\
             'SER', 'THR', 'TRP', 'TYD', 'TYR', 'VAL'):
                test_name = atom.name
                res_test = lines[lines.ResName == RES_NAME]
                atom_test = res_test[res_test.AtomName == test_name]
                #############################################
                ##    Deal with AMOEBA ALA Standard...     ##
                #############################################
                # IGNORE Proline, Glycine, Cysteine Anion (CYM), & N-TERM
                if RES_NAME in ('ALA'):
                    if atom.name in ('HB1', 'HB2', 'HB3'):
                        test_name = 'HB'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('H'):
                        test_name = 'HN'
                        atom_test = res_test[res_test.AtomName == test_name]
                elif RES_NAME in ('ARG'):
                    if atom.name in ('N','CA','C','O','HA'):
                        test_RN = 'ALA'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('NH1', 'NH2'):
                        test_name = 'NH'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('H', 'HN'):
                        test_name = 'HN'
                        test_RN = 'ALA'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('HB2', 'HB3'):
                        test_name = 'HB'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('HG2', 'HG3'):
                        test_name = 'HG'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('HD2', 'HD3'):
                        test_name = 'HD'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('HH11', 'HH12', 'HH21', 'HH22'):
                        test_name = 'HH'
                        atom_test = res_test[res_test.AtomName == test_name]
                elif RES_NAME in ('ASH','TRP'):
                    if atom.name in ('N','CA','C','O','HA'):
                        test_RN = 'ALA'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('H','HN'):
                        test_name = 'HN'
                        test_RN = 'ALA'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('HB2','HB3'):
                        test_name = 'HB'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    else:
                        atom_test = res_test[res_test.AtomName == test_name]
                elif RES_NAME in ('ASN'):
                    if atom.name in ('N','CA','C','O','HA'):
                        test_RN = 'ALA'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('HB2', 'HB3'):
                        test_name = 'HB'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('H', 'HN'):
                        test_RN = 'ALA'
                        test_name = 'HN'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                elif RES_NAME in ('ASP'):
                    if atom.name in ('N','CA','C','O','HA'):
                        test_RN = 'ALA'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('H', 'HN'):
                        test_RN = 'ALA'
                        test_name = 'HN'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('OD1', 'OD2'):
                        test_name = 'OD'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('HB2', 'HB3'):
                        test_name = 'HB'
                        atom_test = res_test[res_test.AtomName == test_name]
                elif RES_NAME in ('CYM'):
                    if atom.name in ('N','CA','C','O','HA'):
                        test_RN = 'ALA'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('H', 'HN'):
                        test_RN = 'ALA'
                        test_name = 'HN'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('HB2', 'HB3'):
                        test_name = 'HB'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('CB', 'HB'):
                        test_RN = 'CYS'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('SG'):
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                elif RES_NAME in ('CYS'):
                    if atom.name in ('N','CA','C','O','HA'):
                        test_RN = 'ALA'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('H', 'HN'):
                        test_RN = 'ALA'
                        test_name = 'HN'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('HB2', 'HB3'):
                        test_name = 'HB'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                elif RES_NAME in ('CYX'):
                    if atom.name in ('N','CA','C','O','HA'):
                        test_RN = 'ALA'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('H', 'HN'):
                        test_RN = 'ALA'
                        test_name = 'HN'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('HB2', 'HB3'):
                        test_name = 'HB'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('CB'):
                        test_RN = 'CYS'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                elif RES_NAME in ('GLH'):
                    if atom.name in ('N','CA','C','O','HA'):
                        test_RN = 'ALA'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('H','HN'):
                        test_name = 'HN'
                        test_RN = 'ALA'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('HB2', 'HB3'):
                        test_name = 'HB'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('HG2', 'HG3'):
                        test_name = 'HG'
                        atom_test = res_test[res_test.AtomName == test_name]
                elif RES_NAME in ('GLN'):
                    if atom.name in ('N','CA','C','O','HA'):
                        test_RN = 'ALA'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('H', 'HN'):
                        test_RN = 'ALA'
                        test_name = 'HN'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('HB2', 'HB3'):
                        test_name = 'HB'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('HG2', 'HG3'):
                        test_name = 'HG'
                        atom_test = res_test[res_test.AtomName == test_name]
                    else:
                        test_RN = RES_NAME
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                elif RES_NAME in ('GLU'):
                    if atom.name in ('N','CA','C','O','HA'):
                        test_RN = 'ALA'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('H', 'HN'):
                        test_RN = 'ALA'
                        test_name = 'HN'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('OE1', 'OE2'):
                        test_name = 'OE'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('HB2', 'HB3'):
                        test_name = 'HB'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('HG2', 'HG3'):
                        test_name = 'HG'
                        atom_test = res_test[res_test.AtomName == test_name]
                elif RES_NAME in ('GLY'):
                    if atom.name in ('H', 'HN'):
                        test_name = 'HN'
                        atom_test = res_test[res_test.AtomName == test_name]
                    if atom.name in ('HA2', 'HA3'):
                        test_name = 'HA'
                        atom_test = res_test[res_test.AtomName == test_name]
                elif RES_NAME in ('HIE','HIP','HID'):
                    if atom.name in ('N','CA','C','O','HA'):
                        test_RN = 'ALA'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('H', 'HN'):
                        test_RN = 'ALA'
                        test_name = 'HN'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('HB2', 'HB3'):
                        test_RN = RES_NAME
                        test_name = 'HB'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    else:
                        test_RN = RES_NAME
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                elif RES_NAME in ('ILE'):
                    if atom.name in ('N','CA','C','O','HA'):
                        test_RN = 'ALA'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('H', 'HN'):
                        test_RN = 'ALA'
                        test_name = 'HN'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('CD1', 'CD2'):
                        test_name = 'CD'
                        atom_test = res_test[res_test.AtomName == test_name]
                    ## Might need condition for HD
                    elif atom.name in ('HD1', 'HD2', 'HD3'):
                        test_name = 'HD'
                        atom_test = res_test[res_test.AtomName == test_name]
                    ## Deal with Amoebabio09 not labeling HG1 and HG2, but
                    ## rather relying on order
                    elif atom.name in ('HG1', 'HG2'):
                        test_RN = 'ILE'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                        if atom_test.empty == True:
                            test_name = 'HG'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == test_name]
                            ## If there's not 1 match, select the first row
                            if len(atom_test) > 1:
                                atom_test = atom_test.iloc[0]
                    elif atom.name in ('CG1', 'CG2'):
                        test_RN = 'ILE'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                        if atom_test.empty == True:
                            test_name = 'CG'
                            atom_test = res_test[res_test.AtomName == test_name]
                            ## If there's not 1 match, select the first row
                            if len(atom_test) > 1:
                                atom_test = atom_test.iloc[0]
                elif RES_NAME in ('LEU'):
                    if atom.name in ('N','CA','C','O','HA'):
                        test_RN = 'ALA'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('H', 'HN'):
                        test_RN = 'ALA'
                        test_name = 'HN'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('CD1', 'CD2'):
                        test_name = 'CD'
                        atom_test = res_test[res_test.AtomName == test_name]
                    ## Might need condition for HB
                    elif atom.name in ('HB2', 'HB3'):
                        test_name = 'HB'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('HD1', 'HD2'):
                        test_name = 'HD'
                        atom_test = res_test[res_test.AtomName == test_name]
                elif RES_NAME in ('LYN', 'LYS'):
                    if atom.name in ('N','CA','C','O','HA'):
                        test_RN = 'ALA'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('H', 'HN'):
                        test_RN = 'ALA'
                        test_name = 'HN'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('HB2', 'HB3'):
                        test_RN = RES_NAME
                        test_name = 'HB'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('HZ1', 'HZ2', 'HZ3'):
                        test_RN = RES_NAME
                        ## Why it's not HZ is yet another thing to ask the
                        ## TINKER developers
                        test_name = 'HN'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('HD2', 'HD3'):
                        test_RN = RES_NAME
                        test_name = 'HD'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('HE2', 'HE3'):
                        test_RN = RES_NAME
                        test_name = 'HE'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('HG2', 'HG3'):
                        test_RN = RES_NAME
                        test_name = 'HG'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                elif RES_NAME in ('MET'):
                    test_RN = RES_NAME
                    if atom.name in ('N','CA','C','O','HA'):
                        test_RN = 'ALA'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('H','HN'):
                        test_name = 'HN'
                        test_RN = 'ALA'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('HB2', 'HB3'):
                        test_name = 'HB'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('HE1', 'HE2', 'HE3'):
                        test_name = 'HE'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('HG2', 'HG3'):
                        test_name = 'HG'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                elif RES_NAME in ('PHE'):
                    if atom.name in ('N','CA','C','O','HA'):
                        test_RN = 'ALA'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('H','HN'):
                        test_name = 'HN'
                        test_RN = 'ALA'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('CE1', 'CE2'):
                        test_name = 'CE'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('CD1', 'CD2'):
                        test_name = 'CD'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('HB2', 'HB3'):
                        test_name = 'HB'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('HE1', 'HE2'):
                        test_name = 'HE'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('HD1', 'HD2'):
                        test_name = 'HD'
                        atom_test = res_test[res_test.AtomName == test_name]
                elif RES_NAME in ('PRO'):
                    if atom.name in ('HB2', 'HB3'):
                        test_name = 'HB'
                        atom_test = res_test[res_test.AtomName == test_name]
                    ## Might need condition for HG
                    elif atom.name in ('HG2', 'HG3'):
                            test_name = 'HG'
                            atom_test = res_test[res_test.AtomName == test_name]
                        ## Might need condition for HD
                    elif atom.name in ('HD2', 'HD3'):
                            test_name = 'HD'
                            atom_test = res_test[res_test.AtomName == test_name]
                elif RES_NAME in ('SER'):
                    if atom.name in ('N','CA','C','O','HA'):
                        test_RN = 'ALA'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('H','HN'):
                        test_name = 'HN'
                        test_RN = 'ALA'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('HB2', 'HB3'):
                        test_name = 'HB'
                        atom_test = res_test[res_test.AtomName == test_name]
                elif RES_NAME in ('THR'):
                    if atom.name in ('N','CA','C','O','HA'):
                        test_RN = 'ALA'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('H','HN'):
                        test_name = 'HN'
                        test_RN = 'ALA'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    ## Named differently between AMOEBA versions
                    elif atom.name in ('OG1'):
                        test_RN = 'THR'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                        if atom_test.empty == True:
                            test_name = 'O'
                            atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('CG2'):
                        test_RN = 'THR'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                        if atom_test.empty == True:
                            test_name = 'CG'
                            atom_test = res_test[res_test.AtomName == test_name]
                    ## Deal with Amoebabio09 not labeling HG1 and HG2, but
                    ## rather relying on order
                    elif atom.name in ('HG1', 'HG2'):
                        test_RN = 'THR'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                        if atom_test.empty == True:
                            test_name = 'HG'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == test_name]
                            ## If there's not 1 match, select the first row
                            if len(atom_test) > 1:
                                atom_test = atom_test.iloc[0]
                elif RES_NAME in ('TYD','TYR'):
                    if atom.name in ('N','CA','C','O','HA'):
                        test_RN = 'ALA'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('H', 'HN'):
                        test_RN = 'ALA'
                        test_name = 'HN'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('CD1', 'CD2'):
                        test_name = 'CD'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('HB1', 'HB2', 'HB3'):
                        test_name = 'HB'
                        atom_test = res_test[res_test.AtomName == test_name]
                    ## Might need condition for HD
                    elif atom.name in ('CE1', 'CE2'):
                        test_name = 'CE'
                        atom_test = res_test[res_test.AtomName == test_name]
                    ## Might need condition for HE
                    elif atom.name in ('HE1', 'HE2'):
                        test_name = 'HE'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('HD1', 'HD2'):
                        test_name = 'HD'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('HD1', 'HD2'):
                        test_name = 'HD'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('OH'):
                        if RES_NAME in ('TYD'):
                            test_name = 'O-'
                            atom_test = res_test[res_test.AtomName == test_name]
                        else:
                            atom_test = res_test[res_test.AtomName == test_name]
                elif RES_NAME in ('VAL'):
                    if atom.name in ('N','CA','C','O','HA'):
                        test_RN = 'ALA'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('H', 'HN'):
                        test_RN = 'ALA'
                        test_name = 'HN'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('CG1', 'CG2'):
                        test_name = 'CG'
                        atom_test = res_test[res_test.AtomName == test_name]
                    ## Might need condition for HG
                    elif atom.name in ('HG1', 'HG2'):
                        test_name = 'HG'
                        atom_test = res_test[res_test.AtomName == test_name]
            #############################################
            ##    Deal with Nucleic Acids, Water...    ##
            #############################################
            elif RES_NAME in ('WAT'):
                res_test = lines[lines.ResName == RES_NAME]
                if atom.name in ('H1','H2','H'):
                    test_name = 'H'
                    atom_test = res_test[res_test.AtomName == test_name]
                else:
                    test_name = 'O'
                    atom_test = res_test[res_test.AtomName == test_name]
            elif RES_NAME in ('DC3', 'DG3', 'DA3', 'DT3'):
                res_test = lines[lines.ResName == RES_NAME]
                if atom.name in ('O3\'T', 'O3T', 'O3\''):
                    test_name = 'O3T'
                    test_RN = 'DX3'
                    res_test = lines[lines.ResName == test_RN]
                    atom_test = res_test[res_test.AtomName == test_name]
                elif atom.name in ('HO\'3', 'H3T', 'HO3\''):
                    test_RN = 'DX'
                    atom_test = res_test[res_test.AtomName == test_name]
                    if atom_test.empty == True:
                        test_name = 'HO3\''
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    if atom_test.empty == True:
                        test_RN = 'DX3'
                        test_name = 'H3T'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                elif atom.name in ('O5\''):
                    test_RN = 'DX'
                    res_test = lines[lines.ResName == test_RN]
                    atom_test = res_test[res_test.AtomName == test_name]
                    if atom_test.empty == True:
                        if RES_NAME in ('DC3', 'DT3'):
                            test_RN = 'DCT'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == test_name]
                        else:
                            test_RN = 'DAG'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == test_name]
                elif atom.name in ('OP1', 'OP2', 'O1P', 'O2P', 'P'):
                    test_RN = 'DRP'
                    if atom.name in ('OP1', 'OP2', 'O1P', 'O2P'):
                        test_name = 'OP'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    else:
                        test_name = 'P'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    if atom_test.empty == True:
                        test_RN = 'DP'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                elif atom.name in ('N9'):
                    if RES_NAME in ('DA3'):
                        test_RN = 'DA'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == atom.name]
                        if atom_test.empty == True:
                            test_RN = 'DRA'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
                    else:
                        test_RN = 'DG'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == atom.name]
                        if atom_test.empty == True:
                            test_RN = 'DRG'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
                else:
                    test_RN = RES_NAME[:-1] #residue.name[:-1]
                    if atom.name in ('C4\'','C3\'','C2\'','C1\'',
                     'H1\'','H3\'','H4\'','C5\'','O3\'','HO\'3','O5\'',
                     'O4\''):
                        test_RN = 'DX'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == atom.name]
                        if atom_test.empty == True:
                            if residue.name in ('DC3', 'DT3'):
                                test_RN = 'DCT'
                                res_test = lines[lines.ResName == test_RN]
                                atom_test = res_test[res_test.AtomName == atom.name]
                            else:
                                test_RN = 'DAG'
                                res_test = lines[lines.ResName == test_RN]
                                atom_test = res_test[res_test.AtomName == atom.name]
                    elif atom.name in ('H2\'1', 'H2\'2'):
                        test_RN = 'DX'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                        if atom_test.empty == True:
                            test_name = 'H2\'1'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == test_name]
                            if atom_test.empty == True:
                                    if RES_NAME in ('DC3', 'DT3'):
                                        test_RN = 'DCT'
                                        res_test = lines[lines.ResName == test_RN]
                                        atom_test = res_test[res_test.AtomName == atom.name]
                                    else:
                                        test_RN = 'DAG'
                                        res_test = lines[lines.ResName == test_RN]
                                        atom_test = res_test[res_test.AtomName == atom.name]
                    elif atom.name in ('H5\'1', 'H5\'2'):
                        test_RN = 'DX'
                        test_name = 'H5\'1'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                        if atom_test.empty == True:
                            test_name = 'H5\''
                            test_RN = 'DX'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == test_name]
                            if atom_test.empty == True:
                                test_name = 'H5\'1'
                                if RES_NAME in ('DC3', 'DT3'):
                                    test_RN = 'DCT'
                                    res_test = lines[lines.ResName == test_RN]
                                    atom_test = res_test[res_test.AtomName == atom.name]
                                else:
                                    test_RN = 'DAG'
                                    res_test = lines[lines.ResName == test_RN]
                                    atom_test = res_test[res_test.AtomName == atom.name]
                    elif atom.name in ('O4'):
                        if test_RN in ('DT'):
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
                        else:
                            test_RN = 'DX'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
                    elif atom.name in ('OP1', 'OP2', 'O1P', 'O2P'):
                        test_RN = 'DRP'
                        test_name = 'OP'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('P'):
                        test_RN = 'DRP'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == atom.name]
                    elif test_RN in ('DA'):
                        test_RN = 'DRA'
                        res_test = lines[lines.ResName == test_RN]
                        if atom.name in ('H61', 'H62'):
                            test_name = 'H61'
                            atom_test = res_test[res_test.AtomName == test_name]
                        elif atom.name in ('N9'):
                            atom_test = res_test[res_test.AtomName == atom.name]
                            if atom_test.empty == True:
                                test_RN = 'DA'
                                res_test = lines[lines.ResName == test_RN]
                                atom_test = res_test[res_test.AtomName == atom.name]
                        else:
                            atom_test = res_test[res_test.AtomName == atom.name]
                    elif test_RN in ('DC'):
                        test_RN = 'DRC'
                        res_test = lines[lines.ResName == test_RN]
                        if atom.name in ('H41', 'H42'):
                            test_name = 'H41'
                            atom_test = res_test[res_test.AtomName == test_name]
                        else:
                            atom_test = res_test[res_test.AtomName == atom.name]
                        if atom.name in ('N1'):
                            atom_test = res_test[res_test.AtomName == atom.name]
                            if atom_test.empty == True:
                                test_RN = 'DC'
                                res_test = lines[lines.ResName == test_RN]
                                atom_test = res_test[res_test.AtomName == atom.name]
                    elif test_RN in ('DG'):
                        if atom.name in ('H21', 'H22'):
                            test_RN = 'DRG'
                            test_name = 'H21'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == test_name]
                        else:
                            test_RN = 'DRG'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
                    elif test_RN in ('DT'):
                        if atom.name in ('H71', 'H72', 'H73'):
                            test_name = 'H7'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == test_name]
                        else:
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
            elif RES_NAME in ('DC5', 'DG5', 'DA5', 'DT5'):
                res_test = lines[lines.ResName == RES_NAME]
                if atom.name in ('O5\'T', 'O5T', 'O5\''):
                    test_RN = 'DX5'
                    test_name = 'O5T'
                    res_test = lines[lines.ResName == test_RN]
                    atom_test = res_test[res_test.AtomName == test_name]
                elif atom.name in ('HO\'5', 'H5T', 'HO5\''):
                    test_RN = 'DX'
                    test_name = 'HO5\''
                    res_test = lines[lines.ResName == test_RN]
                    atom_test = res_test[res_test.AtomName == test_name]
                    if atom_test.empty == True:
                        test_RN = 'DX5'
                        test_name = 'H5T'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                elif atom.name in ('OP1', 'OP2', 'O1P', 'O2P', 'P'):
                    test_RN = 'DRP'
                    if atom.name in ('OP1', 'OP2', 'O1P', 'O2P'):
                        test_name = 'OP'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    else:
                        test_name = 'P'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    if atom_test.empty == True:
                        test_RN = 'DP'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                elif atom.name in ('N9'):
                    if RES_NAME in ('DA5'):
                        test_RN = 'DA'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == atom.name]
                        if atom_test.empty == True:
                            test_RN = 'DRA'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
                    else:
                        test_RN = 'DG'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == atom.name]
                        if atom_test.empty == True:
                            test_RN = 'DRG'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
                elif atom.name in ('H61', 'H62'):
                    test_RN = 'DRA'
                    test_name = 'H61'
                    res_test = lines[lines.ResName == test_RN]
                    atom_test = res_test[res_test.AtomName == test_name]
                else:
                    test_RN = RES_NAME[:-1]
                    if atom.name in ('C4\'','C3\'','C2\'','C1\'',
                     'H1\'','H3\'','H4\'','C5\'','O3\'','HO\'3','O5\'',
                     'O4\''):
                        test_RN = 'DX'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == atom.name]
                        if atom_test.empty == True:
                            if residue.name in ('DC5', 'DT5'):
                                test_RN = 'DCT'
                                res_test = lines[lines.ResName == test_RN]
                                atom_test = res_test[res_test.AtomName == atom.name]
                            else:
                                test_RN = 'DAG'
                                res_test = lines[lines.ResName == test_RN]
                                atom_test = res_test[res_test.AtomName == atom.name]
                    elif atom.name in ('H2\'1', 'H2\'2'):
                        test_RN = 'DX'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                        if atom_test.empty == True:
                            test_name = 'H2\'1'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == test_name]
                            if atom_test.empty == True:
                                    if RES_NAME in ('DC5', 'DT5'):
                                        test_RN = 'DCT'
                                        res_test = lines[lines.ResName == test_RN]
                                        atom_test = res_test[res_test.AtomName == atom.name]
                                    else:
                                        test_RN = 'DAG'
                                        res_test = lines[lines.ResName == test_RN]
                                        atom_test = res_test[res_test.AtomName == atom.name]
                    elif atom.name in ('H5\'1', 'H5\'2'):
                        test_RN = 'DX'
                        test_name = 'H5\'1'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                        if atom_test.empty == True:
                            test_name = 'H5\''
                            test_RN = 'DX'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == test_name]
                            if atom_test.empty == True:
                                if RES_NAME in ('DC5', 'DT5'):
                                    test_RN = 'DCT'
                                    res_test = lines[lines.ResName == test_RN]
                                    atom_test = res_test[res_test.AtomName == atom.name]
                                else:
                                    test_RN = 'DAG'
                                    res_test = lines[lines.ResName == test_RN]
                                    atom_test = res_test[res_test.AtomName == atom.name]
                    elif atom.name in ('O4'):
                        if test_RN in ('DT'):
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
                        else:
                            test_RN = 'DX'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
                    elif atom.name in ('OP1', 'OP2', 'O1P', 'O2P', 'P'):
                        test_RN = 'DRP'
                        if atom.name in ('OP1', 'OP2', 'O1P', 'O2P'):
                            test_name = 'OP'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == test_name]
                        else:
                            test_name = 'P'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == test_name]
                        if atom_test.empty == True:
                            test_RN = 'DP'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('HO5\'', 'H5T'):
                        test_RN = 'DX'
                        test_name = 'HO5\''
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif test_RN in ('DA'):
                        test_RN = 'DRA'
                        res_test = lines[lines.ResName == test_RN]
                        if atom.name in ('H61', 'H62'):
                            test_name = 'H61'
                            atom_test = res_test[res_test.AtomName == test_name]
                        elif atom.name in ('N9'):
                            test_RN = 'DA'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
                            if atom_test.empty == True:
                                test_RN = 'DRA'
                                res_test = lines[lines.ResName == test_RN]
                                atom_test = res_test[res_test.AtomName == atom.name]
                        else:
                            atom_test = res_test[res_test.AtomName == atom.name]
                    elif test_RN in ('DC'):
                        test_RN = 'DRC'
                        res_test = lines[lines.ResName == test_RN]
                        if atom.name in ('N1'):
                            atom_test = res_test[res_test.AtomName == atom.name]
                            if atom_test.empty == True:
                                test_RN = 'DC'
                                res_test = lines[lines.ResName == test_RN]
                                atom_test = res_test[res_test.AtomName == atom.name]
                        elif atom.name in ('H41', 'H42'):
                            test_name = 'H41'
                            atom_test = res_test[res_test.AtomName == test_name]
                        else:
                            atom_test = res_test[res_test.AtomName == atom.name]
                    elif test_RN in ('DG'):
                        test_RN = 'DRG'
                        res_test = lines[lines.ResName == test_RN]
                        if atom.name in ('N9'):
                            atom_test = res_test[res_test.AtomName == atom.name]
                            if atom_test.empty == True:
                                test_RN = 'DG'
                                res_test = lines[lines.ResName == test_RN]
                                atom_test = res_test[res_test.AtomName == atom.name]
                        elif atom.name in ('H21', 'H22'):
                            atom_test = res_test[res_test.AtomName == atom.name]
                            if atom_test.empty == True:
                                test_RN = 'DRG'
                                test_name = 'H21'
                                res_test = lines[lines.ResName == test_RN]
                                atom_test = res_test[res_test.AtomName == test_name]
                        else:
                            atom_test = res_test[res_test.AtomName == atom.name]
                    elif test_RN in ('DT'):
                        if atom.name in ('H71', 'H72', 'H73'):
                            test_name = 'H7'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == test_name]
                        else:
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
            elif RES_NAME in ('DC', 'DG', 'DA', 'DT'):
                test_RN = RES_NAME
                res_test = lines[lines.ResName == test_RN]
                if atom.name in ('C4\'','C3\'','C2\'','C1\'',
                 'H1\'','H3\'','H4\'','C5\'','O3\'','HO\'3','O5\'',
                 'HO5\'', 'O4\'', 'C4\''):
                    test_RN = 'DX'
                    res_test = lines[lines.ResName == test_RN]
                    atom_test = res_test[res_test.AtomName == atom.name]
                    if atom_test.empty == True:
                        if residue.name in ('DC', 'DT'):
                            test_RN = 'DCT'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
                        else:
                            test_RN = 'DAG'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
                elif atom.name in ('H2\'1', 'H2\'2'):
                    test_RN = 'DX'
                    res_test = lines[lines.ResName == test_RN]
                    atom_test = res_test[res_test.AtomName == test_name]
                    if atom_test.empty == True:
                        test_name = 'H2\'1'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                        if atom_test.empty == True:
                                if RES_NAME in ('DC', 'DT'):
                                    test_RN = 'DCT'
                                    res_test = lines[lines.ResName == test_RN]
                                    atom_test = res_test[res_test.AtomName == atom.name]
                                else:
                                    test_RN = 'DAG'
                                    res_test = lines[lines.ResName == test_RN]
                                    atom_test = res_test[res_test.AtomName == atom.name]
                elif atom.name in ('H5\'1', 'H5\'2'):
                    test_RN = 'DX'
                    test_name = 'H5\'1'
                    res_test = lines[lines.ResName == test_RN]
                    atom_test = res_test[res_test.AtomName == test_name]
                    if atom_test.empty == True:
                        test_name = 'H5\''
                        test_RN = 'DX'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                        if atom_test.empty == True:
                            if RES_NAME in ('DC', 'DT'):
                                test_RN = 'DCT'
                                res_test = lines[lines.ResName == test_RN]
                                atom_test = res_test[res_test.AtomName == atom.name]
                            else:
                                test_RN = 'DAG'
                                res_test = lines[lines.ResName == test_RN]
                                atom_test = res_test[res_test.AtomName == atom.name]
                elif atom.name in ('HO\'3', 'H3T', 'HO3\''):
                    test_RN = 'DX'
                    res_test = lines[lines.ResName == test_RN]
                    atom_test = res_test[res_test.AtomName == test_name]
                    if atom_test.empty == True:
                        test_name = 'HO3\''
                        atom_test = res_test[res_test.AtomName == test_name]
                elif atom.name in ('HO\'5', 'H5T', 'HO5\''):
                    test_name = 'HO5\''
                    test_RN = 'DX'
                    res_test = lines[lines.ResName == test_RN]
                    atom_test = res_test[res_test.AtomName == test_name]
                elif atom.name in ('O4'):
                    if test_RN in ('DT'):
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == atom.name]
                    else:
                        test_RN = 'DX'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == atom.name]
                elif atom.name in ('OP1', 'OP2', 'O1P', 'O2P', 'P'):
                    test_RN = 'DRP'
                    if atom.name in ('OP1', 'OP2', 'O1P', 'O2P'):
                        test_name = 'OP'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    else:
                        test_name = 'P'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    if atom_test.empty == True:
                        test_RN = 'DP'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                elif test_RN in ('DA'):
                    test_RN = 'DRA'
                    res_test = lines[lines.ResName == test_RN]
                    if atom.name in ('H61', 'H62'):
                        test_name = 'H61'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('N9'):
                        atom_test = res_test[res_test.AtomName == atom.name]
                        if atom_test.empty == True:
                            test_RN = 'DA'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
                    else:
                        atom_test = res_test[res_test.AtomName == atom.name]
                elif test_RN in ('DC'):
                    test_RN = 'DRC'
                    res_test = lines[lines.ResName == test_RN]
                    if atom.name in ('N1'):
                        atom_test = res_test[res_test.AtomName == atom.name]
                        if atom_test.empty == True:
                            test_RN = 'DC'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
                    elif atom.name in ('H41', 'H42'):
                        test_name = 'H41'
                        atom_test = res_test[res_test.AtomName == test_name]
                    else:
                        atom_test = res_test[res_test.AtomName == atom.name]
                elif test_RN in ('DG'):
                    test_RN = 'DRG'
                    res_test = lines[lines.ResName == test_RN]
                    if atom.name in ('N9'):
                        atom_test = res_test[res_test.AtomName == atom.name]
                        if atom_test.empty == True:
                            test_RN = 'DG'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
                    elif atom.name in ('H21', 'H22'):
                        atom_test = res_test[res_test.AtomName == atom.name]
                        if atom_test.empty == True:
                            test_RN = 'DRG'
                            test_name = 'H21'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == test_name]
                    else:
                        atom_test = res_test[res_test.AtomName == atom.name]
                elif test_RN in ('DT'):
                    res_test = lines[lines.ResName == test_RN]
                    if atom.name in ('H71', 'H72', 'H73'):
                        test_name = 'H7'
                        atom_test = res_test[res_test.AtomName == test_name]
                    else:
                        atom_test = res_test[res_test.AtomName == atom.name]
                else:
                    res_test = lines[lines.ResName == test_RN]
                    atom_test = res_test[res_test.AtomName == atom.name]
            ## Deal with RNA now
            elif RES_NAME in ('RC3', 'RG3', 'RA3', 'RU3'):
                res_test = lines[lines.ResName == RES_NAME]
                if atom.name in ('O3\'T', 'O3T', 'O3\''):
                    test_name = 'O3T'
                    test_RN = 'RX3'
                    res_test = lines[lines.ResName == test_RN]
                    atom_test = res_test[res_test.AtomName == test_name]
                elif atom.name in ('HO\'3', 'H3T', 'HO3\''):
                    test_RN = 'RX'
                    test_name = 'HO3\''
                    res_test = lines[lines.ResName == test_RN]
                    atom_test = res_test[res_test.AtomName == test_name]
                    if atom_test.empty == True:
                        test_name = 'HO3\''
                        atom_test = res_test[res_test.AtomName == test_name]
                    if atom_test.empty == True:
                        test_RN = 'RX3'
                        test_name = 'H3T'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                elif atom.name in ('O5\'', 'HO\'2'):
                    test_RN = 'RX'
                    res_test = lines[lines.ResName == test_RN]
                    atom_test = res_test[res_test.AtomName == test_name]
                    if atom_test.empty == True:
                        if RES_NAME in ('RC3', 'RU3'):
                            test_RN = 'RCU'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == test_name]
                        else:
                            test_RN = 'RAG'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == test_name]
                elif atom.name in ('H5\'1', 'H5\'2'):
                    test_RN = 'RX'
                    test_name = 'H5\'1'
                    res_test = lines[lines.ResName == test_RN]
                    atom_test = res_test[res_test.AtomName == test_name]
                    if atom_test.empty == True:
                        test_name = 'H5\''
                        test_RN = 'RX'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                        if atom_test.empty == True:
                                if RES_NAME in ('RC3', 'RU3'):
                                    test_RN = 'RCU'
                                    res_test = lines[lines.ResName == test_RN]
                                    atom_test = res_test[res_test.AtomName == atom.name]
                                else:
                                    test_RN = 'RAG'
                                    res_test = lines[lines.ResName == test_RN]
                                    atom_test = res_test[res_test.AtomName == atom.name]
                elif atom.name in ('OP1', 'OP2', 'O1P', 'O2P', 'P'):
                    test_RN = 'DRP'
                    if atom.name in ('OP1', 'OP2', 'O1P', 'O2P'):
                        test_name = 'OP'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    else:
                        test_name = 'P'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    if atom_test.empty == True:
                        test_RN = 'RP'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                elif atom.name in ('N9'):
                    if RES_NAME in ('RA3'):
                        test_RN = 'RA'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == atom.name]
                        if atom_test.empty == True:
                            test_RN = 'DRA'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
                    else:
                        test_RN = 'RG'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == atom.name]
                        if atom_test.empty == True:
                            test_RN = 'DRG'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
                elif atom.name in ('H21', 'H22'):
                    test_RN = 'DRG'
                    atom_test = res_test[res_test.AtomName == atom.name]
                    if atom_test.empty == True:
                        test_RN = 'DRG'
                        test_name = 'H21'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                else:
                    test_RN = RES_NAME[:-1] #residue.name[:-1]
                    if atom.name in ('C4\'','C3\'','C2\'','C1\'','H2\'1',
                     'H1\'','H3\'','H4\'','C5\'','O3\'','HO\'3','O5\'',
                     'O4\'', 'O2\''):
                        test_RN = 'RX'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == atom.name]
                        if atom_test.empty == True:
                            if residue.name in ('RC3', 'RU3'):
                                test_RN = 'RCU'
                                res_test = lines[lines.ResName == test_RN]
                                atom_test = res_test[res_test.AtomName == atom.name]
                            else:
                                test_RN = 'RAG'
                                res_test = lines[lines.ResName == test_RN]
                                atom_test = res_test[res_test.AtomName == atom.name]
                    elif atom.name in ('O4'):
                        if test_RN in ('RU'):
                            test_RN = 'RU'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
                        else:
                            test_RN = 'RX'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
                    elif atom.name in ('OP1', 'OP2', 'O1P', 'O2P', 'P'):
                        test_RN = 'DRP'
                        if atom.name in ('OP1', 'OP2', 'O1P', 'O2P'):
                            test_name = 'OP'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == test_name]
                        else:
                            test_name = 'P'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == test_name]
                        if atom_test.empty == True:
                            test_RN = 'RP'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == test_name]
                    elif test_RN in ('RA'):
                        test_RN = 'DRA'
                        res_test = lines[lines.ResName == test_RN]
                        if atom.name in ('H61', 'H62'):
                            test_name = 'H61'
                            atom_test = res_test[res_test.AtomName == test_name]
                        elif atom.name in ('N9'):
                            atom_test = res_test[res_test.AtomName == atom.name]
                            if atom_test.empty == True:
                                test_RN = 'RA'
                                res_test = lines[lines.ResName == test_RN]
                                atom_test = res_test[res_test.AtomName == atom.name]
                        else:
                            atom_test = res_test[res_test.AtomName == atom.name]
                    elif test_RN in ('RC'):
                        test_RN = 'DRC'
                        res_test = lines[lines.ResName == test_RN]
                        if atom.name in ('N1'):
                            atom_test = res_test[res_test.AtomName == atom.name]
                            if atom_test.empty == True:
                                test_RN = 'RC'
                                res_test = lines[lines.ResName == test_RN]
                                atom_test = res_test[res_test.AtomName == atom.name]
                        elif atom.name in ('H41', 'H42'):
                            test_name = 'H41'
                            atom_test = res_test[res_test.AtomName == test_name]
                        else:
                            atom_test = res_test[res_test.AtomName == atom.name]
                    elif test_RN in ('RG'):
                        test_RN = 'DRG'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == atom.name]
                    elif test_RN in ('RU'):
                        test_RN = 'RU'
                        if atom.name in ('H71', 'H72', 'H73'):
                            test_name = 'H7'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == test_name]
                        else:
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
            elif RES_NAME in ('RC5', 'RG5', 'RA5', 'RU5'):
                res_test = lines[lines.ResName == RES_NAME]
                if atom.name in ('O5\'T', 'O5T', 'O5\''):
                    test_RN = 'RX5'
                    test_name = 'O5T'
                    res_test = lines[lines.ResName == test_RN]
                    atom_test = res_test[res_test.AtomName == test_name]
                elif atom.name in ('HO\'5', 'H5T', 'HO5\''):
                    test_RN = 'RX'
                    test_name = 'HO5\''
                    res_test = lines[lines.ResName == test_RN]
                    atom_test = res_test[res_test.AtomName == test_name]
                    if atom_test.empty == True:
                        test_RN = 'RX5'
                        test_name = 'H5T'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                elif atom.name in ('H5\'1', 'H5\'2'):
                    test_RN = 'RX'
                    test_name = 'H5\'1'
                    res_test = lines[lines.ResName == test_RN]
                    atom_test = res_test[res_test.AtomName == test_name]
                    if atom_test.empty == True:
                        test_name = 'H5\''
                        test_RN = 'RX'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                        if atom_test.empty == True:
                                if RES_NAME in ('RC5', 'RU5'):
                                    test_RN = 'RCU'
                                    res_test = lines[lines.ResName == test_RN]
                                    atom_test = res_test[res_test.AtomName == atom.name]
                                else:
                                    test_RN = 'RAG'
                                    res_test = lines[lines.ResName == test_RN]
                                    atom_test = res_test[res_test.AtomName == atom.name]
                elif atom.name in ('OP1', 'OP2', 'O1P', 'O2P', 'P'):
                    test_RN = 'DRP'
                    if atom.name in ('OP1', 'OP2', 'O1P', 'O2P'):
                        test_name = 'OP'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    else:
                        test_name = 'P'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    if atom_test.empty == True:
                        test_RN = 'RP'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                elif atom.name in ('N9'):
                    if RES_NAME in ('RA5'):
                        test_RN = 'RA'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == atom.name]
                        if atom_test.empty == True:
                            test_RN = 'DRA'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == test_name]
                    else:
                        test_RN = 'RG'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == atom.name]
                        if atom_test.empty == True:
                            test_RN = 'DRG'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == test_name]
                elif atom.name in ('H61', 'H62'):
                    test_RN = 'DRA'
                    test_name = 'H61'
                    res_test = lines[lines.ResName == test_RN]
                    atom_test = res_test[res_test.AtomName == test_name]
                elif atom.name in ('H21', 'H22'):
                    test_RN = 'DRG'
                    atom_test = res_test[res_test.AtomName == atom.name]
                    if atom_test.empty == True:
                        test_RN = 'DRG'
                        test_name = 'H21'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                else:
                    test_RN = RES_NAME[:-1]
                    if atom.name in ('C4\'','C3\'','C2\'','C1\'','H2\'1',
                     'H1\'','H3\'','H4\'','C5\'','O3\'','HO\'3','O5\'',
                     'O4\'', 'HO\'2', 'O2\''):
                        test_RN = 'RX'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == atom.name]
                        if atom_test.empty == True:
                            if residue.name in ('RC5', 'RU5'):
                                test_RN = 'RCU'
                                res_test = lines[lines.ResName == test_RN]
                                atom_test = res_test[res_test.AtomName == atom.name]
                            else:
                                test_RN = 'RAG'
                                res_test = lines[lines.ResName == test_RN]
                                atom_test = res_test[res_test.AtomName == atom.name]
                    elif atom.name in ('O4'):
                        if test_RN in ('RU'):
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
                        else:
                            test_RN = 'RX'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
                    elif atom.name in ('OP1', 'OP2', 'O1P', 'O2P', 'P'):
                        test_RN = 'DRP'
                        if atom.name in ('OP1', 'OP2', 'O1P', 'O2P'):
                            test_name = 'OP'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == test_name]
                        else:
                            test_name = 'P'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == test_name]
                        if atom_test.empty == True:
                            test_RN = 'RP'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('HO5\'', 'H5T'):
                        test_RN = 'RX'
                        test_name = 'HO5\''
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif test_RN in ('RA'):
                        test_RN = 'DRA'
                        res_test = lines[lines.ResName == test_RN]
                        if atom.name in ('H61', 'H62'):
                            test_name = 'H61'
                            atom_test = res_test[res_test.AtomName == test_name]
                        elif atom.name in ('N9'):
                            atom_test = res_test[res_test.AtomName == atom.name]
                            if atom_test.empty == True:
                                test_RN = 'RA'
                                res_test = lines[lines.ResName == test_RN]
                                atom_test = res_test[res_test.AtomName == atom.name]
                        else:
                            atom_test = res_test[res_test.AtomName == atom.name]
                    elif test_RN in ('RC'):
                        test_RN = 'DRC'
                        res_test = lines[lines.ResName == test_RN]
                        if atom.name in ('N1'):
                            atom_test = res_test[res_test.AtomName == atom.name]
                            if atom_test.empty == True:
                                test_RN = 'RC'
                                res_test = lines[lines.ResName == test_RN]
                                atom_test = res_test[res_test.AtomName == atom.name]
                        elif atom.name in ('H41', 'H42'):
                            test_name = 'H41'
                            atom_test = res_test[res_test.AtomName == test_name]
                        else:
                            atom_test = res_test[res_test.AtomName == atom.name]
                    elif test_RN in ('RG'):
                        test_RN = 'DRG'
                        res_test = lines[lines.ResName == test_RN]
                        if atom.name in ('N9'):
                            atom_test = res_test[res_test.AtomName == atom.name]
                            if atom_test.empty == True:
                                test_RN = 'RG'
                                res_test = lines[lines.ResName == test_RN]
                                atom_test = res_test[res_test.AtomName == atom.name]
                        else:
                            atom_test = res_test[res_test.AtomName == atom.name]
                    elif test_RN in ('RU'):
                        test_RN = 'RU'
                        if atom.name in ('H71', 'H72', 'H73'):
                            test_name = 'H7'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == test_name]
                        else:
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
            elif RES_NAME in ('RC', 'RG', 'RA', 'RU'):
                test_RN = RES_NAME
                res_test = lines[lines.ResName == test_RN]
                if atom.name in ('C4\'','C3\'','C2\'','C1\'','H2\'1',
                 'H1\'','H3\'','H4\'','C5\'','O3\'','HO\'3','O5\'',
                 'HO5\'', 'O4\'', 'C4\'', 'HO\'2', 'O2\''):
                    test_RN = 'RX'
                    res_test = lines[lines.ResName == test_RN]
                    atom_test = res_test[res_test.AtomName == atom.name]
                    if atom_test.empty == True:
                        if residue.name in ('RC', 'RU'):
                            test_RN = 'RCU'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
                        else:
                            test_RN = 'RAG'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
                elif atom.name in ('HO\'3', 'H3T', 'HO3\''):
                    test_RN = 'RX'
                    res_test = lines[lines.ResName == test_RN]
                    atom_test = res_test[res_test.AtomName == test_name]
                    if atom_test.empty == True:
                        test_name = 'HO3\''
                        atom_test = res_test[res_test.AtomName == test_name]
                elif atom.name in ('HO\'5', 'H5T', 'HO5\''):
                    test_name = 'HO5\''
                    test_RN = 'RX'
                    res_test = lines[lines.ResName == test_RN]
                    atom_test = res_test[res_test.AtomName == test_name]
                elif atom.name in ('O4'):
                    if test_RN in ('RU'):
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == atom.name]
                    else:
                        test_RN = 'RX'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == atom.name]
                elif atom.name in ('H5\'1', 'H5\'2'):
                    test_RN = 'RX'
                    test_name = 'H5\'1'
                    res_test = lines[lines.ResName == test_RN]
                    atom_test = res_test[res_test.AtomName == test_name]
                    if atom_test.empty == True:
                        test_name = 'H5\''
                        test_RN = 'RX'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                        if atom_test.empty == True:
                                if RES_NAME in ('RC', 'RU'):
                                    test_RN = 'RCU'
                                    res_test = lines[lines.ResName == test_RN]
                                    atom_test = res_test[res_test.AtomName == atom.name]
                                else:
                                    test_RN = 'RAG'
                                    res_test = lines[lines.ResName == test_RN]
                                    atom_test = res_test[res_test.AtomName == atom.name]
                elif atom.name in ('OP1', 'OP2', 'O1P', 'O2P', 'P'):
                    test_RN = 'DRP'
                    if atom.name in ('OP1', 'OP2', 'O1P', 'O2P'):
                        test_name = 'OP'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    else:
                        test_name = 'P'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                    if atom_test.empty == True:
                        test_RN = 'RP'
                        res_test = lines[lines.ResName == test_RN]
                        atom_test = res_test[res_test.AtomName == test_name]
                elif test_RN in ('RA'):
                    test_RN = 'DRA'
                    res_test = lines[lines.ResName == test_RN]
                    if atom.name in ('H61', 'H62'):
                        test_name = 'H61'
                        atom_test = res_test[res_test.AtomName == test_name]
                    elif atom.name in ('N9'):
                        atom_test = res_test[res_test.AtomName == atom.name]
                        if atom_test.empty == True:
                            test_RN = 'RA'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
                    else:
                        atom_test = res_test[res_test.AtomName == atom.name]
                elif test_RN in ('RC'):
                    test_RN = 'DRC'
                    res_test = lines[lines.ResName == test_RN]
                    if atom.name in ('N1'):
                        atom_test = res_test[res_test.AtomName == atom.name]
                        if atom_test.empty == True:
                            test_RN = 'RC'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
                    elif atom.name in ('H41', 'H42'):
                        test_name = 'H41'
                        atom_test = res_test[res_test.AtomName == test_name]
                    else:
                        atom_test = res_test[res_test.AtomName == atom.name]
                elif test_RN in ('RG'):
                    test_RN = 'DRG'
                    res_test = lines[lines.ResName == test_RN]
                    if atom.name in ('N9'):
                        atom_test = res_test[res_test.AtomName == atom.name]
                        if atom_test.empty == True:
                            test_RN = 'RG'
                            res_test = lines[lines.ResName == test_RN]
                            atom_test = res_test[res_test.AtomName == atom.name]
                    elif atom.name in ('H21', 'H22'):
                        atom_test = res_test[res_test.AtomName == atom.name]
                        if atom_test.empty == True:
                            test_name = 'H21'
                            atom_test = res_test[res_test.AtomName == test_name]
                    else:
                        atom_test = res_test[res_test.AtomName == atom.name]
                elif test_RN in ('RU'):
                    res_test = lines[lines.ResName == test_RN]
                    atom_test = res_test[res_test.AtomName == atom.name]
                else:
                    res_test = lines[lines.ResName == test_RN]
                    atom_test = res_test[res_test.AtomName == atom.name]

            ## Now that we've renamed the relevant ones, skip these before
            ## getting to things like water and ions
            elif RES_NAME in ('NALA', 'NARG', 'NASH', 'NASN', 'NASP', 'NCYS',
             'NCYM', 'NCYX', 'NGLH', 'NGLN', 'NGLU', 'NGLY', 'NHID', 'NHIE',
             'NHIP', 'NHIS', 'NILE', 'NLEU', 'NLYN', 'NLYS', 'NMET', 'NPHE',
             'NPRO', 'NSER', 'NTHR', 'NTRP', 'NTYD', 'NTYR', 'NVAL',
             'CALA', 'CARG', 'CASN', 'CASH', 'CASP', 'CCYS',
             'CCYM', 'CCYX', 'CGLH', 'CGLN', 'CGLU', 'CGLY', 'CHID', 'CHIE',
             'CHIP', 'CHIS', 'CILE', 'CLEU', 'CLYN', 'CLYS', 'CMET', 'CPHE',
             'CPRO', 'CSER', 'CTHR', 'CTRP', 'CTYD', 'CTYR', 'CVAL'):
                pass
            else:
                test_name = atom.name
                res_test = lines[lines.ResName == RES_NAME]
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

system = load_pdb(infile, AMOEBA)

if AMOEBA == True:
    system = convert_names_AMOEBA(system, lines)
else:
    # system = convert_names(system, lines)
    print("Sorry, Charlie.")

write_xyz(system, outfile)

# ## Write out the XYZ
# system.save(outfile, overwrite=True)

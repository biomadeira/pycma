﻿#!/usr/bin/env python
# -*- coding: UTF-8 -*-                                                      
# Copyright (c) 2012 F. Madeira and L. Krippahl, 2012                                            
# This code is part of PyCMA distribution.                                 
# GNU General Public License - See LICENSE for more details.

#-----------------------------------------------------------------------

"""
Utils Library.
"""

"""
This file is part of PyCMA.

    pycma is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published 
    by the Free Software Foundation, either version 3 of the License, 
    or (at your option) any later version.

    pycma is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this file. If not, see <http://www.gnu.org/licenses/>.
"""

#-----------------------------------------------------------------------

aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'K', 'I', 'L',
      'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'Y', 'W', '-']

#-----------------------------------------------------------------------
# Amino Acid Reduction Alphabets - Adapted from
# http://bio.math-inf.uni-greifswald.de/viscose/html/alphabets.html
# Nomenclature for alphabet names: name_N(_V) (N = number of groups of
# amino acids; V = variation)

# Chan H.S., Dill K.A. (1989) Compact polymers. Macromolecules, 22, 
# 4559-4573
# Lau K.F., Dill K.A. (1989) A lattice statistical mechanics model of
# the conformational and sequence spaces of proteins. Macromolecules, 
# 22, 3986-3997
# AGTSNQDEHRKP - hydrophilic 'A'
# CMFILVWY - hydrophobic 'C'
chemical_2_1 = {'A':'A', 'C':'C', 'D':'A',
                'E':'A', 'F':'C', 'G':'A',
                'H':'A', 'K':'A', 'I':'C',
                'L':'C', 'M':'C', 'N':'A',
                'P':'A', 'Q':'A', 'R':'A',
                'S':'A', 'T':'A', 'V':'C',
                'Y':'C', 'W':'C', '-':'-'}

# Adapted from http://www.russelllab.org/aas/aas.html
# IVL - Aliphatic 'I'
# FYWH - Aromatic 'F'
# KRDE - Charged 'K'
# GACS - Tiny 'G'
# TMQNP - Diverse 'T'
chemical_5 = {'A':'G', 'C':'G', 'D':'K',
              'E':'K', 'F':'F', 'G':'G',
              'H':'F', 'K':'K', 'I':'I',
              'L':'I', 'M':'T', 'N':'T',
              'P':'T', 'Q':'T', 'R':'K',
              'S':'G', 'T':'T', 'V':'I',
              'Y':'F', 'W':'F', '-':'-'}

# Adapted from http://www.russelllab.org/aas/aas.html
# IVL - Aliphatic 'I'
# FYWH - Aromatic 'F'
# KR - Pos. charged 'K'
# DE - Neg. charged 'D'
# GACS - Tiny 'G'
# TMQNP - Diverse 'T'
chemical_6 = {'A':'G', 'C':'G', 'D':'D',
              'E':'D', 'F':'F', 'G':'G',
              'H':'F', 'K':'K', 'I':'I',
              'L':'I', 'M':'T', 'N':'T',
              'P':'T', 'Q':'T', 'R':'K',
              'S':'G', 'T':'T', 'V':'I',
              'Y':'F', 'W':'F', '-':'-'}

# Caporaso, J. G., Smit, S., Easton, B. C., Hunter, L., Huttley, G. a, 
# Knight, R. (2008). Detecting coevolution without phylogenetic trees? 
# Tree-ignorant metrics of coevolution perform as well as tree-aware 
# metrics. BMC evolutionary biology, 8, 327. doi:10.1186/1471-2148-8-327
# CHARGE_2 
chemical_2_2 = {'A':'A', 'C':'A', 'D':'K',
                'E':'K', 'F':'A', 'G':'A',
                'H':'A', 'K':'K', 'I':'A',
                'L':'A', 'M':'A', 'N':'A',
                'P':'A', 'Q':'A', 'R':'K',
                'S':'A', 'T':'A', 'V':'A',
                'Y':'A', 'W':'A', '-':'-'}

# CHARGE_HIS_2
chemical_2_3 = {'A':'A', 'C':'A', 'D':'K',
                'E':'K', 'F':'A', 'G':'A',
                'H':'K', 'K':'K', 'I':'A',
                'L':'A', 'M':'A', 'N':'A',
                'P':'A', 'Q':'A', 'R':'K',
                'S':'A', 'T':'A', 'V':'A',
                'Y':'A', 'W':'A', '-':'-'}

# SIZE_2
chemical_2_4 = {'A':'G', 'C':'G', 'D':'G',
                'E':'M', 'F':'M', 'G':'G',
                'H':'M', 'K':'M', 'I':'G',
                'L':'G', 'M':'M', 'N':'G',
                'P':'G', 'Q':'M', 'R':'M',
                'S':'G', 'T':'G', 'V':'G',
                'Y':'M', 'W':'M', '-':'-'}

# CHARGE_3
chemical_3_1 = {'A':'A', 'C':'A', 'D':'D',
                'E':'D', 'F':'A', 'G':'A',
                'H':'A', 'K':'K', 'I':'A',
                'L':'A', 'M':'A', 'N':'A',
                'P':'A', 'Q':'A', 'R':'K',
                'S':'A', 'T':'A', 'V':'A',
                'Y':'A', 'W':'A', '-':'-'}

# CHARGE_HIS_3
chemical_3_2 = {'A':'A', 'C':'A', 'D':'D',
                'E':'D', 'F':'A', 'G':'A',
                'H':'K', 'K':'K', 'I':'A',
                'L':'A', 'M':'A', 'N':'A',
                'P':'A', 'Q':'A', 'R':'K',
                'S':'A', 'T':'A', 'V':'A',
                'Y':'A', 'W':'A', '-':'-'}

# HYDROPATHY_3 
chemical_3_3 = {'A':'A', 'C':'A', 'D':'D',
                'E':'D', 'F':'A', 'G':'G',
                'H':'D', 'K':'D', 'I':'A',
                'L':'A', 'M':'A', 'N':'D',
                'P':'A', 'Q':'D', 'R':'D',
                'S':'G', 'T':'G', 'V':'A',
                'Y':'G', 'W':'G', '-':'-'}

# POLARITY_HIS_4
chemical_4 = {'A':'A', 'C':'G', 'D':'D',
              'E':'D', 'F':'A', 'G':'G',
              'H':'K', 'K':'K', 'I':'A',
              'L':'A', 'M':'A', 'N':'G',
              'P':'A', 'Q':'G', 'R':'K',
              'S':'G', 'T':'G', 'V':'A',
              'Y':'G', 'W':'A', '-':'-'}

# Murphy L.R., Wallqvist A, Levy RM. (2000) Simplified amino acid 
# alphabets for protein fold recognition and implications for folding. 
# Protein Eng. 13(3):149-152
murphy_2 = {'A':'L', 'C':'L', 'D':'E',
            'E':'E', 'F':'L', 'G':'L',
            'H':'E', 'K':'E', 'I':'L',
            'L':'L', 'M':'L', 'N':'E',
            'P':'L', 'Q':'E', 'R':'E',
            'S':'L', 'T':'L', 'V':'L',
            'Y':'L', 'W':'L', '-':'-'}

murphy_4 = {'A':'A', 'C':'L', 'D':'E',
            'E':'E', 'F':'F', 'G':'A',
            'H':'E', 'K':'E', 'I':'L',
            'L':'L', 'M':'L', 'N':'E',
            'P':'A', 'Q':'E', 'R':'E',
            'S':'A', 'T':'A', 'V':'L',
            'Y':'F', 'W':'F', '-':'-'}
            
murphy_8 = {'A':'A', 'C':'L', 'D':'E',
            'E':'E', 'F':'F', 'G':'A',
            'H':'H', 'K':'K', 'I':'L',
            'L':'L', 'M':'L', 'N':'E',
            'P':'P', 'Q':'E', 'R':'K',
            'S':'S', 'T':'S', 'V':'L',
            'Y':'F', 'W':'F', '-':'-'}
            
murphy_10 = {'A':'A', 'C':'C', 'D':'E',
            'E':'E', 'F':'F', 'G':'G',
            'H':'H', 'K':'K', 'I':'L',
            'L':'L', 'M':'L', 'N':'E',
            'P':'P', 'Q':'E', 'R':'K',
            'S':'S', 'T':'S', 'V':'L',
            'Y':'F', 'W':'F', '-':'-'}
            
murphy_15 = {'A':'A', 'C':'C', 'D':'D',
            'E':'E', 'F':'F', 'G':'G',
            'H':'H', 'K':'K', 'I':'L',
            'L':'L', 'M':'L', 'N':'N',
            'P':'P', 'Q':'Q', 'R':'K',
            'S':'S', 'T':'T', 'V':'L',
            'Y':'F', 'W':'W', '-':'-'}
            
# Wang J., Wang W. (1999) A computational approach to simplifying the 
# protein folding alphabet. Nat Struct Biol. (11):1033-1038.
wang_2 = {'A':'A', 'C':'C', 'D':'A',
          'E':'A', 'F':'C', 'G':'A',
          'H':'A', 'K':'A', 'I':'C',
          'L':'C', 'M':'C', 'N':'A',
          'P':'A', 'Q':'A', 'R':'A',
          'S':'A', 'T':'A', 'V':'C',
          'Y':'C', 'W':'C', '-':'-'}

wang_3 = {'A':'A', 'C':'C', 'D':'D',
          'E':'D', 'F':'C', 'G':'A',
          'H':'A', 'K':'D', 'I':'C',
          'L':'C', 'M':'C', 'N':'D',
          'P':'A', 'Q':'D', 'R':'A',
          'S':'D', 'T':'A', 'V':'C',
          'Y':'C', 'W':'C', '-':'-'}

wang_5_1 = {'A':'A', 'C':'C', 'D':'D',
            'E':'D', 'F':'C', 'G':'G',
            'H':'A', 'K':'S', 'I':'C',
            'L':'C', 'M':'C', 'N':'S',
            'P':'G', 'Q':'S', 'R':'S',
            'S':'S', 'T':'A', 'V':'C',
            'Y':'C', 'W':'C', '-':'-'}
            
wang_5_2 = {'A':'A', 'C':'C', 'D':'N',
            'E':'N', 'F':'C', 'G':'A',
            'H':'H', 'K':'H', 'I':'C',
            'L':'L', 'M':'C', 'N':'N',
            'P':'H', 'Q':'N', 'R':'H',
            'S':'A', 'T':'A', 'V':'L',
            'Y':'L', 'W':'L', '-':'-'}

# Li T., Fan K., Wang J., Wang W. (2003) Reduction of protein sequence 
# complexity by residue grouping. Protein Eng. (5):323-330.
li_3 = {'A':'G', 'C':'C', 'D':'N',
        'E':'N', 'F':'C', 'G':'G',
        'H':'N', 'K':'N', 'I':'C',
        'L':'C', 'M':'C', 'N':'N',
        'P':'G', 'Q':'N', 'R':'N',
        'S':'G', 'T':'G', 'V':'C',
        'Y':'C', 'W':'C', '-':'-'}

li_4 = {'A':'G', 'C':'C', 'D':'N',
        'E':'N', 'F':'C', 'G':'G',
        'H':'N', 'K':'N', 'I':'M',
        'L':'M', 'M':'M', 'N':'N',
        'P':'G', 'Q':'N', 'R':'N',
        'S':'G', 'T':'G', 'V':'M',
        'Y':'C', 'W':'C', '-':'-'}

li_5 = {'A':'P', 'C':'C', 'D':'N',
        'E':'N', 'F':'C', 'G':'G',
        'H':'N', 'K':'N', 'I':'M',
        'L':'M', 'M':'M', 'N':'N',
        'P':'P', 'Q':'N', 'R':'N',
        'S':'P', 'T':'P', 'V':'M',
        'Y':'C', 'W':'C', '-':'-'}

li_10 = {'A':'A', 'C':'C', 'D':'Q',
         'E':'Q', 'F':'F', 'G':'G',
         'H':'N', 'K':'R', 'I':'I',
         'L':'M', 'M':'M', 'N':'N',
         'P':'P', 'Q':'Q', 'R':'R',
         'S':'A', 'T':'A', 'V':'I',
         'Y':'F', 'W':'F', '-':'-'}


# Extended for Alphabets Reduction - Probably no the best
# approach
# B = "Asx";  Aspartic acid (R) or Asparagine (N)
# X = "Xxx";  Unknown or 'other' amino acid
# Z = "Glx";  Glutamic acid (E) or Glutamine (Q)
# J = "Xle";  Leucine (L) or Isoleucine (I), used in mass-spec (NMR)
# U = "Sec";  Selenocysteine
# O = "Pyl";  Pyrrolysine
# . = Gap
# - = Gap
extended = {'A':'A', 'C':'C', 'D':'D',
            'E':'E', 'F':'F', 'G':'G',
            'H':'H', 'K':'K', 'I':'I',
            'L':'L', 'M':'M', 'N':'N',
            'P':'P', 'Q':'Q', 'R':'R',
            'S':'S', 'T':'T', 'V':'V',
            'Y':'Y', 'W':'W', '-':'-',
            '.':'-', 'B':'-', 'X':'-',
            'Z':'-', 'J':'-', 'U':'-',
            'O':'-'}

#-----------------------------------------------------------------------
# Gao, H., Dou, Y., Yang, J., & Wang, J. (2011). New methods to measure 
#residues coevolution in proteins. BMC bioinformatics, 12(1), 206. 
#hydrophobic = ['A', 'G', 'C', 'T', 'I', 'V', 'L', 'K', 'H', 'F', 'Y', 'W', 'M'] 
#aromatic = ['F', 'Y', 'W', 'H']
#aliphatic = ['I', 'V', 'L']
#tiny = ['A', 'S', 'G', 'C']
#small = ['P', 'N', 'D', 'T', 'C', 'A', 'G', 'S', 'V']
#proline = ['P']
#charged = ['K', 'H', 'R', 'D', 'E']
#negative = ['D', 'E']
#polar = ['N', 'Q', 'S', 'D', 'E', 'C', 'T', 'K', 'R', 'H', 'Y', 'W']
#positive = ['K', 'H', 'R']

# converted to
A = ['A', 'G', 'C', 'T', 'I', 'V', 'L', 'K', 'H', 'F', 'Y', 'W', 'M', 'S', 'P', 'N', 'D']
R = ['R', 'K', 'H', 'D', 'E', 'N', 'Q', 'S', 'C', 'T', 'Y', 'W'] 
N = ['N', 'P', 'D', 'T', 'C', 'A', 'G', 'S', 'V', 'Q', 'E', 'K', 'R', 'H', 'Y', 'W']
D = ['D', 'P', 'N', 'T', 'C', 'A', 'G', 'S', 'V', 'E', 'Q', 'K', 'R', 'H', 'Y', 'W']
C = ['C', 'G', 'A', 'T', 'I', 'V', 'L', 'K', 'H', 'F', 'Y', 'W', 'M', 'S', 'N', 'Q', 'D', 'P', 'N', 'E']
Q = ['Q', 'N', 'S', 'D', 'E', 'C', 'T', 'K', 'R', 'H', 'Y', 'W']
E = ['E', 'D', 'K', 'H', 'R', 'N', 'Q', 'S', 'C', 'T', 'R', 'H', 'Y', 'W']
G = ['G', 'A', 'C', 'T', 'I', 'V', 'L', 'K', 'H', 'F', 'Y', 'W', 'M', 'S', 'P', 'N', 'D', 'T', 'V']
H = ['H', 'A', 'G', 'C', 'T', 'I', 'V', 'L', 'K', 'F', 'Y', 'W', 'M', 'R', 'D', 'E', 'S', 'N', 'Q']
I = ['I', 'A', 'G', 'C', 'T', 'V', 'L', 'K', 'H', 'F', 'Y', 'W', 'M']
L = ['L', 'A', 'G', 'C', 'T', 'V', 'I', 'K', 'H', 'F', 'Y', 'W', 'M']
K = ['K', 'G', 'C', 'T', 'I', 'V', 'L', 'A', 'H', 'F', 'Y', 'W', 'M', 'R', 'D', 'E', 'N', 'Q', 'S']
M = ['M', 'A', 'G', 'C', 'T', 'I', 'V', 'L', 'K', 'H', 'F', 'Y', 'W']
F = ['F', 'A', 'G', 'C', 'T', 'I', 'V', 'L', 'K', 'H', 'Y', 'W', 'M']
P = ['P', 'N', 'D', 'T', 'C', 'A', 'G', 'S', 'V']
S = ['S', 'A', 'G', 'C', 'P', 'N', 'D', 'T', 'V', 'Q', 'E', 'K', 'R', 'H', 'Y', 'W']
T = ['T', 'A', 'G', 'C', 'I', 'V', 'L', 'K', 'H', 'F', 'Y', 'W', 'M', 'N', 'D', 'S', 'Q', 'E', 'R']
W = ['W', 'A', 'G', 'C', 'T', 'I', 'V', 'L', 'K', 'H', 'F', 'Y', 'M', 'N', 'Q', 'S', 'D', 'E', 'R']
Y = ['Y', 'A', 'G', 'C', 'T', 'I', 'V', 'L', 'K', 'H', 'F', 'W', 'M', 'N', 'Q', 'S', 'D', 'E', 'R']
V = ['V', 'A', 'G', 'C', 'T', 'I', 'L', 'K', 'H', 'F', 'Y', 'W', 'M', 'P', 'N', 'D', 'S']
gap = ['-']
#-----------------------------------------------------------------------
# Henikoff, S.; Henikoff, J.G.. (1992). "Amino Acid Substitution 
# Matrices from Protein Blocks". PNAS 89 (22) pp. 10915–10919. 
# DOI:10.1073/pnas.89.22.10915. PMID 1438297
BLOSUM62 = [[4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0, -2, -1, 0, -1],
            [-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3, -1, 0, -1, -1],
            [-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3, 3, 0, -1, -1],
            [-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3, 4, 1, -1, -1],
            [0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -1],
            [-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2, 0, 3, -1, -1],
            [-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -1],
            [0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3, -1, -2, -1, -1],
            [-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3, 0, 0, -1, -1],
            [-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3, -3, -3, -1, -1],
            [-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1, -4, -3, -1, -1],
            [-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2, 0, 1, -1, -1],
            [-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1, -3, -1, -1, -1],
            [-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1, -3, -3, -1, -1],
            [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2, -2, -1, -2, -1],
            [1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2, 0, 0, 0, -1],
            [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0, -1, -1, 0, -1],
            [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 1, 2, -3, -4, -3, -2, -1],
            [-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1, -3, -2, -1, -1],
            [0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4, -3, -2, -1, -1],
            [-2, -1, 3, 4, -3, 0, 1, -1, 0, -3, -4, 0, -3, -3, -2, 0, -1, -4, -3, -3, 4, 1, -1, -1],
            [-1, 0, 0, 1, -3, 3, 4, -2, 0, -3, -3, 1, -1, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -1],
            [0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, 0, 0, -2, -1, -1, -1, -1, -1, -1],
            [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]]

# Dayhoff, M.O., Schwartz, R. and Orcutt, B.C. (1978). "A model of 
# Evolutionary Change in Proteins". Atlas of protein sequence and 
# structure (volume 5, supplement 3 ed.). Nat. Biomed. Res. Found.
# pp. 345–358. ISBN 0-912466-07-3
PAM250 = [[2, -2, 0, 0, -2, 0, 0, 1, -1, -1, -2, -1, -1, -3, 1, 1, 1, -6, -3, 0, 0, 0, 0, -8],
          [-2, 6, 0, -1, -4, 1, -1, -3, 2, -2, -3, 3, 0, -4, 0, 0, -1, 2, -4, -2, -1, 0, -1, -8],
          [0, 0, 2, 2, -4, 1, 1, 0, 2, -2, -3, 1, -2, -3, 0, 1, 0, -4, -2, -2, 2, 1, 0, -8],
          [0, -1, 2, 4, -5, 2, 3, 1, 1, -2, -4, 0, -3, -6, -1, 0, 0, -7, -4, -2, 3, 3, -1, -8],
          [-2, -4, -4, -5, 12, -5, -5, -3, -3, -2, -6, -5, -5, -4, -3, 0, -2, -8, 0, -2, -4, -5, -3, -8],
          [0, 1, 1, 2, -5, 4, 2, -1, 3, -2, -2, 1, -1, -5, 0, -1, -1, -5, -4, -2, 1, 3, -1, -8],
          [0, -1, 1, 3, -5, 2, 4, 0, 1, -2, -3, 0, -2, -5, -1, 0, 0, -7, -4, -2, 3, 3, -1, -8],
          [1, -3, 0, 1, -3, -1, 0, 5, -2, -3, -4, -2, -3, -5, 0, 1, 0, -7, -5, -1, 0, 0, -1, -8],
          [-1, 2, 2, 1, -3, 3, 1, -2, 6, -2, -2, 0, -2, -2, 0, -1, -1, -3, 0, -2, 1, 2, -1, -8],
          [-1, -2, -2, -2, -2, -2, -2, -3, -2, 5, 2, -2, 2, 1, -2, -1, 0, -5, -1, 4, -2, -2, -1, -8],
          [-2, -3, -3, -4, -6, -2, -3, -4, -2, 2, 6, -3, 4, 2, -3, -3, -2, -2, -1, 2, -3, -3, -1, -8],
          [-1, 3, 1, 0, -5, 1, 0, -2, 0, -2, -3, 5, 0, -5, -1, 0, 0, -3, -4, -2, 1, 0, -1, -8],
          [-1, 0, -2, -3, -5, -1, -2, -3, -2, 2, 4, 0, 6, 0, -2, -2, -1, -4, -2, 2, -2, -2, -1, -8],
          [-3, -4, -3, -6, -4, -5, -5, -5, -2, 1, 2, -5, 0, 9, -5, -3, -3, 0, 7, -1, -4, -5, -2, -8],
          [1, 0, 0, -1, -3, 0, -1, 0, 0, -2, -3, -1, -2, -5, 6, 1, 0, -6, -5, -1, -1, 0, -1, -8],
          [1, 0, 1, 0, 0, -1, 0, 1, -1, -1, -3, 0, -2, -3, 1, 2, 1, -2, -3, -1, 0, 0, 0, -8],
          [1, -1, 0, 0, -2, -1, 0, 0, -1, 0, -2, 0, -1, -3, 0, 1, 3, -5, -3, 0, 0, -1, 0, -8],
          [-6, 2, -4, -7, -8, -5, -7, -7, -3, -5, -2, -3, -4, 0, -6, -2, -5, 17, 0, -6, -5, -6, -4, -8],
          [-3, -4, -2, -4, 0, -4, -4, -5, 0, -1, -1, -4, -2, 7, -5, -3, -3, 0, 10, -2, -3, -4, -2, -8],
          [0, -2, -2, -2, -2, -2, -2, -1, -2, 4, 2, -2, 2, -1, -1, -1, 0, -6, -2, 4, -2, -2, -1, -8],
          [0, -1, 2, 3, -4, 1, 3, 0, 1, -2, -3, 1, -2, -4, -1, 0, 0, -5, -3, -2, 3, 2, -1, -8],
          [0, 0, 1, 3, -5, 3, 3, 0, 2, -2, -3, 0, -2, -5, 0, 0, -1, -6, -4, -2, 2, 3, -1, -8],
          [0, -1, 0, -1, -3, -1, -1, -1, -1, -1, -1, -1, -1, -2, -1, 0, 0, -4, -2, -1, -1, -1, -1, -8],
          [-8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, 1]]

# McLachlan, A.D. Repeating sequences and gene duplication in 
# proteins. J. Mol. Biol. 64, 417-437 (1972)
MCLACHLAN = [[8.0, 2.0, 3.0, 3.0, 1.0, 3.0, 4.0, 3.0, 3.0, 2.0, 2.0, 3.0, 3.0, 1.0, 4.0, 4.0, 3.0, 1.0, 1.0, 3.0],
             [2.0, 8.0, 3.0, 1.0, 1.0, 5.0, 3.0, 3.0, 5.0, 1.0, 2.0, 5.0, 1.0, 1.0, 3.0, 4.0, 3.0, 3.0, 2.0, 2.0],
             [3.0, 3.0, 8.0, 5.0, 1.0, 4.0, 4.0, 3.0, 4.0, 1.0, 1.0, 4.0, 2.0, 0.0, 1.0, 5.0, 3.0, 0.0, 2.0, 1.0],
             [3.0, 1.0, 5.0, 8.0, 1.0, 4.0, 5.0, 3.0, 4.0, 0.0, 1.0, 3.0, 2.0, 1.0, 3.0, 3.0, 3.0, 0.0, 1.0, 1.0],
             [1.0, 1.0, 1.0, 1.0, 9.0, 0.0, 0.0, 1.0, 3.0, 1.0, 0.0, 0.0, 3.0, 0.0, 0.0, 2.0, 2.0, 2.0, 1.0, 1.0],
             [3.0, 5.0, 4.0, 4.0, 0.0, 8.0, 5.0, 2.0, 4.0, 0.0, 3.0, 4.0, 3.0, 0.0, 3.0, 4.0, 3.0, 2.0, 1.0, 2.0],
             [4.0, 3.0, 4.0, 5.0, 0.0, 5.0, 8.0, 3.0, 2.0, 1.0, 1.0, 4.0, 1.0, 0.0, 4.0, 4.0, 4.0, 1.0, 2.0, 2.0],
             [3.0, 3.0, 3.0, 3.0, 1.0, 2.0, 3.0, 8.0, 2.0, 1.0, 1.0, 3.0, 1.0, 0.0, 3.0, 3.0, 2.0, 1.0, 0.0, 2.0],
             [3.0, 5.0, 4.0, 4.0, 3.0, 4.0, 2.0, 2.0, 8.0, 2.0, 2.0, 4.0, 3.0, 4.0, 3.0, 3.0, 4.0, 3.0, 4.0, 2.0],
             [2.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 2.0, 8.0, 5.0, 1.0, 5.0, 3.0, 1.0, 2.0, 3.0, 3.0, 3.0, 5.0],
             [2.0, 2.0, 1.0, 1.0, 0.0, 3.0, 1.0, 1.0, 2.0, 5.0, 8.0, 2.0, 6.0, 5.0, 1.0, 2.0, 3.0, 3.0, 3.0, 5.0],
             [3.0, 5.0, 4.0, 3.0, 0.0, 4.0, 4.0, 3.0, 4.0, 1.0, 2.0, 8.0, 1.0, 0.0, 3.0, 3.0, 3.0, 1.0, 1.0, 2.0],
             [3.0, 1.0, 2.0, 2.0, 3.0, 3.0, 1.0, 1.0, 3.0, 5.0, 6.0, 1.0, 8.0, 5.0, 1.0, 2.0, 3.0, 1.0, 2.0, 4.0],
             [1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 4.0, 3.0, 5.0, 0.0, 5.0, 9.0, 1.0, 2.0, 1.0, 6.0, 6.0, 3.0],
             [4.0, 3.0, 1.0, 3.0, 0.0, 3.0, 4.0, 3.0, 3.0, 1.0, 1.0, 3.0, 1.0, 1.0, 8.0, 3.0, 3.0, 0.0, 0.0, 2.0],
             [4.0, 4.0, 5.0, 3.0, 2.0, 4.0, 4.0, 3.0, 3.0, 2.0, 2.0, 3.0, 2.0, 2.0, 3.0, 8.0, 5.0, 3.0, 3.0, 2.0],
             [3.0, 3.0, 3.0, 3.0, 2.0, 3.0, 4.0, 2.0, 4.0, 3.0, 3.0, 3.0, 3.0, 1.0, 3.0, 5.0, 8.0, 2.0, 1.0, 3.0],
             [1.0, 3.0, 0.0, 0.0, 2.0, 2.0, 1.0, 1.0, 3.0, 3.0, 3.0, 1.0, 1.0, 6.0, 0.0, 3.0, 2.0, 9.0, 6.0, 2.0],
             [1.0, 2.0, 2.0, 1.0, 1.0, 1.0, 2.0, 0.0, 4.0, 3.0, 3.0, 1.0, 2.0, 6.0, 0.0, 3.0, 1.0, 6.0, 9.0, 3.0],
             [3.0, 2.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 5.0, 5.0, 2.0, 4.0, 3.0, 2.0, 2.0, 3.0, 2.0, 3.0, 8.0]]


# Glaser F, Steinberg DM, Vakser IA, Ben-Tal N. Residue frequencies
# and pairing preferences at protein-protein interfaces. Proteins. 2001
# May 1;43(2):89-102.
CPVN = [[3.89, 4.91, 4.59, 5.33, 1.76, 5.25, 2.84, 0.77, 3.05, 1.00, 6.24, 5.61, 3.27, 3.38, 3.20, 3.60, 2.30, 1.59, 3.23, 3.80],
        [4.91, 3.74, 4.20, 4.69, 2.89, 4.37, 2.57, -0.41, 2.83, 1.42, 2.92, 3.95, 2.90, 3.21, 3.22, 3.22, 1.93, 1.36, 4.45, 4.18],
        [4.59, 4.20, 4.03, 4.86, 2.93, 5.32, 2.77, -0.37, 2.07, 1.41, 5.77, 4.19, 2.50, 4.88, 3.12, 3.46, 1.40, 2.31, 3.15, 4.99],
        [5.33, 4.69, 4.86, 5.34, 3.68, 5.28, 3.00, 0.14, 3.34, 1.75, 5.83, 5.83, 4.25, 3.47, 2.87, 4.25, 0.99, 3.11, 3.57, 4.49],
        [1.76, 2.89, 2.93, 3.68, 7.65, 1.84, 1.46, -0.25, 1.03, 2.48, 2.14, 2.47, 2.74, 4.12, 2.51, 1.33, 0.24, -0.42, 2.05, 2.81],
        [5.25, 4.37, 5.32, 5.28, 1.84, 6.02, 2.30, 0.91, 2.09, 1.61, 4.89, 4.81, 3.38, 4.65, 3.88, 4.18, 0.36, 2.30, 3.93, 3.62],
        [2.84, 2.57, 2.77, 3.00, 1.46, 2.30, -0.52, -1.77, 1.21, 0.39, 3.37, 2.47, 1.22, 2.59, 1.71, 1.72, 1.13, 1.69, 2.13, 1.90],
        [0.77, -0.41, -0.37, 0.14, -0.25, 0.91, -1.77, -4.40, 0.21, -1.53, 1.42, 1.25, -0.51, 1.08, -0.89, 0.70, -0.08, -0.54, 1.33, 1.59],
        [3.05, 2.83, 2.07, 3.34, 1.03, 2.09, 1.21, 0.21, 1.27, 1.91, 5.12, 3.14, 2.65, 2.71, 2.88, 1.82, 3.88, 2.52, 3.67, 3.77],
        [1.00, 1.42, 1.41, 1.75, 2.48, 1.61, 0.39, -1.53, 1.91, -0.09, 2.87, 2.30, 1.33, 0.80, 2.60, 2.00, 2.94, 1.77, 2.74, 2.82],
        [6.24, 2.92, 5.77, 5.83, 2.14, 4.89, 3.37, 1.42, 5.12, 2.87, 5.85, 6.19, 7.87, 6.46, 1.20, 1.37, 2.62, 3.54, 5.76, 8.57],
        [5.61, 3.95, 4.19, 5.83, 2.47, 4.81, 2.47, 1.25, 3.14, 2.30, 6.19, 5.93, 4.22, 6.05, 4.54, 2.05, 1.76, 3.66, 5.26, 5.28],
        [3.27, 2.90, 2.50, 4.25, 2.74, 3.38, 1.22, -0.51, 2.65, 1.33, 7.87, 4.22, 0.60, 2.89, 3.17, 3.50, 1.46, 3.09, 3.75, 3.99],
        [3.38, 3.21, 4.88, 3.47, 4.12, 4.65, 2.59, 1.08, 2.71, 0.80, 6.46, 6.05, 2.89, 5.37, 2.30, 4.00, 5.20, 2.38, 2.72, 4.90],
        [3.20, 3.22, 3.12, 2.87, 2.51, 3.88, 1.71, -0.89, 2.88, 2.60, 1.20, 4.54, 3.17, 2.30, 1.65, 1.95, 0.08, 2.68, 5.32, 5.75],
        [3.60, 3.22, 3.46, 4.25, 1.33, 4.18, 1.72, 0.70, 1.82, 2.00, 1.37, 2.05, 3.50, 4.00, 1.95, 2.83, 3.26, 3.45, 3.50, 4.50],
        [2.30, 1.93, 1.40, 0.99, 0.24, 0.36, 1.13, -0.08, 3.88, 2.94, 2.62, 1.76, 1.46, 5.20, 0.08, 3.26, 0.13, 3.85, 3.90, 4.94],
        [1.59, 1.36, 2.31, 3.11, -0.42, 2.30, 1.69, -0.54, 2.52, 1.77, 3.54, 3.66, 3.09, 2.38, 2.68, 3.45, 3.85, 2.92, 3.17, 3.85],
        [3.23, 4.45, 3.15, 3.57, 2.05, 3.93, 2.13, 1.33, 3.67, 2.74, 5.76, 5.26, 3.75, 2.72, 5.32, 3.50, 3.90, 3.17, 3.24, 2.29],
        [3.80, 4.18, 4.99, 4.49, 2.81, 3.62, 1.90, 1.59, 3.77, 2.82, 8.57, 5.28, 3.99, 4.90, 5.75, 4.50, 4.94, 3.85, 2.29, 2.87]]

# Singer MS, Vriend G, Bywater RP. Prediction of protein residue 
# contacts with a PDB-derived likelihood matrix. Protein Eng.
# 2002 Sep;15(9):721-5.
CLM = [[0.50, 0.90, 0.30, 0.30, 1.40, 0.30, 0.60, 1.30, 0.30, 1.30, 1.10, 0.40, 0.50, 0.50, 0.50, 0.40, 0.60, 1.10, 1.40, 1.30],
       [0.90, 9.60, 0.60, 0.50, 2.60, 0.70, 1.80, 1.80, 0.60, 1.70, 1.50, 0.70, 1.00, 0.80, 0.80, 0.80, 0.80, 1.70, 2.70, 2.20],
       [0.30, 0.60, 0.50, 0.30, 0.50, 0.40, 1.10, 0.50, 1.00, 0.40, 0.50, 0.70, 0.40, 0.50, 1.50, 0.60, 0.50, 0.40, 1.10, 1.20],
       [0.30, 0.40, 0.40, 0.40, 0.80, 0.30, 1.00, 0.60, 1.10, 0.60, 0.70, 0.60, 0.50, 0.60, 1.50, 0.50, 0.60, 0.50, 1.20, 1.30],
       [1.40, 2.40, 0.60, 0.70, 4.30, 0.90, 1.40, 3.30, 0.70, 3.20, 2.80, 0.80, 1.30, 1.10, 1.40, 0.90, 1.10, 2.60, 4.10, 3.10],
       [0.40, 0.60, 0.40, 0.30, 1.00, 0.40, 0.60, 0.60, 0.50, 0.60, 0.70, 0.60, 0.40, 0.50, 0.70, 0.40, 0.50, 0.60, 1.10, 1.00],
       [0.70, 1.50, 1.40, 1.10, 1.20, 0.70, 2.20, 1.10, 0.70, 1.20, 1.60, 0.80, 0.60, 0.90, 1.10, 1.00, 1.00, 0.90, 2.30, 1.90],
       [1.40, 1.70, 0.50, 0.60, 3.40, 0.60, 1.10, 3.50, 0.60, 3.20, 2.60, 0.60, 0.90, 0.80, 1.00, 0.70, 1.10, 2.70, 3.20, 2.80],
       [0.30, 0.60, 1.00, 1.10, 0.80, 0.40, 0.60, 0.60, 0.30, 0.60, 0.60, 0.60, 0.40, 0.70, 0.50, 0.50, 0.50, 0.50, 0.90, 1.10],
       [1.20, 1.60, 0.40, 0.60, 3.70, 0.60, 1.10, 3.30, 0.60, 3.30, 2.70, 0.60, 0.90, 0.90, 1.00, 0.60, 1.00, 2.60, 3.40, 2.50],
       [1.10, 1.60, 0.50, 0.50, 3.40, 0.60, 1.40, 2.40, 0.50, 2.50, 2.90, 0.80, 1.00, 1.00, 1.00, 0.70, 0.80, 2.20, 3.60, 2.80],
       [0.40, 0.70, 0.70, 0.60, 0.80, 0.60, 0.90, 0.60, 0.60, 0.60, 0.90, 0.90, 0.70, 0.90, 1.00, 0.60, 0.80, 0.60, 1.30, 1.20],
       [0.50, 1.00, 0.50, 0.50, 1.40, 0.30, 1.00, 0.80, 0.40, 0.90, 1.10, 0.70, 0.70, 0.70, 0.90, 0.50, 0.60, 0.80, 2.20, 1.70],
       [0.40, 0.80, 0.60, 0.50, 1.10, 0.60, 1.00, 0.80, 0.70, 0.80, 0.90, 0.80, 0.70, 0.80, 1.60, 0.50, 0.80, 0.70, 1.50, 1.20],
       [0.50, 0.90, 1.70, 1.60, 1.50, 0.70, 1.40, 1.10, 0.50, 1.10, 1.10, 1.00, 0.90, 1.00, 1.00, 0.90, 0.90, 1.00, 1.90, 1.90],
       [0.30, 0.60, 0.50, 0.50, 0.80, 0.30, 0.90, 0.60, 0.50, 0.60, 0.60, 0.60, 0.40, 0.70, 0.70, 0.40, 0.60, 0.60, 1.20, 0.90],
       [0.50, 0.90, 0.60, 0.70, 1.10, 0.40, 1.00, 1.10, 0.60, 1.00, 1.00, 0.60, 0.60, 0.80, 0.90, 0.70, 0.80, 0.90, 1.40, 1.20],
       [1.10, 1.40, 0.40, 0.60, 3.00, 0.50, 1.00, 2.60, 0.60, 2.80, 2.20, 0.60, 0.80, 0.80, 0.80, 0.70, 0.90, 2.40, 2.30, 2.00],
       [1.30, 2.40, 1.20, 1.00, 3.60, 1.10, 2.30, 3.00, 1.20, 3.70, 3.70, 1.00, 1.90, 1.60, 1.80, 1.10, 1.10, 2.70, 4.20, 3.40],
       [1.20, 1.80, 1.10, 1.20, 3.30, 0.90, 1.90, 2.20, 1.20, 2.20, 2.60, 1.00, 1.50, 1.30, 1.60, 1.00, 1.10, 2.00, 3.00, 2.50]]


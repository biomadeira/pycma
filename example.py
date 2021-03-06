#!/usr/bin/env python
# -*- coding: UTF-8 -*-                                                         
# Copyright (c) 2012 F. Madeira and L. Krippahl, 2012                                            
# This code is part of PyCMA distribution.                                 
# GNU General Public License - See LICENSE for more details.

#-----------------------------------------------------------------------

"""
Testing the installation of PyCMA.
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

def main():
    "Example for testing purposes."
    
    import pycma
    
    # load alignments
    p = pycma.Alignment('example1.fasta','example2.fasta','fasta')
    p.number_sequences()
    p.length_sequences()
    
    # converting between alignment formats
    p.write_alignment(p.alignment1, 'test1.aln', form='clustal')
    p.write_alignment(p.alignment2, 'test2.aln', form='clustal')
    
    # pre-processing tasks
    pre = pycma.PreProcessing(p)
    
    # consider gaps in the query sequence (1st sequence)
    pre.gaps_in_query(False)
    
    # print alphabet reduction dictionaries
    pre.alphabets()
    
    # alphabet reduction
    pre.alphabet_reduction('murphy_10')
    
    # print alignment formats
    pre.formats()
    
    #re-write processed alignments in ClustalW format
    pre.write_alignment(pre.alignment1, 'test3.aln', form='clustal')
    pre.write_alignment(pre.alignment2, 'test4.aln', form='clustal')
    
    # re-print alignment info
    pre.number_sequences()
    pre.length_sequences()
    
    # coevolution object
    c1 = pycma.Coevolution(pre)
    
    # print methods
    c1.methods()
    
    # print Mutual Information normalisations
    c1.mi_normalizations()
    
    # print matrices for Correlation Coefficient methods 
    c1.cc_matrices()
    
    # Mutual Information without normalisation
    c1.MI('MI')
    
    # results  
    r = pycma.Results(c1)
    
    # normalizing scores [~0:1]
    r.normalization(True)
    
    # 20 best scores
    r.scores(20)
    
    # writing histograms, heatmaps, and coevolution scores
    r.histogram('test1.png')
    r.heatmap('test2.png')
    r.write_file('example1.txt')
    
    # different format
    r.write_file('example2.txt', scores=20, form='plain')
    
    # another example using Statistical Coupling Analysis
    c2 = pycma.Coevolution(pre)
    c2.SCA()
        
    r = pycma.Results(c2)
    
    # sensing neighbour residues influence
    r.grouping(3)

    r.normalization(True)
    r.scores(50)
    r.histogram('test3.png')
    r.heatmap('test4.png')
    r.write_file('example3.txt')
    r.write_file('example4.txt', scores=50, form='basic')
    
    print "Analysis complete!"
    print "Check the output files on your working directory!"
    
    
if __name__ == "__main__":
    main()

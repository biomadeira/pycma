#!/usr/bin/env python
# -*- coding: UTF-8 -*-                                                      
# Copyright (c) 2012 F. Madeira and L. Krippahl, 2012                                            
# This code is part of PyCMA distribution.                                 
# GNU General Public License - See LICENSE for more details.

#TODO:
# - Correlation Measures (CoMap, MirrorTree, CAPS, LnLCorr)
# - Conservation Class and Measures (Valdar, 2002)
# - Interaction Map - NetworkX, Graphviz?
# - Statistical analysis of results: ROC, PCA, Stats Significance
# - Docstring (doctest, unittest, etc.) and _test()
#-----------------------------------------------------------------------

"""
PyCMA - A Python Module for Correlated Mutation Analysis.

    PyCMA is a module for python by Fabio Madeira & Ludwig Krippahl.
    See https://github.com/fmadeira/pycma/ for documentation and the 
    latest versions.
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

"""
 example usage:
    >>> import pycma
    
    # load alignments
    >>> p = pycma.Alignment('example1.fasta','example2.fasta','fasta')
    >>> p.number_sequences()
    >>> p.length_sequences()
    
    # converting between alignment formats
    >>> p.write_alignment(p.alignment1, 'test1.aln', form='clustal')
    >>> p.write_alignment(p.alignment2, 'test2.aln', form='clustal')
    
    # pre-processing tasks
    >>> pre = pycma.PreProcessing(p)
    
    # consider gaps in the query sequence (1st sequence)?
    >>> pre.gaps_in_query(False)
    
    # print alphabet reduction dictionaries
    >>> pre.alphabets()
    
    # alphabet reduction
    >>> pre.alphabet_reduction('murphy_10')
    
    # print alignment formats
    >>> pre.formats()
    
    #re-write processed alignments in ClustalW format
    >>> pre.write_alignment(pre.alignment1, 'test3.aln', form='clustal')
    >>> pre.write_alignment(pre.alignment2, 'test4.aln', form='clustal')
    
    # re-print alignment info
    >>> pre.number_sequences()
    >>> pre.length_sequences()
    
    # coevolution object
    >>> c1 = pycma.Coevolution(pre)
    
    # print methods
    >>> c1.methods()
    
    # print Mutual Information normalizations
    >>> c1.mi_normalizations()
    
    # print matrices for Correlation Coefficient methods 
    >>> c1.cc_matrices()
    
    # Mutual Information without normalization
    >>> c1.MI('MI')
    
    # results  
    >>> r = pycma.Results(c1)
    
    # normalizing scores [~0:1]
    >>> r.normalization(True)
    
    # 20 best scores
    >>> r.scores(20)
    
    # writing histograms, heatmaps, and coevolution scores
    >>> r.histogram('test1.png')
    >>> r.heatmap('test2.png')
    >>> r.write_file('example1.txt')
    
    # different format
    >>> r.write_file('example2.txt', scores=20, form='plain')
    
    # another example using Statistical Coupling Analysis
    >>> c2 = pycma.Coevolution(pre)
    >>> c2.SCA()
        
    >>> r = pycma.Results(c2)
    
    # sensing neighbour residues influence
    >>> r.grouping(3)

    >>> r.normalization(True)
    >>> r.scores(50)
    >>> r.histogram('test3.png')
    >>> r.heatmap('test4.png')
    >>> r.write_file('example3.txt')
    >>> r.write_file('example4.txt', scores=50, form='basic')
"""

__version__ = "0.1.0"
  
#-----------------------------------------------------------------------
import os
import sys
import tempfile
from math import factorial, log, e
from numpy import mean, std, sqrt, zeros
from sets import Set
from matplotlib import pyplot
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
from utils import *
#-----------------------------------------------------------------------
 
class Alignment(object):
    """
    pycma.Alignment object
     
    Inputs the alignment file (MSA) in one of the available formats
    (see Alignformats). The user can access the same functionallity 
    obtained in Biopython (AlignIO).
    
    The default use of pycma needs the input of two 
    prepared alignments. For inter-protein coevolution analysis 
    both alignments must have the same number of sequences (sorted
    by taxa or organism). 
    If the user inputs only one MSA, intra-protein coevolution
    analysis will be performed. In this case the score of one
    position with itself will be automatically ignored.
    
    Attributes: id1, id2, file1, file2, form, alignment1, alignment2,
                sequences1, sequences2, info1, info2
              
    Methods: number_sequences, length_sequences, write_alignment, 
             formats
    
    Usage:
    >>> import pycma
    >>> p = pycma.Alignment('example1.fasta','example2.fasta','fasta')
    >>> p.number_sequences()
    >>> p.length_sequences()
    >>> p.formats()
    >>> p.write_alignment(p.alignment1, 'test1.aln', form='clustal')
    >>> p.write_alignment(p.alignment2, 'test2.aln', form='clustal')
    """     
    def __init__(self, file1=None, file2=None, form=None):
        if file1 is not None and not os.path.exists(file1): 
            raise IOError('Input file1 does not exist!')
            sys.exit()
        elif file2 is not None and not os.path.exists(file2): 
            raise IOError('Input file2 does not exist!')
            sys.exit()       
        elif file1 is not None and file2 is not None:
            # attributes
            self.id1 = file1.split(".")[0]
            self.id2 = file2.split(".")[0]
            self.file1 = file1
            self.file2 = file2
            self.form = Alignformat(form).form
            self.alignment1 = Alignment.ReadAlignment(
                                self.file1, self.form)
            self.alignment2 = Alignment.ReadAlignment(
                                self.file2, self.form)
            self.sequences1 = Alignment.ParseSequences(
                                        self.alignment1)
            self.sequences2 = Alignment.ParseSequences(
                                        self.alignment2)
            self.info1 = Alignment.ParseInfo(
                                self.alignment1)
            self.info2 = Alignment.ParseInfo(
                                self.alignment2)
            self.model = "inter"
            
        elif file1 is not None and file2 is None:
            # attributes
            self.id1 = file1.split(".")[0]
            self.id2 = self.id1
            self.file1 = file1
            self.file2 = self.file1
            self.form = Alignformat(form).form
            self.alignment1 = Alignment.ReadAlignment(
                                self.file1, self.form)
            self.alignment2 = self.alignment1
            self.sequences1 = Alignment.ParseSequences(
                                        self.alignment1)
            self.sequences2 = self.sequences1
            self.info1 = Alignment.ParseInfo(
                                self.alignment1)
            self.info2 = self.info1
            self.model = "intra"
        else:
            self.alignment1 = file1
            self.alignment2 = file2
            print self.__usage__()
            
    def __name__(self):
        return "Alignment Object"
    
    def __len__(self):
        return len(self.alignment1), len(self.alignment2)
    
    def __str__(self):
        return str(self.alignment1), str(self.alignment2)
    
    def __getitem__(self, index):
        return self.alignment1[index], self.alignment2[index]
    
    def __getslice__(self, start, end):
        return self.alignment1[start, end], self.alignment2[start, end]
    
    def __iter__(self):
        return iter(self.alignment1), iter(self.alignment2)
    
    def __usage__(self):
        return """usage pycma.Alignment(file1="align1.format", form="format")
        See https://github.com/fmadeira/pycma/ for documentation."""
    
    def __id__(self):
        return self.id1, self.id2
    
    def __file__(self):
        return self.file1, self.file2
    
    def formats(self):
        return Alignformat(self.form).formats()
    
    @staticmethod
    def ReadAlignment(inp, form):
        """
        Opens the alignment file in memory. Uses the Biopython
        Bio.AlignIO.read() function. Each sequence as a SeqRecord 
        is store in a dictionary. It uses an extended Protein 
        Alphabet including gaps.
        """

        alignment = AlignIO.read(inp, form,
                            alphabet=Gapped(IUPAC.extended_protein))
        # minimum number of sequences
        if len(alignment) < 15:
            raise NumberSequencesError()
        
        return alignment
    
    @staticmethod
    def ParseSequences(alignment):
        """
        Parses sequences and stores them into a list. 
        The sequences keep the SeqRecord attributes.
        """
        sequences = []
        for record in alignment:
            seq = record.seq
            sequences.append(seq)
            
        return sequences
    
    @staticmethod
    def ParseInfo(alignment):
        """
        Parses sequence descriptions and stores them in a list.
        The sequence descriptions keep the SeqRecord attributes.
        """ 
        info = []  
        for record in alignment:
            description = record.id
            info.append(description)
            
        return info
    
    @staticmethod
    def writeToAlignment(alignment, filename, form):
        """
        Writes processed alignments in one of the available formats.
        """
        
        return AlignIO.write(alignment, filename, form)
    
    def number_sequences(self):
        """
        Prints the number of sequences in the inputed alignments.
        """  
        
        return "Alignments have %i sequences." % (len(self.sequences1))
        
    def length_sequences(self):
        """
        Prints the length of the sequences in the inputed alignments.
        """
        
        return "Sequences1 have %i residues and Sequences2 have %i residues." \
            % (len(self.sequences1[0]), len(self.sequences2[0]))
   
    def write_alignment(self, alignment, filename, form=None):
        """
        Writes processed alignments in one of the available formats.
        """
        write_form = Alignformat(form).form
        
        return Alignment.writeToAlignment(alignment, filename,
                                          write_form)


_align_format_builtins = ['clustal', 'fasta', 'nexus',
                           'phylip', 'stockholm']
 
class Alignformat(object):
    """
    pycma.Alignformat object (extends the Alignment object)
        
    formats: 'clustal', 'fasta', 'nexus', 'phylip', 'stockholm'
        default = 'fasta'
        
    Attributes: form
        
    Methods: default, formats
    """
    def __init__(self, form=None):
        if form is not None and form not in _align_format_builtins:
            raise AlignformatError()
        elif form is not None and form in _align_format_builtins:
            self.form = form
        else:
            return self.__default__()
    
    def __call__(self):
        return self.form
    
    def __name__(self):
        return "Alignformat Object"
    
    def __str__(self):
        return str(self.form)
     
    def __repr__(self):
        return repr(self.form)
        
    def __default__(self):
        self.form = "fasta"
        return self.form
    
    def default(self):
        return self.__default__()
    
    def formats(self):
        return """Formats: 'clustal', 'fasta', 'nexus', 'phylip' or 'stockholm'."""

class PreProcessing(Alignment):
    """
    pycma.PreProcessing object
        
    Preprocessing functions to handle gaps and alphabets. 
        
    gaps_per_sequence=100 -> defines the ratio of gaps in each sequence.
    Routine to preprocess the ratio of gaps per sequence. It allows to 
    specify a minimum number of sequences in the alignment, prior to 
    the analysis. This routine updates self.sequences and self.info.
    It ranges from [0:100], i.e. no gaps in a sequence to gaps in all 
    positions of a sequence. For example gaps_per_sequence=50, removes 
    sequences (in both alignments) that have gaps in more than half 
    the sequence length. 
    
    number_sequences=15 -> ensures a minimum of 15 sequences after 
    processing by the previous function.
    
    gaps_in_query=true -> if False positions (columns) with gaps in the 
    site 1 (first sequence) are remove from the alignment. This makes  
    sence when your first sequence is the sequence of interest, and 
    it is not supposed to have gaps (e.g. when you are aligning your 
    sequences with a protein structure).
    
    alphabet=None -> select from the available alphabets or import a 
    custom one as a dictionary. Example {'A','D','G','-'}:
    test = {'A':'A', 'C':'A', 'D':'D',
          'E':'D', 'F':'A', 'G':'G',
          'H':'D', 'K':'D', 'I':'A',
          'L':'A', 'M':'A', 'N':'D',
          'P':'A', 'Q':'D', 'R':'D',
          'S':'G', 'T':'G', 'V':'A',
          'Y':'G', 'W':'G', '-':'-'}
          
    where 'A', 'D', 'G' and '-' correspond to the possible groups of 
    amino acids. Each residue will be updated accordingly. 
        
    Attributes: id1, id2, file1, file2, form, alignment1, alignment2,
                sequences1, sequences2, info1, info2, _gaps_in_query, 
                _gaps_per_sequence, _number_sequences, alpha, 
                model
                
    Methods: gaps_in_query, gaps_per_sequence, alphabet_reduction,
             alphabets
    
    Usage:
    >>> import pycma
    >>> p = pycma.Alignment('example1.fasta','example2.fasta','fasta')
    >>> pre = pycma.PreProcessing(p)
    >>> pre.gaps_in_query(False)
    >>> pre.alphabets()
    >>> pre.alphabet_reduction('murphy_10')
    >>> pre.formats()
    >>> pre.write_alignment(pre.alignment1, 'test3.aln', form='clustal')
    >>> pre.write_alignment(pre.alignment2, 'test4.aln', form='clustal')
    >>> pre.number_sequences()
    >>> pre.length_sequences()
    """
    def __init__(self, alignments):
        try:
            assert len(alignments.alignment1) == len(alignments.alignment2)
        except:
            raise SameNumberSequencesError()

        # attributes
        self.alignment1 = alignments.alignment1
        self.alignment2 = alignments.alignment2
        self.sequences1 = alignments.sequences1
        self.sequences2 = alignments.sequences2
        self.id1 = alignments.id1
        self.id2 = alignments.id2
        self.info1 = alignments.info1
        self.info2 = alignments.info2
        self.file1 = alignments.file1
        self.file2 = alignments.file2
        self.form = alignments.form
        self.model = alignments.model
        self._gaps_in_query = True
        self._gaps_per_sequence = 100
        self._number_sequences = 15
        self.alpha = None
    
    def __name__(self):
        return "PreProcessing Object"
    
    def __str__(self):
        return str(self.alignment1), str(self.alignment2)
    
    def __repr__(self):
        return repr(self.alignment1), repr(self.alignment2)
    
    def alphabets(self):
        return Alphabet(self.alpha).alphabets()
    
    @staticmethod
    def cutMSA(alignment_to_cut, gaps_in_query=True, alpha=None):
        """
        Removes columns that have gaps in the first sequence.
        If alphabet!=None, it converts residues in accordance
        with the alphabet specified.
        """
               
        align = []
        columns = []
        positions = []
        blocks = []
        new_align = []
        new_align_ord = []
        new_align_concate = []
        cut_alignment = []
        aa_red = alpha
        
        alignment = alignment_to_cut
        
        align_length = alignment.get_alignment_length()
        for position in range(0, align_length):
            column = alignment[:, position]
            align.append(column)
            if gaps_in_query == False:
                if column[0] != "-":
                    columns.append(column)
                    positions.append(position)
            else:
                columns.append(column)
                positions.append(position)     
                
        for i in range(0, len(positions), 1):
            beg = int(positions[i])
            end = int(positions[i] + 1)
            block = alignment[:, beg:end]         
            blocks.append(block)
        
        for block in blocks:
            for record in block:
                seq = str(record.seq)
                new_align.append(seq)
        
        numb_blocks = len(new_align) / len(columns[0])
        for i in range(0, len(columns[0])):
            for j in range(0, len(new_align), len(columns[0])):
                new_align_ord.append(new_align[i + j])
    
        for i in range(0, len(new_align_ord), numb_blocks):
            pseudolist = new_align_ord[i:i + numb_blocks]
            lst = ""
            for j in pseudolist:
                lst += j
            new_align_concate.append(lst)
    
        for seq in new_align_concate:
            if aa_red != None:
                red = [Alphabet.alphabetReduct(e, aa_red) for e in seq]
                new_seq = ""
                for i in red:
                    new_seq += str(i)
                cut_alignment.append(new_seq)
            else:
                cut_alignment.append(seq)      
           
        return cut_alignment    
    
    @staticmethod
    def writeTmpAlignment(sequence, info):
        "Writes temporary alignments in FASTA form"
        
        tmp = tempfile.mkstemp()
        out = open(tmp[1], 'w')
        for seq, inf in zip(sequence, info):
            print >> out, ">%s\n%s" % (inf, seq)
        out.close()
        
        return tmp[1]
    
    def updateAlignments(self, sequence1, info1, sequence2, info2,
                         alignment1, alignment2):
        "Temporary alignments to update self.alignments"
        
        tmp_align1 = PreProcessing.writeTmpAlignment(
                        self.sequences1, self.info1)
        tmp_align2 = PreProcessing.writeTmpAlignment(
                        self.sequences2, self.info2)
        
        self.alignment1 = AlignIO.read(tmp_align1, "fasta",
                                alphabet=Gapped(IUPAC.extended_protein))
        self.alignment2 = AlignIO.read(tmp_align2, "fasta",
                                alphabet=Gapped(IUPAC.extended_protein))
        
    
    def gaps_in_query(self, inp):
        """
        Defines if the first sequence has gaps in it, as if it was 
        the protein sequence in stuy. If False columns with '-' in
        first site will be removed.
        """
        if inp is True or inp is False:
                self._gaps_in_query = inp
        else:
            raise BooleanError()      
        
        # gaps in query
        align1 = PreProcessing.cutMSA(self.alignment1, self._gaps_in_query)
        align2 = PreProcessing.cutMSA(self.alignment2, self._gaps_in_query)
        
        self.sequences1 = []
        for seq in align1:
            self.sequences1.append(seq)
        
        self.sequences2 = []
        for seq in align2:
            self.sequences2.append(seq)
        
        self.updateAlignments(self.sequences1, self.info1,
                              self.sequences2, self.info2,
                              self.alignment1, self.alignment2)
        
        return 

    def gaps_per_sequence(self, inp=100, sequences=15):
        """
        Defines the ratio of gaps per sequence. Optionally, it
        defines the minimum number of sequences, after processing.
        """
        try:
            inp = int(inp)
            sequences = int(inp)
        except:
            raise IntError()
        
        self._gaps_per_sequence = inp
        self._number_sequences = sequences
        
        # updating gaps per sequence
        pos = []
        index = -1
        length = len(self.sequences1[0])
        for seq in self.sequences1:
            count = seq.count('-')
            percent = int(count * 100.0 / length)
            index += 1
            if percent > int(self._gaps_per_sequence):
                pos.append(index)
        
        index = -1
        length = len(self.sequences2[0]) 
        for seq in self.sequences2:
            count = seq.count('-')
            percent = int(count * 100.0 / length)
            index += 1
            if percent > int(self._gaps_per_sequence):
                if index not in pos:
                    pos.append(index)
        count = 0
        for p in pos:
            self.sequences1.pop(p - count)
            self.sequences2.pop(p - count)
            self.info1.pop(p - count)
            self.info2.pop(p - count)
            count += 1 
        
        # checking the number of sequences
        nseq = self._number_sequences
        if len(self.sequences1) < nseq:
            raise NumberSequencesError(nseq)
        elif len(self.sequences2) < nseq:
            raise NumberSequencesError(nseq)
        else: 
            pass
        
        self.updateAlignments(self.sequences1, self.info1,
                              self.sequences2, self.info2,
                              self.alignment1, self.alignment2)
        return
    
    def alphabet_reduction(self, inp=None):
        """
        Defines an alphabet for alphabet reduction. User defined
        alphabets can be inputed as a python dictionary. See
        Preprocessing class docstring for more detais.
        """
        self.alpha = Alphabet(inp).alpha
        
        # gaps in query = True
        align1 = PreProcessing.cutMSA(self.alignment1,
                                    alpha=self.alpha)
        align2 = PreProcessing.cutMSA(self.alignment2,
                                    alpha=self.alpha)
        
        self.sequences1 = []
        for seq in align1:
            self.sequences1.append(seq)
        
        self.sequences2 = []
        for seq in align2:
            self.sequences2.append(seq)
        
        self.updateAlignments(self.sequences1, self.info1,
                              self.sequences2, self.info2,
                              self.alignment1, self.alignment2)
        return 

_alphabet_builtins = ['chemical_2_1', 'chemical_2_2', 'chemical_2_3',
                      'chemical_2_4', 'chemical_3_1', 'chemical_3_2',
                      'chemical_3_3', 'chemical_4', 'chemical_5',
                      'chemical_6', 'murphy_2', 'murphy_4', 'murphy_8',
                      'murphy_10', 'murphy_15', 'wang_2', 'wang_3',
                      'wang_5_1', 'wang_5_2', 'li_3', 'li_4', 'li_5',
                      'li_10']

class Alphabet(object):
    """
    pycma.Alphabet object
        
    Methods and attributes of the class Alphabet. Defines the 
    architecture of the alphabets for alphabet reduction operations.
    
    Alphabets: 'chemical_2_1', 'chemical_2_2', 'chemical_2_3', 
        'chemical_2_4', 'chemical_3_1', 'chemical_3_2', 'chemical_3_3', 
        'chemical_4', 'chemical_5', 'chemical_6', 'murphy_2', 'murphy_4', 
        'murphy_8', 'murphy_10', 'murphy_15', 'wang_2', 'wang_3', 'wang_5_1', 
        'wang_5_2', 'li_3', 'li_4', 'li_5' or 'li_10'
    
    Attributes: alpha
    
    Methods: default, alphabets
    """
    def __init__(self, alphabet=None):
        if alphabet is not None and alphabet not in _alphabet_builtins \
        and type(alphabet) != dict:
            raise AlphabetError()
        elif alphabet is not None and alphabet not in _alphabet_builtins \
        and type(alphabet) == dict:
            global cAlphabet
            cAlphabet = alphabet
            self.alpha = "custom"
        elif alphabet is not None and alphabet in _alphabet_builtins \
        and type(alphabet) != dict:
            self.alpha = alphabet
        else:
            return self.__default__()
    
    def __name__(self):
        return "Alphabet Object"
    
    def __str__(self):
        return str(self.alpha)
     
    def __repr__(self):
        return repr(self.alpha)
        
    def __default__(self):
        self.alpha = None
        return self.alpha
    
    def default(self):
        return self.__default__()
    
    def alphabets(self):
        return """Alphabets: 'chemical_2_1', 'chemical_2_2', 'chemical_2_3', 
        'chemical_2_4', 'chemical_3_1', 'chemical_3_2', 'chemical_3_3', 
        'chemical_4', 'chemical_5', 'chemical_6', 'murphy_2', 'murphy_4', 
        'murphy_8', 'murphy_10', 'murphy_15', 'wang_2', 'wang_3', 'wang_5_1', 
        'wang_5_2', 'li_3', 'li_4', 'li_5' or 'li_10'."""
    
    @staticmethod
    def alphabetReduct(aminoacid, method):
        """
        Performs alphabet reduction. Extended conversion is 
        performed to ensure that any alphabet character is
        processed.
        Note: All characters other than the 20 amino acids
        (IUPAC alpahbet) are considered gaps ('-') for 
        alphabet reduction.
        """
        
        if method == None:
            return extended[aminoacid]
        elif method == "chemical_2_1":
            new = extended[aminoacid]
            return chemical_2_1[new]
        elif method == "chemical_2_2":
            new = extended[aminoacid]
            return chemical_2_2[new]
        elif method == "chemical_2_3":
            new = extended[aminoacid]
            return chemical_2_3[new]
        elif method == "chemical_2_4":
            new = extended[aminoacid]
            return chemical_2_4[new]
        elif method == "chemical_3_1":
            new = extended[aminoacid]
            return chemical_3_1[new]
        elif method == "chemical_3_2":
            new = extended[aminoacid]
            return chemical_3_2[new]
        elif method == "chemical_3_3":
            new = extended[aminoacid]
            return chemical_3_3[new]
        elif method == "chemical_4":
            new = extended[aminoacid]
            return chemical_4[new]
        elif method == "chemical_5":
            new = extended[aminoacid]
            return chemical_5[new]
        elif method == "chemical_6":
            new = extended[aminoacid]
            return chemical_6[new]
        elif method == "murphy_2":
            new = extended[aminoacid]
            return murphy_2[new]
        elif method == "murphy_4":
            new = extended[aminoacid]
            return murphy_4[new]
        elif method == "murphy_8":
            new = extended[aminoacid]
            return murphy_8[new]
        elif method == "murphy_10":
            new = extended[aminoacid]
            return murphy_10[new]             
        elif method == "murphy_15":
            new = extended[aminoacid]
            return murphy_15[new]
        elif method == "wang_2":
            new = extended[aminoacid]
            return wang_2[new]
        elif method == "wang_3":
            new = extended[aminoacid]
            return wang_3[new]
        elif method == "wang_5_1":
            new = extended[aminoacid]
            return wang_5_1[new]
        elif method == "wang_5_2":
            new = extended[aminoacid]
            return wang_5_2[new]
        elif method == "li_3":
            new = extended[aminoacid]
            return li_3[new]
        elif method == "li_4":
            new = extended[aminoacid]
            return li_4[new]
        elif method == "li_5":
            new = extended[aminoacid]
            return li_5[new]
        elif method == "li_10":
            new = extended[aminoacid]
            return li_10[new]
        elif method == "custom":
            new = extended[aminoacid]
            return cAlphabet[new]
        else:
            return extended[aminoacid]


class Coevolution(object):
    """
    pycma.Coevolution object
        
    Methods and attributes of the class Coevolution.
    Computes coevolution analysis using one of the available 
    methods.
    
    Coevolution Methods: see Methods below
    
    Attributes: id1, id2, file1, file2, form, alignment1, alignment2,
                sequences1, sequences2, info1, info2, norm, model, 
                groups*, normalized*
    
    Methods: methods, mi_normalizations, cc_matrices,
                        MI, OMES, CC, QUARTETS, SCA, ELSC 
    
    Usage:
    >>> import pycma
    >>> p = pycma.Alignment('example1.fasta','example2.fasta','fasta')
    >>> pre = pycma.PreProcessing(p)
    >>> pre.gaps_in_query(False)
    >>> c1 = pycma.Coevolution(pre)
    >>> c1.methods()
    >>> c1.mi_normalizations()
    >>> c1.cc_matrices()
    >>> c1.MI('MI')
    >>> c2 = pycma.Coevolution(pre)
    >>> c2.SCA()
    """
    def __init__(self, alignments):
        try:
            assert len(alignments.alignment1) == len(alignments.alignment2)
        except:
            raise SameNumberSequencesError()
        
        self.model = alignments.model
        if self.model == "intra":
            print "Intra-protein Coevolution Analysis"
            sys.stdout.flush()
        elif self.model == "inter":
            print "Inter-protein Coevolution Analysis"
            sys.stdout.flush()
        else: pass
        
        # attributes
        self.alignment1 = alignments.alignment1
        self.alignment2 = alignments.alignment2
        self.sequences1 = alignments.sequences1
        self.sequences2 = alignments.sequences2
        self.info1 = alignments.info1
        self.info2 = alignments.info2 
        self.id1 = alignments.id1
        self.id2 = alignments.id2
        self.file1 = alignments.file1
        self.file2 = alignments.file2
        self.form = alignments.form
        self.norm = None
        self.matrix = None
        self.groups = 1
        self.normalized = True
        
    def __name__(self):
        return "Coevolution Object"
    
    def __str__(self):
        return str(self.alignment1), str(self.alignment2)
    
    def __repr__(self):
        return repr(self.alignment1), repr(self.alignment2)
    
    def methods(self):
        return """Methods: MI, OMES, CC, QUARTETS, SCA and ELSC."""
    
    def mi_normalizations(self):
        return normalizationsMI(self.norm).norm
    
    def cc_matrices(self):
        return matricesCC(self.matrix).matrix
    
    @staticmethod
    def mi(i, j, cols1, cols2, pD1, pD2, entropy):
        """
        Mutual Information - MI.
        MI divided by pair Entropy - NMI.
        NMI with Background Distribution - NMIB.
        NMI with Physicochemical Properties - NMIP.
        NMI with Physicochemical Properties - NMIP.
        NMI with Background Distribution and Physicochemical Properties - NMIBP.
        """
        
        col1, col2 = cols1[i], cols2[j]
        n = len(col1)
        assert n == len(col2)
        score = 0
        pairs = [col1[k] + col2[k] for k in range(n) \
                 if col1[k] in aa and col2[k] in aa and col1[k] != '-' and col2[k] != '-']
        pL = sorted(list(set(pairs)))
        for p in pL:
            pXY = pairs.count(p) * 1.0 / n
            pX = pD1[i][p[0]]
            pY = pD2[j][p[1]]
            inside = (pXY * 1.0) / (pX * pY)
            outside = pXY * log20(inside)
            score += outside
        if entropy is False:
            return score
        else:
            mi = score
            entropy = 0
            for p in pL: 
                pXY = pairs.count(p) * 1.0 / n
                inside = pXY
                outside = pXY * log20(inside)
                entropy += outside
            entropy = -entropy
            if entropy == 0.0:
                score = 0.0
            else: 
                score = mi / entropy
            return score
    
    @staticmethod
    def rcw(mi, i_all, all_j, n):
        "Row and Column Weighed MI - RCW"
        
        bottom = (i_all + all_j - 2.0 * mi) / (n - 1)
        if bottom == 0.0:
            score = 0.0
        else: 
            score = mi / bottom
        
        return score
    
    @staticmethod
    def mip(mi, i_all, all_j, mi_all):
        "Normalized MI - MIp"
        
        if mi_all == 0.0:
            apc = 0.0
        else: 
            apc = (i_all * all_j * 1.0) / (mi_all)
            score = mi - apc
        
        return score
    
    @staticmethod
    def omes(column1, column2):
        "Observed Minus Expected Squared - OMES."
     
        assert len(column1) == len(column2)
        
        L = []
        Nvalid = []
        Cxi = []
        Cyj = []
        for i, j in zip(column1, column2):
            if i in aa and j in aa and i != '-' and j != '-':
                value = [i, j]
                Nvalid.append(value)
                Cxi.append(i)
                Cyj.append(j)
                if value not in L:
                    L.append(value)
    
        len_Nvalid = len(Nvalid)
        score = 0.0
        for value in L:
            Nobs = Nvalid.count(value)
            i = value[0]
            j = value[1]
            Ci = Cxi.count(i)
            Cj = Cyj.count(j)
            Nex = Ci * Cj / len_Nvalid    
            top = (Nobs - Nex) ** 2
            score += top * 1.0 / len_Nvalid
        
        return score
    
    @staticmethod
    def cc(d_matrix1, d_matrix2, N, negative_correlation):
        """
        Correlation Coefficient, Pearson's Correlation, McBASC, Spearman's rank
        correlation, etc. - CC.
        """
        
        assert len(d_matrix1) == len(d_matrix2)
        
        sigma_i = std(d_matrix1)
        av_Si = mean(d_matrix1)
        
        sigma_j = std(d_matrix2)
        av_Sj = mean(d_matrix1)
        
        Si = []
        for i in (d_matrix1):
            Si.append(i - av_Si)
            
        Sj = []
        for j in (d_matrix2):
            Sj.append(j - av_Sj)
        
        top = 0.0
        for i, j in zip(Si, Sj):
            top += float(i * j)
    
        bottom = sigma_i * sigma_j
        if bottom == 0.0:
            score = 0.0
        else:
            if negative_correlation is True:
                score = abs((1.0 / N ** 2) * (top / bottom))
            else:
                score = (1.0 / N ** 2) * (top / bottom)   
        
        return score
    
    
    @staticmethod
    def quartets(column1, column2):
        """Quartets method - QUARTETS."""
     
        assert len(column1) == len(column2)
        
        score = 0.0
        
        pairs = []
        for i, j in zip(column1, column2):
            value = [i, j]
            pairs.append(value)
        
        N = len(column1)    
        for i in column1:
            for j in column2:
                if i in aa and j in aa and i != '-' and j != '-':
                    Pix = column1.count(i) * 1.0 / N
                    Piy = column2.count(i) * 1.0 / N
                    Pjx = column1.count(j) * 1.0 / N
                    Pjy = column2.count(j) * 1.0 / N
                    val = [i, j]
                    Dmin = pairs.count(val)
                    Dif = 1.0 * (len(pairs) - Dmin)
                    if Dif != 0.0:
                        DQmin = Dmin * 1.0 / Dif
                    else:
                        DQmin = 0.0
                    try :
                        if ((Pix * Pjy > Piy * Pjx) and ((Pix > Dmin) or (Pjy > Dmin)) or\
                        (Pix * Pjy < Piy * Pjx) and ((Piy > Dmin) or (Pjx > Dmin)))\
                        and\
                       (((Pix * Pjy) * 1.0 / (Piy * Pjx) > DQmin) or\
                        ((Piy * Pjx) * 1.0 / (Pix * Pjy) > DQmin)):
                            score += 1.0
                    except:
                        score += 0.0
        return score / 400.0 # matrix 20 x 20

    @staticmethod
    def sca(column1, column2, pos_j, columns2):
        """Statistical Coupling Analysis - SCA."""
        
        assert len(column1) == len(column2)
        
        new_alignment = subAlignment(column1, columns2)
        i = column1
        j1 = column2
        j2 = new_alignment[pos_j]
        
        inside = 0.0
        for x in i:
            if x in aa and x != '-':
                Pxj1 = j1.count(x) * 1.0 / len(j1)
                Pxj2 = j2.count(x) * 1.0 / len(j2)
                if Pxj2 != 0.0:
                    inside += (ln(Pxj2) - Pxj1) ** 2
                else: 
                    inside += 0.0
                
        score = sqrt(inside)
        return score
    
    @staticmethod
    def elsc(column1, column2, pos_j, columns2):
        """Explicit Likelihood of Subset Covariation - ELSC."""
        
        assert len(column1) == len(column2)
        
        new_alignment = subAlignment2(column1, column2, columns2)
        i = column1
        j1 = column2
        j2 = new_alignment[pos_j]
               
        comb_x = []
        comb_all = []
        for x in i:
            if x in aa and x != '-':
                Nxj1 = j1.count(x)
                nxj2 = j2.count(x)
                Nall = len(j1)
                nall = len(j2)
                mxj = int(round((Nxj1 * 1.0 / Nall) * nall))
                top = long(factorial(Nxj1))
                bott1 = factorial(nxj2) * factorial(Nxj1 - nxj2)
                bott2 = factorial(mxj) * factorial(Nxj1 - mxj)
                comb_x.append(top / bott1)
                comb_all.append(top / bott2)          
        
        product = 1.0
        for k, l in zip(comb_x, comb_all):    
            product *= (k * 1.0 / l) 
            
        if product != 0.0:
            score = -ln(product)
        else: 
            score = 0.0
        
        return score
    
    def MI(self, normalization='MI'):
        """
        Mutual Information (MI) as an interpretation of [1,2,3,4]
            
        Attributes: normalization, groups*, normalized*
        
        Normalizations: 'NMI', 'RCW', 'MIp', 'NMIB', 'NMIP', and 'NMIBP'
            -NMI: MI divided by pair Entropy [5,6]       
            -RCW: Row and Column Weighed MI [7]
            -MIp: Normalized MI [8]
            -NMIB: NMI with Background Distribution (B) [9]
            -NMIP: NMI with Physicochemical Properties (P) [9]
            -NMIBP: NMI with B and P [9]
        
        References:
        [1] Korber, B. T., Farber, R. M., Wolpert, D. H., & Lapedes, a S. (1993). 
            Covariation of mutations in the V3 loop of human immunodeficiency virus 
            type 1 envelope protein: an information theoretic analysis. Proceedings 
            of the National Academy of Sciences of the United States of America, 
            90(15), 7176-80.
        [2] Clarke, N. D. (1995). Covariation of residues in the homeodomain 
            sequence family. Protein science: a publication of the Protein Society, 
            4(11), 2269-78.
        [3] Halperin, I., Wolfson, H., & Nussinov, R. (2006). Correlated Mutations: 
            Advances and Limitations . A Study on Fusion Proteins and on the Cohesin-
            Dockerin Families. Proteins: Structure, Function, and Bioinformatics, 
            63(4), 832-845.
        [4] Horner, D. S., Pirovano, W., & Pesole, G. (2007). Correlated substitution 
            analysis and the prediction of amino acid structural contacts. Briefings 
            in bioinformatics, 9(1), 46-56.
        [5] Gloor, G. B., Martin, L. C., Wahl, L. M., & Dunn, S. D. (2005). Mutual 
            information in protein multiple sequence alignments reveals two classes 
            of coevolving positions. Biochemistry, 44(19), 7156-65. 
        [6] Martin, L. C., Gloor, G. B., Dunn, S. D., & Wahl, L. M. (2005). Using 
            information theory to search for co-evolving residues in proteins. 
            Bioinformatics, 21(22), 4116-24.
        [7] Gouveia-Oliveira, R., & Pedersen, A. G. (2007). Finding coevolving amino 
            acid residues using row and column weighting of mutual information and 
            multi-dimensional amino acid representation. Algorithms for molecular biology:
            AMB, 2, 12.
        [8] Dunn, S. D., Wahl, L. M., & Gloor, G. B. (2008). Mutual information without 
            the influence of phylogeny or entropy dramatically improves residue contact 
            prediction. Bioinformatics, 24(3), 333-40.
        [9] Gao, H., Dou, Y., Yang, J., & Wang, J. (2011). New methods to measure residues 
            coevolution in proteins. BMC bioinformatics, 12(1), 206. 
        """
        # attributes
        self.norm = normalizationsMI(normalization).norm
        
        # preparing the columns for analysis
        self.score = dict()
        alignment1 = [e for e in self.sequences1]
        columns1 = transpose(alignment1)
    
        alignment2 = [e for e in self.sequences2]
        columns2 = transpose(alignment2)
        
        # coevolution measures
        coevolution = self.norm
        
        # Mutual Information
        if coevolution == "MI":
            Flash('Mutual Information')
            
            pD1 = probabilityDict(columns1)
            pD2 = probabilityDict(columns2)
            entropy = False
            
            for i in range(len(columns1)):
                for j in range(len(columns2)):
                    if self.model == "intra":
                        if i != j:
                            self.score[(i, j)] = Coevolution.mi(i, j,
                                        columns1, columns2, pD1, pD2, entropy)
                    else:
                        self.score[(i, j)] = Coevolution.mi(i, j,
                                        columns1, columns2, pD1, pD2, entropy)
            
        # MI divided by pair Entropy
        elif coevolution == "NMI":
            Flash('MI divided by pair Entropy')
            
            pD1 = probabilityDict(columns1)
            pD2 = probabilityDict(columns2)
            entropy = True
            
            for i in range(len(columns1)):
                for j in range(len(columns2)):
                    if self.model == "intra":
                        if i != j:
                            self.score[(i, j)] = Coevolution.mi(i, j,
                                        columns1, columns2, pD1, pD2, entropy)
                    else:
                        self.score[(i, j)] = Coevolution.mi(i, j,
                                        columns1, columns2, pD1, pD2, entropy)

        # Row and Column Weighed MI                  
        elif coevolution == "RCW":
            Flash('Row and Column Weighed MI')
            
            pD1 = probabilityDict(columns1)
            pD2 = probabilityDict(columns2)
            entropy = False
         
            i_all = dict()
            all_j = dict()
            for i in range(len(columns1)):
                v_i = 0
                for j in range(len(columns2)):
                    v_i += Coevolution.mi(i, j, columns1, columns2,
                                             pD1, pD2, entropy)
                    i_all[i] = v_i
    
            for j in range(len(columns2)):
                v_j = 0
                for i in range(len(columns1)):
                    v_j += Coevolution.mi(i, j, columns1, columns2,
                                             pD1, pD2, entropy)
                    all_j[j] = v_j
            
            column = columns1[0]
            n = len(column)
            for i in range(len(columns1)):
                for j in range(len(columns2)):
                    if self.model == "intra":
                        if i != j:
                            mi_score = Coevolution.mi(i, j, columns1, columns2,
                                           pD1, pD2, entropy)    
                            self.score[(i, j)] = Coevolution.rcw(mi_score,
                                                   i_all[i], all_j[j], n)
                    else:
                        mi_score = Coevolution.mi(i, j, columns1, columns2,
                                           pD1, pD2, entropy)    
                        self.score[(i, j)] = Coevolution.rcw(mi_score,
                                                   i_all[i], all_j[j], n)
        # Normalized MI                  
        elif coevolution == "MIp":
            Flash('Normalized MI - Very Slow!')
            
            pD1 = probabilityDict(columns1)
            pD2 = probabilityDict(columns2)
            mi_all = []
            entropy = False
            
            for i in range(len(columns1)):
                for j in range(len(columns2)):
                    if self.model == "intra":
                        if i != j:
                            mi_all.append(Coevolution.mi(i, j,
                                        columns1, columns2, pD1, pD2, entropy))
                    else:
                        mi_all.append(Coevolution.mi(i, j,
                                        columns1, columns2, pD1, pD2, entropy))
            
            mi_all = float(mean(mi_all))
            
            column = columns1[0]
            n = len(column)
            for i in range(len(columns1)):
                for j in range(len(columns2)):
                    if self.model == "intra":
                        if i != j:
                            i_all = []
                            all_j = []
                            for ii in range(len(columns1)):
                                if ii != i:
                                    i_all.append(Coevolution.mi(i, ii, columns1, columns1,
                                                                pD1, pD1, entropy))
                            for jj in range(len(columns2)):
                                if jj != j:
                                    all_j.append(Coevolution.mi(j, jj, columns2, columns2,
                                                                pD2, pD2, entropy))
                            i_all = float(mean(i_all))
                            all_j = float(mean(all_j))
                            mi_score = Coevolution.mi(i, j, columns1, columns2,
                                           pD1, pD2, entropy)    
                            self.score[(i, j)] = Coevolution.mip(mi_score,
                                                   i_all, all_j, mi_all)
                    else:
                        i_all = []
                        all_j = []
                        for ii in range(len(columns1)):
                            if ii != i:
                                i_all.append(Coevolution.mi(i, ii, columns1, columns1,
                                                            pD1, pD1, entropy))
                        for jj in range(len(columns2)):
                            if jj != j:
                                all_j.append(Coevolution.mi(j, jj, columns2, columns2,
                                                            pD2, pD2, entropy))
                        i_all = float(mean(i_all))
                        all_j = float(mean(all_j))
                        mi_score = Coevolution.mi(i, j, columns1, columns2,
                                       pD1, pD2, entropy)    
                        self.score[(i, j)] = Coevolution.mip(mi_score,
                                               i_all, all_j, mi_all)
            
            
                        
        # NMI with Background Distribution                  
        elif coevolution == "NMIB":
            Flash('NMI with Background Distribution')
            
            pD1 = probabilityDictBackDist(columns1)
            pD2 = probabilityDictBackDist(columns2)
            entropy = True
            
            for i in range(len(columns1)):
                for j in range(len(columns2)):
                    if self.model == "intra":
                        if i != j:
                            self.score[(i, j)] = Coevolution.mi(i, j,
                                        columns1, columns2, pD1, pD2, entropy)
                    else:
                        self.score[(i, j)] = Coevolution.mi(i, j,
                                        columns1, columns2, pD1, pD2, entropy)
            
        # NMI with Physicochemical Properties                 
        elif coevolution == "NMIP":
            Flash('NMI with Physicochemical Properties')
            
            pD1 = probabilityDictPhysicProp(columns1)
            pD2 = probabilityDictPhysicProp(columns2)
            entropy = True
         
            for i in range(len(columns1)):
                for j in range(len(columns2)):
                    if self.model == "intra":
                        if i != j:
                            self.score[(i, j)] = Coevolution.mi(i, j,
                                        columns1, columns2, pD1, pD2, entropy)
                    else:
                        self.score[(i, j)] = Coevolution.mi(i, j,
                                        columns1, columns2, pD1, pD2, entropy)
        
        # NMI with Background Distribution and Physicochemical Properties                 
        elif coevolution == "NMIBP":
            Flash('NMI with Background Distribution and Physicochemical Properties')
            
            pD1 = probabilityDictBackPhysic(columns1)
            pD2 = probabilityDictBackPhysic(columns2)
            entropy = True
         
            for i in range(len(columns1)):
                for j in range(len(columns2)):
                    if self.model == "intra":
                        if i != j:
                            self.score[(i, j)] = Coevolution.mi(i, j,
                                        columns1, columns2, pD1, pD2, entropy)
                    else:
                        self.score[(i, j)] = Coevolution.mi(i, j,
                                        columns1, columns2, pD1, pD2, entropy)
        
        else: pass
        
    
    def OMES(self):
        """
        Observed Minus Expected Squared (OMES) as an interpretation of [1,2,3,4,5]
        
        Attributes: 
        
        References:
        [1] Larson, S. M., Di Nardo, a a, & Davidson, a R. (2000). Analysis of covariation
            in an SH3 domain sequence alignment: applications in tertiary contact 
            prediction and the design of compensating hydrophobic core substitutions. 
            Journal of molecular biology, 303(3), 433-46.
        [2] Kass, I., & Horovitz, A. (2002). Mapping pathways of allosteric communication 
            in GroEL by analysis of correlated mutations. Proteins, 48(4), 611-7.
        [3] Fodor, A. a, & Aldrich, R. W. (2004). Influence of conservation on calculations 
            of amino acid covariance in multiple sequence alignments. Proteins, 56(2), 211-21.
        [4] Halperin, I., Wolfson, H., & Nussinov, R. (2006). Correlated Mutations: 
            Advances and Limitations . A Study on Fusion Proteins and on the Cohesin-
            Dockerin Families. Proteins: Structure, Function, and Bioinformatics, 
            63(4), 832-845.
        [5] Horner, D. S., Pirovano, W., & Pesole, G. (2007). Correlated substitution 
            analysis and the prediction of amino acid structural contacts. Briefings 
            in bioinformatics, 9(1), 46-56.
        """
        # preparing the columns for analysis
        self.score = dict()
        alignment1 = [e for e in self.sequences1]
        columns1 = transpose(alignment1)
    
        alignment2 = [e for e in self.sequences2]
        columns2 = transpose(alignment2)
        
        # coevolution measures
        Flash('Observed Minus Expected Squared')
            
        for i in range(len(columns1)):
            for j in range(len(columns2)):
                if self.model == "intra":
                    if i != j:
                        self.score[(i, j)] = Coevolution.omes(columns1[i], columns2[j])
                else:
                    self.score[(i, j)] = Coevolution.omes(columns1[i], columns2[j])
                    
    
    def CC(self, self_comparison=False, rank_order=False,
                    negative_correlation=True, matrix='MCLACHLAN'):
        """
        Correlation Coefficient (CC) as an interpretation of [1,2,3,4,5]
        
        Attributes: self_comparison, rank_order, negative_correlation, matrix
                    
        self_comparison=False [6]
        self_comparison=True [1] (Pearson Correlation)
        
        rank_order=False [1]
        rank_order=True [7] (Spearman Rank Correlation)
        
        negative_correlation=False [1]
        negative_correlation=True [3] (McLachlan Based Substitution Correlation)
        
        matrix='MCLACHLAN' [8] (default), 'BLOSUM62' [9], 'PAM250' [10], 
                                                        'CPVN' [11] or 'CLM' [12]
        
        References:
        [1] Gobel, U., Sander, C., Schneider, R., & Valencia, A. (1994). Correlated Mutations 
            and Residue Contacts in Proteins. PROTEINS: Structure, Function, and Genetics, 18, 
            309-317.
        [2] Neher, E. (1994). How frequent are correlated changes in families of protein 
            sequences? Proceedings of the National Academy of Sciences of the United States 
            of America, 91(1), 98-102. 
        [3] Fodor, A. a, & Aldrich, R. W. (2004). Influence of conservation on calculations 
            of amino acid covariance in multiple sequence alignments. Proteins, 56(2), 211-21.
        [4] Halperin, I., Wolfson, H., & Nussinov, R. (2006). Correlated Mutations: 
            Advances and Limitations . A Study on Fusion Proteins and on the Cohesin-
            Dockerin Families. Proteins: Structure, Function, and Bioinformatics, 
            63(4), 832-845.
        [5] Horner, D. S., Pirovano, W., & Pesole, G. (2007). Correlated substitution 
            analysis and the prediction of amino acid structural contacts. Briefings 
            in bioinformatics, 9(1), 46-56.
        [6] Pollock, D., & Taylor, W. (1997). Effectiveness of correlation analysis in 
            identifying protein residues undergoing correlated evolution. Protein 
            Engineering, 10(6), 647. 
        [7] Pazos, F., Helmer-Citterich, M., Ausiello, G., & Valencia, a. (1997). Correlated 
            mutations contain information about protein-protein interaction. Journal of 
            molecular biology, 271(4), 511-23.
        [8] McLachlan A. D. (1971). Tests for Comparing Related Amino-acid Sequences 
            Cytochrome c and Cytochrome c551. Journal of Molecular Biology, 61, 409-424.
        [9] Henikoff, S., & Henikoff, J. G. (1992). Amino acid substitution matrices from 
            protein blocks. Proceedings of the National Academy of Sciences of the United 
            States of America, 89(22), 10915-9. 
        [10] Dayhoff.M.O. and Eck.R.V. (1968). In Dayhoff.M.O. (ed.), Atlas of Protein 
            Sequence and Structure. National Biomedical Research Foundation, Silver 
            Spring, MD, Vol. 3, pp. 33-41.
        [11] Glaser, F., Steinberg, D. M., Vakser, I. A., & Ben-tal, N. (2001). Residue 
            Frequencies and Pairing Preferences at Protein â Protein Interfaces. PROTEINS: 
            Structure, Function, and Genetics, 102(43), 89-102.
        [12] Singer, M. S., Vriend, G., & Bywater, R. P. (2002). Prediction of protein 
            residue contacts with a PDB-derived likelihood matrix. Protein engineering, 
            15(9), 721-5.
        
        Usage:
            >>> c = Coevolution(p)
            >>> c.CC()
        """
        # attributes
        if self_comparison is True or self_comparison is False:
            self.self_comparison = self_comparison
        else: 
            raise BooleanError()
        
        if rank_order is True or rank_order is False:
            self.rank_order = rank_order
        else:
            raise BooleanError()
        
        if negative_correlation is True or negative_correlation is False:
            self.negative_correlation = negative_correlation
        else:
            raise BooleanError()
        
        self.matrix = matricesCC(matrix).matrix
        
        # preparing the columns for analysis
        self.score = dict()
        alignment1 = [e for e in self.sequences1]
        columns1 = transpose(alignment1)
    
        alignment2 = [e for e in self.sequences2]
        columns2 = transpose(alignment2)
        
        # coevolution measures
        Flash("Correlation Coefficient - Very Slow!")
            
        N = len(columns1[0])
        for i in range(len(columns1)):
            for j in range(len(columns2)):
                if self.model == "intra":
                    if i != j:
                        d_matrix1 = twoDimensionalMatrix(
                                            columns1[i], self.matrix,
                                            self.self_comparison, self.rank_order)
                        d_matrix2 = twoDimensionalMatrix(
                                            columns2[j], self.matrix,
                                            self.self_comparison, self.rank_order)
                        self.score[(i, j)] = Coevolution.cc(
                                            d_matrix1, d_matrix2,
                                            N, self.negative_correlation)
                else:
                    d_matrix1 = twoDimensionalMatrix(
                                            columns1[i], self.matrix,
                                            self.self_comparison, self.rank_order)
                    d_matrix2 = twoDimensionalMatrix(
                                            columns2[j], self.matrix,
                                            self.self_comparison, self.rank_order)
                    self.score[(i, j)] = Coevolution.cc(
                                            d_matrix1, d_matrix2,
                                            N, self.negative_correlation)

   
    def QUARTETS(self):
        """
        QUARTETS as an interpretation of [1,2]
        
        Attributes: groups*, normalized*
        
        References:
        [1] Galitsky, B. (2003). Revealing the set of mutually correlated positions for 
            the protein families of immunoglobulin fold. In silico biology, 3(3), 241-64.
        [2] Halperin, I., Wolfson, H., & Nussinov, R. (2006). Correlated Mutations: 
            Advances and Limitations . A Study on Fusion Proteins and on the Cohesin-
            Dockerin Families. Proteins: Structure, Function, and Bioinformatics, 
            63(4), 832-845.
        """
        # preparing the columns for analysis
        self.score = dict()
        alignment1 = [e for e in self.sequences1]
        columns1 = transpose(alignment1)
    
        alignment2 = [e for e in self.sequences2]
        columns2 = transpose(alignment2)
        
        # coevolution measures  
        Flash('Quartets - Slow!')
            
        for i in range(len(columns1)):
            for j in range(len(columns2)):
                if self.model == "intra":
                    if i != j:
                        self.score[(i, j)] = Coevolution.quartets(
                                                columns1[i], columns2[j])
                else:
                    self.score[(i, j)] = Coevolution.quartets(
                                                columns1[i], columns2[j])
                
    
    def SCA(self):
        """
        Statistical Coupling Analysis (SCA) as an interpretation of [1,2,3,4]
        
        Attributes:
        
        References:
        [1] Lockless, S. W., & Ranganathan, R. (1999). Evolutionarily conserved pathways 
            of energetic connectivity in protein families. Science, 286(5438), 295-9. 
        [2] Süel, G. M., Lockless, S. W., Wall, M. a, & Ranganathan, R. (2003). 
            Evolutionarily conserved networks of residues mediate allosteric communication 
            in proteins. Nature structural biology, 10(1), 59-69. 
        [3] Halabi, N., Rivoire, O., Leibler, S., & Ranganathan, R. (2009). Protein sectors: 
            evolutionary units of three-dimensional structure. Cell, 138(4), 774-86. 
        [4] Halperin, I., Wolfson, H., & Nussinov, R. (2006). Correlated Mutations: 
            Advances and Limitations . A Study on Fusion Proteins and on the Cohesin-
            Dockerin Families. Proteins: Structure, Function, and Bioinformatics, 
            63(4), 832-845.
        """
        # preparing the columns for analysis
        self.score = dict()
        alignment1 = [e for e in self.sequences1]
        columns1 = transpose(alignment1)
    
        alignment2 = [e for e in self.sequences2]
        columns2 = transpose(alignment2)
        
        # coevolution measures
        Flash('Statistical Coupling Analysis')
              
        for i in range(len(columns1)):
            for j in range(len(columns2)):
                if self.model == "intra":
                    if i != j:
                        self.score[(i, j)] = Coevolution.sca(
                                                columns1[i], columns2[j],
                                              j, columns2)
                else:
                    self.score[(i, j)] = Coevolution.sca(
                                                columns1[i], columns2[j],
                                              j, columns2)
        
    
    def ELSC(self):
        """
        Explicit Likelihood of Subset Covariation (ELSC) as an interpretation of [1,2,3]
        
        Attributes:
        
        References:
        [1] Dekker, J. P., Fodor, A., Aldrich, R. W., & Yellen, G. (2004). A perturbation-
            based method for calculating explicit likelihood of evolutionary co-variance 
            in multiple sequence alignments. Bioinformatics, 20(10), 1565-72.
        [2] Halperin, I., Wolfson, H., & Nussinov, R. (2006). Correlated Mutations: 
            Advances and Limitations . A Study on Fusion Proteins and on the Cohesin-
            Dockerin Families. Proteins: Structure, Function, and Bioinformatics, 
            63(4), 832-845.
        [3] Horner, D. S., Pirovano, W., & Pesole, G. (2007). Correlated substitution 
            analysis and the prediction of amino acid structural contacts. Briefings 
            in bioinformatics, 9(1), 46-56.
        """
        # preparing the columns for analysis
        self.score = dict()
        alignment1 = [e for e in self.sequences1]
        columns1 = transpose(alignment1)
    
        alignment2 = [e for e in self.sequences2]
        columns2 = transpose(alignment2)
        
        # coevolution measures
        Flash('Explicit Likelihood of Subset Covariation') 
          
        for i in range(len(columns1)):
            for j in range(len(columns2)):
                if self.model == "intra":
                    if i != j:
                        self.score[(i, j)] = Coevolution.elsc(
                                                columns1[i], columns2[j],
                                               j, columns2)
                else:
                    self.score[(i, j)] = Coevolution.elsc(
                                                columns1[i], columns2[j],
                                               j, columns2)
        


_mi_builtins = ['MI', 'NMI', 'RCW', 'MIp', 'NMIB', 'NMIP', 'NMIBP']

class normalizationsMI(Coevolution):
    """
    pycma.normalizationsMI object
    
    Normalizations: 'MI', 'NMI', 'RCW', 'MIp', 'NMIB', 'NMIP' and 'NMIBP'
    
    Attributes: norm
    
    Methods: mi_normalizations
    """
    def __init__(self, norm):
        if norm not in _mi_builtins:
            raise CoevolutionError()
        else:
            self.norm = norm
    
    def __name__(self):
        return "normalizationsMI Object"
    
    def __str__(self):
        return str(self.norm)
     
    def __repr__(self):
        return repr(self.norm)
    
    def mi_normalizations(self):
        return """Methods: 'MI', 'NMI', 'RCW', 'MIp', 'NMIB', 'NMIP' and 'NMIBP'."""


_matrix_builtins = ['BLOSUM62', 'PAM250', 'MCLACHLAN', 'CPVN', 'CLM']

class matricesCC(Coevolution):
    """
    pycma.matrixCC object
    
    Matrix: 'BLOSUM62', 'PAM250' and 'MCLACHLAN'
    
    Attributes: matrix
    
    Methods: cc_matrices
    """
    def __init__(self, matrix):
        if matrix not in _matrix_builtins:
            raise CoevolutionError()
        else:
            self.matrix = matrix
    
    def __name__(self):
        return "matrixCC Object"
    
    def __str__(self):
        return str(self.matrix)
     
    def __repr__(self):
        return repr(self.matrix)
    
    def cc_matrices(self):
        return """Matrix: 'BLOSUM62', 'PAM250', and 'MCLACHLAN' (default)."""


class Results(object):
    """
    pycma.Results object
        
    Methods and attributes of the class Results. Object responsible
    for the output and analysis of results.    
    
    Attributes: id1, id2, file1, file2, form, alignment1, alignment2,
                sequences1, sequences2, info1, info2, score, 
                model
    
    Methods: pair_scores, residue_scores, stats, histogram, heatmap, 
            networkmap, write_file, formats, grouping, normalization, 
    
    Usage:
    >>> import pycma
    >>> p = pycma.Alignment('example1.fasta','example2.fasta','fasta')
    >>> pre = pycma.PreProcessing(p)
    >>> pre.gaps_in_query(False)
    >>> c1 = pycma.Coevolution(pre)
    >>> c1.methods()
    >>> c1.mi_normalizations()
    >>> c1.cc_matrices()
    >>> c1.MI('MI')
    >>> r = pycma.Results(c1)
    >>> r.normalization(True)
    >>> r.scores(20)
    >>> r.histogram('test1.png')
    >>> r.heatmap('test2.png')
    >>> r.write_file('example1.txt')
    >>> r.write_file('example2.txt', scores=20, form='plain')
    >>> c2 = pycma.Coevolution(pre)
    >>> c2.SCA()
    >>> r = pycma.Results(c2)
    >>> r.grouping(3)
    >>> r.normalization(True)
    >>> r.scores(50)
    >>> r.histogram('test3.png')
    >>> r.heatmap('test4.png')
    >>> r.write_file('example3.txt')
    >>> r.write_file('example4.txt', scores=50, form='basic')
    """
    def __init__(self, alignments):
        # attributes
        self.alignment1 = alignments.alignment1
        self.alignment2 = alignments.alignment2
        self.sequences1 = alignments.sequences1
        self.sequences2 = alignments.sequences2
        self.info1 = alignments.info1
        self.info2 = alignments.info2 
        self.id1 = alignments.id1
        self.id2 = alignments.id2
        self.file1 = alignments.file1
        self.file2 = alignments.file2
        self.form = alignments.form
        self.score = alignments.score
        self.model = alignments.model

    def __name__(self):
        return "Results Object"
    
    def __len__(self):
        return len(self.score)
    
    def __str__(self):
        return "Huge dictionary! Its better not print that out..."
    
    def __repr__(self):
        return "Huge dictionary! Its better not print that out..."
    
    def __getitem__(self, index):
        return self.score[index]
    
    def __getslice__(self, start, end):
        return self.score[start, end]
    
    def __iter__(self):
        return iter(self.score)
    
    def formats(self):
        return """formats: 'plain', 'tab' or 'basic'."""
    
    @staticmethod
    def groups(info, groups, length1, length2):
        """
        Algorithm to average the coevolution score at one position by
        the score of neighboring positions. This is achieved by the 
        ponderation 0.5*pair_score + 0.5*(mean(groups_score)).
        For simplicity the first and last scores are not averaged.
        Besides, groups number must be odd and if a even number is
        provided, it is automatically changed for the following odd
        number. 
        """
        scoreDict = dict()
        for i in range(length1):
            for j in range(length2):
                score = info[(i, j)]
                try:
                    groups_score = []
                    nb = groupsNumber(groups)
                    for l in range(1,nb):
                        pos_neighbor = info[(i+l, j+l)] 
                        neg_neighbor = info[(i-l, j-l)] 
                        groups_score.append(pos_neighbor)
                        groups_score.append(neg_neighbor)
                    mean_groups = mean(groups_score)
                    nscore = (0.5 * score) + (0.5 * mean_groups)
                except: 
                    nscore = score
                scoreDict[(i, j)] = nscore
        
        return scoreDict
    
    @staticmethod
    def normalizes(info, length1, length2): 
        """
        Normalizes Coevolution Scores. Scores range from 0 to 1.
        maximum(score) = 1 and minimum(score) = 0.
        """
        
        max_pos = []
        for i in range(length1):
            for j in range(length2):
                max_pos.append(info[(i, j)])
        max_val = max(max_pos)
        
        scoreDict = dict()       
        for i in range(length1):
            for j in range(length2):
                if info[(i, j)] > 0.0:
                    scoreDict[(i, j)] = info[(i, j)] * 1.0 / max_val
                else:
                    scoreDict[(i, j)] = 0.0
                    
        return scoreDict
    
    @staticmethod
    def bestResults(info, scores):
        "Creates a new list of best scores"
        
        results = dict2List(info)
        
        a = results
        sort = sorted(a, key=lambda a: a[2])
        length = len(sort)
        position = length - scores
        threshold = sort[position]
        
        top = []
        count = 0
        for l in results:
            res1 = l[0]
            res2 = l[1]
            score = float(l[2])
            value = [res1, res2, score]
            if score >= threshold[2]:
                count += 1 
                if count <= scores:
                    top.append(value)
                else: pass
        return top
    
    @staticmethod
    def drawHistogram(info, output):
        "Creates a histogram of scores"  
         
        results = dict2List(info)
           
        data = []
        info = []
        for l in results:
            res1 = l[0]
            res2 = l[1]
            score = float(l[2])
            value = [res1, res2, score]
            data.append(value)
            info.append(score)
    
        maxi = max(info)
        L = [t[2] for t in data]
        X = maxi
        pyplot.hist(L, bins=X * 50)
        ax = pyplot.axes()
        ax.set_xlabel('Score')
        ax.set_ylabel('Frequency')
        ax.set_xlim(0, X)
        pyplot.savefig(output)
    
    @staticmethod                
    def drawHeatmap(id1, id2, info, output): 
        "Creates a heatmap of scores"
                   
        results = dict2List(info)
        
        data = []  
        residue1 = []
        residue2 = [] 
        for l in results:
            res1 = l[0]
            res2 = l[1]
            score = float(l[2])
            value = [res1, res2, score]
            data.append(value)
            if res1 not in residue1:
                residue1.append(res1)
            if res2 not in residue2:
                residue2.append(res2)
           
        startX = int(data[0][0])
        startY = int(data[0][1])
        length = len(data)
        endX = int(data[length - 1][0])
        endY = int(data[length - 1][1])
        
        lenX = len(residue1)
        lenY = len(residue2)    
        heatmap = zeros((lenY + 1, lenX + 1))
        for i in range(length):
            X = int(data[i][0])
            Y = int(data[i][1])
            XY = float(data[i][2])
            heatmap[Y][X] = XY
                
        pyplot.figure()
        pyplot.pcolormesh(heatmap)
        pyplot.colorbar() 
        pyplot.axes().set_xlabel(id1)
        pyplot.axes().set_ylabel(id2)
        pyplot.axes().set_xlim(startX, endX)
        pyplot.axes().set_ylim(startY, endY)
        pyplot.savefig(output)
    
    @staticmethod
    def drawInteractionMap():
        pass
    
    @staticmethod
    def writeToFile(id1, id2, info, output, form=None): 
        "Writes results to file using one of the available formats"
                   
        results = info
        
        try:
            out = open(output, "w")
        except: 
            raise IOError, "Exception: Unable to open filename %s" % (output)
        for l in results:
            res1 = l[0]
            res2 = l[1]
            score = float(l[2])
            
            if form == None or form == 'plain':
                print >> out, "%s:Res %i - %s:Res %i: Score %f" % (id1, res1,
                                                          id2, res2,
                                                          score)
            elif form == 'tab':
                print >> out, "%s:Res\t%i\t-\t%s:Res\t%i:\tScore\t%f" % (id1, res1,
                                                          id2, res2,
                                                          score)
            elif form == 'basic':
                print >> out, res1, res2, score
            else: pass
        out.close()
    
    
    def scores(self, scores, form=None):
        """
        Gets the specified scores. formats can be selected to 
        change the way results are presented. 
        """
        
        results = Results.bestResults(self.score, scores)
        
        print "Best %i scores:" % scores
        if form != 'basic':
            for e in results:
                res1 = e[0]
                res2 = e[1]
                score = e[2]
                print "%s:Res %i - %s:Res %i: Score %f" % (self.id1, res1,
                                                          self.id2, res2,
                                                          score)
        else:
            print "%s:Res %s:Res Score" % (self.id1, self.id2)
            for e in results:
                res1 = e[0]
                res2 = e[1]
                score = e[2]
                print "%i %i %f" % (res1, res2, score)
        return
    
    def sats(self):
        """
        TODO: Defines commun statistical information.
        """
        pass
        
    def histogram(self, filename):
        """
        Writes histograms to disk in png form.
        """
        
        return Results.drawHistogram(self.score, filename)
    
    def heatmap(self, filename):
        """
        Writes heatmaps to disk in png form.
        """
        
        return Results.drawHeatmap(self.id1, self.id2, self.score, filename)
    
    def networkmap(self):
        """
        TODO: Writes interaction maps to disk in png form.
        """
        
        return Results.drawInteractionMap()
    
    def write_file(self, filename, scores=None, form=None):
        """
        Writes files to disk. If scores !=None, only the highest 
        scores are written. formats can be selected to change the
        way results are presented. 
        """
        
        form = Fileform(form).form
        
        if scores != None:
            results = Results.bestResults(self.score, scores)
            return Results.writeToFile(self.id1, self.id2,
                               results, filename, form)
        else:
            form = Fileform(form).form
            results = dict2List(self.score)
            return Results.writeToFile(self.id1, self.id2,
                               results, filename, form)
    
    def grouping(self, group=3):
        """
        if groups = n and n > 1 -> scores pairs as a ponderation
        of neighbor scores. For example position res:5 - res:60 
        for n = 3 would be considered as half the score of residues 
        5 - 60, and half score of the mean of 4 - 59 and 6 - 61.
        This is always simetric, and we recomend the use of odd
        numbers.
        """
        # checks input
        try:
            self.group = int(group)
        except:
            raise IntError()
        
        # preparing the columns for analysis
        alignment1 = [e for e in self.sequences1]
        columns1 = transpose(alignment1)
    
        alignment2 = [e for e in self.sequences2]
        columns2 = transpose(alignment2)
        
        # updates self.score
        Flash("Grouping %i residue(s) for each position" %(self.group))
        if self.group !=1:
            self.score = Results.groups(self.score, self.group,
                                              len(columns1), len(columns2))
        
    def normalization(self, normalized=True):
        """
        if normalized=True -> the scores are normalized [~0:1]
        The max(score) = 1 and min(score)=0 or negative
        else -> scoring schemes specific for each method (this
        include negative scores)
        """
        # checks input
        if normalized is True or normalized is False:
                self.normalized = normalized
        else:
            raise BooleanError()
        
        # preparing the columns for analysis
        alignment1 = [e for e in self.sequences1]
        columns1 = transpose(alignment1)
    
        alignment2 = [e for e in self.sequences2]
        columns2 = transpose(alignment2)
        
        # updates self.score
        Flash("Results normalized [~0:1]")
        if normalized:
            self.score = Results.normalizes(self.score, 
                                              len(columns1), len(columns2))


_file_form_builtins = ['plain', 'tab', 'basic']
    
class Fileform(Results):
    """
    pycma.Fileform object 
    
    formats: 'plain', 'tab' or 'basic'.
    
    Attributes: form
    
    Methods: update, default
    
    Usage:

    """
    def __init__(self, form=None):
        if form is not None and form not in _file_form_builtins:
            raise AlignformatError()
        elif form is not None and form in _file_form_builtins:
            self.form = form
        else:
            return self.__default__()
            
    def __name__(self):
        return "Fileform Object"
    
    def __str__(self):
        return str(self.form)
     
    def __repr__(self):
        return repr(self.form)
        
    def __default__(self):
        self.form = "plain"
        return
    
    def default(self):
        return self.__default__()


#-----------------------------------------------------------------------   

class IntError:
    """
    Input must be an integer.
    """
    def __init__(self):
        print IntError.__doc__
        
class StrError:
    """
    Input must be a string.
    """
    def __init__(self):
        print StrError.__doc__

class BooleanError:
    """
    Input must be a boolean; True or False.
    """
    def __init__(self):
        print BooleanError.__doc__

class DictError:
    """
    Input must be a dictionary with results.
    """
    def __init__(self):
        print DictError.__doc__
    
class AlignformatError:
    """
    See class Alignform for more details.
    """
    def __init__(self):
        print AlignformatError.__doc__

class FileformatError:
    """
    See class Fileform for more details.
    """
    def __init__(self):
        print FileformatError.__doc__

class AlphabetError:
    """ 
    See class Alphabet for more details.
    """
    def __init__(self):
        print AlphabetError.__doc__
        
class CoevolutionError:
    """
    See class Coevolution for more details.
    """
    def __init__(self):
        print CoevolutionError.__doc__
        
class ConservationError:
    """
    See class Conservation for more details.
    """
    def __init__(self):
        print ConservationError.__doc__

class NumberSequencesError:
    """
    Alignments must have at least n sequences.
    """
    def __init__(self, n=15):
        print "Alignments must have at least %s sequences." % n

class SameNumberSequencesError:
    """
    Alignments must have the same number of sequences.
    """
    def __init__(self):
        print SameNumberSequencesError.__doc__

  
#-----------------------------------------------------------------------    
# utility methods
def Flash(message):
    print message
    sys.stdout.flush()

def log20(n):  
    return log(n) * 1.0 / log(20)

def ln(n): 
    return log(n) * 1.0 / log(e)

def groupsNumber(number):
    "Checks if an int is odd and if not gets the next highest odd"
    if number & 1:
        return int(((number - 1) * 1.0 / 2) + 1)
    else:
        return int((number * 1.0 / 2) + 1)

def dict2List(info):
    "Converts a dictionary to a list"
    lst = []
    for i, j in sorted(info.keys()):
        res1 = int(i + 1)
        res2 = int(j + 1)
        score = float(round((info[(i, j)]), 8))
        value = [res1, res2, score]
        lst.append(value)
    return lst

def matchScore(alpha, beta, score_matrix):
    """
    Matches scores from a matrix. It gets the specific 
    scores from a list of lists (matrix).
    """
        
    alphabet = {}    
    alphabet["I"] = 0
    alphabet["V"] = 1
    alphabet["L"] = 2
    alphabet["F"] = 3
    alphabet["C"] = 4
    alphabet["M"] = 5
    alphabet["A"] = 6
    alphabet["G"] = 7
    alphabet["T"] = 8
    alphabet["S"] = 9
    alphabet["W"] = 10
    alphabet["Y"] = 11
    alphabet["P"] = 12
    alphabet["H"] = 13
    alphabet["E"] = 14
    alphabet["Q"] = 15
    alphabet["D"] = 16
    alphabet["N"] = 17
    alphabet["K"] = 18
    alphabet["R"] = 19
    try:
        alphabet["B"] = 20
        alphabet["Z"] = 21
        alphabet["X"] = 22
        alphabet["-"] = 23
    except:
        pass
    lut_x = alphabet[alpha]
    lut_y = alphabet[beta]
    
    return score_matrix[lut_x][lut_y]

def matchScoreMCLACHLAN(alpha, beta, score_matrix):
    """
    Matches scores from a matrix. It gets the specific 
    scores from a list of lists (matrix). Different 
    residue order.
    """
    
    alphabet = {}    
    alphabet["A"] = 0
    alphabet["R"] = 1
    alphabet["N"] = 2
    alphabet["D"] = 3
    alphabet["C"] = 4
    alphabet["Q"] = 5
    alphabet["E"] = 6
    alphabet["G"] = 7
    alphabet["H"] = 8
    alphabet["I"] = 9
    alphabet["L"] = 10
    alphabet["K"] = 11
    alphabet["M"] = 12
    alphabet["F"] = 13
    alphabet["P"] = 14
    alphabet["S"] = 15
    alphabet["T"] = 16
    alphabet["W"] = 17
    alphabet["Y"] = 18
    alphabet["V"] = 19
    lut_x = alphabet[alpha]
    lut_y = alphabet[beta]
    
    return score_matrix[lut_x][lut_y]

def matchAA(alpha):
    if alpha == 'A':
        beta = A
    elif alpha == 'R':
        beta = R
    elif alpha == 'N':
        beta = N
    elif alpha == 'D':
        beta = D
    elif alpha == 'C':
        beta = C
    elif alpha == 'Q':
        beta = Q
    elif alpha == 'E':
        beta = E 
    elif alpha == 'G':
        beta = G 
    elif alpha == 'H':
        beta = H
    elif alpha == 'I':
        beta = I
    elif alpha == 'L':
        beta = L
    elif alpha == 'K':
        beta = K
    elif alpha == 'M':
        beta = M
    elif alpha == 'F':
        beta = F
    elif alpha == 'P':
        beta = P
    elif alpha == 'S':
        beta = S
    elif alpha == 'T':
        beta = T
    elif alpha == 'W':
        beta = W 
    elif alpha == 'Y':
        beta = Y
    elif alpha == 'V':
        beta = V
    elif alpha == 'B':
        beta = gap
    elif alpha == 'Z':
        beta = gap
    elif alpha == 'X':
        beta = gap
    elif alpha == '-':
        beta = gap
    else:
        pass
    return beta
   
def mapMatrix(matrix):
    "Maps a matrix of floats"
    
    if matrix == "BLOSUM62":
        score_matrix = BLOSUM62
    elif matrix == "PAM250":
        score_matrix = PAM250
    elif matrix == "MCLACHLAN":
        score_matrix = MCLACHLAN
    elif matrix == "CPVN":
        score_matrix = CPVN
    elif matrix == "CLM":
        score_matrix = CLM
    else: pass
    
    return score_matrix

def mapRankMatrix(matrix):
    "Maps a matrix of floats"
    
    if matrix == "BLOSUM62":
        matrix = BLOSUM62
        rank_BLOSUM62 = []
        maximum = 0.0
        for i in matrix:
            for j in i:
                if j > maximum:
                    maximum = j
        for i in matrix:
            value = []
            for j in i:
                value.append(j * 1.0 / maximum)
            rank_BLOSUM62.append(value)
        score_matrix = rank_BLOSUM62
    elif matrix == "PAM250":
        matrix = PAM250
        rank_PAM250 = []
        maximum = 0.0
        for i in matrix:
            for j in i:
                if j > maximum:
                    maximum = j
        for i in matrix:
            value = []
            for j in i:
                value.append(j * 1.0 / maximum)
            rank_PAM250.append(value)
        score_matrix = rank_PAM250
    elif matrix == "MCLACHLAN":
        matrix = MCLACHLAN
        rank_MCLACHLAN = []
        maximum = 0.0
        for i in matrix:
            for j in i:
                if j > maximum:
                    maximum = j
        for i in matrix:
            value = []
            for j in i:
                value.append(j * 1.0 / maximum)
            rank_MCLACHLAN.append(value)
        score_matrix = rank_MCLACHLAN
    else: pass
    
    return score_matrix

def mapUserMatrix(filename):
    "TODO: Maps a matrix of floats from user input"
    
    inp = filename
    
    score_matrix = []
    input_matrix = open(inp, 'r')
    for line in input_matrix.readlines():
        score_matrix.append(map(float, line.split()))
    input_matrix.close()
    
    return score_matrix
   
def twoDimensionalMatrix(column, matrix, self_comparison, rank_order):
    """
    For each column in the alignment constructs a two-dimensional matrix.
    Very slow - use numpy
    """
    
    if rank_order is True:
        score_matrix = mapRankMatrix(matrix)
    else:
        score_matrix = mapMatrix(matrix) 
        
    two_d = []
    for i in range(len(column)):
        for j in range(len(column)):               
            if self_comparison is False:
                if i != j:
                    res1 = column[i]
                    res2 = column[j]
                    if res1 in aa and res2 in aa and res1 != '-' and res2 != '-':
                            if matrix is 'MCLACHLAN':
                                s = float(matchScoreMCLACHLAN(res1, res2, score_matrix))
                            else:
                                s = float(matchScore(res1, res2, score_matrix))
                            two_d.append(s)
                    else:
                        s = 0.0
                        two_d.append(s)
            else:
                res1 = column[i]
                res2 = column[j]
                if res1 in aa and res2 in aa and res1 != '-' and res2 != '-':
                        if matrix is 'MCLACHLAN':
                            s = float(matchScoreMCLACHLAN(res1, res2, score_matrix))
                        else:
                            s = float(matchScore(res1, res2, score_matrix))
                        two_d.append(s)
                else:
                    s = 0.0
                    two_d.append(s)
                     
    return two_d

def transpose(L):
    "Transposes row and columns"
    
    R = range(len(L[0]))
    rL = list()
    for i in R:
        rL.append(''.join([item[i] for item in L]))
    return rL
    

def probabilityDict(columns):
    "Caches character probabilities for each column"
    
    n = len(columns[0])
    pD = list()
    for col in columns:
        aa_list = list(set(col))
        values = [col.count(k) * 1.0 / n for k in aa_list]
        pD.append(dict(zip(aa_list, values)))
    return pD

def probabilityDictBackDist(columns):
    "Caches character probabilities for each column, with Background Distribution"
    
    score_matrix = mapMatrix('BLOSUM62')
    
    n = len(columns[0])
    pD = list()
    for col in columns:
        aa_list = list(set(col))
        values = [] 
        pim_total = 0.0
        qim_total = 0.0
        for x in aa_list:
            pim = aa_list.count(x) * 1.0 / n
            qim = float(matchScore(x, x, score_matrix))
            pim_total += pim
            qim_total += qim
        for x in aa_list:
            pix = aa_list.count(x) * 1.0 / n
            qix = float(matchScore(x, x, score_matrix))
            val = (pix * 1.0 / qix) / (pim_total * 1.0 / qim_total)
            values.append(val)
        pD.append(dict(zip(aa_list, values)))
    return pD

def probabilityDictPhysicProp(columns):
    "Caches character probabilities for each column, with Physicochemical Properties"
    
    n = len(columns[0])
    pD = list()
    for col in columns:
        aa_list = list(set(col))
        values = []
        for x in aa_list:
            pix = 0.0
            group = matchAA(x)
            for beta in group:
                pix += aa_list.count(beta) * 1.0 / n
            values.append(pix)
        pD.append(dict(zip(aa_list, values)))
    print pD
    return pD

def probabilityDictBackPhysic(columns):
    """Caches character probabilities for each column, with Background Distribution and
    Physicochemical Properties."""
    
    score_matrix = mapMatrix('BLOSUM62')
    
    n = len(columns[0])
    pD = list()
    for col in columns:
        values = [] 
        aa_list = list(set(col))
        pim_total = 0.0
        qim_total = 0.0
        for x in aa_list:
            pim = 0.0
            qim = 0.0
            group = matchAA(x)
            for aa in group:
                pim += aa_list.count(aa) * 1.0 / n
                qim += float(matchScore(x, x, score_matrix))
            pim_total += pim
            qim_total += qim
        for x in aa_list:
            pix = 0.0
            qix = 0.0
            group = matchAA(x)
            for aa in group:
                pix += aa_list.count(aa) * 1.0 / n
                qix += float(matchScore(x, x, score_matrix))
            val = (pix * 1.0 / qix) / (pim_total * 1.0 / qim_total)
            values.append(val)
        pD.append(dict(zip(aa_list, values)))
    return pD

def subAlignment (column, columns):
    "Creates a sub_alignment based on the most frequent AA in column i"
    
    pD = []
    i = column
    for j in range(len(i)):
        if i[j] in aa:
            freq = i.count(i[j])
            freq_aa = i[j]
            value = [freq_aa, freq]
            pD.append(value)
    
    sort = sorted(pD, key=lambda pD: pD[1])
    aa_x = sort[0][0]
    
    col_positions = []
    pos = -1
    for j in i:
        pos += 1
        if j == aa_x:
            col_positions.append(pos)
    
    sub_align = []
    for col in columns:
        sub_col = []
        for pos in col_positions:
            sub_col.append(col[pos])
        sub_align.append(sub_col)
    return sub_align

def subAlignment2 (column1, column2, columns):
    "Creates a sub_alignment based on AA identity of column i"
    
    i = column1
    j = column2
    
    list_i = []
    for x in i:
        if x in aa and x != '-':
            if x not in list_i:
                list_i.append(x)
    
    col_positions = []
    pos = -1
    for x in j:
        pos += 1
        if x in list_i:
            col_positions.append(pos)
    
    sub_align = []
    for col in columns:
        sub_col = []
        for pos in col_positions:
            sub_col.append(col[pos])
        sub_align.append(sub_col)
    return sub_align


def shannonEntropy(column):
    """    
    Sequence entropy - Shannon's informion theoretic entropy
    Shannon CE. A mathematical theory of communication. The Bell
    System Technical J 1948. 27:379-423, 623-656.
    low entropy => high conservation.
    """
    
    stList = list(column)
    alphabet = list(Set(stList)) 
    
    # calculate the frequency of each aa
    freqList = []
    for symbol in alphabet:
        ctr = 0
        for sym in stList:
            if sym == symbol:
                ctr += 1
        freqList.append(float(ctr) / len(stList))
    
    # Shannon entropy
    entropy = 0.0
    for freq in freqList:
        entropy = entropy + freq * log(freq, 2)
    H = -entropy
    
    return H

#-----------------------------------------------------------------------            
def _test():
    """Run the pycma module's doctests.

    This will try and locate the unit tests directory, and run the doctests
    from there in order that the relative paths used in the examples work.
    
    TODO: Finish this in the End
    """
    import doctest
    if os.path.isdir(os.path.join("..", "Tests")):
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join("..", "Tests"))
        doctest.testmod()
        os.chdir(cur_dir)
        del cur_dir
        print "Done"
    elif os.path.isdir(os.path.join("Tests")):
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join("Tests"))
        doctest.testmod()
        os.chdir(cur_dir)
        del cur_dir
        print "Done"

if __name__ == "__main__":
    #Run the doctests
    #_test()
    
    # print usage message
    Alignment()
    sys.exit()
    
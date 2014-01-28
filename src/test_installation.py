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

#-----------------------------------------------------------------------
import sys, pycma, version
from pycma import *
#-----------------------------------------------------------------------

TESTFILES = ['../Tests/example1.fasta', '../Tests/example2.fasta']

def module_info():
    print('pycma module information:')
    try:
        print('  * Version {0}'.format(version.__version__))
    except:
        try:
            print ('  * Version {0}'.format(pycma.__version__))
        except:
            pass
        pass
    print('  * Installation path {0}\n'.format(pycma.__file__))
    return

def alignment_info():
    align = Alignment(TESTFILES[0], TESTFILES[1], "fasta")
    align.number_sequences()
    align.length_sequences()
    print "Test succeeded"
    return

def main():
    module_info()
    alignment_info()
    return

if (__name__ == '__main__'):
    main()

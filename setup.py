#!/usr/bin/env python
# -*- coding: UTF-8 -*-                                                    
# Copyright (c) 2012 F. Madeira and L. Krippahl, 2012                                            
# This code is part of PyCMA distribution.                                 
# GNU General Public License - See LICENSE for more details.

#-----------------------------------------------------------------------

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
from distutils.core import setup
#-----------------------------------------------------------------------

setup(name='pycma',
      version='0.1.0',
      packages=['pycma'],
      package_dir={'pycma': 'src'},
      description='Inter and intra-protein coevolution analysis',
      long_description='PyCMA - a python module for correlated mutation analysis',
      author='Fábio Madeira & Ludwig Krippahl',
      author_email='fmadeira@campus.fct.unl.pt',
      url='https://github.com/fmadeira/pycma/',
      license='GNU General Public License (GPL)',
      platforms='any that supports python 2.7+',
      classifiers=[
           'Development Status :: Beta',
           'Environment :: Console',
           'Intended Audience :: Education',
           'Intended Audience :: Science/Research',
           'Intended Audience :: Developers',
           'License :: OSI Approved :: GNU General Public License (GPL)',
           'Operating System :: MacOS :: MacOS X',
           'Operating System :: Microsoft :: Windows',
           'Operating System :: POSIX',
           'Operating System :: POSIX :: SunOS/Solaris',
           'Operating System :: Unix'
           'Programming Language :: Python :: 2.7',
           'Topic :: Scientific/Engineering :: Bioinformatics',
           'Topic :: Scientific/Engineering :: Biology',
           'Topic :: Scientific/Engineering :: Computational Biology.',
           'Topic :: Software Development :: Libraries :: Python Modules'
              ],
      )
#-----------------------------------------------------------------------
#!/usr/bin/env python
# -*- coding: UTF-8 -*-                                                    
# Copyright (c) 2012 F. Madeira and L. Krippahl, 2012                                            
# This code is part of PyCMA distribution.                                 
# GNU General Public License - See LICENSE for more details.

#-----------------------------------------------------------------------

"""
PyCMA - A Python Module for Correlated Mutation Analysis.
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

__author__ = "F. Madeira (fmadeira@campus.fct.unl.pt)"

__all__ = ["pycma", "utils", "version", "test_installation"]

# ensure the user is running the version of python we require
import sys

if not hasattr(sys, "version_info") or sys.version_info < (2,7):
    raise RuntimeError("pycoevol requires Python 2.7 or later.")
del sys

# ensure the user has Biopython installed
try: 
    import Bio
    del Bio
except ImportError:
    raise RuntimeError("pycoevol requires Biopython 1.5.7 or later.")

# ensure the user has Biopython installed
try: 
    import numpy
    del numpy
except ImportError:
    raise RuntimeError("pycoevol requires Numpy 1.5.1 or later.")

# ensure the user has MatplotLib installed
try: 
    import matplotlib
    del matplotlib
except ImportError:
    raise RuntimeError("pycoevol requires MatplotLib 1.0.1 or later.")

# imports
from pycma import *

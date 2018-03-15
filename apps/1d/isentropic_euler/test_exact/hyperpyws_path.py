# coding: utf8

#==============================================================================#
# This file is part of HYPERPYWS: Hyperbolic Python WENO Solver
#
#   *** This software is made available "as is" without any assurance that it
#   *** will work for your purposes.  The software may in fact have defects, so
#   *** use the software at your own risk.
#
# License: GPL, see COPYING for details
#
# Copyright (C) 2013 
#
#    David Seal,  seal@math.msu.edu,  Michigan State University
#    Yaman Guclu, guclu@math.msu.edu, Michigan State University
#
#===============================================================================


"""
This module searches for the main library directory LIB_NAME by going up 
MAX_DEPTH levels in the directory tree, starting from this file location.

  * If the required directory is found, its path is added to sys.path
  * If the required directory is not found, the program exits with an error

Importing this module permits the script applications to use the parent library
with neither the need of installing the library itself on the local machine, nor 
the need of adding an environmental variable at runtime.

Usage
-----
>> import <library>_path

Required modules
----------------
  * Built-in: os, sys 

"""
#
# Author: Yaman Güçlü, December 2012 - Michigan State University
#
# Last revision: 25 Mar 2013
#

__all__ = []
__docformat__ = 'reStructuredText'

#===============================================================================

def ImportLibraryPath (LIB_NAME, MAX_DEPTH):
  """
  Add the parent library path to sys.path, so that the calling script can 
  import all relevant library modules even if the library is not properly 
  installed, and no specific environmental variables are set.
  
  Parameters
  ----------
  LIB_NAME : str
    Name of the parent library.
  
  MAX_DEPTH : int
    Maximum number of levels to be traversed in the directory tree.
  
  """
  assert (isinstance(LIB_NAME ,str))
  assert (isinstance(MAX_DEPTH,int))
  import os, sys
  
  lib_found = False
  
  # Determine directory from which program is called, and where this file is
  call_dir = os.path.abspath(os.path.curdir)
  file_dir = '/home/scott/bitbucket/hyperpyws/apps/advection'#os.path.dirname(os.path.abspath(__file__))
  
  # Look for the library by searching recursively in the parent directory
  os.chdir(file_dir)
  for i in range(MAX_DEPTH):
    if os.path.isdir(LIB_NAME):
      sys.path.append(os.path.abspath(os.path.curdir))
      lib_found = True
      break
    else:
      os.chdir(os.path.pardir)
  os.chdir(call_dir)
  
  # Stop execution with error if library search failed
  if not lib_found:
    sys.exit('Error: could not find library directory.')

#===============================================================================

if __file__.endswith('_path.py') or __file__.endswith('_path.pyc'):
  import os.path
  lib_name = os.path.split(__file__)[1].rpartition('_')[0]
  ImportLibraryPath (lib_name, 10)
else:
  import sys
  sys.exit('Error: file name is not in the form <library>_path.py')

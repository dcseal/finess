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

import math

#===============================================================================
# FUNCTION: Load data
#===============================================================================

def LoadData( *file_names ):
  """
  Load data from files.
  
  Returns
  -------
  mx : list of int
    Number of cells in domain for each of the simulations.
  
  data : list of dict
    For each file, dictionary of error norms for convergence analysis.
  
  """
  # Import numpy function for loading .dat files, and import ordered dictionary
  from numpy       import loadtxt
  from collections import OrderedDict
  
  # Load data and store it in a list of dictionaries
  data = []
  for fn in file_names:
    [mx, L1, L2, Li] = loadtxt( fn, unpack=True )
    record = OrderedDict()
    record['mx'] = mx.astype(int).tolist()
    record['L1'] = L1.tolist()
    record['L2'] = L2.tolist()
    record['Li'] = Li.tolist()
    data.append( record )

  # Check compatibility
  mx = data[0]['mx']
  for d in data:
    assert( d['mx'] == mx )
  
  return [mx, data]

#===============================================================================
# FUNCTION: Compute convergence ratios
#===============================================================================

def ComputeConvergence( mx, err ):
  """
  Compute convergence ratios.
  
  Parameters
  ----------
  mx : list of int
    Number of cells in domain for each of the simulations.
  
  err : list of floats
    Error norm for each of the simulations.
  
  """
  order = [None]
  for i in range(1,len(mx)):
    LogRatio_err = math.log(        err[i] /      err[i-1] )
    LogRatio_mx  = math.log( float(mx[i-1]) / float(mx[i]) )
    order.append( LogRatio_err / LogRatio_mx )
  
  return order

#===============================================================================
# FUNCTION: Create LaTeX table
#===============================================================================

def CreateLatexTable( mx, data,
                      file_name ='table.tex', 
                      mx_digits = 4,
                     err_digits = 2,
                     ord_digits = 2
                    ):
  
  # Import HyperPyWS library
# try               :  import hyperpyws
# except ImportError:  import hyperpyws_path
  
  # Number of runs and number of datasets
  NR = len(mx)
  ND = len(data)
  
  # Templates
  mx_fmt  = '{{:{0}d}}' .format(  mx_digits )
  err_fmt = '{:s}'
  ord_fmt = '{{:.{0}f}}'.format( ord_digits )
  
  row_template0 = ('${0}$' + ' & ${1}$ & ---' * ND + r'\\' + '\n') \
                  .format( mx_fmt, err_fmt )
  row_template  = ('${0}$' + ' & ${1}$ & ${2}$' * ND + r'\\' + '\n') \
                  .format( mx_fmt, err_fmt, ord_fmt )
  
  # Precomputed lines
  hline  = r'\hline'+'\n'
  begin  = r'\begin{tabular}{|r|'+'|c|c|'*ND+'}\n'
  end    = r'\end{tabular}'+'\n'
  header = r'\bf{Mesh} & ' + \
    ' & '.join( \
    [r'\bf{{M{} error}} & \bf{{Order}}'.format(i) for i in range(ND)]) + \
     r'\\' + '\n'
  
  # Exponential notation converter
  from output_utilities import Float2LatexConverter
  convert = Float2LatexConverter( err_digits )
  
  # Open output file
  f = open( file_name, 'w' )

  # Write table to file
  f.write( begin )
  f.write( hline )
  f.write( header)
  f.write( hline )
  f.write( hline )
  f.write( row_template0.format( mx[0], *[ convert(d['L1'][0]) for d in data ] ))
  f.write( hline )
  
  for i in range(1,NR):
    
    row_data = [ mx[i] ]
    for d in data:
      row_data += [ convert(d['L1'][i]), d['L1_order'][i] ]
    
    f.write( row_template.format( *row_data ) )
    f.write( hline )

  f.write( r'\end{tabular}' )
  f.write( '\n' )
  
  # Close output file
  f.close()

#===============================================================================
# FUNCTION: Parse input arguments
#===============================================================================

def parse_input():
  
  import argparse, sys
  
  parser = argparse.ArgumentParser (
      prog = 'python '+ sys.argv[0],
      description = '''Load error data from .dat files, estimate
                       convergence rates, and produce LaTeX table.''',
      formatter_class = argparse.ArgumentDefaultsHelpFormatter,
      )
  
  parser.add_argument(metavar = 'FILE',
                      nargs   = '*',
                      dest    = 'files',
                      help    = 'input files')
  
  parser.add_argument('-o', '--output',
                      default = 'table.tex',
                      help    = 'output file name')
  
  parser.add_argument('-d', '--digits',
                      metavar = ('MX','ERR','ORDER'),
                      type    =  int,
                      nargs   =  3,
                      default = [4,2,2],
                      help    = "number of digits for 'mx', 'err', and 'order'")
  
  return parser.parse_args()

#===============================================================================
# FUNCTION: Main script
#===============================================================================

def main():
  
  # Parse input arguments
  args = parse_input()
  print(args)
  print('')
  
  # Load data
  mx, data = LoadData( *args.files )
  
  # Compute convergence ratios
  for d in data:
    d['L1_order'] = ComputeConvergence( mx, d['L1'] )
  
  # Create LaTeX table
  CreateLatexTable( mx, data,
                    file_name = args.output,
                    mx_digits = args.digits[0],
                   err_digits = args.digits[1],
                   ord_digits = args.digits[2])

#===============================================================================
if __name__ == '__main__':
  #Run as main program
  main()

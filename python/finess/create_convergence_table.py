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
# data = []
# for fn in file_names:
#   [mx, L1, L2, Li] = loadtxt( fn, unpack=True )
#   record = OrderedDict()
#   record['mx'] = mx.astype(int).tolist()
#   record['L1'] = L1.tolist()
#   record['L2'] = L2.tolist()
#   record['Li'] = Li.tolist()
#   data.append( record )

  data = []
  for fn in file_names:
    [mx, Lerr] = loadtxt( fn, unpack=True )
    record = OrderedDict()
    record['mx'  ] = mx.astype(int).tolist()
    record['Lerr'] = Lerr.tolist()
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
# CLASS: Convert floats to friendly LaTeX format (This file was formerly found
# in an executable titled 'output_utilities.py'.  The problem with this name is
# that output* files are not only ignored in the git repo, but they are often
# deleted because this is where temporary output data gets dumped - this is the
# only file that is needed from that document.)
#===============================================================================


class Float2LatexConverter( object ):
  """
  Converter from Python number to string in LaTeX exponential notation.
  
  Number is first converted to a string in IEEE 754 floating-point format, 
  and then to a string in LaTeX-friendly exponential notation.
  
  Parameters
  ----------
  digits : int
    Number of digits to be kept in mantissa.
  
  Examples
  --------
  ::
    >>> convert = Float2LatexConverter( 5 )
    >>> convert( 1.25e-7 )
    '1.25000\\times 10^{-07}'
    >>> convert.set_prec( 1 )
    >>> convert( 1.25e-7 )
    '1.2\\times 10^{-07}'
  
  """
  def __init__( self, digits ):
    self.set_prec( digits )
  
  #-----------------------------------------------------------------------------
  def set_prec( self, digits ):
    """ Change number of digits in mantissa. """
    
    self._prec = digits
    self._str  = '{{:.{:d}e}}'.format( digits )
  
  #-----------------------------------------------------------------------------
  def __call__( self, x ):
    """ Convert number to LaTeX exponential notation. """
     
    float_str = self._str.format( x )
    return r'{0}\times 10^{{{1}}}'.format( *float_str.split('e') )


#===============================================================================
# FUNCTION: Create LaTeX table
#===============================================================================

def CreateLatexTable( mx, data,
                      file_name ='table.tex', 
                      mx_digits = 4,
                     err_digits = 2,
                     ord_digits = 2
                    ):
  
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
  convert = Float2LatexConverter( err_digits )
  
  # Open output file
  f = open( file_name, 'w' )

  f.write(r'\documentclass[draft]{amsart}'); f.write('\n');
  f.write('\n');
  f.write(r'\topmargin 0.0in'); f.write('\n'); 
  f.write(r'\headsep 0.2in'); f.write('\n');
  f.write(r'%\headheight  -2in%10.6pt'); f.write('\n');
  f.write(r'\oddsidemargin 0.0in '); f.write('\n');
  f.write(r'\textheight 8.75in '); f.write('\n');
  f.write(r'\textwidth 6.5in'); f.write('\n');
  f.write('\n');
  f.write(r'\usepackage{pslatex}'); f.write('\n');
  f.write(r'\usepackage{amsmath}'); f.write('\n');
  f.write(r'\usepackage{amssymb}'); f.write('\n');
  f.write(r'\usepackage{latexsym,pifont,color,comment}'); f.write('\n');
  f.write(r'\usepackage{stmaryrd}'); f.write('\n');
  f.write(r'\usepackage{esint}'); f.write('\n');
  f.write('\n');
  f.write(r'\pagestyle{plain}'); f.write('\n');
  f.write('\n');
  f.write(r'\begin{document}'); f.write('\n');
  f.write('\n');
  f.write(r'\begin{table}'); f.write('\n');
  f.write(r'\begin{center}'); f.write('\n');

  # ------------------- #
  # Write table to file
  # ------------------- #
  f.write( begin )
  f.write( hline )
  f.write( header)
  f.write( hline )
  f.write( hline )
  f.write( row_template0.format( mx[0], *[ convert(d['Lerr'][0]) for d in data ] ))
  f.write( hline )
  
  for i in range(1,NR):
    
    row_data = [ mx[i] ]
    for d in data:
      row_data += [ convert(d['Lerr'][i]), d['Order'][i] ]
    
    f.write( row_template.format( *row_data ) )
    f.write( hline )

  f.write(r'\end{tabular}' ); f.write('\n')
  f.write(r'\end{center}'); f.write('\n')
  f.write(r'\end{table}')

  # Close output file
  f.write('\n')
  f.write('\n');
  f.write('\end{document}');
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
    d['Order'] = ComputeConvergence( mx, d['Lerr'] )
  
  # Create LaTeX table
  CreateLatexTable( mx, data,
                    file_name = args.output,
                    mx_digits = args.digits[0],
                   err_digits = args.digits[1],
                   ord_digits = args.digits[2])

#===============================================================================
if __name__ == '__main__':
  'Run main program'

  main()

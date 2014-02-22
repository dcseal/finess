#!/usr/bin/env python
from __future__ import with_statement
from contextlib import closing
from subprocess import call, Popen, PIPE
import os
from math import log
import numpy as np

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
  from hyperpyws.output_utilities import Float2LatexConverter
  convert = Float2LatexConverter( err_digits )
  
  # Open output file
  f = open( file_name, 'w' )

  # Write table to file
  f.write( begin )
  f.write( hline )
  f.write( header)
  f.write( hline )
  f.write( hline )
  f.write( row_template0.format( mx[0], *[ convert(d['L2'][0]) for d in data ] ))
  f.write( hline )
  
  for i in range(1,NR):
    
    row_data = [ mx[i] ]
    for d in data:
      row_data += [ convert(d['L2'][i]), d['L2_order'][i] ]
    
    f.write( row_template.format( *row_data ) )
    f.write( hline )

  f.write( r'\end{tabular}' )
  f.write( '\n' )
  
  # Close output file
  f.close()


def main( ):
    '''Write some help documentation here
'''

    print "# leading comments can be given a '#' character"
    my_dictionary = {}
    old_err = i = 0
    while( 1 ):

        directory_num = my_dictionary['dir_num'] =  i

        folder = (os.getcwd() + '/output_%(dir_num)03i/') % my_dictionary

#       print folder
        if( not os.path.exists(folder) ):
            print 'did not find folder: %s' % folder
            break

        my_dictionary['curr_folder'] = folder
        directory_num = i

        try:
            qex  = np.loadtxt(folder + "/q0000.dat")[1:]
            qapp = np.loadtxt(folder + "/q0001.dat")[1:]
        except IOError:
            print('''Did not find the data file.
Please Wait for simulation to finish running.''')
            break
 
        diff = qex - qapp

        new_err = np.sqrt( np.dot(diff,diff) ) / np.sqrt( np.dot(qex,qex) )

        r1 = 'L2-error = %(new).3e;  ' % {'old': old_err, 'new' : new_err}
        
        if( old_err > 0 and new_err > 0 ):
            result = r1 + '   log2(ratio) = %(rat).3f' % \
                {'rat' : log( (old_err/new_err), 2) } 
        else:
            result = r1 + '   log2(ratio) = %(rat).3f' % \
                {'old' : old_err, 'new' : new_err, 'rat' : (old_err/new_err) } 

        # Live user feedback
        print result

        old_err = new_err
        i       = i + 1

# This is exactly the format I want (for my error table):
        err_digits = 2
        convert = Float2LatexConverter( err_digits )

#{\normalsize $25$}   &   {\normalsize $1.747\times 10^{-4}$}   & {\normalsize --} &  {\normalsize $8.292\times 10^{-5}$} & {\normalsize --}  \\

if __name__ == '__main__':
    import optparse
    parser = optparse.OptionParser(
        usage='''%%prog (-h |
    
%s''' % main.__doc__)

    opts, args = parser.parse_args()

    main( )

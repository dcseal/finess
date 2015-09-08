# The purpose of this module is to contain common IO (Input/Output) routines
# used for plotting unkowns in FINESS.
#
# Currently only supports 1D plotting.  
# See $(FINESS)/lib/[1-3]d/Output.cpp to see the order in which variables are saved to file.

#----------------------------------------------------------
def read_qfile(mx, meqn, qfile ):
    """Read solution from file.  The format of this file is identical to that
    used by FINESS for output.
    """

    import string
    import numpy as np

    # solution to be returned
    qsoln = np.zeros((mx, meqn), float);

    # open file
    Rqfile = open(qfile,'r')

    # get time
    linestring = Rqfile.readline()
    linelist   = string.split(linestring)
    time       = float(linelist[0])
    
    # extract all point values
    k = 0
    for me in range(meqn):
        for i in range(mx):

            linestring  = Rqfile.readline()
            linelist    = string.split(linestring)
            qsoln[i,me] = float(linelist[0])

            k = k+1

    # close file
    Rqfile.close()

    # return time
    return time, qsoln
#----------------------------------------------------------


def parse_ini_parameters(parameters_file, params={} ):
    """Parse a parameters.ini file.

    This function parses an input file parameters.ini and returns a dictionary
    containing key-value pairs for each required object in a DoGPack
    parameters.ini file.

    TODO - test for errors on checking these files - read in from the default
    parameters file if that object is missing.

    Input
    -----

        parameters_file - a string pointing to a valid parameters.ini file.

    Returns
    -------

        params - a dictionary containing key-value pairs.  If params already
        exists as a dictionary, this routine will add new keys.

    TODO - this should be rewritten to be consistent with the parameters that
    $FINESS/python/finess/params defines.
    """

    import ConfigParser

    config = ConfigParser.RawConfigParser()
    config.read( parameters_file )

    # TODO - make sure this is consistent with the code generation software as
    # well.  There should be a clear way to do this.  -DS 9/8/2015

    params['ndims']       = config.getint  ('finess', 'ndims'       )
    params['nout']        = config.getint  ('finess', 'nout'        )
    params['tfinal']      = config.getfloat('finess', 'tfinal'      )
    params['dt_init']     = config.getfloat('finess', 'initial_dt'  )
    params['dt_max']      = config.getfloat('finess', 'max_dt'      )
    params['cfl_max']     = config.getfloat('finess', 'max_cfl'     )
    params['cfl_des']     = config.getfloat('finess', 'desired_cfl' )
    params['nv']          = config.getint  ('finess', 'nv'          )

    params['time_stepping_method'] = config.get('finess', 'time_stepping_method' )
    params['space_order'] = config.getint     ('finess', 'space_order' )
    params['time_order']  = config.getint     ('finess', 'time_order'  )
    params['mcapa']       = config.getint     ('finess', 'mcapa'       )
    params['maux']        = config.getint     ('finess', 'maux'        )
    params['source_term'] = config.getboolean ('finess', 'source_term' )
    params['meqn']        = config.getint     ('finess', 'meqn'        )

    if( config.has_option('finess', 'output_dir' ) ):
        params['output_dir']  = config.get        ('finess', 'output_dir' )
    else:
        params['output_dir']  = 'output'

    # Parse grid data  (Every solver at least has this portion!)
    params['mx' ]         = config.getint  ('grid', 'mx'    )
    params['mbc']         = config.getint  ('grid', 'mbc'   )
    params['xlow' ]       = config.getfloat('grid', 'xlow'  )
    params['xhigh']       = config.getfloat('grid', 'xhigh' )

    # 2D "grid" data
    if( params['ndims'] > 1 ):
        params['my' ]         = config.getint  ('grid', 'my'    )
        params['ylow' ]       = config.getfloat('grid', 'ylow'  )
        params['yhigh']       = config.getfloat('grid', 'yhigh' )

    # 3D "grid" data
    if( params['ndims'] > 2 ):
        params['mz' ]         = config.getint  ('grid', 'mz'    )
        params['zlow' ]       = config.getfloat('grid', 'zlow'  )
        params['zhigh']       = config.getfloat('grid', 'zhigh' )

    return params


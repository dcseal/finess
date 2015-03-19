"""
This module provides functions for generating coordinates arrays. 
See also: TODO.
"""

from __future__ import absolute_import


def meshgrid(params):
    """Returns meshgrid (a pair (X, Y)) that can be used for 2D plotting.
    params is what is returned by finess.params.util.read_params."""
    assert(params['finess', 'ndims'] == 3)
    mx    = params['grid', 'mx']
    my    = params['grid', 'my']
    mz    = params['grid', 'mz']
    xlow  = params['grid', 'xlow']
    xhigh = params['grid', 'xhigh']
    ylow  = params['grid', 'ylow']
    yhigh = params['grid', 'yhigh']
    zlow  = params['grid', 'zlow']
    zhigh = params['grid', 'zhigh']
    
    dx = (xhigh-xlow) / float(mx)
    dy = (yhigh-ylow) / float(my)
    dz = (zhigh-zlow) / float(mz)
    
    from pylab import meshgrid, linspace
    X, Y, Z = meshgrid(linspace(xlow + 0.5*dx, xhigh - 0.5*dx, mx),
                       linspace(ylow + 0.5*dy, yhigh - 0.5*dy, my), 
                       linspace(zlow + 0.5*dz, zhigh - 0.5*dz, mz),
                       indexing='ij' )
    return X, Y, Z


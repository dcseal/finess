from __future__ import absolute_import

def meshgrid(params):
    """Returns meshgrid X that can be used for 1D plotting.
    params is what is returned by finess.params.util.read_params."""
    assert(params['finess', 'ndims'] == 1)
    mx    = params['grid', 'mx']
    xlow  = params['grid', 'xlow']
    xhigh = params['grid', 'xhigh']
    
    dx = (xhigh-xlow) / float(mx)
    
    from pylab import linspace
    X = linspace(xlow + 0.5*dx, xhigh - 0.5*dx, mx)
    return X

def ask_which_component_and_which_frame_and_plot_red_dots(parameters_ini_filename):
    """Call this and you get an interactive plot program as the
    function name suggests."""
    from finess.params.util import read_params
    from generate_iniparams import parameter_list
    from finess.viz import ask_which_component, interactive_plot
    
   
    params = read_params(parameters_ini_filename, parameter_list)
   
    component = ask_which_component(params)
    def draw_ith_frame(fig, i):
        import matplotlib.figure 
        assert isinstance(fig, matplotlib.figure.Figure)
        from finess.dim1 import read_qa
        
        import finess.viz.dim1
        X = finess.viz.dim1.meshgrid(params)

        t, q, aux = read_qa(params, i)

        ax = fig.add_subplot(111)
        ax.set_title("$t = %(t)f$" % {"t": t})

        ax.plot(X, q[:, component - 1], 'r.')

   
    interactive_plot(params["finess", "nout"], 
                     draw_ith_frame = draw_ith_frame)
 

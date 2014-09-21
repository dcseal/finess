from __future__ import absolute_import

def meshgrid(params):
    """Returns meshgrid (a pair (X, Y)) that can be used for plotting."""
    assert params['finess', 'ndims'] == 2
    mx = params['grid', 'mx']
    my = params['grid', 'my']
    xlow = params['grid', 'xlow']
    xhigh = params['grid', 'xhigh']
    ylow = params['grid', 'ylow']
    yhigh = params['grid', 'yhigh']
    
    dx = (xhigh-xlow) / float(mx)
    dy = (yhigh-ylow) / float(my)
    
    from pylab import meshgrid, linspace
    X, Y = meshgrid(linspace(xlow + 0.5*dx, xhigh - 0.5*dx, mx), linspace(ylow + 0.5*dy, yhigh - 0.5*dy, my), indexing = "ij")
    return X, Y



def draw_ith_frame_jth_component(params, fig, i, j, 
                                 plotting_method_on_Axes3D):
    from finess.dim2 import read_qa
    from finess.params.util import read_params

    output_dir = params["finess", "output_dir"]
    t, q, aux = read_qa(params, i, output_dir = output_dir)

    from mpl_toolkits.mplot3d import Axes3D
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title("$t = %(t)f$" % {"t": t})

    import finess.viz.dim2
    X, Y = finess.viz.dim2.meshgrid(params)
    plotting_method_on_Axes3D(ax, X, Y, q[:, :, j - 1])
 
def ask_which_component_and_which_frame_and_plot_wireframe():
    from finess.params.util import read_params
    from generate_iniparams import parameter_list
    from finess.viz import ask_which_component, interactive_plot
    from mpl_toolkits.mplot3d import Axes3D
    from finess.viz.dim2 import draw_ith_frame_jth_component
    
    params = read_params("parameters.ini", parameter_list)
    
    component = ask_which_component(params)
    
    draw_ith_frame = \
        lambda fig, i: \
           draw_ith_frame_jth_component( \
             params = params, fig = fig, i = i, j = component,
             plotting_method_on_Axes3D = Axes3D.plot_wireframe)
    
    interactive_plot(params["finess", "nout"], 
                     draw_ith_frame = draw_ith_frame)
    

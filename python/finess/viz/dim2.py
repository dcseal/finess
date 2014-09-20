from __future__ import absolute_import

def draw_ith_frame_jth_component(params, fig, i, j, plot_on_Axes3D):
    from finess.dim2 import read_qa
    from finess.params.util import read_params

    if ("finess", "output_dir") in params and \
       params["finess", "output_dir"] != "":
        output_dir = params["finess", "output_dir"]
    else:
        output_dir = "output"

    t, q, aux = read_qa(params, i, output_dir = output_dir)

    from mpl_toolkits.mplot3d import Axes3D
    ax = fig.add_subplot(111, projection='3d')

    import finess.dim2
    X, Y = finess.dim2.meshgrid(params)
    plot_on_Axes3D(ax, X, Y, q[:, :, j - 1])
 
    

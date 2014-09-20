
import finess.dim2
import finess.params.util

from generate_iniparams import parameter_list

params = finess.params.util.read_params("parameters.ini", parameter_list)



from finess.viz import interactive_plot


from mpl_toolkits.mplot3d import Axes3D

from finess.viz.dim2 import draw_ith_frame_jth_component

interactive_plot(params["finess", "nout"], \
   lambda fig, i: 
       draw_ith_frame_jth_component( \
           params = params, fig = fig, i = i, j = 1,
           plot_on_Axes3D = Axes3D.plot_wireframe))






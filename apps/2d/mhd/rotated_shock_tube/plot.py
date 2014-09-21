
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


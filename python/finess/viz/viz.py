"""Provides basic plotting framework.

All names in this module can be imported via finess.viz.
"""

def interactive_plot(max_i, draw_ith_frame, new_window_for_each_frame = False):
    """Run a loop asking for which frame to plot, and draw that frame
    using draw_ith_frame.
    draw_ith_frame must be callable in this form:
    draw_ith_frame(fig, i)
    where fig is an instance of matplotlib.figure.Figure.
    """
    def get_next_i(current_i):
        """Returns -1 if want to quit."""
        print "Plot which frame ( 0 - %(max_i)d ) [type anything else to quit] ? " % {"max_i": max_i}, 
        i_str = raw_input()
        
        if i_str == "":
            if current_i < max_i:
                return current_i + 1
            else:
                print "End of plots"
                return max_i
            
        try:
            i_entered = int(i_str)
        except ValueError:
            i_entered = -1
        
        if i_entered > max_i:
            print "End of plots"
            return max_i
        else:
            return i_entered
    
    import matplotlib.pyplot as pyplot    
    i = -1
    i = get_next_i(i)
    pyplot.ion()
    if i != -1:
        if not new_window_for_each_frame:
            while i != -1:
                fig = pyplot.figure(1)
                fig.clf()
                draw_ith_frame(fig, i)
                i = get_next_i(i)
                fig.canvas.draw()
        else:
            while i != -1:
                fig = pyplot.figure()
                draw_ith_frame(fig, i)
                i = get_next_i(i)

    pyplot.ioff()            

def ask_which_component(params):
    """Ask for which component to plot.
    params is what is returned from finess.params.util.read_params."""
    def raise_error(component_str, meqn):
        raise RuntimeError( \
          "Invalid input: %(component_str)s is not one of %(valid_range)s"% \
          {"component_str": component_str,
           "valid_range": str(range(1, meqn + 1))})
    meqn = params["finess", "meqn"]

    print 'Plot which component of q  ( 1 - %(meqn)d ) ? ' \
          % {"meqn": meqn},
    component_str = raw_input()

    try:
        component = int(component_str)
    except ValueError:
        raise_error(component_str, meqn)

    if component >= 1 and component <= meqn:
        return component
    else:
        raise_error(component_str, meqn)
        
    

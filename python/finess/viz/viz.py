def interactive_plot(max_i, draw_ith_frame, new_window_for_each_frame = False):
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
 

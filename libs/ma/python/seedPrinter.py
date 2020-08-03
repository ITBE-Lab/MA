from MS import *
from MA import *
from bokeh.plotting import figure, show

class SeedPrinter(Module):
    def __init__(self, parameter_manager, name_a="Data", name_b="Ground Truth", x_range=None, y_range=None, 
                 do_print=True):
        self.name_a = name_a
        self.name_b = name_b
        self.x_range = x_range
        self.y_range = y_range
        self.do_print = do_print

    ##
    # @brief Execute LineSweep for all given seeds.
    # @details
    # Reimplemented from MA.aligner.Module.execute.
    def execute(self, *input):
        assert(len(input) >= 1)

        plot = figure(title="Seeds", x_range=self.x_range, y_range=self.y_range)
        if len(input) > 2:
            helper_ret = input[2]
            l = []
            r = []
            t = []
            b = []
            for rect in helper_ret.rectangles:
                l.append(rect.x_axis.start)
                r.append(rect.x_axis.start+rect.x_axis.size)
                b.append(rect.y_axis.start)
                t.append(rect.y_axis.start+rect.y_axis.size)
            plot.quad(left=l, right=r, bottom=b, top=t, color="lightgrey", alpha=0.2)
            
            if self.do_print:
                for a, b in helper_ret.jump_seeds:
                    print("jump", a.start, a.start_ref, a.size, "forw" if a.on_forward_strand else "rev", "-",
                          b.start, b.start_ref, b.size, "forw" if b.on_forward_strand else "rev")

        def render(seeds, name, dash=(10,0), width=1):
            forw_x = []
            forw_y = []
            rev_x = []
            rev_y = []
            for seed in seeds:
                if self.do_print:
                    print(name, seed.start, seed.start_ref, seed.size, "forw" if seed.on_forward_strand else "rev")
                if seed.on_forward_strand:
                    forw_x.append(seed.start_ref)
                    forw_x.append(seed.start_ref + seed.size)
                    forw_x.append(float("NaN"))
                    forw_y.append(seed.start)
                    forw_y.append(seed.start + seed.size)
                    forw_y.append(float("NaN"))
                else:
                    rev_x.append(seed.start_ref)
                    rev_x.append(seed.start_ref - seed.size)
                    rev_x.append(float("NaN"))
                    rev_y.append(seed.start)
                    rev_y.append(seed.start + seed.size)
                    rev_y.append(float("NaN"))
            plot.line(x=forw_x, y=forw_y,
                        legend_label=name + " - forward", line_dash=dash, line_width=width)
            plot.line(x=rev_x, y=rev_y,
                        legend_label=name + " - reverse", line_color="orange", line_dash=dash,
                        line_width=width)

        if len(input) > 1 and not input[1] is None:
            seeds_b = input[1]
            render(seeds_b, self.name_b, dash=(10,10), width=4)

        seeds_a = input[0]
        render(seeds_a, self.name_a)

        show(plot)


        return None

class SeedPointPrinter(Module):
    def __init__(self, parameter_manager, name_a="Data", name_b="Ground Truth", x_range=None, y_range=None, 
                 do_print=True):
        self.name_a = name_a
        self.name_b = name_b
        self.x_range = x_range
        self.y_range = y_range
        self.do_print = do_print

    ##
    # @brief Execute LineSweep for all given seeds.
    # @details
    # Reimplemented from MA.aligner.Module.execute.
    def execute(self, *input):
        assert(len(input) == 2)

        plot = figure(title="Seeds", x_range=self.x_range, y_range=self.y_range)

        def render(seeds, name, dash=(10,0), width=1):
            forw_x = []
            forw_y = []
            rev_x = []
            rev_y = []
            for seed in seeds:
                if self.do_print:
                    print(name, seed.start, seed.start_ref, seed.size, "forw" if seed.on_forward_strand else "rev")
                if seed.on_forward_strand:
                    forw_x.append(seed.start_ref)
                    forw_x.append(seed.start_ref + seed.size)
                    forw_x.append(float("NaN"))
                    forw_y.append(seed.start)
                    forw_y.append(seed.start + seed.size)
                    forw_y.append(float("NaN"))
                else:
                    rev_x.append(seed.start_ref)
                    rev_x.append(seed.start_ref - seed.size)
                    rev_x.append(float("NaN"))
                    rev_y.append(seed.start)
                    rev_y.append(seed.start + seed.size)
                    rev_y.append(float("NaN"))
            plot.line(x=forw_x, y=forw_y,
                        legend_label=name + " - forward", line_dash=dash, line_width=width)
            plot.line(x=rev_x, y=rev_y,
                        legend_label=name + " - reverse", line_color="orange", line_dash=dash,
                        line_width=width)

        seeds_a = input[0]
        render(seeds_a, self.name_a)

        points = input[1]
        plot.x(x=[p[1] for p in points],y=[p[0] for p in points],legend_label=self.name_b,
               line_color=["blue" if p[2] else "orange" for p in points])

        show(plot)


        return None
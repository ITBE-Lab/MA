from bokeh.plotting import figure
from bokeh.models.tools import HoverTool

class SeedPlot:
    def __init__(self, main_plot):
        self.left_plot = figure(
            width=300,
            height=900,
            y_range=main_plot.plot.y_range,
            tools=["xpan", "xwheel_zoom"],
            active_scroll="xwheel_zoom",
            toolbar_location=None
        )
        self.left_plot.xaxis.major_label_orientation = math.pi/2
        self.left_plot.yaxis.axis_label = "Reference Position"
        self.left_plot.xaxis.axis_label = "Read Id"

        self.bottom_plot = figure(
            width=900,
            height=300,
            x_range=main_plot.plot.x_range,
            y_range=self.left_plot.x_range,
            tools=["ypan", "ywheel_zoom"],
            active_scroll="ywheel_zoom",
            toolbar_location=None
        )
        self.bottom_plot.xaxis.axis_label = "Reference Position"
        self.bottom_plot.yaxis.axis_label = "Read Id"

        # ambigious regions (red rectangles)
        self.ambiguous_regions = ColumnDataSource()
        self.bottom_plot.quad(left="l", bottom="b", right="r", top="t", fill_alpha=0.5,
                            fill_color="red", line_width=0, source=self.ambiguous_regions, name="ambiguous_regions")
        self.left_plot.quad(left="b", bottom="l", right="t", top="r", fill_alpha=0.5,
                            fill_color="red", line_width=0, source=self.ambiguous_regions, name="ambiguous_regions")
                            
        hover_ambiguous_regions = HoverTool(tooltips=[("left", "@l"),
                                                      ("bottom", "@b"),
                                                      ("right", "@r"),
                                                      ("top", "@t"),
                                                      ("fill percentage", "@f"),
                                                      ("additional seed size", "@s")],
                                            names=['ambiguous_regions'],
                                            name="Hover ambiguous regions")
        self.left_plot.add_tools(hover_ambiguous_regions)
        self.bottom_plot.add_tools(hover_ambiguous_regions)

        # seeds
        self.seeds = ColumnDataSource()
        self.left_plot.rect(x="category", y="center", width=1, height="size",
                            fill_color="c", line_width=0, source=self.seeds, name="seeds")
        self.bottom_plot.rect(y="category", x="center", height=1, width="size",
                              fill_color="c", line_width=0, source=self.seeds, name="seeds")

        hover_seeds = HoverTool(tooltips=[("read id", "@r_id"),
                                          ("q, r, l", "@q, @r, @l"),
                                          ("index", "@idx"),
                                          ("reseeding-layer", "@layer"),
                                          ("parlindrome-filtered", "@parlindrome")],
                                names=['seeds'],
                                name="Hover reads")
        self.left_plot.add_tools(self.hover_seeds)
        self.bottom_plot.add_tools(self.hover_seeds)
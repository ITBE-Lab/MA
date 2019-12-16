from bokeh.plotting import figure
from bokeh.models.tools import HoverTool
from bokeh.plotting import ColumnDataSource

class ReadPlotNucs:
    def __init__(self, nuc_plot, read_plot):
        self.left_plot = figure(
            width=40,
            height=900,
            y_range=read_plot.plot.y_range,
            tools=["ypan", "ywheel_zoom"],
            active_scroll="ywheel_zoom",
            toolbar_location=None
        )
        self.left_plot.xaxis.visible = False
        self.left_plot.grid.visible = False
        self.left_plot.yaxis.axis_label = "Read Position"

        self.bottom_plot = figure(
            width=900,
            height=40,
            x_range=read_plot.plot.x_range,
            tools=["xpan", "xwheel_zoom"],
            active_scroll="xwheel_zoom",
            toolbar_location=None
        )
        self.bottom_plot.yaxis.visible = False
        self.bottom_plot.grid.visible = False
        self.bottom_plot.xaxis.axis_label = "Reference Position"

        # the nucleotides from the read
        self.left_nucs = ColumnDataSource({"c":[], "x":[], "y":[]})
        self.left_plot.rect(x=0.5, y="y", width=1, height=1, fill_color="c", line_width=0,
                            source=self.left_nucs, name="nucleotides")

        # the nucleotides from the rendered region on the genome (left and bottom)
        self.bottom_plot.rect(x="x", y=0.5, width=1, height=1, fill_color="c", line_width=0,
                            source=nuc_plot.left_nucs, name="nucleotides")
        self.bottom_plot.rect(x="x", y=0.5, width=1, height=1, fill_color="c", line_width=0,
                            source=nuc_plot.bottom_nucs, name="nucleotides")

        hover_nucleotides = HoverTool(tooltips="@i", names=['nucleotides'], name="Hover nucleotides")
        self.left_plot.add_tools(hover_nucleotides)
        self.bottom_plot.add_tools(hover_nucleotides)

class ReadPlot:
    def __init__(self, nuc_plot):

        self.plot = figure(
            width=900,
            height=400,
            tools=[
                "pan", "box_zoom",
                "wheel_zoom", "reset"
            ],
            active_scroll="wheel_zoom"
        )
        self.plot.axis.visible = False
        self.plot.toolbar.logo = None

        # the ambiguity rectangles
        self.ambiguity_rect = ColumnDataSource({"l":[], "b":[], "r":[], "t":[], "c":[]})
        self.plot.quad(left="l", bottom="b", right="r", top="t", fill_color="c",
                        fill_alpha=0.2, line_width=0, name="ambiguity_rect",
                        source=self.ambiguity_rect)

        self.plot.add_tools(HoverTool(tooltips=[("left", "@l"),
                                                ("bottom", "@b"),
                                                ("right", "@r"),
                                                ("top", "@t"),
                                                ("fill percentage", "@f"),
                                                ("additional seed size", "@s")],
                                     names=['ambiguity_rect'],
                                     name="Hover ambiguity rects"))

        # the seeds
        self.seeds = ColumnDataSource({"x":[], "y":[], "c":[]})
        self.plot.multi_line(xs="x", ys="y", line_color="c", line_width=5, source=self.seeds, name="seeds")

        self.plot.add_tools(HoverTool(tooltips=[("read id", "@r_id"),
                                                ("q, r, l", "@q, @r, @l"),
                                                ("index", "@idx"),
                                                ("reseeding-layer", "@layer"),
                                                ("parlindrome-filtered", "@parlindrome")],
                                      names=['seeds'],
                                      name="Hover seeds"))

        self.nuc_plot = ReadPlotNucs(nuc_plot, self)


from bokeh.plotting import figure
from bokeh.models.tools import HoverTool
from bokeh.plotting import ColumnDataSource
from bokeh.models import BasicTickFormatter
from bokeh.palettes import Plasma256
from bokeh.events import Tap
import math

class ReadPlotNucs:
    def __init__(self, nuc_plot, read_plot):
        self.left_plot = figure(
            width=100,
            height=400,
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
            height=100,
            x_range=read_plot.plot.x_range,
            tools=["xpan", "xwheel_zoom"],
            active_scroll="xwheel_zoom",
            toolbar_location=None
        )
        self.bottom_plot.yaxis.visible = False
        self.bottom_plot.grid.visible = False
        self.bottom_plot.xaxis.axis_label = "Reference Position"
        self.bottom_plot.xaxis.major_label_orientation = math.pi/4

        # the nucleotides from the read
        self.nucs_by_r_id = {} # dict of {"p": [], "c": [], "i": []}
        self.left_nucs = ColumnDataSource({"c":[], "p":[]})
        self.left_plot.rect(x=0.5, y="p", width=1, height=1, fill_color="c", line_width=0,
                            source=self.left_nucs, name="nucleotides")

        # the nucleotides from the rendered region on the genome (left and bottom)
        self.bottom_plot.rect(x="p", y=0.5, width=1, height=1, fill_color="c", line_width=0,
                            source=nuc_plot.left_nucs, name="nucleotides")
        self.bottom_plot.rect(x="p", y=0.5, width=1, height=1, fill_color="c", line_width=0,
                            source=nuc_plot.bottom_nucs, name="nucleotides")

        hover_nucleotides = HoverTool(tooltips="@i", names=['nucleotides'], name="Hover nucleotides")
        self.left_plot.add_tools(hover_nucleotides)
        self.bottom_plot.add_tools(hover_nucleotides)

    def copy_nts(self, renderer):
        if not renderer.selected_read_id is None:
            self.left_nucs.data = self.nucs_by_r_id[renderer.selected_read_id]

    def reset_nts(self):
        self.left_nucs.data = {"c":[], "p":[]}

class ReadPlot:
    def __init__(self, nuc_plot, renderer):

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
                                                ("additional seed size", "@s"),
                                                ("triggered ambiguity filter", "@dp")
                                                ],
                                     names=['ambiguity_rect'],
                                     name="Hover ambiguity rects"))

        # the seeds
        self.seeds = ColumnDataSource({"category":[], "center":[], "size":[], "x":[], "y":[], "c":[]})
        self.plot.multi_line(xs="x", ys="y", line_color="c", line_width=5, source=self.seeds, name="seeds")

        self.plot.add_tools(HoverTool(tooltips=[("read id", "@r_id"),
                                                ("read name", "@r_name"),
                                                ("q, r, l", "@q, @r, @l"),
                                                ("index", "@idx"),
                                                ("reseeding-layer", "@layer"),
                                                ("parlindrome-filtered", "@parlindrome")],
                                      names=['seeds'],
                                      name="Hover seeds"))

        self.nuc_plot = ReadPlotNucs(nuc_plot, self)

        # auto adjust y-range of read plot
        renderer.main_plot.plot.x_range.on_change("start", lambda x,y,z: self.auto_adjust_y_range(renderer))

        # make seeds clickable
        self.plot.on_event(Tap, lambda tap: self.seed_tap(renderer, tap.x, tap.y))

    def auto_adjust_y_range(self, renderer):
        if renderer.widgets.range_link_radio.active == 0:
            self.plot.x_range.start = renderer.main_plot.plot.x_range.start
            self.plot.x_range.end = renderer.main_plot.plot.x_range.end
        else:
            self.plot.x_range.start = renderer.main_plot.plot.y_range.start
            self.plot.x_range.end = renderer.main_plot.plot.y_range.end

        min_seed = 10000000
        max_seed = 0
        def in_x_range(val):
            return self.plot.x_range.start <= val and val <= self.plot.x_range.end
        for idx, _ in enumerate(self.seeds.data["x"]):
            if in_x_range(self.seeds.data["x"][idx][0]) or in_x_range(self.seeds.data["x"][idx][1]):
                min_seed = min(min_seed, *self.seeds.data["y"][idx])
                max_seed = max(max_seed, *self.seeds.data["y"][idx])

        size = max_seed - min_seed
        self.plot.y_range.start = min_seed - size / 20
        self.plot.y_range.end = max_seed + size / 20

    def copy_seeds(self, renderer, condition):
        seed_dict = dict((key, []) for key in renderer.seed_plot.seeds.data.keys())
        found_at_least_one = False
        for idx, _ in enumerate(renderer.seed_plot.seeds.data["center"]):
            if condition(idx):
                for key in renderer.seed_plot.seeds.data.keys():
                    seed_dict[key].append(renderer.seed_plot.seeds.data[key][idx])
                found_at_least_one = True
        if found_at_least_one:
            self.seeds.data = seed_dict
            self.ambiguity_rect.data = renderer.read_plot_rects[renderer.selected_read_id]

    def reset_seeds(self, renderer):
        seed_dict = dict((key, []) for key in renderer.seed_plot.seeds.data.keys())
        self.seeds.data = seed_dict

    def seed_tap(self, renderer, x, y):
        renderer.selected_call_id = set()
        renderer.selected_jump_id = set()
        renderer.selected_seed_id = None
        mouse_delta = x-y
        mouse_add = x+y
        min_dist = 10
        for idx, _ in enumerate(self.seeds.data["center"]):
            if self.seeds.data["f"][idx]:
                seed_delta = self.seeds.data["x"][idx][0] - self.seeds.data["y"][idx][0]
                seed_add_1 = self.seeds.data["x"][idx][0] + self.seeds.data["y"][idx][0]
                seed_add_2 = self.seeds.data["x"][idx][1] + self.seeds.data["y"][idx][1]
                if seed_add_1 <= mouse_add and mouse_add <= seed_add_2 and abs(seed_delta-mouse_delta) <= min_dist:
                    min_dist = abs(seed_delta-mouse_delta)
                    renderer.selected_seed_id = (self.seeds.data["idx"][idx], self.seeds.data["r_id"][idx])
            else:
                seed_delta_1 = self.seeds.data["x"][idx][1] - self.seeds.data["y"][idx][1]
                seed_delta_2 = self.seeds.data["x"][idx][0] - self.seeds.data["y"][idx][0]
                seed_add = self.seeds.data["x"][idx][0] + self.seeds.data["y"][idx][0]
                if seed_delta_1 <= mouse_delta and mouse_delta <= seed_delta_2 and abs(seed_add-mouse_add) <= min_dist:
                    min_dist = abs(seed_add-mouse_add)
                    renderer.selected_seed_id = (self.seeds.data["idx"][idx], self.seeds.data["r_id"][idx])

        renderer.seed_plot.update_selection(renderer)
        self.copy_seeds(renderer, lambda idx: renderer.seed_plot.seeds.data["r_id"][idx] == renderer.selected_read_id)
        renderer.main_plot.update_selection(renderer)

    def reset_cds(self):
        self.ambiguity_rect.data = {"l":[], "b":[], "r":[], "t":[], "c":[]}
        self.seeds.data = {"category":[], "center":[], "size":[], "x":[], "y":[], "c":[]}
        #self.nuc_plot.reset_nts()


from bokeh.plotting import figure
from bokeh.models.tools import HoverTool
from bokeh.plotting import ColumnDataSource
import copy
from bokeh.events import Tap

class MainPlot:
    def __init__(self, renderer):
        self.plot = figure(
            width=900,
            height=900,
            tools=[
                "pan", "box_zoom",
                "wheel_zoom", "save"
            ],
            active_scroll="wheel_zoom"
        )
        self.plot.axis.visible = False

        # the diagonal line & genome outline rectangle
        self.diagonal_line = ColumnDataSource({"xs":[], "ys":[]})
        self.plot.line(x="xs", y="ys", line_color="black", line_width=3, source=self.diagonal_line)

        self.genome_outline = ColumnDataSource({"x":[], "y":[], "w":[], "h":[]})
        self.plot.quad(left="x", bottom="y", right="w", top="h",
                       fill_alpha=0, line_color="black", line_width=3,
                       source=self.genome_outline)

        # the sv jumps
        self.jump_quads = []
        for _ in range(4):
            self.jump_quads.append(ColumnDataSource({"x":[], "y":[], "w":[], "h":[], "c":[], "a":[]}))
            self.plot.quad(left="x", bottom="y", right="w", top="h", fill_color="c",
                           line_color="c", line_width=3, fill_alpha="a",
                           source=self.jump_quads[-1], name="jump_quads")

        self.jump_x = ColumnDataSource({"x":[], "y":[]})
        self.plot.multi_line(xs="x", ys="y", line_width=1.5, line_alpha=0.5, color="black", source=self.jump_x)

        self.plot.add_tools(HoverTool(tooltips=[("supp. nt", "@n"),
                                                ("read id", "@r"),
                                                ("|query|", "@q"),
                                                ("from", "@f"),
                                                ("to", "@t"),
                                                ("fuzziness", "@fuzz nt @f_dir"),
                                                ("id", "@i")],
                                      names=["jump_quads"],
                                      name="Hover jumps"))

        # the quads and x's for the calls:
        self.call_quad = ColumnDataSource({"x":[], "y":[], "w":[], "h":[]})
        self.plot.quad(left="x", bottom="y", right="w", top="h", line_color="magenta", line_width=3,
                       fill_alpha=0, source=self.call_quad, name="call_quad")

        self.ground_truth_quad = ColumnDataSource({"x":[], "y":[], "w":[], "h":[]})
        self.plot.quad(left="x", bottom="y", right="w", top="h", line_color="green",
                       line_width=3, fill_alpha=0, source=self.ground_truth_quad,
                       name="ground_truth_quad")

        self.call_x = ColumnDataSource({"x":[], "y":[]})
        self.plot.x(x="x", y="y", size=20, line_width=3, line_alpha=0.5, color="green",
                    source=self.call_x, name="call_x")

        self.ground_truth_x = ColumnDataSource({"x":[], "y":[]})
        self.plot.x(x="x", y="y", size=20, line_width=3, line_alpha=0.5, color="magenta",
                    source=self.ground_truth_x, name="ground_truth_x")

        self.plot.add_tools(HoverTool(tooltips=[("supp. nt", "@n"),
                                                ("coverage", "@c"),
                                                ("#reads", "@r"),
                                                ("score", "@s"),
                                                ("id", "@idx")],
                                      names=["call_quad", "ground_truth_quad", "call_x", "ground_truth_x"],
                                      name="Hover calls"))
        # the overview plot
        self.overview_quad = ColumnDataSource({"x":[], "y":[], "w":[], "h":[], "c":[]})
        self.plot.quad(left="x", bottom="y", right="w", top="h", color="c", line_width=0,
                       source=self.overview_quad, name="overview_quad")

        self.plot.add_tools(HoverTool(tooltips=[("from", "@f"), ("to", "@t"), ("#calls", "@i")],
                                      names=["overview_quad"],
                                      name="Hover heatmap"))

        # range change callback
        self.plot.y_range.on_change("start", lambda x,y,z: self.range_change_callback(renderer))

        # make jumps and calls clickable
        self.plot.on_event(Tap, lambda tap: self.jump_or_call_tap(renderer, tap.x, tap.y))

    def range_change_callback(self, renderer):
        x_r = self.plot.x_range
        y_r = self.plot.y_range
        # plot range still uninitialized
        if None in [x_r.end, x_r.start, y_r.end, y_r.start]:
            return
        w = renderer.xe - renderer.xs
        h = renderer.ye - renderer.ys
        if x_r.start < renderer.xs - w * renderer.render_area_factor or \
                y_r.start < renderer.ys - h * renderer.render_area_factor or \
                x_r.end > renderer.xe + w * renderer.render_area_factor or \
                y_r.end > renderer.ye + h * renderer.render_area_factor or \
                ( not renderer.redered_everything and ( w/4 > x_r.end - x_r.start or h/4 > y_r.end - y_r.start ) ):
            # there is actually a change in the range
            renderer.xs = x_r.start
            renderer.xe = x_r.end
            renderer.ys = y_r.start
            renderer.ye = y_r.end

            renderer.render()

    def update_selection(self, renderer):
        def highlight_jump(condition):
            for quad_idx, _ in enumerate(self.jump_quads):
                if len(self.jump_quads[quad_idx].data["c"]) > 0:
                    repl_dict = copy.copy(self.jump_quads[quad_idx].data)
                    for idx, _ in enumerate(repl_dict["c"]):
                        if condition(quad_idx, idx):
                            repl_dict["c"][idx] = ["orange", "blue", "lightgreen", "yellow"][quad_idx]
                        else:
                            repl_dict["c"][idx] = "lightgrey"
                    # this copy is inefficient
                    self.jump_quads[quad_idx].data = repl_dict

        if not renderer.selected_seed_id is None:
            # search tabbed seed
            seed_r = None
            seed_size = None
            for idx, _ in enumerate(renderer.seed_plot.seeds.data["center"]):
                if renderer.selected_seed_id == (renderer.seed_plot.seeds.data["idx"][idx],
                                                 renderer.seed_plot.seeds.data["r_id"][idx]):
                    seed_r = renderer.seed_plot.seeds.data["r"][idx]
                    seed_size = renderer.seed_plot.seeds.data["size"][idx]
                    break
            highlight_jump(lambda quad_idx, idx: \
                              ( self.jump_quads[quad_idx].data["f"][idx] == seed_r or \
                                self.jump_quads[quad_idx].data["f"][idx] == seed_r + seed_size - 1 or \
                                self.jump_quads[quad_idx].data["t"][idx] == seed_r or \
                                self.jump_quads[quad_idx].data["t"][idx] == seed_r + seed_size - 1 ) and \
                                self.jump_quads[quad_idx].data["r"][idx] == renderer.selected_read_id
                            )
        elif len(renderer.selected_jump_id) != 0:
            highlight_jump(lambda quad_idx, idx: self.jump_quads[quad_idx].data["i"][idx] in renderer.selected_jump_id)
        elif not renderer.selected_read_id is None:
            highlight_jump(lambda quad_idx, idx: self.jump_quads[quad_idx].data["r"][idx] == renderer.selected_read_id)
        else:
            highlight_jump(lambda quad_idx, idx: True)

    def jump_or_call_tap(self, renderer, x, y):
        renderer.selected_seed_id = None
        renderer.selected_call_id = set()
        renderer.selected_jump_id = set()

        def is_hit(dict, idx):
            return dict["x"][idx] <= x and x <= dict["w"][idx] and dict["y"][idx] <= y and y <= dict["h"][idx]

        # check calls first
        for idx, _ in enumerate(self.call_quad.data["x"]):
            if is_hit(self.call_quad.data, idx):
                renderer.selected_call_id.add(self.call_quad.data["idx"][idx])
                for idx in self.call_quad.data["supporing_jump_ids"][idx]:
                    renderer.selected_jump_id.add(idx)

        # check jumps if no call was hit
        sel_read_id = True
        if len(renderer.selected_call_id) == 0:
            for quad_idx, _ in enumerate(self.jump_quads):
                for idx, _ in enumerate(self.jump_quads[quad_idx].data["c"]):
                    if is_hit(self.jump_quads[quad_idx].data, idx):
                        renderer.selected_jump_id.add(self.jump_quads[quad_idx].data["i"][idx])
                        if sel_read_id == True:
                            sel_read_id = self.jump_quads[quad_idx].data["r"][idx]
                        elif sel_read_id != self.jump_quads[quad_idx].data["r"][idx]:
                            sel_read_id = False

        # if nothing at all was hit do total reset
        if len(renderer.selected_jump_id) == 0:
            renderer.selected_read_id = None

        if sel_read_id != True and sel_read_id != False:
            renderer.selected_read_id = sel_read_id
            renderer.read_plot.nuc_plot.copy_nts(renderer)

        self.update_selection(renderer)
        renderer.seed_plot.update_selection(renderer)

        if len(renderer.selected_jump_id) == 0:
            renderer.read_plot.reset_seeds(renderer)
            renderer.read_plot.nuc_plot.reset_nts(renderer)
        else:
            renderer.read_plot.copy_seeds(renderer, lambda idx: renderer.seed_plot.seeds.data["r_id"][idx] == \
                                                                renderer.selected_read_id)

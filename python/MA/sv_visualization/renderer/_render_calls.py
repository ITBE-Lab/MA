from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import FuncTickFormatter, TapTool, OpenURL, LabelSet, FixedTicker
from bokeh.models.callbacks import CustomJS
from MA import *
import math
from .util import *

def render_calls(self):
    rendered_everything = False
    self.params = ParameterSetManager()
    accepted_boxes_data = {
        "x": [],
        "w": [],
        "y": [],
        "h": [],
        "n": [],
        "c": [],
        "r": [],
        "s": []
    }
    accepted_plus_data = {
        "x": [],
        "y": [],
        "n": [],
        "c": [],
        "r": [],
        "s": []
    }
    ground_plus_data = {
        "x": [],
        "y": [],
        "n": [],
        "c": [],
        "r": [],
        "s": []
    }
    calls_from_db = SvCallsFromDb(self.params, self.sv_db, self.run_id, int(self.xs - self.w), int(self.ys - self.h),
                                  self.w*3, self.h*3, self.min_score)
    while calls_from_db.hasNext():
        def score(jump):
            if jump.coverage == 0:
                return None
            return str(jump.num_supp_nt / jump.coverage)
        jump = calls_from_db.next()
        if jump.num_supp_nt > self.min_score * jump.coverage:
            if jump.from_size == 1 and jump.to_size == 1:
                accepted_plus_data["x"].append(jump.from_start)
                accepted_plus_data["y"].append(jump.to_start)
            else:
                accepted_boxes_data["x"].append(jump.from_start - 0.5)
                accepted_boxes_data["y"].append(jump.to_start - 0.5)
                accepted_boxes_data["w"].append(
                    jump.from_start + jump.from_size + 1)
                accepted_boxes_data["h"].append(
                    jump.to_start + jump.to_size + 1)
                accepted_boxes_data["n"].append(jump.num_supp_nt)
                accepted_boxes_data["c"].append(jump.coverage)
                accepted_boxes_data["r"].append(len(jump.supporing_jump_ids))
                accepted_boxes_data["s"].append(score(jump))
                accepted_plus_data["x"].append(
                    jump.from_start + jump.from_size/2)
                accepted_plus_data["y"].append(jump.to_start + jump.to_size/2)
            accepted_plus_data["n"].append(jump.num_supp_nt)
            accepted_plus_data["c"].append(jump.coverage)
            accepted_plus_data["r"].append(len(jump.supporing_jump_ids))
            accepted_plus_data["s"].append(score(jump))
    calls_from_db = SvCallsFromDb(self.params, self.sv_db, self.ground_truth_id,
                                  int(self.xs - self.w), int(self.ys - self.h), self.w*3, self.h*3, self.min_score)
    while calls_from_db.hasNext():
        def score(jump):
            if jump.coverage == 0:
                return ""
            return " score: " + str(jump.num_supp_nt / jump.coverage)
        jump = calls_from_db.next()
        if jump.num_supp_nt > self.min_score * jump.coverage:
            if jump.from_size == 1 and jump.to_size == 1:
                ground_plus_data["x"].append(
                    jump.from_start + jump.from_size/2)
                ground_plus_data["y"].append(jump.to_start + jump.to_size/2)
                ground_plus_data["n"].append(jump.num_supp_nt)
                ground_plus_data["c"].append(jump.coverage)
                ground_plus_data["r"].append(len(jump.supporing_jump_ids))
                ground_plus_data["s"].append(score(jump))
            else:
                print("ground truth with fuzziness?!?!")

    num_jumps = libMA.get_num_jumps_in_area(self.sv_db, self.pack, self.sv_db.get_run_jump_id(self.run_id),
                                            int(self.xs - self.w), int(self.ys - self.h), self.w*3, self.h*3)
    if num_jumps < self.max_num_ele:
        rendered_everything = self.render_jumps()
    # the sv - boxes
    self.plot.quad(left="x", bottom="y", right="w", top="h", line_color="magenta", line_width=3, fill_alpha=0,
                   source=ColumnDataSource(accepted_boxes_data), name="hover2")
    self.plot.x(x="x", y="y", size=20, line_width=3, line_alpha=0.5, color="green",
                source=ColumnDataSource(ground_plus_data), name="hover2")
    self.plot.x(x="x", y="y", size=20, line_width=3, line_alpha=0.5, color="magenta",
                source=ColumnDataSource(accepted_plus_data), name="hover2")

    return rendered_everything

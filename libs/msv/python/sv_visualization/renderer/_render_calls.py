from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import FuncTickFormatter, TapTool, OpenURL, LabelSet, FixedTicker
from bokeh.models.callbacks import CustomJS
from MA import *
from MSV import *
import math
from .util import *

def render_calls(self):
    accepted_boxes_data = {
        "x": [],
        "w": [],
        "y": [],
        "h": [],
        "n": [],
        "c": [],
        "r": [],
        "idx": [],
        "supporing_jump_ids": [],
        "s": []
    }
    accepted_plus_data = {
        "x": [],
        "y": [],
        "n": [],
        "c": [],
        "r": [],
        "idx": [],
        "s": []
    }
    ground_boxes_data = {
        "x": [],
        "w": [],
        "y": [],
        "h": [],
        "n": [],
        "c": [],
        "r": [],
        "idx": [],
        "supporing_jump_ids": [],
        "s": []
    }
    ground_plus_data = {
        "x": [],
        "y": [],
        "n": [],
        "c": [],
        "r": [],
        "idx": [],
        "s": []
    }
    num_call_jumps = 0
    jump_list = []
    with self.measure("SvCallFromDb(run_id)"):
        calls_from_db = SvCallsFromDb(self.params, self.db_conn, self.get_run_id(), int(self.xs - self.w),
                                      int(self.ys - self.h), self.w*3, self.h*3, self.get_min_score())
    while calls_from_db.hasNext():
        jump = calls_from_db.next()
        if self.do_render_call_jumps_only:
            num_call_jumps += len(jump.supporing_jump_ids)
            for idx in range(len(jump.supporing_jump_ids)):
                jump_list.append(jump.get_jump(idx))
        if jump.x.size == 0 and jump.y.size == 0:
            accepted_plus_data["x"].append(jump.x.start + 0.5)
            accepted_plus_data["y"].append(jump.y.start + 0.5)
        else:
            accepted_boxes_data["x"].append(jump.x.start - 1)
            accepted_boxes_data["y"].append(jump.y.start - 1)
            accepted_boxes_data["w"].append(
                jump.x.start + jump.x.size + 1)
            accepted_boxes_data["h"].append(
                jump.y.start + jump.y.size + 1)
            accepted_boxes_data["n"].append(jump.num_supp_reads)
            accepted_boxes_data["c"].append(jump.reference_ambiguity)
            accepted_boxes_data["r"].append(len(jump.supporing_jump_ids))
            accepted_boxes_data["s"].append(str(jump.get_score()))
            accepted_boxes_data["idx"].append(jump.id)
            accepted_boxes_data["supporing_jump_ids"].append(list(jump.supporing_jump_ids))
            accepted_plus_data["x"].append(
                jump.x.start + jump.x.size/2)
            accepted_plus_data["y"].append(jump.y.start + jump.y.size/2)
        accepted_plus_data["idx"].append(jump.id)
        accepted_plus_data["n"].append(jump.num_supp_reads)
        accepted_plus_data["c"].append(jump.reference_ambiguity)
        accepted_plus_data["r"].append(len(jump.supporing_jump_ids))
        accepted_plus_data["s"].append(str(jump.get_score()))
    with self.measure("SvCallFromDb(run_id)"):
        calls_from_db = SvCallsFromDb(self.params, self.db_conn, self.get_gt_id(),
                                      int(self.xs - self.w), int(self.ys - self.h), self.w*3, self.h*3,
                                      self.get_min_score())
    while calls_from_db.hasNext():
        jump = calls_from_db.next()
        if jump.x.size == 0 and jump.y.size == 0:
            ground_plus_data["x"].append(jump.x.start + 0.5)
            ground_plus_data["y"].append(jump.y.start + 0.5)
        else:
            ground_boxes_data["x"].append(jump.x.start - 1)
            ground_boxes_data["y"].append(jump.y.start - 1)
            ground_boxes_data["w"].append(jump.x.start + jump.x.size + 1)
            ground_boxes_data["h"].append(jump.y.start + jump.y.size + 1)
            ground_boxes_data["s"].append(str(jump.get_score()))
            ground_boxes_data["idx"].append(jump.id)
            ground_boxes_data["supporing_jump_ids"].append(list(jump.supporing_jump_ids))
            ground_boxes_data["n"].append(jump.num_supp_reads)
            ground_boxes_data["c"].append(jump.reference_ambiguity)
            ground_boxes_data["r"].append(len(jump.supporing_jump_ids))
            ground_plus_data["x"].append(jump.x.start + jump.x.size/2)
            ground_plus_data["y"].append(jump.y.start + jump.y.size/2)
        ground_plus_data["idx"].append(jump.id)
        ground_plus_data["n"].append(jump.num_supp_reads)
        ground_plus_data["c"].append(jump.reference_ambiguity)
        ground_plus_data["r"].append(len(jump.supporing_jump_ids))
        ground_plus_data["s"].append(str(jump.get_score()))

    # the sv - boxes
    self.main_plot.call_quad.data = accepted_boxes_data
    self.main_plot.ground_truth_quad.data = ground_boxes_data
    self.main_plot.call_x.data = ground_plus_data
    self.main_plot.ground_truth_x.data = accepted_plus_data

    with self.measure("get_num_jumps_in_area"):
        if self.do_render_call_jumps_only:
            num_jumps = num_call_jumps
        else:
            num_jumps = get_num_jumps_in_area(self.db_conn, self.pack,
                                                SvCallerRunTable(self.db_conn).jump_run_id(self.get_run_id()),
                                                int(self.xs - self.w), int(self.ys - self.h), self.w*3, self.h*3,
                                                self.get_max_num_ele())
    if num_jumps < self.get_max_num_ele():
        with self.measure("render_jumps"):
            self.render_jumps(jump_list)
    else:
        self.analyze.analyze()

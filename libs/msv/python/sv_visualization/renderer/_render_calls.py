from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import FuncTickFormatter, TapTool, OpenURL, LabelSet, FixedTicker
from bokeh.models.callbacks import CustomJS
from MA import *
from MSV import *
import math
from .util import *

def render_calls(self, render_all=False):
    accepted_boxes_data = {
        "x": [],
        "w": [],
        "y": [],
        "h": [],
        "n": [],
        "c": [],
        "col": [],
        "idx": [],
        "desc": [],
        "supporing_jump_ids": [],
        "s": []
    }
    accepted_plus_data = {
        "x": [],
        "y": [],
        "n": [],
        "c": [],
        "col": [],
        "idx": [],
        "desc": [],
        "s": []
    }
    ground_boxes_data = {
        "x": [],
        "w": [],
        "y": [],
        "h": [],
        "n": [],
        "c": [],
        "col": [],
        "idx": [],
        "desc": [],
        "supporing_jump_ids": [],
        "s": []
    }
    ground_plus_data = {
        "x": [],
        "y": [],
        "n": [],
        "c": [],
        "col": [],
        "idx": [],
        "desc": [],
        "s": []
    }
    call_colors = {
        True: {
            True: "darkblue",
            False: "darkturquoise",
        },
        False: {
            True: "darkviolet",
            False: "darkorange",
        }
    }
    desc_table = CallDescTable(self.db_conn)
    num_call_jumps = 0
    jump_list = []
    with self.measure("SvCallFromDb(run_id)"):
        calls_from_db = SvCallsFromDb(self.params, self.db_conn, self.get_run_id(), int(self.xs - self.w),
                                      int(self.ys - self.h), self.w*3, self.h*3, self.get_min_score())
    with self.measure("SvCallFromDb(run_id) extract"):
        while calls_from_db.hasNext():
            call = calls_from_db.next(self.do_render_call_jumps_only)
            if self.do_render_call_jumps_only:
                num_call_jumps += len(call.supporing_jump_ids)
                for idx in range(len(call.supporing_jump_ids)):
                    jump_list.append(call.get_jump(idx))
            if call.x.size == 0 and call.y.size == 0:
                accepted_plus_data["x"].append(call.x.start + 0.5)
                accepted_plus_data["y"].append(call.y.start + 0.5)
            else:
                accepted_boxes_data["x"].append(call.x.start - 1)
                accepted_boxes_data["y"].append(call.y.start - 1)
                accepted_boxes_data["w"].append(
                    call.x.start + call.x.size + 1)
                accepted_boxes_data["h"].append(
                    call.y.start + call.y.size + 1)
                accepted_boxes_data["n"].append(call.num_supp_reads)
                accepted_boxes_data["c"].append(call.reference_ambiguity)
                accepted_boxes_data["col"].append(call_colors[call.from_forward][call.to_forward])
                accepted_boxes_data["s"].append(str(call.get_score()))
                accepted_boxes_data["idx"].append(call.id)
                accepted_boxes_data["desc"].append(desc_table.get_desc(call.id))
                accepted_boxes_data["supporing_jump_ids"].append(list(call.supporing_jump_ids))
                accepted_plus_data["x"].append(call.x.start + call.x.size/2)
                accepted_plus_data["y"].append(call.y.start + call.y.size/2)
            accepted_plus_data["idx"].append(call.id)
            accepted_plus_data["n"].append(call.num_supp_reads)
            accepted_plus_data["c"].append(call.reference_ambiguity)
            accepted_plus_data["col"].append(call_colors[call.from_forward][call.to_forward])
            accepted_plus_data["s"].append(str(call.get_score()))
            accepted_plus_data["desc"].append(desc_table.get_desc(call.id))
    with self.measure("SvCallFromDb(gt_id)"):
        calls_from_db_gt = SvCallsFromDb(self.params, self.db_conn, self.get_gt_id(),
                                    int(self.xs - self.w), int(self.ys - self.h), self.w*3, self.h*3,
                                    self.get_min_score())
    with self.measure("SvCallFromDb(gt_id) extract"):
        while calls_from_db_gt.hasNext():
            call = calls_from_db_gt.next()
            if call.x.size == 0 and call.y.size == 0:
                ground_plus_data["x"].append(call.x.start + 0.5)
                ground_plus_data["y"].append(call.y.start + 0.5)
            else:
                ground_boxes_data["x"].append(call.x.start - 1)
                ground_boxes_data["y"].append(call.y.start - 1)
                ground_boxes_data["w"].append(call.x.start + call.x.size + 1)
                ground_boxes_data["h"].append(call.y.start + call.y.size + 1)
                ground_boxes_data["s"].append(str(call.get_score()))
                ground_boxes_data["idx"].append(call.id)
                ground_boxes_data["supporing_jump_ids"].append(list(call.supporing_jump_ids))
                ground_boxes_data["n"].append(call.num_supp_reads)
                ground_boxes_data["c"].append(call.reference_ambiguity)
                ground_boxes_data["col"].append(call_colors[call.from_forward][call.to_forward])
                ground_plus_data["x"].append(call.x.start + call.x.size/2)
                ground_plus_data["y"].append(call.y.start + call.y.size/2)
                ground_boxes_data["desc"].append(desc_table.get_desc(call.id))
            ground_plus_data["idx"].append(call.id)
            ground_plus_data["n"].append(call.num_supp_reads)
            ground_plus_data["c"].append(call.reference_ambiguity)
            ground_plus_data["col"].append(call_colors[call.from_forward][call.to_forward])
            ground_plus_data["s"].append(str(call.get_score()))
            ground_plus_data["desc"].append(desc_table.get_desc(call.id))

    # the sv - boxes
    def callback():
        self.main_plot.call_quad.data = accepted_boxes_data
        self.main_plot.call_x.data = accepted_plus_data
        self.main_plot.ground_truth_quad.data = ground_boxes_data
        self.main_plot.ground_truth_x.data = ground_plus_data
    self.curdoc.add_next_tick_callback(callback)

    with self.measure("get_num_jumps_in_area"):
        if self.do_render_call_jumps_only:
            num_jumps = num_call_jumps
        else:
            num_jumps = get_num_jumps_in_area(self.db_conn, self.pack,
                                                SvCallerRunTable(self.db_conn).jump_run_id(self.get_run_id()),
                                                int(self.xs - self.w), int(self.ys - self.h), self.w*3, self.h*3,
                                                self.get_max_num_ele())
    if num_jumps < self.get_max_num_ele() or render_all:
        with self.measure("render_jumps"):
            self.render_jumps(jump_list, render_all)
    else:
        self.analyze.analyze()

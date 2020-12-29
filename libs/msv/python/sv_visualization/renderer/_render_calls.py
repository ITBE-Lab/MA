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
    desc_table = CallDescTable(self.db_conn_2)
    jump_list = []
    stats_data_1 = None
    if self.widgets.compute_stats() and self.read_plot.recalc_stat and self.get_run_id() != -1 and \
                    self.get_gt_id() != -1:
        self.read_plot.recalc_stat = False
        with self.measure("count - SvCallFromDb(run_id)"):
            min_blur = 0
            max_blur = 500
            blur_step = 10
            stats_data_1 = {"l":[], "x":[], "y":[], "c":[]}
            stats_data_4 = {"l":[], "x":[], "y":[], "c":[]}
            stats_data_2 = {"l":[], "x":[], "y":[], "c":[]}
            stats_data_3 = {"x":[], "w":[], "t":[], "c":[], "l":[]}
            print("computing stats 1...")
            stats, gt_total = self.count_calls_from_db.count(self.get_run_id(), self.get_gt_id(),
                                                                         self.widgets.get_blur())
            print("done")
            stats_data_1["l"].append("Ground Truth")
            stats_data_1["x"].append([self.get_min_score(), self.get_max_score()])
            stats_data_1["y"].append([gt_total, gt_total])
            stats_data_1["c"].append("black")
            stats_data_1["l"].append("Calls")
            stats_data_1["x"].append([])
            stats_data_1["y"].append([])
            stats_data_1["c"].append("blue")
            stats_data_1["l"].append("true-positives")
            stats_data_1["x"].append([])
            stats_data_1["y"].append([])
            stats_data_1["c"].append("green")

            stats_data_4["x"].append([])
            stats_data_4["y"].append([])
            stats_data_4["c"].append("green")
            stats_data_4["l"].append("Accuracy")
            stats_data_4["x"].append([])
            stats_data_4["y"].append([])
            stats_data_4["c"].append("blue")
            stats_data_4["l"].append("Recall")
            for x, num_calls, num_tp in stats:
                for n in [-1, -2]:
                    stats_data_1["x"][n].append(x)
                    if num_calls > 0 and gt_total > 0:
                        stats_data_4["x"][n].append(x)
                stats_data_1["y"][-1].append(num_tp)
                stats_data_1["y"][-2].append(num_calls)
                if num_calls > 0 and gt_total > 0:
                    stats_data_4["y"][-1].append(num_tp/gt_total) # recall
                    stats_data_4["y"][-2].append(num_tp/num_calls) # accuracy
            if self.print_to_tsv:
                with open(self.widgets.file_input.value + "-" + self.widgets.run_id_dropdown.value + ".tsv",
                          "w") as out_file:
                    out_file.write("//|ground truth| = " + str(gt_total) + "\n")
                    out_file.write("//#supporting reads\t#true positives\t#num entries\trecall\taccuracy\n")
                    for x, num_calls, num_tp in stats:
                        out_file.write(str(x) + "\t" + str(num_tp) + "\t" + str(num_calls) +
                                       "\t" + str(num_tp/gt_total) + "\t" + str(num_tp/num_calls) + "\n")
            stats_data_2["l"].append("Ground Truth")
            stats_data_2["x"].append([min_blur, max_blur])
            stats_data_2["y"].append([gt_total, gt_total])
            stats_data_2["c"].append("black")
            stats_data_2["l"].append("true-positives")
            stats_data_2["x"].append([])
            stats_data_2["y"].append([])
            stats_data_2["c"].append("green")
            if False: # deprecated ....
                for x, y in blur_stats:
                    stats_data_2["x"][-1].append(x)
                    stats_data_2["y"][-1].append(y)
            if False: # deprecated ....
                print("computing stats 2...")
                tp_bars, fp_bars, fn_bars, bar_width = self.count_calls_from_db.count_by_supp_nt(self.get_run_id(),
                                                                self.get_gt_id(),
                                                                self.widgets.get_blur(), 50,
                                                                self.get_min_score(), self.get_max_score(),
                                                                10000)
                print("done")
                for legend, color, bar in [("false-positive", "red", fp_bars),
                                        ("true-positive", "green", tp_bars),
                                        ("false-negative", "orange", fn_bars)]:
                    for c, t in bar:
                        stats_data_3["x"].append(c)
                        stats_data_3["w"].append(bar_width)
                        stats_data_3["t"].append(t)
                        stats_data_3["c"].append(color)
                        stats_data_3["l"].append(legend)

    with self.measure("SvCallFromDb(run_id)"):
        default_args = [self.get_run_id(),
                        int(self.xs - self.w),
                        int(self.ys - self.h),
                        self.w*3,
                        self.h*3]
        if self.widgets.get_render_t_p() and self.widgets.get_render_f_p():
            self.calls_from_db.init(*default_args,
                                        self.get_min_score(),
                                        self.get_max_score())
        elif self.widgets.get_render_t_p() or self.widgets.get_render_f_p():
            self.calls_from_db.init(*default_args,
                                        self.get_gt_id(),
                                        self.widgets.get_render_t_p(),
                                        self.widgets.get_blur(),
                                        self.get_min_score(),
                                        self.get_max_score())

    with self.measure("SvCallFromDb(run_id) extract"):
        while (self.widgets.get_render_t_p() or self.widgets.get_render_f_p()) and self.calls_from_db.hasNext():
            call = self.calls_from_db.next(self.do_render_call_jumps_only)
            if self.do_render_call_jumps_only:
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
        default_args = [self.get_gt_id(),
                        int(self.xs - self.w),
                        int(self.ys - self.h),
                        self.w*3,
                        self.h*3]
        if self.widgets.get_render_t_p() and self.widgets.get_render_f_n():
            self.calls_from_db_gt.init(*default_args)
        elif self.widgets.get_render_t_p() or self.widgets.get_render_f_n():
            self.calls_from_db_gt.init(self.get_min_score(),
                                        self.get_max_score(),
                                        *default_args,
                                        self.get_run_id(),
                                        self.widgets.get_render_t_p(),
                                        self.widgets.get_blur())
    with self.measure("SvCallFromDb(gt_id) extract"):
        while (self.widgets.get_render_t_p() or self.widgets.get_render_f_n()) and self.calls_from_db_gt.hasNext():
            call = self.calls_from_db_gt.next(False)
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
        if not stats_data_1 is None:
            self.read_plot.stat_lines_1.data = stats_data_1
            self.read_plot.stat_lines_2.data = stats_data_2
            self.read_plot.stat_lines_4.data = stats_data_4
            self.read_plot.stat_lines_3.data = stats_data_3
    self.do_callback(callback)

    with self.measure("render_jumps"):
        self.render_jumps(jump_list, render_all)

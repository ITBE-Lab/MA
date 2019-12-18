from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import FuncTickFormatter, TapTool, OpenURL, LabelSet, FixedTicker
from bokeh.models.callbacks import CustomJS
from MA import *
import math
from .util import *

def render_jumps(self):
    self.read_ids = set()
    out_dicts = []
    patch = {
        "x": [],
        "y": []
    }
    for _ in range(4):
        out_dicts.append({
            "x": [],
            "y": [],
            "w": [],
            "h": [],
            "a": [],
            "n": [],
            "r": [],
            "q": [],
            "f": [],
            "t": [],
            "c": [],
            "fuzz": [],
            "f_dir": [],
            "i": []
        })
    with self.measure("SortedSvJumpFromSql"):
        sweeper = SortedSvJumpFromSql(self.params, self.sv_db, self.sv_db.get_run_jump_id(self.get_run_id()),
                                        int(self.xs - self.w), int(self.ys - self.h), self.w*3, self.h*3)
    with self.measure("render jumps"):
        while sweeper.has_next_start():
            jump = sweeper.get_next_start()
            idx = None
            if jump.switch_strand_known():
                if jump.does_switch_strand():
                    idx = 0
                else:
                    idx = 1
            else:
                if jump.from_known():
                    idx = 2
                else:
                    idx = 3

            f_dir = ""
            if not jump.from_fuzziness_is_rightwards():
                f_dir = "left-"
            else:
                f_dir = "right-"
            if not jump.to_fuzziness_is_downwards():
                f_dir += "up"
            else:
                f_dir += "down"

            out_dicts[idx]["c"].append("lightgrey")
            out_dicts[idx]["f_dir"].append(f_dir)
            out_dicts[idx]["f"].append(jump.from_pos if jump.from_known() else "unknown")
            out_dicts[idx]["t"].append(jump.to_pos if jump.to_known() else "unknown")
            out_dicts[idx]["x"].append(jump.from_start_same_strand())
            out_dicts[idx]["y"].append(jump.to_start())
            out_dicts[idx]["w"].append(
                jump.from_start_same_strand() + jump.from_size() + 1.0)
            out_dicts[idx]["h"].append(jump.to_start() + jump.to_size() + 1.0)
            out_dicts[idx]["a"].append(jump.num_supp_nt() / 1000)
            out_dicts[idx]["n"].append(jump.num_supp_nt())
            out_dicts[idx]["r"].append(jump.read_id)
            out_dicts[idx]["q"].append(jump.query_distance())
            out_dicts[idx]["i"].append(jump.id)
            if jump.from_known() and jump.to_known():
                out_dicts[idx]["fuzz"].append(jump.fuzziness())
            else:
                out_dicts[idx]["fuzz"].append(jump.query_distance())
            self.read_ids.add(jump.read_id)

            f = jump.from_pos
            t = jump.to_pos
            if not jump.from_known():
                f = t
            if not jump.to_known():
                t = f
            patch["x"].append([f + 0.25, f + 0.75])
            patch["x"].append([f + 0.25, f + 0.75])
            patch["y"].append([t + 0.25, t + 0.75])
            patch["y"].append([t + 0.75, t + 0.25])
        for idx in range(4):
            self.main_plot.jump_quads[idx].data = out_dicts[idx]
        self.main_plot.jump_x.data = patch
        self.main_plot.update_selection(self)

    if self.w*3+self.h*3 < self.get_max_num_ele():
        # render nucs in read plot
        # render nucs in sv plot
        with self.measure("render_nucs"):
            self.render_nucs()

    # render the seeds in the main seed plot
    if len(self.read_ids)*self.read_penalty_factor < self.get_max_num_ele():
        with self.measure("render_reads"):
            self.render_reads()

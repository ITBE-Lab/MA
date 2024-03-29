from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import FuncTickFormatter, TapTool, OpenURL, LabelSet, FixedTicker
from bokeh.models.callbacks import CustomJS
from MA import *
from MSV import *
import math
from .util import *

def render_jumps(self, jump_list=[], render_all=False):
    with self.measure("get_num_jumps_in_area"):
        num_jumps = 0
        if True:
            if self.do_render_call_jumps_only:
                num_jumps = len(jump_list)
            else:
                num_jumps = get_num_jumps_in_area(self.db_conn, self.pack,
                                                    SvCallerRunTable(self.db_conn).jump_run_id(self.get_run_id()),
                                                    int(self.xs - self.w), int(self.ys - self.h), self.w*3, self.h*3,
                                                    self.get_max_num_ele())
    self.read_ids = set()
    if num_jumps < self.get_max_num_ele() or render_all:
        out_dicts = []
        patch = {
            "x": [],
            "y": []
        }
        for _ in range(6):
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
                "fs": [],
                "t": [],
                "ts": [],
                "m": [],
                "c": [],
                "fuzz": [],
                "f_dir": [],
                "ref_ambig": [],
                "i": []
            })
        if not self.do_render_call_jumps_only:
            with self.measure("SvJumpFromSql"):
                sweeper = SvJumpFromSql(self.db_conn, 
                                            SvCallerRunTable(self.db_conn).jump_run_id(self.get_run_id()),
                                            int(self.xs - self.w), int(self.ys - self.h), self.w*3, self.h*3,
                                            self.get_max_num_ele())
        with self.measure("render jumps"):
            if not self.do_render_call_jumps_only:
                while sweeper.has_next():
                    jump_list.append(sweeper.get_next())
            max_supp_nt = 1
            for jump in jump_list:
                max_supp_nt = max(jump.num_supp_nt(), max_supp_nt)
            for jump in jump_list:
                idx = None
                if not jump.to_known():
                    idx = 0
                elif not jump.from_known():
                    idx = 1
                else:
                    idx = 2
                    if not jump.from_forward:
                        idx += 2
                    if not jump.to_forward:
                        idx += 1

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
                out_dicts[idx]["ref_ambig"].append(jump.ref_ambiguity(5, 20, self.pack))
                out_dicts[idx]["f"].append(jump.from_pos if jump.from_known() else "unknown")
                out_dicts[idx]["fs"].append("forw" if jump.from_forward else "rev")
                out_dicts[idx]["t"].append(jump.to_pos if jump.to_known() else "unknown")
                out_dicts[idx]["ts"].append("forw" if jump.to_forward else "rev")
                out_dicts[idx]["m"].append("yes" if jump.was_mirrored else "no")
                out_dicts[idx]["x"].append(jump.from_start_same_strand())
                out_dicts[idx]["y"].append(jump.to_start())
                out_dicts[idx]["w"].append(
                    jump.from_start_same_strand() + jump.from_size() + 1.0)
                out_dicts[idx]["h"].append(jump.to_start() + jump.to_size() + 1.0)
                if not jump.from_known():
                    out_dicts[idx]["y"][-1] += jump.from_size() + 1
                    out_dicts[idx]["h"][-1] += jump.from_size() + 1
                out_dicts[idx]["a"].append(0.1)
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
                    t += 1
                if not jump.to_known():
                    t = f
                    f -= 1
                patch["x"].append([f + 0.25, f + 0.75])
                patch["x"].append([f + 0.25, f + 0.75])
                patch["y"].append([t + 0.25, t + 0.75])
                patch["y"].append([t + 0.75, t + 0.25])

            def callback():
                for idx in range(6):
                    self.main_plot.jump_quads[idx].data = out_dicts[idx]
                self.main_plot.jump_x.data = patch
                self.main_plot.update_selection(self)
            self.do_callback(callback)

    if self.w*3+self.h*3 < self.get_max_num_ele() or render_all:
        # render nucs in read plot
        # render nucs in sv plot
        with self.measure("render_nucs"):
            self.render_nucs()

    # render the seeds in the main seed plot
    if len(self.read_ids)*self.read_penalty_factor >= self.get_max_num_ele() and not render_all:
        print("not rendering", len(self.read_ids), "reads")
        self.read_ids = set()
    with self.measure("render_reads"):
        self.render_reads(render_all)

    self.analyze.analyze()

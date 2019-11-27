from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import FuncTickFormatter, TapTool, OpenURL, LabelSet, FixedTicker
from bokeh.models.callbacks import CustomJS
from MA import *
import math
from .util import *

PATCH = False

def render_jumps(self):
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
    sweeper = SortedSvJumpFromSql(self.params, self.sv_db, self.sv_db.get_run_jump_id(self.run_id),
                                    int(self.xs - self.w), int(self.ys - self.h), self.w*3, self.h*3)
    while sweeper.has_next_start():
        jump = sweeper.get_next_start()
        idx = None
        if jump.switch_strand_known():
            if jump.does_switch_strand():
                idx = 0
                out_dicts[idx]["c"].append("orange")
            else:
                idx = 1
                out_dicts[idx]["c"].append("blue")
        else:
            if jump.from_known():
                idx = 2
                out_dicts[idx]["c"].append("lightgreen")
            else:
                idx = 3
                out_dicts[idx]["c"].append("yellow")

        if self.selected_read_id != -1 and self.selected_read_id != jump.read_id:
            out_dicts[idx]["c"][-1] = "lightgrey"

        f_dir = ""
        if not jump.from_fuzziness_is_rightwards():
            f_dir = "left-"
        else:
            f_dir = "right-"
        if not jump.to_fuzziness_is_downwards():
            f_dir += "up"
        else:
            f_dir += "down"

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
        if PATCH:
            if not jump.from_fuzziness_is_rightwards():
                if not jump.to_fuzziness_is_downwards():
                    patch["x"].extend([f - 2.5, f + .5, f + .5, float("NaN")])
                    patch["y"].extend([t - .5, t + 2.5, t - .5, float("NaN")])
                else:
                    patch["x"].extend([f - 2.5, f + .5, f + .5, float("NaN")])
                    patch["y"].extend([t + .5, t - 2.5, t + .5, float("NaN")])
            else:
                if not jump.to_fuzziness_is_downwards():
                    patch["x"].extend([f + 2.5, f - .5, f - .5, float("NaN")])
                    patch["y"].extend([t - .5, t + 2.5, t - .5, float("NaN")])
                else:
                    patch["x"].extend([f + 2.5, f - .5, f - .5, float("NaN")])
                    patch["y"].extend([t + .5, t - 2.5, t + .5, float("NaN")])
        else:
            patch["x"].append([f + 0.25, f + 0.75])
            patch["x"].append([f + 0.25, f + 0.75])
            patch["y"].append([t + 0.25, t + 0.75])
            patch["y"].append([t + 0.75, t + 0.25])
    self.quads.append(self.plot.quad(left="x", bottom="y", right="w", top="h", fill_color="c", line_color="c",
                                line_width=3, fill_alpha="a", source=ColumnDataSource(out_dicts[0]), name="hover3"))
    self.quads.append(self.plot.quad(left="x", bottom="y", right="w", top="h", fill_color="c", line_color="c",
                                line_width=3, fill_alpha="a", source=ColumnDataSource(out_dicts[1]), name="hover3"))
    self.quads.append(self.plot.quad(left="x", bottom="y", right="w", top="h", fill_color="c", line_color="c",
                                line_width=3, fill_alpha="a", source=ColumnDataSource(out_dicts[2]), name="hover3"))
    self.quads.append(self.plot.quad(left="x", bottom="y", right="w", top="h", fill_color="c", line_color="c",
                                line_width=3, fill_alpha="a", source=ColumnDataSource(out_dicts[3]), name="hover3"))
    if PATCH:
        self.plot.patch(x="x", y="y", line_width=1, color="black",
                        source=ColumnDataSource(patch))
    else:
        self.plot.multi_line(xs="x", ys="y", line_width=1.5, line_alpha=0.5, color="black",
                        source=ColumnDataSource(patch))
    if len(self.read_ids) < self.max_num_ele:
        return self.render_reads()
    return False
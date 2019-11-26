from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import FuncTickFormatter, TapTool, OpenURL, LabelSet, FixedTicker
from bokeh.models.callbacks import CustomJS
from MA import *
import math
from .util import *

def render_jumps(self):
    read_ids = set()
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

        out_dicts[idx]["f"].append(jump.from_pos)
        out_dicts[idx]["t"].append(jump.to_pos)
        out_dicts[idx]["x"].append(jump.from_start_same_strand() - 0.5)
        out_dicts[idx]["y"].append(jump.to_start() - 0.5)
        out_dicts[idx]["w"].append(
            jump.from_start_same_strand() + jump.from_size() + 1)
        out_dicts[idx]["h"].append(jump.to_start() + jump.to_size() + 1)
        out_dicts[idx]["a"].append(jump.num_supp_nt() / 1000)
        out_dicts[idx]["n"].append(jump.num_supp_nt())
        out_dicts[idx]["r"].append(jump.read_id)
        out_dicts[idx]["q"].append(jump.query_distance())
        out_dicts[idx]["i"].append(jump.id)
        read_ids.add(jump.read_id)

        f = jump.from_pos
        t = jump.to_pos
        if not jump.from_known():
            f = t
        if not jump.to_known():
            t = f
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
    quads = []
    quads.append(self.plot.quad(left="x", bottom="y", right="w", top="h", fill_color="c", line_color="c",
                                line_width=3, fill_alpha="a", source=ColumnDataSource(out_dicts[0]), name="hover3"))
    quads.append(self.plot.quad(left="x", bottom="y", right="w", top="h", fill_color="c", line_color="c",
                                line_width=3, fill_alpha="a", source=ColumnDataSource(out_dicts[1]), name="hover3"))
    quads.append(self.plot.quad(left="x", bottom="y", right="w", top="h", fill_color="c", line_color="c",
                                line_width=3, fill_alpha="a", source=ColumnDataSource(out_dicts[2]), name="hover3"))
    quads.append(self.plot.quad(left="x", bottom="y", right="w", top="h", fill_color="c", line_color="c",
                                line_width=3, fill_alpha="a", source=ColumnDataSource(out_dicts[3]), name="hover3"))
    self.plot.patch(x="x", y="y", line_width=1, color="black",
                    source=ColumnDataSource(patch))
    if len(read_ids) < self.max_num_ele:
        seeder = BinarySeeding(self.params)
        jumps_from_seeds = libMA.SvJumpsFromSeeds(
            self.params, -1, self.sv_db, self.pack)
        read_dict = {
            "center": [],
            "r_id": [],
            "size": [],
            "q": [],
            "r": [],
            "l": [],
            "idx": [],
            "c": [],
            "f": [],
            "layer": [],
            "x": [],
            "y": [],
            "category": []
        }
        # create a column data source for the read plot...
        read_plot_dict = {
            "center": [],
            "r_id": [],
            "size": [],
            "q": [],
            "r": [],
            "l": [],
            "idx": [],
            "c": [],
            "f": [],
            "layer": [],
            "x": [],
            "y": [],
            "category": []
        }
        read_id_n_cols = {}
        col_ids = []
        all_col_ids = []
        category_counter = 0
        for read_id in sorted(read_ids, reverse=True):
            read = self.sv_db.get_read(read_id)
            segments = seeder.execute(self.fm_index, read)
            seeds = libMA.Seeds()
            layer_of_seeds = jumps_from_seeds.cpp_module.execute_helper(
                segments, self.pack, self.fm_index, read, seeds)
            end_column = []
            seeds_n_idx = list(enumerate(sorted([(x, y) for x, y in zip(seeds, layer_of_seeds)],
                                                key=lambda x: x[0].start)))
            for idx, (seed, layer) in sorted(seeds_n_idx, key=lambda x: x[1][0].start_ref):
                seed_size = seed.size - 1
                if seed.on_forward_strand:
                    read_dict["center"].append(
                        seed.start_ref + seed_size/2)
                    if self.selected_read_id != -1 and self.selected_read_id != read_id:
                        read_dict["c"].append("lightgrey")
                    else:
                        read_dict["c"].append("green")
                    read_dict["r"].append(seed.start_ref)
                    read_dict["x"].append(
                        [seed.start_ref, seed.start_ref+seed_size])
                    curr_end = seed.start_ref + seed_size
                    curr_start = seed.start_ref
                else:
                    read_dict["center"].append(
                        seed.start_ref - seed_size/2 + 1)
                    if self.selected_read_id != -1 and self.selected_read_id != read_id:
                        read_dict["c"].append("lightgrey")
                    else:
                        read_dict["c"].append("purple")
                    read_dict["r"].append(seed.start_ref - seed.size + 1)
                    read_dict["x"].append(
                        [seed.start_ref + 1, seed.start_ref - seed_size + 1])
                    curr_end = seed.start_ref
                    curr_start = seed.start_ref - seed.size
                curr_column = 0
                while curr_column < len(end_column):
                    if curr_start >= end_column[curr_column]:
                        break
                    else:
                        curr_column += 1
                if curr_column >= len(end_column):
                    end_column.append(0)
                    all_col_ids.append(curr_column + category_counter)
                end_column[curr_column] = curr_end
                read_dict["r_id"].append(read_id)
                read_dict["size"].append(seed_size)
                read_dict["l"].append(seed.size)
                read_dict["q"].append(seed.start)
                read_dict["y"].append([seed.start, seed.start+seed.size])
                read_dict["idx"].append(idx)
                read_dict["layer"].append(layer)
                read_dict["f"].append(seed.on_forward_strand)
                read_dict["category"].append(
                    category_counter + curr_column)
                if self.selected_read_id == read_id:
                    for key in read_dict.keys():
                        read_plot_dict[key].append(read_dict[key][-1])
            if (len(end_column)-1) % 2 == 0:
                # prevent forming of float if possible (to stay javascript compatible)
                curr_col_id = category_counter + (len(end_column)-1)//2
            else:
                curr_col_id = category_counter + (len(end_column)-1)/2
            col_ids.append(curr_col_id)
            category_counter += len(end_column) + 2
            read_id_n_cols[curr_col_id] = read_id
        if len(read_dict["c"]) < self.max_num_ele:
            read_source = ColumnDataSource(read_dict)
            self.l_plot[1].rect(x="category", y="center", width=1, height="size",
                                fill_color="c", line_width=0, source=read_source, name="hover5")
            self.d_plot[1].rect(y="category", x="center", height=1, width="size",
                                fill_color="c", line_width=0, source=read_source, name="hover5")
            self.l_plot[1].xaxis.ticker = FixedTicker(ticks=col_ids)
            self.l_plot[1].xaxis.formatter = FuncTickFormatter(
                args={"read_id_n_cols": read_id_n_cols},
                code="""
                        if(!tick in read_id_n_cols)
                            return "";
                        return read_id_n_cols[tick];
                    """)
            self.l_plot[1].xgrid.ticker = FixedTicker(ticks=all_col_ids)
            self.d_plot[1].yaxis.ticker = FixedTicker(ticks=col_ids)
            self.d_plot[1].yaxis.formatter = FuncTickFormatter(
                args={"read_id_n_cols": read_id_n_cols},
                code="""
                        if(!tick in read_id_n_cols)
                            return "";
                        return read_id_n_cols[tick];
                    """)
            self.d_plot[1].ygrid.ticker = FixedTicker(ticks=all_col_ids)

            num_nt = self.w*3+self.h*3
            if num_nt < self.max_num_ele:
                l_plot_nucs = {"y": [], "c": [], "i": []}

                def append_nuc_type(dict_, nuc):
                    if nuc == "A" or nuc == "a":
                        dict_["c"].append("blue")
                        dict_["i"].append("A")
                    elif nuc == "C" or nuc == "c":
                        dict_["c"].append("red")
                        dict_["i"].append("C")
                    elif nuc == "G" or nuc == "g":
                        dict_["c"].append("green")
                        dict_["i"].append("G")
                    elif nuc == "T" or nuc == "t":
                        dict_["c"].append("yellow")
                        dict_["i"].append("T")
                    else:
                        dict_["c"].append("lightgreen")
                        dict_["i"].append(nuc)
                nuc_seq = self.pack.extract_from_to(
                    int(self.ys - self.h), int(self.ye + self.h + 1))
                for y_add, nuc in enumerate(str(nuc_seq)):
                    l_plot_nucs["y"].append(
                        int(self.ys - self.h) + y_add + 0.5)
                    append_nuc_type(l_plot_nucs, nuc)
                self.l_plot[0].rect(x=0.5, y="y", width=1, height=1, fill_color="c", line_width=0,
                                    source=ColumnDataSource(l_plot_nucs), name="hover4")
                d_plot_nucs = {"x": [], "c": [], "i": []}
                nuc_seq = self.pack.extract_from_to(
                    int(self.xs - self.w), int(self.xe + self.w + 1))
                for x_add, nuc in enumerate(str(nuc_seq)):
                    d_plot_nucs["x"].append(
                        int(self.xs - self.w) + x_add + 0.5)
                    append_nuc_type(d_plot_nucs, nuc)
                self.d_plot[0].rect(x="x", y=0.5, width=1, height=1, fill_color="c", line_width=0,
                                    source=ColumnDataSource(d_plot_nucs), name="hover4")

                read_plot_line = self.read_plot.multi_line(xs="x", ys="y", line_color="c", line_width=5,
                                                            source=ColumnDataSource(read_plot_dict), name="hover5")
                # auto adjust y-range of read plot
                js_auto_adjust_y_range = js_file("auto_adjust")
                self.plot.x_range.js_on_change('start', CustomJS(args=dict(checkbox_group=self.checkbox_group,
                                                                            read_plot=self.read_plot, plot=self.plot,
                                                                            read_plot_line=read_plot_line.data_source),
                                                                    code=js_auto_adjust_y_range+"auto_adjust();"))

                # the tapping callback on jumps
                self.plot.js_on_event("tap", CustomJS(args=dict(srcs=[x.data_source for x in quads],
                                                                read_source=read_source),
                                                        code=js_file("jump_tap")))
                # the tapping callback on seeds
                code = js_auto_adjust_y_range+js_file("seed_tap")
                self.l_plot[1].js_on_event("tap", CustomJS(args=dict(srcs=[x.data_source for x in quads],
                                                                        checkbox_group=self.checkbox_group,
                                                                        plot=self.plot,
                                                                        read_source=read_source,
                                                                        range=self.d_plot[1].y_range,
                                                                        read_plot_line=read_plot_line.data_source,
                                                                        read_plot=self.read_plot),
                                                            code="""
                                                        var curr_x = cb_obj.y;
                                                        var curr_y = cb_obj.x;
                                                        """ + code))
                self.d_plot[1].js_on_event("tap", CustomJS(args=dict(srcs=[x.data_source for x in quads],
                                                                        checkbox_group=self.checkbox_group,
                                                                        plot=self.plot,
                                                                        read_source=read_source,
                                                                        range=self.d_plot[1].y_range,
                                                                        read_plot_line=read_plot_line.data_source,
                                                                        read_plot=self.read_plot),
                                                            code="""
                                                        var curr_x = cb_obj.x;
                                                        var curr_y = cb_obj.y;
                                                        """ + code))

                rendered_everything = True
    return rendered_everything
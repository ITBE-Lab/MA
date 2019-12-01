from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import FuncTickFormatter, TapTool, OpenURL, LabelSet, FixedTicker
from bokeh.models.callbacks import CustomJS
from bokeh.palettes import Plasma, Plasma256
from MA import *
import math
from .util import *


def render_reads(self):
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
    l_plot_nucs = {}  # dict of {"y": [], "c": [], "i": []}
    read_plot_rects = {}  # dict of {"l": [], "b": [], "t": [], "r": []}
    initial_l_data = {"y": [0], "c": ["black"], "i": [""]}
    initial_rect_plot_data = {"l": [], "b": [], "t": [], "r": []}
    for read_id in sorted(self.read_ids, reverse=True):
        l_plot_nucs[read_id] = {"y": [], "c": [], "i": []}
        read_plot_rects[read_id] = {"l": [], "b": [], "t": [], "r": []}
        read = self.sv_db.get_read(read_id)
        if self.selected_read_id == read_id:
            for value in initial_l_data.values():
                value.clear()
            for value in initial_rect_plot_data.values():
                value.clear()
        for y, nuc in enumerate(str(read)):
            append_nuc_type(l_plot_nucs[read_id], nuc, y, "y")
            if self.selected_read_id == read_id:
                append_nuc_type(initial_l_data, nuc, y, "y")
        segments = seeder.execute(self.fm_index, read)
        seeds = libMA.Seeds()
        layer_of_seeds, rectangles = jumps_from_seeds.cpp_module.execute_helper(
            segments, self.pack, self.fm_index, read, seeds)
        if self.render_mems == 1 and self.selected_read_id == read_id:
            hash_map_seeder = HashMapSeeding(self.params)
            hash_map_seeder.cpp_module.seed_size = 9
            query_section = NucSeq(str(read)[int(self.seed_plot_y_s):int(self.seed_plot_y_e)])
            ref_section = self.pack.extract_from_to(int(self.xs), int(self.xe))
            all_k_mers = hash_map_seeder.execute(query_section, ref_section)
            all_mems = SeedLumping(self.params).execute(all_k_mers)
            filter_module = FilterToUnique(self.params)
            filter_module.cpp_module.num_mm = 0
            filtered_mems = filter_module.execute(all_mems, query_section, ref_section)
            #filtered_mems = all_mems
            for idx in range(len(filtered_mems)):
                filtered_mems[idx].start += int(self.seed_plot_y_s)
                filtered_mems[idx].start_ref += int(self.xs)
                layer_of_seeds.append(-1)
            seeds.extend(filtered_mems)
        for rectangle in rectangles:
            if self.selected_read_id == read_id:
                initial_rect_plot_data["l"].append(rectangle.x_axis.start)
                initial_rect_plot_data["b"].append(rectangle.y_axis.start)
                initial_rect_plot_data["r"].append(rectangle.x_axis.start + rectangle.x_axis.size)
                initial_rect_plot_data["t"].append(rectangle.y_axis.start + rectangle.y_axis.size)
            read_plot_rects[read_id]["l"].append(rectangle.x_axis.start)
            read_plot_rects[read_id]["b"].append(rectangle.y_axis.start)
            read_plot_rects[read_id]["r"].append(rectangle.x_axis.start + rectangle.x_axis.size)
            read_plot_rects[read_id]["t"].append(rectangle.y_axis.start + rectangle.y_axis.size)
        end_column = []
        seeds_n_idx = list(enumerate(sorted([(x, y) for x, y in zip(seeds, layer_of_seeds)],
                                            key=lambda x: x[0].start)))
        max_seed_size = max(seed.size for seed in seeds)
        for idx, (seed, layer) in sorted(seeds_n_idx, key=lambda x: x[1][0].start_ref):
            seed_size = seed.size - 1
            if seed.on_forward_strand:
                read_dict["center"].append(seed.start_ref + seed.size/2)
                if self.selected_read_id != -1 and self.selected_read_id != read_id:
                    read_dict["c"].append("lightgrey")
                else:
                    if layer == -1:
                        read_dict["c"].append(Plasma256[ (255 * seed_size) // max_seed_size])
                    else:
                        read_dict["c"].append("green")
                read_dict["r"].append(seed.start_ref)
                read_dict["x"].append(
                    [seed.start_ref, seed.start_ref+seed.size])
                curr_end = seed.start_ref + seed_size + 1
                curr_start = seed.start_ref
            else:
                read_dict["center"].append(seed.start_ref - seed.size/2 + 1)
                if self.selected_read_id != -1 and self.selected_read_id != read_id:
                    read_dict["c"].append("lightgrey")
                else:
                    if layer == -1:
                        read_dict["c"].append(Plasma256[ (255 * seed_size) // max_seed_size])
                    else:
                        read_dict["c"].append("purple")
                read_dict["r"].append(seed.start_ref - seed.size + 1)
                read_dict["x"].append(
                    [seed.start_ref + 1, seed.start_ref - seed.size + 1])
                curr_end = seed.start_ref + 1
                curr_start = seed.start_ref - seed.size
            curr_column = 0
            while curr_column < len(end_column):
                if curr_start > end_column[curr_column]:
                    break
                else:
                    curr_column += 1
            if curr_column >= len(end_column):
                end_column.append(0)
                all_col_ids.append(curr_column + category_counter)
            end_column[curr_column] = curr_end
            read_dict["y"].append([seed.start, seed.start+seed.size])
            read_dict["r_id"].append(read_id)
            read_dict["size"].append(seed.size)
            read_dict["l"].append(seed.size)
            read_dict["q"].append(seed.start)
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
    if len(read_dict["c"]) < self.max_num_ele or self.render_mems == 1:
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

        rect_read_plot_data = self.read_plot.quad(left="l", bottom="b", right="r", top="t", fill_color="grey",
                                                  fill_alpha=0.2, line_width=0, name="hover6",
                                                  source=ColumnDataSource(initial_rect_plot_data))
        read_plot_line = self.read_plot.multi_line(xs="x", ys="y", line_color="c", line_width=5,
                                                   line_alpha=0.5 if self.render_mems == 1 else 1,
                                                   source=ColumnDataSource(read_plot_dict), name="hover5")
        l_read_plot_data = self.l_read_plot.rect(x=0.5, y="y", width=1, height=1, fill_color="c", line_width=0,
                                                 source=ColumnDataSource(initial_l_data), name="hover4")
        # auto adjust y-range of read plot
        js_auto_adjust_y_range = js_file("auto_adjust")
        self.plot.x_range.js_on_change('start', CustomJS(args=dict(radio_group=self.radio_group,
                                                                   read_plot=self.read_plot, plot=self.plot,
                                                                   read_plot_line=read_plot_line.data_source),
                                                         code=js_auto_adjust_y_range+"auto_adjust();"))

        # the tapping callback on jumps
        self.plot.js_on_event("tap", CustomJS(args=dict(srcs=[x.data_source for x in self.quads],
                                                        read_source=read_source,
                                                        l_plot_nucs=l_plot_nucs,
                                                        read_plot_line=read_plot_line.data_source,
                                                        l_read_plot_data=l_read_plot_data.data_source,
                                                        rect_read_plot_data=rect_read_plot_data.data_source,
                                                        read_plot_rects=read_plot_rects),
                                              code=js_file("jump_tap")))
        # the tapping callback on seeds
        code = js_auto_adjust_y_range+js_file("seed_tap")
        self.l_plot[1].js_on_event("tap", CustomJS(args=dict(srcs=[x.data_source for x in self.quads],
                                                             radio_group=self.radio_group,
                                                             plot=self.plot,
                                                             read_source=read_source,
                                                             range=self.d_plot[1].y_range,
                                                             read_plot_line=read_plot_line.data_source,
                                                             read_plot=self.read_plot,
                                                             l_plot_nucs=l_plot_nucs,
                                                             l_read_plot_data=l_read_plot_data.data_source,
                                                             rect_read_plot_data=rect_read_plot_data.data_source,
                                                             read_plot_rects=read_plot_rects),
                                                   code="""
                                                var curr_x = cb_obj.y;
                                                var curr_y = cb_obj.x;
                                                """ + code))
        self.d_plot[1].js_on_event("tap", CustomJS(args=dict(srcs=[x.data_source for x in self.quads],
                                                             radio_group=self.radio_group,
                                                             plot=self.plot,
                                                             read_source=read_source,
                                                             range=self.d_plot[1].y_range,
                                                             read_plot_line=read_plot_line.data_source,
                                                             read_plot=self.read_plot,
                                                             l_plot_nucs=l_plot_nucs,
                                                             l_read_plot_data=l_read_plot_data.data_source,
                                                             rect_read_plot_data=rect_read_plot_data.data_source,
                                                             read_plot_rects=read_plot_rects),
                                                   code="""
                                                var curr_x = cb_obj.x;
                                                var curr_y = cb_obj.y;
                                                    """ + code))

        num_nt = self.w*3+self.h*3
        if num_nt < self.max_num_ele:
            # render nucs in read plot

            # render nucs in sv plot
            return self.render_nucs()
    return False

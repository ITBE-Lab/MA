from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import FuncTickFormatter, TapTool, OpenURL, LabelSet, FixedTicker
from bokeh.models.callbacks import CustomJS
from bokeh.palettes import Plasma, Plasma256
from MA import *
from MSV import *
import math
from .util import *

def add_rectangle(self, seed_sample_size, read_id, rectangle, fill, read_ambiguous_reg_dict, end_column,
                  category_counter, use_dp):
    
    if rectangle.x_axis.size != 0:
        seed_sample_size /= rectangle.x_axis.size
    # if
    self.read_plot_rects[read_id]["l"].append(rectangle.x_axis.start)
    self.read_plot_rects[read_id]["b"].append(rectangle.y_axis.start)
    self.read_plot_rects[read_id]["r"].append(rectangle.x_axis.start + rectangle.x_axis.size)
    self.read_plot_rects[read_id]["t"].append(rectangle.y_axis.start + rectangle.y_axis.size)
    self.read_plot_rects[read_id]["f"].append(fill)
    self.read_plot_rects[read_id]["c"].append("lightgrey")
    self.read_plot_rects[read_id]["s"].append(seed_sample_size)
    self.read_plot_rects[read_id]["dp"].append(use_dp)
    if use_dp and len(self.read_ids) <= self.do_compressed_seeds:
        read_ambiguous_reg_dict["l"].append(rectangle.x_axis.start)
        read_ambiguous_reg_dict["b"].append(category_counter - 0.5)
        read_ambiguous_reg_dict["r"].append(rectangle.x_axis.start + rectangle.x_axis.size)
        read_ambiguous_reg_dict["t"].append(category_counter + len(end_column) - 0.5)
        read_ambiguous_reg_dict["f"].append("lightgrey")
        read_ambiguous_reg_dict["s"].append(seed_sample_size)

def render_reads(self, render_all=False):
    seeder = BinarySeeding(self.params)
    with self.measure("SvJumpsFromSeeds"):
        jumps_from_seeds = SvJumpsFromSeeds(self.params, self.pack)
    read_dict = {
        "center": [],
        "r_id": [],
        "r_name": [],
        "size": [],
        "q": [],
        "r": [],
        "l": [],
        "idx": [],
        "c": [],
        "f": [],
        "layer": [],
        "parlindrome": [],
        "x": [],
        "y": [],
        "category": []
    }
    # create a column data source for the read plot...
    read_plot_dict = {
        "center": [],
        "r_id": [],
        "r_name": [],
        "size": [],
        "q": [],
        "r": [],
        "l": [],
        "idx": [],
        "c": [],
        "f": [],
        "layer": [],
        "parlindrome": [],
        "x": [],
        "y": [],
        "category": []
    }
    read_ambiguous_reg_dict = {
        "l": [],
        "r": [],
        "t": [],
        "b": [],
        "f": [],
        "s": []
    }
    read_id_n_cols = {}
    col_ids = []
    all_col_ids = []
    category_counter = 0
    end_column = None

    with self.measure("computing seeds"):
        if self.do_render_seeds:
            try:
                read_table = ReadTable(self.db_conn)
                all_seeds = []
                for read_id in sorted(self.read_ids, reverse=True):
                    self.read_plot.nuc_plot.nucs_by_r_id[read_id] = {"p": [], "c": [], "i": []}
                    self.read_plot_rects[read_id] = {"l": [], "b": [], "t": [], "r": [], "f":[], "s":[], "c":[], "dp": []}
                    read = read_table.get_read(read_id)
                    for y, nuc in enumerate(str(read)):
                        append_nuc_type(self.read_plot.nuc_plot.nucs_by_r_id[read_id], nuc, y, "p")
                    with self.measure("seeder.execute"):
                        segments = seeder.execute(self.fm_index, read)
                    # execute_helper is not threadsave 
                    with self.measure("jumps_from_seeds.cpp_module.execute_helper"):
                        helper_ret = jumps_from_seeds.cpp_module.execute_helper(segments, self.pack, self.fm_index, read)
                    seeds = helper_ret.seeds
                    layer_of_seeds = helper_ret.layer_of_seeds
                    rectangles = helper_ret.rectangles
                    parlindromes = helper_ret.parlindrome
                    fill_of_rectangles = helper_ret.rectangles_fill
                    seed_sample_sizes = helper_ret.rectangle_ambiguity
                    rectangle_used_dp = helper_ret.rectangle_used_dp
                    # for
                    with self.measure("main seed loop"):
                        seeds_n_idx = list(enumerate(sorted([(x, y, z, read_id, read.name) for x, y, z in zip(seeds,
                                                                                                        layer_of_seeds,
                                                                                                        parlindromes)],
                                                            key=lambda x: x[0].start)))
                        if len(self.read_ids) <= self.do_compressed_seeds:
                            end_column = []
                            max_seed_size = max(seed.size for seed in seeds)
                            sorted_for_main_loop = sorted(seeds_n_idx, key=lambda x: x[1][0].start_ref)
                            for idx, (seed, layer, parlindrome, read_id, read_name) in sorted_for_main_loop:
                                add_seed(seed, read_dict, max_seed_size, end_column, all_col_ids, category_counter,
                                            parlindrome, layer, read_id, idx, read_name)
                            if (len(end_column)-1) % 2 == 0:
                                # prevent forming of float if possible (to stay javascript compatible)
                                curr_col_id = category_counter + (len(end_column)-1)//2
                            else:
                                curr_col_id = category_counter + (len(end_column)-1)/2
                            col_ids.append(curr_col_id)
                        else:
                            all_seeds.extend(seeds_n_idx)

                        for rectangle, fill, seed_sample_size, use_dp in zip(rectangles, fill_of_rectangles,
                                                                    seed_sample_sizes, rectangle_used_dp):
                            self.add_rectangle(seed_sample_size, read_id, rectangle, fill, read_ambiguous_reg_dict,
                                            end_column, category_counter, use_dp)

                        if len(self.read_ids) <= self.do_compressed_seeds:
                            category_counter += len(end_column) + 2
                            read_id_n_cols[curr_col_id] = read_id

                if len(self.read_ids) > self.do_compressed_seeds:
                    end_column = []
                    if len(all_seeds) > 0:
                        max_seed_size = max(seed[1][0].size for seed in all_seeds)
                        sorted_for_main_loop = sorted(all_seeds, key=lambda x: (x[1][0].start_ref, x[1][3]))
                        for idx, (seed, layer, parlindrome, read_id, read_name) in sorted_for_main_loop:
                            add_seed(seed, read_dict, max_seed_size, end_column, all_col_ids, category_counter,
                                        parlindrome, layer, read_id, idx, read_name)
                        category_counter += len(end_column)
            except Exception as e:
                print(e)

    def callback():
        if len(read_dict["c"]) < self.get_max_num_ele() or render_all:
            with self.measure("rendering seeds"):
                # render ambiguous regions on top and left
                self.seed_plot.ambiguous_regions.data = read_ambiguous_reg_dict

                # render seeds on top and left
                self.seed_plot.seeds.data = read_dict
                if self.seed_plot.left_plot.x_range.start == 0 and self.seed_plot.left_plot.x_range.end == 0:
                    self.seed_plot.left_plot.x_range.start = -1
                    self.seed_plot.left_plot.x_range.end = category_counter

                if len(self.read_ids) <= self.do_compressed_seeds:
                    self.seed_plot.left_plot.xaxis.ticker = FixedTicker(ticks=col_ids)
                    self.seed_plot.bottom_plot.yaxis.ticker = FixedTicker(ticks=col_ids)
                    self.seed_plot.left_plot.xaxis.formatter = FuncTickFormatter(
                        args={"read_id_n_cols": read_id_n_cols},
                        code="""
                                if(!tick in read_id_n_cols)
                                    return "";
                                return read_id_n_cols[tick];
                            """)
                    self.seed_plot.bottom_plot.yaxis.formatter = FuncTickFormatter(
                        args={"read_id_n_cols": read_id_n_cols},
                        code="""
                                if(!tick in read_id_n_cols)
                                    return "";
                                return read_id_n_cols[tick];
                            """)
                    self.seed_plot.left_plot.xaxis.axis_label = "Read Id"
                    self.seed_plot.bottom_plot.yaxis.axis_label = "Read Id"
                else:
                    self.seed_plot.left_plot.xaxis.ticker = []
                    self.seed_plot.left_plot.xaxis.axis_label = "compressed seeds"
                    self.seed_plot.bottom_plot.yaxis.ticker = []
                    self.seed_plot.bottom_plot.yaxis.axis_label = "compressed seeds"
                self.seed_plot.left_plot.xgrid.ticker = FixedTicker(ticks=all_col_ids)
                self.seed_plot.bottom_plot.ygrid.ticker = FixedTicker(ticks=all_col_ids)

                self.seed_plot.update_selection(self)

    self.curdoc.add_next_tick_callback(callback)
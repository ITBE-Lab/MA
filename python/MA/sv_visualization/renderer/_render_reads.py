from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import FuncTickFormatter, TapTool, OpenURL, LabelSet, FixedTicker
from bokeh.models.callbacks import CustomJS
from bokeh.palettes import Plasma, Plasma256
from MA import *
import math
from .util import *


def add_seed(self, seed, read_dict, max_seed_size, end_column, all_col_ids, category_counter, parlindrome, layer,
             read_id, read_plot_dict, idx):
    seed_size = seed.size - 1
    if seed.on_forward_strand:
        read_dict["center"].append(seed.start_ref + seed.size/2)
        read_dict["r"].append(seed.start_ref)
        read_dict["x"].append(
            [seed.start_ref, seed.start_ref+seed.size])
        curr_end = seed.start_ref + seed_size + 1
        curr_start = seed.start_ref
    else:
        read_dict["center"].append(seed.start_ref - seed.size/2 + 1)
        read_dict["r"].append(seed.start_ref - seed.size + 1)
        read_dict["x"].append(
            [seed.start_ref + 1, seed.start_ref - seed.size + 1])
        curr_end = seed.start_ref
        curr_start = seed.start_ref - seed.size
    curr_column = 0
    while curr_column < len(end_column):
        add_dist = 100
        if end_column[curr_column][1] == read_id:
            add_dist = 0
        if curr_start > end_column[curr_column][0] + add_dist:
            break
        else:
            curr_column += 1
    if curr_column >= len(end_column):
        end_column.append( None )
        all_col_ids.append(curr_column + category_counter)
    end_column[curr_column] = (curr_end, read_id)
    read_dict["y"].append([seed.start, seed.start+seed.size])
    read_dict["c"].append("lightgrey")
    read_dict["r_id"].append(read_id)
    read_dict["size"].append(seed.size)
    read_dict["l"].append(seed.size)
    read_dict["q"].append(seed.start)
    read_dict["idx"].append(idx)
    read_dict["layer"].append(layer)
    read_dict["parlindrome"].append(parlindrome)
    read_dict["f"].append(seed.on_forward_strand)
    read_dict["category"].append(
        category_counter + curr_column)

def add_rectangle(self, seed_sample_size, read_id, rectangle, fill, read_plot_rects,
                  read_ambiguous_reg_dict, end_column, category_counter):
    color = Plasma256[min(seed_sample_size, 255)]
    # if
    read_plot_rects[read_id]["l"].append(rectangle.x_axis.start)
    read_plot_rects[read_id]["b"].append(rectangle.y_axis.start)
    read_plot_rects[read_id]["r"].append(rectangle.x_axis.start + rectangle.x_axis.size)
    read_plot_rects[read_id]["t"].append(rectangle.y_axis.start + rectangle.y_axis.size)
    read_plot_rects[read_id]["f"].append(fill)
    read_plot_rects[read_id]["c"].append(color)
    read_plot_rects[read_id]["s"].append(seed_sample_size)
    if seed_sample_size > 10 and not self.do_compressed_seeds:
        read_ambiguous_reg_dict["l"].append(rectangle.x_axis.start)
        read_ambiguous_reg_dict["b"].append(category_counter - 0.5)
        read_ambiguous_reg_dict["r"].append(rectangle.x_axis.start + rectangle.x_axis.size)
        read_ambiguous_reg_dict["t"].append(category_counter + len(end_column) - 0.5)
        read_ambiguous_reg_dict["f"].append(fill)
        read_ambiguous_reg_dict["s"].append(seed_sample_size)

def render_reads(self):
    seeder = BinarySeeding(self.params)
    with self.measure("SvJumpsFromSeeds"):
        jumps_from_seeds = libMA.SvJumpsFromSeeds(self.params, -1, self.sv_db, self.pack)
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
        "parlindrome": [],
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
    read_plot_rects = {}  # dict of {"l": [], "b": [], "t": [], "r": [], "f":[], "s":[], "c":[]}
    end_column = None

    with self.measure("computing seeds"):
        if self.do_render_seeds:
            all_seeds = []
            for read_id in sorted(self.read_ids, reverse=True):
                self.read_plot.nuc_plot.nucs_by_r_id[read_id] = {"p": [], "c": [], "i": []}
                read_plot_rects[read_id] = {"l": [], "b": [], "t": [], "r": [], "f":[], "s":[], "c":[]}
                read = self.sv_db.get_read(read_id)
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
                if self.render_mems and self.selected_read_id == read_id:
                    hash_map_seeder = HashMapSeeding(self.params)
                    hash_map_seeder.cpp_module.seed_size = 9
                    query_section = NucSeq(str(read)[int(self.seed_plot_y_s):int(self.seed_plot_y_e)])
                    ref_section = self.pack.extract_from_to(int(self.xs), int(self.xe))
                    all_k_mers = hash_map_seeder.execute(query_section, ref_section)
                    all_mems = SeedLumping(self.params).execute(all_k_mers)
                    filter_module = FilterToUnique(self.params)
                    filter_module.cpp_module.num_mm = 0
                    #filtered_mems = filter_module.execute(all_mems, query_section, ref_section)
                    filtered_mems = all_mems
                    for idx in range(len(filtered_mems)):
                        filtered_mems[idx].start += int(self.seed_plot_y_s)
                        filtered_mems[idx].start_ref += int(self.xs)
                        layer_of_seeds.append(-1)
                    seeds.extend(filtered_mems)
                # for
                with self.measure("main seed loop"):
                    seeds_n_idx = list(enumerate(sorted([(x, y, z, read_id) for x, y, z in zip(seeds, layer_of_seeds,
                                                                                      parlindromes)],
                                                        key=lambda x: x[0].start)))
                    if not self.do_compressed_seeds:
                        end_column = []
                        max_seed_size = max(seed.size for seed in seeds)
                        sorted_for_main_loop = sorted(seeds_n_idx, key=lambda x: x[1][0].start_ref)
                        for idx, (seed, layer, parlindrome, read_id) in sorted_for_main_loop:
                            self.add_seed(seed, read_dict, max_seed_size, end_column, all_col_ids, category_counter,
                                        parlindrome, layer, read_id, read_plot_dict, idx)
                        if (len(end_column)-1) % 2 == 0:
                            # prevent forming of float if possible (to stay javascript compatible)
                            curr_col_id = category_counter + (len(end_column)-1)//2
                        else:
                            curr_col_id = category_counter + (len(end_column)-1)/2
                        col_ids.append(curr_col_id)
                    else:
                        all_seeds.extend(seeds_n_idx)

                    for rectangle, fill, seed_sample_size in zip(rectangles, fill_of_rectangles, seed_sample_sizes):
                        add_rectangle(self, seed_sample_size, read_id, rectangle, fill,
                                    read_plot_rects, read_ambiguous_reg_dict, end_column, category_counter)

                    if not self.do_compressed_seeds:
                        category_counter += len(end_column) + 2
                        read_id_n_cols[curr_col_id] = read_id

            if self.do_compressed_seeds:
                end_column = []
                if len(all_seeds) > 0:
                    max_seed_size = max(seed[1][0].size for seed in all_seeds)
                    sorted_for_main_loop = sorted(all_seeds, key=lambda x: (x[1][0].start_ref, x[1][3]))
                    for idx, (seed, layer, parlindrome, read_id) in sorted_for_main_loop:
                        self.add_seed(seed, read_dict, max_seed_size, end_column, all_col_ids, category_counter,
                                    parlindrome, layer, read_id, read_plot_dict, idx)
                    category_counter += len(end_column)

    with self.measure("rendering seeds"):
        if len(read_dict["c"]) < self.get_max_num_ele() or self.render_mems:
            # render ambiguous regions on top and left
            self.seed_plot.ambiguous_regions.data = read_ambiguous_reg_dict

            # render seeds on top and left
            self.seed_plot.seeds.data = read_dict
            if self.seed_plot.left_plot.x_range.start == 0 and self.seed_plot.left_plot.x_range.end == 0:
                self.seed_plot.left_plot.x_range.start = -1
                self.seed_plot.left_plot.x_range.end = category_counter

            if not self.do_compressed_seeds:
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
            else:
                self.seed_plot.left_plot.xaxis.ticker = []
                self.seed_plot.left_plot.xaxis.axis_label = "compressed seeds"
                self.seed_plot.bottom_plot.yaxis.ticker = []
                self.seed_plot.bottom_plot.yaxis.axis_label = "compressed seeds"
            self.seed_plot.left_plot.xgrid.ticker = FixedTicker(ticks=all_col_ids)
            self.seed_plot.bottom_plot.ygrid.ticker = FixedTicker(ticks=all_col_ids)

            self.seed_plot.update_selection(self)
    return False

from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import FuncTickFormatter, TapTool, OpenURL, LabelSet, FixedTicker
from bokeh.models.callbacks import CustomJS
from bokeh.palettes import Plasma, Plasma256
from MA import *
from MSV import *
import math
from .util import *
import sys, traceback

def add_rectangle(self, seed_sample_size, read_id, rectangle, fill, read_ambiguous_reg_dict, end_column_len,
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
        read_ambiguous_reg_dict["t"].append(category_counter + end_column_len - 0.5)
        read_ambiguous_reg_dict["f"].append("lightgrey")
        read_ambiguous_reg_dict["s"].append(seed_sample_size)

def render_reads(self, render_all=False):
    print("rendering reads")
    self.params.by_name("Fixed SoC Width").set(50)
    self.params.by_name("Max Size Reseed").set(2000)
    """
    seeding_module = MMFilteredSeeding(self.params, 300)
    seed_lumper = SeedLumping(self.params)
    soc_module = StripOfConsiderationSeeds(self.params)
    soc_filter = GetAllFeasibleSoCs(self.params, 100)
    with self.measure("SvJumpsFromSeeds"):
        jumps_from_seeds = SvJumpsFromSeeds(self.params, self.pack)
    """
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
                seed_info_list, rectangles_info_list, reads_n_cols, reads = seedDisplaysForReadIds(self.params, 
                                                        self.db_pool, list(self.read_ids), self.pack,
                                                        self.mm_index, self.mm_counter,
                                                        len(self.read_ids) > self.do_compressed_seeds, 3)

    with self.measure("render seeds"):
                for seed_info in seed_info_list:
                    read_dict["r_id"].append(seed_info.iReadId)
                    read_dict["r_name"].append(seed_info.sReadName)
                    read_dict["size"].append(seed_info.uiSize)
                    read_dict["l"].append(seed_info.uiL)
                    read_dict["q"].append(seed_info.uiQ)
                    read_dict["idx"].append(seed_info.uiSeedOrderOnQuery)
                    read_dict["layer"].append(seed_info.uiLayer)
                    read_dict["parlindrome"].append(seed_info.bParlindrome)
                    read_dict["f"].append(seed_info.bOnForward)
                    read_dict["category"].append(seed_info.uiCategory)
                    
                    read_dict["center"].append(seed_info.fCenter)
                    read_dict["r"].append(seed_info.uiR)
                    read_dict["x"].append([*seed_info.xX])
                    read_dict["y"].append([*seed_info.xY])
                    read_dict["c"].append("lightgrey")

                for read in reads:
                    self.read_plot.nuc_plot.nucs_by_r_id[read.id] = {"p": [], "c": [], "i": []}
                    self.read_plot_rects[read.id] = {"l": [], "b": [], "t": [], "r": [], "f":[], "s":[], "c":[], "dp": []}
                    for y, nuc in enumerate(str(read)):
                        append_nuc_type(self.read_plot.nuc_plot.nucs_by_r_id[read.id], nuc, y, "p")
                for x, y in reads_n_cols:
                    read_id_n_cols[x] = y
                for r_i in rectangles_info_list:
                    for rectangle, fill, seed_sample_size, use_dp in zip(r_i.vRectangles, r_i.vRectangleFillPercentage,
                                                                r_i.vRectangleReferenceAmbiguity, r_i.vRectangleUsedDp):
                        self.add_rectangle(seed_sample_size, r_i.iReadId, rectangle, fill, read_ambiguous_reg_dict,
                                        r_i.uiEndColumnSize, r_i.uiCategory, use_dp)

    def callback():
        if len(read_dict["c"]) < self.get_max_num_ele() or render_all or True:
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

    print("done rendering reads")
    self.curdoc.add_next_tick_callback(callback)
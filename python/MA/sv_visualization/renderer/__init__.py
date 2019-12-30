from MA import AnalyzeRuntimes, ParameterSetManager
import datetime
from .visual_elements.main_plot import *
from .visual_elements.seed_plot import *
from .visual_elements.nuc_plot import *
from .visual_elements.read_plot import *
from .visual_elements.widgets import *

class Renderer():
    def __init__(self):
        self.main_plot = MainPlot(self)
        self.nuc_plot = NucPlot(self.main_plot)
        self.seed_plot = SeedPlot(self.main_plot, self)
        self.read_plot = ReadPlot(self.nuc_plot, self)
        self.widgets = Widgets(self)
        self.pack = None
        self.fm_index = None
        self.sv_db = None
        self.w = None
        self.h = None
        self.params = ParameterSetManager()
        self.quads = []
        self.read_ids = set() # @todo make two dicts? read ids new and read ids old?
        self.give_up_factor = 1000
        self.analyze = AnalyzeRuntimes()
        self.do_render_seeds = True
        self.do_compressed_seeds = 30
        self.xs = 0
        self.ys = 0
        self.xe = 0
        self.ye = 0
        self.redered_everything = False
        self.do_render_call_jumps_only = True
        self.render_area_factor = 1 # @todo make this adjustable 
        # number of reads needs to be smaller than max_num_elements / read_penalty_factor to be rendered
        self.read_penalty_factor = 10
        self.selected_read_id = None
        self.selected_seed_id = None
        self.selected_call_id = set()
        self.selected_jump_id = set()
        self.read_plot_rects = {}  # dict of {"l": [], "b": [], "t": [], "r": [], "f":[], "s":[], "c":[]}

    def get_run_id(self):
        if self.widgets.run_id_dropdown.value is None:
            return -1
        return int(self.widgets.run_id_dropdown.value)

    def get_gt_id(self):
        if self.widgets.ground_truth_id_dropdown.value is None:
            return -1
        return int(self.widgets.ground_truth_id_dropdown.value)

    def get_min_score(self):
        return self.widgets.score_slider.value

    def get_max_num_ele(self):
        return self.widgets.max_elements_slider.value

    def measure(self, name):
        class MeasureHelper:
            def __init__(self, name, analyze):
                self.name = name
                self.start = None
                self.analyze = analyze
            def __enter__(self):
                self.start = datetime.datetime.now()
            def __exit__(self, exc_type, exc_val, exc_tb):
                end = datetime.datetime.now()
                delta = end - self.start
                self.analyze.register(name, delta.total_seconds(), False, lambda x: x)
        return MeasureHelper(name, self.analyze)

    def reset_runtimes(self):
        self.analyze.reset()

    def reset_cds(self):
        self.main_plot.reset_cds()
        self.nuc_plot.reset_cds()
        self.seed_plot.reset_cds()
        self.read_plot.reset_cds()

    # imported methdos
    from ._render import render
    from ._setup import setup
    from ._render_overview import render_overview
    from ._render_calls import render_calls
    from ._render_jumps import render_jumps
    from ._render_reads import add_seed
    from ._render_reads import add_rectangle
    from ._render_reads import render_reads
    from ._render_nucs import render_nucs

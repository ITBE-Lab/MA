from MA import AnalyzeRuntimes, ParameterSetManager
from MSV import *
import datetime
from .visual_elements.main_plot import *
from .visual_elements.seed_plot import *
from .visual_elements.nuc_plot import *
from .visual_elements.read_plot import *
from .visual_elements.widgets import *
from bokeh.plotting import curdoc

class Renderer():
    def __init__(self):
        self.main_plot = MainPlot(self)
        self.nuc_plot = NucPlot(self.main_plot)
        self.seed_plot = SeedPlot(self.main_plot, self)
        self.read_plot = ReadPlot(self.nuc_plot, self)
        self.widgets = Widgets(self)
        self.pack = None
        self.fm_index = None
        self.mm_counter = HashCounters()
        self.db_conn = None
        self.db_conn_2 = None
        self.db_pool = None
        self.w = None
        self.h = None
        self.params = ParameterSetManager()
        self.params.by_name("Min Size Edge").set(10)
        self.params.by_name("Maximal Ambiguity SV").set(100)
        self.params.by_name("Min NT in SoC").set(25)
        self.params.by_name("Rectangular SoC").set(False)
        self.params.by_name("Fixed SoC Width").set(100)
        self.params.by_name("Max Size Reseed").set(2000)
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
        self.do_render_call_jumps_only = False
        self.render_area_factor = 1 # @todo make this adjustable 
        # number of reads needs to be smaller than max_num_elements / read_penalty_factor to be rendered
        self.read_penalty_factor = 100
        self.selected_read_id = None
        self.selected_seed_id = None
        self.selected_call_id = set()
        self.selected_jump_id = set()
        self.read_plot_rects = {}  # dict of {"l": [], "b": [], "t": [], "r": [], "f":[], "s":[], "c":[]}
        self.cached_global_overview = None
        self.cached_overview_min_score = None
        self.cached_overview_max_render = None
        self.global_overview_threshold = 0.2
        self._do_overview_cache = False
        self.curdoc = curdoc()

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
        self.main_plot.reset_cds(self)
        self.nuc_plot.reset_cds()
        self.seed_plot.reset_cds()
        self.read_plot.reset_cds()

    def get_global_overview(self):
        if self.cached_global_overview is None or self.cached_overview_min_score != self.get_min_score() or \
                self.cached_overview_max_render != self.get_max_num_ele():
            div = int(math.sqrt(self.get_max_num_ele()/10))
            self.cached_global_overview = get_call_overview(self.db_pool, self.pack,
                                        self.get_run_id(), self.get_min_score(),
                                        0, 0, self.pack.unpacked_size_single_strand,
                                        self.pack.unpacked_size_single_strand,
                                        self.pack.unpacked_size_single_strand//div,
                                        self.pack.unpacked_size_single_strand//div, self.give_up_factor)
            self.cached_overview_max_render = self.get_max_num_ele()
            self.cached_overview_min_score = self.get_min_score()
        return self.cached_global_overview

    def do_overview_cache(self):
        p_s = self.pack.unpacked_size_single_strand
        return self._do_overview_cache and self.w*3 * self.h*3 >= p_s * p_s * self.global_overview_threshold

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

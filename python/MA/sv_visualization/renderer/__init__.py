from MA import AnalyzeRuntimes
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
        self.seed_plot = SeedPlot(self.main_plot)
        self.read_plot = ReadPlot(self.nuc_plot)
        self.widgets = Widgets()
        self.pack = None
        self.fm_index = None
        self.sv_db = None
        self.w = None
        self.h = None
        self.params = None
        self.quads = []
        self.read_ids = set()
        self.give_up_factor = 1000
        self.analyze = AnalyzeRuntimes()
        self.do_render_seeds = True
        self.do_compressed_seeds = True
        self.xs = 0
        self.ys = 0
        self.xe = 0
        self.ye = 0
        self.redered_everything = False
        self.render_area_factor = 2

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
                self.analyze.register(name, delta.total_seconds(), lambda x: x)
        return MeasureHelper(name, self.analyze)


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

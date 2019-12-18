from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import FuncTickFormatter, TapTool, OpenURL, LabelSet, FixedTicker
from bokeh.models.callbacks import CustomJS
from MA import *
import math
from .util import *

def render_nucs(self):
    l_plot_nucs = {"p": [], "c": [], "i": []}

    nuc_seq = self.pack.extract_from_to(max(int(self.ys - self.h), 0),
                                        min(int(self.ye + self.h + 1), self.pack.unpacked_size_single_strand))
    for y_add, nuc in enumerate(str(nuc_seq)):
        append_nuc_type(l_plot_nucs, nuc, int(self.ys - self.h) + y_add, "p")
    self.nuc_plot.left_nucs.data = l_plot_nucs

    d_plot_nucs = {"p": [], "c": [], "i": []}
    nuc_seq = self.pack.extract_from_to(max(int(self.xs - self.w), 0),
                                        min(int(self.xe + self.w + 1), self.pack.unpacked_size_single_strand))
    for x_add, nuc in enumerate(str(nuc_seq)):
        append_nuc_type(d_plot_nucs, nuc, int(self.xs - self.w) + x_add, "p")
    self.nuc_plot.bottom_nucs.data = d_plot_nucs

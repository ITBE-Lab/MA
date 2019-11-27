from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import FuncTickFormatter, TapTool, OpenURL, LabelSet, FixedTicker
from bokeh.models.callbacks import CustomJS
from MA import *
import math
from .util import *

def render_nucs(self):
    l_plot_nucs = {"y": [], "c": [], "i": []}

    nuc_seq = self.pack.extract_from_to(
        int(self.ys - self.h), int(self.ye + self.h + 1))
    for y_add, nuc in enumerate(str(nuc_seq)):
        append_nuc_type(l_plot_nucs, nuc, int(self.ys - self.h) + y_add, "y")
    self.l_plot[0].rect(x=0.5, y="y", width=1, height=1, fill_color="c", line_width=0,
                        source=ColumnDataSource(l_plot_nucs), name="hover4")
    self.d_read_plot.rect(x="y", y=0.5, width=1, height=1, fill_color="c", line_width=0,
                        source=ColumnDataSource(l_plot_nucs), name="hover4")
    d_plot_nucs = {"x": [], "c": [], "i": []}
    nuc_seq = self.pack.extract_from_to(
        int(self.xs - self.w), int(self.xe + self.w + 1))
    for x_add, nuc in enumerate(str(nuc_seq)):
        append_nuc_type(d_plot_nucs, nuc, int(self.xs - self.w) + x_add, "x")
    self.d_plot[0].rect(x="x", y=0.5, width=1, height=1, fill_color="c", line_width=0,
                        source=ColumnDataSource(d_plot_nucs), name="hover4")
    self.d_read_plot.rect(x="x", y=0.5, width=1, height=1, fill_color="c", line_width=0,
                        source=ColumnDataSource(d_plot_nucs), name="hover4")
    return True
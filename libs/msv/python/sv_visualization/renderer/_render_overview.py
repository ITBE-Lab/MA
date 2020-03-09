from bokeh.models import TapTool, OpenURL
from bokeh.plotting import ColumnDataSource
import math
from MA import *
from .util import *


def render_overview(self):
    #self.plot.grid.visible = False
    with self.measure("get_call_overview"):
        if self.do_overview_cache():
            rect_vec = self.get_global_overview()
        else:
            div = int(math.sqrt(self.get_max_num_ele()/10))
            rect_vec = get_call_overview(self.db_pool, self.pack, self.get_run_id(), self.get_min_score(),
                                            int(self.xs - self.w),
                                            int(self.ys - self.h),
                                            self.w*3, self.h*3,
                                            self.w//div, self.h//div, self.give_up_factor)

    cds = {
        'x': [],
        'y': [],
        'w': [],
        'h': [],
        'c': [],
        'f': [],
        't': [],
        'i': []
    }
    max_ = max(*[r.c for r in rect_vec], 0)
    names = self.pack.contigNames()
    for rect in rect_vec:
        cds["x"].append(rect.x)
        cds["y"].append(rect.y)
        cds["w"].append(rect.x + rect.w)
        cds["h"].append(rect.y + rect.h)
        cds["c"].append(format(light_spec_approximation(rect.c/max_)))
        cds["f"].append(names[rect.i])
        cds["t"].append(names[rect.j])
        cds["i"].append(str(rect.c))
    self.main_plot.overview_quad.data = cds

    self.analyze.analyze()

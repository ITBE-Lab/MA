from MA import *
from MSV import *
from threading import Thread
from bokeh.models.tickers import FixedTicker

def render(self, render_all=False):
    print("rendering")
    self.widgets.show_spinner(self)
    self.reset_runtimes()
    self.reset_cds()
    if self.xs < 0:
        self.xs = 0
    if self.ys < 0:
        self.ys = 0
    if self.xe < 0:
        self.xe = 0
    if self.ye < 0:
        self.ye = 0
    self.w = int(self.xe - self.xs)
    self.h = int(self.ye - self.ys)

    s = max(min(self.xs - self.w, self.ys - self.h), 0)
    e = min(max(self.xe + self.w, self.ye + self.h), self.pack.unpacked_size_single_strand)


    # plot diagonal; we need s and e since too large lines sometimes do not render...
    self.main_plot.diagonal_line.data = {"xs":[s, e, float("NaN")], "ys":[s, e, float("NaN")]}

    def blocking_task():
        if not self.widgets.run_id_dropdown.value is None:
            with self.measure("get_call_overview_area"):
                if self.do_overview_cache():
                    num_ele = self.get_max_num_ele() + 1
                else:
                    num_ele = get_call_overview_area(self.db_conn, self.pack, self.get_run_id(),
                                                    self.get_min_score(),
                                                    int(self.xs - self.w), int(self.ys - self.h), self.w*3, self.h*3,
                                                    self.get_max_num_ele() + 1)
            if num_ele > self.get_max_num_ele() and not render_all:
                self.render_overview()
            else:
                self.render_calls(render_all)

        print("done")
        self.widgets.hide_spinner(self)

    thread = Thread(target=blocking_task)
    thread.start()

from MA import *


def render(self):
    # genome outline
    self.plot.quad(left=0, bottom=0, right=self.pack.unpacked_size_single_strand,
                   top=self.pack.unpacked_size_single_strand,
                   fill_alpha=0, line_color="black", line_width=3)

    if not self.sv_db.run_exists(self.run_id):
        return True

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
    e = min(max(self.xe + self.w, self.ye + self.h),
            self.pack.unpacked_size_single_strand)
    # plot diagonal; we need s and e since too large lines sometimes do not render...
    self.plot.line(x=[s, e], y=[s, e], line_color="black", line_width=3)

    if libMA.get_call_overview_area(self.sv_db, self.pack, self.run_id, self.min_score, int(self.xs - self.w), int(self.ys - self.h), self.w*3, self.h*3) > self.max_num_ele:
        return self.render_overview()
    else:
        return self.render_calls()

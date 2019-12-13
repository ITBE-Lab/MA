from MA import AnalyzeRuntimes
import datetime

class Renderer():
    def __init__(self, plot, l_plot, d_plot, xs, xe, ys, ye, pack, fm_index, sv_db, run_id, ground_truth_id, min_score,
                 max_num_ele, dataset_name, active_tools, radio_group, read_plot, selected_read_id, l_read_plot,
                 d_read_plot, render_mems, seed_plot_y_s, seed_plot_y_e, index_prefix):
        self.plot = plot
        self.l_plot = l_plot
        self.d_plot = d_plot
        self.xs = xs
        self.xe = xe
        self.ys = ys
        self.ye = ye
        self.pack = pack
        self.fm_index = fm_index
        self.sv_db = sv_db
        self.run_id = run_id
        self.ground_truth_id = ground_truth_id
        self.min_score = min_score
        self.max_num_ele = max_num_ele
        self.dataset_name = dataset_name
        self.active_tools = active_tools
        self.radio_group = radio_group
        self.read_plot = read_plot
        self.selected_read_id = selected_read_id
        self.index_prefix = index_prefix
        self.render_mems = render_mems
        self.seed_plot_y_s = seed_plot_y_s
        self.seed_plot_y_e = seed_plot_y_e
        self.w = None
        self.h = None
        self.params = None
        self.l_read_plot = l_read_plot
        self.d_read_plot = d_read_plot
        self.quads = []
        self.read_ids = set()
        self.give_up_factor = 1000
        self.analyze = AnalyzeRuntimes()
        self.do_render_seeds = True
        self.do_compressed_seeds = True

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
    from ._render_overview import render_overview
    from ._render_calls import render_calls
    from ._render_jumps import render_jumps
    from ._render_reads import add_seed
    from ._render_reads import add_rectangle
    from ._render_reads import render_reads
    from ._render_nucs import render_nucs

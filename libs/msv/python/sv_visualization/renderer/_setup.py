from MSV import *
from renderer.util import *
from bokeh.models.tickers import CompositeTicker, FixedTicker
from bokeh.models import FuncTickFormatter
import os.path
import json

JSON_PREFIX = "/MAdata/sv_datasets/"


def setup(self):
    dataset_name = self.widgets.file_input.value
    if not dataset_name is None and os.path.isfile(JSON_PREFIX + dataset_name + "/info.json"):
        with open(JSON_PREFIX + dataset_name + "/info.json", "r") as json_file:
            json_info_file = json.loads(json_file.read(), object_hook=decode)
        ref_genome = json_info_file["reference_path"] + "/ma/genome"
        print("Genome:", json_info_file["reference_path"] + "/ma/genome")

        self.pack = Pack()
        self.pack.load(ref_genome)
        self.fm_index = FMIndex()
        self.fm_index.load(ref_genome)
        self.mm_index = MinimizerIndex(ParameterSetManager(), self.pack.contigSeqs(), self.pack.contigNames())
        self.mm_index.set_max_occ(2)
        self.db_conn = DbConn(dataset_name)
        print("NUM THREADS", self.params.get_num_threads())
        self.db_pool = PoolContainer(self.params.get_num_threads() + 1, dataset_name)

        # chromosome lines
        xs = [*self.pack.contigStarts(), self.pack.unpacked_size_single_strand]
        self.main_plot.plot.x_range.end = self.pack.unpacked_size_single_strand
        self.main_plot.plot.y_range.end = self.pack.unpacked_size_single_strand
        self.main_plot.plot.x_range.start = 0
        self.main_plot.plot.y_range.start = 0
        self.main_plot.plot.grid.ticker = FixedTicker(ticks=xs)
        self.seed_plot.bottom_plot.grid.ticker = FixedTicker(ticks=xs)
        self.seed_plot.left_plot.grid.ticker = FixedTicker(ticks=xs)
        self.read_plot.plot.xgrid.ticker = FixedTicker(ticks=xs)

        self.main_plot.plot.grid.bounds = (0, self.pack.unpacked_size_single_strand)
        self.seed_plot.bottom_plot.grid.bounds = (0, self.pack.unpacked_size_single_strand)
        self.seed_plot.left_plot.grid.bounds = (0, self.pack.unpacked_size_single_strand)
        self.read_plot.plot.xgrid.bounds = (0, self.pack.unpacked_size_single_strand)

        #contig_ends = [x+s for x,s in zip(self.pack.contigStarts(), self.pack.contigLengths())]
        #self.seed_plot.left_plot.yaxis[0].formatter = FuncTickFormatter(
        #                args={"contig_starts": [*self.pack.contigStarts(), self.pack.unpacked_size_single_strand],
        #                    "genome_end":self.pack.unpacked_size_single_strand,
        #                    "contig_names": [*self.pack.contigNames()]},
        #                code="""
        #                        if(tick < 0 || tick >= genome_end)
        #                            return "n/a";
        #                        idx = 0;
        #                        while(contig_starts[idx + 1] < tick)
        #                            idx += 1;
        #                        return contig_names[idx] + ": " + (tick - contig_starts[idx]);
        #                    """)
        #reference_ticks_center = [x+s/2 for x,s in zip(self.pack.contigStarts(), self.pack.contigLengths())]
        #self.seed_plot.left_plot.yaxis.ticker.max_interval = min(*self.pack.contigLengths())
        #self.seed_plot.left_plot.yaxis.ticker.min_interval = 1
        #self.seed_plot.left_plot.yaxis[0].ticker = CompositeTicker(tickers=[
        #                                                        self.seed_plot.left_plot.yaxis.ticker,
        #                                                        f_ticker])

        self.seed_plot.bottom_plot.xaxis[0].formatter = self.seed_plot.left_plot.yaxis[0].formatter
        self.seed_plot.bottom_plot.xaxis[0].ticker = self.seed_plot.left_plot.yaxis[0].ticker

        self.read_plot.nuc_plot.bottom_plot.xaxis[0].formatter = self.seed_plot.left_plot.yaxis[0].formatter
        self.read_plot.nuc_plot.bottom_plot.xaxis[0].ticker = self.seed_plot.left_plot.yaxis[0].ticker


        self.xs = 0
        self.xe = self.pack.unpacked_size_single_strand
        self.ys = 0
        self.ye = self.pack.unpacked_size_single_strand

        menu = []
        self.widgets.run_id_dropdown.label="select run id here"
        self.widgets.ground_truth_id_dropdown.label="select ground truth id here"
        if not self.db_conn is None:
            run_table = SvCallerRunTable(self.db_conn)
            for run_id in run_table.getIds():
                text = run_table.getName(run_id) + " - " + run_table.getDate(run_id) + " - " + run_table.getDesc(run_id)
                text += " - " + str(SvCallTable(self.db_conn).num_calls(run_id, self.widgets.score_slider.value))
                menu.append((text, str(run_id)))
        self.widgets.run_id_dropdown.menu = menu
        self.widgets.ground_truth_id_dropdown.menu = menu

        self.render()
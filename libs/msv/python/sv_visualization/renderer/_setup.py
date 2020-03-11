from MSV import *
from renderer.util import *
import os.path
import json

JSON_PREFIX = "/MAdata/sv_datasets/"

def setup(self):
    dataset_name = self.widgets.file_input.value
    if not dataset_name is None and os.path.isfile(JSON_PREFIX + dataset_name + "/info.json"):
        with open(JSON_PREFIX + dataset_name + "/info.json", "r") as json_file:
            json_info_file = json.loads(json_file.read(), object_hook=decode)
        ref_genome = json_info_file["reference_path"] + "/ma/genome"

        self.pack = Pack()
        self.pack.load(ref_genome)
        self.fm_index = FMIndex()
        self.fm_index.load(ref_genome)
        self.db_conn = DbConn(dataset_name)
        print("NUM THREADS", self.params.get_num_threads())
        self.db_pool = PoolContainer(self.params.get_num_threads() + 1, dataset_name)

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
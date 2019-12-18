from MA import *
from renderer.util import *
import os.path
import json

def setup(self):
    dataset_name = self.widgets.file_input.value
    if os.path.isfile(dataset_name + "/info.json"):
        with open(dataset_name + "/info.json", "r") as json_file:
            json_info_file = json.loads(json_file.read(), object_hook=decode)
        ref_genome = json_info_file["reference_path"] + "/ma/genome"

        self.pack = Pack()
        self.pack.load(ref_genome)
        self.fm_index = FMIndex()
        self.fm_index.load(ref_genome)
        self.sv_db = SV_DB(dataset_name + "/svs.db", "open")

        self.xs = 0
        self.xe = self.pack.unpacked_size_single_strand
        self.ys = 0
        self.ye = self.pack.unpacked_size_single_strand

        menu = []
        self.widgets.run_id_dropdown.label="select run id here"
        self.widgets.ground_truth_id_dropdown.label="select ground truth id here"
        if not self.sv_db is None:
            jump_fetcher = libMA.SvCallerRunsFromDb(self.sv_db)
            while not jump_fetcher.eof():
                text = jump_fetcher.name() + " - " + \
                    self.sv_db.get_run_date(jump_fetcher.id()) + " - " + \
                    jump_fetcher.desc()
                text += " - " + \
                    str(self.sv_db.get_num_calls(jump_fetcher.id(), self.widgets.score_slider.value))
                menu.append((text, str(jump_fetcher.id())))
                jump_fetcher.next()
        self.widgets.run_id_dropdown.menu = menu
        self.widgets.ground_truth_id_dropdown.menu = menu

        self.render()
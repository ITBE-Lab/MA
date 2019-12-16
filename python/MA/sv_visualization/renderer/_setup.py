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

        self.render()
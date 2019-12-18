from bokeh.plotting import figure
from bokeh.models import Button, Slider
from bokeh.models.widgets import Dropdown, TextInput, RadioGroup, CheckboxGroup
from bokeh.events import ButtonClick
from ..util import *
from MA import *
import copy

class Widgets:
    def __init__(self, renderer):
        self.file_input = TextInput(value="/MAdata/sv_datasets/minimal", title="Dataset Name:")
        self.run_id_dropdown = Dropdown(label="select run id here", menu=[])
        self.ground_truth_id_dropdown = Dropdown(label="select ground truth id here", menu=[])
        self.score_slider = Slider(start=0, end=1, value=0, step=.1, callback_policy='mouseup', title="min score")
        self.max_elements_slider = Slider(start=1000, end=100000, value=25000, step=1000,
                                          callback_policy='mouseup', title="max render")
        self.range_link_radio = RadioGroup(labels=["Link read plot to x-range", "Link read plot to y-range"],
                                           active=0)
        self.render_mems_button = Button(label="render MEMs")
        self.delete_button = Button(label="Delete Dataset")

        self.file_input.on_change("value", lambda x,y,z: self.file_input_change(renderer))
        self.run_id_dropdown.on_change("value", lambda x,y,z: self.run_id_change(renderer))
        self.ground_truth_id_dropdown.on_change("value", lambda x,y,z: self.ground_id_change(renderer))
        self.score_slider.on_change("value_throttled", lambda x,y,z: self.slider_change(renderer))
        self.max_elements_slider.on_change("value_throttled", lambda x,y,z: self.slider_change(renderer))

        self.render_mems_button.on_event(ButtonClick, lambda x: self.render_mems_button_event(renderer))
        self.delete_button.on_event(ButtonClick, lambda x: self.delete_button_event(renderer))

    def file_input_change(self, renderer):
        renderer.setup()

    def run_id_change(self, renderer):
        self.run_id_dropdown.label = "Selected run: " + renderer.sv_db.get_run_name(int(self.run_id_dropdown.value)) + \
                                     " - " + self.run_id_dropdown.value
        self.score_slider.end = renderer.sv_db.get_max_score(int(self.run_id_dropdown.value))
        renderer.render()

    def ground_id_change(self, renderer):
        self.ground_truth_id_dropdown.label = "Selected ground truth: " + \
                                     renderer.sv_db.get_run_name(int(self.ground_truth_id_dropdown.value)) + \
                                     " - " + self.ground_truth_id_dropdown.value
        renderer.render()

    def slider_change(self, renderer):
        renderer.render()

    def render_mems_button_event(self, renderer):
        if not renderer.selected_read_id is None:
            read = renderer.sv_db.get_read(renderer.selected_read_id)

            seed_plot_y_s = max(renderer.read_plot.plot.y_range.start, 0)
            seed_plot_y_e = min(renderer.read_plot.plot.y_range.end, len(read))
            seed_plot_x_s = max(renderer.read_plot.plot.x_range.start, 0)
            seed_plot_x_e = min(renderer.read_plot.plot.x_range.end, renderer.pack.unpacked_size_single_strand)

            hash_map_seeder = HashMapSeeding(renderer.params)
            hash_map_seeder.cpp_module.seed_size = 9
            query_section = NucSeq(str(read)[int(seed_plot_y_s):int(seed_plot_y_e)])
            ref_section = renderer.pack.extract_from_to(int(seed_plot_x_s), int(seed_plot_x_e))
            all_k_mers = hash_map_seeder.execute(query_section, ref_section)
            all_mems = SeedLumping(renderer.params).execute(all_k_mers)
            filter_module = FilterToUnique(renderer.params)
            filter_module.cpp_module.num_mm = 0
            #filtered_mems = filter_module.execute(all_mems, query_section, ref_section)
            filtered_mems = all_mems
            seed_data_new = dict((key, []) for key in renderer.read_plot.seeds.data.keys())
            if len(filtered_mems) > 0:
                max_seed_size = max( seed.size for seed in filtered_mems )
                for idx in range(len(filtered_mems)):
                    filtered_mems[idx].start += int(seed_plot_y_s)
                    filtered_mems[idx].start_ref += int(seed_plot_x_s)

                    add_seed(filtered_mems[idx], seed_data_new, max_seed_size, [], [],
                                0, False, -1, renderer.selected_read_id, idx)
            renderer.read_plot.seeds.data = seed_data_new

    def delete_button_event(self, renderer):
        renderer.sv_db.delete_run(renderer.get_run_id())
        renderer.setup()

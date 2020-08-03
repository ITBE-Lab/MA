from bokeh.plotting import figure
from bokeh.models import Button, Slider
from bokeh.models.widgets import Dropdown, TextInput, RadioGroup, CheckboxGroup, Div
from bokeh.events import ButtonClick
from ..util import *
from MA import *
from MSV import *
import copy

class Widgets:
    def __init__(self, renderer):
        self.file_input = Dropdown(label="Select dataset", menu=[])
        for dataset_name in SQLDBInformer(DbConn({})).get_all_schemas():
            if dataset_name in set(["information_schema", "performance_schema", "sys", "defaultDB", "mysql"]):
                continue
            self.file_input.menu.append(dataset_name)
        self.run_id_dropdown = Dropdown(label="select run id here", menu=[])
        self.ground_truth_id_dropdown = Dropdown(label="select ground truth id here", menu=[])
        self.score_slider = Slider(start=0, end=1, value=0, step=.1, callback_policy='mouseup', title="min score")
        self.max_elements_slider = Slider(start=1000, end=100000, value=25000, step=1000,
                                          callback_policy='mouseup', title="max render")
        self.range_link_radio = RadioGroup(labels=["Link read plot to x-range", "Link read plot to y-range"],
                                           active=0)
        self.full_render_button = Button(label="render everything")
        self.render_mems_button = Button(label="render MEMs")
        self.delete_button = Button(label="Delete Dataset")

        self.file_input.on_change("value", lambda x,y,z: self.file_input_change(renderer))
        self.run_id_dropdown.on_change("value", lambda x,y,z: self.run_id_change(renderer))
        self.ground_truth_id_dropdown.on_change("value", lambda x,y,z: self.ground_id_change(renderer))
        self.score_slider.on_change("value_throttled", lambda x,y,z: self.slider_change(renderer))
        self.max_elements_slider.on_change("value_throttled", lambda x,y,z: self.slider_change(renderer))

        self.full_render_button.on_event(ButtonClick, lambda x: self.full_render(renderer))
        self.render_mems_button.on_event(ButtonClick, lambda x: self.render_mems_button_event(renderer))
        self.delete_button.on_event(ButtonClick, lambda x: self.delete_button_event(renderer))

        self.spinner_div = Div(text=html_file("spinner"), sizing_mode="scale_both", visible=False)

    def show_spinner(self, renderer):
        def callback():
            self.spinner_div.visible = True
        renderer.curdoc.add_next_tick_callback(callback)

    def hide_spinner(self, renderer):
        def callback():
            self.spinner_div.visible = False
        renderer.curdoc.add_next_tick_callback(callback)

    def file_input_change(self, renderer):
        self.file_input.label = "Selected dataset: " + self.file_input.value
        renderer.setup()

    def run_id_change(self, renderer):
        print("new run_id:", self.run_id_dropdown.value)
        renderer.cached_global_overview = None
        run_table = SvCallerRunTable(renderer.db_conn)
        self.run_id_dropdown.label = "Selected run: " + run_table.getName(int(self.run_id_dropdown.value)) + \
                                     " - " + self.run_id_dropdown.value
        call_table = SvCallTable(renderer.db_conn)
        self.score_slider.end = 0
        if call_table.num_calls(int(self.run_id_dropdown.value), 0) > 0:
            self.score_slider.end = call_table.max_score(int(self.run_id_dropdown.value))
        renderer.k_mer_counter = KMerFilterTable(renderer.db_conn).get_counter(int(self.run_id_dropdown.value), 1)
        renderer.render()

    def ground_id_change(self, renderer):
        run_table = SvCallerRunTable(renderer.db_conn)
        self.ground_truth_id_dropdown.label = "Selected ground truth: " + \
                                            run_table.getName(int(self.ground_truth_id_dropdown.value)) + \
                                            " - " + self.ground_truth_id_dropdown.value
        renderer.render()

    def slider_change(self, renderer):
        renderer.render()

    def full_render(self, renderer):
        renderer.render(render_all=True)

    def render_mems_button_event(self, renderer):
        if not renderer.selected_read_id is None:
            read = ReadTable(renderer.db_conn).get_read(renderer.selected_read_id)

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
        print("unimplemented at the moment...")
        #renderer.sv_db.delete_run(renderer.get_run_id())
        renderer.setup()
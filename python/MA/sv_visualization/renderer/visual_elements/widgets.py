from bokeh.plotting import figure
from bokeh.models import Button, Slider
from bokeh.models.widgets import Dropdown, TextInput, RadioGroup, CheckboxGroup

class Widgets:
    def __init__(self):
        self.file_input = TextInput(value="/MAdata/sv_datasets/minimal", title="Dataset Name:")
        self.run_id_dropdown = Dropdown(label="select run id here", menu=[])
        self.ground_truth_id_dropdown = Dropdown(label="select ground truth id here", menu=[])
        self.score_slider = Slider(start=0, end=1, value=0, step=.1, callback_policy='mouseup', title="min score")
        self.max_elements_slider = Slider(start=1000, end=100000, value=25000, step=1000,
                                          callback_policy='mouseup', title="max render")
        self.range_link_radio = RadioGroup(labels=["Link read plot to x-range", "Link read plot to y-range"],
                                           active=0)
        self.render_mems_button = Button(label="render MEMs")
        self.reset_button = Button(label="Reset")
        self.delete_button = Button(label="Delete Dataset")
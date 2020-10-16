from renderer import Renderer
from renderer.util import *
import os.path
import json
from bokeh.layouts import column, row, grid
from bokeh.models import Button, Slider, FuncTickFormatter
from bokeh.plotting import figure, curdoc
from bokeh.models.callbacks import CustomJS
from bokeh.events import ButtonClick
from bokeh.models.widgets import Dropdown, TextInput, RadioGroup, CheckboxGroup
from bokeh.models.tools import HoverTool
from bokeh.models.axes import LinearAxis
from MA import *
import math

renderer = Renderer()

with renderer.measure("render"):
    redered_everything = renderer.setup()

# render this document
curdoc().add_root(row(
    grid([[renderer.seed_plot.left_plot, renderer.nuc_plot.left_plot, renderer.main_plot.plot],
          [None,                         None,                        renderer.nuc_plot.bottom_plot],
          [None,                         None,                        renderer.seed_plot.bottom_plot]]),
    column(row(renderer.widgets.file_input,
               renderer.widgets.run_id_dropdown,
               renderer.widgets.ground_truth_id_dropdown),
           renderer.widgets.score_slider,
           renderer.widgets.max_elements_slider,
           renderer.widgets.range_link_radio,
           row(renderer.widgets.full_render_button,
               renderer.widgets.render_mems_button,
               renderer.widgets.delete_button),
           renderer.widgets.force_read_id,
           grid([[renderer.read_plot.nuc_plot.left_plot, renderer.read_plot.plot], 
                 [None,                                  renderer.read_plot.nuc_plot.bottom_plot]]),
            renderer.widgets.subset_buttons,
            renderer.widgets.blur_slider,
            renderer.read_plot.stat_plot ),
    renderer.widgets.spinner_div
))
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
    redered_everything = renderer.render()
renderer.analyze.analyze()

# render this document
curdoc().add_root(row(
    grid([[render.seed_plot.left_plot, render.nuc_plot.left_plot, render.main_plot.plot],
          [None,                       None,                      render.nuc_plot.bottom_plot],
          [None,                       None,                      render.seed_plot.bottom_plot]]),
    column(renderer.widgets.file_input,
           renderer.widgets.run_id_dropdown,
           renderer.widgets.ground_truth_id_dropdown,
           renderer.widgets.score_slider,
           renderer.widgets.max_elements_slider,
           renderer.widgets.range_link_radio,
           renderer.widgets.render_mems_button,
           renderer.widgets.reset_button,
           renderer.widgets.delete_button,
           grid([[render.read_plot.nuc_plot.left_plot, render.read_plot.plot], 
                 [None,                                render.read_plot.nuc_plot.bottom_plot]]))
))

if os.path.isfile(dataset_name + "/info.json"):
    with open(dataset_name + "/info.json", "r") as json_file:
        json_info_file = json.loads(json_file.read(), object_hook=decode)
    ref_genome = json_info_file["reference_path"] + "/ma/genome"

    if server_context.ref_genome == ref_genome:
        print("using cached pack")
    else:
        pack = Pack()
        pack.load(ref_genome)
        fm_index = FMIndex()
        fm_index.load(ref_genome)
        server_context.pack = pack
        server_context.ref_genome = ref_genome
        server_context.fm_index = fm_index
    if server_context.dataset_name == dataset_name:
        print("using cached sv_db")
    else:
        sv_db = SV_DB(dataset_name + "/svs.db", "open")
        server_context.sv_db = sv_db
        server_context.dataset_name = dataset_name

    xe = args_get("xe", float, server_context.pack.unpacked_size_single_strand)
    ye = args_get("ye", float, server_context.pack.unpacked_size_single_strand)

    plot.x_range.start = xs
    plot.x_range.end = xe
    plot.y_range.start = ys
    plot.y_range.end = ye

    renderer = Renderer(plot, [l_plot, l2_plot], [d_plot, d2_plot], xs, xe, ys, ye,
                        server_context.pack, server_context.fm_index, server_context.sv_db, run_id,
                        ground_truth_id, min_score, max_elements, dataset_name, active_tools,
                        radio_group, read_plot, selected_read_id, l_read_plot, d_read_plot, render_mems, seed_plot_y_s, 
                        seed_plot_y_e, json_info_file["reference_path"] + "/ma/genome")

    with renderer.measure("render"):
        redered_everything = renderer.render()
    renderer.analyze.analyze()
else:
    xe = 0
    ye = 0

make_js_callback_file = js_file("reload_callback")


def make_js_callback(condition):
    return CustomJS(args=dict(xr=plot.x_range, yr=plot.y_range, xs=xs, xe=xe, ys=ys, ye=ye,
                              run_id=run_id, min_score=min_score, plot=plot, max_elements=max_elements,
                              ground_truth_id=ground_truth_id, radio_group=radio_group,
                              read_plot_range=l2_plot.x_range, dataset_name=dataset_name,
                              read_plot_y_range=read_plot.y_range, render_mems=0,
                              selected_read_id=selected_read_id),
                    code=condition + make_js_callback_file)

render_mems_button.js_on_event(ButtonClick, make_js_callback("""
        render_mems = 1;
    """))


# callback if area changed too much
plot.y_range.js_on_change('start', make_js_callback("""
        var w = xe - xs;
        var h = ye - ys;
        if(
                xr.start < xs - w * """ + str(render_area_factor) + """ ||
                yr.start < ys - h * """ + str(render_area_factor) + """ ||
                xr.end > xe + w * """ + str(render_area_factor) + """ ||
                yr.end > ye + h * """ + str(render_area_factor) +
    (" " if redered_everything else " || w/4 > xr.end - xr.start || h/4 > yr.end - yr.start ") + """
            )
    """))

menu = []
label = "select run id here"
label_2 = "select ground truth id here"
if not server_context.sv_db is None:
    for idx in range(100):  # sv_db.get_num_runs() + 1
        if server_context.sv_db.run_exists(idx):
            text = server_context.sv_db.get_run_name(idx) + " - " + \
                server_context.sv_db.get_run_date(idx) + " - " + \
                server_context.sv_db.get_run_desc(idx)
            text += " - " + \
                str(server_context.sv_db.get_num_calls(idx, min_score))
            if idx == run_id:
                label = text
            if idx == ground_truth_id:
                label_2 = text
            menu.append((text, str(idx)))

run_id_dropdown = Dropdown(label=label, menu=menu)
run_id_dropdown.js_on_change("value", make_js_callback("""
        run_id = cb_obj.value;
    """))
ground_truth_id_dropdown = Dropdown(label=label_2, menu=menu)
ground_truth_id_dropdown.js_on_change("value", make_js_callback("""
        ground_truth_id = cb_obj.value;
    """))

max_score_slider = 100
try:
    if not server_context.sv_db is None and server_context.sv_db.run_exists(run_id):
            max_score_slider = server_context.sv_db.get_max_score(run_id)
except:
    pass

score_slider.callback = make_js_callback("""
        min_score = cb_obj.value;
    """)

max_elements_slider.callback = make_js_callback("""
        max_elements = cb_obj.value;
    """)

file_input.js_on_change("value", make_js_callback("""
        dataset_name = cb_obj.value;
    """))

# render this document
curdoc().add_root(row(
    grid([[l2_plot, l_plot, plot], [None, None, d_plot], [None, None, d2_plot]]),
    column(file_input, run_id_dropdown, ground_truth_id_dropdown, score_slider,
           max_elements_slider, radio_group, render_mems_button, reset_button, delete_button,
           grid([[l_read_plot, read_plot], [None, d_read_plot]]))
))

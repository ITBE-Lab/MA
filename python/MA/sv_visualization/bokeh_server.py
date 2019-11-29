# from render_genome_overview import render_overview
from render_region import render_region
import os.path
import json
from bokeh.layouts import column, row
from bokeh.models import Button, Slider
from bokeh.plotting import figure, curdoc
from bokeh.models.callbacks import CustomJS
from bokeh.events import ButtonClick
from bokeh.models.widgets import Dropdown, TextInput
from bokeh.models.tools import HoverTool
from MA import *

# AKFIX
"""Markus @ Zeus""" 
# svdb_dir = "/MAdata/sv_datasets/" # AKFIX

"""Arne @ home """
svdb_dir = "C:/Users/Markus/Desktop/MA-Database/sv_datasets/" 


def _decode(o):
    if isinstance(o, str):
        try:
            return float(o)
        except ValueError:
            return o
    elif isinstance(o, dict):
        return {_decode(k): _decode(v) for k, v in o.items()}
    elif isinstance(o, list):
        return [_decode(v) for v in o]
    else:
        return o

args = curdoc().session_context.request.arguments

def args_get(name, convert, default):
    if not args.get(name) is None:
        if convert is None:
            return args.get(name)[0]
        else:
            return convert(args.get(name)[0])
    else:
        return default

dataset_name = args_get("dataset_name", None, b'minimal-2').decode()
json_info_file = None # noop

render_area_factor = 1

xs = args_get("xs", float, 0)
ys = args_get("ys", float, 0)
run_id = args_get("run_id", int, -1)
ground_truth_id = args_get("ground_truth_id", int, -1)
min_score = args_get("min_score", float, 0)
max_elements = args_get("max_elements", int, 2000)
active_tools = args_get("active_tools", None, b"pan.wheel_zoom.tap.hover1,hover2,hover3").decode().split(".")
print(active_tools)

def text_to_none(text):
    if text == "None":
        return None
    return text

sv_db = None

plot = figure(
            width=900,
            height=900,
            tools=[
                "pan", "box_zoom", 
                "wheel_zoom", "save",
                "tap"
            ],
            active_drag=text_to_none(active_tools[0]),
            active_scroll=text_to_none(active_tools[1]),
            active_tap=text_to_none(active_tools[2])
        )
hover1 = HoverTool(tooltips=[("from", "@f"), ("to", "@t"), ("#calls", "@i")], names=['hover1'],
                   name="Hover heatmap")
hover2 = HoverTool(tooltips=[("supp. nt", "@n"), ("coverage", "@c"), ("#reads", "@r"), ("score", "@s")],
                   names=['hover2'], name="Hover calls")
hover3 = HoverTool(tooltips=[("supp. nt", "@n"), ("read id", "@r"), ("|query|:", "@q")],
                   names=['hover3'], name="Hover jumps")
plot.add_tools(hover1)
plot.add_tools(hover2)
plot.add_tools(hover3)

plot.tools[0].name = "pan"
plot.tools[1].name = "box_zoom"
plot.tools[2].name = "wheel_zoom"
plot.tools[3].name = "tap"
plot.tools[4].name = "save"
plot.tools[5].name = "hover1"
plot.tools[6].name = "hover2"
plot.tools[7].name = "hover3"

l = []
for x in active_tools[3].split(","):
    if(x == "hover1"):
        l.append(hover1)
    if(x == "hover2"):
        l.append(hover2)
    if(x == "hover3"):
        l.append(hover3)
plot.toolbar.active_inspect = l

redered_everything = True
print("File: ", svdb_dir + dataset_name + "/info.json")
if os.path.isfile(svdb_dir + dataset_name + "/info.json"):
    with open(svdb_dir + dataset_name + "/info.json", "r") as json_file:
        json_info_file = json.loads(json_file.read(), object_hook=_decode)
    ref_genome = json_info_file["reference_path"] + "/ma/genome"

    pack = Pack()
    pack.load(ref_genome)

    xe = args_get("xe", float, pack.unpacked_size_single_strand)
    ye = args_get("ye", float, pack.unpacked_size_single_strand)

    plot.x_range.start = xs
    plot.x_range.end = xe
    plot.y_range.start = ys
    plot.y_range.end = ye

    sv_db = SV_DB(svdb_dir + dataset_name + "/svs.db", "open")

    redered_everything = render_region(plot, xs, xe, ys, ye, pack, sv_db, run_id, ground_truth_id, min_score,
                                             max_elements, dataset_name, active_tools)

# reset button
reset_button = Button(label="Reset")
reset_button.js_on_event(ButtonClick, CustomJS(code="""
    document.location.href = document.location.href.split("?")[0];
"""))

def make_js_callback(condition):
    return CustomJS(args=dict(xr=plot.x_range, yr=plot.y_range, xs=xs, xe=xe, ys=ys, ye=ye,
                              run_id=run_id, min_score=min_score, plot=plot, max_elements=max_elements,
                              ground_truth_id=ground_truth_id, dataset_name=dataset_name),
                    code=condition + """
            {
                if (typeof window.left_page_already !== 'undefined')
                    return;
                else
                    window.left_page_already = 1;
                var active_drag = "None";
                if(plot.toolbar.tools[0].active)
                    active_drag = "pan";
                if(plot.toolbar.tools[1].active)
                    active_drag = "box_zoom";
                var active_scroll = "None";
                if(plot.toolbar.tools[2].active)
                    active_scroll = "wheel_zoom";
                var active_tap = "tap"; // hmm the active does not work...?
                if(plot.toolbar.tools[3].active)
                    active_tap = "tap";
                var active_inspect = "";
                for(var x = 5; x < 8; x++)
                    if(plot.toolbar.tools[x].active)
                        active_inspect += "," + plot.toolbar.tools[x].name;

                var active_tools = active_drag + "." + active_scroll + "." + active_tap + "." + active_inspect;
                
                plot.toolbar.active_drag = null;
                plot.toolbar.active_scroll = null;
                plot.toolbar.active_tap = null;
                plot.toolbar.tools = [];
                s = document.location.href.split("?")[0] + "?xs=" + xr.start +
                                                        "&xe=" + xr.end +
                                                        "&ys=" + yr.start +
                                                        "&ye=" + yr.end +
                                                        "&run_id=" + run_id +
                                                        "&max_elements=" + max_elements +
                                                        "&dataset_name=" + dataset_name +
                                                        "&min_score=" + min_score +
                                                        "&ground_truth_id=" + ground_truth_id +
                                                        "&active_tools=" + active_tools;
                //alert(s);
                document.location.href = s;
            } // if or scope
        """)

# callback if area changed too much
plot.y_range.js_on_change('start', make_js_callback("""
        var w = xe - xs;
        var h = ye - ys;
        if(
                xr.start < xs - w * """ + str(render_area_factor) + """ ||
                yr.start < ys - h * """ + str(render_area_factor) + """ ||
                xr.end > xe + w * """ + str(render_area_factor) + """ ||
                yr.end > ye + h * """ + str(render_area_factor) +
                ( " " if redered_everything else " || w/4 > xr.end - xr.start || h/4 > yr.end - yr.start " ) + """
            )
    """))

menu = []
label = "select run id here"
label_2 = "select ground truth id here"
if not sv_db is None:
    for idx in range(sv_db.get_num_runs() + 1):
        if sv_db.run_exists(idx):
            text = sv_db.get_run_name(idx) + " - " + sv_db.get_run_date(idx) + " - " + sv_db.get_run_desc(idx)
            text += " - " + str(sv_db.get_num_calls(idx, min_score))
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

if not sv_db is None and sv_db.run_exists(run_id):
    max_score_slider = sv_db.get_max_score(run_id)
else:
    max_score_slider = 100
score_slider = Slider(start=0, end=max_score_slider, 
                      value=min_score, step=.1, callback_policy = 'mouseup', title="min score")
score_slider.callback = make_js_callback("""
        min_score = cb_obj.value;
    """)

max_elements_slider = Slider(start=1000, end=100000, value=max_elements, step=1000, callback_policy = 'mouseup',
                             title="max render")
max_elements_slider.callback = make_js_callback("""
        max_elements = cb_obj.value;
    """)

file_input = TextInput(value=dataset_name, title="Databaset Name:")
file_input.js_on_change("value", make_js_callback("""
        dataset_name = cb_obj.value;
    """))

# render this document
curdoc().add_root(row(plot, column(file_input, run_id_dropdown, ground_truth_id_dropdown, score_slider,
                                   max_elements_slider, reset_button)))

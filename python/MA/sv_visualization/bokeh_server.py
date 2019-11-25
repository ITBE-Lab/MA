from render_region import render_region
import os.path
import json
from bokeh.layouts import column, row, grid
from bokeh.models import Button, Slider, FuncTickFormatter
from bokeh.plotting import figure, curdoc
from bokeh.models.callbacks import CustomJS
from bokeh.events import ButtonClick
from bokeh.models.widgets import Dropdown, TextInput, CheckboxGroup
from bokeh.models.tools import HoverTool
from bokeh.models.tickers import FixedTicker
from bokeh.models.axes import LinearAxis
from MA import *
import math

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

server_context = curdoc().session_context.server_context

args = curdoc().session_context.request.arguments

def args_get(name, convert, default):
    if not args.get(name) is None:
        if convert is None:
            return args.get(name)[0]
        else:
            return convert(args.get(name)[0])
    else:
        return default

dataset_name = args_get("dataset_name", None, b'minimal').decode()
json_info_file = None # noop

render_area_factor = 1

xs = args_get("xs", float, 0)
ys = args_get("ys", float, 0)
render_read_indices = args_get("render_read_indices", None, b'false').decode()
x_range_link = args_get("x_range_link", None, b'false').decode()
y_range_link = args_get("y_range_link", None, b'false').decode()
run_id = args_get("run_id", int, -1)
ground_truth_id = args_get("ground_truth_id", int, -1)
min_score = args_get("min_score", float, 0)
max_elements = args_get("max_elements", int, 2000)
active_tools = args_get("active_tools", None, b"pan.wheel_zoom.hover1,hover2,hover3,hover4").decode().split(".")

def text_to_none(text):
    if text == "None":
        return None
    return text

plot = figure(
            width=900,
            height=900,
            tools=[
                "pan", "box_zoom", 
                "wheel_zoom", "save"
            ],
            active_drag=text_to_none(active_tools[0]),
            active_scroll=text_to_none(active_tools[1])
        )
plot.axis.visible = False
l_plot = figure(
            width=40,
            height=900,
            y_range=plot.y_range,
            tools=["ypan"],
            toolbar_location=None
        )
l_plot.axis.visible = False
l_plot.grid.visible = False

d_plot = figure(
            width=900,
            height=40,
            x_range=plot.x_range,
            tools=["xpan"],
            toolbar_location=None
        )
d_plot.axis.visible = False
d_plot.grid.visible = False

l2_plot = figure(
            width=300,
            height=900,
            y_range=plot.y_range,
            x_range=list(),
            tools=["xpan", "xwheel_zoom"],
            active_scroll="xwheel_zoom",
            toolbar_location=None
        )
l2_plot.xaxis.major_label_orientation = math.pi/2
l2_plot.yaxis.axis_label = "Reference Position"
l2_plot.xaxis.axis_label = "Read Id"


# Remove default axis
#l2_plot.xaxis.visible = False
# Add custom axis
# ticker = FixedTicker()
# ticker.ticks = [1,2,3]
# xaxis = LinearAxis(ticker=ticker)
# l2_plot.add_layout(xaxis, 'below')
# xaxis.formatter = FuncTickFormatter(code="return Math.floor(tick);")

d2_plot = figure(
            width=900,
            height=300,
            x_range=plot.x_range,
            y_range=l2_plot.x_range,
            tools=["ypan", "ywheel_zoom"],
            active_scroll="ywheel_zoom",
            toolbar_location=None
        )
d2_plot.xaxis.axis_label = "Reference Position"
d2_plot.yaxis.axis_label = "Read Id"


hover1 = HoverTool(tooltips=[("from", "@f"), ("to", "@t"), ("#calls", "@i")], names=['hover1'],
                   name="Hover heatmap")
hover2 = HoverTool(tooltips=[("supp. nt", "@n"), ("coverage", "@c"), ("#reads", "@r"), ("score", "@s")],
                   names=['hover2'], name="Hover calls")
hover3 = HoverTool(tooltips=[("supp. nt", "@n"), ("read id", "@r"), ("|query|", "@q"), ("from", "@f"), ("to", "@t")],
                   names=['hover3'], name="Hover jumps")
hover4 = HoverTool(tooltips="@i", names=['hover4'], name="Hover simple")
plot.add_tools(hover1)
plot.add_tools(hover2)
plot.add_tools(hover3)
plot.add_tools(hover4)

plot.tools[0].name = "pan"
plot.tools[1].name = "box_zoom"
plot.tools[2].name = "wheel_zoom"
plot.tools[3].name = "save"
plot.tools[4].name = "hover1"
plot.tools[5].name = "hover2"
plot.tools[6].name = "hover3"
plot.tools[7].name = "hover4"


read_plot = figure(
            width=900,
            height=400,
            tools=[
                "pan", "box_zoom", 
                "wheel_zoom", "save", "reset"
            ]
        )

l = []
for x in active_tools[2].split(","):
    if(x == "hover1"):
        l.append(hover1)
    if(x == "hover2"):
        l.append(hover2)
    if(x == "hover3"):
        l.append(hover3)
    if(x == "hover4"):
        l.append(hover4)
plot.toolbar.active_inspect = l

hover5 = HoverTool(tooltips=[("read id", "@r_id"), ("q, r, l", "@q, @r, @size"), ("index", "@idx"),
                             ("reseeding-layer", "@layer")],
                   names=['hover5'], name="Hover reads")
l2_plot.add_tools(hover5)
d2_plot.add_tools(hover5)
read_plot.add_tools(hover5)

checkbox_group = CheckboxGroup(labels=["Render read indices", "Link read plot to x-range", "Link read plot to y-range"],
                               active=[])
if render_read_indices == "true":
    checkbox_group.active.append(0) # Read Indices
if x_range_link == "true":
    # @todo this requires a reload from the server
    checkbox_group.active.append(1) # x-range link
    read_plot.x_range = plot.x_range
elif y_range_link == "true":
    # @todo this requires a reload from the server
    # @todo y-range of the read plot should auto update if x-range is changed from other plot...
    checkbox_group.active.append(2) # y-range-link
    read_plot.x_range = plot.y_range

redered_everything = True

if os.path.isfile("/MAdata/sv_datasets/" + dataset_name + "/info.json"):
    with open("/MAdata/sv_datasets/" + dataset_name + "/info.json", "r") as json_file:
        json_info_file = json.loads(json_file.read(), object_hook=_decode)
    ref_genome = json_info_file["reference_path"] + "/ma/genome"

    if not hasattr(server_context, "ref_genome"):
        setattr(server_context, "ref_genome", None)
        setattr(server_context, "pack", None)
        setattr(server_context, "sv_db", None)
        setattr(server_context, "dataset_name", None)
    if server_context.ref_genome == ref_genome:
        print("using cached pack")
    else:
        pack = Pack()
        pack.load(ref_genome)
        server_context.pack = pack
        server_context.ref_genome = ref_genome
    if server_context.dataset_name == dataset_name:
        print("using cached sv_db")
    else:
        sv_db = SV_DB("/MAdata/sv_datasets/" + dataset_name + "/svs.db", "open")
        server_context.sv_db = sv_db
        server_context.dataset_name = dataset_name

    xe = args_get("xe", float, server_context.pack.unpacked_size_single_strand)
    ye = args_get("ye", float, server_context.pack.unpacked_size_single_strand)

    plot.x_range.start = xs
    plot.x_range.end = xe
    plot.y_range.start = ys
    plot.y_range.end = ye

    redered_everything = render_region(plot, [l_plot, l2_plot], [d_plot, d2_plot], xs, xe, ys, ye,
                                        server_context.pack, server_context.sv_db, run_id,
                                        ground_truth_id, min_score, max_elements, dataset_name, active_tools, 
                                        checkbox_group, read_plot, json_info_file["reference_path"] + "/ma/genome")
else:
    xe = 0
    ye = 0

# reset button
reset_button = Button(label="Reset")
reset_button.js_on_event(ButtonClick, CustomJS(code="""
    document.location.href = document.location.href.split("?")[0];
"""))


def make_js_callback(condition):
    return CustomJS(args=dict(xr=plot.x_range, yr=plot.y_range, xs=xs, xe=xe, ys=ys, ye=ye,
                              run_id=run_id, min_score=min_score, plot=plot, max_elements=max_elements,
                              ground_truth_id=ground_truth_id, checkbox_group=checkbox_group, 
                              dataset_name=dataset_name),
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
                var active_inspect = "";
                for(var x = 4; x < 8; x++)
                    if(plot.toolbar.tools[x].active)
                        active_inspect += "," + plot.toolbar.tools[x].name;

                var active_tools = active_drag + "." + active_scroll + "." + active_inspect;
                
                plot.toolbar.active_drag = null;
                plot.toolbar.active_scroll = null;
                plot.toolbar.active_tap = null;
                plot.toolbar.tools = [];
                debugger;
                s = document.location.href.split("?")[0] + "?xs=" + xr.start +
                                                        "&xe=" + xr.end +
                                                        "&ys=" + yr.start +
                                                        "&ye=" + yr.end +
                                                        "&run_id=" + run_id +
                                                        "&max_elements=" + max_elements +
                                                        "&dataset_name=" + dataset_name +
                                                        "&min_score=" + min_score +
                                                        "&ground_truth_id=" + ground_truth_id +
                                                        "&render_read_indices=" + checkbox_group.active.includes(0) +
                                                        "&x_range_link=" + checkbox_group.active.includes(1) +
                                                        "&y_range_link=" + checkbox_group.active.includes(2) +
                                                        "&active_tools=" + active_tools;
                //alert(s);
                document.location.href = s;
            } // if or scope
        """)

# reset button
delete_button = Button(label="Delete Dataset")
delete_button.js_on_event(ButtonClick, make_js_callback("""
        s = document.location.href.split("bokeh_server")
        //alert(s);
        document.location.href = s[0] + "delete_run" + s[1];
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
                ( " " if redered_everything else " || w/4 > xr.end - xr.start || h/4 > yr.end - yr.start " ) + """
            )
    """))

menu = []
label = "select run id here"
label_2 = "select ground truth id here"
if not server_context.sv_db is None:
    for idx in range(100): #sv_db.get_num_runs() + 1
        if server_context.sv_db.run_exists(idx):
            text = server_context.sv_db.get_run_name(idx) + " - " + \
                    server_context.sv_db.get_run_date(idx) + " - " + \
                    server_context.sv_db.get_run_desc(idx)
            text += " - " + str(server_context.sv_db.get_num_calls(idx, min_score))
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

if not server_context.sv_db is None and server_context.sv_db.run_exists(run_id):
    max_score_slider = server_context.sv_db.get_max_score(run_id)
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
curdoc().add_root(row(
                        grid([[l2_plot, l_plot, plot], [None, None, d_plot], [None, None, d2_plot]]),
                        column(file_input, run_id_dropdown, ground_truth_id_dropdown, score_slider,
                                   max_elements_slider, checkbox_group, reset_button, delete_button, read_plot)
                        ))

from render_genome_overview import render_overview
from render_region import render_region
import json
from bokeh.layouts import column, row
from bokeh.models import Button, Slider
from bokeh.plotting import figure, curdoc
from bokeh.models.callbacks import CustomJS
from bokeh.events import ButtonClick
from bokeh.models.widgets import Dropdown
from MA import *

print(" --- GOT REQUEST")

args = curdoc().session_context.request.arguments

def args_get(name, convert, default):
    if not args.get(name) is None:
        return convert(args.get(name)[0])
    else:
        return default

dataset_name = "minimal-z"
json_info_file = None # noop
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
with open("/MAdata/sv_datasets/" + dataset_name + "/info.json", "r") as json_file:
    json_info_file = json.loads(json_file.read(), object_hook=_decode)
db_path = "/MAdata/sv_datasets/" + dataset_name + "/svs.db"
ref_genome = json_info_file["reference_path"] + "/ma/genome"


pack = Pack()
pack.load(ref_genome)
sv_db = SV_DB(db_path, "open")

xs = args_get("xs", float, 0)
ys = args_get("ys", float, 0)
xe = args_get("xe", float, pack.unpacked_size_single_strand)
ye = args_get("ye", float, pack.unpacked_size_single_strand)
run_id = args_get("run_id", int, -1)
ground_truth_id = args_get("ground_truth_id", int, -1)
min_score = args_get("min_score", float, 0)
max_elements = args_get("max_elements", int, 2000)

plot, redered_everything = render_region(xs, xe, ys, ye, pack, sv_db, run_id, ground_truth_id, min_score, max_elements)

# reset button
reset_button = Button(label="Reset")
reset_button.js_on_event(ButtonClick, CustomJS(code="""
    document.location.href = document.location.href.split("?")[0];
"""))

# callback if area changed too much
callback = CustomJS(args=dict(xr=plot.x_range, yr=plot.y_range, xs=xs, xe=xe, ys=ys, ye=ye,
                              run_id=run_id, min_score=min_score, plot=plot, max_elements=max_elements,
                              ground_truth_id=ground_truth_id),
    code="""
        function dist(a, b){
            return Math.abs( a - b );
        }
        var w = xe - xs;
        var h = ye - ys;
        if(dist(xs, xr.start) > w || dist(ys, yr.start) > h || """ +
           ( " " if redered_everything else "w/4 > xr.end - xr.start || h/4 > yr.end - yr.start || " ) + """
           w*4 < xr.end - xr.start || h*4 < yr.end - yr.start)
        {
            if (typeof window.left_page_already !== 'undefined')
                return;
            else
                window.left_page_already = 1;
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
                                                    "&min_score=" + min_score +
                                                    "&ground_truth_id=" + ground_truth_id;
            //alert(s);
            document.location.href = s;
        }
    """)
#plot.x_range.js_on_change('start', callback)
plot.y_range.js_on_change('start', callback)


menu = []
label = "select run id here"
label_2 = "select run id here"
for idx in range(sv_db.get_num_runs()):
    if sv_db.run_exists(idx):
        text = sv_db.get_run_name(idx) + " - " + sv_db.get_run_date(idx) + " - " + sv_db.get_run_desc(idx)
        text += " - " + str(sv_db.get_num_calls(idx, min_score))
        if idx == run_id:
            label = text
        if idx == ground_truth_id:
            label_2 = text
        menu.append((text, str(idx)))
run_id_dropdown = Dropdown(label=label, menu=menu)
run_id_dropdown.js_on_change("value", CustomJS(args=dict(plot=plot, xr=plot.x_range, yr=plot.y_range,
                                                         min_score=min_score, max_elements=max_elements,
                                                         ground_truth_id=ground_truth_id),
    code="""
        if (typeof window.left_page_already !== 'undefined')
            return;
        else
            window.left_page_already = 1;
        plot.toolbar.active_drag = null;
        plot.toolbar.active_scroll = null;
        plot.toolbar.active_tap = null;
        plot.toolbar.tools = [];
        s = document.location.href.split("?")[0] + "?xs=" + xr.start +
                                                "&xe=" + xr.end +
                                                "&ys=" + yr.start +
                                                "&ye=" + yr.end +
                                                "&run_id=" + cb_obj.value +
                                                "&max_elements=" + max_elements +
                                                "&min_score=" + min_score +
                                                "&ground_truth_id=" + ground_truth_id;
        //alert(s);
        document.location.href = s;
    """))
ground_truth_id_dropdown = Dropdown(label=label_2, menu=menu)
ground_truth_id_dropdown.js_on_change("value", CustomJS(args=dict(plot=plot, run_id=run_id, xr=plot.x_range,
                                                                  yr=plot.y_range, min_score=min_score,
                                                                  max_elements=max_elements),
    code="""
        if (typeof window.left_page_already !== 'undefined')
            return;
        else
            window.left_page_already = 1;
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
                                                "&min_score=" + min_score +
                                                "&ground_truth_id=" + cb_obj.value;
        //alert(s);
        document.location.href = s;
    """))

score_slider = Slider(start=0, end=sv_db.get_max_score(run_id) if sv_db.run_exists(run_id) else 100, 
                      value=min_score, step=.1, callback_policy = 'mouseup', title="min score")
score_slider.callback = CustomJS(args=dict(plot=plot, xr=plot.x_range, yr=plot.y_range, 
                                           run_id=run_id, max_elements=max_elements, ground_truth_id=ground_truth_id),
    code="""
        if (typeof window.left_page_already !== 'undefined')
            return;
        else
            window.left_page_already = 1;
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
                                                "&min_score=" + cb_obj.value +
                                                "&ground_truth_id=" + ground_truth_id;
        //alert(s);
        document.location.href = s;
    """)

max_elements_slider = Slider(start=1000, end=100000, value=max_elements, step=1000, callback_policy = 'mouseup',
                             title="max render")
max_elements_slider.callback = CustomJS(args=dict(plot=plot, xr=plot.x_range, yr=plot.y_range, run_id=run_id,
                                                  min_score=min_score, ground_truth_id=ground_truth_id),
    code="""
        if (typeof window.left_page_already !== 'undefined')
            return;
        else
            window.left_page_already = 1;
        plot.toolbar.active_drag = null;
        plot.toolbar.active_scroll = null;
        plot.toolbar.active_tap = null;
        plot.toolbar.tools = [];
        s = document.location.href.split("?")[0] + "?xs=" + xr.start +
                                                "&xe=" + xr.end +
                                                "&ys=" + yr.start +
                                                "&ye=" + yr.end +
                                                "&run_id=" + run_id +
                                                "&min_score=" + min_score +
                                                "&max_elements=" + cb_obj.value +
                                                "&ground_truth_id=" + ground_truth_id;
        //alert(s);
        document.location.href = s;
    """)

print(" --- ANSWERED REQUEST")

# render this document
curdoc().add_root(row(plot, column(run_id_dropdown, ground_truth_id_dropdown, score_slider, max_elements_slider,
                                   reset_button)))

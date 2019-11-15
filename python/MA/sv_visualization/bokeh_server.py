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
min_score = args_get("min_score", float, 0)

plot = None
if xe - xs > 100000000 or ye - ys > 100000000:
    plot = render_overview(xs, xe, ys, ye, pack, sv_db, run_id, min_score)
else:
    plot = render_region(xs, xe, ys, ye, pack, sv_db, run_id, min_score)

# reset button
reset_button = Button(label="Reset")
reset_button.js_on_event(ButtonClick, CustomJS(code="""
    document.location.href = document.location.href.split("?")[0];
"""))

# callback if area changed too much
callback = CustomJS(args=dict(xr=plot.x_range, yr=plot.y_range, xs=xs, xe=xe, ys=ys, ye=ye,
                              run_id=run_id, min_score=min_score),
    code="""
        function dist(a, b){
            return Math.abs( a - b );
        }
        var w = xe - xs;
        var h = ye - ys;
        var max_d = 0.3;
        if(dist(xs, xr.start) > max_d*w || dist(ys, yr.start) > max_d*h ||
           dist(w, xr.end - xr.start) > max_d*w || dist(h, yr.end - yr.start) > max_d*h)
        {
            if (typeof window.left_page_already !== 'undefined')
                return;
            else
                window.left_page_already = 1;
            s = document.location.href.split("?")[0] + "?xs=" + xr.start +
                                                    "&xe=" + xr.end +
                                                    "&ys=" + yr.start +
                                                    "&ye=" + yr.end +
                                                    "&run_id=" + run_id +
                                                    "&min_score=" + min_score;
            //alert(s);
            document.location.href = s;
        }
    """)
#plot.x_range.js_on_change('start', callback)
plot.y_range.js_on_change('start', callback)


menu = []
label = "select run id here"
for idx in range(sv_db.get_num_runs()):
    if sv_db.run_exists(idx):
        text = sv_db.get_run_name(idx) + " - " + sv_db.get_run_date(idx) + " - " + sv_db.get_run_desc(idx)
        text += " - " + str(sv_db.get_num_calls(idx, min_score))
        if idx == run_id:
            label = text
        menu.append((text, str(idx)))
run_id_dropdown = Dropdown(label=label, menu=menu)
run_id_dropdown.js_on_change("value", CustomJS(args=dict(xr=plot.x_range, yr=plot.y_range, ),
    code="""
        if (typeof window.left_page_already !== 'undefined')
            return;
        else
            window.left_page_already = 1;
        s = document.location.href.split("?")[0] + "?xs=" + xr.start +
                                                "&xe=" + xr.end +
                                                "&ys=" + yr.start +
                                                "&ye=" + yr.end +
                                                "&run_id=" + cb_obj.value
        //alert(s);
        document.location.href = s;
    """))

score_slider = Slider(start=0, end=sv_db.get_max_score(run_id) if sv_db.run_exists(run_id) else 100, 
                      value=min_score, step=.1, callback_policy = 'mouseup', title="min score")
score_slider.callback = CustomJS(args=dict(xr=plot.x_range, yr=plot.y_range, run_id=run_id),
    code="""
        if (typeof window.left_page_already !== 'undefined')
            return;
        else
            window.left_page_already = 1;
        s = document.location.href.split("?")[0] + "?xs=" + xr.start +
                                                "&xe=" + xr.end +
                                                "&ys=" + yr.start +
                                                "&ye=" + yr.end +
                                                "&run_id=" + run_id +
                                                "&min_score=" + cb_obj.value;
        //alert(s);
        document.location.href = s;
    """)

print(" --- ANSWERED REQUEST")

# render this document
curdoc().add_root(row(plot, column(run_id_dropdown, score_slider, reset_button)))

from render_genome_overview import render_overview

import json
from bokeh.layouts import layout
from bokeh.models import Button
from bokeh.plotting import figure, curdoc
from bokeh.models.callbacks import CustomJS
from bokeh.events import ButtonClick

print(" --- GOT REQUEST")

args = curdoc().session_context.request.arguments

def args_get(name, convert, default):
    if not args.get(name) is None:
        return convert(args.get(name)[0])
    else:
        return default


xs = args_get("xs", float, 0)
ys = args_get("ys", float, 0)
xe = args_get("xe", float, 3*10**9)
ye = args_get("ye", float, 3*10**9)

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

print("plot dimensions:", xs, xe, ys, ye)
plot = None
if xe-xs > 1000 or ye - ys > 1000:
    plot = render_overview(xs, xe, ys, ye, db_path, ref_genome)
else:
    plot = figure(x_range=(xs, xe), y_range=(ys, ye), toolbar_location=None)

# reset button
reset_button = Button(label="Reset")
reset_button.js_on_event(ButtonClick, CustomJS(code="""
    document.location.href = document.location.href.split("?")[0];
"""))

# callback if area changed too much
callback = CustomJS(args=dict(xr=plot.x_range, yr=plot.y_range, xs=xs, xe=xe, ys=ys, ye=ye),
    code="""
        function dist(a, b){
            return Math.abs( a - b );
        }
        return;
        var max_d = 100;
        if(dist(xs, xr.start) > max_d || dist(xe, xr.end) > max_d ||
           dist(ys, yr.start) > max_d || dist(ye, yr.end) > max_d)
        {
            s = document.location.href.split("?")[0] + "?xs=" + xr.start +
                                                    "&xe=" + xr.end +
                                                    "&ys=" + yr.start +
                                                    "&ye=" + yr.end;
            //alert(s);
            document.location.href = s;
        }
    """)
plot.x_range.js_on_change('start', callback)
plot.y_range.js_on_change('start', callback)


print(" --- ANSWERED REQUEST")

# render this document
curdoc().add_root(layout([[plot, reset_button]]))

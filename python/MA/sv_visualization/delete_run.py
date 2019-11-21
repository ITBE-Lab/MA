import os.path
from bokeh.layouts import column
from bokeh.models import Button
from bokeh.plotting import curdoc
from bokeh.models.callbacks import CustomJS
from bokeh.events import ButtonClick
from MA import *


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

run_id = args_get("run_id", int, -1)

def text_to_none(text):
    if text == "None":
        return None
    return text

sv_db = None
if os.path.isfile("/MAdata/sv_datasets/" + dataset_name + "/info.json"):
    sv_db = SV_DB("/MAdata/sv_datasets/" + dataset_name + "/svs.db", "open")
    if sv_db.run_exists(run_id):
        sv_db.delete_run(run_id)
    else:
        print("tried deleting run that does not exist...")


# reset button
done_button = Button(label="Back")
done_button.js_on_event(ButtonClick, CustomJS(code="""
        s = document.location.href.split("delete_run")
        //alert(s);
        document.location.href = s[0] + "bokeh_server" + s[1];
    """))

# render this document
curdoc().add_root(column(done_button))

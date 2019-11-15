from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import FuncTickFormatter, TapTool, OpenURL
from MA import *

def render_region(xs, xe, ys, ye, pack, sv_db, run_id, min_score):
    plot = figure(
            x_range=(xs, xe),
            y_range=(ys, ye),
            width=900,
            height=900,
            tooltips="@i",
            tools=[
                "pan", #"wheel_zoom",
                "box_zoom", "save",
                "reset", "hover", "tap"
            ],
            active_drag="box_zoom",
            #active_scroll="wheel_zoom"
            )

    if xs < 0:
        xs = 0
    if ys < 0:
        ys = 0
    if xe < 0:
        xe = 0
    if ye < 0:
        ye = 0

    params = ParameterSetManager()
    calls_from_db = SvCallsFromDb(params, sv_db, run_id, int(xs), int(ys), int(xe-xs), int(ye-ys), min_score)
    accepted_boxes_data = {
        "x": [],
        "w": [],
        "y": [],
        "h": [],
        "i": []
    }
    accepted_plus_data = {
        "x": [],
        "y": [],
        "i": []
    }
    max_render = 1000
    cnt_render = 0
    while calls_from_db.hasNext():
        def score(jump):
            if jump.coverage == 0:
                return ""
            return " score: " + str(jump.num_supp_nt / jump.coverage)
        jump = calls_from_db.next()
        if jump.num_supp_nt > min_score * jump.coverage:
            cnt_render += 1
            if cnt_render >= max_render:
                print("hit max_render you wont see the full picture")
                break
            if jump.from_size == 1 and jump.to_size == 1:
                accepted_plus_data["x"].append(jump.from_start)
                accepted_plus_data["y"].append(jump.to_start)
                accepted_plus_data["i"].append(" suppNt: " + str(jump.num_supp_nt) + " cov: " +
                                            str(jump.coverage) + " #reads: " + str(len(jump.supporing_jump_ids)) + 
                                            score(jump))
            else:
                accepted_boxes_data["x"].append(jump.from_start - 0.5)
                accepted_boxes_data["y"].append(jump.to_start - 0.5)
                accepted_boxes_data["w"].append(jump.from_start + jump.from_size + 1)
                accepted_boxes_data["h"].append(jump.to_start + jump.to_size + 1)
                accepted_boxes_data["i"].append(" suppNt: " + str(jump.num_supp_nt) + " cov: " +
                                            str(jump.coverage) + " #reads: " + str(len(jump.supporing_jump_ids)) + 
                                            score(jump))
                accepted_plus_data["x"].append(jump.from_start)
                accepted_plus_data["y"].append(jump.to_start)
                accepted_plus_data["i"].append(" suppNt: " + str(jump.num_supp_nt) + " cov: " +
                                            str(jump.coverage) + " #reads: " + str(len(jump.supporing_jump_ids)) + 
                                            score(jump))
    plot.quad(left="x", bottom="y", right="w", top="h",
              source=ColumnDataSource(accepted_boxes_data))
    plot.x(x="x", y="y", source=ColumnDataSource(accepted_plus_data))

    return plot

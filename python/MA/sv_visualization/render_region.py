from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import FuncTickFormatter, TapTool, OpenURL
from MA import *
import math

def light_spec_approximation(x):
    #map input [0, 1] to wavelength [350, 645]
    w = 370 + x * (645-370)
    r = 0.0
    g = 0.0
    b = 0.0
    if w < 440:
        r = -(w - 440.) / (440. - 380.)
        b = 1.0
    elif w >= 440 and w < 490:
        g = (w - 440.) / (490. - 440.)
        b = 1.0
    elif w >= 490 and w < 510:
        g = 1.0
        b = -(w - 510.) / (510. - 490.)
    elif w >= 510 and w < 580:
        r = (w - 510.) / (580. - 510.)
        g = 1.0
    elif w >= 580 and w < 645:
        r = 1.0
        g = -(w - 645.) / (645. - 580.)
    elif w >= 645:
        r = 1.0

    #intensity
    i = 1.0
    if w > 650:
        i = .3 + .7*(780-w)/(780-650)
    elif w < 420:
        i = .3 + .7*(w-380)/(420-380)

    #gamma
    m = .8

    return (i*r**m, i*g**m, i*b**m)

def format(rgb):
        def clamp(x):
            return max(0, min(x, 255))
        red, green, blue = rgb
        return "#{0:02x}{1:02x}{2:02x}".format(clamp(int(red * 255)), clamp(int(green * 255)),
                                               clamp(int(blue * 255)))

def render_region(plot, xs, xe, ys, ye, pack, sv_db, run_id, ground_truth_id, min_score, max_num_ele, dataset_name, 
                  active_tools):
    plot.quad(left=0, bottom=0, right=pack.unpacked_size_single_strand, top=pack.unpacked_size_single_strand, 
              fill_alpha=0, line_color="black", line_width=3)
    lengths = pack.contigLengths()
    names = pack.contigNames()
    plot.axis.formatter = FuncTickFormatter(
        args={"lengths":[x for x in lengths], "names":[x for x in names]}, 
        code="""
                var i = 0;
                while(tick > lengths[i])
                {
                    tick -= lengths[i];
                    i += 1;
                    if(i >= lengths.length)
                        return tick + " - " + names[lengths.length-1];
                }
                return tick + " - " + names[i];
            """)
    if not sv_db.run_exists(run_id):
        return plot, True
    
    rendered_everything = False

    if xs < 0:
        xs = 0
    if ys < 0:
        ys = 0
    if xe < 0:
        xe = 0
    if ye < 0:
        ye = 0
    w = int(xe - xs)
    h = int(ye - ys)

    s = max(min(xs - w, ys - h), 0)
    e = min(max(xe + w, ye + h), pack.unpacked_size_single_strand)
    plot.line(x=[s,e], y=[s,e], line_color="black", line_width=3)

    give_up_factor = 1000
    if libMA.get_call_overview_area(sv_db, pack, run_id, min_score, int(xs - w), int(ys - h), w*3, h*3) > max_num_ele:
        plot.grid.visible = False
        div = int(math.sqrt(max_num_ele))
        rect_vec = libMA.get_call_overview(sv_db, pack, run_id, min_score, int(xs - w), int(ys - h), w*3, h*3,
                                           w//div, h//div, give_up_factor)

        cds = {
            'x': [],
            'y': [],
            'w': [],
            'h': [],
            'c': [],
            'f': [],
            't': [],
            'i': []
        }
        max_ = max(*[r.c for r in rect_vec], 0)
        for rect in rect_vec:
            cds["x"].append(rect.x)
            cds["y"].append(rect.y)
            cds["w"].append(rect.x + rect.w)
            cds["h"].append(rect.y + rect.h)
            cds["c"].append(format(light_spec_approximation(rect.c/max_)))
            cds["f"].append(names[rect.i])
            cds["t"].append(names[rect.j])
            cds["i"].append(str(rect.c))
        plot.quad(left="x", bottom="y", right="w", top="h", color="c", line_width=0, source=ColumnDataSource(cds),
                  name="hover1")

        url = "http://localhost:5006/bokeh_server?xs=@x&ys=@y&xe=@w&ye=@h&run_id=" + str(run_id) + \
            "&min_score=" + str(min_score) + "&ground_truth_id=" + str(ground_truth_id) + "&dataset_name=" + \
            dataset_name
        taptool = plot.select(type=TapTool)
        taptool.callback = OpenURL(url=url, same_tab=True)

        return plot, rendered_everything

    else:
        params = ParameterSetManager()
        accepted_boxes_data = {
            "x": [],
            "w": [],
            "y": [],
            "h": [],
            "n": [],
            "c": [],
            "r": [],
            "s": []
        }
        accepted_plus_data = {
            "x": [],
            "y": [],
            "n": [],
            "c": [],
            "r": [],
            "s": []
        }
        ground_plus_data = {
            "x": [],
            "y": [],
            "n": [],
            "c": [],
            "r": [],
            "s": []
        }
        calls_from_db = SvCallsFromDb(params, sv_db, run_id, int(xs - w), int(ys - h), w*3, h*3, min_score)
        while calls_from_db.hasNext():
            def score(jump):
                if jump.coverage == 0:
                    return None
                return str(jump.num_supp_nt / jump.coverage)
            jump = calls_from_db.next()
            if jump.num_supp_nt > min_score * jump.coverage:
                if jump.from_size == 1 and jump.to_size == 1:
                    accepted_plus_data["x"].append(jump.from_start)
                    accepted_plus_data["y"].append(jump.to_start)
                else:
                    accepted_boxes_data["x"].append(jump.from_start - 0.5)
                    accepted_boxes_data["y"].append(jump.to_start - 0.5)
                    accepted_boxes_data["w"].append(jump.from_start + jump.from_size + 1)
                    accepted_boxes_data["h"].append(jump.to_start + jump.to_size + 1)
                    accepted_boxes_data["n"].append(jump.num_supp_nt)
                    accepted_boxes_data["c"].append(jump.coverage)
                    accepted_boxes_data["r"].append(len(jump.supporing_jump_ids))
                    accepted_boxes_data["s"].append(score(jump))
                    accepted_plus_data["x"].append(jump.from_start + jump.from_size/2)
                    accepted_plus_data["y"].append(jump.to_start + jump.to_size/2)
                accepted_plus_data["n"].append(jump.num_supp_nt)
                accepted_plus_data["c"].append(jump.coverage)
                accepted_plus_data["r"].append(len(jump.supporing_jump_ids))
                accepted_plus_data["s"].append(score(jump))
        calls_from_db = SvCallsFromDb(params, sv_db, ground_truth_id, int(xs - w), int(ys - h), w*3, h*3, min_score)
        while calls_from_db.hasNext():
            def score(jump):
                if jump.coverage == 0:
                    return ""
                return " score: " + str(jump.num_supp_nt / jump.coverage)
            jump = calls_from_db.next()
            if jump.num_supp_nt > min_score * jump.coverage:
                if jump.from_size == 1 and jump.to_size == 1:
                    ground_plus_data["x"].append(jump.from_start + jump.from_size/2)
                    ground_plus_data["y"].append(jump.to_start + jump.to_size/2)
                    ground_plus_data["n"].append(jump.num_supp_nt)
                    ground_plus_data["c"].append(jump.coverage)
                    ground_plus_data["r"].append(len(jump.supporing_jump_ids))
                    ground_plus_data["s"].append(score(jump))
                else:
                    print("ground truth with fuzziness?!?!")
        
        num_jumps = libMA.get_num_jumps_in_area(sv_db, pack, sv_db.get_run_jump_id(run_id), int(xs - w), int(ys - h),
                                                w*3, h*3)
        if num_jumps < max_num_ele:
            out_dicts = []
            patch = {
                    "x": [],
                    "y": []
                }
            for _ in range(4):
                out_dicts.append({
                    "x": [],
                    "y": [],
                    "w": [],
                    "h": [],
                    "a": [],
                    "n": [],
                    "r": [],
                    "q": []
                })
            sweeper = SortedSvJumpFromSql(params, sv_db, sv_db.get_run_jump_id(run_id), int(xs - w), int(ys - h),
                                          w*3, h*3)
            while sweeper.has_next_start():
                jump = sweeper.get_next_start()
                idx = None
                if jump.switch_strand_known():
                    if jump.does_switch_strand():
                        idx = 0
                    else:
                        idx = 1
                else:
                    if jump.from_known():
                        idx = 2
                    else:
                        idx = 3
                
                out_dicts[idx]["x"].append( jump.from_start_same_strand() - 0.5 )
                out_dicts[idx]["y"].append( jump.to_start() - 0.5 )
                out_dicts[idx]["w"].append( jump.from_start_same_strand() + jump.from_size() + 1 )
                out_dicts[idx]["h"].append( jump.to_start() + jump.to_size() + 1 )
                out_dicts[idx]["a"].append( jump.num_supp_nt() / 1000 )
                out_dicts[idx]["n"].append(jump.num_supp_nt())
                out_dicts[idx]["r"].append(jump.read_id)
                out_dicts[idx]["q"].append(jump.query_distance())
                
                f = jump.from_pos
                t = jump.to_pos
                if not jump.from_known():
                    f = t
                if not jump.to_known():
                    t = f
                if not jump.from_fuzziness_is_rightwards():
                    if not jump.to_fuzziness_is_downwards():
                        patch["x"].extend([f - 2.5, f + .5, f + .5, float("NaN")])
                        patch["y"].extend([t - .5, t + 2.5, t - .5, float("NaN")])
                    else:
                        patch["x"].extend([f - 2.5, f + .5, f + .5, float("NaN")])
                        patch["y"].extend([t + .5, t - 2.5, t + .5, float("NaN")])
                else:
                    if not jump.to_fuzziness_is_downwards():
                        patch["x"].extend([f + 2.5, f - .5, f - .5, float("NaN")])
                        patch["y"].extend([t - .5, t + 2.5, t - .5, float("NaN")])
                    else:
                        patch["x"].extend([f + 2.5, f - .5, f - .5, float("NaN")])
                        patch["y"].extend([t + .5, t - 2.5, t + .5, float("NaN")])
                
            plot.quad(left="x", bottom="y", right="w", top="h", fill_color="orange", line_color="orange", line_width=3,
                      fill_alpha="a", source=ColumnDataSource(out_dicts[0]), name="hover3")
            plot.quad(left="x", bottom="y", right="w", top="h", fill_color="blue", line_color="blue", line_width=3,
                      fill_alpha="a", source=ColumnDataSource(out_dicts[1]), name="hover3")
            plot.quad(left="x", bottom="y", right="w", top="h", fill_color="grey", line_color="grey", line_width=3,
                      fill_alpha="a", source=ColumnDataSource(out_dicts[2]), name="hover3")
            plot.quad(left="x", bottom="y", right="w", top="h", fill_color="yellow", line_color="yellow", line_width=3,
                      fill_alpha="a", source=ColumnDataSource(out_dicts[3]), name="hover3")
            plot.patch(x="x", y="y", line_width=1, color="black", source=ColumnDataSource(patch))
            rendered_everything = True
        # the sv - boxes
        plot.quad(left="x", bottom="y", right="w", top="h", line_color="magenta", line_width=3, fill_alpha=0,
                  source=ColumnDataSource(accepted_boxes_data), name="hover2")
        plot.x(x="x", y="y", size=20, line_width=3, line_alpha=0.5, color="green",
               source=ColumnDataSource(ground_plus_data), name="hover2")
        plot.x(x="x", y="y", size=20, line_width=3, line_alpha=0.5, color="magenta",
               source=ColumnDataSource(accepted_plus_data), name="hover2")

    return rendered_everything

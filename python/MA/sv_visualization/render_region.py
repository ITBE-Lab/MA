from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import FuncTickFormatter, TapTool, OpenURL, LabelSet
from bokeh.models.callbacks import CustomJS
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

def render_region(plot, l_plot, d_plot, xs, xe, ys, ye, pack, sv_db, run_id, ground_truth_id, min_score, max_num_ele,
                  dataset_name, active_tools, index_prefix):
    plot.quad(left=0, bottom=0, right=pack.unpacked_size_single_strand, top=pack.unpacked_size_single_strand, 
              fill_alpha=0, line_color="black", line_width=3)
    lengths = pack.contigLengths()
    names = pack.contigNames()
    #plot.axis.formatter = FuncTickFormatter(
    #    args={"lengths":[x for x in lengths], "names":[x for x in names]}, 
    #    code="""
    #            var i = 0;
    #            while(tick > lengths[i])
    #            {
    #                tick -= lengths[i];
    #                i += 1;
    #                if(i >= lengths.length)
    #                    return tick + " - " + names[lengths.length-1];
    #            }
    #            return tick + " - " + names[i];
    #        """)
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
            read_ids = set()
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
                    "q": [],
                    "f": [],
                    "t": [],
                    "c": [],
                    "i": []
                })
            sweeper = SortedSvJumpFromSql(params, sv_db, sv_db.get_run_jump_id(run_id), int(xs - w), int(ys - h),
                                          w*3, h*3)
            while sweeper.has_next_start():
                jump = sweeper.get_next_start()
                idx = None
                if jump.switch_strand_known():
                    if jump.does_switch_strand():
                        idx = 0
                        out_dicts[idx]["c"].append( "orange" )
                    else:
                        idx = 1
                        out_dicts[idx]["c"].append( "blue" )
                else:
                    if jump.from_known():
                        idx = 2
                        out_dicts[idx]["c"].append( "lightgreen" )
                    else:
                        idx = 3
                        out_dicts[idx]["c"].append( "yellow" )
                
                out_dicts[idx]["f"].append( jump.from_pos )
                out_dicts[idx]["t"].append( jump.to_pos )
                out_dicts[idx]["x"].append( jump.from_start_same_strand() - 0.5 )
                out_dicts[idx]["y"].append( jump.to_start() - 0.5 )
                out_dicts[idx]["w"].append( jump.from_start_same_strand() + jump.from_size() + 1 )
                out_dicts[idx]["h"].append( jump.to_start() + jump.to_size() + 1 )
                out_dicts[idx]["a"].append( jump.num_supp_nt() / 1000 )
                out_dicts[idx]["n"].append( jump.num_supp_nt() )
                out_dicts[idx]["r"].append( jump.read_id )
                out_dicts[idx]["q"].append( jump.query_distance() )
                out_dicts[idx]["i"].append( jump.id )
                read_ids.add( jump.read_id )

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
            quads = []
            quads.append(plot.quad(left="x", bottom="y", right="w", top="h", fill_color="c", line_color="c",
                             line_width=3, fill_alpha="a", source=ColumnDataSource(out_dicts[0]), name="hover3"))
            quads.append(plot.quad(left="x", bottom="y", right="w", top="h", fill_color="c", line_color="c",
                                    line_width=3, fill_alpha="a", source=ColumnDataSource(out_dicts[1]), name="hover3"))
            quads.append(plot.quad(left="x", bottom="y", right="w", top="h", fill_color="c", line_color="c",
                                line_width=3, fill_alpha="a", source=ColumnDataSource(out_dicts[2]), name="hover3"))
            quads.append(plot.quad(left="x", bottom="y", right="w", top="h", fill_color="c", line_color="c",
                                    line_width=3, fill_alpha="a", source=ColumnDataSource(out_dicts[3]), name="hover3"))
            plot.patch(x="x", y="y", line_width=1, color="black", source=ColumnDataSource(patch))
            if len(read_ids) < max_num_ele:
                params = libMA.ParameterSetManager()
                seeder = BinarySeeding(params)
                jumps_from_seeds = libMA.SvJumpsFromSeeds(params, -1, sv_db, pack)
                fm_index = FMIndex()
                fm_index.load(index_prefix)
                read_dict = {
                    "center": [],
                    "r_id": [],
                    "size": [],
                    "q": [],
                    "r": [],
                    "idx": [],
                    "c": [],
                    "f": []
                }
                read_id_n_cols = []
                for read_id in read_ids:
                    read = sv_db.get_read(read_id)
                    segments = seeder.execute(fm_index, read)
                    seeds = libMA.Seeds()
                    jumps_from_seeds.cpp_module.execute_helper(segments, pack, fm_index, read, seeds)
                    end_column = []
                    seeds_n_idx = list(enumerate(sorted([x for x in seeds], key=lambda x: x.start)))
                    for idx, seed in sorted(seeds_n_idx, key=lambda x: x[1].start_ref):
                        if seed.on_forward_strand:
                            read_dict["center"].append(seed.start_ref + seed.size/2)
                            read_dict["c"].append("green")
                            read_dict["r"].append(seed.start_ref)
                            curr_end = seed.start_ref + seed.size + 3
                            curr_start = seed.start_ref
                        else:
                            read_dict["center"].append(seed.start_ref - seed.size/2 + 1)
                            read_dict["c"].append("purple")
                            read_dict["r"].append(seed.start_ref - seed.size + 1)
                            curr_end = seed.start_ref + 3
                            curr_start = seed.start_ref - seed.size
                        curr_column = 0
                        while curr_column < len(end_column):
                            if curr_start > end_column[curr_column]:
                                break
                            else:
                                curr_column += 1
                        if curr_column >= len(end_column):
                            read_id_n_cols.append( (read_id, curr_column) )
                            end_column.append(0)
                        end_column[curr_column] = curr_end
                        read_dict["r_id"].append( (str(read_id), str(curr_column)) )
                        read_dict["size"].append(seed.size)
                        read_dict["q"].append(seed.start)
                        read_dict["idx"].append(idx)
                        read_dict["f"].append(seed.on_forward_strand)
                if len(read_dict["c"]) < max_num_ele:
                    sorted_r_ids = [(str(x), str(y)) for x, y in sorted(read_id_n_cols, reverse=True)]
                    l_plot[1].x_range.factors = sorted_r_ids
                    l_plot[1].x_range.bounds = "auto"
                    d_plot[1].y_range.factors = sorted_r_ids
                    d_plot[1].y_range.bounds = "auto"
                    read_source = ColumnDataSource(read_dict)
                    l_plot[1].rect(x='r_id', y="center", width=1, height="size",
                                    fill_color="c", line_width=0, source=read_source, name="hover5")
                    labels = LabelSet(x='r_id', y='center', text='idx',
                                    source=read_source, render_mode='canvas', text_align="center",
                                    text_baseline="middle", text_color="white")
                    l_plot[1].add_layout(labels)
                    d_plot[1].rect(y='r_id', x="center", height=1, width="size",
                                  fill_color="c", line_width=0, source=read_source, name="hover5")
                    labels = LabelSet(x='center', y='r_id', text='idx',
                                    source=read_source, render_mode='canvas', text_align="center",
                                    text_baseline="middle", text_color="white")
                    d_plot[1].add_layout(labels)

                    num_nt = w*3+h*3
                    if num_nt < max_num_ele:
                        l_plot_nucs = {"y":[], "c":[], "i":[]}
                        def append_nuc_type(dict_, nuc):
                            if nuc == "A" or nuc =="a":
                                dict_["c"].append("blue")
                                dict_["i"].append("A")
                            elif nuc == "C" or nuc =="c":
                                dict_["c"].append("red")
                                dict_["i"].append("C")
                            elif nuc == "G" or nuc =="g":
                                dict_["c"].append("green")
                                dict_["i"].append("G")
                            elif nuc == "T" or nuc =="t":
                                dict_["c"].append("yellow")
                                dict_["i"].append("T")
                            else:
                                dict_["c"].append("lightgreen")
                                dict_["i"].append(nuc)
                        nuc_seq = pack.extract_from_to(int(ys - h), int(ye + h + 1))
                        for y_add, nuc in enumerate(str(nuc_seq)):
                            l_plot_nucs["y"].append(int(ys - h) + y_add + 0.5)
                            append_nuc_type(l_plot_nucs, nuc)
                        l_plot[0].rect(x=0.5, y="y", width=1, height=1, fill_color="c", line_width=0,
                                source=ColumnDataSource(l_plot_nucs), name="hover4")
                        d_plot_nucs = {"x":[], "c":[], "i":[]}
                        nuc_seq = pack.extract_from_to(int(xs - w), int(xe + w + 1))
                        for x_add, nuc in enumerate(str(nuc_seq)):
                            d_plot_nucs["x"].append(int(xs - w) + x_add + 0.5)
                            append_nuc_type(d_plot_nucs, nuc)
                        d_plot[0].rect(x="x", y=0.5, width=1, height=1, fill_color="c", line_width=0,
                                source=ColumnDataSource(d_plot_nucs), name="hover4")

                        # the tapping callback on jumps
                        plot.js_on_event("tap", CustomJS(args=dict(srcs=[x.data_source for x in quads],
                                                                   read_source=read_source),
                                                         code="""
                    var foud_one = false;
                    for(var i = 0; i < srcs.length; i++)
                    {
                        src = srcs[i];
                        for(var idx = 0; idx < src.data.a.length; idx++)
                        {
                            if(src.data.x[idx] <= cb_obj.x && src.data.w[idx] >= cb_obj.x &&
                                src.data.y[idx] <= cb_obj.y && src.data.h[idx] >= cb_obj.y)
                            {
                                src.data.c[idx] = ["orange", "blue", "grey", "yellow"][i];
                                for(var j = 0; j < read_source.data.r_id.length; j++)
                                    if(read_source.data.r_id[j][0] == src.data.r[idx].toString() &&
                                        (   read_source.data.r[j] == src.data.f[idx] || 
                                            read_source.data.r[j] == src.data.t[idx] || 
                                            read_source.data.r[j] + read_source.data.size[j] - 1 == src.data.f[idx] || 
                                            read_source.data.r[j] + read_source.data.size[j] - 1 == src.data.t[idx]
                                        )
                                            )
                                        read_source.data.c[j] = read_source.data.f[j] ? "green" : "purple";
                                    else
                                        read_source.data.c[j] = "lightgrey";
                                foud_one = true;
                            }
                            else
                                src.data.c[idx] = "lightgrey";
                        }
                    }
                    if(!foud_one)
                    {
                        for(var j = 0; j < read_source.data.r_id.length; j++)
                            read_source.data.c[j] = read_source.data.f[j] ? "green" : "purple";
                        for(var i = 0; i < srcs.length; i++)
                        {
                            src = srcs[i];
                            for(var idx = 0; idx < src.data.a.length; idx++)
                                src.data.c[idx] = ["orange", "blue", "lightgreen", "yellow"][i];
                        }
                    }
                    read_source.change.emit();
                    for(var i = 0; i < srcs.length; i++)
                        srcs[i].change.emit();
                                """))
                        # the tapping callback on seeds
                        code = """
        for(var i = 0; i < srcs.length; i++)
            for(var idx = 0; idx < srcs[i].data.a.length; idx++)
                srcs[i].data.c[idx] = ["orange", "blue", "lightgreen", "green"][i];
        for(var read_id in range._mapping)
        {
            for(var col in range._mapping[read_id].mapping)
            {
                pos = range._mapping[read_id].mapping[col].value;
                if(Math.abs(pos - cb_obj.y) <= 0.5)
                {
                    for(var j = 0; j < read_source.data.r_id.length; j++)
                    {
                        if(read_source.data.r_id[j][0] == read_id &&
                            Math.abs(read_source.data.center[j] - cb_obj.x) <= read_source.data.size[j]/2)
                        {
                            //console.log(read_id + " " + read_source.data.idx[j]);
                            
                            // find the appropriate jump
                            
                            for(var i = 0; i < srcs.length; i++)
                            {
                                src = srcs[i];
                                for(var idx = 0; idx < src.data.a.length; idx++)
                                {
                                    if(read_id == src.data.r[idx].toString() &&
                                        (   read_source.data.r[j] == src.data.f[idx] || 
                                            read_source.data.r[j] == src.data.t[idx] || 
                                            read_source.data.r[j] + read_source.data.size[j] - 1 == src.data.f[idx] || 
                                            read_source.data.r[j] + read_source.data.size[j] - 1 == src.data.t[idx]
                                        )
                                      )
                                        continue;
                                    else
                                        src.data.c[idx] = "lightgrey";
                                }
                                src.change.emit();
                            }
                            
                            for(var j2 = 0; j2 < read_source.data.r_id.length; j2++)
                            {
                                if(j2 == j)
                                    read_source.data.c[j2] = read_source.data.f[j] ? "green" : "purple";
                                else
                                    read_source.data.c[j2] = "lightgrey";
                            }
                            read_source.change.emit();
                            return;
                        }
                    }
                    for(var j = 0; j < read_source.data.r_id.length; j++)
                    {
                        if(read_source.data.r_id[j][0] == read_id)
                            read_source.data.c[j] = read_source.data.f[j] ? "green" : "purple";
                        else
                            read_source.data.c[j] = "lightgrey";
                    }
                    // correct column but no single seed matches...
                    for(var i = 0; i < srcs.length; i++)
                    {
                        src = srcs[i];
                        for(var idx = 0; idx < src.data.a.length; idx++)
                            if(read_id != src.data.r[idx].toString())
                                src.data.c[idx] = "lightgrey";
                        src.change.emit();
                    }
                    read_source.change.emit();
                    return;
                }
            }
        }
        for(var i = 0; i < srcs.length; i++)
            srcs[i].change.emit();
        return;
                        """
                        l_plot[1].js_on_event("tap", CustomJS(args=dict(srcs=[x.data_source for x in quads],
                                                                        read_source=read_source),
                                                              code=code))
                        d_plot[1].js_on_event("tap", CustomJS(args=dict(srcs=[x.data_source for x in quads],
                                                                        read_source=read_source,
                                                                        range=d_plot[1].y_range),
                                                              code=code))

                        rendered_everything = True
        # the sv - boxes
        plot.quad(left="x", bottom="y", right="w", top="h", line_color="magenta", line_width=3, fill_alpha=0,
                  source=ColumnDataSource(accepted_boxes_data), name="hover2")
        plot.x(x="x", y="y", size=20, line_width=3, line_alpha=0.5, color="green",
               source=ColumnDataSource(ground_plus_data), name="hover2")
        plot.x(x="x", y="y", size=20, line_width=3, line_alpha=0.5, color="magenta",
               source=ColumnDataSource(accepted_plus_data), name="hover2")


    return rendered_everything

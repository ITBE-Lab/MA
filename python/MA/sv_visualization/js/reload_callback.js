{
    if (typeof window.left_page_already !== 'undefined')
        return;
    else
        window.left_page_already = 1;
    if (typeof window.selected_read_id == 'undefined')
        window.selected_read_id = selected_read_id;
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
                                            "&range_link=" + checkbox_group.active +
                                            "&read_plot_start=" + read_plot_range.start +
                                            "&read_plot_end=" + read_plot_range.end +
                                            "&selected_read_id=" + window.selected_read_id +
                                            "&active_tools=" + active_tools;
    //alert(s);
    document.location.href = s;
} // if or scope
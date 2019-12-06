for (var j = 0; j < read_source.data.r_id.length; j++)
    read_source.data.c[j] = "lightgrey";

// reset read plot
for (var data_list_name in read_plot_line.data)
    read_plot_line.data[data_list_name] = [];
// reset y-axis nucleotides in read plot
for (var data_list_name in l_read_plot_data.data)
    l_read_plot_data.data[data_list_name] = [];
// reset rect data in read plot
for (var data_list_name in rect_read_plot_data.data)
    rect_read_plot_data.data[data_list_name] = [];
window.selected_read_id = -1;

var found_one = false;
for (var i = 0; i < srcs.length; i++)
{
    src = srcs[i];
    for (var idx = 0; idx < src.data.a.length; idx++)
    {
        if (src.data.x[idx] <= cb_obj.x && src.data.w[idx] >= cb_obj.x &&
            src.data.y[idx] <= cb_obj.y && src.data.h[idx] >= cb_obj.y)
        {
            src.data.c[idx] = ["orange", "blue", "lightgreen", "yellow"][i];
            for (var j = 0; j < read_source.data.r_id.length; j++)
                if (read_source.data.r_id[j] == src.data.r[idx] &&
                    (read_source.data.r[j] == src.data.f[idx] ||
                        read_source.data.r[j] == src.data.t[idx] ||
                        read_source.data.r[j] + read_source.data.size[j] - 1 == src.data.f[idx] ||
                        read_source.data.r[j] + read_source.data.size[j] - 1 == src.data.t[idx]
                    )
                )
                    read_source.data.c[j] = read_source.data.parlindrome[j] ? "red" : 
                                                                        (read_source.data.f[j] ? "green" : "purple");
            window.selected_read_id = src.data.r[idx];
            found_one = true;
        }
        else
            src.data.c[idx] = "lightgrey";
    }
}
if (!found_one)
{
    for (var j = 0; j < read_source.data.r_id.length; j++)
        read_source.data.c[j] = read_source.data.parlindrome[j] ? "red" : 
                                                                   (read_source.data.f[j] ? "green" : "purple");
    for (var i = 0; i < srcs.length; i++)
    {
        src = srcs[i];
        for (var idx = 0; idx < src.data.a.length; idx++)
            src.data.c[idx] = ["orange", "blue", "lightgreen", "yellow"][i];
    }
}
else
{
    for (var j = 0; j < read_source.data.r_id.length; j++)
        if (read_source.data.r_id[j] == window.selected_read_id)
            for (var data_list_name in read_source.data)
                read_plot_line.data[data_list_name].push(read_source.data[data_list_name][j]);
    // copy nucleotides over to read plot
    for (var data_list_name in l_read_plot_data.data)
        l_read_plot_data.data[data_list_name] =
            l_plot_nucs[window.selected_read_id][data_list_name];
    // copy rect data over
    for (var data_list_name in rect_read_plot_data.data)
        rect_read_plot_data.data[data_list_name] = read_plot_rects[window.selected_read_id][data_list_name];
}
read_source.change.emit();
for (var i = 0; i < srcs.length; i++)
    srcs[i].change.emit();
read_plot_line.change.emit();
l_read_plot_data.change.emit();
rect_read_plot_data.change.emit();
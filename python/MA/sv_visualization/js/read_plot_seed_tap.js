



// no seed hit?
// just mark everything that has to do with the read
if(read_plot_line.data.r_id.length > 0)
{
    r_id = read_plot_line.data.r_id[0];

    // set sv jumps back to correct colors
    for (var i = 0; i < srcs.length; i++)
        for (var idx = 0; idx < srcs[i].data.r.length; idx++)
            if(srcs[i].data.r[idx] == r_id)
                srcs[i].data.c[idx] = ["orange", "blue", "lightgreen", "green"][i];
            else
                srcs[i].data.c[idx] = "lightgrey";

    // set seeds to correct colors
    for (var j = 0; j < read_source.data.r_id.length; j++)
    {
        if (read_source.data.r_id[j] == r_id)
            read_source.data.c[j] = read_source.data.parlindrome[j] ? "red" :
                (read_source.data.f[j] ? "green" : "purple");
        else
            read_source.data.c[j] = "lightgrey";
    }
    for (var j = 0; j < read_plot_line.data.r_id.length; j++)
        if (read_plot_line.data.r_id[j] == r_id)
            read_plot_line.data.c[j] = read_plot_line.data.parlindrome[j] ? "red" :
                (read_plot_line.data.f[j] ? "green" : "purple");
} // if

for (var i = 0; i < srcs.length; i++)
    srcs[i].change.emit();
read_source.change.emit();
read_plot_line.change.emit();
rect_read_plot_data.change.emit();
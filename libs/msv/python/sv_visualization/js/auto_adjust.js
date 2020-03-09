function auto_adjust()
{
    if (radio_group.active == 0)
    {
        read_plot.x_range.start = plot.x_range.start;
        read_plot.x_range.end = plot.x_range.end;
    }
    if (radio_group.active == 1)
    {
        read_plot.x_range.start = plot.y_range.start;
        read_plot.x_range.end = plot.y_range.end;
    }
    read_plot.x_range.change.emit();

    var min_seed = 10000000;
    var max_seed = 0;
    for (var i = 0; i < read_plot_line.data.x.length; i++)
    {
        if (
            (read_plot_line.data.x[i][0] >= read_plot.x_range.start &&
                read_plot_line.data.x[i][0] <= read_plot.x_range.end)
            ||
            (read_plot_line.data.x[i][1] >= read_plot.x_range.start &&
                read_plot_line.data.x[i][1] <= read_plot.x_range.end)
        )
        {
            min_seed = Math.min(min_seed, read_plot_line.data.y[i][0]);
            min_seed = Math.min(min_seed, read_plot_line.data.y[i][1]);
            max_seed = Math.max(max_seed, read_plot_line.data.y[i][0]);
            max_seed = Math.max(max_seed, read_plot_line.data.y[i][1]);
        } // if
    } // for

    var size = max_seed - min_seed;

    read_plot.y_range.start = min_seed - size / 20;
    read_plot.y_range.end = max_seed + size / 20;
    read_plot.y_range.change.emit();
}
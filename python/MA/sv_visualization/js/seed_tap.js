for (var data_list_name in read_plot_line.data)
    read_plot_line.data[data_list_name] = [];
for (var i = 0; i < srcs.length; i++)
    for (var idx = 0; idx < srcs[i].data.a.length; idx++)
        srcs[i].data.c[idx] = ["orange", "blue", "lightgreen", "green"][i];
for (var j = 0; j < read_source.data.r_id.length; j++) {
    if (Math.abs(read_source.data.category[j] - curr_y) <= 1 / 2 &&
        Math.abs(read_source.data.center[j] - curr_x) <= read_source.data.size[j] / 2) {
        //console.log(read_id + " " + read_source.data.idx[j]);

        // find the appropriate jump

        for (var i = 0; i < srcs.length; i++) {
            src = srcs[i];
            for (var idx = 0; idx < src.data.a.length; idx++) {
                if (read_source.data.r_id[j] == src.data.r[idx] &&
                    (read_source.data.r[j] == src.data.f[idx] ||
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

        for (var j2 = 0; j2 < read_source.data.r_id.length; j2++) {
            if (j2 == j)
                read_source.data.c[j2] = read_source.data.f[j] ? "green" : "purple";
            else
                read_source.data.c[j2] = "lightgrey";
            // copy over to read viewer
            if (read_source.data.r_id[j] == read_source.data.r_id[j2])
                for (var data_list_name in read_source.data)
                    read_plot_line.data[data_list_name].push(read_source.data[data_list_name][j2]);
        }
        read_source.change.emit();
        read_plot_line.change.emit();
        window.selected_read_id = read_source.data.r_id[j];
        auto_adjust();
        return;
    }
}
for (var outer_j = 0; outer_j < read_source.data.r_id.length; outer_j++) {
    // correct column but no single seed matches...
    if (Math.abs(read_source.data.category[outer_j] - curr_y) <= 1 / 2) {
        for (var j = 0; j < read_source.data.r_id.length; j++) {
            if (read_source.data.r_id[j] == read_source.data.r_id[outer_j]) {
                read_source.data.c[j] = read_source.data.f[j] ? "green" : "purple";
                // copy over to read viewer
                for (var data_list_name in read_source.data)
                    read_plot_line.data[data_list_name].push(read_source.data[data_list_name][j]);
            }
            else
                read_source.data.c[j] = "lightgrey";
        }
        for (var i = 0; i < srcs.length; i++) {
            src = srcs[i];
            for (var idx = 0; idx < src.data.a.length; idx++) {
                if (read_source.data.r_id[outer_j] != src.data.r[idx])
                    src.data.c[idx] = "lightgrey";
            }
            src.change.emit();
        }
        read_source.change.emit();
        //read_plot.reset.emit();
        read_plot_line.change.emit();
        auto_adjust();
        window.selected_read_id = read_source.data.r_id[outer_j];
        return;
    }
}
for (var i = 0; i < srcs.length; i++)
    srcs[i].change.emit();
read_source.change.emit();
read_plot_line.change.emit();
auto_adjust();
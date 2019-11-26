for (var j = 0; j < read_source.data.r_id.length; j++)
    read_source.data.c[j] = "lightgrey";
var found_one = false;
for (var i = 0; i < srcs.length; i++) {
    src = srcs[i];
    for (var idx = 0; idx < src.data.a.length; idx++) {
        if (src.data.x[idx] <= cb_obj.x && src.data.w[idx] >= cb_obj.x &&
            src.data.y[idx] <= cb_obj.y && src.data.h[idx] >= cb_obj.y) {
            src.data.c[idx] = ["orange", "blue", "grey", "yellow"][i];
            for (var j = 0; j < read_source.data.r_id.length; j++)
                if (read_source.data.r_id[j] == src.data.r[idx] &&
                    (read_source.data.r[j] == src.data.f[idx] ||
                        read_source.data.r[j] == src.data.t[idx] ||
                        read_source.data.r[j] + read_source.data.size[j] - 1 == src.data.f[idx] ||
                        read_source.data.r[j] + read_source.data.size[j] - 1 == src.data.t[idx]
                    )
                )
                    read_source.data.c[j] = read_source.data.f[j] ? "green" : "purple";
            found_one = true;
        }
        else
            src.data.c[idx] = "lightgrey";
    }
}
if (!found_one) {
    for (var j = 0; j < read_source.data.r_id.length; j++)
        read_source.data.c[j] = read_source.data.f[j] ? "green" : "purple";
    for (var i = 0; i < srcs.length; i++) {
        src = srcs[i];
        for (var idx = 0; idx < src.data.a.length; idx++)
            src.data.c[idx] = ["orange", "blue", "lightgreen", "yellow"][i];
    }
}
read_source.change.emit();
for (var i = 0; i < srcs.length; i++)
    srcs[i].change.emit();
window.selected_read_id = -1;
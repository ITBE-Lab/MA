from bokeh.plotting import figure, output_file, show
from bokeh.layouts import row, column
from bokeh.palettes import d3
from bokeh.io import export_png


resolution = 800
min_x = -1
min_y = -1
max_x = 11
max_y = 11


plot = figure(
            title="seed shadows",
            plot_width=resolution, plot_height=resolution,
            x_range=(min_x+1,max_x-1),
            y_range=(min_y+1,max_y-1)
        )
plot.axis.visible = False
plot.grid.grid_line_color = None

# x y size draw_shadow?
seeds = [(2,3,2, True), (5,1,3, False), (6,8,0, True)]

seeds_x = []
seeds_y = []
patch_x = []
patch_y = []
for x, y, size, shadow in seeds:
    seeds_x.append(x)
    seeds_x.append(x + size)
    seeds_x.append(float('nan'))
    seeds_y.append(y)
    seeds_y.append(y + size)
    seeds_y.append(float('nan'))

    if shadow:
        patch_x.append(min_x)
        patch_y.append(y)
        patch_x.append(x)
        patch_y.append(y)
        patch_x.append(x + size)
        patch_y.append(y + size)
        patch_x.append(x + size)
        patch_y.append(max_y)
        patch_x.append(min_x)
        patch_y.append(max_y)
        patch_x.append(float('nan'))
        patch_y.append(float('nan'))
        patch_x.append(x)
        patch_y.append(min_y)
        patch_x.append(x)
        patch_y.append(y)
        patch_x.append(x + size)
        patch_y.append(y + size)
        patch_x.append(max_x)
        patch_y.append(y + size)
        patch_x.append(max_x)
        patch_y.append(min_y)
        patch_x.append(float('nan'))
        patch_y.append(float('nan'))


plot.patch(
        patch_x,
        patch_y,
        fill_color="#C0C0C0",
        fill_alpha=.5,
        line_color="#555555",
        line_width=2,
        line_dash=[10,10],
        legend="shadows"
    )

plot.line(
        seeds_x,
        seeds_y,
        legend="seeds",
        color="black",
        line_width=5
    )

plot.x(
        [seeds[-1][0]],
        [seeds[-1][1]],
        color="black",
        size=10,
        line_width=3,
        legend="P"
    )

show(plot)
from bokeh.plotting import figure, output_file, show
from bokeh.layouts import row, column
from bokeh.models import Arrow, OpenHead, NormalHead, VeeHead
from bokeh.palettes import d3
from bokeh.io import export_png, export_svgs
from bokeh.models import FuncTickFormatter, FixedTicker, Label
import math

font = "Times New Roman"

greys = [
        "#acacac",
        "#b9b9b9",
        "#c6c6c6",
        "#d3d3d3",
        "#dfdfdf",
        "#ececec",
        "#f9f9f9"
    ]

def save(plot, name):
    export_png(plot, filename="paperGraphics/" + name + ".png")
    plot.output_backend = "svg"
    export_svgs(plot, filename="paperGraphics/" + name + ".svg")
    #show(plot)

resolution = 300
min_x = 0
min_y = 0
max_x = 10
max_y = 10


plot = figure(
            title="Figure 2",
            plot_width=resolution, plot_height=resolution,
            x_range=(min_x,max_x),
            y_range=(min_y,max_y)
        )
plot.axis.visible = False
plot.grid.grid_line_color = None
plot.toolbar.logo = None
plot.toolbar_location = None

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
        fill_color=greys[2],
        fill_alpha=.5,
        line_color=None,
        #line_width=2,
        #line_dash=[2,2],
        legend="seed-shadow"
    )

plot.line(
        seeds_x,
        seeds_y,
        legend="seed",
        color="black",
        line_width=5
    )

plot.line(
        [4,6,6,8],
        [0,2,8,10],
        legend="path",
        color="black",
        line_dash=[2,2],
        line_width=2
    )

plot.x(
        [seeds[-1][0]],
        [seeds[-1][1]],
        color="black",
        size=10,
        line_width=3,
    )

plot.add_layout(
        Label(x=6.3, y=7.5, text='P', text_font=font, text_color="black") 
    )
"""
plot.add_layout(Arrow(end=OpenHead(line_color="black",size=10),
    x_start=seeds[0][0] + seeds[0][2], y_start=seeds[0][1] + seeds[0][2], 
    x_end=seeds[0][0] + seeds[0][2]+1, y_end=seeds[0][1] + seeds[0][2]))

plot.add_layout(Arrow(end=OpenHead(line_color="black",size=10),
    x_start=seeds[0][0] + seeds[0][2], y_start=seeds[0][1] + seeds[0][2], 
    x_end=seeds[0][0] + seeds[0][2], y_end=seeds[0][1] + seeds[0][2]+1))
"""

plot.title.text_font=font
plot.legend.label_text_font=font
plot.legend.label_text_baseline="bottom"
plot.axis.axis_label_text_font=font
plot.axis.major_label_text_font=font
plot.legend.location = "top_left"
save(plot, "shadows")


I = "I"
M = "M"
MM = "\=M"
D = "D"


query =     [                         "C", "A", "C", "A", "T", "A", "T", "T" ]
reference = ["A", "G", "G", "A", "G", "C", "A",           "T", "T", "T", "T", "C", "A"]
"""
for index in range(len(query)):
    query[index] = query[index] + " (" + str(index) + ")"
for index in range(len(reference)):
    reference[index] = reference[index] + " (" + str(index) + ")"
"""

alignment = [D,D,D,D,D,M,M,I,I,M,MM,M,M,D,D]


plot = figure(
            title="Figure 1",
            plot_width=resolution, plot_height=resolution,
            x_axis_label = "reference", y_axis_label = "query"
        )

plot.xaxis.major_tick_line_color = None
plot.yaxis.major_tick_line_color = None

cur_x = -.5
cur_y = -.5
m_x = []
m_y = []
mm_x = []
mm_y = []
i_x = []
i_y = []
d_x = []
d_y = []

c_alignment = []
for symbol in alignment:
    if len(c_alignment) > 0 and c_alignment[-1][0] == symbol:
        c_alignment[-1] = (symbol, c_alignment[-1][1]+1)
    else:
        c_alignment.append( (symbol,1) )

for symbol, amount in c_alignment:
    if symbol == I:
        i_x.append([cur_x, cur_x])
        i_y.append([cur_y, cur_y+amount])
        cur_y += amount
    elif symbol == D:
        d_x.append([cur_x, cur_x+amount])
        d_y.append([cur_y, cur_y])
        cur_x += amount
    elif symbol == M:
        m_x.append([cur_x, cur_x+amount])
        m_y.append([cur_y, cur_y+amount])
        cur_x += amount
        cur_y += amount
    elif symbol == MM:
        mm_x.append([cur_x, cur_x+amount])
        mm_y.append([cur_y, cur_y+amount])
        cur_x += amount
        cur_y += amount

plot.multi_line(
        m_x,
        m_y,
        legend="match",
        color="black",
        line_width=5
    )

plot.multi_line(
        mm_x,
        mm_y,
        legend="missmatch",
        color=greys[0],
        line_width=5
    )

plot.multi_line(
        i_x,
        i_y,
        legend="insertion",
        color="black",
        line_width=2,
        line_dash=[2,2]
    )

plot.multi_line(
        d_x,
        d_y,
        legend="deletion",
        color="black",
        line_width=2,
        line_dash=[5,5]
    )

plot.xaxis.ticker = FixedTicker(ticks=range(len(reference)))
plot.legend.location = "top_left"
plot.toolbar.logo = None
plot.toolbar_location = None
grid = []
for p in range(-1,len(reference)):
    grid.append(p+.5)
plot.xgrid.ticker = FixedTicker(ticks=grid)
plot.xgrid.band_fill_color = greys[3]
plot.xgrid.band_fill_alpha = 0.2
plot.xgrid.grid_line_color = greys[0]
plot.xaxis.formatter = FuncTickFormatter(code="""
    var labels = %s;
    return labels[tick];
""" % reference)
plot.yaxis.ticker = FixedTicker(ticks=range(len(query)))
grid = []
for p in range(-1,len(query)):
    grid.append(p+.5)
plot.ygrid.ticker = FixedTicker(ticks=grid)
plot.ygrid.band_fill_color = greys[3]
plot.ygrid.band_fill_alpha = 0.2
plot.ygrid.grid_line_color = greys[0]
plot.yaxis.formatter = FuncTickFormatter(code="""
    var labels = %s;
    return labels[tick];
""" % query)


plot.title.text_font=font
plot.legend.label_text_font=font
plot.legend.label_text_baseline="bottom"
plot.axis.axis_label_text_font=font
plot.axis.major_label_text_font=font
save(plot, "alignment")


plot = figure(
            title="Figure 3: Strip of consideration",
            plot_width=resolution, plot_height=resolution,
            x_axis_label = "reference", y_axis_label = "query"
        )

plot.patch(
        [-.5,7.5,12.5,12.5,5.5],
        [-.5,7.5,7.5,6.5,-.5],
        fill_color=greys[2],
        fill_alpha=.75,
        line_color=None,
        #line_width=2,
        #line_dash=[2,2],
        legend="strip of consideration"
    )


plot.line(
    [-.5,7.5],
    [-.5,7.5],
    color="black",
    legend="strip of consideration",
    line_width=1,
    line_dash=[2,2]
)

plot.line(
    [5.5,12.5],
    [-.5,6.5],
    color="black",
    legend="strip of consideration",
    line_width=1,
    line_dash=[2,2]
)

plot.line(
    [5.5,7.5],
    [2.5,4.5],
    color="black",
    legend="anchor",
    line_width=5
)


plot.xaxis.major_tick_line_color = None
plot.yaxis.major_tick_line_color = None
plot.xaxis.ticker = FixedTicker(ticks=range(len(reference)))
plot.legend.location = "top_left"
plot.toolbar.logo = None
plot.toolbar_location = None
grid = []
for p in range(-1,len(reference)):
    grid.append(p+.5)
plot.xgrid.ticker = FixedTicker(ticks=grid)
plot.xgrid.band_fill_color = greys[3]
plot.xgrid.band_fill_alpha = 0.2
plot.xgrid.grid_line_color = greys[0]
plot.xaxis.formatter = FuncTickFormatter(code="""
    var labels = %s;
    return labels[tick];
""" % reference)
plot.yaxis.ticker = FixedTicker(ticks=range(len(query)))
grid = []
for p in range(-1,len(query)):
    grid.append(p+.5)
plot.ygrid.ticker = FixedTicker(ticks=grid)
plot.ygrid.band_fill_color = greys[3]
plot.ygrid.band_fill_alpha = 0.2
plot.ygrid.grid_line_color = greys[0]
plot.yaxis.formatter = FuncTickFormatter(code="""
    var labels = %s;
    return labels[tick];
""" % query)

plot.title.text_font=font
plot.legend.label_text_font=font
plot.legend.label_text_baseline="hanging"
plot.axis.axis_label_text_font=font
plot.axis.major_label_text_font=font
save(plot, "stripOfConsideration")
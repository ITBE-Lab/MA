from bokeh.plotting import figure, output_file, show
from bokeh.layouts import row, column, gridplot
from bokeh.models import Arrow, OpenHead, NormalHead, VeeHead
from bokeh.palettes import d3
from bokeh.io import export_png, export_svgs
from bokeh.models import FuncTickFormatter, FixedTicker, Label, ColorBar, FactorRange
from bokeh.models import LinearAxis, Range1d, LogColorMapper, FixedTicker, LinearColorMapper
from bokeh.models import ColumnDataSource, CompositeTicker
from bokeh.transform import dodge
from bokeh.core.properties import value
import math
import random
import numpy as np


font = "Helvetica"

dark_greys = [
        "#9f9f9f",
        "#939393",
        "#868686",
]

greys = [
        "#acacac",
        "#b9b9b9",
        "#c6c6c6",
        "#d3d3d3",
        "#dfdfdf",
        "#ececec",
        "#f9f9f9"
    ]

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

    return (i*r**m,i*g**m,i*b**m)
"""
for x in np.linspace(0, 1, 100):
    def format(rgb):
        def clamp(x):
            return int(max(0, min(x*255, 255)))
        r, g, b = rgb
        return (clamp(r),clamp(g),clamp(b))
    print(format(light_spec_approximation(x)))
"""

def heatmap_palette(scheme, num_colors):
    def format(rgb):
        def clamp(x):
            return max(0, min(x, 255))
        red, green, blue = rgb
        return "#{0:02x}{1:02x}{2:02x}".format(clamp(int(red * 255)), clamp(int(green * 255)),
                                               clamp(int(blue * 255)))
    return [format(scheme(x)) for x in np.linspace(0, 1, num_colors)]

def simulate_max_length(q_len, mutation_amount, indel_amount, indel_size, sim_amount):
    match_lens = []
    if mutation_amount + indel_amount/2 * indel_size >= q_len:
        return None
    for _ in range(sim_amount):
        q = []
        for index in range(q_len):
            q.append(index)
        ##
        # helper function
        # returns a list (amount elements) of random indices that are part of [0, interval_length]
        # and at least min_distance apart from each other
        #
        # if it is not possible to fit enough indices into [0, interval_length] the last indices 
        # are behind interval_length
        def get_random_spots(amount, interval_length, min_distance):
            spots = []
            for index in range(interval_length - min_distance):
                spots.append(index)
            random.shuffle(spots)
            spots = spots[:amount]
            spots.sort()
            for index in range(1,len(spots)):
                if spots[index-1] + min_distance + 1 > spots[index]:
                    spots[index] = spots[index-1] + min_distance + 1
            if len(spots) >=1 and spots[-1] > interval_length:
                start = spots[0]
                for index in range(len(spots)):
                    spots[index] -= start
            return spots

        deletion_amount = int(indel_amount/2)
        insertion_amount = int( (indel_amount+1) /2)

        # deletion
        deletion_spots = get_random_spots(deletion_amount, len(q), indel_size + 1)
        for pos in reversed(deletion_spots):
            q = q[:pos] + q[pos + indel_size:]

        # mutations
        mutation_spots = []
        for index in range(len(q)):
            mutation_spots.append(index)
        random.shuffle(mutation_spots)
        for pos in mutation_spots[:mutation_amount]:
            l = len(q)
            q = q[:pos] + [-2] + q[pos+1:]
            if not len(q) == l:
                print("ERROR: mutation changed query length:" + str(l) + " != " + str(len(q)))
                print(q)

        # insertion
        insertion_spots = get_random_spots(insertion_amount, len(q), 1)
        # pos are sorted in order to we need to reverse them in order 
        # to not insert twice at the same location
        for pos in reversed(insertion_spots):
            for _ in range(indel_size):
                q = q[:pos] + [-2] + q[pos:]

        # get the results
        matches = []
        last = -2
        match_len = 1
        for num in q:
            if num == last + 1:
                match_len += 1
            elif not num == -2:
                matches.append(match_len)
                match_len = 1
            last = num
        matches = matches[1:]
        matches.append(match_len)

        match_lens.append(matches)

    return match_lens

def only_max(li):
    ret = []
    if li is None:
        return None
    for l in li:
        ret.append(0)
        for x in l:
            if x > ret[-1]:
                ret[-1] = x
    return ret

def mean(li):
    if li is None:
        return float('NaN')
    return sorted(li)[int(len(li)/2)]

def avg(li):
    if li is None:
        return float('NaN')
    avg = 0
    for x in li:
        avg += x
    return avg / len(li)

def save(plot, name, grid=False):
    if grid:
        export_png(gridplot(plot), filename="paperGraphics/" + name + ".png")
        for r in plot:
            for p in r:
                p.output_backend = "svg"
        export_svgs(gridplot(plot), filename="paperGraphics/" + name + ".svg")
    else:
        export_png(plot, filename="paperGraphics/" + name + ".png")
        plot.output_backend = "svg"
        export_svgs(plot, filename="paperGraphics/" + name + ".svg")
    #show(plot)

resolution = 300
min_x = 0
min_y = 0
max_x = 10
max_y = 10

def ambiguity_per_length():
    data1 = [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1e-05, 3.0000000000000004e-05, 4e-05, 6e-05], [0.00035000000000000027, 0.0020700000000000046, 0.006419999999999907, 0.012899999999999643, 0.01866999999999941, 0.024689999999999164, 0.029309999999998976, 0.029509999999998968, 0.029589999999998964, 0.02994999999999895], [0.049500000000004485, 0.09312999999999659, 0.10157999999999331, 0.08540999999999958, 0.07046000000000538, 0.0524800000000054, 0.03988000000000154, 0.030989999999998907, 0.023459999999999214, 0.019469999999999377], [0.2835600000001428, 0.2028500000000621, 0.10839999999999067, 0.06604000000000709, 0.04710000000000375, 0.039120000000001307, 0.03392999999999972, 0.02893999999999899, 0.024359999999999177, 0.021269999999999303], [0.584869999999973, 0.19182000000005106, 0.08475999999999984, 0.047250000000003796, 0.02921999999999898, 0.018079999999999433, 0.011219999999999711, 0.00755999999999986, 0.005549999999999942, 0.0035300000000000084], [0.8230899999988889, 0.1222399999999853, 0.032189999999999185, 0.010929999999999723, 0.004100000000000001, 0.0025300000000000058, 0.0013400000000000029, 0.0008200000000000015, 0.0005400000000000008, 0.0004700000000000006], [0.9413299999983508, 0.049270000000004414, 0.006239999999999914, 0.0015000000000000033, 0.0007100000000000012, 0.00031000000000000016, 0.0001, 0.00013000000000000002, 7.000000000000001e-05, 2e-05], [0.9838899999981571, 0.014619999999999573, 0.000990000000000002, 0.00026000000000000003, 0.00011, 4e-05, 3.0000000000000004e-05, 0.0, 0.0, 0.0]]
    data2 = [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.06633000000000698, 0.1979700000000572, 0.15913000000001837, 0.10464999999999212, 0.24937000000010862, 0.16641000000002565, 0.04387000000000276], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0029000000000000067, 0.09215999999999697, 0.21072000000006996, 0.15462000000001386, 0.10385999999999243, 0.224770000000084, 0.15348000000001272, 0.04309000000000252, 0.01024999999999975, 0.003020000000000007], [0.00016, 0.0007000000000000012, 0.006189999999999916, 0.046100000000003444, 0.11030999999998993, 0.19924000000005848, 0.13397999999999322, 0.10748999999999102, 0.19883000000005807, 0.1395699999999988, 0.04340000000000262, 0.009899999999999765, 0.0025000000000000057, 0.000990000000000002, 0.0003700000000000003], [0.028599999999999005, 0.053940000000005844, 0.09152999999999721, 0.1085899999999906, 0.08126000000000119, 0.10789999999999086, 0.16968000000002892, 0.1186499999999867, 0.04092000000000186, 0.01098999999999972, 0.0028200000000000065, 0.001010000000000002, 0.0004300000000000005, 0.00016, 4e-05], [0.016279999999999507, 0.028739999999999, 0.05114000000000499, 0.08732999999999884, 0.11884999999998662, 0.08246000000000073, 0.0335399999999996, 0.010659999999999734, 0.002950000000000007, 0.001030000000000002, 0.00035000000000000027, 0.00013000000000000002, 8e-05, 4e-05, 0.0], [0.017829999999999444, 0.02891999999999899, 0.0358400000000003, 0.030989999999998907, 0.01866999999999941, 0.007619999999999858, 0.0028600000000000066, 0.0009400000000000018, 0.0004100000000000004, 0.00018, 8e-05, 5e-05, 1e-05, 3.0000000000000004e-05, 0.0], [0.0031800000000000075, 0.003710000000000009, 0.003800000000000009, 0.002590000000000006, 0.0014700000000000032, 0.0007200000000000012, 0.00035000000000000027, 0.00021, 5e-05, 5e-05, 0.0, 0.0, 0.0, 1e-05, 0.0], [0.00024, 0.00035000000000000027, 0.00038000000000000035, 0.0003300000000000002, 0.00022, 9e-05, 6e-05, 6e-05, 0.0, 1e-05, 1e-05, 0.0, 0.0, 0.0, 0.0], [1e-05, 4e-05, 6e-05, 0.0001, 5e-05, 4e-05, 1e-05, 0.0, 0.0, 0.0, 1e-05, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 2e-05, 0.0, 1e-05, 0.0, 2e-05, 1e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]

    r1max = 10
    r2size = 15
    min_len=10
    max_len=20

    color_mapper = LinearColorMapper(
                    palette=heatmap_palette(light_spec_approximation, 256),
                    low=0,
                    high=1
                )

    plot = figure(title="ambiguity on human genome",
            x_range=(0,r1max), y_range=(min_len, max_len),
            x_axis_label='ambiguity', y_axis_label='sequence length',
            plot_width=resolution*2, plot_height=resolution,
            min_border_bottom=10, min_border_top=10,
            min_border_left=10, min_border_right=15
        )
    plot.image(image=[data1], color_mapper=color_mapper,
            dh=[max_len - min_len], dw=[r1max], x=[0], y=[min_len])

    plot2 = figure(x_range=(r1max,2**r2size+r1max), y_range=(min_len, max_len),
            min_border_bottom=10, min_border_top=10,
            min_border_left=20, min_border_right=15,
            plot_width=resolution, plot_height=resolution,tools=[],
            x_axis_type="log"
        )
    plot2.image(image=[data2], color_mapper=color_mapper,
            dh=[max_len - min_len], dw=[2**r2size+r1max], x=[r1max], y=[min_len])

    color_bar = ColorBar(color_mapper=color_mapper, border_line_color=None, location=(0,0))
    color_bar.major_label_text_font=font
    plot.add_layout(color_bar, 'left')

    plot.legend.label_text_font=font
    plot.axis.axis_label_text_font=font
    plot.axis.major_label_text_font=font
    plot.xaxis.major_label_standoff = 15

    plot2.legend.label_text_font=font
    plot2.yaxis.visible = False
    plot2.axis.axis_label_text_font=font
    plot2.axis.major_label_text_font=font
    plot2.xaxis.major_label_standoff = 15

    save([[plot, plot2]], "ambiguityPerQueryLen", True)

def theoretical_max_acc():

    def binomial_no_success(n, p):
        return (1-p)**n

    """
    def prob_x_exists(x_len, ref_len):
        p = 0.25**x_len # (1/4)^x_len
        n = ref_len - x_len
        return 1-binomial_no_success(n,p)
    """

    def prob_non_enclosed(x_len, ref_len):
        expected_num_matches = (0.25**x_len) * (ref_len - x_len)
        prob_extendable = 1.0 - 3.0/4.0 ** 2
        return binomial_no_success(expected_num_matches, prob_extendable)

    def prob_all_non_enclosed(li, ref_len):
        ret = 1.0
        # compute weather all are enclosed
        for x in li:
            ret *= 1-prob_non_enclosed(x, ref_len)
        # return opposite
        return 1-ret

    def exp_match_amount(x_len, ref_len):
        p = 0.25**x_len
        n = ref_len - x_len
        return n*p

    ref_len = 3000000000 # three billion => human genome length
    q_len = 1000
    indel_size = 3
    prob_refindable = []
    quality = 100

    print("creating query length matrix...")
    max_indels = int(q_len/indel_size)*2
    for num_mut in range(0, q_len, max(10,int(q_len/quality))):
        prob_refindable.append( [] )
        for num_indel in range(0,max_indels, max(2,int(max_indels/quality))):
            q_len_e = simulate_max_length(q_len, num_mut, num_indel, indel_size, 512)
            if q_len_e is None:
                prob_refindable[-1].append(float('NaN'))
            else:
                probs = []
                for x in q_len_e:
                    probs.append(prob_all_non_enclosed(x, ref_len))
                prob_refindable[-1].append(avg(probs))
            #prob_refindable[-1].append(q_len_e)
        if num_mut % 100 == 0:
            print(num_mut, "/", q_len)
    print("done")

    w = q_len
    h = max_indels

    color_mapper = LinearColorMapper(
                        palette=heatmap_palette(light_spec_approximation, 127),
                        low=0,
                        high=1
                    )

    tick_formater = FuncTickFormatter(code="""
        return Math.max(Math.floor( (tick+1)/2),0) + '; ' +
                Math.max(Math.floor( (tick)/2),0)"""
        )
    #tick_formater = FuncTickFormatter(code="return 'a')

    plot = figure(title="theoretical max accuracy",
            x_range=(0,h), y_range=(0,w),
            x_axis_label='num ' + str(indel_size) + ' nt insertions; num ' + str(indel_size) + ' nt deletions', y_axis_label='num mutations'
        )
    plot.xaxis.formatter = tick_formater
    plot.image(image=[prob_refindable], color_mapper=color_mapper,
            dh=[w], dw=[h], x=[0], y=[0])

    color_bar = ColorBar(color_mapper=color_mapper, border_line_color=None, location=(0,0))

    plot.add_layout(color_bar, 'left')

    plot.title.text_font=font
    plot.legend.label_text_font=font
    plot.axis.axis_label_text_font=font
    plot.axis.major_label_text_font=font
    plot.legend.label_text_baseline="bottom"
    save(plot, "upperBoundShortIndels")


#static stuff
I = "I"
M = "M"
MM = "\=M"
D = "D"
query =     [                         "T", "A", "C", "A", "T", "T", "C", "T" ]
reference = ["T", "T", "C", "A", "G",           "C", "A", "T", "A", "C", "T", "C", "A"]
l_alignment = [D,D,I,I,M,M,D,D,M,MM,M,M,D,D]

def seed_shadows():
    min_x = -1.0
    min_y = -1.0
    max_x = 13
    max_y = 8
    plot = figure(
                title="Figure X: Shadows",
                plot_width=resolution, plot_height=resolution,
                x_range=[-1,13], y_range=[-1,8]
            )
    # x y size draw_shadow?
    seeds = [(1.5,1.5,2, False), (5.5,2.5,2, True), (1.5,5.5,1, True)]

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
            fill_color=greys[0],
            fill_alpha=.5,
            line_color=None,
            #line_width=2,
            #line_dash=[2,2],
        )
    """
    plot.patch(
            [seeds[-2][0]+seeds[-2][2], max_x, max_x, seeds[-2][0], seeds[-2][0]],
            [seeds[-2][1]+seeds[-2][2], seeds[-2][1]+seeds[-2][2], min_y, min_y, seeds[-2][1]],
            fill_color=dark_greys[2],
            fill_alpha=.5,
            line_color=None,
            #line_width=2,
            #line_dash=[2,2],
        )
    """

    plot.line(
            seeds_x,
            seeds_y,
            color="black",
            line_width=5
        )

    plot.line(
            [min_x, seeds[0][0], seeds[0][0]+1, seeds[1][0], seeds[1][0]+5.5],
            [min_y, seeds[0][1], seeds[0][1]+1, seeds[1][1], max_y],
            color="black",
            line_dash=[2,2],
            line_width=2
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
    save(plot, "shadows")

def alignment():
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
    for symbol in l_alignment:
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
            legend="ma",
            color="black",
            line_width=5
        )

    plot.multi_line(
            mm_x,
            mm_y,
            legend="mis",
            color=greys[0],
            line_width=5
        )

    plot.multi_line(
            i_x,
            i_y,
            legend="ins_____",
            color="black",
            line_width=2,
            line_dash=[2,2]
        )

    plot.multi_line(
            d_x,
            d_y,
            legend="del",
            color="black",
            line_width=2,
            line_dash=[5,5]
        )

    plot.xaxis.ticker = FixedTicker(ticks=range(len(reference)))
    plot.legend.location = "bottom_right"
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

def stripOfConsideration():
    plot = figure(
                title="Figure X: Strip of consideration",
                plot_width=resolution, plot_height=resolution,
                x_range=[-1,13], y_range=[-1,8]
            )
    plot.patch(
            [-.5,7.5,12.5,12.5,5.5],
            [-.5,7.5,7.5,6.5,-.5],
            fill_color=greys[2],
            fill_alpha=.75,
            line_color=None,
            #line_width=2,
            #line_dash=[2,2],
        )


    plot.line(
        [-.5,7.5],
        [-.5,7.5],
        color="black",
        line_width=1,
        line_dash=[2,2]
    )

    plot.line(
        [5.5,12.5],
        [-.5,6.5],
        color="black",
        line_width=1,
        line_dash=[2,2]
    )

    plot.line(
        [4.5,7.5],
        [1.5,4.5],
        color="black",
        line_width=6
    )
    plot.line(
        [2.5,4.5],
        [-.5,1.5],
        color="black",
        line_width=1,
        line_dash=[8,2]
    )
    #lines to right and left
    plot.line(
        [3,6],
        [3,3],
        color="black",
        line_width=1
    )
    plot.line(
        [6,9],
        [3,3],
        color="black",
        line_width=1
    )


    plot.line(
        [10.5,11.5],
        [1.5,2.5],
        color=dark_greys[2],
        line_width=3
    )
    plot.line(
        [8.5,10.5],
        [-.5,1.5],
        color=dark_greys[2],
        line_width=1,
        line_dash=[8,2]
    )

    plot.line(
        [1.5,2.5],
        [1.5,2.5],
        color=dark_greys[2],
        line_width=3
    )
    plot.line(
        [-.5,1.5],
        [-.5,1.5],
        color=dark_greys[2],
        line_width=1,
        line_dash=[8,2]
    )

    plot.line(
        [-.5,2.5],
        [3.5,6.5],
        color=dark_greys[2],
        line_width=3
    )
    plot.line(
        [-1.0,-.5],
        [3.0,3.5],
        color=dark_greys[2],
        line_width=1,
        line_dash=[8,2]
    )

    plot.line(
        [0.5,1.5],
        [-.5,0.5],
        color="black",
        line_width=3
    )

    plot.line(
        [9.5,11.5],
        [4.5,6.5],
        color="black",
        line_width=3
    )
    plot.line(
        [4.5,9.5],
        [-.5,4.5],
        color="black",
        line_width=1,
        line_dash=[8,2]
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


def optimal_matching():
    min_x = -1.0
    min_y = -1.0
    max_x = 13
    max_y = 8
    plot = figure(
                title="Figure X: Optimal matching",
                plot_width=resolution, plot_height=resolution,
                x_range=[-1,13], y_range=[-1,8]
            )
    plot.patch(
            [1.5,1.5,4.5,4.5],
            [.5,1.5,1.5,0.5],
            fill_color=greys[2],
            fill_alpha=.75,
            line_color=None,
            #line_width=2,
            #line_dash=[2,2],
        )
    plot.patch(
            [10.5,10.5,max_x,max_x],
            [5.5,max_y,max_y,5.5],
            fill_color=greys[2],
            fill_alpha=.75,
            line_color=None,
            #line_width=2,
            #line_dash=[2,2],
        )
    plot.patch(
            [min_x,min_x,0.5,0.5],
            [min_y,-.5,-.5,min_y],
            fill_color=greys[2],
            fill_alpha=.75,
            line_color=None,
            #line_width=2,
            #line_dash=[2,2],
        )

    plot.line(
        [0.5,1.5],
        [-.5,0.5],
        color="black",
        line_width=3
    )
    plot.line(
        [4.5,7.5],
        [1.5,4.5],
        color="black",
        line_width=3
    )
    plot.line(
        [9.5,10.5],
        [4.5,5.5],
        color="black",
        line_width=3
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
    save(plot, "optimalMatching")

def unrelated_non_enclosed_seeds():
    plot = figure(
                title="Supplementary Figure X: unrelated Seeds",
                plot_width=resolution, plot_height=resolution/2,
                y_range=["reference", "1", "2", "3", "query1", "query2", "query3", "query4"]
            )

    seeds = [
        ("query1", 0, 500, 20),
        ("query1", 24, 525, 13),
        ("query1", 47, 547, 24),

        ("query2", 22, 320, 5),
        ("query2", 36, 1000, 12),
        
        ("query3", 19, 1500, 4),
        ("query3", 35, 320, 3),
        ("query3", 46, 800, 5),

        ("query4", 20, 800, 4),
        ("query4", 45, 1340, 4)
    ]

    max_q = 0
    max_r = 0
    for seed in seeds:
        if max_q < seed[1] + seed[-1]:
            max_q = seed[1] + seed[-1]
        if max_r < seed[2] + seed[-1]:
            max_r = seed[2] + seed[-1]

    r_fac = max_q/float(max_r)
    for seed in seeds:
        plot.patch(
            [seed[1]+seed[-1], (seed[2]+seed[-1])*r_fac, seed[2]*r_fac, seed[1]],
            [seed[0], "reference", "reference", seed[0]],
            fill_color=greys[2],
            fill_alpha=.5,
            line_color=None,
        )
    plot.line([0, max_q], ["reference", "reference"], color="black", line_width=1)
    for seed in seeds:
        plot.line([seed[1], seed[1] + seed[-1]], [seed[0], seed[0]], color="black", line_width=1)



    plot.title.text_font=font
    plot.legend.label_text_font=font
    plot.legend.label_text_baseline="hanging"
    plot.axis.axis_label_text_font=font
    plot.axis.major_label_text_font=font
    plot.xaxis.major_tick_line_color = None
    plot.yaxis.major_tick_line_color = None
    plot.xaxis.minor_tick_line_color = None
    plot.yaxis.minor_tick_line_color = None
    plot.xaxis.major_label_text_alpha = 0
    plot.toolbar.logo = None
    plot.toolbar_location = None
    plot.grid.grid_line_color = None
    #plot.xgrid.ticker = FixedTicker(ticks=[-5,3,6,21,25,40])
    save(plot, "unrelatedNonEnclosedSeeds")

def required_nmw_band_size():
    data = ([0.429135585783255, 0.0013163668275559457, 0.00153576129881527, 0.0013163668275559457, 0.00021939447125932427, 0.0, 0.0006581834137779728, 0.002193944712593243, 0.0013163668275559457, 0.0026327336551118918, 0.0026327336551118918, 0.0035103115401491896, 0.0024133391838525673, 0.0039491004826678385, 0.005265467310223785, 0.004607283896445812, 0.0063624396665204076, 0.004826678367705136, 0.005704256252742434, 0.006581834137779732, 0.007020623080298381, 0.00548486178148311, 0.005923650724001759, 0.006581834137779732, 0.007678806494076354, 0.007678806494076354, 0.0068012286090390565, 0.008336989907854328, 0.008775778850372977, 0.006581834137779732, 0.008775778850372977, 0.009214567792891626, 0.009653356735410274, 0.013163668275559466, 0.010969723562966221, 0.013163668275559466, 0.012286090390522168, 0.011627906976744195, 0.01382185168933744, 0.011189118034225546, 0.009653356735410274, 0.013163668275559466, 0.011847301448003519, 0.016235190873190006, 0.014918824045634061, 0.01842913558578325, 0.010969723562966221, 0.015796401930671358, 0.011627906976744195, 0.20645019745502155], [0.9785914757525838, 0.00019580967299784609, 0.00019580967299784609, 9.790483649892306e-05, 9.790483649892306e-05, 9.790483649892306e-05, 9.790483649892306e-05, 0.00013053978199856407, 0.00019580967299784609, 0.00016317472749820508, 6.526989099928203e-05, 9.790483649892306e-05, 0.00019580967299784609, 9.790483649892306e-05, 0.0002284446184974871, 0.00019580967299784609, 0.00019580967299784609, 0.0002284446184974871, 3.2634945499641017e-05, 0.00016317472749820508, 0.00026107956399712813, 0.0002284446184974871, 0.00032634945499641015, 0.00013053978199856407, 0.0004568892369949742, 0.00029371450949676914, 0.00026107956399712813, 0.00032634945499641015, 0.00035898440049605116, 0.0005221591279942563, 0.0004242542914953332, 0.0005547940734938973, 0.00039161934599569217, 0.00039161934599569217, 0.00039161934599569217, 0.0005221591279942563, 0.0006200639644931795, 0.0004242542914953332, 0.0006526989099928205, 0.00032634945499641015, 0.0005221591279942563, 0.0005221591279942563, 0.0005547940734938973, 0.0006853338554924616, 0.00039161934599569217, 0.0006853338554924616, 0.0005221591279942563, 0.0004895241824946152, 0.0004242542914953332, 0.005972195026434325])

    num_buckets = len(data[0])
    if len(data[1]) != num_buckets:
        print("WARNING LENGTHS DIFFER")
    w = 1.0 / (num_buckets)
    buckets_x = []
    for index in range(num_buckets):
        buckets_x.append(index/float(num_buckets))
    print(buckets_x)
    print(w)

    plot = figure(
            title="Figure X: required DP band size",
            plot_width=resolution, plot_height=resolution,
            x_axis_label='relative size required', y_axis_label='relative amount'
        )

    #plot.vbar(x=buckets_x, bottom=0, top=data[1], width=w, color=greys[0], legend=value("inaccurate"))
    #plot.vbar(x=buckets_x, bottom=0, top=data[0], width=w, color="black", legend=value("accurate"))
    for i in range(num_buckets):
        l = buckets_x[i]
        r = l + w
        if data[0][i] < data[1][i]:
            plot.quad(left=l, bottom=-0.001, top=data[1][i], right=r, color=greys[0], line_width=0, legend=value("inaccurate"))
        plot.quad(left=l, bottom=-0.001, top=data[0][i], right=r, fill_color="black", line_width=0, legend=value("accurate"))
        if data[0][i] > data[1][i]:
            plot.quad(left=l, bottom=-0.001, top=data[1][i], right=r, color=greys[0], line_width=0, legend=value("inaccurate"))


    plot.title.text_font=font
    plot.legend.label_text_font=font
    plot.legend.label_text_baseline="bottom"
    plot.axis.axis_label_text_font=font
    plot.axis.axis_label_text_baseline="bottom"
    plot.axis.major_label_text_font=font
    plot.xaxis.major_label_standoff = 10
    plot.xgrid.grid_line_color = None
    plot.toolbar.logo = None
    plot.toolbar_location = None
    plot.legend.location = "top_center"
    save(plot, "nmwBandSize")

# actually call the functions that create the pictures

#unrelated_non_enclosed_seeds()
ambiguity_per_length()
#theoretical_max_acc()
#seed_shadows()
#alignment()
#stripOfConsideration()
#optimal_matching()
#required_nmw_band_size()

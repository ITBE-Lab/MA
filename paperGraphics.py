from bokeh.plotting import figure, output_file, show
from bokeh.layouts import row, column
from bokeh.models import Arrow, OpenHead, NormalHead, VeeHead
from bokeh.palettes import d3
from bokeh.io import export_png, export_svgs
from bokeh.models import FuncTickFormatter, FixedTicker, Label, ColorBar, FactorRange
from bokeh.models import LinearAxis, Range1d, LogColorMapper, FixedTicker, LinearColorMapper
from bokeh.models import ColumnDataSource
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

def ambiguity_per_length():
    data = [[8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21], [5150, 1165, 263, 46, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0], [9993, 2136, 456, 104, 24, 5, 1, 0, 0, 0, 0, 0, 0, 0], [19405, 4381, 993, 228, 53, 12, 3, 0, 0, 0, 0, 0, 0, 0], [78004, 17143, 3642, 728, 140, 29, 7, 2, 0, 0, 0, 0, 0, 0], [128484, 31693, 7638, 1840, 439, 104, 24, 6, 1, 0, 0, 0, 0, 0], [252549, 66284, 17353, 4411, 1122, 284, 72, 19, 5, 2, 1, 0, 0, 0], [9156386, 7272230, 5947033, 1072039, 485588, 243296, 132622, 74758, 16771, 25234, 1798, 413, 306, 49]]
    desc1 = ['index', 'min', 'quantile5%', 'quantile25%', 'median', 'quantile75%', 'quantile95%', 'max']

    """
    ref_len = 3000000000 # three billion => human genome length

    expected_ambiguity = []
    q_lens = list(range(10,21,1))
    for q_len in q_lens:
        expected_ambiguity.append( (0.25**q_len) * (ref_len - q_len) )

    plot = figure(title="Figure X: ambiguity for query length",
            x_axis_label='query length', y_axis_label='expected ambiguity',
            y_range=(-1,101)
        )

    plot.line(q_lens, expected_ambiguity, color="black")
    for i in range(len(data[0])):
        if data[0][i] < 10:
            data[0][i] = 10
        if data[0][i] > 20:
            data[0][i] = 20
    """
    for ii in range(1, len(data)):
        for i in range(len(data[ii])):
            if data[ii][i] > 1000:
                data[ii][i] = 1000

    plot = figure(
            title="Figure X: ambiguity on the human genome",
            plot_width=resolution, plot_height=resolution,
            y_range=(-1,101), x_range=(9.9, 20.1),
            x_axis_label='sequence length', y_axis_label='ambiguity'
        )

    plot.patch(list(reversed(data[0])) + data[0], list(reversed(data[7])) + data[1], legend="0-100%", color=greys[4])
    plot.patch(data[0] + list(reversed(data[0])), data[2] + list(reversed(data[6])), legend="5-95%", color=greys[2])
    plot.patch(data[0] + list(reversed(data[0])), data[3] + list(reversed(data[5])), legend="25-75%", color=greys[0])
    plot.line(data[0], data[4], legend="median", color="black")

    plot.title.text_font=font
    plot.legend.label_text_font=font
    plot.legend.label_text_baseline="bottom"
    plot.axis.axis_label_text_font=font
    plot.axis.axis_label_text_baseline="bottom"
    plot.axis.major_label_text_font=font
    plot.xaxis.major_label_standoff = 10
    plot.legend.location = "top_right"
    save(plot, "ambiguityPerQueryLen")

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
    data = ([0.29852045256745025, 0.0008703220191470844, 0.0017406440382941688, 0.0, 0.0, 0.0008703220191470844, 0.0, 0.0, 0.0008703220191470844, 0.0008703220191470844, 0.0017406440382941688, 0.0008703220191470844, 0.0, 0.0017406440382941688, 0.0008703220191470844, 0.0, 0.0, 0.0, 0.0, 0.0008703220191470844, 0.0, 0.0, 0.0, 0.0, 0.0008703220191470844, 0.0008703220191470844, 0.0, 0.0, 0.0, 0.0, 0.0017406440382941688, 0.0, 0.0, 0.0008703220191470844, 0.0008703220191470844, 0.0026109660574412533, 0.0, 0.0008703220191470844, 0.0017406440382941688, 0.0, 0.0017406440382941688, 0.0026109660574412533, 0.0026109660574412533, 0.0026109660574412533, 0.005221932114882507, 0.00783289817232376, 0.014795474325500439, 0.010443864229765015, 0.01566579634464752, 0.616187989556137], [0.04914390275780972, 0.007580708404130186, 0.009410534570644362, 0.010586851391974902, 0.011763168213305442, 0.011371062606195262, 0.013592994379819616, 0.016337733629590886, 0.018298261665141813, 0.018167559796105084, 0.011763168213305442, 0.018298261665141813, 0.02862370931904336, 0.018167559796105084, 0.022480721474317122, 0.018690367272251998, 0.03829564762776126, 0.0177754541889949, 0.021565808391060023, 0.0283623055809699, 0.02104300091491311, 0.021435106522023295, 0.02091229904587638, 0.018036857927068356, 0.0577702261142338, 0.015161416808260336, 0.014115801855966522, 0.019997385962619282, 0.016337733629590886, 0.022350019605280394, 0.016468435498627615, 0.011763168213305442, 0.005750882237616, 0.03450529342569614, 0.010194745784864722, 0.010848255130048355, 0.009541236439681088, 0.022219317736243666, 0.008887727094497455, 0.012024571951378895, 0.008495621487387275, 0.011109658868121809, 0.011109658868121809, 0.00522807476146909, 0.010717553261011628, 0.009018428963534181, 0.0039210560711018146, 0.0048359691543589075, 0.002221931773624362, 0.17370278394980673])

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
#ambiguity_per_length()
#theoretical_max_acc()
#seed_shadows()
#alignment()
#stripOfConsideration()
#optimal_matching()
required_nmw_band_size()

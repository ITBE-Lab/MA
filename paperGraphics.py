from bokeh.plotting import figure, output_file, show
from bokeh.layouts import row, column, gridplot
from bokeh.models import Arrow, OpenHead, NormalHead, VeeHead, AdaptiveTicker, FixedTicker
from bokeh.palettes import d3
from bokeh.io import export_png, export_svgs
from bokeh.models import FuncTickFormatter, FixedTicker, Label, ColorBar, FactorRange
from bokeh.models import LinearAxis, Range1d, LogColorMapper, FixedTicker, LinearColorMapper
from bokeh.models import ColumnDataSource, CompositeTicker, SingleIntervalTicker, BasicTickFormatter
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

green = "#77933C"
blue = "#376092"
red = "#953735"
purple = "#604A7B"
orange = "#E46C0A"

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

def simulate_max_length(q_len, mutation_amount, indel_amount,
                        indel_size, sim_amount, max_missmatches
                       ):
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

        #get the starting positions for every match
        starts = []
        last = -1
        s_ind = 0
        for ind, num in enumerate(q):
            if num != last + 1:
                starts.append(s_ind)
                s_ind = ind
            last = num
        starts.append(s_ind)

        # get the results
        matches = []
        for start in starts:
            length = 0
            num_mm = 0
            while start+length+1 < len(q) and num_mm <= max_missmatches:
                if q[start] + length + 1 != q[start+length+1]:
                    num_mm += 1
                length += 1
            matches.append(length)

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

def genome_on_genome():
    return ( )


def random_on_random():
    return (  )



def ambiguity_per_length():
    high = 99
    low = 1
    data1_, data2_ = genome_on_genome()

    data1 = []
    for row in data1_:
        data1.append( [] )
        for ele in row:
            data1[-1].append(ele * (high-low) + low)
    data2 = []
    for row in data2_:
        data2.append( [] )
        for ele in row:
            data2[-1].append(ele * (high-low) + low)

    r1max = 10
    r2size = 15
    min_len=10
    max_len=20


    color_mapper = LogColorMapper(
                    palette=heatmap_palette(light_spec_approximation, 256),
                    low=low,
                    high=high
                )

    plot = figure(title="ambiguity on human genome",
            x_range=(0,r1max), y_range=(min_len, max_len),
            x_axis_label='ambiguity', y_axis_label='sequence length',
            plot_width=resolution, plot_height=resolution,
            min_border_bottom=10, min_border_top=10,
            min_border_left=10, min_border_right=15,
            x_axis_type=None, y_axis_type=None
        )
    for index, row in enumerate(data1):
        plot.image(image=[[row]], color_mapper=color_mapper,
                dh=[.6], dw=[r1max], x=[0], y=[min_len + index + 0.2])

    plot2 = figure(x_range=(r1max,2**r2size+r1max), y_range=(min_len, max_len),
            min_border_bottom=10, min_border_top=10,
            min_border_left=20, min_border_right=15,
            plot_width=resolution*3/4, plot_height=resolution,tools=[],
            x_axis_type="log"
        )
    for index, row in enumerate(data2):
        plot2.image(image=[[row]], color_mapper=color_mapper,
            dh=[.6], dw=[2**r2size+r1max], x=[r1max], y=[min_len + index + 0.2])

    ticks = []
    num_ticks = 6
    for tick in range(num_ticks):
        ticks.append( math.exp( (tick/float(num_ticks-1)) * math.log(high)) )#(*(high-low)) / float(math.exp(num_ticks)) + low )

    print(ticks)

    size = "12pt"
    #color_bar = ColorBar(color_mapper=color_mapper, border_line_color=None, location=(0,0), ticker=FixedTicker(ticks=ticks))
    color_bar = ColorBar(color_mapper=color_mapper, border_line_color=None, location=(0,0))
    color_bar.major_label_text_font=font
    color_bar.major_label_text_font_size=size
    #color_bar.label_standoff = 20

    ticker = SingleIntervalTicker(interval=1, num_minor_ticks=0)
    plot.add_layout(LinearAxis(ticker=ticker), 'left')
    plot.add_layout(LinearAxis(ticker=ticker), 'below')
    plot.add_layout(color_bar, 'left')
    plot.legend.label_text_font=font
    plot.background_fill_color = dark_greys[2]
    plot.background_fill_alpha = 1
    plot.axis.axis_label_text_font=font
    plot.axis.major_label_text_font=font
    plot.axis.major_label_text_font_size=size
    plot.xaxis.major_label_standoff = 15
    plot.xaxis.minor_tick_line_color = None
    plot.yaxis.minor_tick_line_color = None
    plot.grid.grid_line_color = None

    plot2.legend.label_text_font=font
    plot2.yaxis.visible = False
    plot2.background_fill_color = dark_greys[2]
    plot2.background_fill_alpha = 1
    plot2.axis.axis_label_text_font=font
    plot2.axis.major_label_text_font_size=size
    plot2.axis.major_label_text_font=font
    plot2.xaxis.major_label_standoff = 15
    plot2.grid.grid_line_color = None

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
    indel_size = 100
    prob_refindable = []
    quality = 20
    depth = 16
    max_missmatches = 0
    check_for_min_size_instead = None

    print("creating query length matrix...")
    max_indels = int(q_len/indel_size)*2
    for num_mut in range(0, q_len, max(10,int(q_len/quality))):
        prob_refindable.append( [] )
        for num_indel in range(0,max_indels, max(2,int(max_indels/quality))):
            q_len_e = simulate_max_length(q_len, num_mut, num_indel, indel_size, depth, max_missmatches)
            if q_len_e is None:
                prob_refindable[-1].append(float('NaN'))
            else:
                probs = []
                for x in q_len_e:
                    if check_for_min_size_instead is None:
                        probs.append(prob_all_non_enclosed(x, ref_len))
                    else:
                        p = 0
                        for e in x:
                            if e >= check_for_min_size_instead:
                                p = 1
                        probs.append(p)
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
    save(plot, "upperBound")


#static stuff
I = "I"
M = "M"
MM = "\=M"
D = "D"
query =     [                         "T", "A", "C", "A", "T", "T", "C", "T" ]
reference = ["T", "T", "C", "A", "G",           "C", "A", "T", "A", "C", "T", "C", "A"]
l_alignment = [D,D,I,I,M,M,D,D,D,M,MM,M,M,D]

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
    seeds = [
        (1.5,1.5,2, orange),
        (5.5,2.5,2, green),
        (-.5,3.5,3, blue)
    ]

    for x, y, size, color in seeds:
        patch_x = []
        patch_y = []

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
                fill_color=color,
                fill_alpha=.3,
                line_color=None,
                #line_width=2,
                #line_dash=[2,2],
            )

    for x, y, size, color in seeds:
        seeds_x = []
        seeds_y = []
        seeds_x.append(x)
        seeds_x.append(x + size)
        seeds_x.append(float('nan'))
        seeds_y.append(y)
        seeds_y.append(y + size)
        seeds_y.append(float('nan'))

        plot.line(
                seeds_x,
                seeds_y,
                color=color,
                line_width=5
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
            color=green,
            line_width=5
        )

    plot.multi_line(
            mm_x,
            mm_y,
            legend="mm",
            color=orange,
            line_width=5
        )

    plot.multi_line(
            i_x,
            i_y,
            legend="ins_____",
            color=blue,
            line_width=5
        )

    plot.multi_line(
            d_x,
            d_y,
            legend="del",
            color=purple,
            line_width=5
        )

    plot.xaxis.ticker = FixedTicker(ticks=range(len(reference)))
    plot.legend.location = None
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
            fill_color=green,
            fill_alpha=.3,
            line_color=None,
            #line_width=2,
            #line_dash=[2,2],
        )


    plot.line(
        [-.5,7.5],
        [-.5,7.5],
        color=green,
        line_width=1,
        line_dash=[2,2]
    )

    plot.line(
        [5.5,12.5],
        [-.5,6.5],
        color=green,
        line_width=1,
        line_dash=[2,2]
    )

    plot.line(
        [4.5,7.5],
        [1.5,4.5],
        color=blue,
        line_width=3
    )
    plot.line(
        [2.5,4.5],
        [-.5,1.5],
        color=blue,
        line_alpha=.3,
        line_width=1,
        line_dash=[8,2]
    )


    plot.line(
        [10.5,11.5],
        [1.5,2.5],
        color=orange,
        line_width=3
    )
    plot.line(
        [8.5,10.5],
        [-.5,1.5],
        color=orange,
        line_alpha=.3,
        line_width=1,
        line_dash=[8,2]
    )

    plot.line(
        [1.5,2.5],
        [1.5,2.5],
        color=blue,
        line_width=3
    )
    plot.line(
        [-.5,1.5],
        [-.5,1.5],
        color=blue,
        line_alpha=.3,
        line_width=1,
        line_dash=[8,2]
    )

    plot.line(
        [-.5,2.5],
        [3.5,6.5],
        color=orange,
        line_width=3
    )
    plot.line(
        [-1.0,-.5],
        [3.0,3.5],
        color=orange,
        line_alpha=.3,
        line_width=1,
        line_dash=[8,2]
    )

    plot.line(
        [0.5,1.5],
        [-.5,0.5],
        color=blue,
        line_width=3
    )

    plot.line(
        [9.5,11.5],
        [4.5,6.5],
        color=blue,
        line_width=3
    )
    plot.line(
        [4.5,9.5],
        [-.5,4.5],
        color=blue,
        line_alpha=.3,
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


def forced_gap():
    plot = figure(
                #title="Figure X: Forced Gaps",
                plot_width=resolution, plot_height=resolution,
                x_range=[-1,13], y_range=[-1,8]
            )
    plot.line(
        [-1,9],
        [-1,9],
        color=blue,
        line_alpha=.3,
        line_width=1,
        line_dash=[8,2]
    )
    plot.line(
        [1.5,3.5],
        [1.5,3.5],
        color=blue,
        line_width=3
    )

    plot.line(
        [2,12],
        [-1,9],
        color=blue,
        line_alpha=.3,
        line_width=1,
        line_dash=[8,2]
    )
    plot.line(
        [8.5,10.5],
        [5.5,7.5],
        color=blue,
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
    save(plot, "forcedGap")


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
            fill_color=green,
            fill_alpha=.3,
            line_color=None,
            #line_width=2,
            #line_dash=[2,2],
        )
    plot.patch(
            [10.5,10.5,max_x,max_x],
            [5.5,max_y,max_y,5.5],
            fill_color=green,
            fill_alpha=.3,
            line_color=None,
            #line_width=2,
            #line_dash=[2,2],
        )
    plot.patch(
            [min_x,min_x,0.5,0.5],
            [min_y,-.5,-.5,min_y],
            fill_color=green,
            fill_alpha=.3,
            line_color=None,
            #line_width=2,
            #line_dash=[2,2],
        )

    plot.line(
        [0.5,1.5],
        [-.5,0.5],
        color=blue,
        line_width=3
    )
    plot.line(
        [4.5,7.5],
        [1.5,4.5],
        color=blue,
        line_width=3
    )
    plot.line(
        [9.5,10.5],
        [4.5,5.5],
        color=blue,
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
                plot_width=resolution, plot_height=resolution/2,
                y_range=["reference", "1", "2", "3", "query1", "query2", "query3", "query4"]
            )

    seeds = [
        ("query1", 0, [450], green, 7),
        ("query1", 7, [500], green, 13),
        ("query1", 24, [525], green, 13),
        ("query1", 47, [547], green, 24),

        ("query2", 5, [270], orange, 3),
        ("query2", 22, [150], blue, 5),
        ("query2", 36, [1200], blue, 20),
        
        ("query3", 6, [20,1340], orange, 2),
        ("query3", 19, [1500], blue, 4),
        ("query3", 35, [320, 1150, 1477], orange, 3),
        ("query3", 46, [1000], orange, 13),

        ("query4", 20, [30, 750], orange, 4),
        ("query4", 33, [800], orange, 4),
    ]

    max_q = 0
    max_r = 0
    for seed in seeds:
        if max_q < seed[1] + seed[-1]:
            max_q = seed[1] + seed[-1]
        if max_r < seed[2][-1] + seed[-1]:
            max_r = seed[2][-1] + seed[-1]

    r_fac = max_q/float(max_r)
    for seed in seeds:
        for end in seed[2]:
            plot.patch(
                [seed[1]+seed[-1], (end+seed[-1])*r_fac, end*r_fac, seed[1]],
                [seed[0], "reference", "reference", seed[0]],
                fill_color=seed[-2],
                fill_alpha=.3,
                line_color=None,
            )
    plot.line([0, max_q], ["reference", "reference"], color="black", line_width=3)
    for seed in seeds:
        plot.line([seed[1], seed[1] + seed[-1]], [seed[0], seed[0]], color=seed[-2], line_width=3)



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
#forced_gap()
ambiguity_per_length()
#theoretical_max_acc()
#seed_shadows()
#alignment()
#stripOfConsideration()
#optimal_matching()
#required_nmw_band_size()

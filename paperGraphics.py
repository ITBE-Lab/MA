from bokeh.plotting import figure, output_file, show
from bokeh.layouts import row, column
from bokeh.models import Arrow, OpenHead, NormalHead, VeeHead
from bokeh.palettes import d3
from bokeh.io import export_png, export_svgs
from bokeh.models import FuncTickFormatter, FixedTicker, Label, ColorBar, FactorRange
from bokeh.models import LinearAxis, Range1d, LogColorMapper, FixedTicker, LinearColorMapper
import math
import random
import numpy as np


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

    plot.title.text_font=font
    plot.legend.label_text_font=font
    plot.legend.label_text_baseline="bottom"
    plot.axis.axis_label_text_font=font
    plot.axis.major_label_text_font=font
    plot.legend.location = "top_left"
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

def seed_shadows():
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

def alignment():
    I = "I"
    M = "M"
    MM = "\=M"
    D = "D"


    query =     [                         "C", "A", "C", "A", "T", "A", "T", "T" ]
    reference = ["A", "G", "G", "A", "G", "C", "A",           "T", "T", "T", "T", "C", "A"]

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

def stripOfConsideration():
    plot = figure(
                title="Figure 3: Strip of consideration",
                plot_width=resolution, plot_height=resolution,
                x_axis_label = "reference", y_axis_label = "query"
            )
    
    s_gap = 16.0
    s_extend = 1.0
    s_match = 8.0
    s_missmatch = 2.0
    anchor = ((3,5),(106,8))

    q_len = 1000

    bottom_left_x = anchor[0][0] + s_gap / s_extend - (1 + s_match/s_extend) * (anchor[0][0]-1)
    top_left_x = anchor[0][0] + anchor[1][0] + s_gap / s_extend - (1 + s_match/s_extend) * q_len

    bottom_right_x = anchor[0][1] + anchor[1][1] - s_gap / s_extend + (s_match / s_extend) * q_len + (2 + s_match / s_extend) * anchor[0][0] - (1 + s_match/s_extend) * (anchor[0][0]-1)
    top_right_x = anchor[0][1] + anchor[1][1] - s_gap / s_extend + (s_match / s_extend) * q_len + (2 + s_match / s_extend) * anchor[0][0] - (1 + s_match/s_extend) * q_len

    print(bottom_left_x, bottom_right_x, top_left_x, top_right_x)
    return

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

def unrelated_non_enclosed_seeds():
    plot = figure(
                title="Supplementary Figure X: unrelated Seeds",
                plot_width=resolution, plot_height=resolution/2,
                x_axis_label = "query",
                y_range=["related seeds", "unrelated seeds 1", "unrelated seeds 2"]
            )

    related = ["related seeds", "related seeds"]
    unrelated = ["unrelated seeds 1", "unrelated seeds 1"]
    unrelated2 = ["unrelated seeds 2", "unrelated seeds 2"]

    plot.line([-5,3], related, color="black", line_width=10)
    plot.line([6,21], related, color="black", line_width=10)
    plot.line([25,40], related, color="black", line_width=10)

    plot.line([2,7], unrelated, color="black", line_width=10)
    plot.line([20,27], unrelated, color="black", line_width=10)

    plot.line([1,4], unrelated2, color="black", line_width=10)
    plot.line([5,8], unrelated2, color="black", line_width=10)
    plot.line([19,22], unrelated2, color="black", line_width=10)
    plot.line([23,29], unrelated2, color="black", line_width=10)


    plot.title.text_font=font
    plot.legend.label_text_font=font
    plot.legend.label_text_baseline="hanging"
    plot.axis.axis_label_text_font=font
    plot.axis.major_label_text_font=font
    plot.xaxis.major_tick_line_color = None
    plot.yaxis.major_tick_line_color = None
    plot.xaxis.minor_tick_line_color = None
    plot.xaxis.major_label_text_alpha = 0
    plot.toolbar.logo = None
    plot.toolbar_location = None
    plot.ygrid.grid_line_color = None
    plot.xgrid.ticker = FixedTicker(ticks=[-5,3,6,21,25,40])
    save(plot, "unrelatedNonEnclosedSeeds")


# actually call the functions that create the pictures

#unrelated_non_enclosed_seeds()
#ambiguity_per_length()
#theoretical_max_acc()
#seed_shadows()
#alignment()
stripOfConsideration()

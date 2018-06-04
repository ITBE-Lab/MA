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
import json
import os.path

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
def format(rgb):
    def clamp(x):
        return max(0, min(x, 255))
    red, green, blue = rgb
    return "#{0:02x}{1:02x}{2:02x}".format(clamp(int(red * 255)), clamp(int(green * 255)),
                                            clamp(int(blue * 255)))

def heatmap_palette(scheme, num_colors):
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
        show(plot)
        export_png(plot, filename="paperGraphics/" + name + ".png")
        plot.output_backend = "svg"
        export_svgs(plot, filename="paperGraphics/" + name + ".svg")

resolution = 300
min_x = 0
min_y = 0
max_x = 10
max_y = 10

def human():
    return (10, 20, [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2e-05, 0.0], [0.0, 4e-05, 6e-05, 0.00013, 0.00013, 0.00016, 0.00016, 0.0002, 0.00019, 0.00021], [0.0, 0.00041, 0.0009, 0.00102, 0.00142, 0.00146, 0.00207, 0.00257, 0.00242, 0.00264], [0.0, 0.0047, 0.00821, 0.00913, 0.00887, 0.0084, 0.00824, 0.00712, 0.00691, 0.00657], [0.0, 0.02985, 0.03135, 0.03083, 0.03405, 0.03523, 0.03673, 0.03636, 0.03484, 0.03492], [0.0, 0.12077, 0.11714, 0.10139, 0.08106, 0.06348, 0.04928, 0.03868, 0.03112, 0.02517], [0.0, 0.32592, 0.18053, 0.09592, 0.05587, 0.03392, 0.02353, 0.01728, 0.01269, 0.00984], [0.0, 0.52612, 0.14507, 0.05049, 0.02688, 0.01548, 0.01089, 0.0086, 0.00692, 0.00544], [0.0, 0.64777, 0.09116, 0.0288, 0.01576, 0.0105, 0.00727, 0.00595, 0.00476, 0.00399]], [[0.0, 0.0, 0.0, 0.0, 1e-05, 0.00013, 0.00062, 0.00134, 0.00872, 0.02704, 0.02922, 0.07054, 0.29214, 0.33081, 0.15811], [0.0, 2e-05, 8e-05, 0.00029, 0.00061, 0.00202, 0.00909, 0.02523, 0.02721, 0.06647, 0.24793, 0.30554, 0.16821, 0.07773, 0.03757], [0.00015, 0.0004, 0.00114, 0.00323, 0.01093, 0.02215, 0.02588, 0.06355, 0.21827, 0.27593, 0.17325, 0.08584, 0.04475, 0.03016, 0.01967], [0.00284, 0.00618, 0.01069, 0.01539, 0.0237, 0.07261, 0.19663, 0.24542, 0.1636, 0.08782, 0.0481, 0.03107, 0.02395, 0.02122, 0.016], [0.00676, 0.01387, 0.03277, 0.08466, 0.16403, 0.19931, 0.1494, 0.08671, 0.04892, 0.03162, 0.0237, 0.02054, 0.01909, 0.01927, 0.01375], [0.03148, 0.05795, 0.09006, 0.11482, 0.10302, 0.07422, 0.04895, 0.03154, 0.02375, 0.02041, 0.01817, 0.01738, 0.01831, 0.01725, 0.0134], [0.02068, 0.03121, 0.04063, 0.04282, 0.03597, 0.02952, 0.02392, 0.01999, 0.01822, 0.01759, 0.01665, 0.01735, 0.01739, 0.01574, 0.01145], [0.0084, 0.0126, 0.01669, 0.01903, 0.01885, 0.01828, 0.01667, 0.0165, 0.0159, 0.01539, 0.01633, 0.01615, 0.01703, 0.01495, 0.01037], [0.00491, 0.00742, 0.01034, 0.01258, 0.01365, 0.01535, 0.01528, 0.0146, 0.01457, 0.01423, 0.01462, 0.01599, 0.01639, 0.01465, 0.00933], [0.00341, 0.00576, 0.00831, 0.0113, 0.01223, 0.01367, 0.01372, 0.01361, 0.01357, 0.01357, 0.01402, 0.0151, 0.01508, 0.01323, 0.00835]])

def random():
    return (10, 20, [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 1e-05, 1e-05, 4e-05, 0.00018, 0.00028], [0.0, 0.00389, 0.02021, 0.05774, 0.10963, 0.15351, 0.17115, 0.15822, 0.12546, 0.08887], [0.0, 0.24808, 0.34334, 0.24043, 0.11374, 0.03933, 0.01143, 0.00292, 0.00062, 0.00011], [0.0, 0.70461, 0.24673, 0.04322, 0.00507, 0.00032, 5e-05, 0.0, 0.0, 0.0], [0.0, 0.9169, 0.07942, 0.0036, 7e-05, 1e-05, 0.0, 0.0, 0.0, 0.0], [0.0, 0.97892, 0.02081, 0.00027, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]], [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.02647, 0.97352, 1e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.00086, 0.0051, 0.06046, 0.54301, 0.38955, 0.0005, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.05331, 0.04551, 0.01218, 0.00032, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]])

def plasmodium():
    return (7, 17, [[0.0, 0.0, 0.0006103515625, 0.0003662109375, 0.000732421875, 0.000732421875, 0.0003662109375, 0.0006103515625, 0.0015869140625, 0.0008544921875], [0.0, 0.008062044019373429, 0.009441481209000062, 0.01149530991355527, 0.011265403715284164, 0.012231009748022806, 0.012660167984795537, 0.010545030960701368, 0.011571945312978971, 0.011403347434246827], [0.0, 0.00237, 0.00477, 0.00633, 0.00695, 0.00776, 0.00872, 0.00823, 0.0086, 0.00911], [0.0, 0.01046, 0.01326, 0.01504, 0.01476, 0.01487, 0.01464, 0.01473, 0.01415, 0.01367], [0.0, 0.03502, 0.03565, 0.03222, 0.0287, 0.02671, 0.02369, 0.02135, 0.02025, 0.01853], [0.0, 0.09254, 0.06767, 0.05184, 0.04241, 0.03619, 0.03063, 0.02605, 0.02457, 0.02125], [0.0, 0.18386, 0.10034, 0.06672, 0.05128, 0.03853, 0.03196, 0.0284, 0.02445, 0.02134], [0.0, 0.29991, 0.12096, 0.07134, 0.04915, 0.0377, 0.0297, 0.02453, 0.02142, 0.01879], [0.0, 0.42895, 0.12538, 0.06683, 0.04359, 0.03125, 0.02397, 0.02098, 0.01858, 0.01508], [0.0, 0.54107, 0.11882, 0.05623, 0.03647, 0.02437, 0.01874, 0.01554, 0.01417, 0.01175]], [[0.0018310546875, 0.0035400390625, 0.0084228515625, 0.0162353515625, 0.0355224609375, 0.07470703125, 0.1068115234375, 0.1400146484375, 0.162841796875, 0.157470703125, 0.1221923828125, 0.0838623046875, 0.04150390625, 0.0263671875, 0.0059814453125], [0.010682974679664031, 0.02131996811967384, 0.0383176997118509, 0.06792961804916928, 0.09824658206118571, 0.12735270676230764, 0.14323156152289865, 0.13990558518791, 0.10286003310649255, 0.07341671264790632, 0.041153209490527864, 0.02090613696278585, 0.010713628839433512, 0.002651584820060082, 0.0016553246275519589], [0.00875, 0.01758, 0.03488, 0.06692, 0.11473, 0.16784, 0.18784, 0.14941, 0.09773, 0.04868, 0.02592, 0.01094, 0.00316, 0.00196, 0.00066], [0.01285, 0.02457, 0.0455, 0.07602, 0.11725, 0.14881, 0.16104, 0.13525, 0.08107, 0.04459, 0.01717, 0.00633, 0.00265, 0.00088, 0.00038], [0.01797, 0.03295, 0.05456, 0.08493, 0.11277, 0.13122, 0.12098, 0.0963, 0.06192, 0.0263, 0.01172, 0.00414, 0.00124, 0.0007, 0.00014], [0.01908, 0.03476, 0.05541, 0.0833, 0.09977, 0.09939, 0.08932, 0.06301, 0.03414, 0.01846, 0.00683, 0.0018, 0.00104, 0.0004, 0.0001], [0.01893, 0.03375, 0.05126, 0.06676, 0.07547, 0.07247, 0.0562, 0.03819, 0.02293, 0.01115, 0.00371, 0.00108, 0.0008, 0.00036, 2e-05], [0.01617, 0.0277, 0.04238, 0.05228, 0.05393, 0.04795, 0.03649, 0.02387, 0.01463, 0.00696, 0.00216, 0.00098, 0.00072, 0.00022, 2e-05], [0.01233, 0.02171, 0.02982, 0.03606, 0.03807, 0.03034, 0.02317, 0.01556, 0.01036, 0.0048, 0.00147, 0.00096, 0.00056, 0.00012, 2e-05], [0.00908, 0.01562, 0.0218, 0.0273, 0.02691, 0.01995, 0.01675, 0.01143, 0.00802, 0.00321, 0.00119, 0.00102, 0.00038, 0.00012, 2e-05]])

def mouse():
    return (10, 20, [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2e-05, 2e-05, 1e-05], [0.0, 1e-05, 3e-05, 7e-05, 6e-05, 0.00017, 0.00022, 0.00017, 0.00015, 0.00021], [0.0, 0.00055, 0.00121, 0.00128, 0.00138, 0.002, 0.00216, 0.00253, 0.00261, 0.00289], [0.0, 0.0063, 0.00924, 0.00992, 0.00984, 0.00835, 0.0078, 0.0076, 0.00716, 0.00701], [0.0, 0.0338, 0.03425, 0.03441, 0.03691, 0.03902, 0.04283, 0.04177, 0.04091, 0.03957], [0.0, 0.1351, 0.1303, 0.1104, 0.08666, 0.06607, 0.04834, 0.03531, 0.02882, 0.02148], [0.0, 0.3585, 0.18855, 0.09019, 0.04575, 0.0264, 0.01824, 0.01152, 0.0088, 0.00738], [0.0, 0.55785, 0.13091, 0.04064, 0.01994, 0.01174, 0.00845, 0.00627, 0.0051, 0.00419], [0.0, 0.66305, 0.07425, 0.02247, 0.01204, 0.00764, 0.00629, 0.00484, 0.00419, 0.0039]], [[0.0, 0.0, 0.0, 0.0, 2e-05, 0.00022, 0.00083, 0.00167, 0.01034, 0.02986, 0.02595, 0.06818, 0.34316, 0.33024, 0.12399], [1e-05, 4e-05, 0.00012, 0.00027, 0.00079, 0.00222, 0.01057, 0.02693, 0.02522, 0.06636, 0.29185, 0.32358, 0.13109, 0.05386, 0.04064], [0.00024, 0.00051, 0.0012, 0.00354, 0.01176, 0.0239, 0.02401, 0.06896, 0.25392, 0.30159, 0.14742, 0.05858, 0.03209, 0.02511, 0.02775], [0.00332, 0.00538, 0.01112, 0.01667, 0.0241, 0.07884, 0.22551, 0.25921, 0.15153, 0.07004, 0.0359, 0.02286, 0.01997, 0.02019, 0.02314], [0.00725, 0.01482, 0.03885, 0.09458, 0.18756, 0.20175, 0.12898, 0.07399, 0.04611, 0.02768, 0.01945, 0.01712, 0.01542, 0.01876, 0.02141], [0.03578, 0.06267, 0.09454, 0.10615, 0.08358, 0.05625, 0.0398, 0.03797, 0.02741, 0.01821, 0.01666, 0.0149, 0.01406, 0.01754, 0.02025], [0.01735, 0.02518, 0.03106, 0.03105, 0.02878, 0.0245, 0.02703, 0.03039, 0.02115, 0.0165, 0.01445, 0.01307, 0.0128, 0.01714, 0.01778], [0.00606, 0.00923, 0.0139, 0.0166, 0.01798, 0.01807, 0.02367, 0.02553, 0.01915, 0.01411, 0.01392, 0.01213, 0.01202, 0.01744, 0.01676], [0.00402, 0.00647, 0.00991, 0.01313, 0.01481, 0.0163, 0.0207, 0.02389, 0.01647, 0.0133, 0.01227, 0.01238, 0.01183, 0.01735, 0.01526], [0.00327, 0.00537, 0.0089, 0.01102, 0.01357, 0.01477, 0.02035, 0.02339, 0.01559, 0.0122, 0.01216, 0.01114, 0.01196, 0.01764, 0.01387]])

def ambiguity_per_length():
    high = 99
    low = 1
    min_len, max_len, data1_, data2_ = plasmodium()

    data1 = []
    for row in data1_:
        data1.append( [] )
        for ele in row[1:]:
            data1[-1].append(ele * (high-low) + low)
    data2 = []
    for row in data2_:
        data2.append( [] )
        for ele in row:
            data2[-1].append(ele * (high-low) + low)

    r1max = 10
    r2size = 15

    color_mapper = LogColorMapper(
                    palette=heatmap_palette(light_spec_approximation, 256),
                    low=low,
                    high=high
                )

    plot = figure(title="ambiguity on human genome",
            x_range=(1,r1max), y_range=(min_len, max_len),
            x_axis_label='ambiguity', y_axis_label='sequence length',
            plot_width=resolution, plot_height=resolution,
            min_border_bottom=10, min_border_top=10,
            min_border_left=10, min_border_right=15,
            x_axis_type=None, y_axis_type=None
        )
    for index, row in enumerate(data1):
        plot.image(image=[[row]], color_mapper=color_mapper,
                dh=[.6], dw=[r1max-1], x=[1], y=[min_len + index + 0.2])

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
    query = ["G", "C", "A", "C", "A", "G", "T" ]
    reference = ["G", "C", "G", "A", "T", "A", "G", "A", "C", "G", "T" ]
    min_x = -1.0
    min_y = -1.0
    max_x = 13
    max_y = 8
    plot = figure(
                plot_width=resolution, plot_height=resolution*len(query)/len(reference),
                x_range=[-1,11], y_range=[-1,7]
            )
    # x y size draw_shadow?
    seeds = [
        #(1.5,1.5,2, orange),
        #(5.5,2.5,2, green),
        #(-.5,3.5,3, blue)
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

    #plot.line(
    #        [min_x, seeds[0][0], seeds[0][0]+1, seeds[1][0], seeds[1][0]+5.5],
    #        [min_y, seeds[0][1], seeds[0][1]+1, seeds[1][1], max_y],
    #        color="black",
    #        line_dash=[2,2],
    #        line_width=2
    #    )


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
    #plot.legend.label_text_baseline="hanging"
    plot.axis.axis_label_text_font=font
    plot.axis.major_label_text_font=font
    save(plot, "shadows")

def alignment():
    plot = figure(
                plot_width=resolution*13/8, plot_height=resolution
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
            legend="match",
            color=green,
            line_width=5
        )

    plot.multi_line(
            mm_x,
            mm_y,
            legend="missmatch",
            color=orange,
            line_width=5
        )

    plot.multi_line(
            i_x,
            i_y,
            legend="insertion",
            color=blue,
            line_width=5
        )

    plot.multi_line(
            d_x,
            d_y,
            legend="deletion",
            color=purple,
            line_width=5
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
    plot.legend.label_text_font_size="11pt"
    #plot.legend.label_text_baseline="bottom"
    plot.axis.major_label_text_font=font
    plot.axis.major_label_text_font_size="11pt"
    print("saving")
    save(plot, "alignment")

def stripOfConsideration():
    query =     [                         "A", "C", "A", "T", "T", "C" ]
    reference = ["T", "T", "C", "A", "G",      "C", "A", "T", "C", "A", "T", "C"]
    plot = figure(
                plot_width=resolution, plot_height=resolution*6/12,
                x_range=[-1,12], y_range=[-1,6]
            )
    plot.patch(
            [-.5,5.5,11.5,5.5],
            [-.5,5.5,5.5,-.5],
            fill_color=green,
            fill_alpha=.3,
            line_color=None,
            #line_width=2,
            #line_dash=[2,2],
        )


    plot.line(
        [-.5,5.5],
        [-.5,5.5],
        color=green,
        line_width=1,
        line_dash=[2,2]
    )

    plot.line(
        [5.5,11.5],
        [-.5,5.5],
        color=green,
        line_width=1,
        line_dash=[2,2]
    )

    plot.line(
        [4.5,7.5],
        [0.5,3.5],
        color=blue,
        line_width=3
    )
    plot.line(
        [3.5,4.5],
        [-.5,0.5],
        color=blue,
        line_alpha=.3,
        line_width=1,
        line_dash=[8,2]
    )


    plot.line(
        [7.5,9.5],
        [0.5,2.5],
        color=orange,
        line_width=3
    )
    plot.line(
        [6.5,7.5],
        [-0.5,0.5],
        color=orange,
        line_alpha=.3,
        line_width=1,
        line_dash=[8,2]
    )

    plot.line(
        [1.5,3.5],
        [0.5,2.5],
        color=blue,
        line_width=3
    )
    plot.line(
        [0.5,1.5],
        [-.5,0.5],
        color=blue,
        line_alpha=.3,
        line_width=1,
        line_dash=[8,2]
    )

    plot.line(
        [-.5,2.5],
        [2.5,5.5],
        color=orange,
        line_width=3
    )
    plot.line(
        [-1.0,-.5],
        [2.0,2.5],
        color=orange,
        line_alpha=.3,
        line_width=1,
        line_dash=[8,2]
    )

    """
    plot.line(
        [0.5,1.5],
        [-.5,0.5],
        color=blue,
        line_width=3
    )
    """

    plot.line(
        [9.5,11.5],
        [3.5,5.5],
        color=blue,
        line_width=3
    )
    plot.line(
        [5.5,9.5],
        [-.5,3.5],
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
    #plot.legend.label_text_baseline="hanging"
    plot.axis.axis_label_text_font=font
    plot.axis.major_label_text_font=font
    plot.axis.major_label_text_font_size="11pt"
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

def contradicting_seeds():
    min_x = -1.0
    min_y = -1.0
    max_x = 13
    max_y = 8
    reference = ["G", "C", "G", "A", "T", "A", "T", "A", "C", "G", "G"]
    query = ["G", "C", "A", "T", "T", "G", "G"]
    plot = figure(
                plot_width=resolution, plot_height=resolution*len(query)/len(reference)
            )

    plot.line(
        [-.5,1.5],
        [-.5,1.5],
        color=blue,
        line_width=3
    )
    plot.line(
        [2.5,4.5],
        [1.5,3.5],
        color=blue,
        line_width=3
    )
    plot.line(
        [4.5,6.5],
        [1.5,3.5],
        color=blue,
        line_width=3
    )
    plot.line(
        [8.5,10.5],
        [4.5,6.5],
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
    save(plot, "contradicting_seeds")

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

nan = float("nan")


def accuracy_pics():
    show_coverage = True
    #plots = [ [], [], [] ]
    plots = [ [] ]

    def _decode(o):
        if isinstance(o, str) or isinstance(o, unicode):
            try:
                return float(o)
            except ValueError:
                return o
        elif isinstance(o, dict):
            return {_decode(k): _decode(v) for k, v in o.items()}
        elif isinstance(o, list):
            return [_decode(v) for v in o]
        else:
            return o

    def makePicFromDict(d, title, print_out=False,desc2=None, set_max=None, inner=None):
        x_hover = []
        y_hover = []
        c_inner = []
        desc1_hover = []
        desc2_hover = []
        desc3_hover = []
        desc4_hover = []
        c_hover = []
        min_ = 0.0
        max_ = 1.0
        if set_max != None:
            max_ = set_max
        for x, value in d.items():
            for y, value in value.items():
                max_ = max(max_, value)
                min_ = min(min_, value)

        colors = heatmap_palette(light_spec_approximation, 127)
        w = 0
        h = 0
        w_keys = sorted(d.keys())
        width = 0
        height = 0
        if len(w_keys) > 1:
            width = w_keys[1] - w_keys[0]
            h_keys = sorted(d[w_keys[0]].keys())
            height = h_keys[1] - h_keys[0]
            for x, value in d.items():
                for y, value in value.items():
                    x_hover.append(x + width/2)
                    y_hover.append(y + height/2)
                    desc1_hover.append(str(x))
                    desc3_hover.append(str(y))
                    if inner != None:
                        if inner[x][y] == None or max_ == 0:
                            c_inner.append("#AAAAAA")
                            desc2_hover.append("")
                        else:
                            c_inner.append( colors[int(126*(inner[x][y]-min_) /max_)] )
                            desc2_hover.append(str(inner[x][y]))
                    else:
                        desc2_hover.append(str(d[x][y]))
                    if desc2 != None:
                        desc4_hover.append(desc2[x][y])
                    if d[x][y] == None or max_ == 0:
                        c_hover.append("#AAAAAA")
                    else:
                        c_hover.append( colors[int(126*(d[x][y]-min_) /max_)] )

        color_mapper = LinearColorMapper(
                palette=heatmap_palette(light_spec_approximation, 127),
                low=min_,
                high=max_
            )

        tick_formater = FuncTickFormatter(code="""
            return Math.max(Math.floor( (tick+1)/2),0) + '; ' +
                    Math.max(Math.floor( (tick)/2),0)"""
            )
        #tick_formater = FuncTickFormatter(code="return 'a')

        source = ColumnDataSource(data=dict(
            x=x_hover,
            y=y_hover,
            c_outer=c_hover,
            c=c_hover,
            desc1=desc1_hover,
            desc2=desc2_hover,
            desc3=desc3_hover,
        ))
        # hover = HoverTool(tooltips=[
        #     ("#indel", "@desc1"),
        #     ("#mut", "@desc3"),
        #     ("val", "@desc2"),
        # ])
        if desc2 != None:
            source = ColumnDataSource(data=dict(
                x=x_hover,
                y=y_hover,
                c_outer=c_hover,
                c=c_inner,
                desc1=desc1_hover,
                desc2=desc2_hover,
                desc3=desc3_hover,
                desc4=desc4_hover,
            ))
            # hover = HoverTool(tooltips=[
            #     ("#indel", "@desc1"),
            #     ("#mut", "@desc3"),
            #     ("val", "@desc2"),
            #     ("fails", "@desc4"),
            # ])

        plot = figure(plot_width=200, plot_height=200) # , tools=[hover]
        plot.axis.visible = False
        plot.toolbar.logo = None
        plot.toolbar_location = None
        plot.grid.grid_line_color = None

        plot.xaxis.formatter = tick_formater
        plot.xaxis.ticker = FixedTicker(ticks=[0, 4, 8, 12, 16, 20])
        if show_coverage:
            plot.rect('x', 'y', width, height, color='c', line_alpha=0, source=source)
        else:
            plot.rect('x', 'y', width, height, color='c', fill_alpha =0, line_alpha=0, source=source)
        plot.circle('x', 'y', radius=width*2/7, color='c_outer', line_alpha=0, source=source)

        #color_bar = ColorBar(color_mapper=color_mapper, border_line_color=None, location=(0,0))

        #plot.add_layout(color_bar, 'left')

        return plot

    files = [
        "paperGraphics/sw_human_200.db.html.json", #[done]
        "paperGraphics/sw_human_1000.db.html.json", #[done]
        "paperGraphics/human_30000.db.html.json",#[done]

        "paperGraphics/sw_zebrafish_200.db.html.json", #[done] # 3
        "paperGraphics/sw_zebrafish_1000.db.html.json", #[done]
        "paperGraphics/zebrafish_30000.db.html.json",

        "paperGraphics/sw_plasmodium_200.db.html.json", #[done] # 6
        "paperGraphics/sw_plasmodium_1000.db.html.json", #[done]
        "paperGraphics/plasmodium_30000.db.html.json", #[done]

        "paperGraphics/sw_human_1000_10.db.html.json", #[done] # 9
        "paperGraphics/human_30000_10.db.html.json", # #[done] 

        "paperGraphics/sw_eColi_200.db.html.json", #[done] # 11
        "paperGraphics/sw_eColi_1000.db.html.json", # #[done] 
        "paperGraphics/eColi_30000.db.html.json", # #[done] 
    ]

    for file in files[9:]:
        if not os.path.isfile(file):
            continue
        print file
        with open(file, "r") as f:
            json_file = json.loads(f.read(), object_hook=_decode)
            for approach, accuracy, coverage, runtime, alignments, fails, runtime_tup in json_file:
                if approach != "BWA MEM ont2d":
                    continue
                tot_runtime = ""
                if len(runtime_tup) > 0:
                    tot_runtime = str( runtime_tup[0][0] * 1000 )[:7] + "ms"
                avg_acc = 0
                count = 0
                for key, value in accuracy.items():
                    for key, value in value.items():
                        avg_acc += value
                        count += 1
                avg_acc = str(avg_acc * 100 / count)[:7] + "%"
                name = approach + "                  "
                print name[:15], tot_runtime, "\t", avg_acc
                plots[0].append(
                        makePicFromDict(accuracy, approach + tot_runtime, desc2=fails, inner=coverage, set_max=1)
                    )
                #plots[1].append(makePicFromDict(runtime, None)) #, set_max=50
                #plots[2].append(makePicFromDict(alignments, None, set_max=500))

    save(plots, "accuracyPics", True)
    #code from python test:
    """
        
            plots = [ [], [], [] ]

            def makePicFromDict(d, w, h, divideBy, title, ignore_max_n = 0, log = False, set_max = None, set_min=None, xax=True, yax=True):
                pic = []
                min_ = 10000.0
                max_ = 0.0
                if len(d.keys()) == 0:
                    return None

                for x in range(0,401,20):
                    pic.append( [] )
                    for y in range(0,20,2):
                        if x not in d or y not in d[x] or divideBy[x][y] == 0:
                            pic[-1].append( float("nan") )
                        else:
                            pic[-1].append( d[x][y] / divideBy[x][y] )
                            if pic[-1][-1] < min_:
                                min_ = pic[-1][-1]
                            if pic[-1][-1] > max_:
                                max_ = pic[-1][-1]

                for _ in range(ignore_max_n):
                    max_x = 0
                    max_y = 0
                    for x, row in enumerate(pic):
                        for y, p in enumerate(row):
                            if p > pic[max_x][max_y]:
                                max_x = x
                                max_y = y
                    pic[max_x][max_y] = float('nan')
                if ignore_max_n > 0:
                    max_ = 0
                    for row in pic:
                        for p in row:
                            if p > max_:
                                max_ = p

                if set_max is not None:
                    max_ = set_max
                if set_min is not None:
                    min_ = set_min

                color_mapper = LinearColorMapper(
                        palette=heatmap_palette(light_spec_approximation, 256),
                        low=min_,
                        high=max_
                    )
                if log:
                    color_mapper = LogColorMapper(
                            palette=heatmap_palette(light_spec_approximation, 256),
                            low=min_,
                            high=max_
                        )

                plot_width=500
                if yax:
                    plot_width += 150
                plot_height=500
                if xax:
                    plot_height += 40
                plot = figure(title=title,
                        x_range=(0,h), y_range=(0,w),
                        x_axis_label=str(indel_size) + 'nt indels', y_axis_label='mutations',
                        plot_width=plot_width, plot_height=plot_height,
                        min_border_bottom=10, min_border_top=10,
                        min_border_left=10, min_border_right=15,
                        tools=["save"]
                    )
                #plot.xaxis.formatter = tick_formater
                plot.image(image=[pic], color_mapper=color_mapper,
                        dh=[w], dw=[h], x=[0], y=[0])

                font = "Helvetica"
                font_size = '15pt'
                if yax:
                    color_bar = ColorBar(color_mapper=color_mapper, border_line_color=None, location=(0,0))
                    color_bar.major_label_text_font=font
                    color_bar.major_label_text_font_size=font_size
                    plot.add_layout(color_bar, 'left')

                #if not xax:
                #    plot.xaxis.visible = False
                #if not yax:
                #    plot.yaxis.visible = False
                if not title is None:
                    plot.title.text_font=font
                    plot.title.text_font_size=font_size
                plot.legend.label_text_font=font
                plot.legend.label_text_font_size=font_size
                plot.axis.axis_label_text_font=font
                plot.axis.major_label_text_font=font
                plot.axis.axis_label_text_font_size=font_size
                plot.axis.major_label_text_font_size=font_size

                return plot
            if indel_size > 0:
                avg_hits = makePicFromDict(hits, 400, max_indel, tries, approach, set_max=1, set_min=0, yax=False)
                avg_runtime = makePicFromDict(run_times, 400, max_indel, tries, None, 0, False, 0.2, 0, yax=False)
                #avg_score = makePicFromDict(scores, max_mut, max_indel, tries, "score " + approach, yax=yax)
                #avg_seeds = makePicFromDict(nums_seeds, 400, max_indel, tries, approach, yax=yax)
                #avg_rel = makePicFromDict(hits, 400, max_indel, nums_seeds, approach, yax=yax, set_min=0, set_max=0.01)
                yax = False

                if not avg_hits is None:
                    plots[0].append(avg_hits)
                if not avg_runtime is None:
                    plots[1].append(avg_runtime)
                #if not avg_seeds is None:
                #    plots[2].append(avg_seeds)
                #if not avg_rel is None:
                #    plots[3].append(avg_rel)

        plot = figure(title="BWA-pic",
                x_axis_label='#wrong / #mapped', y_axis_label='#mapped / total',
                x_axis_type="log", y_range=(0,0.4),
                plot_width=650, plot_height=500,
                min_border_bottom=10, min_border_top=10,
                min_border_left=10, min_border_right=15
            )


        c_palette = heatmap_palette(light_spec_approximation, len(approaches))

        for index, approach_, in enumerate(approaches):
            approach = approach_[0]
            data = mapping_qual[index]

            total_amount_1 = len(data[0])
            total_amount_2 = len(data[1])
            all_data = []
            min_val = 1
            max_val = 0
            for x in data[0]:
                if x is None:
                    continue
                all_data.append(x)
                if x < min_val:
                    min_val = x
                if x > max_val:
                    max_val = x
            for x in data[1]:
                if x is None:
                    continue
                all_data.append(x)
                if x < min_val:
                    min_val = x
                if x > max_val:
                    max_val = x

            values = []
            amount = 1000
            step = len(all_data) / amount
            c = 0
            for v in sorted(all_data):
                if c >= step:
                    c = 0
                    values.append(v)
                c += 1

            line_x = []
            line_y = []
            for val in values:
                mapped = 0
                wrong = 0
                total = total_amount_1 + total_amount_2
                for ele in data[0]:
                    if not ele is None and ele > val:
                        mapped += 1
                for ele in data[1]:
                    if not ele is None and ele > val:
                        mapped += 1
                        wrong += 1

                if mapped > 0:
                    line_x.append( wrong/mapped )
                    line_y.append( mapped/total )

            plot.line(line_x, line_y, legend=approach, color=c_palette[index])
            plot.x(line_x, line_y, legend=approach, color=c_palette[index])

        plot.legend.location = "top_left"
        plots[-1].append(plot)

        save(gridplot(plots))
    """


# actually call the functions that create the pictures

#seed_relevance_average()
#seed_relevance_pics()
accuracy_pics()
#unrelated_non_enclosed_seeds()
#forced_gap()
#ambiguity_per_length()
#theoretical_max_acc()
#seed_shadows()
#alignment()
#stripOfConsideration()
#optimal_matching()
#contradicting_seeds()
#required_nmw_band_size()

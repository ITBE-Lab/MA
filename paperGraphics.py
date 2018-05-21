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
        max_ = 0.0
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
        "paperGraphics/human_30000.db.html.json", # [missing bowtie & blasr pic]

        "paperGraphics/sw_zebrafish_200.db.html.json", #[done] # 3
        "paperGraphics/sw_zebrafish_1000.db.html.json", #[done]
        "paperGraphics/zebrafish_30000.db.html.json",

        "paperGraphics/sw_plasmodium_200.db.html.json", #[done] # 6
        "paperGraphics/sw_plasmodium_1000.db.html.json", #[done]
        "paperGraphics/plasmodium_30000.db.html.json", # [missing bowtie & blasr pic]

        "paperGraphics/sw_human_1000_10.db.html.json", #[done] # 9
    ]

    with open(files[5], "r") as f:
        json_file = json.loads(f.read(), object_hook=_decode)
        for approach, accuracy, coverage, runtime, alignments, fails, runtime_tup in json_file:
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


def smems():
    return [[0.9719981015662079, 0.08672167950231531, 0.05403315879230418, 0.041792016199016485, 0.031350716930571204, 0.022762294802288127, 0.01682020292602066, 0.012207627796844098, 0.008715509852365261, 0.0060908817865408735, 0.0042696243911077136, 0.0030579080900885345, 0.002121632601699336, 0.0014566951593409582, 0.000948664072616158, 0.0006253146513137031, 0.00041226271362989475, 0.0002498569367539554, 0.00017258663298167276, 9.818901870985691e-05], [0.051225022023357845, 0.06126125744225825, 0.047492970720413245, 0.03514241881252127, 0.02632680435853358, 0.01887796816560842, 0.013749493456845778, 0.009739768989055578, 0.006709526575511839, 0.0047251664451157535, 0.003270067619394893, 0.002272923671163612, 0.0014789889046291761, 0.0009817094411494595, 0.0006080374596535891, 0.00040360798044188203, 0.00024270600521354732, 0.0001814712401465365, 9.594282531707897e-05, 5.8392455694724244e-05], [0.04007890685142417, 0.05173044986921481, 0.040743432847211186, 0.030342501706215743, 0.021241494914630867, 0.015589457251649886, 0.01104737872999926, 0.007729150596197675, 0.005247162299188104, 0.003655686317641883, 0.0023799557655346357, 0.0015668173355379048, 0.001001834136954437, 0.0006174521820734938, 0.00045469360741841016, 0.00026758775773733686, 0.00015441901812515935, 9.369065507785334e-05, 5.584355250541504e-05, 2.876998339560958e-05], [0.035861773066127964, 0.04365774594798341, 0.03467713820776304, 0.025657788777588342, 0.018325117736841517, 0.01272805647195045, 0.008792940225299695, 0.005940276749680226, 0.004162449046589883, 0.002617183706875425, 0.0017107624824050296, 0.0010845938660615686, 0.0007123972762531531, 0.0004940959874564175, 0.00028225353651541315, 0.00013969499524615523, 8.149324171192912e-05, 5.144772113206275e-05, 2.2266586116314788e-05, 1.6220243095100443e-05], [0.030693654666824322, 0.037260036513176095, 0.02979288008984773, 0.021555814647108266, 0.014615294179841578, 0.010259592387539084, 0.006694200787827646, 0.004359350623061346, 0.002996360046232588, 0.0018172642686866522, 0.0011884765816379103, 0.0007202496039556133, 0.0004286799784372131, 0.00023780306275822508, 0.00014625878970553032, 9.026495418906083e-05, 4.5689489542964525e-05, 2.5529845840346788e-05, 1.324761168600626e-05, 5.734375689917075e-06], [0.028229697995162135, 0.03255061208761893, 0.02435141140754387, 0.017458567380245708, 0.011673248673031286, 0.007812654382849706, 0.005146792759801828, 0.0034161631405059626, 0.0021155659122793187, 0.001276596548279333, 0.0007994734904791372, 0.00043139584051308, 0.0002527151839706369, 0.00013221793804338893, 6.798442151760147e-05, 4.7229577143683456e-05, 2.0184707373646615e-05, 6.8796215290876135e-06, 4.5511564204016965e-06, 2.8337538624065144e-06], [0.02655455973390003, 0.027552153426679774, 0.02018972123041317, 0.013977585605149882, 0.009339918147106569, 0.005896245622607705, 0.0037990317453966585, 0.002304129783323549, 0.0014193596667373435, 0.0008670460018456131, 0.00044564759902827086, 0.0002595024847063945, 0.0001325374773384581, 7.946979621840637e-05, 3.4708021660119384e-05, 1.552461110849173e-05, 6.865939105986069e-06, 4.533004237225711e-06, 2.26182808862521e-06, 1.6798328230374513e-06], [0.023588213509971347, 0.02256752033880249, 0.01663694381365193, 0.011320743545488725, 0.007194493972930179, 0.0044877977380663, 0.002718426252303697, 0.001608902468794347, 0.0009094843345901308, 0.0004905468905647067, 0.0002744170415946297, 0.00012798384350746182, 6.033872069148174e-05, 3.6255111683008325e-05, 1.3155260674350102e-05, 6.223771752789522e-06, 2.810847397408286e-06, 5.600025984120566e-07, 0.0, 0.0], [0.02175835776710551, 0.019538230959280593, 0.013439937437909626, 0.008820484783224081, 0.005357064608514404, 0.003224505483113225, 0.0018435607058203845, 0.0009978550947690528, 0.0004960083563381333, 0.0002655133099634313, 0.00012438504961452258, 6.396835468205134e-05, 3.2092113535020805e-05, 9.662351952875005e-06, 5.650477295817177e-06, 1.6808351733809495e-06, 1.1136700774668905e-06, 1.1079117083001422e-06, 0.0, 0.0], [0.018544421989716304, 0.01615860039694077, 0.010786562704529077, 0.006496404036814282, 0.003882486848361465, 0.0022037279987536947, 0.0010958231815319513, 0.0005681909074050265, 0.00026243198296066624, 0.0001451770002614343, 5.655465607628366e-05, 1.9308703000345286e-05, 1.1242700676585727e-05, 2.2491685105162686e-06, 1.6729084740621118e-06, 0.0, 1.1028598811006761e-06, 0.0, 0.0, 5.467878129745093e-07], [0.01646773105798122, 0.013211357461672813, 0.008297784598327463, 0.004838108702140591, 0.0025837549792394427, 0.001310124951079123, 0.0006553926600720227, 0.00029842793052737394, 0.00013813358033069523, 5.533421868083223e-05, 2.543884840029215e-05, 1.6896239629228917e-06, 2.2265107292768905e-06, 0.0, 0.0, 0.0, 0.0, 0.0, 5.419798632801597e-07, 0.0], [0.014299619697263123, 0.010542207821860152, 0.00614373334690163, 0.0032866381286704156, 0.0016276680961579892, 0.0007721565044241058, 0.00034123967994945966, 0.00012530507514032745, 5.984075584519987e-05, 1.856423111455143e-05, 1.0058130405895852e-05, 1.6664935365048186e-06, 1.65716283779193e-06, 0.0, 0.0, 0.0, 5.429574348518948e-07, 0.0, 5.409063101589237e-07, 5.404269048292008e-07], [0.011857564749757738, 0.008090297604773759, 0.004362972814000808, 0.0020273041056754456, 0.0009165748739057337, 0.0003646188102093267, 0.0001455048312679731, 5.849637322486006e-05, 1.392438058785393e-05, 6.071300246218822e-06, 5.501527499110128e-07, 1.096377623151987e-06, 0.0, 0.0, 0.0, 5.432400028900368e-07, 0.0, 0.0, 5.397600766459309e-07, 1.077157897344698e-06], [0.009751127724455811, 0.005852686983842518, 0.00275068604581048, 0.001207902031972392, 0.000444046221483024, 0.00015615791090939314, 4.709540476284139e-05, 2.0385293057871094e-05, 3.8388136968872705e-06, 2.730879474272931e-06, 1.083969713886194e-06, 0.0, 5.406594856273785e-07, 0.0, 5.382638300162825e-07, 0.0, 0.0, 5.382852707386512e-07, nan, nan], [0.007632016716502434, 0.004010537583831787, 0.0016676469610301461, 0.0006323664639677814, 0.0002076625822647247, 5.950308314123388e-05, 1.643251249281762e-05, 4.8999564992750786e-06, 1.0870746820306556e-06, 5.416302420328899e-07, 0.0, 0.0, 0.0, 0.0, 0.0, nan, nan, nan, nan, nan], [0.005588692828096903, 0.0025497031650361136, 0.0008836640894524517, 0.0002538848233049746, 6.21304532416564e-05, 2.2741633921421156e-05, 4.314761879078799e-06, 2.6857851597478156e-06, 0.0, 0.0, 0.0, 0.0, 5.329782987226109e-07, nan, nan, nan, nan, nan, nan, nan], [0.003920403817226153, 0.0013924070940292456, 0.00038824024655436783, 0.00011371089301852841, 2.744762884778083e-05, 7.496556938491821e-06, 1.0698843080603479e-06, 0.0, 0.0, 5.310874706375015e-07, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan], [0.002493895694638646, 0.0006781615315156884, 0.0001772635947028085, 3.7723819138196696e-05, 1.003993251051683e-05, 2.6462545178180255e-06, 5.272971171612011e-07, 5.270108890989905e-07, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan], [0.0013540130642521345, 0.0003484041825343194, 7.293391085797039e-05, 2.3001532738499757e-05, 6.2667897742702326e-06, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan], [0.0008747457526810827, 0.00018794701959399314, 5.679356931577689e-05, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan]]

def max_spannning():
    return [[0.9719981015662079, 0.21808112988372227, 0.13878809738503156, 0.11205134968785721, 0.08331148512221946, 0.05942213217865703, 0.04277537539130896, 0.0301382086350052, 0.02037343768666234, 0.01345214851475749, 0.008586467041270834, 0.0059832270924993864, 0.004077208763296715, 0.0025403087222247126, 0.001569411176911434, 0.0010587750697631244, 0.0005943859107678151, 0.0003534068419564603, 0.0002573332997681203, 0.00015413749064165235], [0.23995359831903346, 0.15184096787924428, 0.13203538901125067, 0.09711268960211035, 0.07350921726103783, 0.04974524484346074, 0.03394809858005549, 0.022802247485298428, 0.015091050583657587, 0.010179826216042964, 0.0068506219724621, 0.004539088419650669, 0.0027157567491915263, 0.0016193788587901634, 0.0010862929322170407, 0.0006329349566253397, 0.00043119483242585137, 0.0002738551334006445, 0.00014950580027132534, 7.623424605691976e-05], [0.18945562644670827, 0.15380968246239013, 0.1184828767747462, 0.0889690453968098, 0.058067660109261586, 0.04061877530039196, 0.027606979841762146, 0.01793138331338901, 0.011613754851679392, 0.007661227802960517, 0.004710941318643753, 0.003010149682785596, 0.0017878288998642521, 0.0011421178127442029, 0.0007663680376453554, 0.00041763718260274857, 0.00021290656166962984, 0.00015064696021824635, 6.50944547661346e-05, 3.212696576871797e-05], [0.17804452670356463, 0.13153637568175516, 0.10796350840034859, 0.07448061683443992, 0.05165169371135877, 0.033360401241939334, 0.021455588887439675, 0.013548010648772838, 0.009044908670112924, 0.005166515394819521, 0.0031643615253634165, 0.0019656033798017706, 0.001197335714732633, 0.0008286778923395212, 0.0004050054231890571, 0.0002137841437945041, 0.0001008394200370653, 7.829923266751986e-05, 2.1346710605312662e-05, 1.8485899220159136e-05], [0.13419290080121407, 0.11964873765093303, 0.09423543450139214, 0.06497874078166663, 0.03947053590600939, 0.026572799116540705, 0.016476162310321382, 0.009711212881902893, 0.006104446638563549, 0.003405936059990269, 0.0022184993110073016, 0.001168351645894418, 0.0006939972018032823, 0.00036566589684372596, 0.00022205252354569137, 0.00012081023396915312, 5.315727643378464e-05, 4.2175852215813837e-05, 2.874674639097666e-05, 1.0445418769895259e-05], [0.12544038858429402, 0.10850829609565676, 0.07317574706898407, 0.050519339954047435, 0.03120830604274313, 0.01888785060591228, 0.011677182960092715, 0.007156007317855625, 0.00408443060752335, 0.0024668438025162373, 0.001374530206207412, 0.0007620789513793629, 0.000432544897065265, 0.00017701235818586844, 0.00011817124625007722, 7.951295663628388e-05, 1.3175404155522471e-05, 1.0415229148060294e-05, 1.2963342260755037e-05, 5.161423520606983e-06], [0.1265599704911986, 0.0963974037115954, 0.062933807618962, 0.04013943365835065, 0.024531630931731073, 0.01423307353022955, 0.008139214966984594, 0.004677417519261628, 0.0026310717846878325, 0.0015622615331252168, 0.0007456510516957416, 0.0003988046785798183, 0.00022162581454243728, 0.00011471499991996628, 5.53371349066647e-05, 2.089563934126497e-05, 7.773490841532207e-06, 2.5695051133151754e-06, 2.5634715557196175e-06, 2.547692809391815e-06], [0.10127498077606807, 0.07430094928826443, 0.052212328033325144, 0.032502108328049686, 0.018371126943925013, 0.010206807190229653, 0.005798354698958014, 0.003079634136108081, 0.001631164315736197, 0.0008551211649318362, 0.0004175590307242629, 0.00017673237897848686, 0.00010517210099676858, 6.55498716533513e-05, 5.209744305749474e-06, 1.0285235287613749e-05, 0.0, 0.0, 0.0, 0.0], [0.10361608883215025, 0.0683712379578924, 0.04155565478777126, 0.024139896252529717, 0.012829148546997467, 0.0070146563210282935, 0.0036677961413456883, 0.001881592701701036, 0.0008905323917582313, 0.00045617622355855023, 0.00021036079505862493, 9.682872829858839e-05, 4.951823965264258e-05, 1.2912922995075011e-05, 5.129467766424556e-06, 0.0, 2.5313059260403034e-06, 2.5161220520485006e-06, 0.0, 0.0], [0.08023706485749987, 0.05641273204892698, 0.03279433112567054, 0.017048903734399625, 0.009184435596637334, 0.0045849703365488026, 0.0020867202758289597, 0.0010037061985567026, 0.00041065476915903203, 0.00021665323592744758, 6.770886230882012e-05, 3.622288458354032e-05, 1.020866511494957e-05, 7.678623990580887e-06, 0.0, 0.0, 4.990891622788411e-06, 0.0, 0.0, 0.0], [0.0721026763940476, 0.04468444184701259, 0.024467917839939395, 0.012214130573655295, 0.005650895937337904, 0.0026390187794941254, 0.0011225613551603241, 0.0005034145210050367, 0.00023741256825611338, 7.035537280530323e-05, 3.582248423810693e-05, 0.0, 7.573099341392794e-06, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.06099096737900571, 0.03529159345661616, 0.016815048191977516, 0.007625468409942563, 0.003398386635329852, 0.0014510011119084607, 0.0006285174811314358, 0.00015270652911000563, 8.72976748007559e-05, 2.0505722378151154e-05, 1.0171775864918817e-05, 0.0, 5.028802466124729e-06, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.049096093834068655, 0.02587199439342785, 0.01134350188610253, 0.00458315203776222, 0.0017576035457741097, 0.0006163817093279098, 0.00025223731928354306, 9.481855060876072e-05, 2.530428401528379e-05, 5.016919561220215e-06, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.038769424836196435, 0.01806058912639125, 0.006694435398304462, 0.0025239466189149345, 0.0007756760285387339, 0.0002540127666816534, 6.81044267877412e-05, 4.2587197286443994e-05, 0.0, 7.430819074417176e-06, 4.932912391475927e-06, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, nan, nan], [0.029360861968227243, 0.011189671877864195, 0.0037810197950532945, 0.0012616691548872334, 0.0003359992669106904, 0.00010517195614830247, 2.739616852130426e-05, 7.4331941674203096e-06, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, nan, nan, nan, nan, nan], [0.019684044396765143, 0.006606951796192856, 0.001703042955652559, 0.00043497933848142214, 9.869670994517398e-05, 3.4667795846302805e-05, 7.361927058026709e-06, 2.44179374168264e-06, 0.0, 0.0, 0.0, 0.0, 0.0, nan, nan, nan, nan, nan, nan, nan], [0.01290818403136795, 0.003430514577186578, 0.0006876816255906297, 0.0001930030294146389, 6.609611844424425e-05, 1.9368349695553752e-05, 0.0, 0.0, 0.0, 2.426789818160639e-06, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan], [0.007764319077830224, 0.0015671098227620862, 0.00038005306125432126, 6.284853744201618e-05, 1.9203763937731795e-05, 4.812400593850234e-06, 0.0, 0.0, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan], [0.0041353056551271025, 0.0007542608485434547, 0.0001655057279373669, 4.5574696928265424e-05, 1.4324561131258341e-05, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan], [0.002530909709618875, 0.0004473346705784656, 9.869998026000394e-05, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan]]

def kmer():
    return [[0.03421827096000926, 0.02840721119255558, 0.02241906668237961, 0.018315835771568128, 0.01574332534095572, 0.012455221556962166, 0.011063162015927813, 0.008991200311309348, 0.007438158023613808, 0.005780514702182143, 0.00486702110739045, 0.003933853561797531, 0.003167231646968023, 0.002577190465976076, 0.0019312016930517206, 0.001439888027814221, 0.001015922928633557, 0.0006564973157892624, 0.000490829739758453, 0.0003039225976109881], [0.03512124861492885, 0.027782765909246863, 0.021407403496180887, 0.017331213258294773, 0.014671231927164697, 0.012203688057257352, 0.00993651862517369, 0.008350284204322947, 0.006314664099000141, 0.005416196296489374, 0.0044845553863438236, 0.0034096429907544946, 0.0026341847626683815, 0.0018781133168093308, 0.0014411967042311763, 0.0010499868269706531, 0.000751371664076617, 0.0005290008791097199, 0.0003378822913312364, 0.0001867104651974729], [0.03507034857300553, 0.027877752375115285, 0.02051300285275768, 0.017420075263483645, 0.013427528856472588, 0.01107101006915953, 0.009146576433795446, 0.007604315429036893, 0.005798600812507547, 0.00460120301263743, 0.0037515471620370514, 0.002823286091932814, 0.0021575697388916525, 0.0014592338618151562, 0.001100169095402535, 0.0008143443283063333, 0.0004606267692393931, 0.00030336435056411797, 0.00021499149842012634, 0.00010237272446510252], [0.03466681081520236, 0.026904927991947596, 0.020573844378341963, 0.01684192823858025, 0.013466416566190492, 0.010854133760041998, 0.008540076704802993, 0.006639919492998363, 0.00539063735127829, 0.003994144220336255, 0.0030034761756240333, 0.002271065025739962, 0.001777086456301444, 0.0013187368714645562, 0.000792343284046627, 0.0004953839801130498, 0.0002769004149965299, 0.00020175164216452925, 0.00011896339622471465, 4.666586945806342e-05], [0.03527566727771075, 0.02622736187191452, 0.020691115145031344, 0.016657037376986466, 0.01241639555651072, 0.010367488264316842, 0.0073490905902362244, 0.00598700358647372, 0.004492407777851579, 0.0032078403200948924, 0.0025161342152478375, 0.0017638649262161897, 0.0012284082971983214, 0.0007636410662545464, 0.0005185411186084242, 0.0003135194489784244, 0.0001471942669935781, 8.221481159875476e-05, 4.770804336868568e-05, 1.4859582254686871e-05], [0.03546698543342871, 0.025763038439761677, 0.018973490958435642, 0.015554933271447463, 0.01126421170596498, 0.008496191040580971, 0.006457176567316086, 0.005253008959699594, 0.003943015340797381, 0.0028249201356273655, 0.001977928080776825, 0.0012406519400130175, 0.0008366428040653528, 0.00047434272939832813, 0.00028911074228047827, 0.00019280725643881555, 9.32058800155801e-05, 2.8163930766799525e-05, 2.4982625719386063e-05, 1.0935721046767218e-05], [0.03341839855181572, 0.0251383466911308, 0.017815914930558525, 0.013720495487097278, 0.010190001095899657, 0.007833526523016313, 0.005707135809541997, 0.004112595521666822, 0.0030883652578391355, 0.0021117003552561368, 0.001331451681269417, 0.0008774178558722489, 0.0005389315832746856, 0.00028682506476946436, 0.0001345576533145262, 5.9831812774291463e-05, 3.13172578609017e-05, 1.3044411477659074e-05, 1.2257350426617082e-05, 7.8904905123427e-06], [0.03061865547872329, 0.022372864942496695, 0.01670159472892683, 0.013139141936141043, 0.009544482716937266, 0.007237892751869816, 0.005197739577102439, 0.003801855479530981, 0.0024987223247036573, 0.0015059438911710515, 0.0009971196300211354, 0.0005092296142551172, 0.0002404333688803983, 0.00015650060753341435, 4.8848150436020626e-05, 2.2741537873757176e-05, 2.449134537082958e-06, 1.2954609638747755e-06, 0.0, 0.0], [0.03317621745390254, 0.02261470397350489, 0.01632400946874741, 0.012440848240702414, 0.008579423827839869, 0.005738339134609664, 0.004217967967362133, 0.002621938300577761, 0.0016432530530720276, 0.0009271010643811226, 0.0005234220709165789, 0.00032224364738229013, 0.00011134682211328523, 3.643329119583998e-05, 1.1814186475119323e-05, 6.445790898543251e-06, 1.3130545193366973e-06, 1.405449772036047e-06, 0.0, 0.0], [0.02905251582832726, 0.021351153369096185, 0.014903624715164753, 0.010190208735157297, 0.007511883552129491, 0.005278474076036175, 0.003061015886280917, 0.0019412316789720525, 0.001018602377935648, 0.0005635178395640681, 0.0002241112013415033, 7.853624753809263e-05, 5.5712840104148904e-05, 1.0796311979827692e-05, 1.2927512849947773e-06, 0.0, 1.4213186425838436e-06, 0.0, 0.0, 0.0], [0.032171734345669514, 0.01955113328245462, 0.013624409274373396, 0.009055559674063685, 0.005750784821717742, 0.003604157889749023, 0.0023022303251730852, 0.0011424248512461503, 0.0005349967471585473, 0.000276182215946118, 0.00010770982763290399, 1.4346670556208461e-05, 6.209398296886235e-06, 0.0, 0.0, 0.0, 0.0, 0.0, 3.179746288043677e-06, 0.0], [0.03168226814392149, 0.01988819166512157, 0.012356110283247521, 0.007830426425710425, 0.004447326867313941, 0.002749095590202518, 0.0013824775729812508, 0.0006105852655126579, 0.00027854323807474754, 9.685720386460243e-05, 4.7720888060610553e-05, 1.566841456535818e-05, 8.270283905062654e-06, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.029534223162692275, 0.016733187707642982, 0.010886475536532464, 0.0060005735899810175, 0.0032624417986717954, 0.0015076267153306647, 0.0006615693120470505, 0.0003496024405675851, 6.450225757901526e-05, 1.1605505651881252e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.744198547789609e-06], [0.028173148414812138, 0.015316953270069494, 0.008377038373783622, 0.004642610358333277, 0.002094668452708899, 0.0008482349644115537, 0.00024053729272295173, 8.017871218186012e-05, 2.523046033638843e-05, 7.032744458197367e-06, 4.396789757238581e-06, 0.0, 4.737532787674835e-06, 0.0, 0.0, 0.0, 0.0, 0.0, nan, nan], [0.026201686301426875, 0.013060784416900091, 0.006340209187734565, 0.003023568228083526, 0.001144761558110664, 0.00042553600964833407, 0.00010144922090755722, 3.483265588283366e-05, 2.881707353973082e-06, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, nan, nan, nan, nan, nan], [0.02364020805178355, 0.01090423376506431, 0.00428500283936953, 0.0014724838477537109, 0.00043187807976472126, 0.00018508353479126762, 2.7487063900551804e-05, 8.568968053458936e-06, 0.0, 0.0, 0.0, 0.0, 0.0, nan, nan, nan, nan, nan, nan, nan], [0.02098161696702209, 0.008185795455325435, 0.0023232444018945436, 0.0008810645151009523, 0.0002248165336502176, 4.7874040583105814e-05, 1.4740739867215415e-06, 0.0, 0.0, 0.0, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan], [0.017550879838695608, 0.004978024465528048, 0.0017812578182021185, 0.0003022597035356891, 8.233685492111184e-05, 3.1777763654327264e-05, 0.0, 4.6559631247720515e-06, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan], [0.014686352118777333, 0.0037034705323636245, 0.0007651245763003851, 0.00019472526035992853, 8.341519048146385e-05, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan], [0.010886950227798642, 0.001994101445803058, 0.000625910889027691, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan]]


def seed_relevance_pics():
    h = 400
    w = 20


    data = []
    high = 99
    low = 1
    for row in kmer():
        data.append( [] )
        for ele in row:
            data[-1].append( math.log(ele * high + low) )

    color_mapper = LinearColorMapper(
                    palette=heatmap_palette(light_spec_approximation, 256),
                    low=math.log(low),
                    high=math.log(high+low)
                )

    plot = figure(
            x_range=(0,w), y_range=(0,h),
            plot_width=resolution, plot_height=resolution
        )
    plot.axis.visible = False
    plot.toolbar.logo = None
    plot.toolbar_location = None

    plot.image(image=[data], color_mapper=color_mapper,
                dh=[h], dw=[w], x=[0], y=[0])

    print("color scale:")
    for x_ in [1.0, 0.8, 0.6, 0.4, 0.2, 0.0]:
        #adjust x_ to the linear scale
        x = x_ * ( math.log(high+low) - math.log(low) ) + math.log(low)
        #    x = log( y * high + low )
        # => e^x = y * high + low
        # => (e^x - low) / high = y
        print(x_, "entspricht", 100 * (math.e**x - low) / high, "%" )

    save(plot, "relevancePics")


# actually call the functions that create the pictures

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

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
        export_png(plot, filename="paperGraphics/" + name + ".png")
        plot.output_backend = "svg"
        export_svgs(plot, filename="paperGraphics/" + name + ".svg")
    #show(plot)

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

nan = float("nan")
def get_accuracy_data():
    return [(400, 20.0, 'MA Accurate', [[1.0, 0.96875, 1.0, 0.9375, 0.96875, 0.96875, 0.96875, 0.9375, 0.90625, 0.59375], [1.0, 0.9375, 0.96875, 0.96875, 0.90625, 0.96875, 0.9375, 0.84375, 0.5625, 0.1875], [0.90625, 0.9375, 0.9375, 1.0, 0.84375, 0.96875, 0.84375, 0.375, 0.125, 0.0], [1.0, 0.96875, 0.9375, 0.875, 0.875, 0.875, 0.4375, 0.21875, 0.03125, 0.0], [0.9375, 0.9375, 0.96875, 0.96875, 0.75, 0.625, 0.1875, 0.0625, 0.0, 0.0], [0.90625, 0.96875, 0.875, 0.875, 0.65625, 0.375, 0.15625, 0.0, 0.0, nan], [0.96875, 0.9375, 0.8125, 0.6875, 0.25, 0.15625, 0.125, 0.0, 0.0, nan], [0.90625, 0.9375, 0.71875, 0.53125, 0.25, 0.09375, 0.0, 0.0, 0.0, nan], [0.9375, 0.78125, 0.53125, 0.15625, 0.03125, 0.03125, 0.0, 0.0, 0.0, nan], [0.78125, 0.5625, 0.375, 0.15625, 0.15625, 0.0, 0.0, 0.0, 0.0, nan], [0.5625, 0.3125, 0.15625, 0.09375, 0.0, 0.0, 0.0, 0.0, nan, nan], [0.59375, 0.34375, 0.0625, 0.0625, 0.0, 0.0, 0.0, 0.0, nan, nan], [0.375, 0.125, 0.0625, 0.0, 0.0, 0.0, 0.0, 0.0, nan, nan], [0.40625, 0.09375, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, nan, nan], [0.09375, 0.03125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, nan, nan], [0.125, 0.03125, 0.03125, 0.0, 0.0, 0.0, 0.0, nan, nan, nan], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, nan, nan, nan], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, nan, nan, nan], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, nan, nan, nan], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, nan, nan, nan]], [[3.9218478202819824, 3.94997239112854, 4.056169509887695, 4.057829856872559, 4.080075025558472, 4.257440567016602, 4.240171194076538, 4.467890977859497, 4.418534755706787, 4.238884925842285], [4.216655731201172, 4.091922760009766, 4.109312534332275, 4.415005683898926, 4.354067087173462, 4.493427276611328, 4.560627460479736, 4.421206712722778, 4.301068544387817, 4.258890390396118], [4.297042369842529, 4.299480676651001, 4.486332416534424, 4.318108797073364, 4.370590925216675, 4.736164331436157, 4.440006256103516, 4.3853371143341064, 4.32783579826355, 4.361993074417114], [4.184879541397095, 4.297570466995239, 4.45026707649231, 4.43636679649353, 4.421149730682373, 4.482136964797974, 4.2955322265625, 4.272418260574341, 4.3274548053741455, 4.329248905181885], [4.278815269470215, 4.439603567123413, 4.564253330230713, 4.5324015617370605, 4.483762502670288, 4.313724994659424, 4.300660133361816, 4.3334643840789795, 4.23433780670166, 4.478815317153931], [4.484564542770386, 4.582150459289551, 4.606868743896484, 4.680434465408325, 4.570631742477417, 4.474807262420654, 4.304866313934326, 4.211396217346191, 4.262079954147339, nan], [4.468990087509155, 4.5177083015441895, 4.506632328033447, 4.524904012680054, 4.213464975357056, 4.358732461929321, 4.2207207679748535, 4.215402603149414, 4.236968040466309, nan], [4.600229024887085, 4.537827491760254, 4.54267954826355, 4.448955774307251, 4.32643723487854, 4.305803537368774, 4.275101661682129, 4.329741477966309, 4.365920305252075, nan], [4.503157615661621, 4.645795822143555, 4.305184602737427, 4.303131580352783, 4.341109752655029, 4.213728904724121, 4.264341831207275, 4.261372089385986, 4.490002393722534, nan], [4.490204095840454, 4.49869704246521, 4.617083787918091, 4.302822589874268, 4.284229278564453, 4.300840616226196, 4.311132907867432, 4.695278167724609, 4.325638294219971, nan], [4.350025415420532, 4.223368883132935, 4.1948912143707275, 4.185719013214111, 4.226929426193237, 4.301138162612915, 4.231539249420166, 4.222170352935791, nan, nan], [4.358469009399414, 4.234894037246704, 4.17148232460022, 4.2356116771698, 4.206381559371948, 4.216463088989258, 4.286026477813721, 4.3179895877838135, nan, nan], [4.274434328079224, 4.298497915267944, 4.229882478713989, 4.244692087173462, 4.219872951507568, 4.2518017292022705, 4.236212253570557, 4.354483127593994, nan, nan], [4.290680646896362, 4.4145684242248535, 4.318278074264526, 4.1998841762542725, 4.3009796142578125, 4.297672510147095, 4.311764240264893, 4.264423131942749, nan, nan], [4.200058221817017, 4.365768671035767, 4.366191625595093, 4.195141315460205, 4.275837659835815, 4.406377553939819, 4.23265814781189, 4.341414928436279, nan, nan], [4.204110622406006, 4.341148376464844, 4.241948366165161, 4.260713815689087, 4.479329347610474, 4.229283571243286, 4.2236247062683105, nan, nan, nan], [4.190190076828003, 4.3049116134643555, 4.321344614028931, 4.285593509674072, 4.217620134353638, 4.454849004745483, 4.303012371063232, nan, nan, nan], [4.283354997634888, 4.304572343826294, 4.202645301818848, 4.22454833984375, 4.32361102104187, 4.246096134185791, 4.314420223236084, nan, nan, nan], [4.382885932922363, 4.268249034881592, 4.2031190395355225, 4.226048946380615, 4.23385214805603, 4.321706295013428, 4.228805065155029, nan, nan, nan], [4.170299768447876, 4.250947713851929, 4.308153867721558, 4.234647989273071, 4.253759860992432, 4.341252565383911, 4.297029256820679, nan, nan, nan]], ([None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None], [None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None])), (400, 20.0, 'MA Fast', [[1.0, 1.0, 1.0, 0.96875, 1.0, 1.0, 0.96875, 0.96875, 0.9375, 0.65625], [1.0, 1.0, 0.96875, 1.0, 0.9375, 1.0, 0.90625, 0.84375, 0.4375, 0.1875], [0.9375, 0.9375, 0.9375, 0.96875, 0.875, 0.9375, 0.875, 0.375, 0.09375, 0.0], [1.0, 0.96875, 1.0, 0.90625, 0.8125, 0.90625, 0.40625, 0.1875, 0.03125, 0.0], [0.9375, 0.90625, 0.9375, 0.9375, 0.78125, 0.5625, 0.15625, 0.0625, 0.0, 0.0], [0.90625, 0.96875, 0.78125, 0.71875, 0.40625, 0.34375, 0.09375, 0.0, 0.0, nan], [0.9375, 0.75, 0.71875, 0.625, 0.1875, 0.25, 0.125, 0.0, 0.0, nan], [0.75, 0.78125, 0.5625, 0.375, 0.125, 0.0625, 0.0, 0.0, 0.0, nan], [0.78125, 0.625, 0.40625, 0.125, 0.0, 0.0625, 0.0, 0.0, 0.0, nan], [0.5625, 0.40625, 0.28125, 0.0625, 0.125, 0.0, 0.0, 0.0, 0.0, nan], [0.46875, 0.28125, 0.09375, 0.03125, 0.0, 0.0, 0.0, 0.0, nan, nan], [0.375, 0.25, 0.03125, 0.09375, 0.0, 0.0, 0.0, 0.0, nan, nan], [0.3125, 0.09375, 0.0625, 0.0, 0.0, 0.0, 0.0, 0.0, nan, nan], [0.28125, 0.09375, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, nan, nan], [0.09375, 0.0625, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, nan, nan], [0.03125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, nan, nan, nan], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, nan, nan, nan], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, nan, nan, nan], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, nan, nan, nan], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, nan, nan, nan]], [[3.8617701530456543, 3.8627562522888184, 3.897467851638794, 3.94268798828125, 3.9755184650421143, 4.066900014877319, 4.174433708190918, 4.168971538543701, 4.068197011947632, 3.913459300994873], [3.9196135997772217, 3.976292133331299, 4.021225452423096, 4.062448024749756, 4.155272483825684, 4.186021089553833, 4.459471702575684, 4.140146732330322, 3.928741455078125, 3.8986704349517822], [3.95206356048584, 4.054339408874512, 4.168201684951782, 4.152499437332153, 4.216059923171997, 4.284073829650879, 3.9055826663970947, 4.099047899246216, 3.8781213760375977, 3.9173076152801514], [4.088428497314453, 4.138639688491821, 4.268385648727417, 4.331273555755615, 4.216109991073608, 4.15776515007019, 3.9110870361328125, 3.924645185470581, 3.973381280899048, 3.962174892425537], [4.159902095794678, 4.248851299285889, 4.371466398239136, 4.343582391738892, 4.035197734832764, 4.063980579376221, 3.93902325630188, 3.8902547359466553, 3.959409236907959, 3.9109935760498047], [4.368624448776245, 4.3749706745147705, 4.279805421829224, 4.115503549575806, 3.9841761589050293, 3.9648728370666504, 4.015820026397705, 3.9702038764953613, 3.9713144302368164, nan], [4.279850721359253, 4.344956636428833, 4.242412328720093, 4.010035276412964, 3.963073492050171, 3.966193675994873, 3.902120590209961, 3.9641339778900146, 3.925208806991577, nan], [4.289578914642334, 4.213010549545288, 3.9485597610473633, 3.965155839920044, 3.9145991802215576, 3.9280714988708496, 3.963259696960449, 3.9258382320404053, 3.8931751251220703, nan], [3.9948768615722656, 4.0842955112457275, 3.8930695056915283, 3.922520875930786, 3.9424421787261963, 3.912259817123413, 3.906968355178833, 3.9459869861602783, 3.9172849655151367, nan], [3.998366355895996, 4.004817962646484, 3.9731104373931885, 3.9531126022338867, 3.908541440963745, 3.890634536743164, 3.934767007827759, 3.9509117603302, 3.9135208129882812, nan], [3.9038119316101074, 3.9158589839935303, 3.8892221450805664, 3.9071216583251953, 3.8962562084198, 3.9563965797424316, 3.93165922164917, 4.028118133544922, nan, nan], [3.9203667640686035, 3.9279403686523438, 3.9132080078125, 3.9068009853363037, 3.888261318206787, 3.9080748558044434, 3.9628865718841553, 3.8951351642608643, nan, nan], [3.896059989929199, 3.9299676418304443, 3.92405104637146, 3.9567909240722656, 3.912309408187866, 3.8923375606536865, 3.9018540382385254, 3.934128761291504, nan, nan], [3.9971001148223877, 3.891047477722168, 3.8804919719696045, 3.8868253231048584, 3.9071640968322754, 3.9590206146240234, 3.9294674396514893, 3.9028544425964355, nan, nan], [3.9328298568725586, 3.9610538482666016, 3.907522439956665, 3.9147276878356934, 3.925855875015259, 3.9033000469207764, 3.9451510906219482, 3.939100980758667, nan, nan], [3.931523084640503, 3.933715581893921, 3.9236276149749756, 3.966701030731201, 3.8926219940185547, 3.9657015800476074, 3.904157876968384, nan, nan, nan], [3.901181221008301, 3.9346792697906494, 3.900826930999756, 3.959944725036621, 3.9868242740631104, 3.8977408409118652, 3.9738049507141113, nan, nan, nan], [3.9247257709503174, 3.954254388809204, 4.039129257202148, 3.9534945487976074, 3.8998773097991943, 3.936549425125122, 3.9527087211608887, nan, nan, nan], [3.894724130630493, 3.9030275344848633, 3.917459726333618, 3.9152777194976807, 3.9238998889923096, 3.909499406814575, 3.9124739170074463, nan, nan, nan], [3.930255651473999, 3.9127771854400635, 3.9634478092193604, 3.9370079040527344, 3.897848129272461, 3.9063377380371094, 3.8870956897735596, nan, nan, nan]], ([None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None], [None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None])), (400, 20.0, 'MA Fast PY', [[1.0, 1.0, 1.0, 0.96875, 1.0, 1.0, 0.96875, 0.96875, 0.9375, 0.65625], [1.0, 1.0, 0.96875, 1.0, 0.9375, 1.0, 0.90625, 0.84375, 0.4375, 0.1875], [0.9375, 0.9375, 0.9375, 0.96875, 0.875, 0.9375, 0.875, 0.375, 0.09375, 0.0], [1.0, 0.96875, 1.0, 0.90625, 0.8125, 0.90625, 0.40625, 0.1875, 0.03125, 0.0], [0.967741935483871, 0.90625, 0.9375, 0.9375, 0.78125, 0.5625, 0.15625, 0.0625, 0.0, 0.0], [0.90625, 0.96875, 0.78125, 0.71875, 0.40625, 0.34375, 0.09375, 0.0, 0.0, nan], [0.9375, 0.75, 0.71875, 0.625, 0.1875, 0.25, 0.125, 0.0, 0.0, nan], [0.75, 0.78125, 0.5625, 0.375, 0.125, 0.0625, 0.0, 0.0, 0.0, nan], [0.78125, 0.625, 0.40625, 0.125, 0.0, 0.0625, 0.0, 0.0, 0.0, nan], [0.5625, 0.40625, 0.28125, 0.0625, 0.125, 0.0, 0.0, 0.0, 0.0, nan], [0.46875, 0.28125, 0.09375, 0.03125, 0.0, 0.0, 0.0, 0.0, nan, nan], [0.375, 0.25, 0.03125, 0.09375, 0.0, 0.0, 0.0, 0.0, nan, nan], [0.3125, 0.09375, 0.0625, 0.0, 0.0, 0.0, 0.0, 0.0, nan, nan], [0.28125, 0.09375, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, nan, nan], [0.09375, 0.0625, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, nan, nan], [0.03125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, nan, nan, nan], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, nan, nan, nan], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, nan, nan, nan], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, nan, nan, nan], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, nan, nan, nan]], [[0.0015136059999999998, 0.001522984, 0.001682003, 0.008274776, 0.002962329, 0.0038887, 0.0051478389999999995, 0.124941599, 0.03834592, 0.003115126], [0.0017308019999999998, 0.00284123, 0.006090379999999999, 0.011636356, 0.001954774, 0.017210396, 0.0031669449999999996, 0.0030475100000000002, 0.002737034, 0.0035592609999999998], [0.003043372, 0.016784805, 0.019664129000000002, 0.01006351, 0.01136131, 0.024623692, 0.002439496, 0.003555404, 0.00345213, 0.003662714], [0.019280899, 0.024736663000000002, 0.035284156, 0.043082425, 0.0032964570000000005, 0.031056285, 0.0035815540000000002, 0.003637489, 0.00330293, 0.003592695], [0.014394078999999999, 0.014527202, 0.080629323, 0.036650322, 0.002886918, 0.0030291140000000003, 0.002652948, 0.003463021, 0.003310153, 0.003483579], [0.031711728, 0.00912208, 0.002359396, 0.022938846000000002, 0.003241764, 0.002790205, 0.0043382659999999995, 0.003641385, 0.003141285, nan], [0.049836443, 0.071025118, 0.016307316999999998, 0.011022664, 0.0036095560000000003, 0.0029104219999999997, 0.0040801380000000005, 0.002795156, 0.0039314599999999995, nan], [0.031242587, 0.003043053, 0.051419719999999995, 0.002793882, 0.003205035, 0.003604848, 0.0036682660000000003, 0.0032497999999999997, 0.0034641339999999994, nan], [0.023518787, 0.002802037, 0.002582374, 0.0031874309999999993, 0.003000109, 0.002491033, 0.003599858, 0.003423217, 0.00421279, nan], [0.012390113000000001, 0.003620267, 0.0035406360000000002, 0.0029347559999999997, 0.002705216, 0.0033838719999999996, 0.0040380319999999996, 0.003404762, 0.003440628, nan], [0.002410523, 0.0026585380000000002, 0.003379032, 0.0035695510000000002, 0.003019645, 0.0038771190000000005, 0.003587734, 0.003660001, nan, nan], [0.02504284, 0.0029000319999999994, 0.0024978859999999995, 0.0029237969999999997, 0.0034266619999999996, 0.0034104730000000002, 0.003804753, 0.003785647, nan, nan], [0.003486786, 0.002927173, 0.003520337, 0.003175831, 0.0029503069999999997, 0.0035908809999999998, 0.0038760360000000002, 0.003789735, nan, nan], [0.003099998, 0.003400895, 0.0028864370000000004, 0.003174988, 0.00420155, 0.003654721, 0.003701519, 0.004147557, nan, nan], [0.0030547320000000005, 0.0033111049999999995, 0.004109016, 0.0030423299999999998, 0.003502927, 0.0037695060000000006, 0.003597193, 0.003550105, nan, nan], [0.00318592, 0.003429657, 0.003502125, 0.0032812690000000003, 0.0032791659999999996, 0.003836223, 0.003831292, nan, nan, nan], [0.003836603, 0.00380285, 0.003451949999999999, 0.0035710420000000004, 0.0033023489999999996, 0.004447441, 0.003788773, nan, nan, nan], [0.0030913009999999994, 0.0032595679999999996, 0.0032297830000000004, 0.003896977, 0.0036348920000000002, 0.0037536380000000005, 0.003266952, nan, nan, nan], [0.003100239, 0.0032018100000000002, 0.003608686, 0.003921302, 0.0033148720000000004, 0.0036884050000000005, 0.003863633, nan, nan, nan], [0.0034543850000000004, 0.0038549579999999995, 0.00335166, 0.004237467, 0.002790396, 0.0039046330000000002, 0.0026832849999999997, nan, nan, nan]], ([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.9735154725397268, 0.9566865186789388, 0.9491094147582697, 0.9738143729997091, 0.9744130042143287, 0.973349126443589, 0.9720176730486009, 0.9300699300699301, 0.0, 0.9726704841228527, 0.976840565577767, 0.9759919130654536, 0.9478351591027647, 0.9684431977559608, 0.9767008387698043, 0.9454649827784156, 0.9647335423197492, 0.9517766497461929, 0.9767576990122022, 0.9726609963547995, 0.9757207890743551, 0.9770525242223356, 0.9501743896362731, 0.9537465309898242, 0.6902533039647577, 0.9763903462749213, 0.9658886894075404, 0.9730266893810335, 0.9310216256524981, 0.9780594831789371, 0.0, 0.972936400541272, 0.8175965665236051, 0.9732367758186398, 0.9679075738125802, 0.9680462568472307, 0.9724593775819333, 0.8761811665037471, 0.9516351911561493, 0.8715144887916895, 0.8854564387365498, 0.7671798723244461, 0.9353796445880452, 0.9551368326603858, 0.9313138429024305, 0.9747849302877485, 0.9599325179249262, 0.9626865671641791, 0.9588194921070693, 0.9747340425531915, 0.9660835415922885, 0.9701351776171016, 0.9676823638042474, 0.975, 0.9745444801714899, 0.9703557312252964, 0.9460431654676259, 0.9717593900028241, 0.9762783342119136, 0.9656357388316151, 0.969751966122202, 0.964801778436458, 0.9490373725934315, 0.9637909319899244, 0.9511400651465798, 0.9412423055400112, 0.9761833879130694, 0.9441441441441442, 0.0, 0.9635435654393001, 0.9537712895377128, 0.8865814696485623, 0.9701907790143084, 0.9683639650497138, 0.9716949716949717, 0.9593147751605996, 0.9591503267973857, 0.9255079006772009, 0.9602313810556761, 0.9682539682539683, 0.9525566684238271, 0.6198830409356725, 0.9610757254069356, 0.8823529411764706, 0.963744232036915, 0.9704229517894114, 0.9698029243483789, 0.9637418419144308, 0.8452611218568665, 0.7158955223880598, 0.9017241379310345, 0.8993610223642172, 0.9701403404001194, 0.9662379421221865, 0.902896642527979, 0.9614760746147607, 0.9265785609397944, 0.940517844646606, 0.9606164383561644, 0.9620909816440543, 0.9289267945984364, 0.9600840336134454, 0.9317561419472248, 0.9306839186691312, 0.9543874891398784, 0.962241653418124, 0.924074074074074, 0.9587706146926537, 0.918195339613287, 0.9604496253122398, 0.8737904922170804, 0.9638597759306108, 0.9555160142348754, 0.9499165275459098, 0.9394306480920654, 0.9611992945326279, 0.9435665914221218, 0.9302473050095117, 0.9676099556767814, 0.9390986601705238, 0.9625836943678614, 0.9569672131147541, 0.0, 0.9441715051362215, 0.9454905847373637, 0.9589659417316373, 0.9582784365393061, 0.9240368963646229, 0.8004750593824228, 0.9445740956826137, 0.9075678798382438, 0.9249482401656315, 0.9559981472904122, 0.8007518796992481, 0.9352750809061489, 0.2682555780933063, 0.936268829663963, 0.960456942003515, 0.8981571290009699, 0.9519923187710033, 0.9440089585666294, 0.0, 0.9539058709364386, 0.9455394190871369, 0.949708840656432, 0.9502982107355865, 0.9380387931034483, 0.8991257565568258, 0.9459167117360735, 0.9397186872069658, 0.8792535675082327, 0.9025641025641026, 0.9073461283917935, 0.9524714828897338, 0.9579596412556054, 0.9551166965888689, 0.8164556962025317, 0.9459193706981318, 0.8587404355503238, 0.9543170397441755, 0.9448529411764706, 0.8426707597851113, 0.8146604046242775, 0.9425861208187718, 0.9419354838709677, 0.933920704845815, 0.9350086655112652, 0.92003046458492, 0.9254287844891872, 0.9234889058913542, 0.9541284403669725, 0.9292149292149292, 0.9358974358974359, 0.9657198824681684, 0.9436475409836066, 0.90625, 0.9262422360248447, 0.9334875650665124, 0.9160671462829736, 1.0, 0.6817275747508306, 0.9252669039145908, 0.9464184997179921, 0.9423177766124803, 0.9172185430463576, 0.953117674636662, 0.9464094319399786, 0.9292929292929293, 0.9315494710640946, 0.9228187919463087, 0.7514204003813155, 0.5305147058823529, 0.6958976237973103, 0.3378747948883696, 0.8480929678188319, 0.8957345971563981, 0.9041811846689896, 0.9022801302931596, 0.9380657455931396, 0.7460256410256411, 0.8548927294398092, 0.8948475289169295, 0.54, 0.8796498905908097, 0.8779134295227525, 0.9182114255205552, 0.8539226519337015, 0.9555837563451777, 0.5519713261648745, 0.5867114093959731, 0.0, 0.9210734017363852, 0.8556227140986834, 0.6638405733764531, 0.9265426052889324, 0.7650753768844221, 0.5793526785714286, 0.8300293139822402, 0.6636157337367625, 0.8217025862068966, 0.5734962406015038, 0.105, 0.01, 0.42518128548808404, 0.6145795039249184, 0.075, 0.4355246913580247, 0.5470706253791001, 0.6602875196786799, 0.315, 0.3355291294910039, 0.4522122614227877, 0.5937167573697864, 0.105, 0.115, 0.6467439785905441, 0.2482960077896787, 0.4821653543307087, 0.35514018691588783, 0.045000000000000005, 0.7615629319006765, 0.3612423129196776, 0.19646031580048104, 0.21194829393268894, 0.31628771122710786, 0.750858959489158, 0.32, 0.24114173228346455, 0.5934026305966608, 0.338062015503876, 0.4915625946283229, 0.07, 0.02, 0.05073800738007381, 0.03982300884955752, 0.055, 0.025, 0.014218009478672985, 0.21853146853146851, 0.065, 0.2947643360977481, 0.015, 0.0697211155378486, 0.0, 0.05876591576885406, 0.0776255707762557, 0.05, 0.08720930232558138, 0.2745651838685853, 0.11916583912611718, 0.05639097744360902, 0.13270142180094788, 0.9350012263919548, 0.958586416344561, 0.9684897591717309, 0.9716681776971895, 0.9745934959349594, 0.9723695844385499, 0.9672719594594594, 0.7983392645314353, 0.9831838565022422, 0.9646188485043422, 0.9425691514299109, 0.9592391304347826, 0.9797958315610379, 0.9675098630772802, 0.961181798576666, 0.9804856895056374, 0.9697775628626693, 0.933920704845815, 0.9788512911843277, 0.9677659310686834, 0.9285384970032273, 0.9800399201596807, 0.9363131079203335, 0.9343268502147006, 0.9540461779869984, 0.9327601410934744, 0.977079240340537, 0.9461749885478699, 0.9392835458409229, 0.976321036889332, 0.9538538288850087, 0.9792213473315835, 0.564618966977138, 0.011818620356970575, 0.9568384138117075, 0.9447278911564626, 0.9738288406176394, 0.9733570159857904, 0.9757340451346761, 0.974974974974975, 0.9775336994508238, 0.9767498776309349, 0.0, 0.9235971943887775, 0.9430441105186621, 0.9722958323295592, 0.9431818181818182, 0.9642194540153723, 0.975018735948039, 0.9607396149949341, 0.9461206896551724, 0.9739854318418314, 0.7981132075471699, 0.9685838569357177, 0.9754601226993865, 0.9756157034869544, 0.9285125103106956, 0.9712643678160919, 0.9583843329253366, 0.9732815156667476, 0.9433819828887771, 0.977459559798462, 0.9613501674826076, 0.969824193125164, 0.9123346662565365, 0.9120326735783851, 0.8788115715402658, 0.961890243902439, 0.968895800933126, 0.9719730941704036, 0.9709746410021387, 0.9652432969215492, 0.9346219442865265, 0.8987854251012146, 0.9573820395738204, 0.9576629974597799, 0.9720588235294118, 0.9747876857749469, 0.9741453605285837, 0.9674267100977199, 0.9631050767414404, 0.9568965517241379, 0.9122023809523809, 0.9705622608183692, 0.9759229534510433, 0.9448428019856592, 0.3302477926516662, 0.9700956937799043, 0.9736842105263158, 0.9727011494252874, 0.9653350940904589, 0.9753154141524959, 0.9568841621155505, 0.9648975007020499, 0.9729583558680368, 0.966044142614601, 0.9702380952380952, 0.9389473684210526, 0.9646583641871423, 0.4046659953004364, 0.9602922490470139, 0.9698986058301647, 0.0, 0.917864476386037, 0.9616858237547893, 0.9470809149880506, 0.95724973281083, 0.9677523379554982, 0.9619508820477344, 0.9586503473370823, 0.9094460227272727, 0.9117276166456494, 0.9707773232028054, 0.9638678596008259, 0.9588073754413495, 0.9667931688804554, 0.9486125385405961, 0.9636929460580913, 0.9474009900990099, 0.965034965034965, 0.9668141592920354, 0.9626068376068376, 0.1839557399723375, 0.94679186228482, 0.9646142958244869, 0.9437052200614124, 0.9061406299713649, 0.9495677233429395, 0.9113197704747, 0.9602423324498296, 0.9471992653810836, 0.9272953884503531, 0.9643393393393394, 0.9595572584078331, 0.9434762129062647, 0.0, 0.9363057324840764, 0.9152910512597741, 0.9689681923972071, 0.9585921325051759, 0.8722826086956522, 0.9658273381294964, 0.9537648612945839, 0.9650227352221056, 0.9164467897977133, 0.9610743479953289, 0.703509843856076, 0.8010323574730355, 0.9696739954510993, 0.9439601494396015, 0.9488553336580614, 0.9129930394431555, 0.5726275016139445, 0.9488708990200255, 0.9610136452241715, 0.9649859943977591, 0.9366359447004609, 0.922509225092251, 1.0, 0.9243379571248423, 0.8850574712643678, 0.45956639566395663, 0.6972147651006712, 0.6731537102473498, 0.9325567136725935, 0.015, 0.7807378904481665, 0.9252275682704811, 0.46595419847328245, 0.9341692789968652, 0.9580536912751678, 0.9369171695008228, 0.781517094017094, 0.9445324881141046, 0.8973701090442592, 0.9609375, 0.93993993993994, 0.44549105094079855, 0.9403710247349824, 0.956140350877193, 0.946319018404908, 0.9497991967871486, 0.9611769513690233, 0.9511480214948705, 0.8999884792626728, 0.9250340754202635, 0.8922742200328406, 0.5736054421768708, 0.9518348623853211, 0.04, 0.3448291925465839, 0.9347705252422234, 0.684901185770751, 0.9424986931521171, 0.38077847188851516, 0.895, 0.8229123711340206, 0.9442586399108138, 0.8548710601719198, 0.8432, 0.6283280757097791, 0.7651662635439404, 0.5608699531121442, 0.9022801302931596, 0.9560527367159408, 0.8373027927486526, 0.571665485471297, 0.2576633165829146, 0.48429539295392954, 0.40475667120939024, 0.7646586345381526, 0.271727078891258, 0.6174812030075189, 0.7479452054794521, 0.6438653603034135, 0.8873764906303236, 0.6742052023121388, 0.8324255667031912, 0.2766236003445306, 0.18000000000000002, 0.025, 0.03, 0.16421641791044778, 0.23690387941685198, 0.6731279872543813, 0.4252647678668291, 0.0499001996007984, 0.04, 0.3502893245158243, 0.2713368983957219, 0.11, 0.18316901408450706, 0.11499999999999999, 0.048828125, 0.03, 0.20069686411149826, 0.4798148302610727, 0.3055229142185664, 0.5424521195173985, 0.41531697341513296, 0.40346323362890446, 0.25271739130434784, 0.045, 0.5429781879194631, 0.02, 0.4703908891480093, 0.060000000000000005, 0.3270662599208914, 0.03, 0.09460737937559129, 0.065, 0.17359066236208182, 0.03, 0.06756756756756757, 0.10858835143139191, 0.03710575139146567, 0.055, 0.10480349344978165, 0.015, 0.06298449612403101, 0.04488330341113106, 0.04633920296570899, 0.0044014084507042256, 0.004500450045004501, 0.10793255436308805, 0.9718492854049372, 0.21607014172471775, 0.9744671403197158, 0.9780268072950999, 0.9751007825468342, 0.945792462570986, 0.9671398527865405, 0.9566395663956639, 0.9617203332582752, 0.795, 0.9451342766387525, 0.9654406319427302, 0.972972972972973, 0.9691710436172642, 0.9679982343853454, 0.9692552949214301, 0.9746660525103639, 0.9657456040191824, 0.9421450823319982, 0.9732555312424022, 0.9520692073883563, 0.9262131065532766, 0.976910644193027, 0.9605003291639236, 0.9536546610169492, 0.9768250289687138, 0.9692089057318807, 0.9653243847874721, 0.9476148016049933, 0.9779298168174796, 0.959733407386837, 0.9703743603555077, 0.972, 0.971042471042471, 0.973159509202454, 0.9624648107600876, 0.951481772882245, 0.95603517186251, 0.9698976520168573, 0.9095022624434389, 0.9041095890410958, 0.8512669849430775, 0.9561623246492986, 0.9505806770447245, 0.9341486359360301, 0.2519751017476658, 0.9111405835543767, 0.04835164835164835, 0.9708609271523179, 0.9655787863335034, 0.9708576186511241, 0.9729510413849067, 0.8197820620284996, 0.9708994708994709, 0.9736578023080783, 0.376953125, 0.9487836107554417, 0.9384010484927916, 0.9714518760195758, 0.9353448275862069, 0.8976, 0.9478487614080835, 0.9616368286445013, 0.9333116460637606, 0.9565670604586518, 0.971238268240993, 0.9567474048442907, 0.9597747385358005, 0.9562142135399124, 0.9363057324840764, 0.9665246500304321, 0.9701046337817638, 0.9250681198910081, 0.9074916590840157, 0.12410841654778887, 0.9714983713355049, 0.9492242595204513, 0.9676624576532183, 0.9653121902874133, 0.9708051628764598, 0.9498229043683589, 0.8847926267281107, 0.8615569823434992, 0.9723936342968497, 0.9662596401028277, 0.9033613445378151, 0.9393442622950819, 0.9735138316656857, 0.9637812386816371, 0.9377408834865105, 0.94384765625, 0.945424013434089, 0.8779795686719637, 0.9721932568647897, 0.9620587264929066, 0.9629368160960113, 0.8652151373768792, 0.9578754578754579, 0.74, 0.9612903225806452, 0.9563227953410982, 0.004068190623789229, 0.961281239000352, 0.9440820130475303, 0.9415086728519564, 0.9076992345790185, 0.9603174603174603, 0.9680817108203, 0.9717613836921991, 0.8675877950489349, 0.9285714285714286, 0.9113233287858117, 0.9644212523719166, 0.962432915921288, 0.588495575221239, 0.9605055292259084, 0.9619047619047619, 0.8614179104477612, 0.8989637305699482, 0.9666444296197465, 0.9493321050207277, 0.910958904109589, 0.9619952494061758, 0.9695387293298521, 0.7647648902821317, 0.7767000911577029, 0.7076661514683154, 0.38039748953974895, 0.9419104991394148, 0.9545889101338432, 0.6199763733018311, 0.5467395264116576, 0.06, 0.6794809010773751, 0.9515905947441217, 0.2970972354623451, 0.9564183835182251, 0.034999999999999996, 0.9270699270699271, 0.015000000000000001, 0.9514759401536595, 0.6317816365366319, 0.9075, 0.906630891273915, 0.9522445081184336, 0.9542809642560266, 0.729590527873705, 0.9603505843071787, 0.9366041896361632, 0.08, 0.19151162790697673, 0.6769841033579423, 0.3551672104404568, 0.43577279752704795, 0.6702239469641847, 0.6586961972419557, 0.41231410701876303, 0.2532961309523809, 0.5854503641241855, 0.5997814207650274, 0.3505422446406053, 0.1607607776838546, 0.33228488792480115, 0.33150812064965196, 0.8734906114588349, 0.37757575757575756, 0.030000000000000002, 0.12594850948509487, 0.18321428571428572, 0.26628212450028554, 0.3015555555555555, 0.43342412451361867, 0.14, 0.2140205552998926, 0.21277460770328102, 0.8230587720531933, 0.44958715596330273, 0.37398373983739835, 0.6745253164556961, 0.015, 0.0, 0.075, 0.24276744853482782, 0.045000000000000005, 0.15000000000000002, 0.3015241882041087, 0.01, 0.005, 0.04, 0.035, 0.039999999999999994, 0.085, 0.25665768194070077, 0.15037974683544303, 0.01, 0.015000000000000001, 0.18886177114025218, 0.0, 0.32845564074479733, 0.024999999999999998, 0.16, 0.0, 0.024999999999999998, 0.030000000000000002, 0.30676804655965456, 0.005, 0.10500000000000001, 0.09, 0.3268806619644329, 0.09877704609595485, 0.08499999999999999, 0.05293551491819056, 0.08490566037735849, 0.015, 0.15355642474286543, 0.085, 0.034999999999999996, 0.2459665991902834, 0.07440476190476192, 0.03255813953488372, 0.042293233082706765, 0.06126295947219604, 0.9322459222082811, 0.976806640625, 0.9758787043418332, 0.9437751004016064, 0.9414634146341463, 0.9786830885836096, 0.967948717948718, 0.9764521193092621, 0.9652605459057072, 0.9771376314586191, 0.9769585253456221, 0.9224035017906884, 0.9547613662067406, 0.9446808510638298, 0.9420625724217845, 0.9778451492537313, 0.9721633085896076, 0.9511555873242793, 0.9239728243286962, 0.943426724137931, 0.9537156081254822, 0.9724077839093813, 0.9780701754385965, 0.9649824912456229, 0.9746063991874048, 0.9394325788970354, 0.9736842105263158, 0.9758103531688437, 0.9405045216563541, 0.9706916764361079, 0.9732490272373541, 0.940337224383917, 0.9607258938244854, 0.9688524590163935, 0.981651376146789, 0.8952590959206174, 0.9303849303849304, 0.9680974477958236, 0.9716752090639331, 0.9638360567907849, 0.9653425753132499, 0.967208814270724, 0.9519959575543203, 0.9706427688504327, 0.9679393762751385, 0.9312, 0.9747340425531915, 0.8602698650674663, 0.94441431670282, 0.9755725190839695, 0.9176125961186379, 0.9649957994959395, 0.9694189602446484, 0.9620462046204621, 0.9790794979079498, 0.94199207696661, 0.9751504054407534, 0.9744, 0.9724972497249725, 0.9628308058281296, 0.9655742219774167, 0.9720506031185643, 0.9303530936183998, 0.9644507643085674, 0.9567949725058916, 0.761079478054567, 0.7803824969400245, 0.9568965517241379, 0.967632552404439, 0.9624701467076083, 0.5403414264036419, 0.9616122840690979, 0.9733664772727273, 0.9696582561481955, 0.9669867947178872, 0.44727860220201054, 0.0, 0.949305294780323, 0.741015625, 0.3987640449438203, 0.9109707064905227, 0.9721244131455399, 0.6329333333333333, 0.9667497921862012, 0.9720181881776845, 0.9291459211030256, 0.9322493224932249, 0.9584717607973422, 0.9640780020526856, 0.9556737588652482, 0.9686912961803381, 0.017346053772766695, 0.810451320457233, 0.9685732243871779, 0.9691758598312784, 0.7033836003770028, 0.025, 0.8988621444201312, 0.9395373291272345, 0.30718166383701184, 0.9335748792270532, 0.9548440065681445, 0.35082767978290363, 0.9439834024896265, 0.9561595791319597, 0.8199476987447698, 0.6586005183265458, 0.7989043583535109, 0.637062937062937, 0.39881435257410297, 0.9475982532751092, 0.9722016308376575, 0.9682971014492754, 0.965034965034965, 0.9517019319227231, 0.5289064976228209, 0.6082926829268293, 0.7917989792169396, 0.7437619263970923, 0.8253049907578558, 0.9582463465553236, 0.5692817679558011, 0.8191298145506419, 0.5909480986639261, 0.5786350148367952, 0.5247796610169491, 0.7603547209181012, 0.5830214152700186, 0.42766520056182367, 0.7705393835616439, 0.5119793459552496, 0.8819753086419753, 0.5032578209277239, 0.009999999999999998, 0.8115447154471545, 0.21063157894736842, 0.8333568904593639, 0.5875177076625886, 0.335390625, 0.4190732662546847, 0.09689453125, 0.395, 0.699228650137741, 0.049999999999999996, 0.4789905362776026, 0.7783589059372915, 0.4616470588235294, 0.7769141427102656, 0.7047163120567376, 0.0, 0.18026077097505672, 0.2776441351888668, 0.10058823529411764, 0.12549861495844875, 0.125, 0.19319227230910763, 0.025, 0.02, 0.005, 0.21025751072961374, 0.085, 0.31001901140684407, 0.2341159056416406, 0.2750554323725055, 0.346315331010453, 0.09999999999999999, 0.01, 0.015, 0.06999999999999999, 0.10348071495766697, 0.1982891061452514, 0.14466346153846155, 0.015, 0.03, 0.08999999999999998, 0.30092065106815874, 0.045, 0.005, 0.34683778234086243, 0.015, 0.03, 0.105, 0.3160982478097622, 0.1672316384180791, 0.03, 0.004999999999999999, 0.085, 0.08, 0.030000000000000002, 0.0348605577689243, 0.009737098344693282, 0.05500000000000001, 0.03349282296650718, 0.07272049226179376, 0.03234750462107209, 0.05769230769230769, 0.0, 0.019157088122605363, 0.01833180568285976, 0.9451679232350926, 0.9512069851052902, 0.9728331177231565, 0.9653441682600382, 0.9705056179775281, 0.9577729121050985, 0.9703054298642534, 0.9748364368394564, 0.9749874937468734, 0.9747048903878583, 0.9549168779938011, 0.9667379425041578, 0.9536560247167868, 0.9386411254115534, 0.9303944315545244, 0.9657869012707723, 0.9819971195391263, 0.914833215046132, 0.9756157034869544, 0.9807247494217425, 0.9619660620245758, 0.9742810381108253, 0.972074130489972, 0.9591451917033312, 0.9367769706752758, 0.9647373540856031, 0.9311995213879749, 0.9686648501362398, 0.972244250594766, 0.975023786869648, 0.9715524248171228, 0.9461704091048908, 0.9599615631005766, 0.9619771863117871, 0.9707357859531772, 0.9511599511599511, 0.979050279329609, 0.9535165723524657, 0.004419391206313416, 0.9676793794440853, 0.5346593255333792, 0.9419376244193762, 0.9437265527303043, 0.9494485294117647, 0.9808743169398907, 0.9593991067803491, 0.964562942963213, 0.8531013229104029, 0.9668874172185431, 0.3570435684647303, 0.9699781659388647, 0.9685804055984004, 0.37091487669053297, 0.9722870478413069, 1.0, 0.9452221545952526, 0.9310360941024816, 0.7436147592245153, 0.9253403601229688, 0.9656419529837251, 0.33763856544014903, 0.6292862241256246, 0.9626998223801065, 0.5861114003123373, 0.7419637883008356, 0.747741577193632, 0.06103092783505155, 0.9542961608775137, 0.5888626172034244, 0.8485802469135802, 0.7945565552699229, 0.9481417458945549, 0.6301987883377509, 0.6578324761204996, 0.9700070571630205, 0.4400062266500623, 0.7677660695468914, 0.6350173731758166, 0.9521040027369141, 0.807536231884058, 0.8359903201787043, 0.36434272300469484, 0.9481621112158342, 0.3556369982547993, 0.38237205523964257, 0.9105754276827371, 0.34221982758620695, 0.7464371772805508, 0.588821913380737, 0.038420074349442376, 0.4068507333908542, 0.2776912181303116, 0.8075404376784016, 0.7025975975975975, 0.7726109570041608, 0.768094038623006, 0.41990310077519377, 0.3580859375, 0.863167119071537, 0.6988767688679245, 0.41086393088552914, 0.5410029069767441, 0.517783752640676, 0.43891050583657587, 0.44105046343975285, 0.37879341864716637, 0.44253389434315094, 0.6437303958177745, 0.3862548262548262, 0.3903481331987891, 0.3168411927877947, 0.31360927152317886, 0.2811336032388664, 0.5524137931034483, 0.4665189873417721, 0.05365853658536585, 0.3818921389396709, 0.049999999999999996, 0.7432995714842229, 0.025, 0.25970149253731345, 0.034999999999999996, 0.125, 0.4332963149924281, 0.3622332616433087, 0.18108886107634542, 0.024999999999999998, 0.02, 0.09, 0.06999999999999999, 0.060000000000000005, 0.020000000000000004, 0.2891699092088197, 0.030000000000000002, 0.27991568296795954, 0.025, 0.4640625, 0.03, 0.0817391304347826, 0.7519055944055943, 0.4885610640870617, 0.25600311041990664, 0.47893599334995846, 0.015, 0.009999999999999998, 0.025, 0.015, 0.02, 0.075, 0.075, 0.2957479627473807, 0.2320915841584158, 0.015, 0.055, 0.2299009900990099, 0.019999999999999997, 0.049999999999999996, 0.19787253983130274, 0.16748768472906403, 0.25097256857855366, 0.03, 0.0, 0.02, 0.08, 0.01, 0.085, 0.08828996282527882, 0.03381234150464919, 0.009389671361502348, 0.9706047032474804, 0.710512978986403, 0.9622703412073491, 0.5827870493991989, 0.8668131206740897, 0.9662775616083009, 0.9706572769953051, 0.4844933920704846, 0.9013320177602367, 0.9256823127137084, 0.8953460287970274, 0.8907118974644749, 0.7354669358909642, 0.9648123324396782, 0.9749812359269452, 0.9586939417781275, 0.9152573529411765, 0.9713384924047005, 0.9687304565353346, 0.33995196266813066, 0.7025, 0.9657086871325931, 0.8135477582846004, 0.9691539365452408, 0.8078127233738385, 0.6830748175182482, 0.6307687028140014, 0.9188673242125358, 0.9808673469387755, 0.6251478226128646, 0.33548821548821545, 0.7895732838589982, 0.7935833968012186, 0.5233426834969612, 0.6153340132585416, 0.9625407166123778, 0.9481431917029107, 0.4750118990956687, 0.10874999999999999, 0.9498997995991983, 0.5709510086455332, 0.3685453359425962, 0.6891460055096419, 0.7844131695667007, 0.6576546906187626, 0.7885933147632311, 0.6698395356777057, 0.5765056213916743, 0.840214899713467, 0.4629888739042481, 0.9738863287250384, 0.5279616506877867, 0.24670731707317076, 0.9290518191841235, 0.889807774455361, 0.9670482136663198, 0.8108306364617045, 0.9408450704225352, 0.46714285714285714, 0.959830866807611, 0.3005460750853242, 0.9733244414804935, 0.4004867256637168, 0.328430566967954, 0.3838668555240793, 0.4770978380359106, 0.12728197674418604, 0.1874074074074074, 0.6219222462203023, 0.1, 0.5542413620430645, 0.6300681302043906, 0.7977827447474295, 0.9436323479442745, 0.19003134796238247, 0.5078426395939086, 0.5709996492458786, 0.009999999999999998, 0.4845443037974684, 0.6961065573770492, 0.40598988720342283, 0.14260387811634348, 0.6468232484076434, 0.41497107969151675, 0.06334841628959276, 0.08331360946745563, 0.22888321995464853, 0.025, 0.20825581395348838, 0.18658536585365854, 0.025, 0.19281425891181989, 0.22937573616018847, 0.20717781402936378, 0.7189145658263305, 0.0, 0.20866946778711482, 0.30938084112149533, 0.0, 0.07999999999999999, 0.08821428571428572, 0.01, 0.015, 0.16794585987261146, 0.725, 0.04, 0.29579268292682925, 0.4757170099160946, 0.02, 0.030000000000000006, 0.01, 0.1752510849349039, 0.039999999999999994, 0.009999999999999998, 0.009066183136899365, 0.049999999999999996, 0.06999999999999999, 0.0, 0.01, 0.025, 0.024999999999999998, 0.030000000000000002, 0.004999999999999999, 0.03, 0.039999999999999994, 0.061611374407582936, 0.0, 0.009999999999999998, 0.0181653042688465, 0.04, 0.004725897920604915, 0.039999999999999994, 0.015, 0.034482758620689655, 0.035, 0.916045820958846, 0.7114236853356135, 0.5753512458721105, 0.3416157989228007, 0.9726775956284153, 0.4961876947040499, 0.5426599823061045, 0.8848183556405355, 0.39151168572849326, 0.22544802867383515, 0.7481302638966872, 0.4433528493364559, 0.7663507109004739, 0.5539468008626887, 0.35337726523887975, 0.5925609756097561, 0.555883476599809, 0.9339049660593068, 0.2061423220973783, 0.46765845070422535, 0.511588785046729, 0.42628465804066545, 0.9661213720316623, 0.577065981611682, 0.7806417112299465, 0.8883213773314205, 0.8042691415313225, 0.4948828125, 0.5025911458333333, 1.0, 0.6940740740740741, 0.317880496054115, 0.0, 0.5954956169925826, 0.3163551051051051, 0.46746789727126803, 0.3275, 0.968501326259947, 0.2967919075144509, 0.5527230046948357, 0.5516365688487584, 0.6641342377722241, 0.16499999999999998, 0.6708110367892977, 0.4311827956989247, 0.48710017574692444, 0.554185303514377, 0.8230837004405286, 0.5164610913787537, 0.32425851703406816, 0.5008657810352726, 0.29923423423423423, 0.9148034793814432, 0.5057849829351536, 0.11839357429718876, 0.015, 0.005, 0.2697244094488189, 0.01, 0.3983333333333333, 0.17387900355871885, 0.20953281423804226, 0.06, 0.035, 0.3377860696517413, 0.01, 0.3420722433460076, 0.2799226305609284, 0.24643979057591622, 0.43405156537753226, 0.08267506709557074, 0.0, 0.03, 0.015, 0.3229409647228222, 0.29805097451274365, 0.511304347826087, 0.049999999999999996, 0.18820754716981133, 0.08, 0.04, 0.025, 0.009999999999999998, 0.03839378238341969, 0.004999999999999999, 0.035, 0.003154639175257732, 0.06202247191011236, 0.20088442756375474, 0.13006711409395974, 0.3441560798548094, 0.31414012738853503, 0.15893258426966295, 0.005, 0.26159105851413544, 0.0, 0.3238847583643123, 0.0, 0.065, 0.0, 0.06, 0.029850746268656712, 0.04, 0.019999999999999997, 0.023741690408357073, 0.08, 0.004999999999999999, 0.05, 0.01488095238095238, 0.055, 0.015, 0.005, 0.02, 0.015, 0.1050228310502283, 0.5096358169061144, 0.06, 0.20217967599410896, 0.24741197183098593, 0.7024932091579356, 0.40753007217321574, 0.24619537275064265, 0.6997428666224286, 0.5096858638743456, 0.5499152542372882, 0.35177435554067626, 0.23113207547169812, 0.5020668316831682, 0.43095137420718815, 0.23463931718061673, 0.25414370078740156, 0.3418265259293254, 0.24447426586675086, 0.3968503937007874, 0.005, 0.37215474209650584, 0.6858291457286432, 0.5034146341463415, 0.5659697439875873, 0.38175771971496436, 0.3763327032136106, 0.07863013698630136, 0.13489583333333333, 0.5348971312662997, 0.7285950413223141, 0.20093406593406593, 0.5459595107549556, 0.21561215370866846, 0.4630424339471577, 0.23191550925925927, 0.43679892400806997, 0.015, 0.6756424747174301, 0.3746614327772326, 0.7424801236749117, 0.33556594948550045, 0.06000000000000001, 0.2333131067961165, 0.7453709394205443, 0.35105250709555347, 0.0, 0.46375, 0.3717171717171717, 0.04, 0.16686046511627908, 0.045000000000000005, 0.08499999999999999, 0.11461956521739129, 0.075, 0.030000000000000002, 0.2655806938159879, 0.025, 0.32381740523224767, 0.075, 0.19729041916167664, 0.18511848341232226, 0.02, 0.065, 0.1, 0.025, 0.025, 0.18393617021276598, 0.015, 0.02, 0.27676009892827697, 0.014258555133079848, 0.045, 0.055, 0.015, 0.045, 0.045, 0.2684182908545727, 0.02, 0.035, 0.005, 0.055, 0.039999999999999994, 0.02, 0.04113345521023766, 0.06445672191528544, 0.4578689397975492, 0.20783316378433367, 0.3263176895306859, 0.30080876158382475, 0.009999999999999998, 0.2274113475177305, 0.20265486725663714, 0.2, 0.23, 0.07, 0.3320930232558139, 0.025, 0.4991228540772532, 0.035, 0.09999999999999999, 0.2762560220233998, 0.5591549295774648, 0.12258169934640523, 0.215, 0.16454391891891892, 0.18508009153318075, 0.3977371048252912, 0.01, 0.45070005385029616, 0.4538512553292278, 0.0, 0.0, 0.3169834710743802, 0.19694371727748688, 0.185, 0.025, 0.005000000000000001, 0.15703703703703703, 0.035, 0.0, 0.363926812104152, 0.21352480417754569, 0.18896455223880596, 0.1953125, 0.034999999999999996, 0.09000000000000001, 0.3930584192439862, 0.41570302233902756, 0.2096177442189712, 0.05499999999999999, 0.06000000000000001, 0.115, 0.16753863134657837, 0.06999999999999999, 0.0, 0.065, 0.01, 0.05500000000000001, 0.075, 0.18313953488372092, 0.03, 0.005000000000000001, 0.055, 0.01, 0.01, 0.045, 0.03, 0.024999999999999998, 0.005, 0.09499999999999999, 0.024999999999999998, 0.2504025044722719, 0.16, 0.01, 0.005, 0.251680506993007, 0.005, 0.020000000000000004, 0.3021287779237845, 0.01, 0.01, 0.18878048780487805, 0.1, 0.015, 0.22879028491359177, 0.21338350391330527, 0.4386015325670498, 0.045000000000000005, 0.1901082543978349, 0.095, 0.30164499717354437, 0.17, 0.18310055865921787, 0.28045364891518737, 0.1620039292730845, 0.17639817629179333, 0.01, 0.025, 0.015, 0.42129753914988816, 0.02, 0.2623423423423423, 0.025, 0.02, 0.02, 0.4474113475177305, 0.06, 0.29051724137931034, 0.01, 0.045, 0.015, 0.01, 0.009999999999999998, 0.034999999999999996, 0.05, 0.030000000000000002, 0.055, 0.23138443935926775, 0.035, 0.0, 0.035, 0.01, 0.04, 0.11499999999999999, 0.005, 0.09000000000000001, 0.005, 0.01, 0.0, 0.06863849765258216, 0.004999999999999999, 0.08, 0.009999999999999998, 0.045, 0.030000000000000002, 0.04, 0.05500000000000001, 0.035, 0.0, 0.024999999999999998, 0.0, 0.015, 0.009999999999999998, 0.03, 0.22410564225690274, 0.05, 0.02, 0.09000000000000001, 0.10914893617021278, 0.1406528189910979, 0.065, 0.55378146101903, 0.08499999999999999, 0.004999999999999999, 0.2725323686214775, 0.09500000000000001, 0.01, 0.18488964346349746, 0.039999999999999994, 0.02, 0.049999999999999996, 0.015, 0.0, 0.015, 0.020000000000000004, 0.05, 0.02, 0.055, 0.05500000000000001, 0.009999999999999998, 0.005, 0.01, 0.034999999999999996, 0.0, 0.004999999999999999, 0.16624716553287983, 0.01, 0.004999999999999999, 0.005, 0.055, 0.009999999999999998, 0.060000000000000005, 0.04, 0.06, 0.03, 0.065, 0.009999999999999998, 0.039999999999999994, 0.030000000000000002, 0.02, 0.019999999999999997, 0.025, 0.024999999999999998, 0.045, 0.045, 0.0, 0.035, 0.015, 0.075, 0.01], [0.6656891495601173, 0.060000000000000005, 0.0, 0.0, 0.003, 0.005000000000000001, 0.0, 0.02446183953033268, 0.005000000000000001, 0.0, 0.01824817518248175, 0.01, 0.0, 0.00493583415597236, 0.005000000000000001, 0.004892367906066536, 1.0, 0.23, 0.275, 0.0, 0.005, 0.00463821892393321, 0.08, 0.12499999999999999, 0.019743336623889437, 0.0, 0.0, 0.03910068426197458, 0.0, 0.0, 0.01, 0.020000000000000004, 0.02, 0.020000000000000004, 0.004965243296921549, 0.07115749525616699, 0.004999999999999999, 0.005000000000000001, 0.015, 0.0, 0.015, 0.004999999999999999, 0.030000000000000002, 0.0, 0.05, 0.01, 0.01, 0.03, 0.025, 0.019627085377821395, 0.004708097928436911, 0.005000000000000001, 0.0, 0.005000000000000001, 0.0, 0.03710575139146567, 0.03, 0.0, 0.005000000000000001, 0.0, 0.0048828125, 0.004940711462450593, 0.0, 0.09250693802035152, 0.004999999999999999, 0.015, 0.004999999999999999, 0.004703668861712135, 0.024390243902439022, 0.0, 0.014619883040935672, 0.0, 0.215, 0.0, 0.7105831533477321, 0.0, 0.0, 0.04, 0.30635782747603835, 0.05, 0.18738574040219377, 0.07, 0.0, 0.010000000000000002, 0.07, 0.005, 0.015000000000000003, 0.05763688760806916, 0.005000000000000001, 0.03827751196172249, 0.00983284169124877, 0.009999999999999998, 0.005000000000000001, 0.004873294346978557, 0.005, 0.0, 0.03, 0.0, 0.15458937198067632, 0.018674136321195144, 0.01, 0.075, 0.0, 0.07500000000000001, 0.06529850746268656, 0.06, 0.0, 0.025, 0.02, 0.0, 0.03, 0.015, 0.009442870632672332, 0.0, 0.0, 0.014164305949008499, 0.09250693802035152, 0.02, 0.015, 0.0, 0.03, 0.0, 0.015, 0.0, 0.0, 0.023854961832061067, 0.00983284169124877, 0.020000000000000004, 0.004897159647404505, 0.005, 0.025, 0.009871668311944718, 0.004999999999999999, 0.004999999999999999, 0.0, 0.01932367149758454, 0.005, 0.004784688995215312, 0.0, 0.018921475875118263, 0.018921475875118256, 0.009999999999999998, 0.023540489642184557, 0.024826216484607744, 0.013736263736263736, 0.0046641791044776115, 0.005, 0.02334267040149393, 0.033206831119544596, 0.0, 0.015, 0.012733446519524618, 0.004344048653344918, 0.0, 0.0, 0.004625346901017576, 0.005, 0.02332089552238806, 0.09596928982725528, 0.0049900199600798395, 0.017953321364452424, 0.019157088122605363, 0.009861932938856016, 0.009182736455463728, 0.03, 0.022301516503122214, 0.0, 0.0, 0.017985611510791366, 0.13, 0.2, 0.024999999999999998, 0.0, 0.024366471734892786, 0.015, 0.020000000000000004, 0.035, 0.0, 0.024999999999999998, 0.030000000000000002, 0.00942507068803016, 0.06999999999999999, 0.0, 0.004887585532746823, 0.015000000000000001, 0.050597976080956765, 0.01, 0.005000000000000001, 0.005, 0.053554040895813046, 0.1, 0.009794319294809012, 0.115, 0.01, 0.005000000000000001, 0.0, 0.045000000000000005, 0.005, 0.045, 0.005, 0.03, 0.0, 0.02973240832507433, 0.025, 0.005000000000000001, 0.0, 0.0, 0.03, 0.0048828125, 0.02360717658168083, 0.01, 0.005, 0.01, 0.009451795841209828, 0.00949667616334283, 0.004664179104477612, 0.05752636625119846, 0.020000000000000004, 0.0, 0.08, 0.005, 0.03451676528599606, 0.009999999999999998, 0.005, 0.005, 0.005000000000000001, 0.02321262766945218, 0.0, 0.004999999999999999, 0.0, 0.023299161230195712, 0.005000000000000001, 0.009505703422053232, 0.018416206261510127, 0.0, 0.07, 0.05499999999999999, 0.024975024975024976, 0.0, 0.00922509225092251, 0.0, 0.025, 0.009999999999999998, 0.0, 0.01374885426214482, 0.02, 0.015, 0.038204393505253106, 0.005, 0.01, 0.0, 0.005000000000000001, 0.019249278152069296, 0.01949317738791423, 0.015, 0.014231499051233396, 0.004743833017077799, 0.009910802775024777, 0.004721435316336167, 0.0, 0.014548981571290007, 0.015, 0.014177693761814745, 0.005000000000000001, 0.01, 0.004999999999999999, 0.004708097928436912, 0.005000000000000001, 0.015, 0.014111006585136405, 0.014124293785310733, 0.0, 0.0, 0.08662175168431184, 0.01, 0.0049504950495049506, 0.0, 0.005000000000000001, 0.00909090909090909, 0.004766444232602479, 0.014436958614051972, 0.018298261665141813, 0.0, 0.004570383912248629, 0.015, 0.004424778761061947, 0.024366471734892786, 0.004629629629629629, 0.0, 0.004633920296570899, 0.0, 0.06999999999999999, 0.23, 0.09000000000000001, 0.0, 0.06999999999999999, 0.03, 0.035, 0.04000000000000001, 0.015000000000000003, 0.01, 0.0, 0.005, 0.030000000000000002, 0.004999999999999999, 0.034999999999999996, 0.015, 0.004999999999999999, 0.009930486593843098, 0.005, 0.005, 0.04, 0.034999999999999996, 0.020000000000000004, 0.004999999999999999, 0.105, 0.045, 0.060000000000000005, 0.005, 0.02906976744186046, 0.05223171889838556, 0.01, 0.01, 0.039999999999999994, 0.025, 0.025, 0.015, 0.005000000000000001, 0.005000000000000001, 0.0, 0.0, 0.0, 0.015, 0.015, 0.0, 0.0, 0.025, 0.0, 0.06737247353224254, 0.07, 0.005000000000000001, 0.014895729890764646, 0.03, 0.010000000000000002, 0.0, 0.0, 0.024999999999999998, 0.00947867298578199, 0.0, 0.01450676982591876, 0.034999999999999996, 0.01, 0.018484288354898334, 0.0, 0.009999999999999998, 0.025, 0.024999999999999998, 0.015, 0.004826254826254826, 0.0, 0.004830917874396135, 0.004553734061930783, 0.009606147934678193, 0.015, 0.0, 0.009633911368015413, 0.019175455417066153, 0.0099601593625498, 0.0, 0.0, 0.0, 0.0, 0.029325513196480937, 0.004999999999999999, 0.0, 0.0, 0.1, 0.075, 0.005, 0.004999999999999999, 0.015, 0.004868549172346641, 0.004766444232602478, 0.02, 0.0, 0.0, 0.0249003984063745, 0.005000000000000001, 0.0, 0.004999999999999999, 0.024999999999999998, 0.008880994671403197, 0.004868549172346641, 0.014520813165537268, 0.0, 0.015, 0.0, 0.009389671361502348, 0.004826254826254827, 0.019455252918287938, 0.027675276752767528, 0.004873294346978557, 0.020000000000000004, 0.015, 0.005, 0.005000000000000001, 0.01, 0.01, 0.004770992366412213, 0.009765625, 0.00949667616334283, 0.014124293785310733, 0.019665683382497544, 0.0, 0.004574565416285453, 0.024999999999999998, 0.0, 0.004366812227074236, 0.01, 0.005000000000000001, 0.015, 0.020000000000000004, 0.0244140625, 0.09803921568627451, 0.02, 0.03358925143953935, 0.0049751243781094535, 0.005000000000000001, 0.0, 0.004633920296570899, 0.020000000000000004, 0.004873294346978557, 0.0, 0.03, 0.03, 0.0, 0.025, 0.02, 0.024999999999999998, 0.004180602006688963, 0.09000000000000001, 0.0, 0.03, 0.014999999999999998, 0.04, 0.1411214953271028, 0.04, 0.235, 0.0, 0.02, 0.04, 0.07, 0.005, 0.06999999999999999, 0.015000000000000001, 0.085, 0.015, 0.095, 0.01, 0.005, 0.019999999999999997, 0.005, 0.015, 0.005, 0.055, 0.0, 0.0, 0.05, 0.005, 0.004999999999999999, 0.055, 0.0, 0.04, 0.005, 0.005000000000000001, 0.03, 0.025, 0.025, 0.00499001996007984, 0.005, 0.020000000000000004, 0.025, 0.014285714285714285, 0.005, 0.015, 0.0, 0.005000000000000001, 0.019999999999999997, 0.01, 0.00911577028258888, 0.060000000000000005, 0.015000000000000001, 0.055, 0.0, 0.0, 0.015, 0.013661202185792348, 0.005000000000000001, 0.004999999999999999, 0.0, 0.01, 0.0, 0.0, 0.01, 0.035, 0.005, 0.004940711462450593, 0.005000000000000001, 0.004999999999999999, 0.03283302063789869, 0.024999999999999998, 0.015, 0.005, 0.009999999999999998, 0.03394762366634335, 0.0, 0.0, 0.005, 0.02360717658168083, 0.025, 0.005, 0.025, 0.034999999999999996, 0.03, 0.014619883040935672, 0.015, 0.005000000000000001, 0.0, 0.009416195856873824, 0.020000000000000004, 0.0, 0.00484027105517909, 0.004775549188156638, 0.034999999999999996, 0.014044943820224717, 0.0, 0.0, 0.009319664492078284, 0.0, 0.005000000000000001, 0.014450867052023121, 0.009999999999999998, 0.005, 0.014191106906338694, 0.03, 0.0, 0.0, 0.010000000000000002, 0.020000000000000004, 0.004739336492890996, 0.009416195856873822, 0.014018691588785047, 0.005000000000000001, 0.0, 0.004999999999999999, 0.025, 0.005000000000000001, 0.034999999999999996, 0.024975024975024976, 0.03, 0.025, 0.004708097928436911, 0.0, 0.0, 0.005000000000000001, 0.03314393939393939, 0.005000000000000001, 0.015, 0.09174311926605505, 0.019398642095053344, 0.00940733772342427, 0.01, 0.0, 0.005000000000000001, 0.009861932938856016, 0.02365930599369085, 0.0, 0.024154589371980676, 0.034999999999999996, 0.005, 0.0, 0.01, 0.0, 0.0, 0.015, 0.005000000000000001, 0.0, 0.009999999999999998, 0.004574565416285453, 0.034999999999999996, 0.0, 0.005000000000000001, 0.0, 0.01, 0.0, 0.015000000000000001, 0.01, 0.015, 0.0, 0.02, 0.01, 0.055, 0.04, 0.02, 0.055, 0.06999999999999999, 0.005, 0.08, 0.01, 0.065, 0.005, 0.024999999999999998, 0.02, 0.005000000000000001, 0.005, 0.004999999999999999, 0.034999999999999996, 0.04499999999999999, 0.0, 0.0, 0.004999999999999999, 0.0, 0.009999999999999998, 0.0, 0.0, 0.005, 0.02, 0.01, 0.049999999999999996, 0.051210428305400374, 0.004999999999999999, 0.004844961240310078, 0.01, 0.009999999999999998, 0.035, 0.004999999999999999, 0.005, 0.045000000000000005, 0.005000000000000001, 0.005000000000000001, 0.105, 0.005, 0.030000000000000002, 0.025, 0.0, 0.0, 0.009999999999999998, 0.009999999999999998, 0.005000000000000001, 0.020000000000000004, 0.03, 0.030000000000000002, 0.004999999999999999, 0.004999999999999999, 0.005, 0.0, 0.005, 0.005, 0.0, 0.0, 0.0, 0.0, 0.024999999999999998, 0.01483679525222552, 0.005, 0.01, 0.0, 0.005000000000000001, 0.015, 0.0, 0.015, 0.004686035613870666, 0.0, 0.0, 0.005000000000000001, 0.020000000000000004, 0.105, 0.019999999999999997, 0.005, 0.004999999999999999, 0.004999999999999999, 0.020000000000000004, 0.005, 0.025, 0.015, 0.0, 0.005, 0.04, 0.02, 0.009746588693957114, 0.0, 0.05, 0.009999999999999998, 0.01, 0.0, 0.0, 0.024999999999999998, 0.013901760889712697, 0.005, 0.03, 0.01, 0.0, 0.022810218978102186, 0.01, 0.0, 0.00916590284142988, 0.015, 0.005000000000000001, 0.0, 0.0, 0.004999999999999999, 0.00984251968503937, 0.005000000000000001, 0.02, 0.005, 0.004743833017077799, 0.004999999999999999, 0.015, 0.0, 0.0, 0.0, 0.01, 0.004999999999999999, 0.005000000000000001, 0.005000000000000001, 0.01404494382022472, 0.01, 0.01932367149758454, 0.025, 0.004721435316336167, 0.0, 0.020000000000000004, 0.005000000000000001, 0.019999999999999997, 0.008802816901408451, 0.09442870632672333, 0.005, 0.004999999999999999, 0.0, 0.0, 0.005, 0.03, 0.032740879326473335, 0.005000000000000001, 0.025, 0.009541984732824428, 0.035, 0.0, 0.004816955684007707, 0.02, 0.020000000000000004, 0.034999999999999996, 0.004655493482309125, 0.01921229586935639, 0.014619883040935672, 0.0, 0.0, 0.0, 0.005, 0.0, 0.0, 0.015, 0.0, 0.075, 0.0, 0.0044286979627989375, 0.07500000000000001, 0.025, 0.045000000000000005, 0.025, 0.0, 0.0, 0.13999999999999999, 0.06999999999999999, 0.009999999999999998, 0.06, 0.0, 0.005000000000000001, 0.009999999999999998, 0.015000000000000001, 0.020000000000000004, 0.005000000000000001, 0.034999999999999996, 0.0, 0.0, 0.004999999999999999, 0.035, 0.035, 0.09499999999999999, 0.005, 0.04, 0.01, 0.024999999999999998, 0.015, 0.005, 0.015, 0.0, 0.004999999999999999, 0.025, 0.005000000000000001, 0.039999999999999994, 0.015, 0.035, 0.05, 0.009999999999999998, 0.004999999999999999, 0.0, 0.025, 0.005, 0.005, 0.0, 0.0, 0.005, 0.025, 0.0, 0.004999999999999999, 0.004911591355599214, 0.0, 0.015, 0.0, 0.005, 0.0048543689320388345, 0.0, 0.009999999999999998, 0.1, 0.060000000000000005, 0.004999999999999999, 0.03, 0.0, 0.005000000000000001, 0.025, 0.02, 0.015, 0.025, 0.025, 0.034999999999999996, 0.020000000000000004, 0.065, 0.02, 0.0, 0.005000000000000001, 0.015, 0.005000000000000001, 0.009999999999999998, 0.004999999999999999, 0.004999999999999999, 0.005000000000000001, 0.014807502467917077, 0.014285714285714285, 0.01, 0.0, 0.045, 0.004999999999999999, 0.005000000000000001, 0.034999999999999996, 0.053868756121449556, 0.004999999999999999, 0.025, 0.009999999999999998, 0.005, 0.0, 0.005, 0.039999999999999994, 0.005000000000000001, 0.034999999999999996, 0.01, 0.005000000000000001, 0.0, 0.004999999999999999, 0.005, 0.03, 0.039999999999999994, 0.0, 0.009813542688910697, 0.0, 0.015, 0.024999999999999998, 0.0, 0.004668534080298786, 0.0, 0.0, 0.0, 0.0, 0.01, 0.025, 0.015, 0.0, 0.0, 0.0, 0.03, 0.03, 0.014436958614051972, 0.03, 0.004995004995004996, 0.004999999999999999, 0.004604051565377532, 0.03, 0.015, 0.0, 0.0, 0.01, 0.0, 0.0, 0.004999999999999999, 0.0, 0.02, 0.019999999999999997, 0.02, 0.01, 0.0, 0.028957528957528955, 0.0, 0.0, 0.01, 0.005, 0.0, 0.0047169811320754715, 0.025, 0.015, 0.0, 0.015, 0.0, 0.0, 0.004655493482309124, 0.0, 0.005000000000000001, 0.0, 0.004999999999999999, 0.025, 0.0, 0.019999999999999997, 0.019999999999999997, 0.020000000000000004, 0.019980019980019983, 0.005000000000000001, 0.004999999999999999, 0.013799448022079117, 0.005000000000000001, 0.01874414245548266, 0.005, 0.005, 0.02332089552238806, 0.019512195121951216, 0.019267822736030827, 0.014822134387351778, 0.034999999999999996, 0.01440922190201729, 0.025, 0.009999999999999998, 0.01, 0.019743336623889437, 0.005000000000000001, 0.029585798816568046, 0.015, 0.005000000000000001, 0.004347826086956522, 0.009523809523809525, 0.004873294346978557, 0.0, 0.013901760889712697, 0.09560229445506692, 0.0, 0.013636363636363636, 0.0, 0.0, 0.0048828125, 0.0, 0.0, 0.005000000000000001, 0.01, 0.019999999999999997, 0.005000000000000001, 0.025, 0.05500000000000001, 0.0, 0.075, 0.039999999999999994, 0.005, 0.015, 0.02, 0.11, 0.0, 0.055, 0.005, 0.06999999999999999, 0.005, 0.015, 0.0, 0.0, 0.02, 0.01, 0.045, 0.009999999999999998, 0.0, 0.019999999999999997, 0.009999999999999998, 0.01, 0.019999999999999997, 0.04, 0.1, 0.1, 0.145, 0.04429133858267716, 0.04703668861712135, 0.015, 0.0, 0.12, 0.004999999999999999, 0.0, 0.075, 0.045, 0.004999999999999999, 0.004999999999999999, 0.0, 0.009999999999999998, 0.009999999999999998, 0.015, 0.055, 0.01, 0.01, 0.0, 0.025, 0.01, 0.01, 0.004999999999999999, 0.005, 0.005, 0.01, 0.0, 0.025, 0.005, 0.0, 0.0, 0.019999999999999997, 0.0, 0.05, 0.0048543689320388345, 0.0, 0.004999999999999999, 0.045, 0.035, 0.155, 0.005000000000000001, 0.1, 0.015, 0.009999999999999998, 0.019999999999999997, 0.0, 0.0, 0.01818181818181818, 0.035, 0.0, 0.034999999999999996, 0.0, 0.005000000000000001, 0.005, 0.01, 0.02, 0.01, 0.01, 0.03, 0.005, 0.009999999999999998, 0.004999999999999999, 0.0, 0.0, 0.009999999999999998, 0.009861932938856016, 0.02, 0.020000000000000004, 0.03, 0.004999999999999999, 0.06999999999999999, 0.07, 0.015, 0.015000000000000003, 0.024999999999999998, 0.025, 0.03, 0.015, 0.0, 0.0, 0.020000000000000004, 0.029097963142580018, 0.0, 0.0, 0.0, 0.0, 0.015, 0.00927643784786642, 0.0, 0.005000000000000001, 0.0, 0.0, 0.004999999999999999, 0.0, 0.005000000000000001, 0.0, 0.025, 0.01, 0.009596928982725527, 0.009999999999999998, 0.025, 0.027932960893854747, 0.02, 0.0, 0.019999999999999997, 0.02806361085126286, 0.02, 0.0, 0.005000000000000001, 0.024679170779861797, 0.015, 0.014792899408284023, 0.009999999999999998, 0.020000000000000004, 0.009999999999999998, 0.004570383912248629, 0.005, 0.005000000000000001, 0.0, 0.024999999999999998, 0.015, 0.005000000000000001, 0.020000000000000004, 0.005, 0.01, 0.025, 0.025, 0.0, 0.005, 0.015, 0.014677103718199608, 0.005000000000000001, 0.005, 0.0, 0.004812319538017324, 0.025, 0.005000000000000001, 0.02, 0.020000000000000004, 0.0, 0.005, 0.0, 0.005000000000000001, 0.0, 0.005000000000000001, 0.0, 0.005000000000000001, 0.020000000000000004, 0.015, 0.0048543689320388345, 0.03, 0.0, 0.015, 0.0, 0.024999999999999998, 0.004557885141294439, 0.015, 0.035, 0.015, 0.004651162790697674, 0.0, 0.004859086491739554, 0.005, 0.022977941176470586, 0.019588638589618023, 0.0, 0.004826254826254827, 0.0045662100456621, 0.01460564751703992, 0.013876040703052728, 0.0, 0.005, 0.015, 0.015, 0.025, 0.0, 0.013876040703052728, 0.005, 0.004803073967339097, 0.005, 0.020000000000000004, 0.028708133971291863, 0.07, 0.0, 0.015, 0.004659832246039142, 0.02, 0.02, 0.105, 0.005000000000000001, 0.023299161230195712, 0.015, 0.02, 0.015000000000000001, 0.01, 0.045, 0.0, 0.02, 0.024999999999999998, 0.0, 0.019999999999999997, 0.03, 0.049999999999999996, 0.03, 0.024999999999999998, 0.06, 0.005, 0.005, 0.075, 0.03, 0.075, 0.019999999999999997, 0.024999999999999998, 0.005, 0.009999999999999998, 0.06999999999999999, 0.0, 0.004999999999999999, 0.034999999999999996, 0.02, 0.05500000000000001, 0.020000000000000004, 0.01, 0.01, 0.06, 0.0, 0.01, 0.01, 0.0, 0.005, 0.020000000000000004, 0.11, 0.025, 0.0, 0.0, 0.005000000000000001, 0.005000000000000001, 0.034999999999999996, 0.005, 0.005000000000000001, 0.024999999999999998, 0.034999999999999996, 0.0, 0.005, 0.004999999999999999, 0.005, 0.03, 0.055, 0.01, 0.0, 0.005000000000000001, 0.06999999999999999, 0.0, 0.01, 0.03, 0.004999999999999999, 0.005, 0.034999999999999996, 0.034999999999999996, 0.02, 0.004999999999999999, 0.0, 0.0, 0.0, 0.075, 0.005000000000000001, 0.0, 0.005, 0.024999999999999998, 0.0, 0.015, 0.005000000000000001, 0.030000000000000002, 0.0, 0.015, 0.015, 0.01483679525222552, 0.019999999999999997, 0.1, 0.005000000000000001, 0.025, 0.005000000000000001, 0.0047892720306513415, 0.005, 0.004608294930875576, 0.035, 0.02811621368322399, 0.01, 0.0, 0.004999999999999999, 0.0, 0.01, 0.005000000000000001, 0.005000000000000001, 0.005000000000000001, 0.0, 0.01, 0.004999999999999999, 0.0, 0.04, 0.009999999999999998, 0.025, 0.005000000000000001, 0.023786869647954328, 0.0, 0.015, 0.0, 0.009398496240601503, 0.015, 0.015, 0.009871668311944718, 0.0, 0.0, 0.025, 0.005000000000000001, 0.0, 0.025, 0.0, 0.0, 0.015, 0.025, 0.0, 0.0, 0.0, 0.034999999999999996, 0.015, 0.01, 0.025, 0.0, 0.01, 0.0, 0.01, 0.0, 0.005000000000000001, 0.004999999999999999, 0.02, 0.0, 0.005000000000000001, 0.00945179584120983, 0.034999999999999996, 0.004651162790697674, 0.015, 0.01372369624885636, 0.03, 0.115, 0.005, 0.00947867298578199, 0.0, 0.02, 0.01, 0.0, 0.0, 0.009999999999999998, 0.0, 0.02, 0.005000000000000001, 0.10486177311725453, 0.009999999999999998, 0.01996007984031936, 0.0, 0.015, 0.024999999999999998, 0.09871668311944719, 0.025, 0.005000000000000001, 0.025, 0.005000000000000001, 0.0, 0.0, 0.0, 0.019704433497536943, 0.004965243296921549, 0.004574565416285453, 0.02327746741154562, 0.020000000000000004, 0.07360157016683022, 0.02, 0.0, 0.02, 0.004770992366412214, 0.0047892720306513415, 0.004999999999999999, 0.0, 0.004739336492890996, 0.024999999999999998, 0.020000000000000004, 0.02321262766945218, 0.0, 0.015, 0.020000000000000004, 0.02, 0.028708133971291863, 0.013204225352112674, 0.025, 0.015, 0.0, 0.0, 0.004789272030651341, 0.034999999999999996, 0.0, 0.0, 0.0042194092827004225, 0.023674242424242424, 0.013736263736263736, 0.004906771344455349, 0.0, 0.03, 0.0, 0.0, 0.004389815627743635, 0.03, 0.024999999999999998, 0.004930966469428008, 0.02, 0.018921475875118263, 0.03, 0.010000000000000002, 0.0, 0.010000000000000002, 0.0, 0.0, 0.005000000000000001, 0.009523809523809523, 0.004999999999999999, 0.004935834155972359, 0.024999999999999998, 0.028901734104046242, 0.0, 0.0, 0.005000000000000001, 0.0, 0.0, 0.024999999999999998, 0.025, 0.020000000000000004, 0.005000000000000001, 0.065, 0.0, 0.0, 0.025, 0.0, 0.004999999999999999, 0.005, 0.01, 0.01, 0.015, 0.005, 0.045, 0.0, 0.005, 0.015, 0.04, 0.035, 0.0, 0.005, 0.065, 0.0, 0.0, 0.015, 0.015, 0.03, 0.004999999999999999, 0.034999999999999996, 0.009999999999999998, 0.015, 0.0, 0.025, 0.11, 0.049999999999999996, 0.0, 0.075, 0.01, 0.019999999999999997, 0.0, 0.039999999999999994, 0.05, 0.03, 0.01, 0.0, 0.0, 0.015, 0.019999999999999997, 0.005, 0.005, 0.009999999999999998, 0.085, 0.0, 0.01, 0.0, 0.03, 0.03, 0.085, 0.004999999999999999, 0.01, 0.005000000000000001, 0.015, 0.02, 0.01, 0.04, 0.025, 0.025, 0.0, 0.005, 0.005000000000000001, 0.005000000000000001, 0.020000000000000004, 0.015, 0.0, 0.009999999999999998, 0.0, 0.0, 0.004999999999999999, 0.039999999999999994, 0.015, 0.05, 0.015000000000000003, 0.0, 0.01, 0.0, 0.0, 0.025, 0.08, 0.03, 0.0, 0.02, 0.015, 0.01, 0.020000000000000004, 0.004999999999999999, 0.004999999999999999, 0.0, 0.0, 0.0, 0.015000000000000001, 0.02, 0.02, 0.0, 0.004999999999999999, 0.0, 0.005000000000000001, 0.019999999999999997, 0.005, 0.01, 0.0, 0.025, 0.0, 0.03, 0.0, 0.07, 0.0, 0.01, 0.03, 0.025, 0.015, 0.024999999999999998, 0.005, 0.004999999999999999, 0.004999999999999999, 0.03, 0.01, 0.015, 0.01, 0.020000000000000004, 0.005000000000000001, 0.019999999999999997, 0.005000000000000001, 0.024999999999999998, 0.0, 0.0, 0.04, 0.015, 0.024999999999999998, 0.005000000000000001, 0.015, 0.0, 0.01, 0.0, 0.11, 0.005, 0.020000000000000004, 0.0, 0.020000000000000004, 0.025, 0.0, 0.025, 0.02, 0.01, 0.024999999999999998, 0.024999999999999998, 0.005000000000000001, 0.015, 0.02, 0.0, 0.0, 0.024999999999999998, 0.00922509225092251, 0.005, 0.014822134387351778, 0.014492753623188404, 0.01, 0.0, 0.020000000000000004, 0.004999999999999999, 0.020000000000000004, 0.0, 0.005, 0.0, 0.095, 0.0049800796812749, 0.01, 0.0, 0.024999999999999998, 0.027700831024930747, 0.00484027105517909, 0.004955401387512388, 0.034999999999999996, 0.025, 0.009920634920634922, 0.015, 0.0, 0.0, 0.005, 0.005, 0.004816955684007707, 0.005, 0.004999999999999999, 0.0, 0.020000000000000004, 0.009999999999999998, 0.005000000000000001, 0.0, 0.005, 0.0, 0.005, 0.0, 0.005, 0.015, 0.005000000000000001, 0.005, 0.0, 0.0, 0.02, 0.005000000000000001, 0.019999999999999997, 0.0, 0.0, 0.005000000000000001, 0.08804448563484708, 0.015, 0.0, 0.027002700270027002, 0.01, 0.005, 0.005000000000000001, 0.01988071570576541, 0.0, 0.024999999999999998, 0.0, 0.01, 0.034999999999999996, 0.0, 0.04, 0.0, 0.004999999999999999, 0.005, 0.015, 0.015, 0.035, 0.039999999999999994, 0.0, 0.005, 0.004999999999999999, 0.01, 0.025, 0.005, 0.020000000000000004, 0.055, 0.0, 0.04, 0.020000000000000004, 0.024999999999999998, 0.0, 0.004999999999999999, 0.035, 0.02, 0.01, 0.08, 0.01, 0.005000000000000001, 0.02, 0.03, 0.0, 0.025, 0.065, 0.004999999999999999, 0.015, 0.004999999999999999, 0.005, 0.01, 0.0, 0.004999999999999999, 0.004999999999999999, 0.005000000000000001, 0.0, 0.004999999999999999, 0.005000000000000001, 0.0, 0.04, 0.0, 0.0, 0.005000000000000001, 0.025, 0.004999999999999999, 0.005, 0.04, 0.004999999999999999, 0.045, 0.010000000000000002, 0.03, 0.005000000000000001, 0.0, 0.015, 0.025, 0.004955401387512388, 0.0, 0.015, 0.02, 0.005, 0.0, 0.0, 0.034999999999999996, 0.005, 0.005000000000000001, 0.004999999999999999, 0.025, 0.025, 0.0, 0.0, 0.020000000000000004, 0.01, 0.004999999999999999, 0.03, 0.005000000000000001, 0.004999999999999999, 0.0, 0.025, 0.004999999999999999, 0.0, 0.015, 0.015, 0.03, 0.004999999999999999, 0.015, 0.0, 0.0, 0.060000000000000005, 0.015, 0.0, 0.015, 0.01, 0.005000000000000001, 0.0, 0.005000000000000001, 0.0, 0.024999999999999998, 0.01, 0.004999999999999999, 0.020000000000000004, 0.009999999999999998, 0.04, 0.009999999999999998, 0.01, 0.009999999999999998, 0.004999999999999999, 0.004999999999999999, 0.020000000000000004, 0.02, 0.004999999999999999, 0.0, 0.005, 0.024999999999999998, 0.004999999999999999, 0.02, 0.004999999999999999, 0.014691478942213515, 0.01, 0.0, 0.015, 0.0, 0.034999999999999996, 0.02976190476190476, 0.004830917874396135, 0.009999999999999998, 0.015, 0.005000000000000001, 0.0, 0.004921259842519685, 0.01, 0.0, 0.034999999999999996, 0.020000000000000004, 0.004999999999999999, 0.02, 0.0, 0.01, 0.019999999999999997, 0.015, 0.005000000000000001, 0.0, 0.0048216007714561235, 0.005, 0.009999999999999998, 0.014749262536873156, 0.004999999999999999, 0.025, 0.005000000000000001, 0.005, 0.005000000000000001, 0.010000000000000002, 0.01, 0.015, 0.020000000000000004, 0.01, 0.11, 0.024999999999999998, 0.0, 0.02, 0.005000000000000001, 0.0048543689320388345, 0.095, 0.004999999999999999, 0.0, 0.005000000000000001, 0.029154518950437316, 0.004703668861712135, 0.009999999999999998, 0.014312977099236639, 0.020000000000000004, 0.005000000000000001, 0.005, 0.025, 0.03, 0.010000000000000002, 0.0, 0.009115770282588878, 0.014423076923076922, 0.00937207122774133, 0.024975024975024976, 0.0, 0.028544243577545196, 0.0, 0.029527559055118106, 0.0, 0.020000000000000004, 0.005000000000000001, 0.020000000000000004, 0.02446183953033268, 0.01, 0.004999999999999999, 0.0, 0.03, 0.0, 0.004686035613870665, 0.014563106796116504, 0.013837638376383764, 0.004803073967339097, 0.005, 0.0, 0.08490566037735849, 0.0, 0.020000000000000004, 0.005, 0.08, 0.025, 0.005, 0.004734848484848485, 0.005000000000000001, 0.004999999999999999, 0.005, 0.020000000000000004, 0.015, 0.01450676982591876, 0.005000000000000001, 0.009999999999999998, 0.015, 0.0, 0.014044943820224719, 0.009319664492078284, 0.015, 0.0, 0.005, 0.020000000000000004, 0.004935834155972359, 0.005, 0.0, 0.04, 0.025, 0.025, 0.004999999999999999, 0.03, 0.01, 0.019999999999999997, 0.005000000000000001, 0.03, 0.075, 0.005, 0.0, 0.034999999999999996, 0.005000000000000001, 0.004999999999999999, 0.04, 0.03, 0.0, 0.035, 0.05, 0.015, 0.0, 0.035, 0.0, 0.005, 0.039999999999999994, 0.015, 0.004999999999999999, 0.005000000000000001, 0.005, 0.005000000000000001, 0.0, 0.0, 0.004999999999999999, 0.11, 0.005, 0.01, 0.025, 0.004999999999999999, 0.009999999999999998, 0.009999999999999998, 0.004999999999999999, 0.005000000000000001, 0.015, 0.005000000000000001, 0.015, 0.085, 0.004999999999999999, 0.005000000000000001, 0.004999999999999999, 0.0, 0.005000000000000001, 0.1, 0.025, 0.015, 0.004999999999999999, 0.02, 0.0, 0.1, 0.019999999999999997, 0.02, 0.0, 0.0, 0.035, 0.02, 0.03, 0.03, 0.005, 0.0, 0.035, 0.0, 0.020000000000000004, 0.115, 0.01, 0.004999999999999999, 0.0, 0.025, 0.0, 0.024999999999999998, 0.034999999999999996, 0.105, 0.0, 0.0, 0.01, 0.0, 0.0, 0.015000000000000003, 0.0, 0.03, 0.004901960784313725, 0.005000000000000001, 0.005, 0.005000000000000001, 0.03, 0.005000000000000001, 0.004999999999999999, 0.0, 0.0, 0.025, 0.005000000000000001, 0.05, 0.025, 0.015, 0.03, 0.03, 0.009999999999999998, 0.01, 0.0, 0.04, 0.005000000000000001, 0.0, 0.0, 0.004999999999999999, 0.005000000000000001, 0.0, 0.004999999999999999, 0.020000000000000004, 0.03, 0.005000000000000001, 0.095, 0.0, 0.005000000000000001, 0.0, 0.015, 0.005000000000000001, 0.005, 0.01881467544684854, 0.015, 0.0, 0.005000000000000001, 0.0, 0.024999999999999998, 0.015, 0.020000000000000004, 0.0, 0.005000000000000001, 0.024999999999999998, 0.005000000000000001, 0.0, 0.009999999999999998, 0.025, 0.019999999999999997, 0.005, 0.005, 0.005000000000000001, 0.0, 0.0, 0.0, 0.0, 0.024999999999999998, 0.005000000000000001, 0.00999000999000999, 0.020000000000000004, 0.0, 0.009999999999999998, 0.0, 0.005, 0.018433179723502304, 0.01, 0.03, 0.0, 0.009999999999999998, 0.004999999999999999, 0.03, 0.0, 0.02, 0.024999999999999998, 0.023255813953488372, 0.004999999999999999, 0.04, 0.060000000000000005, 0.005000000000000001, 0.03, 0.024999999999999998, 0.04, 0.005000000000000001, 0.0, 0.02, 0.0, 0.0, 0.009999999999999998, 0.014734774066797641, 0.0, 0.0, 0.0, 0.03, 0.004999999999999999, 0.005, 0.03, 0.0, 0.03, 0.004999999999999999, 0.005, 0.005, 0.009823182711198428, 0.004999999999999999, 0.0, 0.025, 0.0, 0.03, 0.015, 0.005000000000000001, 0.005000000000000001, 0.005, 0.0, 0.005, 0.005000000000000001, 0.024975024975024976, 0.0, 0.03, 0.00463821892393321, 0.004830917874396135, 0.01, 0.0, 0.015, 0.020000000000000004, 0.02, 0.0, 0.024999999999999998, 0.0, 0.0, 0.03475670307845084, 0.0047984644913627635, 0.034999999999999996, 0.005000000000000001, 0.024999999999999998, 0.015, 0.0, 0.005000000000000001, 0.004807692307692308, 0.023607176581680833, 0.0, 0.0, 0.02, 0.00903342366757001, 0.01, 0.0, 0.019685039370078743, 0.004999999999999999, 0.005000000000000001, 0.11, 0.0, 0.024999999999999998, 0.005000000000000001, 0.0, 0.05, 0.009999999999999998, 0.005, 0.025, 0.015, 0.020000000000000004, 0.03, 0.0, 0.01, 0.0, 0.004999999999999999, 0.025, 0.02, 0.020000000000000004, 0.024999999999999998, 0.024999999999999998, 0.005, 0.03, 0.025, 0.04, 0.004999999999999999, 0.0, 0.024999999999999998, 0.0, 0.025, 0.005000000000000001, 0.004999999999999999, 0.010000000000000002, 0.03, 0.005000000000000001, 0.035, 0.004999999999999999, 0.019999999999999997, 0.005000000000000001, 0.020000000000000004, 0.0, 0.004999999999999999, 0.025, 0.03, 0.02, 0.0, 0.0, 0.015, 0.0, 0.0, 0.005000000000000001, 0.005, 0.004999999999999999, 0.03, 0.005, 0.0, 0.02, 0.01, 0.005, 0.0, 0.005000000000000001, 0.019999999999999997, 0.0, 0.009999999999999998, 0.0, 0.019999999999999997, 0.0, 0.0, 0.005, 0.015, 0.0, 0.005, 0.0, 0.01, 0.020000000000000004, 0.03, 0.01, 0.019999999999999997, 0.009999999999999998, 0.005, 0.0, 0.005, 0.0, 0.020000000000000004, 0.005, 0.009999999999999998, 0.02, 0.034999999999999996, 0.0, 0.004844961240310078, 0.020000000000000004, 0.020000000000000004, 0.0, 0.009999999999999998, 0.0, 0.035, 0.004999999999999999, 0.025, 0.0, 0.004999999999999999, 0.04, 0.004999999999999999, 0.01, 0.004999999999999999, 0.0, 0.025, 0.004999999999999999, 0.005000000000000001, 0.0, 0.004930966469428008, 0.005, 0.015, 0.025, 0.025, 0.004999999999999999, 0.029880478087649404, 0.03, 0.0, 0.009999999999999998, 0.03, 0.0, 0.02, 0.04, 0.020000000000000004, 0.009999999999999998, 0.020000000000000004, 0.005, 0.023105360443622918, 0.01, 0.005000000000000001, 0.0, 0.0, 0.005000000000000001, 0.024999999999999998, 0.004999999999999999, 0.025, 0.0, 0.0, 0.015, 0.03, 0.005, 0.015, 0.0, 0.0, 0.020000000000000004, 0.0, 0.005, 0.0, 0.004999999999999999, 0.025, 0.028248587570621465, 0.0, 0.013786764705882353, 0.020000000000000004, 0.004999999999999999, 0.005000000000000001, 0.01, 0.03, 0.004999999999999999, 0.009999999999999998, 0.020000000000000004, 0.0, 0.0, 0.02387774594078319, 0.0, 0.020000000000000004, 0.005000000000000001, 0.005000000000000001, 0.020000000000000004, 0.01, 0.015, 0.024999999999999998, 0.005000000000000001, 0.024999999999999998, 0.004999999999999999, 0.005000000000000001, 0.02, 0.005000000000000001, 0.025, 0.005000000000000001, 0.0, 0.009999999999999998, 0.004803073967339097, 0.00499001996007984, 0.03, 0.0, 0.03, 0.0, 0.0, 0.03, 0.005000000000000001, 0.034999999999999996, 0.024999999999999998, 0.0, 0.015, 0.015, 0.03202195791399817, 0.005000000000000001, 0.014436958614051972, 0.010000000000000002, 0.005000000000000001, 0.024999999999999998, 0.02, 0.02, 0.004739336492890996, 0.010000000000000002, 0.03, 0.02, 0.005, 0.004655493482309124, 0.0, 0.034999999999999996, 0.005, 0.0, 0.0, 0.004766444232602479, 0.004999999999999999, 0.005000000000000001, 0.0, 0.013927576601671309, 0.004557885141294439, 0.015, 0.0, 0.005000000000000001, 0.0, 0.024999999999999998, 0.019588638589618023, 0.020000000000000004, 0.010000000000000002, 0.025, 0.0, 0.0, 0.025, 0.003952569169960475, 0.020000000000000004, 0.019999999999999997, 0.013824884792626729, 0.0, 0.004757373929590866, 0.0, 0.02, 0.025, 0.004557885141294439, 0.015, 0.015, 0.024999999999999998, 0.0, 0.034999999999999996, 0.024999999999999998, 0.004999999999999999, 0.004999999999999999, 0.005000000000000001, 0.045, 0.005, 0.005, 0.03, 0.01, 0.005, 0.01, 0.034999999999999996, 0.01, 0.0, 0.0, 0.0, 0.005, 0.015, 0.024999999999999998, 0.015000000000000001, 0.0, 0.01, 0.0, 0.03, 0.0, 0.0, 0.0, 0.004999999999999999, 0.01, 0.015, 0.0, 0.03, 0.01, 0.04, 0.01, 0.0, 0.0, 0.03, 0.0, 0.009999999999999998, 0.0, 0.13, 0.0, 0.0, 0.024999999999999998, 0.0, 0.030000000000000002, 0.025, 0.03, 0.035, 0.005, 0.0, 0.004999999999999999, 0.020000000000000004, 0.01, 0.025, 0.034999999999999996, 0.0, 0.005, 0.005000000000000001, 0.024999999999999998, 0.034999999999999996, 0.005000000000000001, 0.020000000000000004, 0.02, 0.005000000000000001, 0.005000000000000001, 0.0, 0.0, 0.015, 0.005, 0.03, 0.024999999999999998, 0.0, 0.005000000000000001, 0.004999999999999999, 0.0, 0.0, 0.035, 0.005000000000000001, 0.0, 0.0, 0.0, 0.005, 0.019999999999999997, 0.02849002849002849, 0.01, 0.0, 0.0, 0.009999999999999998, 0.0, 0.0, 0.0, 0.025, 0.020000000000000004, 0.0, 0.005000000000000001, 0.0, 0.0, 0.0, 0.025, 0.0, 0.015, 0.0, 0.005, 0.03, 0.0, 0.0, 0.09, 0.005000000000000001, 0.005, 0.024999999999999998, 0.0, 0.005000000000000001, 0.024999999999999998, 0.004999999999999999, 0.024999999999999998, 0.024999999999999998, 0.0, 0.04, 0.01, 0.0, 0.0, 0.02, 0.005000000000000001, 0.005, 0.015, 0.0, 0.005000000000000001, 0.005, 0.005000000000000001, 0.020000000000000004, 0.02, 0.0, 0.015, 0.0, 0.0, 0.01855287569573284, 0.004999999999999999, 0.025, 0.02340823970037453, 0.005000000000000001, 0.005, 0.02, 0.004999999999999999, 0.015, 0.005, 0.004999999999999999, 0.0, 0.005, 0.02, 0.005, 0.02, 0.024999999999999998, 0.0, 0.010000000000000002, 0.004849660523763336, 0.019999999999999997, 0.025, 0.020000000000000004, 0.004999999999999999, 0.024999999999999998, 0.0, 0.03, 0.009999999999999998, 0.025, 0.03, 0.010000000000000002, 0.005, 0.01, 0.015, 0.02, 0.020000000000000004, 0.004999999999999999, 0.020000000000000004, 0.0, 0.02, 0.015, 0.0, 0.005000000000000001, 0.0, 0.024999999999999998, 0.015, 0.009852216748768471, 0.020000000000000004, 0.024999999999999998, 0.0, 0.0, 0.0, 0.004780114722753346, 0.01865671641791045, 0.03, 0.005, 0.024485798237022523, 0.005, 0.0, 0.02358490566037736, 0.01, 0.0, 0.015, 0.0, 0.0, 0.0, 0.005, 0.01926782273603083, 0.023629489603024575, 0.028653295128939826, 0.020000000000000004, 0.0, 0.0, 0.005000000000000001, 0.005000000000000001, 0.004999999999999999, 0.08964143426294821, 0.0, 0.015, 0.005000000000000001, 0.0, 0.004770992366412214, 0.01, 0.019685039370078743, 0.005, 0.014662756598240468, 0.03, 0.013927576601671309, 0.022789425706472195, 0.0, 0.005, 0.03, 0.009099181073703366, 0.009999999999999998, 0.00946073793755913, 0.025, 0.015, 0.02, 0.034999999999999996, 0.0, 0.019474196689386564, 0.005000000000000001, 0.004599816007359705, 0.0046210720887245845, 0.005, 0.024999999999999998, 0.015, 0.005, 0.020000000000000004, 0.020000000000000004, 0.0, 0.03, 0.03, 0.0, 0.025, 0.04, 0.019999999999999997, 0.03, 0.025, 0.04, 0.015, 0.009999999999999998, 0.0, 0.005, 0.0, 0.004999999999999999, 0.009999999999999998, 0.005, 0.019999999999999997, 0.02, 0.005000000000000001, 0.005, 0.01, 0.0, 0.005, 0.015, 0.02, 0.024999999999999998, 0.03, 0.004999999999999999, 0.005, 0.04, 0.004999999999999999, 0.005, 0.0, 0.03, 0.020000000000000004, 0.024999999999999998, 0.039999999999999994, 0.005, 0.005000000000000001, 0.015, 0.0, 0.0, 0.03, 0.0, 0.025, 0.0, 0.015, 0.004999999999999999, 0.005000000000000001, 0.005, 0.004999999999999999, 0.005, 0.0, 0.019999999999999997, 0.025, 0.0, 0.004999999999999999, 0.01, 0.009999999999999998, 0.020000000000000004, 0.0, 0.005, 0.075, 0.004999999999999999, 0.005000000000000001, 0.02, 0.005, 0.005000000000000001, 0.04, 0.004999999999999999, 0.005000000000000001, 0.009999999999999998, 0.015, 0.025, 0.0, 0.015, 0.025, 0.019999999999999997, 0.004999999999999999, 0.015, 0.0, 0.0, 0.025, 0.004999999999999999, 0.020000000000000004, 0.015, 0.1, 0.010000000000000002, 0.005, 0.020000000000000004, 0.03, 0.0, 0.005, 0.005000000000000001, 0.005, 0.015, 0.0, 0.020000000000000004, 0.015, 0.0, 0.004999999999999999, 0.009999999999999998, 0.024999999999999998, 0.005000000000000001, 0.0, 0.005000000000000001, 0.034999999999999996, 0.005000000000000001, 0.01, 0.034999999999999996, 0.004999999999999999, 0.020000000000000004, 0.019999999999999997, 0.0, 0.01, 0.03, 0.01928640308582449, 0.020000000000000004, 0.019999999999999997, 0.03, 0.03, 0.0, 0.035, 0.004999999999999999, 0.0, 0.015, 0.024999999999999998, 0.009999999999999998, 0.0, 0.015, 0.0, 0.0, 0.02, 0.024999999999999998, 0.03305004721435316, 0.020000000000000004, 0.03, 0.019999999999999997, 0.0, 0.004999999999999999, 0.009999999999999998, 0.010000000000000002, 0.0, 0.020000000000000004, 0.0, 0.004999999999999999, 0.020000000000000004, 0.03, 0.01, 0.02, 0.020000000000000004, 0.004999999999999999, 0.024999999999999998, 0.0, 0.005000000000000001, 0.0, 0.0, 0.015, 0.02, 0.005, 0.015, 0.005, 0.025, 0.01, 0.03, 0.034999999999999996, 0.024999999999999998, 0.0, 0.009999999999999998, 0.023809523809523808, 0.015, 0.005000000000000001, 0.0, 0.005, 0.095, 0.01, 0.039999999999999994, 0.01, 0.005000000000000001, 0.004975124378109453, 0.025, 0.004999999999999999, 0.03, 0.023809523809523808, 0.023854961832061067, 0.0, 0.0, 0.01489572989076465, 0.004999999999999999, 0.01, 0.025, 0.0, 0.0, 0.005000000000000001, 0.004668534080298786, 0.01, 0.025, 0.005, 0.005000000000000001, 0.0, 0.009999999999999998, 0.015, 0.0, 0.0, 0.0, 0.005, 0.03, 0.09597806215722121, 0.01439539347408829, 0.005000000000000001, 0.009999999999999998, 0.004734848484848485, 0.009999999999999998, 0.005000000000000001, 0.00940733772342427, 0.020000000000000004, 0.014866204162537165, 0.020000000000000004, 0.005, 0.045000000000000005, 0.01, 0.019999999999999997, 0.005000000000000001, 0.015000000000000001, 0.004999999999999999, 0.005, 0.0, 0.020000000000000004, 0.01, 0.024999999999999998, 0.020000000000000004, 0.01, 0.015, 0.039999999999999994, 0.015, 0.025, 0.03, 0.0, 0.015, 0.020000000000000004, 0.03, 0.01, 0.0, 0.0, 0.005, 0.005, 0.025, 0.024999999999999998, 0.03, 0.015, 0.009999999999999998, 0.005000000000000001, 0.004999999999999999, 0.03, 0.03, 0.0, 0.005, 0.004999999999999999, 0.0, 0.005000000000000001, 0.005, 0.01, 0.034999999999999996, 0.01, 0.025, 0.04, 0.034999999999999996, 0.005, 0.01, 0.0, 0.015, 0.075, 0.0, 0.039999999999999994, 0.020000000000000004, 0.005000000000000001, 0.0, 0.005, 0.024999999999999998, 0.015, 0.004999999999999999, 0.12, 0.015, 0.0, 0.0, 0.004999999999999999, 0.01, 0.015, 0.0, 0.03, 0.0, 0.010000000000000002, 0.0, 0.009999999999999998, 0.0, 0.015, 0.015, 0.005000000000000001, 0.005, 0.0, 0.024999999999999998, 0.004999999999999999, 0.0, 0.034999999999999996, 0.025, 0.025, 0.0, 0.004999999999999999, 0.01, 0.0, 0.0, 0.01, 0.034999999999999996, 0.01, 0.03, 0.06, 0.075, 0.004999999999999999, 0.025, 0.01, 0.02, 0.005, 0.034999999999999996, 0.015, 0.01, 0.005000000000000001, 0.01, 0.005, 0.025, 0.0, 0.020000000000000004, 0.03, 0.0, 0.020000000000000004, 0.005, 0.009999999999999998, 0.0, 0.01, 0.0, 0.095, 0.01, 0.004999999999999999, 0.024999999999999998, 0.01, 0.014955134596211365, 0.0, 0.005000000000000001, 0.03, 0.03, 0.02, 0.019999999999999997, 0.02, 0.0, 0.015, 0.0, 0.01, 0.025, 0.009999999999999998, 0.014018691588785045, 0.024999999999999998, 0.0, 0.004906771344455349, 0.0, 0.025, 0.005000000000000001, 0.02, 0.024999999999999998, 0.005, 0.005000000000000001, 0.020000000000000004, 0.034584980237154145, 0.020000000000000004, 0.029940119760479042, 0.009999999999999998, 0.024999999999999998, 0.04, 0.019999999999999997, 0.005, 0.0, 0.009999999999999998, 0.0, 0.015, 0.004995004995004995, 0.0, 0.01, 0.004999999999999999, 0.0, 0.005000000000000001, 0.005000000000000001, 0.010000000000000002, 0.00496031746031746, 0.005000000000000001, 0.0, 0.03, 0.020000000000000004, 0.01, 0.1, 0.0, 0.0, 0.009999999999999998, 0.025, 0.005000000000000001, 0.0, 0.025, 0.015, 0.015, 0.005000000000000001, 0.005000000000000001, 0.039999999999999994, 0.0, 0.014395393474088292, 0.0, 0.025, 0.03, 0.019512195121951223, 0.015, 0.005, 0.005, 0.024999999999999998, 0.005, 0.018535681186283598, 0.0, 0.0, 0.013914656771799629, 0.015, 0.024038461538461536, 0.0, 0.015, 0.020000000000000004, 0.034999999999999996, 0.0, 0.0, 0.0, 0.005000000000000001, 0.014778325123152709, 0.004625346901017576, 0.025, 0.03, 0.0, 0.005000000000000001, 0.025, 0.004999999999999999, 0.005, 0.0, 0.004999999999999999, 0.04, 0.004999999999999999, 0.0, 0.009999999999999998, 0.005000000000000001, 0.005, 0.03, 0.004999999999999999, 0.025, 0.009999999999999998, 0.004999999999999999, 0.01, 0.020000000000000004, 0.005000000000000001, 0.005, 0.01, 0.045, 0.025, 0.0, 0.005000000000000001, 0.0, 0.01, 0.020000000000000004, 0.020000000000000004, 0.004999999999999999, 0.0, 0.0, 0.019999999999999997, 0.004999999999999999, 0.01, 0.01, 0.01, 0.025, 0.035, 0.0, 0.024999999999999998, 0.015, 0.005000000000000001, 0.005, 0.01, 0.015, 0.095, 0.0, 0.0, 0.005, 0.025, 0.0, 0.0, 0.020000000000000004, 0.005, 0.024999999999999998, 0.0, 0.004999999999999999, 0.0, 0.01, 0.0, 0.020000000000000004, 0.0, 0.02, 0.004999999999999999, 0.024999999999999998, 0.0, 0.0, 0.0, 0.005000000000000001, 0.1, 0.025, 0.024999999999999998, 0.03, 0.0, 0.015, 0.03, 0.020000000000000004, 0.005000000000000001, 0.024999999999999998, 0.004999999999999999, 0.0, 0.005000000000000001, 0.024999999999999998, 0.025, 0.01, 0.01, 0.03, 0.0, 0.015, 0.009999999999999998, 0.034999999999999996, 0.0, 0.0, 0.0, 0.025, 0.025, 0.004999999999999999, 0.0, 0.020000000000000004, 0.005000000000000001, 0.005000000000000001, 0.0, 0.0, 0.005000000000000001, 0.005, 0.03, 0.015, 0.02, 0.005000000000000001, 0.01904761904761905, 0.009999999999999998, 0.005000000000000001, 0.009999999999999998, 0.0, 0.009999999999999998, 0.0, 0.015, 0.075, 0.03, 0.004999999999999999, 0.01, 0.0, 0.005000000000000001, 0.03, 0.034999999999999996, 0.0, 0.009871668311944718, 0.005000000000000001, 0.03, 0.0, 0.024999999999999998, 0.0, 0.025, 0.01, 0.009999999999999998, 0.009999999999999998, 0.005000000000000001, 0.020000000000000004, 0.005000000000000001, 0.045000000000000005, 0.005, 0.01, 0.005, 0.0, 0.0, 0.025, 0.005, 0.005000000000000001, 0.020000000000000004, 0.015, 0.009999999999999998, 0.02, 0.020000000000000004, 0.034999999999999996, 0.0, 0.005000000000000001, 0.020000000000000004, 0.015, 0.020000000000000004, 0.0, 0.005000000000000001, 0.025, 0.004999999999999999, 0.010000000000000002, 0.020000000000000004, 0.009999999999999998, 0.005, 0.015, 0.0, 0.009999999999999998, 0.004887585532746823, 0.015, 0.015, 0.005000000000000001, 0.020000000000000004, 0.015, 0.08, 0.020000000000000004, 0.024999999999999998, 0.025, 0.024999999999999998, 0.024999999999999998, 0.0, 0.005000000000000001, 0.03, 0.01, 0.020000000000000004, 0.03, 0.02, 0.005000000000000001, 0.03, 0.005, 0.0, 0.0, 0.01, 0.005000000000000001, 0.005, 0.004849660523763337, 0.009999999999999998, 0.0, 0.0, 0.03, 0.009643201542912245, 0.04, 0.0, 0.005, 0.03, 0.004999999999999999, 0.0, 0.0049019607843137246, 0.0, 0.005, 0.005000000000000001, 0.0, 0.0, 0.00967117988394584, 0.020000000000000004, 0.0, 0.015, 0.00463821892393321, 0.005, 0.0, 0.025, 0.03, 0.005000000000000001, 0.0, 0.020000000000000004, 0.005000000000000001, 0.0, 0.03, 0.0, 0.004999999999999999, 0.005, 0.0, 0.015, 0.0, 0.004999999999999999, 0.025, 0.0, 0.004999999999999999, 0.005, 0.03, 0.0, 0.0, 0.0, 0.0, 0.005000000000000001, 0.019999999999999997, 0.0, 0.01, 0.03, 0.01, 0.0, 0.005000000000000001, 0.004999999999999999, 0.005, 0.025, 0.005000000000000001, 0.009999999999999998, 0.020000000000000004, 0.0, 0.005000000000000001, 0.0, 0.0, 0.0, 0.01, 0.025, 0.0, 0.0, 0.004999999999999999, 0.01, 0.005, 0.009999999999999998, 0.025, 0.005000000000000001, 0.0, 0.020000000000000004, 0.02, 0.0, 0.03, 0.0, 0.020000000000000004, 0.004999999999999999, 0.0, 0.009999999999999998, 0.03, 0.009999999999999998, 0.025, 0.005000000000000001, 0.020000000000000004, 0.01, 0.0, 0.0, 0.005000000000000001, 0.020000000000000004, 0.0, 0.03, 0.03, 0.0, 0.0, 0.009999999999999998, 0.0, 0.0, 0.004999999999999999, 0.03, 0.005000000000000001, 0.005000000000000001, 0.04, 0.0, 0.0, 0.01, 0.005, 0.004999999999999999, 0.005000000000000001, 0.0, 0.015, 0.0, 0.0, 0.024999999999999998, 0.0, 0.085, 0.019999999999999997, 0.0, 0.03, 0.015, 0.0, 0.03, 0.005000000000000001, 0.08, 0.0, 0.0, 0.005, 0.02, 0.015, 0.005000000000000001, 0.020000000000000004, 0.004999999999999999, 0.0, 0.005000000000000001, 0.005, 0.025, 0.025, 0.0, 0.034999999999999996, 0.005000000000000001, 0.0, 0.009999999999999998, 0.03, 0.009999999999999998, 0.0, 0.0, 0.005000000000000001, 0.019999999999999997, 0.020000000000000004, 0.0, 0.015, 0.005, 0.015, 0.005000000000000001, 0.005, 0.005000000000000001, 0.0, 0.019999999999999997, 0.015, 0.0, 0.025, 0.005, 0.005, 0.024999999999999998, 0.0, 0.024999999999999998, 0.02, 0.0, 0.02, 0.009999999999999998, 0.02, 0.095, 0.005000000000000001, 0.015, 0.0, 0.005000000000000001, 0.0, 0.015, 0.025, 0.005, 0.0, 0.0, 0.005000000000000001, 0.03, 0.015, 0.019665683382497544, 0.005000000000000001, 0.0, 0.01, 0.024999999999999998, 0.0, 0.085, 0.025, 0.005000000000000001, 0.005, 0.01, 0.015, 0.020000000000000004, 0.01, 0.005000000000000001, 0.015000000000000003, 0.01, 0.01, 0.0, 0.0, 0.02321262766945218, 0.005, 0.0, 0.03, 0.0, 0.01, 0.0, 0.0, 0.018744142455482664, 0.0, 0.01, 0.0, 0.005, 0.0, 0.020000000000000004, 0.010000000000000002, 0.0, 0.02, 0.02, 0.004985044865403789, 0.024999999999999998, 0.03, 0.02, 0.00463821892393321, 0.03, 0.005, 0.02, 0.0, 0.0, 0.004999999999999999, 0.004816955684007708, 0.005000000000000001, 0.024999999999999998, 0.015, 0.0, 0.0, 0.019455252918287938, 0.0348605577689243, 0.009737098344693282, 0.015, 0.005000000000000001, 0.034999999999999996, 0.025, 0.015, 0.009999999999999998, 0.005, 0.015, 0.025, 0.035, 0.004999999999999999, 0.009999999999999998, 0.005, 0.0, 0.005000000000000001, 0.02, 0.025, 0.0, 0.0, 0.020000000000000004, 0.0, 0.0, 0.04, 0.005, 0.024999999999999998, 0.004999999999999999, 0.005000000000000001, 0.0, 0.005000000000000001, 0.005, 0.0, 0.020000000000000004, 0.005, 0.0, 0.0, 0.045, 0.005000000000000001, 0.005, 0.005000000000000001, 0.015, 0.03, 0.0, 0.024999999999999998, 0.015, 0.005000000000000001, 0.020000000000000004, 0.095, 0.020000000000000004, 0.0, 0.0, 0.020000000000000004, 0.004999999999999999, 0.0, 0.005000000000000001, 0.0, 0.005, 0.0, 0.004999999999999999, 0.02, 0.005, 0.005, 0.020000000000000004, 0.004999999999999999, 0.025, 0.03, 0.0, 0.0, 0.024999999999999998, 0.020000000000000004, 0.024999999999999998, 0.004999999999999999, 0.025, 0.009999999999999998, 0.02, 0.03, 0.019999999999999997, 0.004999999999999999, 0.005, 0.005000000000000001, 0.0, 0.009999999999999998, 0.02, 0.02, 0.03, 0.005, 0.004999999999999999, 0.0, 0.020000000000000004, 0.005000000000000001, 0.055, 0.03, 0.0, 0.01, 0.004999999999999999, 0.005, 0.015, 0.005000000000000001, 0.004999999999999999, 0.005000000000000001, 0.0, 0.020000000000000004, 0.005000000000000001, 0.0, 0.01, 0.024999999999999998, 0.024999999999999998, 0.0, 0.015, 0.015, 0.0, 0.004999999999999999, 0.004633920296570899, 0.03, 0.0, 0.0, 0.0, 0.03, 0.025, 0.005, 0.005, 0.005000000000000001, 0.005, 0.02, 0.024999999999999998, 0.02, 0.01, 0.015, 0.02, 0.03, 0.0, 0.02, 0.015000000000000003, 0.1, 0.005, 0.015, 0.019999999999999997, 0.005000000000000001, 0.005000000000000001, 0.01, 0.0, 0.005, 0.0, 0.005, 0.009999999999999998, 0.005, 0.0, 0.015, 0.01, 0.009999999999999998, 0.005000000000000001, 0.095, 0.03, 0.009999999999999998, 0.0, 0.009999999999999998, 0.0, 0.009970089730807577, 0.02, 0.005000000000000001, 0.015, 0.005, 0.015, 0.019665683382497544, 0.03, 0.009398496240601503, 0.01, 0.02, 0.005, 0.005, 0.004999999999999999, 0.023786869647954328, 0.013850415512465372, 0.005, 0.020000000000000004, 0.0, 0.0, 0.009523809523809525, 0.005000000000000001, 0.01, 0.034999999999999996, 0.005, 0.0, 0.015000000000000003, 0.024999999999999998, 0.0, 0.0, 0.0, 0.015, 0.0, 0.020000000000000004, 0.01, 0.03, 0.005000000000000001, 0.015, 0.0, 0.020000000000000004, 0.0, 0.025, 0.015, 0.075, 0.02, 0.0048828125, 0.004734848484848484, 0.0, 0.02, 0.034999999999999996, 0.0, 0.025, 0.034999999999999996, 0.0, 0.025, 0.02, 0.0, 0.02, 0.0, 0.020000000000000004, 0.005, 0.01, 0.0, 0.03, 0.0, 0.025, 0.005000000000000001, 0.020000000000000004, 0.03]))]



def accuracy_pics():

    data = get_accuracy_data()

    max_runtime = 0
    #get the maximal runtime
    for _, _, _, _, runtime_mat, _ in data:
        for row in runtime_mat:
            for ele in row:
                if ele > max_runtime:
                    max_runtime = ele

    #make all runtimes percentages and invert them so that the fastest is RED
    for i0, x in enumerate(data):
        for i1, row in enumerate(x[4]):
            for i2, ele in enumerate(row):
                if not ele is nan:
                    data[i0][4][i1][i2] = 1 - ele / max_runtime

    plots = [[], []]

    
    color_mapper = LinearColorMapper(
                    palette=heatmap_palette(light_spec_approximation, 256),
                    low=0,
                    high=1
                )
    color_mapper_fake = LinearColorMapper(
                    palette=heatmap_palette(light_spec_approximation, 256),
                    low=0,
                    high=max_runtime
                )
    plot = figure(
        plot_width=resolution, plot_height=resolution
        )
    color_bar = ColorBar(color_mapper=color_mapper, border_line_color=None, location=(0,0))
    color_bar.major_label_text_font=font
    color_bar.major_label_text_font_size="12pt"
    plot.add_layout(color_bar, 'left')
    plots[0].append(plot)
    plot1 = figure(
        plot_width=resolution, plot_height=resolution
        )
    color_bar1 = ColorBar(color_mapper=color_mapper_fake, border_line_color=None, location=(0,0))
    color_bar1.major_label_text_font=font
    color_bar1.major_label_text_font_size="12pt"
    plot1.add_layout(color_bar1, 'left')
    plots[1].append(plot1)

    for w, h, approach, qual_mat, runtime_mat, mapping_qual in data:
        def draw_matrix(plot, matrix, w, h, palette):
            m_h = w/float(len(matrix))
            m_w = h/float(len(matrix[0]))
            xs = []
            ys = []
            xs2 = []
            ys2 = []
            cs = []
            for y, row in enumerate(matrix):
                for x, ele in enumerate(row):
                    xs.append(x*m_w)
                    ys.append(y*m_h)
                    xs2.append((x+1)*m_w)
                    ys2.append((y+1)*m_h)
                    if ele is nan:
                        cs.append("#AAAAAA")
                    else:
                        cs.append(format(palette(ele)))
            plot.quad(left=xs, bottom=ys, right=xs2, top=ys2, color=cs)

        print(approach)
        plot = figure(
            x_range=(0,h), y_range=(0,w),
            plot_width=resolution, plot_height=resolution
            )
        plot.axis.visible = False

        draw_matrix(plot, qual_mat, w, h, light_spec_approximation)
        plots[0].append(plot)

        plot2 = figure(
            x_range=(0,h), y_range=(0,w),
            plot_width=resolution, plot_height=resolution
            )
        plot2.axis.visible = False

        draw_matrix(plot2, runtime_mat, w, h, light_spec_approximation)
        plots[1].append(plot2)

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
    return [[1.0, 0.07402376910016978, 0.061113527620704655, 0.04037948878783127, 0.028959882941104742, 0.025674827369742623, 0.01832806664751508, 0.010107752455420997, 0.009546188867675135, 0.006400620666246424, 0.003928045351069053, 0.0031586440074739745, 0.002542298588585956, 0.0012496229585900805, 0.0010348967173076126, 0.0006357009662654687, 0.00045295449866172534, 0.0003645643456069996, 0.00015926737009755126, 0.00015804030027657053], [0.05301204819277108, 0.08266381003392251, 0.05579969340827798, 0.04077892325315006, 0.02601795238714713, 0.019736460207813316, 0.013416176389728928, 0.010553693792827489, 0.010101010101010102, 0.004962002682163612, 0.0036750229688935557, 0.001931662087912088, 0.0017648831303000302, 0.0008607615690453744, 0.0004555997349237906, 0.0004743833017077799, 0.00024132244700961268, 0.00011516314779270633, 0.00015541826941757004, 0.00011264644037248422], [0.03238224637681159, 0.0630995811376841, 0.039045719545872964, 0.03429524361948956, 0.021113812941106695, 0.017413365623017932, 0.012169860176074573, 0.008339615529589145, 0.004773988215560621, 0.003953107960741549, 0.0018450958191806098, 0.0014699706005879883, 0.001004016064257028, 0.0005988263004511158, 0.0006027243138988226, 0.00023719165085388995, 7.781192856864957e-05, 0.0, 7.604562737642586e-05, 0.0], [0.03527131782945737, 0.04793028322440087, 0.03616523772408418, 0.022574998490976037, 0.020766773162939296, 0.012110645115687291, 0.009251014148609875, 0.006366139022711631, 0.004147071398033347, 0.0026106362373152555, 0.001585711901185111, 0.0009175842509175843, 0.0009302701828183142, 0.000281135788585887, 0.0001561706945691641, 0.0003270705382127412, 7.4827895839569e-05, 0.0, 3.796795504594123e-05, 3.763218304293832e-05], [0.030835734870317003, 0.035834358343583436, 0.029659882163899302, 0.019268009295120063, 0.01489247311827957, 0.01131401349741961, 0.006973246433617428, 0.004963878916810708, 0.0029681361850720208, 0.00209213774634922, 0.0008897516783952115, 0.001251160350324898, 0.0004701457451810061, 0.000304843196280913, 7.44324525493115e-05, 7.39125614398167e-05, 7.699711260827719e-05, 7.677248474146866e-05, 0.0, 0.0], [0.028275318238143516, 0.03268895710611466, 0.020122528266175225, 0.019185478775346085, 0.012622088655146507, 0.007957796852646639, 0.005089278012593806, 0.0028968470548721608, 0.0022951760318045823, 0.0009748869131180783, 0.00047818290496114764, 0.0003892868265337901, 0.0002753953890943426, 0.0002204261572373255, 0.00015312763188117295, 7.303268212525105e-05, 0.0003715814506539834, 0.0, 0.0, 0.0], [0.02882882882882883, 0.02873644175235949, 0.019257373692656742, 0.014036975447520297, 0.009246287475483329, 0.005674873406670159, 0.003791795932437091, 0.002881844380403458, 0.0013808340237503453, 0.0008981217540708345, 0.0005094443138176973, 0.00019507627482345596, 7.435220640172497e-05, 0.00011008366358432408, 7.252157516861267e-05, 3.737060428267125e-05, 3.589504289457626e-05, 0.0, 0.0, 0.0], [0.018033614657721995, 0.02521339461588969, 0.017266034637995867, 0.01061434544869733, 0.007503069437497158, 0.005324444829228172, 0.002579430068784802, 0.0019298034012784947, 0.0012729026036644166, 0.0006324875362750204, 0.00022814555686528004, 0.00014631648255175947, 0.0, 0.0, 7.220216606498195e-05, 0.0, 0.0, 0.0, 0.0, 0.0], [0.022877396487836314, 0.017508100763039616, 0.01140366326171187, 0.008701169219613886, 0.004718664498988857, 0.0028190568241311265, 0.0019323290480321791, 0.0007686985932815743, 0.00034324942791762013, 0.0002669615956675947, 0.00011167777240069985, 7.11111111111111e-05, 0.0, 0.0, 0.0, 3.5249744439352814e-05, 0.0, 0.0, 0.0, 0.0], [0.017027577271885988, 0.016422287390029325, 0.009567421702510426, 0.00611865980251498, 0.004253386190948134, 0.0017884914463452567, 0.0009669683607952347, 0.0005364190198858193, 0.0001789100797938956, 7.507507507507507e-05, 0.00011318619128466328, 3.620040544454098e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.01823840647370059, 0.014235728322868235, 0.00928753712923064, 0.0046764459011151525, 0.002071329040033141, 0.0016102442203734233, 0.00038082181347347575, 0.0002964500111168754, 0.00025643843645821885, 3.6360991927859794e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.015089163237311385, 0.010778778405114087, 0.005698244437926845, 0.0031037123624491113, 0.0013949163050216986, 0.0008353267266583134, 0.00029847405141215536, 0.0002895927601809955, 0.00017471521420085262, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.5061884225658285e-05, 0.0], [0.012, 0.007449190492593868, 0.003663869518969389, 0.002426559799992647, 0.0010510292838503915, 0.00021689621516104545, 0.00015015015015015014, 7.091192738618636e-05, 3.590148632153371e-05, 3.484563384207959e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.00959507850518777, 0.006396053676685757, 0.0023310914093549375, 0.001149467907597612, 0.0004889607703012751, 0.00010739600486861888, 0.00010636034886194426, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.007257006098973211, 0.0033747450398664935, 0.0015317286652078775, 0.0004296301600372346, 3.72300819061802e-05, 6.90894016857814e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.005176181894865994, 0.002389129460952665, 0.000460878505335555, 0.00010399694942281693, 3.447324875896305e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.004314423040532869, 0.0013797109868564374, 0.00020668274199104375, 0.0001343499143519296, 3.474876641879213e-05, 0.0, 3.4693311129614214e-05, 0.0, 0.0, 0.0], [0.002300007187522461, 0.0006520247083047358, 0.00013566220111921315, 3.5052052297662026e-05, 3.36802398033074e-05, 0.0, 0.0, 0.0], [0.0008107834194790716, 0.00033796343235661903, 0.00013373900832525326, 0.0, 0.0], [0.0013540290620871862, 9.883050568275408e-05, 6.6711140760507e-05]]

def max_spannning():
    return [[1.0, 0.15651006711409396, 0.2147962830593281, 0.10887096774193548, 0.07350290406569197, 0.07815321124580861, 0.04777441155907714, 0.02458120903131828, 0.022973271482217804, 0.014291632145816073, 0.00854052002277472, 0.0060287825748735905, 0.005654281098546042, 0.0027569909413154787, 0.0013235016071090943, 0.0011531808571977706, 0.0007518796992481203, 0.0003601656762110571, 0.0, 0.00035549235691432633], [0.3488372093023256, 0.3284313725490196, 0.1845748187211602, 0.13802476593174268, 0.07839632277834525, 0.05513051305130513, 0.032708142726440986, 0.026631739794804626, 0.021042084168336674, 0.012104095218882388, 0.006867299535447385, 0.0029767441860465114, 0.003945885005636978, 0.0011689070718877848, 0.0009416195856873823, 0.0010883366588064574, 0.0003677822728944465, 0.00017937219730941703, 0.00018443378827001107, 0.0], [0.10634648370497427, 0.281936127744511, 0.09745923397800531, 0.10787011557512383, 0.05796178343949045, 0.04922644163150492, 0.02983081032947462, 0.01596688350088705, 0.010555830446973445, 0.007722813608849927, 0.0032567396417586395, 0.002226758211170904, 0.0018740629685157421, 0.0010212765957446808, 0.0011084426380934787, 0.0006943239021003298, 0.00018198362147406734, 0.0, 0.0, 0.0], [0.16910935738444194, 0.2046142208774584, 0.12311135982092893, 0.05888144938952344, 0.06381909547738693, 0.029878273699741793, 0.023508137432188065, 0.013416971597839345, 0.00865265760197775, 0.00546448087431694, 0.0021798365122615805, 0.0011755485893416929, 0.002616333395626986, 0.00035797386790764276, 0.00018011527377521613, 0.00029073993312981537, 0.0, 0.0, 0.0, 0.0], [0.16609672691744015, 0.0931809376210771, 0.10190615835777127, 0.044613259668508286, 0.0441976427923844, 0.03154643723803221, 0.016170380595543286, 0.012082310741929394, 0.007731958762886598, 0.0043079228319910096, 0.001281112737920937, 0.002971508477538892, 0.0005362888809438684, 0.000544464609800363, 0.0, 0.0, 0.00017809439002671417, 0.0, 0.0, 0.0], [0.13922651933701657, 0.11764705882352941, 0.05609679446888749, 0.05824800910125142, 0.03738923708666522, 0.01832508704416346, 0.010671936758893281, 0.004820825968182549, 0.003784465669489998, 0.0013971358714634998, 0.0011169024571854058, 0.00017500875043752187, 0.0, 0.0003359650596337981, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.1710863986313088, 0.10288486086290528, 0.05854579792256846, 0.04124561765312436, 0.0258952053408861, 0.014627396718719115, 0.00783435926133184, 0.006051437216338881, 0.0030775149067128293, 0.0012711094970038134, 0.001028630207440425, 0.00017543859649122806, 0.0, 0.0001639344262295082, 0.0, 0.0, 0.00016366612111292964, 0.0, 0.0, 0.0], [0.0650280627056319, 0.08002507312996239, 0.055896805896805894, 0.02930914166085136, 0.01875901875901876, 0.01284796573875803, 0.00417298937784522, 0.003922267783918702, 0.001890359168241966, 0.001722356183258698, 0.00016559032952475575, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.12654424040066778, 0.054383651944627555, 0.02849426063470628, 0.021734949028659356, 0.010596235963941167, 0.005853301627949515, 0.004690117252931323, 0.0013468013468013469, 0.0008880994671403197, 0.0003459010722933241, 0.00033400133600534405, 0.00033573946617424877, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.06887280248190279, 0.06423059464366773, 0.024919848440687845, 0.017079889807162536, 0.011227682679919013, 0.003192204301075269, 0.0019798416126709864, 0.0010070493454179255, 0.00016414970453053183, 0.0, 0.0, 0.00016257519102584944, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.08761408083441982, 0.05218422613164657, 0.024485951969863445, 0.011973875181422351, 0.0035864649929830033, 0.002099370188943317, 0.00032035880185808104, 0.0006830601092896175, 0.0003526093088857546, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.07487562189054726, 0.03662358010184097, 0.01687605475342209, 0.007557117750439367, 0.003493013972055888, 0.0014349775784753362, 0.0005036937541974479, 0.0008494733265375467, 0.00031660598385309483, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.057583659278574534, 0.023912634155526268, 0.008285004142502071, 0.0047017337643255955, 0.001043115438108484, 0.000671591672263264, 0.0003427004797806717, 0.0, 0.0001682935038707506, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.04020496649586126, 0.02172722293226219, 0.006035523366097603, 0.003212715590125127, 0.0006821282401091405, 0.0003318400530944085, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, nan, nan], [0.02130445132111402, 0.01032448377581121, 0.0027057138309724655, 0.0004935022207599934, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, nan, nan, nan, nan, nan], [0.018244191486043972, 0.006955312119631368, 0.0006335128286347799, 0.00032346757237586933, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, nan, nan, nan, nan, nan, nan, nan], [0.016657710908113917, 0.003234152652005175, 0.0008163265306122449, 0.00031279324366593683, 0.00015664160401002505, 0.0, 0.0, 0.0, 0.0, 0.0, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan], [0.007569386038687973, 0.001562011871290222, 0.00029282576866764275, 0.0001539408866995074, 0.00014947683109118088, 0.0, 0.0, 0.0, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan], [0.0017994438082774415, 0.0009453285016543249, 0.00046425255338904364, 0.0, 0.0, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan], [0.0038627247065814887, 0.00014947683109118088, 0.0, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan]]

def kmer():
    return [[0.03880134360864002, 0.03501559695380165, 0.027522354229295325, 0.017103798550756687, 0.01687007657030651, 0.0168277709883955, 0.012313454008965557, 0.009873850083956429, 0.008478648103876623, 0.005417321357379048, 0.0042886586491419215, 0.003728041426324069, 0.004236382102956113, 0.0022698685789466076, 0.0033150372262377046, 0.002206658351722153, 0.0006985679357317499, 0.0007738500588126045, 0.000385549600956163, 0.0005779575868047837], [0.030024953714883684, 0.036845353879344754, 0.02734769600326002, 0.023177453061538334, 0.012973716019244899, 0.010514073160749532, 0.008207925125408208, 0.009807747065585353, 0.009597031818212259, 0.005166028198566228, 0.006252649427723612, 0.003130710196126203, 0.003001742536981214, 0.0013456199219881108, 0.0009738236210657526, 0.0012834723761740855, 0.0007545034424219561, 0.0005469462169553327, 0.0004976609933313427, 0.0003824416776441593], [0.03994314622811423, 0.03369219025558685, 0.02432155314071308, 0.016662226689560532, 0.011515103714804006, 0.01188790934858943, 0.007430403404755738, 0.0051736260926538505, 0.006515403512519645, 0.004383936828626307, 0.00292936849980506, 0.0025077381634758685, 0.0030658003111558524, 0.0010164446219189023, 0.002036726782957196, 0.0006042462028619297, 0.000160926939169617, 0.00013969732246798603, 0.00025786487880350697, 0.0], [0.023539591088761726, 0.03156602889902687, 0.021965048247889897, 0.013458186291506167, 0.018190956878678086, 0.012063342913995136, 0.00944886566100487, 0.007694981882362379, 0.005608487622217599, 0.0052287771034602625, 0.003899740589943136, 0.0024129214405305146, 0.003119574652777778, 0.0007210768080333298, 0.0002932004158114988, 0.0006564993435006565, 5.356760231412042e-05, 0.0, 0.0, 0.0], [0.05169884597727836, 0.022900849432458844, 0.019926402361359046, 0.03587552918004272, 0.012684965505795554, 0.015303682448589193, 0.005782772705312831, 0.007412947180696751, 0.005456908814440575, 0.0036407607481303314, 0.002246214042460691, 0.0037464974970878064, 0.0007622463849344243, 0.0014788972735201547, 0.0003522012578616352, 0.00023969319271332693, 0.0002549232043846791, 0.0001184342985728667, 0.0, 0.0], [0.0363159924767483, 0.018336463431218464, 0.01368390963934398, 0.017677931335407654, 0.015901569896302117, 0.010519432119426493, 0.007394736842105263, 0.0046084574599405135, 0.004796982411064493, 0.002707526924851086, 0.001940396866076005, 0.0007667004674399624, 0.0006734460233012324, 0.0009191598878624937, 0.00033899581905156504, 5.438329345225147e-05, 0.0, 0.0, 0.0, 0.0], [0.03170467006568163, 0.020136324106986014, 0.016862172242812928, 0.014082975131823956, 0.009631414192617833, 0.007443363595657281, 0.00809456163392211, 0.005908547467641606, 0.0025372237430727113, 0.0018304960644334614, 0.001520300482918977, 0.00036007725293790303, 0.00020703933747412008, 0.0004046476674952309, 0.0003014954172696575, 0.0001988862370723946, 0.0, 0.0, 0.0, 0.0], [0.02722471169541234, 0.02403699055254794, 0.020235986171072775, 0.012795055166544738, 0.012189838686385606, 0.00792501144827715, 0.00512710673217226, 0.0035722710881946553, 0.0028473973009046417, 0.0012033379548486672, 0.0005663797009515179, 0.00040986966144765965, 0.0, 0.0, 5.676979846721544e-05, 0.0, 0.0, 0.0, 0.0, 0.0], [0.04305219834062847, 0.021434968653597373, 0.012995869378848102, 0.014231159377051137, 0.006344288793103448, 0.0036320561618660234, 0.006492041174122099, 0.0020760687722975834, 0.0011180818684390335, 0.0018347936793229986, 0.0005102040816326531, 0.0003538153084090105, 0.0, 0.0, 0.0, 0.0001501726986033939, 0.0, 0.0, 0.0, 0.0], [0.030558660122778886, 0.02277772972349566, 0.013358631013317681, 0.009441898948041745, 0.011940254449921776, 0.005011249744323993, 0.003351711072790421, 0.0013869014859658778, 0.0007842844338213325, 0.0004465780953444234, 0.0002176041780002176, 6.210408644888833e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.028930434007908498, 0.02048808483881722, 0.016699985372275586, 0.008283959732752255, 0.003472340040491783, 0.0041551893178953265, 0.0009830298013245033, 0.0010488855590934632, 0.000440392830404721, 4.930723337113555e-05, 0.0, 0.0, 0.00024238506907974468, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.034781416394341735, 0.028573609485522365, 0.009199290194043523, 0.0069684499314128946, 0.0042893856194649365, 0.002562316533160854, 0.0008425951931950407, 0.0013784954706577393, 0.0012140255659501536, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0375548454183469, 0.012934841793211931, 0.008673006719836306, 0.007150153217568948, 0.003885095361431599, 0.0021428238691945777, 0.0006481931615621455, 0.0003562776115148924, 0.000245378701128742, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.03074687935481958, 0.01828909164178179, 0.005210308093100622, 0.00474774053311978, 0.001553291607656814, 0.0009315530319833207, 0.0008330556481172943, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.022148370186947014, 0.010415035238841033, 0.005568379702749767, 0.00330729996812241, 0.0006100994931481134, 0.00032615786040443573, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.02346136511592969, 0.007606399867474703, 0.001796986591715431, 0.000748652425633859, 0.0008168451624614268, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.023213084633430123, 0.006791171477079796, 0.00185592817989555, 0.0007124198527665638, 0.00025119316754584274, 0.0, 0.00021119324181626187, 0.0, 0.0, 0.0], [0.015131180188760932, 0.006465200568153989, 0.0012043356081894822, 0.0, 7.160246312473149e-05, 0.0, 0.0, 0.0], [0.003983741481466502, 0.0037062937062937065, 0.0014239048496523998, 0.0, 0.0], [0.019125991057002217, 0.0009356287425149701, 0.0003507233669443227]]

def seed_relevance_pics():
    h = 400
    w = 20


    data = []
    high = 99
    low = 1
    for row in max_spannning():
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

seed_relevance_pics()
#accuracy_pics()
#unrelated_non_enclosed_seeds()
#forced_gap()
#ambiguity_per_length()
#theoretical_max_acc()
#seed_shadows()
#alignment()
#stripOfConsideration()
#optimal_matching()
#required_nmw_band_size()

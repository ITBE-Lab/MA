from bokeh.palettes import Plasma, Plasma256

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

    return (i*r**m, i*g**m, i*b**m)

def format(rgb):
        def clamp(x):
            return max(0, min(x, 255))
        red, green, blue = rgb
        return "#{0:02x}{1:02x}{2:02x}".format(clamp(int(red * 255)), clamp(int(green * 255)),
                                               clamp(int(blue * 255)))

def html_file(name):
    with open("html/" + name + ".html", "r") as html_file:
        return html_file.read()

def css_file(name):
    with open("css/" + name + ".css", "r") as css_file:
        return css_file.read()

def decode(o):
    if isinstance(o, str):
        try:
            return float(o)
        except ValueError:
            return o
    elif isinstance(o, dict):
        return {decode(k): decode(v) for k, v in o.items()}
    elif isinstance(o, list):
        return [decode(v) for v in o]
    else:
        return o

def append_nuc_type(dict_, nuc, pos, pos_key):
    if nuc == "A" or nuc == "a":
        dict_["c"].append("blue")
        dict_["i"].append("A @" + str(pos))
    elif nuc == "C" or nuc == "c":
        dict_["c"].append("red")
        dict_["i"].append("C @" + str(pos))
    elif nuc == "G" or nuc == "g":
        dict_["c"].append("green")
        dict_["i"].append("G @" + str(pos))
    elif nuc == "T" or nuc == "t":
        dict_["c"].append("yellow")
        dict_["i"].append("T @" + str(pos))
    else:
        dict_["c"].append("lightgray")
        dict_["i"].append(nuc + " @" + str(pos))
    dict_[pos_key].append(pos + 0.5)

def add_seed(seed, read_dict, max_seed_size, end_column, all_col_ids, category_counter, parlindrome, layer,
             read_id, idx, r_name):
    seed_size = seed.size - 1
    if seed.on_forward_strand:
        read_dict["center"].append(seed.start_ref + seed.size/2)
        read_dict["r"].append(seed.start_ref)
        read_dict["x"].append(
            [seed.start_ref, seed.start_ref+seed.size])
        curr_end = seed.start_ref + seed_size + 1
        curr_start = seed.start_ref
    else:
        read_dict["center"].append(seed.start_ref - seed.size/2 + 1)
        read_dict["r"].append(seed.start_ref - seed.size + 1)
        read_dict["x"].append(
            [seed.start_ref + 1, seed.start_ref - seed.size + 1])
        curr_end = seed.start_ref
        curr_start = seed.start_ref - seed.size
    curr_column = 0
    while curr_column < len(end_column):
        add_dist = 100
        if end_column[curr_column][1] == read_id:
            add_dist = 0
        if curr_start > end_column[curr_column][0] + add_dist:
            break
        else:
            curr_column += 1
    if curr_column >= len(end_column):
        end_column.append( None )
        all_col_ids.append(curr_column + category_counter)
    end_column[curr_column] = (curr_end, read_id)
    read_dict["y"].append([seed.start, seed.start+seed.size])
    if layer == -1:
        read_dict["c"].append( Plasma256[ (255 * seed.size) // max_seed_size ] )
    else:
        read_dict["c"].append("lightgrey")
    read_dict["r_id"].append(read_id)
    read_dict["r_name"].append(r_name)
    read_dict["size"].append(seed.size)
    read_dict["l"].append(seed.size)
    read_dict["q"].append(seed.start)
    read_dict["idx"].append(idx)
    read_dict["layer"].append(layer)
    read_dict["parlindrome"].append(parlindrome)
    read_dict["f"].append(seed.on_forward_strand)
    read_dict["category"].append(category_counter + curr_column)
    read_dict["overlapping"].append(False)
    read_dict["in_soc_reseed"].append(False)
    read_dict["soc_nt"].append(0)
    read_dict["soc_id"].append(-1)
    read_dict["max_filter"].append(0)
    read_dict["min_filter"].append(0)
    read_dict["oc"].append("red")
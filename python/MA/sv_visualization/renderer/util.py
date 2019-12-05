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
def js_file(name):
    with open("js/" + name + ".js", "r") as js_file:
        return js_file.read()

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
        dict_["i"].append(nuc)
    dict_[pos_key].append(pos + 0.5)
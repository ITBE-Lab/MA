from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import FuncTickFormatter#, TapTool, OpenURL
from MA import *

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

def render_overview(xs, xe, ys, ye, db_path, ref_genome):
    pack = Pack()
    pack.load(ref_genome)
    sv_db = SV_DB(db_path, "open")

    run_id = 0
    min_score = 0

    plot = figure(
            x_range=(xs, xe),
            y_range=(ys, ye),
            width=900,
            height=900,
            tooltips=[("from:", "@f"), ("to:", "@t"), ("num calls:", "@i")],
            tools=[
                "pan", "wheel_zoom", "box_zoom", "save",
                "reset", "hover"#, "tap"
            ],
            active_drag="pan",
            active_scroll="wheel_zoom")
    plot.background_fill_color = "black"
    plot.grid.visible = False
    #plot.ygrid = None

    max_ = pack.unpacked_size_single_strand
    cds = {
            'x': [],
            'y': [],
            #'xe': [],
            #'ye': [],
            'w': [],
            'h': [],
            'c': [],
            'f': [],
            't': [],
            'i': []
        }
    starts = pack.contigStarts()
    lengths = pack.contigLengths()
    names = pack.contigNames()
    num_contigs = pack.num_contigs()
    num_calls_list = libMA.get_call_overview(sv_db, pack)
    max_ = max(num_calls_list)
    for idx, num_calls in enumerate(num_calls_list):
        if num_calls > 0:
            i = idx // num_contigs
            j = idx % num_contigs
            cds["x"].append(starts[i])
            cds["y"].append(starts[j])
            #cds["xe"].append(starts[i] + lengths[i])
            #cds["ye"].append(starts[j] + lengths[j])
            cds["w"].append(lengths[i])
            cds["h"].append(lengths[j])
            cds["c"].append(format(light_spec_approximation(num_calls/max_)))
            cds["f"].append(names[i])
            cds["t"].append(names[j])
            cds["i"].append(str(num_calls))
    plot.rect(x="x", y="y", width="w", height="h", color="c", line_width=0, source=ColumnDataSource(cds))

    plot.axis.formatter = FuncTickFormatter(
            args={"lengths":[x for x in lengths], "names":[x for x in names]}, 
            code="""
                    var i = 0;
                    while(tick > lengths[i])
                    {
                        tick -= lengths[i];
                        i += 1;
                        if(i > lengths.length)
                            return names[lengths.length-1];
                    }
                    return names[i];
                """)

    #url = "http://localhost:5006/bokeh_server?xs=@x&ys=@y&xe=@xe&ye=@ye"
    #taptool = plot.select(type=TapTool)
    #taptool.callback = OpenURL(url=url)

    return plot

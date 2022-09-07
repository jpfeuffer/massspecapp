import json

from bokeh.models import HoverTool, EdgesAndLinkedNodes, NodesAndLinkedEdges, CustomJS, ColumnDataSource, CustomJSHover

from bokeh.plotting import figure, show
import numpy as np
import pandas as pd

# prepare some data
x = np.array([1., 2., 3., 4., 5.])
y = np.array([6., 7., 2., 4., 5.])
annotList = ["1.+", "2.+", "3.+", "4.+", "5.+"]
annot = pd.DataFrame(annotList, columns=["data"])
annotSource = ColumnDataSource(annot)


# def createMultilineData(a,b,annot):
#    assert a.size == b.size and b.size == annot.size
#    c = np.empty(a.size, dtype=[('mz',a.dtype),('ity',b.dtype),('annot',) ...
#    c[0::3] = a
#    c[1::3] = 0
#    c[2::3] = np.nan
#    return c


def spaceWithZeroAndNan(a):
    c = np.empty((a.size * 3,), dtype=a.dtype)
    c[0::3] = a
    c[1::3] = 0
    c[2::3] = np.nan
    return c


def spaceWithPrevAndNan(a):
    c = np.empty((a.size * 3,), dtype=a.dtype)
    c[0::3] = a
    c[1::3] = a
    c[2::3] = np.nan
    return c


# create a new plot with a title and axis labels
p = figure(title="Simple line example", x_axis_label='m/z', y_axis_label='Intensity')

# add a line renderer with legend and line thickness to the plot
line = p.line(spaceWithPrevAndNan(x), spaceWithZeroAndNan(y), legend_label="Spectrum", line_width=2, name="spec")
p.scatter(x, y, name="pnts")

code = f"return ({json.dumps(annotList)})[special_vars.indices[0] / 3];"
formatter = CustomJSHover(code=code)
callback = CustomJS(args=dict(source=line.data_source, annot=annotSource),
                    code="""
                        console.log("Indices: "+ annot.data.data[cb_data['index'].line_indices / 3]);
                        """)

hover = HoverTool(
    tooltips=[
        # ("Series", "@columns"),
        ("m/z", "@x"),
        ("Intensity", "@y"),
        ("Annot", "$index{custom}")
    ],
    formatters={"$index": formatter}
)

hover.renderers = p.select("spec")
hover.line_policy = 'prev'
p.add_tools(hover)

# show the results
show(p)

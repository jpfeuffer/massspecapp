from functools import reduce
from pathlib import Path
from time import perf_counter

import bokeh.layouts
import datashader
import numpy as np
from bokeh.models import ColumnDataSource, Label, CustomJS, Div
from bokeh.plotting import figure
from bokeh.layouts import column, row, layout
import panel as pn
from bokeh.palettes import RdYlBu3
from holoviews.plotting.bokeh.callbacks import LinkCallback
from bokeh.plotting import figure, curdoc
from holoviews.plotting.links import Link
from holoviews.plotting.util import process_cmap
from numpy import log
from pyopenms import *
import pandas as pd
import os

import holoviews as hv
import holoviews.operation.datashader as hd
from holoviews import opts, dim, streams

pn.extension()
hv.extension('bokeh')
renderer = hv.renderer('bokeh').instance(mode='server')


def modify_doc(doc):
    """Add a plotted function to the document.

    Arguments:
        doc: A bokeh document to which elements can be added.
    """
    print("OMP enabled: ", OpenMSBuildInfo.isOpenMPEnabled())
    print("MaxThreads: ", OpenMSBuildInfo.getOpenMPMaxNumThreads())
    # Start the stopwatch / counter
    t1_start = perf_counter()
    exp = MSExperiment()
    loader = MzMLFile()
    opts = loader.getOptions()  # type: PeakFileOptions
    opts.setMSLevels([1])
    opts.setIntensity32Bit(True)
    loader.setOptions(opts)
    #loader.load("/Volumes/Data/UPS1/mzML/UPS1_250amol_R1.mzML", exp)
    loader.load(str(Path(__file__).resolve().parent) + "/static/data/BSA1.mzML", exp)
    load_stop = perf_counter()
    print("Time for loading mzML:",
          load_stop - t1_start)
    cnt = sum(spec.size() for spec in exp if spec.getMSLevel() == 1)

    cols = ["RT", "mzarray", "intarray"]
    expandcols = ["RT", "mz", "inty"]
    spectraarr = np.fromiter(((spec.getRT(), point[0], point[1]) for spec in exp if spec.getMSLevel() == 1
                              for point in zip(*spec.get_peaks())), dtype=[('RT', 'f'), ('mz', 'f'), ('inty', 'f')],
                             count=cnt)
    np_stop = perf_counter()
    print("Time for loading and creating numpy array:",
          np_stop - load_stop)
    # spectradfwide = pd.DataFrame(data=((spec.getRT(), *spec.get_peaks()) for spec in exp if spec.getMSLevel() == 1), columns=cols)
    # spectradf = pd.DataFrame(data=((spec.getRT(), point[0], point[1]) for spec in exp if spec.getMSLevel() == 1 for point in zip(*spec.get_peaks())), columns=expandcols)
    spectradf = pd.DataFrame(data=spectraarr, columns=expandcols)

    # Stop the stopwatch / counter
    df_stop = perf_counter()
    print("Time for loading and creating DF:",
          df_stop - np_stop)

    points = hv.Points(spectradf, kdims=['RT', 'mz'], vdims=['inty'], label="test").opts(
        fontsize={'title': 16, 'labels': 14, 'xticks': 6, 'yticks': 12},
        color=log(dim('int')),
        colorbar=True,
        cmap='Magma',
        width=1000,
        height=1000,
        tools=['hover'])

    # Unfortunately datashade cannot transfer interactiveness (e.g. hover) or colorbars
    # see https://anaconda.org/jbednar/datashade_vs_rasterize/notebook
    # shade = hd.datashade(points,cmap=process_cmap("blues", provider="bokeh"), alpha=250, min_alpha=10).opts(
    #    plot=dict(
    #        width=1000,
    #        height=1000))#, color_key=colors)
    raster = hd.rasterize(points, cmap=process_cmap("blues", provider="bokeh"), aggregator=datashader.sum('inty'),
                          cnorm='log', alpha=50, min_alpha=10).opts(
        tools=['hover']).opts(
        plot=dict(
            width=1000,
            height=1000,
            xlabel="loc",
            ylabel="time")
    )


    # selection = hv.streams.Selection1D(source=points)
    t2_stop = perf_counter()
    print("Time for creating plot objects:",
          t2_stop - df_stop)

    # 2. Instead of Jupyter's automatic rich display, render the object as a bokeh document
    # hv.renderer('bokeh').server_doc(hd.dynspread(shade), doc=doc)
    # Create HoloViews plot and attach the document
    finalplot = hd.dynspread(raster)

    jsupdateinfo = '''
    var d = [ {
  "time" : 1675.3155,
  "m/z array" : [ 359.2407531738281, 359.3271789550781, 360.2444763183594, 360.3111572265625, 360.3471374511719, 361.051025390625, 363.0480651855469 ],
  "intensity array" : [ 266890.28125, 6482.91064453125, 53019.37109375, 29689.6953125, 7935.1279296875, 42180.21484375, 27541.478515625 ]
}, {
  "time" : 1675.55064,
  "m/z array" : [ 359.2407531738281, 359.327392578125, 360.2445983886719, 360.3111267089844, 361.05078125, 363.04833984375 ],
  "intensity array" : [ 264650.5, 10177.9462890625, 57445.171875, 33437.2890625, 44480.46875, 25241.8046875 ]
}];
    console.log("xrange changes");
    renderAll(d);
        '''

    hvplot = renderer.get_plot(finalplot, doc)
    hvplot.state.x_range.js_on_change('start', CustomJS(code=jsupdateinfo))

    doc.title = 'HoloViews Bokeh App'

    layout = bokeh.layouts.layout([hvplot.state])

    doc.add_root(layout)
    return doc


def main():
    modify_doc(curdoc())


main()

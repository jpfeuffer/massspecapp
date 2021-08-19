from functools import reduce
from time import perf_counter

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

import holoviews as hv
import holoviews.operation.datashader as hd
from holoviews import opts, dim, streams

pn.extension()
hv.extension('bokeh')
renderer = hv.renderer('bokeh').instance(mode='server')


class TapLink(Link):
    _requires_target = True


class TapCallback(LinkCallback):
    source_model = ['xrange', 'yrange']
    source_handles = ['xrange', 'yrange']
    on_source_changes = ['xrange', 'yrange']

    target_model = 'text'

    source_code = """
                alert("Foo")
            """

    oldsource_code = """
        var inds = source_selected.indices
        var d = source_cds.data
        var vm = 0
        if (inds.length == 0)
            return
        for (var i = 0; i < inds.length; i++)
            vm += d[column][inds[i]]
        vm /= inds.length
        target_glyph.location = vm
    """

    def validate(self):
        assert self.link.column in self.source_plot.handles['cds'].data


TapLink.register_callback('bokeh', TapCallback)


def modify_doc(doc):
    """Add a plotted function to the document.

    Arguments:
        doc: A bokeh document to which elements can be added.
    """
    print("OMP enabled: ", OpenMSBuildInfo.isOpenMPEnabled())
    print("MaxThreads: ", OpenMSBuildInfo.getOpenMPMaxNumThreads())
    # Start the stopwatch / counter
    t1_start = perf_counter()
    exp = MSExperiment() # type: PeakMap
    loader = MzMLFile()
    opts = loader.getOptions()  # type: PeakFileOptions
    opts.setMSLevels([1])
    opts.setIntensity32Bit(True)
    loader.setOptions(opts)
    #loader.load("/Volumes/Data/UPS1/mzML/UPS1_250amol_R1.mzML", exp)
    loader.load("/Users/pfeuffer/git/OpenMS-fixes-src/share/OpenMS/examples/BSA/BSA1.mzML", exp)
    load_stop = perf_counter()
    print("Time for loading mzML:",
          load_stop - t1_start)
    cnt = sum(spec.size() for spec in exp if spec.getMSLevel() == 1)

    cols = ["RT", "mzarray", "intarray"]
    expandcols = ["RT", "mz", "inty"]
    #spectraarr = np.fromiter(((spec.getRT(), point[0], point[1]) for spec in exp if spec.getMSLevel() == 1
    #                          for point in zip(*spec.get_peaks())), dtype=[('RT', 'f'), ('mz', 'f'), ('inty', 'f')],
    #                         count=cnt)

    def createSpectrum(zippedMzIntarray):
        return hv.Points([zippedMzIntarray[0],zippedMzIntarray[1]], 'mz', 'inty')

    spectraptsdict = {spec.getRT(): createSpectrum(spec.get_peaks()) for spec in exp}

    np_stop = perf_counter()
    print("Time for loading and creating numpy array:",
          np_stop - load_stop)
    # spectradfwide = pd.DataFrame(data=((spec.getRT(), *spec.get_peaks()) for spec in exp if spec.getMSLevel() == 1), columns=cols)
    # spectradirectdf = pd.DataFrame(data=((spec.getRT(), point[0], point[1]) for spec in exp if spec.getMSLevel() == 1 for point in zip(*spec.get_peaks())), columns=expandcols)
    #spectradf = pd.DataFrame(data=spectraarr, columns=expandcols) # from np ndarray

    kdims = [hv.Dimension(('rt', 'RT'))]
    holomap = hv.HoloMap(spectraptsdict, kdims=kdims)

    # Stop the stopwatch / counter
    df_stop = perf_counter()
    print("Time for loading and creating DF:",
          df_stop - np_stop)

    points = holomap.overlay('rt', label="MS1").opts(
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
    #raster = hd.rasterize(points, cmap=process_cmap("blues", provider="bokeh"), aggregator=datashader.sum('inty'),
    #                      cnorm='log', alpha=50, min_alpha=10).opts(
    #    tools=['hover']).opts(
    #    plot=dict(
    #        width=1000,
    #        height=1000)
    #)

    # selection = hv.streams.Selection1D(source=points)
    t2_stop = perf_counter()
    print("Time for creating plot objects:",
          t2_stop - df_stop)

    # 2. Instead of Jupyter's automatic rich display, render the object as a bokeh document
    # hv.renderer('bokeh').server_doc(hd.dynspread(shade), doc=doc)
    # Create HoloViews plot and attach the document
    hvplot = renderer.get_plot(points, doc)


    doc.title = 'HoloViews Bokeh App'

    plot = layout([hvplot.state], sizing_mode='fixed')

    doc.add_root(plot)
    return doc


def main():
    modify_doc(curdoc())


main()

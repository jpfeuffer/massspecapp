import random
from functools import reduce, partial
from pathlib import Path
from time import perf_counter
import sys

import bokeh.layouts
import datashader
import numpy as np
from bokeh.models import ColumnDataSource, CustomJS, Button, Select, Div
import panel as pn
from bokeh.plotting import curdoc
from holoviews.plotting.util import process_cmap
from numpy import log
from pyopenms import MSExperiment, MzMLFile, PeakFileOptions
import pandas as pd
import os
import holoviews as hv
import holoviews.operation.datashader as hd
from holoviews import opts, dim
from tornado import gen
from tornado.ioloop import IOLoop

pn.extension()
hv.extension('bokeh')
renderer = hv.renderer('bokeh').instance(mode='server')


def modify_doc(doc):
    """Add a plotted function to the document.

    Arguments:
        doc: A bokeh document to which elements can be added.
    """

    exp = MSExperiment()
    loader = MzMLFile()
    opts = loader.getOptions()  # type: PeakFileOptions
    opts.setMSLevels([1])
    opts.setIntensity32Bit(True)
    loader.setOptions(opts)

    if len(sys.argv) > 1:
        file = sys.argv[1]
    else:
        file = "/Volumes/Data/UPS1/mzML/UPS1_250amol_R1.mzML"
        #file = str(Path(__file__).resolve().parent) + "/static/data/BSA1.mzML"

    jsupdateinfo = '''
        renderAllCDS(data, xr.start, xr.end, yr.start, yr.end);
        '''

    # Instead of listening to the range, we require a button press. Much easier and performant.
    #hvplot.state.x_range.js_on_change('start', CustomJS(code=jsupdateinfo))

    # If enabled via args, creates a button to stop the server
    def stopbutton_callback():
        sys.exit()  # Stop the server
    stopbutton = Button(label="Stop server", button_type="danger")
    stopbutton.on_click(stopbutton_callback)

    # Fake selectors that we use to trigger both a python function (for data filetering)
    # and a JS function to update the 3D View with the filtered data
    invisText = Select(title="Option:", value="foo", options=["foo", "bar", "baz", "quux"], visible=False, id="InvisibleText")
    invisText2 = Select(title="Option:", value="bar", options=["foo", "bar", "baz", "quux"], visible=False, id="InvisibleText2")


    js_on_invisText2 = '''
    console.log("js_on_invisText2 triggered")
    Bokeh.documents[0].get_model_by_id("InvisibleText").value = "bar" + new Date().timeNow();
    '''

    invisText2.js_on_change("value", CustomJS(code=js_on_invisText2))

    ## Layout
    doc.title = 'Mass-spectrometry Viewer'

    dynamic_col = bokeh.layouts.column(Div(text="Loading spectra..."))
    if len(sys.argv) > 2 and sys.argv[2] == "stoppable":
        layout = bokeh.layouts.layout(dynamic_col, [invisText,invisText2])
    else:
        layout = bokeh.layouts.layout(dynamic_col,[invisText,invisText2,stopbutton])

    def init_2D_plot_and_update_btn():
        # Start the stopwatch / counter
        t1_start = perf_counter()
        loader.load(file, exp)
        load_stop = perf_counter()
        print("Time for loading mzML:",
              load_stop - t1_start)
        cnt = sum(spec.size() for spec in exp)

        cols = ["RT", "mzarray", "intarray"]
        expandcols = ["RT", "mz", "inty"]
        spectraarr = np.fromiter(((spec.getRT(), point[0], point[1]) for spec in exp
                                  for point in zip(*spec.get_peaks())), dtype=[('RT', 'f'), ('mz', 'f'), ('inty', 'f')],
                                 count=cnt)
        np_stop = perf_counter()
        print("Time for loading and creating numpy array:",
              np_stop - load_stop)
        # spectradfwide = pd.DataFrame(data=((spec.getRT(), *spec.get_peaks()) for spec in exp), columns=cols)
        # spectradf = pd.DataFrame(data=((spec.getRT(), point[0], point[1]) for spec in exp for point in zip(*spec.get_peaks())), columns=expandcols)

        # Initial tests showed that loading into numpy array first is faster than direct construction from iter.
        spectradf = pd.DataFrame(data=spectraarr, columns=expandcols)
        # spectracds = ColumnDataSource(spectradf)
        maxrt = spectradf["RT"].max()
        minrt = spectradf["RT"].min()
        maxmz = spectradf["mz"].max()
        minmz = spectradf["mz"].min()

        # Stop the stopwatch / counter
        df_stop = perf_counter()
        print("Time for loading and creating DF:",
              df_stop - np_stop)

        def new_bounds_hook(plot, elem):
            x_range = plot.state.x_range
            y_range = plot.state.y_range
            x_range.bounds = minrt, maxrt
            y_range.bounds = minmz, maxmz

        points = hv.Points(spectradf, kdims=['RT', 'mz'], vdims=['inty'], label="MS1 survey scans").opts(
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
                              cnorm='log', alpha=10, min_alpha=0).opts(
            active_tools=['box_zoom'],
            tools=['hover'],
            hooks=[new_bounds_hook],
        ).opts(
            plot=dict(
                width=1000,
                height=1000,
                xlabel="Retention time (s)",
                ylabel="mass/charge (Da)"
            )
        )

        t2_stop = perf_counter()
        print("Time for creating plot objects:",
              t2_stop - df_stop)

        # Dynspread dynamically increases the size of the points when zooming in.
        finalplot = hd.dynspread(raster, threshold=0.7, how="add", shape="square")
        hvplot = renderer.get_plot(finalplot, doc)

        # updates the JS function to be triggered by invisText with the new filtered data.
        # then triggers a value change in invisText2. This will trigger a js side change for invisText,
        # eventually rendering the new data in 3D.
        # TODO with the new filtering mechanism, we can stop earlier if we see that the points are too many.
        #   We also do not need much of the binary search in the JS code anymore (and do not need to pass the RT range to JS)
        #   We could even call pyOpenMS functions on callback (e.g. extract range or noise filtering)
        def onbuttonclick():
            print("Foo")
            invisText.js_on_change("value", CustomJS(code=jsupdateinfo, args=dict(xr=hvplot.state.x_range,
                                                                                  yr=hvplot.state.y_range,
                                                                                  data=ColumnDataSource(
                                                                                      spectradf[spectradf['RT'].between(
                                                                                          hvplot.state.x_range.start,
                                                                                          hvplot.state.x_range.end)]))))
            invisText2.value = "baz"  # + str(random.randint(0, 10000))

        # The update 3D button
        bt = Button(label='Update 3D View')
        bt.on_click(onbuttonclick)
        dynamic_col.children = []
        dynamic_col.children.append(hvplot.state)
        dynamic_col.children.append(bt)

    @gen.coroutine
    def load_data(event):
            print("Starting load!")
            doc.add_timeout_callback(init_2D_plot_and_update_btn, 0)  # we do this instead

    doc.add_root(layout)
    doc.on_event("document_ready", load_data)
    return doc

def main():
    modify_doc(curdoc())

main()

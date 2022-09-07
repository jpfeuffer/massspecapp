import random
from functools import reduce, partial
from pathlib import Path
from time import perf_counter
import sys

import colorcet as cc
from numpy import linspace
from pandas import CategoricalDtype
from scipy.stats.kde import gaussian_kde

from bokeh.models import ColumnDataSource, FixedTicker, PrintfTickFormatter, FactorRange
from bokeh.plotting import figure
import bokeh.layouts
import numpy as np
from bokeh.models import ColumnDataSource, CustomJS, Button, Select, Div, MultiSelect
from bokeh.plotting import curdoc
from pyopenms import IdXMLFile, ProteinIdentification, PeptideIdentification, ControlledVocabulary, File, CVTerm, \
    DataValue
import pandas as pd
import os

spectradf = pd.DataFrame()

def modify_doc(doc):
    """Add a plotted function to the document.

    Arguments:
        doc: A bokeh document to which elements can be added.
    """

    prots = []  # type: list[ProteinIdentification]
    peps = []  # type: list[PeptideIdentification]

    if len(sys.argv) > 1:
        file = sys.argv[1]
    else:
        #file = "/Volumes/Data/iPRG2016/idXMLs/B1.idXML"
        #file = str(Path(__file__).resolve().parent) + "/static/data/BSA1_F1.idXML"
        file = "C:/git/OpenMS/src/tests/topp/MzTabExporter_2_input.idXML"

    switchDictT = {bool: np.bool_, int: np.intc, float: np.float, str: np.byte}
    # TODO find a heuristic for strings
    #  Especially sequences and spectrum_ids
    switchDict = {bool: '?', int: 'i', float: 'f', str: 'U100'}
    IdXMLFile().load(file, prots, peps)
    metavals = []
    types = []
    for pep in peps:
        hits = pep.getHits()
        if not len(hits) == 0:
            hits[0].getKeys(metavals)
            for k in metavals:
                if k == b"target_decoy":
                    types.append('?')
                else:
                    mv = hits[0].getMetaValue(k)
                    types.append(switchDict[type(mv)])
            break

    # TODO get score type name
    decodedMVs = [m.decode("utf-8") for m in metavals]
    cv = ControlledVocabulary()
    cv.loadFromOBO("psims", File.getOpenMSDataPath() + "/CV/psi-ms.obo")
    clearMVs = [cv.getTerm(m).name if m.startswith("MS:") else m for m in decodedMVs]
    cols = ["id", "RT", "mz", "score", "charge"] + decodedMVs
    clearcols = ["id", "RT", "mz", "score", "charge"] + clearMVs
    coltypes = ['S', 'f', 'f', 'f', 'i'] + types
    dt = list(zip(cols, coltypes))
    def extract(pep):
        hits = pep.getHits()
        if not hits:
            return tuple([pep.getIdentifier(), pep.getRT(), pep.getMZ(), np.NA, np.NA] + [np.NA]*len(metavals))
        else:
            besthit = hits[0]
            ret = [pep.getIdentifier(), pep.getRT(), pep.getMZ(), besthit.getScore(), besthit.getCharge()]
            for k in metavals:
                if besthit.metaValueExists(k):
                    val = besthit.getMetaValue(k)
                    if k == b"target_decoy":
                        if not val or val[0] == 't':
                            ret.append(True)
                        else:
                            ret.append(False)
                    else:
                        ret.append(val)
                else:
                    ret.append(np.nan)
            return tuple(ret)

    #TODO implement hasHits function in C++
    psmarr = np.fromiter((extract(pep) for pep in peps), dtype=dt, count=len(peps))
    #TODO make spectrum_ref the index, if available?
    psmdf = pd.DataFrame(psmarr)
    # rename cols
    psmdf.columns = clearcols
    psmdf['target_decoy'] = psmdf['target_decoy'].astype(CategoricalDtype(categories=[True, False])).cat.rename_categories(["target","decoy"])

    def ridge(category, data, scale=20):
        return list(zip([category] * len(data), scale * data))


    palette = [cc.rainbow[i * 15] for i in range(17)]
    plotcol = bokeh.layouts.column([], name="plotcol")

    def doPlotVal(attr, old, new):
        doPlot(new, selectCats.value)

    def doPlotCats(attr, old, new):
        doPlot(selectScore.value, new)

    def doPlot(scoretype, aggs):
        scoremin = psmdf[scoretype].min()
        scoremax = psmdf[scoretype].max()
        x = linspace(scoremin - 0.2, scoremax + 0.2, 300)
        source = ColumnDataSource(data=dict(x=x))
        gb = psmdf.groupby(aggs)
        nr_grps = len(gb)
        #TODO what if groups (i.e. nr of uniq vals in cat usually ints) are too many? Maybe check unique values before offering? Or fail here.
        #TODO maybe this can be combined with next loop?
        y_labs = []
        height = 900
        # we have to pad a lot, since the kernel dens. est. can have long tails.
        # which actually makes me think if histograms would be better or more exact.
        span = scoremax - scoremin
        p = figure(y_range=FactorRange(), plot_height=height, plot_width=700, x_range=(scoremin - span * 0.5,
                                                                                       scoremax + span * 0.5),
                   name="myfig")
        i = 0
        for name, grp in gb:
            print(name, grp, scoretype)
            if type(name) != str:
                try:
                    iter(name)
                    grpstr = "_".join([str(g) for g in name])
                except TypeError:
                    grpstr = str(name)
            else:
                grpstr = name
            y_labs.append(grpstr)
            p.y_range.factors = y_labs
            if len(grp) > 0:
                try:
                    # TODO play with bandwidth/estimation or let it be configurable
                    # TODO KDE only works with at least two points. What to do else?
                    pdf = gaussian_kde(grp[scoretype], 0.1)
                    p_x = pdf(x)
                    max_y = np.max(p_x)
                    print(max_y)
                    scale = 1/max_y
                    y = ridge(grpstr, pdf(x), scale/2)
                    source.add(y, grpstr)
                    p.patch('x', grpstr, color=palette[i], alpha=0.6, line_color="black", line_width=2, source=source)
                except np.linalg.LinAlgError:
                    print("Warning: sing. matrix")
            else:
                p.patch('x', grpstr, color=palette[i], alpha=0.6, line_color="black", line_width=2, source={grpstr: []})
            i+=1
        #p.y_range = y_labs
        p.y_range.range_padding = -0.2
        p.outline_line_color = None
        p.background_fill_color = "#efefef"

        #p.xaxis.ticker = FixedTicker(ticks=list(range(scoremin, scoremax+1, 10)))
        #p.xaxis.formatter = PrintfTickFormatter(format="%d%")

        p.ygrid.grid_line_color = None
        p.xgrid.grid_line_color = "#dddddd"
        p.xgrid.ticker = p.xaxis[0].ticker

        p.axis.minor_tick_line_color = None
        p.axis.major_tick_line_color = None
        p.axis.axis_line_color = None

        p.y_range.range_padding = 0 #0.12
        listOfSubLayouts = plotcol.children
        try:
            plotToRemove = doc.get_model_by_name('myfig')
            listOfSubLayouts.remove(plotToRemove)
        except Exception:
            pass
        listOfSubLayouts.append(p)

    # TODO load mass error, isotope error?
    scores = ["score"]
    for i in range(len(clearMVs)):
        if types[i] == 'f':
            scores.append(clearMVs[i])
    selectScore = Select(options=scores, value="score")
    selectScore.on_change("value", partial(doPlotVal))

    categories = ["charge"]
    for i in range(len(clearMVs)):
        if types[i] == 'U100' or types[i] == '?' or types[i] == 'i':
            categories.append(clearMVs[i])
    selectCats = MultiSelect(options=categories, value=["target_decoy"])
    selectCats.on_change("value", partial(doPlotCats))

    layout = bokeh.layouts.layout(bokeh.layouts.column([selectScore, selectCats]), plotcol)
    doc.add_root(layout)
    # this is done as initial change event of selectScore
    doPlot("score", ["target_decoy"])
    return doc

def main():
    modify_doc(curdoc())

main()

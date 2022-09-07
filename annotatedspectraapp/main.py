import rdkit.Chem
from bokeh.models import HoverTool, ColumnDataSource, BoxZoomTool
from rdkit.Chem.Draw import rdMolDraw2D
from bokeh.plotting import figure, show
import pandas as pd


def mol2svg(smiles):
    if smiles == "":
        return ""
    mol = rdkit.Chem.MolFromSmiles(smiles)
    d2d = rdMolDraw2D.MolDraw2DSVG(200, 100)
    d2d.DrawMolecule(mol)
    d2d.FinishDrawing()
    return d2d.GetDrawingText()


# prepare some data
mzs = [1., 2., 3.5, 4., 5.]
itys = [6., 7., 2., 4., 5.]
annotList = ["", "Cc1cc(-n2ncc(=O)[nH]c2=O)cc(C)c1C(O)c1ccc(Cl)cc1", "", "", ""]
annotList = [mol2svg(i) for i in annotList]
columns = ["mz", "ity", "annot"]

# create a new plot with a title and axis labels
p = figure(title="Simple line example", x_axis_label='m/z', y_axis_label='Intensity',
           tools=['xwheel_zoom', 'xpan'], active_scroll='xwheel_zoom')

# add a line renderer with legend and line thickness to the plot
for mz, ity, annot in zip(mzs, itys, annotList):
    source = ColumnDataSource(pd.DataFrame(
        {
            "mz": [mz, mz],
            "ity": [0, ity],
            "annot": ["", annot]
        }
    ))
    p.line(x="mz", y="ity", source=source, line_width=2, hover_line_width=3, color=("black" if annot == "" else "red"))

hover = HoverTool(
    tooltips="""
        <div>
        <div>
            @annot
        <div>
            <b>m/z:</b> <span style="font-size: 12px; color: #696;">@mz{0.000}</span><br>
            <b>Intensity:</b> <span style="font-size: 12px; color: #796;">@ity{0.000}</span>
        </div>
    </div>
    """,
)

hover.line_policy = 'next'
boxzoom = BoxZoomTool(dimensions="width")
boxzoom.overlay.fill_color = "grey"
boxzoom.overlay.fill_alpha = 0.1
boxzoom.overlay.line_width = 0.5
p.add_tools(hover)
p.add_tools(boxzoom)

# show the results
show(p)

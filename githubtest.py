import matplotlib.pyplot
import pyopenms.Plotting
from pyopenms import *

annotated_spectra = MSExperiment()
orig_spectra = MSExperiment()
loader = MzMLFile()
loader.load("/Users/pfeuffer/git/OpenMS-fixes-src/cmake-build-debug/src/tests/topp/spec.mzML", annotated_spectra)
loader.load("/Users/pfeuffer/git/OpenMS-fixes-src/src/tests/topp/THIRDPARTY/SiriusAdapter_1_input.mzML", orig_spectra)
spec = annotated_spectra[10]  # type: MSSpectrum
lookup = SpectrumLookup()
lookup.readSpectra(orig_spectra, ".?<SCAN>")
arrs = spec.getStringDataArrays()
foo = arrs[0]
print(spec.getNativeID())
orig = orig_spectra[lookup.findByNativeID(spec.getNativeID())]
pyopenms.Plotting.mirror_plot_spectrum(orig, spec, spectrum_top_kws={"color_ions": False},
                                       spectrum_bottom_kws={"color_ions": False})
matplotlib.pyplot.show()

from pyopenms import *
exp = MSExperiment()  # type: PeakMap
loader = MzMLFile()
loader.load("/Users/pfeuffer/PycharmProjects/pyopenms_test/bokehappfolder/static/data/BSA1.mzML", exp)
spec = exp[0]  # type: PeakSpectrum
metavaluekeys = []
spec.getKeys(metavaluekeys)
print(metavaluekeys) # metavaluekeys[0] == "base peak m/z"
print(spec.getMetaValue("base peak m/z"))
#for spec in exp:  # type: PeakSpectrum
#    instsettings = spec.getInstrumentSettings()  # type: InstrumentSettings
#    spec.metaRegistry()
#    if instsettings.getPolarity() == IonSource().Polarity.POSITIVE:
#        print("positive")
from time import perf_counter
from collections import defaultdict
from functools import reduce
from pathlib import Path
import sys
from pandas import CategoricalDtype
import numpy as np
from pyopenms import ConsensusXMLFile, ProteinIdentification, PeptideIdentification, ControlledVocabulary, File, CVTerm, \
    ConsensusMap, ConsensusFeature, FeatureHandle, ColumnHeader,PeakFileOptions, MzMLFile, PeakSpectrum, DRange1, DPosition1
import pandas as pd
import os
import queue
import threading


def SpectrumGenerator(file_path):
    q = queue.Queue(maxsize=1)
    JOB_DONE = object()
    TIMEOUT = 30

    def task():
        loader = MzMLFile()
        loadopts = loader.getOptions()  # type: PeakFileOptions
        loadopts.setMSLevels([1])
        loadopts.setSkipXMLChecks(True)
        loadopts.setIntensity32Bit(True)
        loadopts.setIntensityRange(DRange1(DPosition1(10000), DPosition1(sys.maxsize)))
        loader.setOptions(loadopts)
        loader.transform(file_path, MSCallback(q), True, True)
        q.put(JOB_DONE)

    t = threading.Thread(target=task)
    t.start()

    while True:
        # better set a timeout, or if task in sub-threading fails, it will result in a deadlock
        chunk = q.get(timeout=TIMEOUT)
        if chunk is JOB_DONE:
            break
        yield from chunk

    t.join()

class MSCallback:
    def __init__(self, q):
        self.q = q

    def setExperimentalSettings(self, s):
        pass

    def setExpectedSize(self, a, b):
        pass

    def consumeChromatogram(self, c):
        pass

    def consumeSpectrum(self, s: PeakSpectrum):
        #if s.getMSLevel() == 1:
            self.q.put((s.getRT(), point[0], point[1]) for point in zip(*s.get_peaks()))

class CntCallback:
    def __init__(self):
        self.cnt = 0

    def getCnt(self):
        return self.cnt

    def setExperimentalSettings(self, s):
        pass

    def setExpectedSize(self, a, b):
        pass

    def consumeChromatogram(self, c):
        pass

    def consumeSpectrum(self, s: PeakSpectrum):
            self.cnt += s.size()

def countMS1Peaks(file):
    cb = CntCallback()
    loader = MzMLFile()
    loadopts = loader.getOptions()  # type: PeakFileOptions
    loadopts.setMSLevels([1])
    loadopts.setIntensityRange(DRange1(DPosition1(10000), DPosition1(sys.maxsize)))
    loadopts.setSkipXMLChecks(True)
    loadopts.setFillData(True)
    loader.setOptions(loadopts)
    loader.transform(file, cb, True, True)
    return cb.getCnt()


if __name__ == '__main__':
    #"/Users/pfeuffer/git/OpenMS-fixes-src/share/OpenMS/examples/BSA/BSA1.mzML"
    start = perf_counter()
    file = "/Volumes/Data/UPS1/mzML/UPS1_250amol_R1.mzML"
    spectraarr = np.fromiter(SpectrumGenerator(file), dtype=[('RT', 'f'), ('mz', 'f'), ('inty', 'f')])
    print(perf_counter() - start)
    print(spectraarr[0])
    exit(0)

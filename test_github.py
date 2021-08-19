from pyopenms import *

if __name__ == '__main__':
    edf = ExperimentalDesignFile().load("/Users/pfeuffer/Downloads/MSV000084264_ExperimentalDesign.tsv".encode(),False)
    #Prepare data loading (save memory by only
    # loading MS1 spectra into memory)
    fh = MzMLFile()

    # Load data
    input_map = MSExperiment()
    fh.load("/Users/pfeuffer/Downloads/180308_Q_QC1X_01_01_ORIGINAL.mzML", input_map)
    input_map.updateRanges()

    ic = InternalCalibration()
    lm = [InternalCalibration_LockMass(100, 1, 2)]
    cd = CalibrationData()
    ic.fillCalibrants(input_map, lm, 0.1, True, False, cd, False)



    exit(0)
    ff = FeatureFinder()
    ff.setLogType(LogType.CMD)
    #ff.setLogType(LogType.GUI)

    # Run the feature finder
    name = "centroided"
    features = FeatureMap()
    seeds = FeatureMap()
    params = ff.getParameters(name)

    for k, v in params.asDict().items():
        print(k, v)

    params.setValue("mass_trace:min_spectra", 5)

    for k, v in params.asDict().items():
        print(k, v)

    ff.run(name, input_map, features, params, seeds)

    features.setUniqueIds()
    fh = FeatureXMLFile()
    fh.store("test.featureXML", features)
    print("Found", features.size(), "features")

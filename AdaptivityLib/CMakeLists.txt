# @file CMakeList.txt  -- First file to read for CMake

#To include hp-refinement strategies
	add_subdirectory(StrategyRef)

#To include error estimators
	add_subdirectory(ErrorEstimators)

#To include computing wavelet coefficients
	add_subdirectory(WaveletCoef)

target_sources(OPS_Reliability
    PRIVATE
        BetaRV.cpp
        ChiSquareRV.cpp
        ExponentialRV.cpp
        GammaRV.cpp
        GumbelRV.cpp
        LaplaceRV.cpp
        LognormalRV.cpp
        NormalRV.cpp
        ParetoRV.cpp
	#  PythonRV.cpp
        RayleighRV.cpp
        ShiftedExponentialRV.cpp
        ShiftedRayleighRV.cpp
        Type1LargestValueRV.cpp
        Type1SmallestValueRV.cpp
        Type2LargestValueRV.cpp
        Type3SmallestValueRV.cpp
        UniformRV.cpp
        UserDefinedRV.cpp
        WeibullRV.cpp
    PUBLIC
        BetaRV.h
        ChiSquareRV.h
        ExponentialRV.h
        GammaRV.h
        GumbelRV.h
        LaplaceRV.h
        LognormalRV.h
        NormalRV.h
        ParetoRV.h
	# PythonRV.h
        RayleighRV.h
        ShiftedExponentialRV.h
        ShiftedRayleighRV.h
        Type1LargestValueRV.h
        Type1SmallestValueRV.h
        Type2LargestValueRV.h
        Type3SmallestValueRV.h
        UniformRV.h
        UserDefinedRV.h
        WeibullRV.h
)

target_include_directories(OPS_Reliability PUBLIC $(CMAKE_CURRENT_LIST_DIR))
target_include_directories(OpenSeesPy PUBLIC $(CMAKE_CURRENT_LIST_DIR))

target_sources(OpenSeesPy
    PRIVATE
	PythonRV.cpp
    PUBLIC
        PythonRV.h
)

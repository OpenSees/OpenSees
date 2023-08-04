target_sources(OPS_Reliability
    PRIVATE
        JonswapSpectrum.cpp
        NarrowBandSpectrum.cpp
        PointsSpectrum.cpp
        Spectrum.cpp
        SpectrumIter.cpp
    PUBLIC
        JonswapSpectrum.h
        NarrowBandSpectrum.h
        PointsSpectrum.h
        Spectrum.h
        SpectrumIter.h
)
target_include_directories(OPS_Reliability PUBLIC ${CMAKE_CURRENT_LIST_DIR})


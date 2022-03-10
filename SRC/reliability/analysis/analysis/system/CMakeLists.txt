target_sources(OPS_Reliability
    PRIVATE
        IPCM.cpp
        MVNcdf.cpp
        PCM.cpp
        SCIS.cpp
    PUBLIC
        IPCM.h
        MVNcdf.h
        PCM.h
        SCIS.h
)
target_include_directories(OPS_Reliability PUBLIC $(CMAKE_CURRENT_LIST_DIR))


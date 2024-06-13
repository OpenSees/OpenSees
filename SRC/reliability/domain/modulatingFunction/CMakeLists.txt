target_sources(OPS_Reliability
    PRIVATE
        ConstantModulatingFunction.cpp
        GammaModulatingFunction.cpp
        KooModulatingFunction.cpp
        ModulatingFunction.cpp
        ModulatingFunctionIter.cpp
        TrapezoidalModulatingFunction.cpp
    PUBLIC
        ConstantModulatingFunction.h
        GammaModulatingFunction.h
        KooModulatingFunction.h
        ModulatingFunction.h
        ModulatingFunctionIter.h
        TrapezoidalModulatingFunction.h
)
target_include_directories(OPS_Renderer PUBLIC ${CMAKE_CURRENT_LIST_DIR})


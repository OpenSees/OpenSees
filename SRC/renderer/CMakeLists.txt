
target_sources(OPS_Renderer
    PRIVATE
        AGL_Device.cpp
        Clipping.cpp
        db.cpp
        Device.cpp
        DofColorMap.cpp
        gMatrix.cpp
        main.cpp
        mainTriangle.cpp
        OpenGlDevice.cpp
        OpenGlRenderer.cpp
        PlainMap.cpp
        Projection.cpp
        Renderer.cpp
        Scan.cpp
        simple.c
        View.cpp
        Viewport.cpp
        VrmlViewer.cpp
        WindowDevice.cpp
        WindowRenderer.cpp
        X11Device.cpp
        X11Renderer.cpp
    PUBLIC
        AGL_Device.h
        Clipping.h
        ColorMap.h
        Device.h
        DofColorMap.h
        OpenGlDevice.h
        OpenGLRenderer.h
        PlainMap.h
        Projection.h
        Renderer.h
        Scan.h
        View.h
        Viewport.h
        VrmlViewer.h
        WindowDevice.h
        WindowRenderer.h
        X11Device.h
        X11Renderer.h
)

target_include_directories(OPS_Renderer PUBLIC $(CMAKE_CURRENT_LIST_DIR))


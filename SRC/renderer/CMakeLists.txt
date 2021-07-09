target_sources(OPS_Renderer
    PRIVATE
        #AGL_Device.cpp
        Clipping.cpp
        db.cpp
        Device.cpp
        DofColorMap.cpp
        gMatrix.cpp
	#OpenGlDevice.cpp
	#OpenGlRenderer.cpp
        PlainMap.cpp
        Projection.cpp
        Renderer.cpp
        Scan.cpp
        View.cpp
        Viewport.cpp
        VrmlViewer.cpp
        WindowDevice.cpp
        WindowRenderer.cpp
    PUBLIC
        #AGL_Device.h
        Clipping.h
        ColorMap.h
        Device.h
        DofColorMap.h
        OpenGlDevice.h
	#OpenGLRenderer.h
        PlainMap.h
        Projection.h
        Renderer.h
        Scan.h
        View.h
        Viewport.h
        VrmlViewer.h
        WindowDevice.h
        WindowRenderer.h
)
add_library(OPS_Renderer_X11 OBJECT
        X11Device.cpp
        X11Renderer.cpp
        #X11Device.h
        #X11Renderer.h
)

target_include_directories(OPS_Renderer PUBLIC $(CMAKE_CURRENT_LIST_DIR))


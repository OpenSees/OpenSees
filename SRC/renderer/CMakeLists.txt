#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Renderer
  PRIVATE
    PlainMap.cpp
    Renderer.cpp
  PUBLIC
    PlainMap.h
    Renderer.h
    #ColorMap.h
)

target_sources(OPS_Graphics
  PRIVATE
    DofColorMap.cpp
  PUBLIC
    DofColorMap.h
)

target_sources(OPS_Graphics_Default
  PRIVATE
    Clipping.cpp
    db.cpp
    Device.cpp
    gMatrix.cpp
    Projection.cpp
    Scan.h
    View.h
    Viewport.cpp
    VrmlViewer.cpp
    WindowDevice.cpp
    WindowRenderer.cpp
    X11Device.cpp
    X11Renderer.cpp
  PUBLIC
    Clipping.h
    Device.h
    Projection.h
    Scan.cpp
    View.cpp
    Viewport.h
    VrmlViewer.h
    WindowDevice.h
    WindowRenderer.h
    X11Device.h
    X11Renderer.h
)

target_sources(OPS_Graphics_GL
  PRIVATE
    AGL_Device.cpp
    OpenGlDevice.cpp
    OpenGlRenderer.cpp
  PUBLIC
    AGL_Device.h
    OpenGlDevice.h
    OpenGLRenderer.h
)

target_include_directories(OPS_Renderer PUBLIC ${CMAKE_CURRENT_LIST_DIR})


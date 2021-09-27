#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Renderer
  PRIVATE
    PlainMap.cpp   # currently added in libOpenSees
    #AGL_Device.cpp
    Clipping.cpp
    db.cpp
    Device.cpp
    DofColorMap.cpp
    gMatrix.cpp
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

add_library(OPS_Renderer_X11 OBJECT EXCLUDE_FROM_ALL)
target_sources(OPS_Renderer_X11
  PRIVATE
    X11Device.cpp
    X11Renderer.cpp
  PUBLIC
    X11Device.h
    X11Renderer.h
)

add_library(OPS_Renderer_OpenGL OBJECT EXCLUDE_FROM_ALL
    OpenGlDevice.cpp
    OpenGlRenderer.cpp
)

target_include_directories(OPS_Renderer PUBLIC ${CMAKE_CURRENT_LIST_DIR})


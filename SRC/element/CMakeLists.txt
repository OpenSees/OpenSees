#add_library(OPS_Element OBJECT)

target_sources(G3_Element
    PRIVATE
      Element.cpp
      ElementalLoad.cpp
      #Information.cpp
    PUBLIC
      Element.h
      ElementalLoad.h
      #Information.h
)

target_sources(G3_Utilities
    PRIVATE
      Information.cpp
    PUBLIC
      Information.h
)

#target_sources(OPS_Tcl
#    PRIVATE
#        TclElementCommands.cpp
#    #PUBLIC
#)


add_subdirectory(zeroLength)

add_subdirectory(brick)

add_subdirectory(truss)

add_subdirectory(beam2d)
add_subdirectory(beam3d)
#add_subdirectory(beamWithHinges)

add_subdirectory(catenaryCable)
add_subdirectory(LHMYS)
add_subdirectory(PML)
add_subdirectory(RockingBC)
add_subdirectory(absorbentBoundaries)
add_subdirectory(adapter)
add_subdirectory(componentElement)

add_subdirectory(dispBeamColumn)
add_subdirectory(dispBeamColumnInt)
add_subdirectory(forceBeamColumn)
add_subdirectory(mixedBeamColumn)
add_subdirectory(elasticBeamColumn)
add_subdirectory(gradientInelasticBeamColumn)
add_subdirectory(updatedLagrangianBeamColumn)

add_subdirectory(fourNodeQuad)

add_subdirectory(frictionBearing)
add_subdirectory(elastomericBearing)
add_subdirectory(feap)
add_subdirectory(generic)
add_subdirectory(joint)
add_subdirectory(mvlem)
add_subdirectory(shell)
add_subdirectory(surfaceLoad)
add_subdirectory(PFEMElement)

add_subdirectory(pyMacro)
add_subdirectory(UP-ucsd)
add_subdirectory(UWelements)
add_subdirectory(HUelements)
add_subdirectory(XMUelements)

target_include_directories(OPS_Element PUBLIC $(CMAKE_CURRENT_LIST_DIR))


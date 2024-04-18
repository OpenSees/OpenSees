#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#add_library(OPS_Element OBJECT)

target_sources(OPS_Element
    PRIVATE
      Element.cpp
      ElementalLoad.cpp
      WrapperElement.cpp
      #Information.cpp
    PUBLIC
      Element.h
      ElementalLoad.h
      WrapperElement.h
      #Information.h
)

target_sources(OPS_Utilities
    PRIVATE
      Information.cpp
    PUBLIC
      Information.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

add_subdirectory(fourNodeQuad)
add_subdirectory(twoNodeLink)
add_subdirectory(tetrahedron)
add_subdirectory(triangle)
add_subdirectory(shell)
add_subdirectory(brick)


add_subdirectory(zeroLength)

add_subdirectory(truss)
add_subdirectory(beam2d)
add_subdirectory(beam3d)
#add_subdirectory(beamWithHinges)

add_subdirectory(dispBeamColumn)
add_subdirectory(dispBeamColumnInt)
add_subdirectory(forceBeamColumn)
add_subdirectory(mixedBeamColumn)
add_subdirectory(elasticBeamColumn)
add_subdirectory(gradientInelasticBeamColumn)
add_subdirectory(updatedLagrangianBeamColumn)
add_subdirectory(nonlinearBeamColumn)

add_subdirectory(catenaryCable)
add_subdirectory(LHMYS)
add_subdirectory(PML)
add_subdirectory(RockingBC)
add_subdirectory(absorbentBoundaries)
add_subdirectory(adapter)
add_subdirectory(componentElement)

add_subdirectory(frictionBearing)
add_subdirectory(elastomericBearing)
add_subdirectory(generic)
add_subdirectory(joint)
add_subdirectory(surfaceLoad)

add_subdirectory(feap)

add_subdirectory(PFEMElement)

add_subdirectory(mvlem)
add_subdirectory(pyMacro)
add_subdirectory(UP-ucsd)
add_subdirectory(UWelements)
add_subdirectory(HUelements)
add_subdirectory(XMUelements)
add_subdirectory(masonry)
add_subdirectory(CEqElement)

add_subdirectory(IGA)

add_subdirectory(mefi)

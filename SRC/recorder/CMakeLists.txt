#add_library(recorder OBJECT)

target_sources(G3_Recorder
    PRIVATE
      DamageRecorder.cpp
      DatastoreRecorder.cpp
      DriftRecorder.cpp
      ElementRecorder.cpp
      ElementRecorderRMS.cpp
      EnvelopeDriftRecorder.cpp
      EnvelopeElementRecorder.cpp
      EnvelopeNodeRecorder.cpp
      GSA_Recorder.cpp
      GmshRecorder.cpp
      MaxNodeDispRecorder.cpp
      NodeRecorder.cpp
      NodeRecorderRMS.cpp
      NormElementRecorder.cpp
      NormEnvelopeElementRecorder.cpp
      PatternRecorder.cpp
      Recorder.cpp
      RemoveRecorder.cpp
      VTK_Recorder.cpp
    PUBLIC
      DamageRecorder.h
      DatastoreRecorder.h
      DriftRecorder.h
      ElementRecorder.h
      ElementRecorderRMS.h
      EnvelopeDriftRecorder.h
      EnvelopeElementRecorder.h
      EnvelopeNodeRecorder.h
      GSA_Recorder.h
      GmshRecorder.h
      MaxNodeDispRecorder.h
      NodeRecorder.h
      NodeRecorderRMS.h
      NormElementRecorder.h
      NormEnvelopeElementRecorder.h
      PatternRecorder.h
      Recorder.h
      RemoveRecorder.h
      VTK_Recorder.h
)
target_sources(OPS_Paraview
    PRIVATE
      PVDRecorder.cpp
    PUBLIC
      PVDRecorder.h
)

target_sources(OPS_ASDEA
    PRIVATE
      MPCORecorder.cpp
    PUBLIC
      MPCORecorder.h
)


#if (OPS_Use_Graphics)
target_sources(OPS_Renderer
    PRIVATE
      AlgorithmIncrements.cpp
      FilePlotter.cpp
    PUBLIC
      AlgorithmIncrements.h
      FilePlotter.h
)
#endif()
add_subdirectory(response)
#target_include_directories(recorder PUBLIC $(CMAKE_CURRENT_LIST_DIR))


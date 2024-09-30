#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#

target_sources(OPS_Runtime PRIVATE
    "commands.cpp"
    "domain.cpp"
    "element.cpp"
    "response.cpp"
    "region.cpp"
    "nodes.cpp"
    "runtime.cpp"
    "rigid_links.cpp"
    "recorder.cpp"
    "TclUpdateMaterialStageCommand.cpp"
    "TclUpdateMaterialCommand.cpp"
    "parameter.cpp"
    "sensitivity.cpp"

# LOADS & PATTERNS
    "loading/groundMotion.cpp"
    "loading/element_load.cpp"
    "loading/pattern.cpp"
    "loading/series.cpp"
    "loading/TclSeriesIntegratorCommand.cpp"
    #"domain/pattern/drm/TclPatternCommand.cpp"
)

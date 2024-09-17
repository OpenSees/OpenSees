#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#

target_sources(OPS_Runtime PRIVATE
# Interpreter
    "interpreter.cpp"
    "pragma.cpp"
    "packages.cpp"
    "parallel/sequential.cpp"

# Modeling
    "modeling/model.cpp"
    "modeling/nodes.cpp"
    "modeling/constraint.cpp"
    "modeling/geomTransf.cpp"
    "modeling/element.cpp"
    "modeling/nDMaterial.cpp"
    "modeling/section.cpp"
    "modeling/uniaxialMaterial.cpp"
    "modeling/uniaxial.cpp"
    "modeling/printing.cpp"
    "modeling/blockND.cpp"
    "modeling/Block2D.cpp"
    "modeling/Block3D.cpp"
    "modeling/invoking/invoke.cpp"
    "modeling/invoking/invoke_uniaxial.cpp"
    "modeling/invoking/invoke_section.cpp"
    "modeling/invoking/invoke_stress.cpp"

# Analysis
    "analysis/algorithm.cpp"
    "analysis/integrator.cpp"
    "analysis/transient.cpp"
    "analysis/analysis.cpp"
    "analysis/numberer.cpp"
    "analysis/ctest.cpp"
    "analysis/solver.cpp"
    "analysis/solver.hpp"

# Utilities
    "utilities/utilities.cpp"
    "utilities/progress.cpp"
    "utilities/formats.cpp"
)

add_subdirectory(domain)


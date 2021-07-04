# Claudio Perez
function (opensees_add_cxx_flag)
    cmake_parse_arguments(
        PARSE_ARGV 0
        OPS_COMPILE_OPT_ARG # prefix of output variables
        "PUBLIC;PRIVATE;INTERFACE" # list of names of the boolean arguments (only defined ones will be true)
        "" # list of names of mono-valued arguments
        "TARGETS;MSVC;GNU" # list of names of multi-valued arguments (output variables are lists)
        #${ARGN} # arguments of the function to parse, here we take the all original ones
    )
    set(flg "PUBLIC")
    if (NOT ${OPS_COMPILE_OPT_ARG_TARGETS})
        add_compile_options(${flg} $<$<COMPILE_LANGUAGE:CXX>:${${OPS_COMPILE_OPT_ARG_}${CXX_COMPILER_ID}}>)
    endif()
    foreach(src ${OPS_COMPILE_OPT_ARG_TARGETS})
        if(MSVC)
            foreach(trgt ${OPS_COMPILE_OPT_ARG_TARGETS})
                    target_compile_options(${trgt} ${flg} $<$<COMPILE_LANGUAGE:CXX>:${OPS_COMPILE_OPT_ARG_MSVC}>)
            endforeach()
        else()
            foreach(trgt ${OPS_COMPILE_OPT_ARG_TARGETS})
                    target_compile_options(${trgt} ${flg} $<$<COMPILE_LANGUAGE:CXX>:${OPS_COMPILE_OPT_ARG_GNU}>)
            endforeach()
        endif()
    endforeach(src)
endfunction()


function (opensees_load lib_name)
    cmake_parse_arguments(
        PARSE_ARGV 1
        OPS_LOAD_ARG # prefix of output variables
        "SEARCH;BUILD;FIND" # list of names of the boolean arguments (only defined ones will be true)
	"LIBRARY;INCLUDE;VENDOR" # list of names of mono-valued arguments
        "PATHS" # list of names of multi-valued arguments (output variables are lists)
        #${ARGN} # arguments of the function to parse, here we take the all original ones
    )
    set(OPS_PKG_FOUND_VAR "${lib_name}_FOUND")# PARENT_SCOPE) 
    set(${OPS_PKG_FOUND_VAR} FALSE)# PARENT_SCOPE)
    if(OPS_LOAD_ARG_BUILD)
	message("OPS >>> Using OpenSees bundled library '${lib_name}'")
        set(${OPS_PKG_FOUND_VAR} TRUE PARENT_SCOPE)
        opensees_build(${lib_name})
        return()
    elseif(OPS_LOAD_ARG_LIBRARY)
	    message("OPS >>> ${lib_name}")
	    set("${lib_name}_LIBRARIES"    ${OPS_LOAD_ARG_LIBRARY} PARENT_SCOPE)
	    set("${lib_name}_INCLUDE_DIRS" ${OPS_LOAD_ARG_INCLUDE} PARENT_SCOPE)
    elseif(OPS_LOAD_ARG_PATHS)
        message("Provided ${lib_name} paths are:")
        foreach(src ${OPS_LOAD_ARG_PATHS})
            message("- ${src}")
        endforeach(src)
        find_package(${lib_name} PATHS ${OPS_LOAD_ARG_PATHS})
    else()
        message("OPS>>> find_package(${lib_name})")
        find_package(${lib_name})
        message("OPS>>> status: ${OPS_PKG_FOUND_VAR}=${${OPS_PKG_FOUND_VAR}}")
    endif()

    if(OPS_LOAD_ARG_SEARCH)
        if(NOT ${${OPS_PKG_FOUND_VAR}})
            message("CMake failed to find ${lib_name}; building OpenSees version")
            set(${OPS_PKG_FOUND_VAR} TRUE PARENT_SCOPE)
            opensees_build(${lib_name})
            set("${lib_name}_LIBRARIES" ${lib_name} PARENT_SCOPE)
        endif()
    endif()
endfunction()

function (opensees_build lib_name)
    add_subdirectory("${OPS_EXTERNALS_DIR}/${lib_name}")
    include_directories("${OPS_EXTERNALS_DIR}/${lib_name}") 
    message(${${lib_name}_LIBRARIES})
endfunction()





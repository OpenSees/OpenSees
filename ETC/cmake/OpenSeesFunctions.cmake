# Claudio Perez
#include(conan.cmake)

function (opensees_library elemlib)
  # opensees_library(<lib_name> [REQUIRES <requirement>] <sources>...)
  cmake_parse_arguments(PARSE_ARGV 1 OPS_LIB_ARG "EXTERNAL" "" "LINK;SOURCES")
  if (OPS_LIB_ARG_EXTERNAL)
    set(OPS_LIB_NEW_TYP "SHARED")
  endif()
  add_library(${elemlib} EXCLUDE_FROM_ALL ${OPS_LIB_NEW_TYP})
  set(private_sources ${OPS_LIB_ARG_SOURCES})
  set(public_sources  ${OPS_LIB_ARG_SOURCES})
  list(FILTER private_sources INCLUDE REGEX ".*\.cpp")
  list(FILTER public_sources INCLUDE REGEX ".*\.h")
  target_sources(${elemlib} PRIVATE ${private_sources} PUBLIC ${public_sources})

  if (OPS_LIB_ARG_LINK)
    target_link_libraries(${elemlib} INTERFACE ${OPS_LIB_ARG_LINK}) 
  endif()
  set_target_properties(${elemlib} PROPERTIES LINKER_LANGUAGE CXX)
  #foreach(lib ${OPS_LIB_ARG_REQUIRE_NUMLIB})
  #  list(APPEND OPS_Numlib_List ${lib})
  #endforeach()
  #set(OPS_Numlib_List ${OPS_Numlib_List} PARENT_SCOPE)
endfunction()

function (opensees_elements elemlib)
  opensees_library(${elemlib} ${ARGN})
endfunction()

function (opensees_uniaxials unilib)
  opensees_library(${unilib} ${ARGN})
endfunction()


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
        add_compile_options($<$<COMPILE_LANGUAGE:CXX>:${${OPS_COMPILE_OPT_ARG_}${CXX_COMPILER_ID}}>)
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
#
#
#
    cmake_parse_arguments(
      PARSE_ARGV 1
      OPS_LOAD_ARG # prefix of output variables
      "BUILD;FIND" # boolean arguments (only defined ones will be true)
      "BUNDLED;LIBRARY;CONAN;INCLUDE" # mono-valued arguments
      "PATHS;AS"   # multi-valued arguments (output variables are lists)
    )
    set(OPS_PKG_FOUND_VAR "${lib_name}_FOUND")
    set(${OPS_PKG_FOUND_VAR} FALSE)

    if(OPS_LOAD_ARG_CONAN)
      if(NOT ${${OPS_PKG_FOUND_VAR}})
        #message(">>> >>> ${OPS_LOAD_ARG_CONAN}")
        conan_cmake_configure(REQUIRES ${OPS_LOAD_ARG_CONAN}
          #BASIC_SETUP CMAKE_TARGETS 
          #BUILD missing
        )
        conan_cmake_install(PATH_OR_REFERENCE . BUILD missing )#REMOTE conan-center)
      endif()
    endif()

    if(OPS_LOAD_ARG_LIBRARY)
      if(NOT ${${OPS_PKG_FOUND_VAR}})
        message("OPS >>> ${lib_name}")
        set("${lib_name}_LIBRARIES"    ${OPS_LOAD_ARG_LIBRARY})
        set("${lib_name}_INCLUDE_DIRS" ${OPS_LOAD_ARG_INCLUDE})
      endif()
    endif()

    if(OPS_LOAD_ARG_PATHS)
      message("Provided ${lib_name} paths are:")
      foreach(src ${OPS_LOAD_ARG_PATHS})
        message("- ${src}")
      endforeach(src)
      find_package(${lib_name} PATHS ${OPS_LOAD_ARG_PATHS})
    endif()

    if(OPS_LOAD_ARG_FIND)
        if(NOT ${${OPS_PKG_FOUND_VAR}})
          message("OPS >>> find_package(${lib_name})")
          find_package(${lib_name})
          message("OPS >>> status: ${OPS_PKG_FOUND_VAR}=${${OPS_PKG_FOUND_VAR}}")
        endif()
    endif()

    if(OPS_LOAD_ARG_BUNDLED)
      if(NOT ${${OPS_PKG_FOUND_VAR}})
        #message("OPS >>> Building OpenSees bundled ${lib_name}")
        set(${OPS_PKG_FOUND_VAR} TRUE PARENT_SCOPE)
        add_subdirectory("${OPS_LOAD_ARG_BUNDLED}")
        include_directories("${OPS_LOAD_ARG_BUNDLED}")
        set("${lib_name}_LIBRARIES" "${lib_name}") 
      endif()
    endif()

    set("${lib_name}_LIBRARIES"    ${${lib_name}_LIBRARIES}    PARENT_SCOPE)
    set("${lib_name}_INCLUDE_DIRS" ${${lib_name}_INCLUDE_DIRS} PARENT_SCOPE)
    target_compile_definitions(OPS_External_packages INTERFACE "OPSDEF_${lib_name}")
    #message("    status: ${lib_name}_LIBRARIES    =${${lib_name}_LIBRARIES}")
    #message("    status: ${lib_name}_INCLUDE_DIRS =${${lib_name}_INCLUDE_DIRS}\n")
endfunction()

function (opensees_build lib_name)
    add_subdirectory("${OPS_BUNDLED_DIR}/${lib_name}")
    include_directories("${OPS_BUNDLED_DIR}/${lib_name}") 
    message(${${lib_name}_LIBRARIES})
endfunction()





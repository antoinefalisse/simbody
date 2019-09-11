## This whole directory exists just so I could define this extra preprocessor value.

add_definitions(-DSimTK_SIMMATH_BUILDING_SHARED_LIBRARY)

#
# Set up file groups for better browsing in Visual Studio.
#

# Override Api headers to put them in the Api group.
source_group("Header Files\\Api" FILES ${API_ABS_INCLUDE_FILES})

# Now put things in their specified groups.
if(SOURCE_GROUPS)
    list(LENGTH SOURCE_GROUPS NGROUPS)
    math(EXPR lastgrpnum ${NGROUPS}-1)
    foreach( grpnum RANGE 0 ${lastgrpnum} )
        list(GET SOURCE_GROUPS ${grpnum} grp)
        list(GET SOURCE_GROUP_FILES ${grpnum} grpfile)
        source_group("${grp}" FILES "${grpfile}")
    endforeach()
endif()

# These are defaults for files that don't otherwise get mentioned.
source_group("Source Files\\Others" FILES ${SOURCE_FILES})
source_group("Header Files\\Others" FILES ${SOURCE_INCLUDE_FILES})

# Do not use IpOpt with Recorder: c/cpp
set(SOURCE_FILES_RECORDER ${SOURCE_FILES})
file(GLOB ipopt_files_to_remove INCLUDE_DIRECTORIES false ${CMAKE_SOURCE_DIR}/SimTKmath/Optimizers/src/IpOpt/*.cpp)
file(GLOB cpodes_files_to_remove1 INCLUDE_DIRECTORIES false ${CMAKE_SOURCE_DIR}/SimTKmath/Integrators/src/CPodes/*.cpp)
file(GLOB cpodes_files_to_remove2 INCLUDE_DIRECTORIES false ${CMAKE_SOURCE_DIR}/SimTKmath/Integrators/src/CPodes/sundials/src/cpodes/*.c)
file(GLOB cpodes_files_to_remove3 INCLUDE_DIRECTORIES false ${CMAKE_SOURCE_DIR}/SimTKmath/Integrators/src/CPodes/sundials/src/nvec_ser/*.c)
file(GLOB cpodes_files_to_remove4 INCLUDE_DIRECTORIES false ${CMAKE_SOURCE_DIR}/SimTKmath/Integrators/src/CPodes/sundials/src/sundials/*.c)
file(GLOB cmaes_files_to_remove INCLUDE_DIRECTORIES false ${CMAKE_SOURCE_DIR}/SimTKmath/Optimizers/src/c-cmaes/*.c)
file(GLOB optimizers_files_to_remove INCLUDE_DIRECTORIES false ${CMAKE_SOURCE_DIR}/SimTKmath/Optimizers/src/*.cpp)
list(REMOVE_ITEM SOURCE_FILES_RECORDER ${ipopt_files_to_remove})
list(REMOVE_ITEM SOURCE_FILES_RECORDER ${cpodes_files_to_remove1})
list(REMOVE_ITEM SOURCE_FILES_RECORDER ${cpodes_files_to_remove2})
list(REMOVE_ITEM SOURCE_FILES_RECORDER ${cpodes_files_to_remove3})
list(REMOVE_ITEM SOURCE_FILES_RECORDER ${cpodes_files_to_remove4})
list(REMOVE_ITEM SOURCE_FILES_RECORDER ${cmaes_files_to_remove})
list(REMOVE_ITEM SOURCE_FILES_RECORDER ${optimizers_files_to_remove})

# Do not use IpOpt with Recorder: h
set(SOURCE_INCLUDE_FILES_RECORDER ${SOURCE_INCLUDE_FILES})
file(GLOB ipopt_files_to_remove_h INCLUDE_DIRECTORIES false ${CMAKE_SOURCE_DIR}/SimTKmath/Optimizers/src/IpOpt/*.h)
file(GLOB ipopt_files_to_remove_hpp INCLUDE_DIRECTORIES false ${CMAKE_SOURCE_DIR}/SimTKmath/Optimizers/src/IpOpt/*.hpp)
file(GLOB cpodes_files_to_remove1_h INCLUDE_DIRECTORIES false ${CMAKE_SOURCE_DIR}/SimTKmath/Integrators/src/CPodes/*.h)
file(GLOB cpodes_files_to_remove2_h INCLUDE_DIRECTORIES false ${CMAKE_SOURCE_DIR}/SimTKmath/Integrators/src/CPodes/sundials/include/cpodes/*.h)
file(GLOB cpodes_files_to_remove3_h INCLUDE_DIRECTORIES false ${CMAKE_SOURCE_DIR}/SimTKmath/Integrators/src/CPodes/sundials/include/nvector/*.h)
file(GLOB cpodes_files_to_remove4_h INCLUDE_DIRECTORIES false ${CMAKE_SOURCE_DIR}/SimTKmath/Integrators/src/CPodes/sundials/include/sundials/*.h)
file(GLOB cmaes_files_to_remove_h INCLUDE_DIRECTORIES false ${CMAKE_SOURCE_DIR}/SimTKmath/Optimizers/src/c-cmaes/*.h)
file(GLOB optimizers_files_to_remove_h INCLUDE_DIRECTORIES false ${CMAKE_SOURCE_DIR}/SimTKmath/Optimizers/src/*.h)
list(REMOVE_ITEM SOURCE_INCLUDE_FILES_RECORDER ${ipopt_files_to_remove_h})
list(REMOVE_ITEM SOURCE_INCLUDE_FILES_RECORDER ${ipopt_files_to_remove_hpp})
list(REMOVE_ITEM SOURCE_INCLUDE_FILES_RECORDER ${cpodes_files_to_remove1_h})
list(REMOVE_ITEM SOURCE_INCLUDE_FILES_RECORDER ${cpodes_files_to_remove2_h})
list(REMOVE_ITEM SOURCE_INCLUDE_FILES_RECORDER ${cpodes_files_to_remove3_h})
list(REMOVE_ITEM SOURCE_INCLUDE_FILES_RECORDER ${cpodes_files_to_remove4_h})
list(REMOVE_ITEM SOURCE_INCLUDE_FILES_RECORDER ${cmaes_files_to_remove_h})
list(REMOVE_ITEM SOURCE_INCLUDE_FILES_RECORDER ${optimizers_files_to_remove_h})

if(BUILD_UNVERSIONED_LIBRARIES)

    add_library(${RECORDER_TARGET} SHARED 
                ${SOURCE_FILES_RECORDER} ${SOURCE_INCLUDE_FILES_RECORDER} 
                ${API_ABS_INCLUDE_FILES})
    
    target_link_libraries(${RECORDER_TARGET} 
                 ${SimTKCOMMON_RECORDER_LIBRARY}
                 ${MATH_LIBS_TO_USE}
				 ${RECORDER_LIBRARIES})
    
    set_target_properties(${RECORDER_TARGET} PROPERTIES
        PROJECT_LABEL "Code - ${RECORDER_TARGET}"
        SOVERSION ${SIMBODY_SONAME_VERSION})
    
    # install library; on Windows both .lib and .dll go in the lib directory.
    install(TARGETS ${RECORDER_TARGET} EXPORT SimbodyTargets
                     PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE 
                     GROUP_READ GROUP_WRITE GROUP_EXECUTE 
                     WORLD_READ WORLD_EXECUTE
                     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                     ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
                     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

endif(BUILD_UNVERSIONED_LIBRARIES)


if(BUILD_VERSIONED_LIBRARIES)

    add_library(${RECORDER_TARGET_VN} SHARED 
                ${SOURCE_FILES_RECORDER} ${SOURCE_INCLUDE_FILES_RECORDER} 
                ${API_ABS_INCLUDE_FILES})
    
    target_link_libraries(${RECORDER_TARGET_VN} 
                 ${SimTKCOMMON_RECORDER_LIBRARY_VN}
                 ${MATH_LIBS_TO_USE_VN}
				 ${RECORDER_LIBRARIES})
    
    set_target_properties(${RECORDER_TARGET_VN} PROPERTIES
        PROJECT_LABEL "Code - ${RECORDER_TARGET_VN}" 
        SOVERSION ${SIMBODY_SONAME_VERSION})
    
    # install library; on Windows both .lib and .dll go in the lib directory.
    install(TARGETS ${RECORDER_TARGET_VN} EXPORT SimbodyTargets
                     PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE 
                     GROUP_READ GROUP_WRITE GROUP_EXECUTE 
                     WORLD_READ WORLD_EXECUTE
                     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                     ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
                     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

endif(BUILD_VERSIONED_LIBRARIES)

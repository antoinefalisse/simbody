# This directory defines a target for building SimTKcommon (as a shared library)
# with SimTK::Real defined as Recorder's recorder.
# We could have put this target definition in the parent directory but we want to
# define the following SHARED_LIBRARY definition, and using add_definitions() to
# do this requires that we create this recorderTarget directory.

add_definitions(-DSimTK_SimTKCOMMON_BUILDING_SHARED_LIBRARY)

if(BUILD_UNVERSIONED_LIBRARIES)
    add_library(${RECORDER_TARGET} SHARED
            ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES}
            ${API_ABS_INCLUDE_FILES})

    add_dependencies(${RECORDER_TARGET} PlatformFiles)

    target_link_libraries(${RECORDER_TARGET} ${MATH_LIBS_TO_USE}
            ${RECORDER_LIBRARIES} ${RECORDER_LIBRARIES})

    set_target_properties(${RECORDER_TARGET} PROPERTIES
            PROJECT_LABEL "Code - ${RECORDER_TARGET}"
            SOVERSION ${SIMBODY_SONAME_VERSION})

    # "PUBLIC" means this definition will be passed to the compiler when
    # compiling any target that links to RECORDER_TARGET.
    target_compile_definitions(${RECORDER_TARGET} PUBLIC SimTK_REAL_IS_ADOUBLE)
    
    # Warning 4996 is not propagated to downstream targets  
    if(WIN32 AND NOT ${CMAKE_VERSION} VERSION_LESS 3.0.2) 
        target_compile_options(${RECORDER_TARGET} INTERFACE 
            "$<$<CONFIG:Debug>:/wd4996>")
    endif()

    # install library; on Windows .dll's go in bin directory

    install(TARGETS ${RECORDER_TARGET} EXPORT SimbodyTargets
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
            GROUP_READ GROUP_WRITE GROUP_EXECUTE
            WORLD_READ WORLD_EXECUTE
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} )
endif(BUILD_UNVERSIONED_LIBRARIES)

if(BUILD_VERSIONED_LIBRARIES)
    add_library(${RECORDER_TARGET_VN} SHARED
            ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES}
            ${API_ABS_INCLUDE_FILES})

    add_dependencies(${RECORDER_TARGET_VN} PlatformFiles)

    target_link_libraries(${RECORDER_TARGET_VN} ${MATH_LIBS_TO_USE_VN}
            ${RECORDER_LIBRARIES}  ${RECORDER_LIBRARIES})

    # "PUBLIC" means this definition will be passed to the compiler when
    # compiling any target that links to RECORDER_TARGET.
    target_compile_definitions(${RECORDER_TARGET_VN} PUBLIC SimTK_REAL_IS_ADOUBLE)

    set_target_properties(${RECORDER_TARGET_VN} PROPERTIES
            PROJECT_LABEL "Code - ${RECORDER_TARGET_VN}"
            SOVERSION ${SIMBODY_SONAME_VERSION})
            
    # Warning 4996 is not propagated to downstream targets
    if(WIN32 AND NOT ${CMAKE_VERSION} VERSION_LESS 3.0.2) 
        target_compile_options(${RECORDER_TARGET_VN} INTERFACE 
            "$<$<CONFIG:Debug>:/wd4996>")
    endif()

    # install library; on Windows .dll's go in bin directory

    install(TARGETS ${RECORDER_TARGET_VN} EXPORT SimbodyTargets
            PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
            GROUP_READ GROUP_WRITE GROUP_EXECUTE
            WORLD_READ WORLD_EXECUTE
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} )
endif(BUILD_VERSIONED_LIBRARIES)


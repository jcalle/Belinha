#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "pz" for configuration "Debug"
set_property(TARGET pz APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(pz PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/pzlib/lib/libpz.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS pz )
list(APPEND _IMPORT_CHECK_FILES_FOR_pz "${_IMPORT_PREFIX}/pzlib/lib/libpz.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)

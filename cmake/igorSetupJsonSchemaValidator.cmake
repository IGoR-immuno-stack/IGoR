### igorSetupJsonSchemaValidator.cmake ---

# json-schema-validator, which igor::Model links privately (ModelSchema.cpp).
#
# conda-forge's win-64 build (2.4.0 and 2.4.1) is a Debug DLL: its recipe runs CMake with NMake
# and no build type, so the DLL links the debug runtime (MSVCP140D, ucrtbased), whose STL layout
# is not a Release build's, and every call that hands it a json or a string crashes. On Windows
# the same release is therefore built from source, as a static library inside igor::Model, until
# the feedstock builds Release. IGOR_BUILD_JSON_SCHEMA_VALIDATOR=ON does so on any platform.

option(IGOR_BUILD_JSON_SCHEMA_VALIDATOR
       "Build json-schema-validator from source rather than use the installed package" ${WIN32})

if(IGOR_BUILD_JSON_SCHEMA_VALIDATOR)
  include(FetchContent)
  # The archive conda-forge's recipe builds, at its hash.
  FetchContent_Declare(nlohmann_json_schema_validator
    URL      https://github.com/pboettch/json-schema-validator/archive/2.4.1.tar.gz
    URL_HASH SHA256=0c8bf7ba59bd50e312d82da4bab1d49cf3e2211cfbd659472f9c90505d2a4e4c
    EXCLUDE_FROM_ALL)
  # A subproject neither installs nor builds its tests and examples; it uses nlohmann_json as found.
  set(JSON_VALIDATOR_SHARED_LIBS OFF)
  FetchContent_MakeAvailable(nlohmann_json_schema_validator)
  set_target_properties(nlohmann_json_schema_validator PROPERTIES POSITION_INDEPENDENT_CODE ON)
else()
  find_package(nlohmann_json_schema_validator REQUIRED)
endif()

######################################################################
### igorSetupJsonSchemaValidator.cmake ends here

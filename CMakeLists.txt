cmake_minimum_required(VERSION 3.14)

include(ExternalProject)
include(FetchContent)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(PANSV VERSION 1.0.0)
add_executable(pansv graph.c sketch.cpp main_sketch.cpp main_kan.cpp main_search.cpp main_map.cpp main.cpp)
set_target_properties(pansv PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")
target_compile_options(pansv PRIVATE) # -Wall -Wextra -Wpedantic -Werror)

# gfatools
############
# message(STATUS "gfatools will be built from source")
# ExternalProject_Add(gfatools
#     GIT_REPOSITORY https://github.com/lh3/gfatools.git
#     GIT_TAG c31be8a62efc6bdea4576029f7fbe84f345a6eed
#     BUILD_IN_SOURCE 1
#     UPDATE_COMMAND ""
#     CONFIGURE_COMMAND ""
#     BUILD_COMMAND "make"
#     INSTALL_COMMAND ""
# )
# ExternalProject_Get_Property(gfatools SOURCE_DIR)
# SET(GFA_INCLUDE_DIR ${SOURCE_DIR})

find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)

# ksw2
########
message(STATUS "gfatools will be built from source")
ExternalProject_Add(ksw2
    GIT_REPOSITORY https://github.com/lh3/ksw2.git
    GIT_TAG ""
    BUILD_IN_SOURCE 1
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make
    INSTALL_COMMAND ""
)
ExternalProject_Get_Property(ksw2 SOURCE_DIR)
SET(KSW2_INCLUDE_DIR ${SOURCE_DIR})

# # sdsl-lite (vg fork)
# #######################
# message(STATUS "sdsl-lite (vg fork) will be built from source")
# ExternalProject_Add(sdsl
#     GIT_REPOSITORY https://github.com/vgteam/sdsl-lite.git
#     GIT_TAG 03288600910bd9cec2b70cd52db12066361b2ca9
#     BUILD_IN_SOURCE 1
#     UPDATE_COMMAND ""
#     CONFIGURE_COMMAND ""
#     BUILD_COMMAND ""
#     INSTALL_COMMAND bash install.sh .
# )
# ExternalProject_Get_Property(sdsl SOURCE_DIR)
# SET(SDSL_SOURCE_DIR ${SOURCE_DIR})
# SET(SDSL_INCLUDE_DIR ${SDSL_SOURCE_DIR}/include)
# add_library(SDSL STATIC IMPORTED)
# set_target_properties(SDSL PROPERTIES IMPORTED_LOCATION ${SDSL_SOURCE_DIR}/lib/libsdsl.a)
# add_dependencies(SDSL sdsl)
#
# # gbwt
# ########
# message(STATUS "gbwt will be built from source")
# ExternalProject_Add(gbwt
#     GIT_REPOSITORY https://github.com/jltsiren/gbwt.git
#     GIT_TAG f198496a890410b6f37b55428ca8bd8f73d4bf52
#     BUILD_IN_SOURCE 1
#     UPDATE_COMMAND ""
#     CONFIGURE_COMMAND ""
#     BUILD_COMMAND make SDSL_DIR=${SDSL_SOURCE_DIR}
#     INSTALL_COMMAND ""
# )
# ExternalProject_Get_Property(gbwt SOURCE_DIR)
# SET(GBWT_SOURCE_DIR ${SOURCE_DIR})
# SET(GBWT_INCLUDE_DIR ${GBWT_SOURCE_DIR}/include)
# add_library(GBWT STATIC IMPORTED)
# set_target_properties(GBWT PROPERTIES IMPORTED_LOCATION ${GBWT_SOURCE_DIR}/lib/libgbwt.a)
# add_dependencies(GBWT gbwt)
# add_dependencies(GBWT SDSL)

# ropebwt3
############
message(STATUS "ropebwt3 will be built from source")
ExternalProject_Add(rb3
    GIT_REPOSITORY https://github.com/lh3/ropebwt3.git
    GIT_TAG 999c17f459d224e623ff4d517e027e50f82fb67b
    BUILD_IN_SOURCE 1
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND "make"
    INSTALL_COMMAND ""
)
ExternalProject_Get_Property(rb3 SOURCE_DIR)
SET(RB3_SOURCE_DIR ${SOURCE_DIR})
SET(RB3_INCLUDE_DIR ${RB3_SOURCE_DIR})
set(rb3-src ${RB3_SOURCE_DIR}/fm-index.c ${RB3_SOURCE_DIR}/rld0.c ${RB3_SOURCE_DIR}/mrope.c ${RB3_SOURCE_DIR}/rope.c ${RB3_SOURCE_DIR}/io.c ${RB3_SOURCE_DIR}/rle.c ${RB3_SOURCE_DIR}/kthread.c ${RB3_SOURCE_DIR}/kalloc.c ${RB3_SOURCE_DIR}/misc.c ${RB3_SOURCE_DIR}/ssa.c ${RB3_SOURCE_DIR}/dawg.c ${RB3_SOURCE_DIR}/libsais16.c)
set_source_files_properties(${rb3-src} PROPERTIES GENERATED TRUE)
add_library(RB3 OBJECT ${rb3-src})
add_dependencies(RB3 rb3)

# abPOA
########
SET(PATCH_CMD "")
ExternalProject_Add(abpoa
  GIT_REPOSITORY https://github.com/yangao07/abPOA.git
  GIT_TAG 138407dcf01d0a147a5b98931a96c6cfe6694775 # v1.5.3
  PATCH_COMMAND ${PATCH_CMD}
  INSTALL_COMMAND ""
)
ExternalProject_Get_Property(abpoa SOURCE_DIR BINARY_DIR)
SET(ABPOA_SOURCE_DIR ${SOURCE_DIR})
SET(ABPOA_INCLUDE_DIR ${ABPOA_SOURCE_DIR}/include)
# NOTE: abpoa library will be set manually
# add_library(ABPOA STATIC IMPORTED)
# set_target_properties(ABPOA PROPERTIES IMPORTED_LOCATION ${ABPOA_SOURCE_DIR}/lib/libabpoa.a)

############################################
### final setup (includes and libraries) ###
############################################

add_dependencies(pansv RB3 abpoa ksw2) # gbwt)

target_include_directories(pansv
    PRIVATE ${RB3_INCLUDE_DIR}
    # PRIVATE ${SDSL_INCLUDE_DIR}
    # PRIVATE ${GBWT_INCLUDE_DIR}
    PRIVATE ${KSW2_INCLUDE_DIR}
    PRIVATE ${ABPOA_INCLUDE_DIR}
    PRIVATE ${PROJECT_SOURCE_DIR}
)
# target_link_libraries(GBWT
#   PUBLIC SDSL)

target_link_libraries(pansv
    PUBLIC ${KSW2_INCLUDE_DIR}/ksw2_extz.o
    PUBLIC ${KSW2_INCLUDE_DIR}/ksw2_extz2_sse.o
    PUBLIC ${BINARY_DIR}/lib/libabpoa.a
    # PUBLIC GBWT
    # PUBLIC SDSL
    PUBLIC RB3
    PUBLIC z
    PUBLIC pthread
    PUBLIC OpenMP::OpenMP_CXX
)

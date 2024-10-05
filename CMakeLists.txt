cmake_minimum_required(VERSION 3.14)

include(ExternalProject)
include(FetchContent)

#set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(PANSV VERSION 1.0.0)
add_executable(pansv gsketch.cpp main_call.cpp main.cpp)
# add_executable(pansv kalloc.c gfa-aug.c gfa-base.c  gfa-bbl.c  gfa-ed.c gfa-io.c index.c kthread.c options.c sys.c misc.c sketch.c main.cpp)
set_target_properties(pansv PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")

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


# # KMC3
# #######
# message(STATUS "KMC3 will be built from source")
# ExternalProject_Add(kmc
#     GIT_REPOSITORY https://github.com/refresh-bio/KMC.git
#     GIT_TAG d6c24dc88010508336d824bbffbf1ef3f49fddf8
#     UPDATE_COMMAND ""
#     BUILD_IN_SOURCE 1
#     PATCH_COMMAND ""
#     CONFIGURE_COMMAND ""
#     BUILD_COMMAND make kmc
#     INSTALL_COMMAND ""
# )
# ExternalProject_Get_Property(kmc SOURCE_DIR)
# SET(KMC_SOURCE_DIR ${SOURCE_DIR})
# SET(KMC_INCLUDE_DIR ${KMC_SOURCE_DIR})
# add_library(KMC STATIC IMPORTED)
# set_target_properties(KMC PROPERTIES IMPORTED_LOCATION ${KMC_INCLUDE_DIR}/bin/libkmc_core.a)
# add_dependencies(KMC kmc)

# abPOA
########
SET(PATCH_CMD "")
# if(DEFINED CONDAPREFIX OR HOLYBUILD)
#   # we have to remove the -native flag
#   message(STATUS "abpoa will be build from source (no -native)")
#   SET(PATCH_CMD sed -i "/native/d" CMakeLists.txt)
# else ()
#   message(STATUS "abpoa will be build from source (with -native)")
# endif ()
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

add_dependencies(pansv RB3 abpoa) # KMC

target_include_directories(pansv
    PRIVATE ${RB3_INCLUDE_DIR}
    # PRIVATE ${KMC_INCLUDE_DIR}
    PRIVATE ${ABPOA_INCLUDE_DIR}
    PRIVATE ${PROJECT_SOURCE_DIR}
)

target_link_libraries(pansv
    # PUBLIC ${GFA_INCLUDE_DIR}/kalloc.o
    # PUBLIC ${GFA_INCLUDE_DIR}/gfa-base.o
    # PUBLIC ${GFA_INCLUDE_DIR}/gfa-io.o
    PUBLIC ${BINARY_DIR}/lib/libabpoa.a
    PUBLIC RB3
    # PUBLIC KMC
    PUBLIC z
    PUBLIC pthread
    PUBLIC OpenMP::OpenMP_CXX
)

target_link_libraries(
    pansv
)

cmake_minimum_required(VERSION 3.14)

include(ExternalProject)
include(FetchContent)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(PALSS VERSION 0.0.1)
add_executable(palss misc.c kmer.c sketch.c graph.c main_sketch.c main.cpp)
# add_executable(palss sfs.cpp graph.c sketch.cpp main_sketch.cpp main_dump.cpp main_search.cpp main_kan.cpp main_call.cpp main.cpp) # main_map.cpp
set_target_properties(palss PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")
target_compile_options(palss PRIVATE) # -Wno-unused-function -Wno-misleading-indentation) # -Wall -Wextra -Wpedantic -Werror)

add_definitions("-Wp,-w")

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
  GIT_TAG e6bb6fdfa40d573558e5e2b545bdf5769631eaaf # v1.5.3 with bug fixes
  BUILD_IN_SOURCE 1
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND make
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

add_dependencies(palss RB3 abpoa ksw2) # gbwt)

target_include_directories(palss
    PRIVATE ${RB3_INCLUDE_DIR}
    PRIVATE ${KSW2_INCLUDE_DIR}
    PRIVATE ${ABPOA_INCLUDE_DIR}
    PRIVATE ${PROJECT_SOURCE_DIR}
)

target_link_libraries(palss
    PUBLIC ${KSW2_INCLUDE_DIR}/ksw2_extz2_sse.o
    PUBLIC ${KSW2_INCLUDE_DIR}/ksw2_extd2_sse.o
    PUBLIC ${BINARY_DIR}/lib/libabpoa.a
    PUBLIC RB3
    PUBLIC z
    PUBLIC pthread
    PUBLIC OpenMP::OpenMP_CXX
)

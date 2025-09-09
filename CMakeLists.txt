cmake_minimum_required(VERSION 3.14)

include(ExternalProject)
include(FetchContent)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(PALSS VERSION 0.0.1)

add_executable(palss misc.cpp kmer.cpp graph.cpp sketch.cpp main_sketch.cpp main_kan.cpp main_stats.cpp main_augment.cpp main.cpp)
# palss misc.cpp kmer.cpp segments.cpp graph.cpp  main_sketch.cpp main_search.cpp main.cpp

set_target_properties(palss PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")
target_compile_options(palss PRIVATE -Wall -Wno-char-subscripts) # -Wno-unused-function -Wno-misleading-indentation) # -Wall -Wextra -Wpedantic -Werror)

# add_definitions("-Wp,-w")

find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)

# ksw2
########
message(STATUS "ksw2 will be built from source")
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
    GIT_TAG 0ea3919ed21f10d857d508b0ea728a27abbd8a35 # 3.9 (r259)
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
  GIT_TAG f01e0da57ca9398c5a7900026e18d7222b55ee9f # v1.5.4
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
    PUBLIC OpenMP::OpenMP_C
    PUBLIC OpenMP::OpenMP_CXX
)

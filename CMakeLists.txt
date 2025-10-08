cmake_minimum_required(VERSION 3.14)

include(ExternalProject)
include(FetchContent)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(PALSS VERSION 0.0.1)

add_executable(palss misc.cpp kmer.cpp sketch.cpp main_sketch.cpp main_stats.cpp main.cpp)
#  main_index.cpp main_search.cpp main_kan.cpp main_map.cpp main_call.cpp graph.cpp sfs.c misc.cpp sketch.c kmer.c
# misc.c kmer.c sfs.c sketch.c segment.c path.c graph.c main_sketch.c main_dump.c main_kan.c main_search.c main_test.c main_call.cpp main.c)
set_target_properties(palss PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")
target_compile_options(palss PRIVATE -Wall) # -Wno-char-subscripts -Wno-unused-function -Wno-misleading-indentation) # -Wall -Wextra -Wpedantic -Werror)
# add_definitions("-Wp,-w")
find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)

SET(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)
SET(LIB_DIR ${CMAKE_SOURCE_DIR}/lib)

# libhandlegraph
message(STATUS "libhandlegraph will be built from source")
ExternalProject_Add(lhg
    GIT_REPOSITORY https://github.com/vgteam/libhandlegraph.git
    GIT_TAG 0e70dadb5054568d8071e280b3b7b11b5658937f
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR} -DCMAKE_INSTALL_LIBDIR=${LIB_DIR}
)

# sdsl-lite
message(STATUS "sdsl-lite will be built from source")
ExternalProject_Add(sdsl
    GIT_REPOSITORY https://github.com/vgteam/sdsl-lite.git
    GIT_TAG 8abd7c9c183055d74f9b57915d5e36967a453062
    BUILD_IN_SOURCE 1
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND BUILD_PORTABLE=1 ./install.sh ${CMAKE_SOURCE_DIR}
)
ExternalProject_Get_Property(sdsl SOURCE_DIR)
SET(SDSL_SRC_DIR ${SOURCE_DIR})

# gbwt
message(STATUS "gbwt will be built from source")
ExternalProject_Add(gbwt
    GIT_REPOSITORY https://github.com/jltsiren/gbwt.git
    GIT_TAG 32187dc2c866a17ef474529895965303034ed9fd
    BUILD_IN_SOURCE 1
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make SDSL_DIR=${SDSL_SRC_DIR}
    INSTALL_COMMAND cp -r include/gbwt ${INCLUDE_DIR}
    COMMAND cp lib/libgbwt.a ${LIB_DIR}
)

# gbwtgraph
message(STATUS "gbwtgraph will be built from source")
ExternalProject_Add(gbwtgraph
    GIT_REPOSITORY https://github.com/jltsiren/gbwtgraph.git
    GIT_TAG 197e1c3f454cebb919a9f5564021ce46273a162d
    BUILD_IN_SOURCE 1
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make SDSL_DIR=${SDSL_SRC_DIR}
    INSTALL_COMMAND cp -r include/gbwtgraph ${INCLUDE_DIR}
    COMMAND cp lib/libgbwtgraph.a ${LIB_DIR}
)

# # ksw2
# message(STATUS "ksw2 will be built from source")
# ExternalProject_Add(ksw2
#     GIT_REPOSITORY https://github.com/lh3/ksw2.git
#     GIT_TAG ""
#     BUILD_IN_SOURCE 1
#     UPDATE_COMMAND ""
#     CONFIGURE_COMMAND ""
#     BUILD_COMMAND make
#     INSTALL_COMMAND ""
# )
# ExternalProject_Get_Property(ksw2 SOURCE_DIR)
# SET(KSW2_INCLUDE_DIR ${SOURCE_DIR})

# ropebwt3
message(STATUS "ropebwt3 will be built from source")
ExternalProject_Add(rb3
    GIT_REPOSITORY https://github.com/lh3/ropebwt3.git
    GIT_TAG dd8d28575cb889650ed45a3a8ffed3688a855293
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

# # abPOA
# SET(PATCH_CMD "")
# ExternalProject_Add(abpoa
#   GIT_REPOSITORY https://github.com/yangao07/abPOA.git
#   GIT_TAG f01e0da57ca9398c5a7900026e18d7222b55ee9f # v1.5.4
#   BUILD_IN_SOURCE 1
#   UPDATE_COMMAND ""
#   CONFIGURE_COMMAND ""
#   BUILD_COMMAND make
#   INSTALL_COMMAND ""
# )
# ExternalProject_Get_Property(abpoa SOURCE_DIR BINARY_DIR)
# SET(ABPOA_SOURCE_DIR ${SOURCE_DIR})
# SET(ABPOA_INCLUDE_DIR ${ABPOA_SOURCE_DIR}/include)
# # NOTE: abpoa library will be set manually
# # add_library(ABPOA STATIC IMPORTED)
# # set_target_properties(ABPOA PROPERTIES IMPORTED_LOCATION ${ABPOA_SOURCE_DIR}/lib/libabpoa.a)


############################################
### final setup (includes and libraries) ###
############################################

add_dependencies(gbwt sdsl)
add_dependencies(gbwtgraph gbwt lhg)
add_dependencies(palss RB3 gbwtgraph) # abpoa ksw2)

target_include_directories(palss
  PUBLIC ${CMAKE_SOURCE_DIR}/include
  PUBLIC ${RB3_INCLUDE_DIR}
  # PUBLIC ${KSW2_INCLUDE_DIR}
  # PUBLIC ${ABPOA_INCLUDE_DIR}
  PUBLIC ${PROJECT_SOURCE_DIR}
)

target_link_libraries(palss
    # PUBLIC ${KSW2_INCLUDE_DIR}/ksw2_extz2_sse.o
    # PUBLIC ${KSW2_INCLUDE_DIR}/ksw2_extd2_sse.o
    # PUBLIC ${ABPOA_SOURCE_DIR}/lib/libabpoa.a

    PUBLIC ${CMAKE_SOURCE_DIR}/lib/libgbwtgraph.a
    PUBLIC ${CMAKE_SOURCE_DIR}/lib/libgbwt.a
    PUBLIC ${CMAKE_SOURCE_DIR}/lib/libsdsl.a
    PUBLIC ${CMAKE_SOURCE_DIR}/lib/libhandlegraph.a

    PUBLIC RB3
    PUBLIC z
    PUBLIC pthread
    PUBLIC OpenMP::OpenMP_C
    PUBLIC OpenMP::OpenMP_CXX
)
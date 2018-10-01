project(faldoi_qtcreator)
cmake_minimum_required(VERSION 2.8.1)
aux_source_directory(. SRC_LIST)

message("Build type: ${CMAKE_BUILD_TYPE}")

include_directories(./lib)

find_library(ZLIB libz.so PATHS lib)

find_library(PNG libpng.so PATHS lib)

#find_library(JPEG libjpeg.so PATHS lib)
#include_directories(${JPEG_INCLUDE_DIR})

#find_library(TIFF libtiff.so PATHS lib)
#include_directories(${TIFF_INCLUDE_DIR})

find_package(OpenMP)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -std=c99 -march=native -mtune=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -std=c++11")

set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)

if (CMAKE_CXX_COMPILER EQUAL clang++)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fvectorize")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvectorize")
elseif ()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ftree-vectorize -ftree-loop-vecorize")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftree-vectorize -ftree-loop-vecorize")
endif ()

set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS} -O3")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3")

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS} -O0 -ggdb -DNDEBUG -Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -O0 -ggdb -DNDEBUG -Wall -Wextra")

#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=undefined -fsanitize=thread")

#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address")

SET(SHARED_C_SRC iio.c mask.c xmalloc.c bicubic_interpolation.c elap_recsep.c)
SET(SHARED_CPP_SRC tvl2_model.cpp nltv_model.cpp tvcsad_model.cpp nltvcsad_model.cpp tvl2w_model.cpp nltvcsadw_model.cpp nltvw_model.cpp tvcsadw_model.cpp aux_energy_model.cpp energy_model.cpp tvl2_model_occ.cpp utils.cpp utils_preprocess.cpp)


add_executable(sparse_flow ${SHARED_C_SRC} ${SHARED_CPP_SRC} sparse_flow.cpp)
add_executable(local_faldoi ${SHARED_C_SRC} ${SHARED_CPP_SRC} local_faldoi.cpp)
add_executable(global_faldoi ${SHARED_C_SRC} ${SHARED_CPP_SRC} global_faldoi.cpp)
target_link_libraries(sparse_flow ${PNG} ${ZLIB}) # ${JPEG} ${TIFF} ${ZLIB})
target_link_libraries(local_faldoi ${PNG} ${ZLIB}) # ${JPEG} ${TIFF} ${ZLIB})
target_link_libraries(global_faldoi ${PNG} ${ZLIB}) # ${JPEG} ${TIFF} ${ZLIB})



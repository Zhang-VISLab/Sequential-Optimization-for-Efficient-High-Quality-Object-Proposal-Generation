#
set(INC_DIRS ${BING++_SOURCE_DIR}/LibLinear ${BING++_SOURCE_DIR}/LibLinear/blas ${BING++_SOURCE_DIR}/../ext ${OPENCV_PATH}include ${BING++_SOURCE_DIR}/BINGpp/FelzenSegment )
set(LIB_DIRS ${BING++_BINARY_DIR}/LibLinear ${BING++_SOURCE_DIR}/../ext ${OPENCV_PATH}lib)

add_subdirectory(ORUtils)

#OPENCV_DIR should be set as an environment variable in the shell if it is not
#present in the standard directories

include_directories(${INC_DIRS})
link_directories(${LIB_DIRS})

if (${CMAKE_COMPILER_IS_GNUCXX})
  set(OPENMP_FLAG "-fopenmp -std=c++11 -ftree-vectorize")
endif (${CMAKE_COMPILER_IS_GNUCXX})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPENMP_FLAG} -O3") # -O0 -pg

set(SRC_FILES CmFile.cpp CmShow.cpp DataSetVOC.cpp FilterBING.cpp SegmentRecursive.cpp BINGPP.cpp stdafx.cpp)
add_library(BINGPP ${SRC_FILES})
add_executable(${PROJECT_NAME} ${SRC_FILES} main.cpp)
target_link_libraries(${PROJECT_NAME} opencv_highgui opencv_core opencv_imgproc opencv_ts opencv_imgcodecs liblinear blas vl)

PROJECT(HighOrderCCD)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#COMPILE SETTINGS
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
ENDIF()
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/")
SET(CUDA_SEPARABLE_COMPILATION ON)
IF(CMAKE_BUILD_TYPE STREQUAL "Debug")
  MESSAGE(STATUS "Build type: Debug")
  IF(MSVC)
    SET(CMAKE_CXX_FLAGS "/openmp /bigobj -Wall -std=c++11 -frounding-math ${CMAKE_CXX_FLAGS}")
  ELSE(MSVC)
    SET(CMAKE_CXX_FLAGS "-fopenmp -O0 -Wall -fPIC -std=c++11 -Wno-unused-function -frounding-math ${CMAKE_CXX_FLAGS}")
  ENDIF(MSVC)
ELSE()
  MESSAGE(STATUS "Build type: Release")
  IF(MSVC)
    SET(CMAKE_CXX_FLAGS "/openmp /bigobj -Wall -std=c++11 -frounding-math ${CMAKE_CXX_FLAGS}")
  ELSE(MSVC)
    SET(CMAKE_CXX_FLAGS "-fopenmp -O3 -Wall -fPIC -std=c++11 -Wno-unused-function -frounding-math ${CMAKE_CXX_FLAGS}")
  ENDIF(MSVC)
ENDIF()
set (CMAKE_CXX_STANDARD 11)
SET(CMAKE_CONFIGURATION_TYPES ${CMAKE_BUILD_TYPE})

#LIBRARY
add_subdirectory(lib/opengjk)
#add_subdirectory(Exact-CDD)
INCLUDE("${PROJECT_SOURCE_DIR}/cmake/LoadLibs.cmake")
#THIRD PARTY HEADERS
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR})

#link_directories( ${PROJECT_SOURCE_DIR}/mosek/lib/mosek8_1 )
#include_directories(${PROJECT_SOURCE_DIR}/mosek/include)
add_definitions(/D_USE_MATH_DEFINES)
add_definitions(/D_CRT_SECURE_NO_WARNINGS)
add_definitions(/D_CRT_NONSTDC_NO_DEPRECATE)

SET(OPENMESH_PATH "D:/programfiles/OpenMesh_9.0/include")
SET(OPENMESH_LIB_PATH "D:/programfiles/OpenMesh_9.0/lib")
include_directories(${OPENMESH_PATH})
link_directories(${OPENMESH_LIB_PATH})

#CGAL
SET(BOOST_PATH "D:/programfiles/boost_1_80_0")
SET(BOOST_LIB_PATH "D:/programfiles/boost_1_80_0/stage/lib")
SET(CGAL_PATH "D:/programfiles/CGAL-5.5.1/include")
SET(GMP_PATH "D:/programfiles/CGAL-5.5.1/auxiliary/gmp/include")
SET(GMP_LIB_PATH "D:/programfiles/CGAL-5.5.1/auxiliary/gmp/lib")
include_directories(${BOOST_PATH})
link_directories(${BOOST_LIB_PATH})
include_directories(${CGAL_PATH})
include_directories(${GMP_PATH})
link_directories(${GMP_LIB_PATH})

#OCCT
SET(OPENCASCADE_PATH "D:/programfiles/OpenCASCADE-7.7.0-vc14-64/opencascade-7.7.0/inc")
SET(OPENCASCADE_LIB_PATH "D:/programfiles/OpenCASCADE-7.7.0-vc14-64/opencascade-7.7.0/win64/vc14/lib")
include_directories(${OPENCASCADE_PATH})
link_directories(${OPENCASCADE_LIB_PATH})

MACRO(SEARCH_SOURCE NAME DIRPRJ)
  FILE(GLOB_RECURSE header${NAME} ${PROJECT_SOURCE_DIR}/${DIRPRJ}/*.h ${PROJECT_SOURCE_DIR}/${DIRPRJ}/*.hpp ${PROJECT_SOURCE_DIR}/${DIRPRJ}/*.hh)
  FILE(GLOB_RECURSE source${NAME} ${PROJECT_SOURCE_DIR}/${DIRPRJ}/*.cpp ${PROJECT_SOURCE_DIR}/${DIRPRJ}/*.cpp ${PROJECT_SOURCE_DIR}/${DIRPRJ}/*.cc ${PROJECT_SOURCE_DIR}/${DIRPRJ}/*.c)
  FILE(GLOB_RECURSE cuda_source${NAME} ${PROJECT_SOURCE_DIR}/${DIRPRJ}/*.cu)
ENDMACRO(SEARCH_SOURCE)
SEARCH_SOURCE(HighOrderCCD HighOrderCCD)

#BUILD LIBRARY
ADD_LIBRARY(HighOrderCCD SHARED
  ${sourceHighOrderCCD} ${headerHighOrderCCD})
#target_link_libraries(HighOrderCCD m )
target_link_libraries(HighOrderCCD openGJKlib )
#target_link_libraries(HighOrderCCD RationalCCD )
target_link_libraries(HighOrderCCD igl::core igl::opengl igl::opengl_glfw)
#target_link_libraries(HighOrderCCD mosek64)
target_link_libraries(HighOrderCCD ${GMP_LIBRARIES})
target_link_libraries(${PROJECT_NAME}
    libgmp-10
    libboost_filesystem-vc143-mt-gd-x64-1_80
    optimized ${OPENMESH_LIB_PATH}/OpenMeshCore.lib
    debug ${OPENMESH_LIB_PATH}/OpenMeshCored.lib
    optimized ${OPENMESH_LIB_PATH}/OpenMeshTools.lib
    debug ${OPENMESH_LIB_PATH}/OpenMeshToolsd.lib
    #OpenCASCADE
    TKernel.lib
    TKMath.lib
    TKG2d.lib
    TKG3d.lib
    TKBRep.lib
    TKGeomBase.lib
    TKGeomAlgo.lib
    TKTopAlgo.lib
    TKPrim.lib
    TKBO.lib
    TKBool.lib
    TKOffset.lib
    TKService.lib
    TKV3d.lib
    TKOpenGl.lib
    TKFillet.lib
    TKXSBase.lib
    TKXDEIGES.lib
    TKSTEP.lib
    TKIGES.lib
    TKXMESH.lib
    TKMESH.lib
    TKMESHVS.lib
)

TARGET_LINK_LIBRARIES(HighOrderCCD ${ALL_LIBRARIES} ${ALL_STATIC_LIBRARIES})

#EXE
MACRO(ADD_EXE NAME)
  ADD_EXECUTABLE(${NAME} ${NAME}.cpp)
  #target_link_libraries(${NAME} m )
  target_link_libraries(${NAME} openGJKlib )
  #target_link_libraries(${NAME} RationalCCD )
  target_link_libraries(${NAME} igl::core igl::opengl igl::opengl_glfw)
  #target_link_libraries(${NAME} mosek64) 
  target_link_libraries(${NAME} ${GMP_LIBRARIES})
  IF(MSVC)
    TARGET_LINK_LIBRARIES(${NAME} HighOrderCCD ${ALL_STATIC_LIBRARIES})
  ELSE(MSVC)
    TARGET_LINK_LIBRARIES(${NAME} HighOrderCCD gomp ${ALL_STATIC_LIBRARIES})
  ENDIF(MSVC)
ENDMACRO(ADD_EXE)
#GEOMETRY
ADD_EXE(debugPathPlanning3D)



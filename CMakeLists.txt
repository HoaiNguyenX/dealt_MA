# Set the name of the project and target:
set_property(GLOBAL PROPERTY CXX_STANDARD 20)
SET(TARGET "tsplines")
SET(CMAKE_C_COMPILER, "mpicc")
SET(CMAKE_CXX_COMPILER, "mpicxx")

# Declare all source files the target consists of. Here, this is only
# the one step-X.cc file, but as you expand your project you may wish
# to add other source files as well. If your project becomes much larger,
# you may want to either replace the following statement by something like
#  FILE(GLOB_RECURSE TARGET_SRC  "source/*.cc")
#  FILE(GLOB_RECURSE TARGET_INC  "include/*.h")
#  SET(TARGET_SRC ${TARGET_SRC}  ${TARGET_INC})
# or switch altogether to the large project CMakeLists.txt file discussed
# in the "CMake in user projects" page accessible from the "User info"
# page of the documentation.

FILE(GLOB_RECURSE TARGET_SRC  "src/*.cc")
INCLUDE_DIRECTORIES(inc)
SET(TARGET_SRC 
  ${TARGET_SRC}  
  main.cc
  )

# Usually, you will not need to modify anything beyond this point...
CMAKE_MINIMUM_REQUIRED(VERSION 3.1.0)

set (DEAL_II_DIR /home/ifam-studenten/nguyen/Software/dealii-9.5.1/dealii-9.5.1)

FIND_PACKAGE(deal.II 9.5.1 QUIET
#FIND_PACKAGE(deal.II EXACT 9.2.0 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )

IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()


DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})
DEAL_II_INVOKE_AUTOPILOT()


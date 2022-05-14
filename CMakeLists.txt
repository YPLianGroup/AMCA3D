# For optimized version, run with option -DCMAKE_BUILD_TYPE=Release
cmake_minimum_required(VERSION 2.8)

PROJECT(AMCA3D)

find_package(Boost COMPONENTS system filesystem program_options REQUIRED)
find_package(MPI REQUIRED)

set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Source and header files
set(SOURCES
        ./CA/src/CafeParsing.cpp
        ./CA/src/Cafe.cpp
        ./CA/src/CellularAutomataSolver.cpp
        ./CA/src/Simulation.cpp
        ./CA/src/CafeEnv.cpp
        ./CA/src/CafeSolver.cpp
        ./CA/src/CellularAutomataManager_3D.cpp
        ./CA/src/CellularAutomataManager_2D.cpp
        ./CA/src/CellularAutomataManager.cpp
        ./CA/src/Envelope.cpp
        ./CA/src/Grain.cpp
        ./CA/src/Octahedron.cpp
        ./CA/src/Orientation.cpp
        ./CA/src/ParallelCommManager_2D.cpp
        ./CA/src/ParallelCommManager_3D.cpp
        ./CA/src/ParallelCommManager.cpp
        ./CA/src/ProblemPhysics.cpp
        ./CA/src/Timer.cpp
        ./CA/src/VtrManager.cpp
        ./CA/src/VtuManager.cpp
        ./CA/src/Solver.cpp
        ./CA/src/Simulation.cpp
        ./CA/src/Node.cpp
        ./CA/src/MapVoxelManager.cpp
        ./CA/src/LoadDataFromFile.cpp
        ./CA/src/FiniteElementMethodSolver.cpp
        ./CA/src/FiniteElementManager.cpp
        ./CA/src/Element.cpp
        ./CA/src/VtkManager.cpp
        ./CA/src/CafeMacro.cpp
        ./CA/src/CellularAutomataManager_3DRemelting.cpp
        )


include_directories ( ${CMAKE_SOURCE_DIR}/./CA/include ${MPI_INCLUDE_PATH} ${Boost_INCLUDE_DIR} ${YAML_INCLUDE_DIR})
link_directories( ${Boost_LIBRARY_DIR} ${YAML_LIBRARY_DIR})

add_executable(AMCA3D ${SOURCES})

target_link_libraries( AMCA3D ${Boost_LIBRARIES} yaml-cpp -lpthread)


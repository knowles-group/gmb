cmake_minimum_required(VERSION 3.15)

project(gmb VERSION 1.0 LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)

option(GMB_COMPILE_LIBTENSOR "whether to build libtensor if not found" ON)

if (MPI)
    find_package(MPI REQUIRED C)
    set(GMB_MPI ON)
endif ()

find_package(OpenMP REQUIRED C)
if (OPENMP_FOUND)
    set(GMB_OMP ON)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include(FetchContent)
FetchContent_Declare(
        LibraryManager
        GIT_REPOSITORY https://gitlab+deploy-token-178345:RQZa2rmdwnV7Hdwg8XEq@gitlab.com/molpro/LibraryManager.git
        GIT_TAG 0.7.2
)
FetchContent_MakeAvailable(LibraryManager)

add_subdirectory(dependencies)

add_subdirectory(src)
if (CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    enable_testing()
    add_subdirectory(test)
endif ()
#add_subdirectory(docs)



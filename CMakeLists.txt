cmake_minimum_required(VERSION 3.15)

project(gmb VERSION 1.0 LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)

# option(GMB_COMPILE_LIBTENSOR "whether to build libtensor if not found" ON)
include(FetchContent)
FetchContent_Declare(libtensor GIT_REPOSITORY https://github.com/epifanovsky/libtensor
                               GIT_TAG 759566bfe96a4dea4f352be6e75cf126a66acc29)

set(LIBTENSOR_STANDALONE 1)
add_definitions(-DUSE_PTHREADS)
FetchContent_MakeAvailable(libtensor)

if (MPI)
    find_package(MPI REQUIRED C)
    set(GMB_MPI ON)
endif ()

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



cmake_minimum_required(VERSION 3.15)

project(gmb VERSION 1.0 LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)

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
enable_testing()
add_subdirectory(test)
#add_subdirectory(docs)



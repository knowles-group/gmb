cmake_minimum_required(VERSION 3.15)

project(gmb VERSION 1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)

########## Compiler Flags ##########
if (MPI)
    find_package(MPI REQUIRED C)
    set(GMB_MPI ON)
endif ()

########## Dependencies ##########
# library manager
include(FetchContent)
FetchContent_Declare(
        LibraryManager
        GIT_REPOSITORY https://gitlab+deploy-token-178345:RQZa2rmdwnV7Hdwg8XEq@gitlab.com/molpro/LibraryManager.git
        GIT_TAG 0.7.0
)
FetchContent_MakeAvailable(LibraryManager)

# libtensor
FetchContent_Declare(libtensor
#        GIT_REPOSITORY https://github.com/epifanovsky/libtensor
        GIT_REPOSITORY https://github.com/pjknowles/libtensor
        GIT_TAG gmb)
set(LIBTENSOR_STANDALONE 1)
add_definitions(-DUSE_PTHREADS)
FetchContent_MakeAvailable(libtensor)

add_subdirectory(dependencies)

########## Targets ##########
add_subdirectory(src)
enable_testing()
add_subdirectory(test)
#add_subdirectory(docs)



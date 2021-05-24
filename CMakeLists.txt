cmake_minimum_required(VERSION 3.15)

project(gmb VERSION 1.0 LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)
include(FetchContent)
FetchContent_Declare(
        LibraryManager
        GIT_REPOSITORY https://gitlab+deploy-token-178345:RQZa2rmdwnV7Hdwg8XEq@gitlab.com/molpro/LibraryManager.git
        GIT_TAG 0.7.0
)
FetchContent_MakeAvailable(LibraryManager)

include(FetchContent)
FetchContent_Declare(libtensor GIT_REPOSITORY https://github.com/epifanovsky/libtensor
        GIT_TAG 759566bfe96a4dea4f352be6e75cf126a66acc29)

set(LIBTENSOR_STANDALONE 1)
add_definitions(-DUSE_PTHREADS)
FetchContent_MakeAvailable(libtensor)

add_subdirectory(dependencies)

add_subdirectory(src)
enable_testing()
add_subdirectory(test)


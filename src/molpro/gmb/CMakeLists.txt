LibraryManager_Append(${PROJECT_NAME}
        PUBLIC_HEADER gmb.h
        SOURCES gmb.cpp init.cpp utils.cpp libtensor_utils.cpp get_integral.cpp
        )

DependencyManager_Populate(Profiler)
DependencyManager_Populate(utilities)
set(LINEARALGEBRA_ARRAY_HDF5 OFF)
DependencyManager_Populate(LinearAlgebra)
DependencyManager_Populate(FCIdump)


FetchContent_Declare(libtensor
        #        GIT_REPOSITORY https://github.com/epifanovsky/libtensor
        GIT_REPOSITORY https://github.com/pjknowles/libtensor
        GIT_TAG gmb)
set(LIBTENSOR_STANDALONE 1)
FetchContent_MakeAvailable(libtensor)
target_compile_definitions(tensor PUBLIC USE_PTHREADS)

target_link_libraries(${PROJECT_NAME} PUBLIC tensor molpro::utilities molpro::Profiler molpro::FCIdump molpro::LinearAlgebra)

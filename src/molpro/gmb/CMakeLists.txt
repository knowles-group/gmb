LibraryManager_Append(${PROJECT_NAME}
        PUBLIC_HEADER gmb.h
        SOURCES gmb.cpp init.cpp utils.cpp libtensor_utils.cpp get_integral.cpp
        )

DependencyManager_Populate(Profiler)
DependencyManager_Populate(utilities)
DependencyManager_Populate(LinearAlgebra)
DependencyManager_Populate(FCIdump)

find_package(libtensor)
if (libtensor_FOUND)
    message(STATUS "using installed libtensor")
    target_compile_definitions(libtensor::tensor INTERFACE USE_PTHREADS)
else ()
    message(STATUS "building libtensor")
    FetchContent_Declare(libtensor
            #        GIT_REPOSITORY https://github.com/epifanovsky/libtensor
            GIT_REPOSITORY https://github.com/pjknowles/libtensor
            GIT_TAG ec5a2f7ae3643c5511e6c14537c481699670dac5)
    set(LIBTENSOR_STANDALONE 1)
    set(BUILD_SHARED_LIBS ON)
    FetchContent_MakeAvailable(libtensor)
    target_compile_definitions(tensor PUBLIC USE_PTHREADS)
    add_library(libtensor::tensor ALIAS tensor)
endif ()

target_link_libraries(${PROJECT_NAME} PUBLIC libtensor::tensor molpro::utilities molpro::Profiler molpro::FCIdump molpro::LinearAlgebra)

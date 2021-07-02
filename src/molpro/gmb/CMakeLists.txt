LibraryManager_Append(${PROJECT_NAME}
        PUBLIC_HEADER gmb.h
        SOURCES gmb.cpp init.cpp utils.cpp libtensor_utils.cpp get_integral.cpp
        )
LibraryManager_LAPACK(${PROJECT_NAME} PRIVATE Intel10_64lp)

DependencyManager_Populate(Profiler)
DependencyManager_Populate(utilities)
DependencyManager_Populate(LinearAlgebra)
DependencyManager_Populate(FCIdump)

find_package(libtensor QUIET)
if (libtensor_FOUND)
    message(STATUS "using installed libtensor")
    target_compile_definitions(libtensor::tensor INTERFACE USE_PTHREADS)
else ()
    if (GMB_COMPILE_LIBTENSOR)
        message(STATUS "building libtensor")
        FetchContent_Declare(libtensor
                #        GIT_REPOSITORY https://github.com/epifanovsky/libtensor
                GIT_REPOSITORY https://github.com/pjknowles/libtensor
                GIT_TAG 5a74fe3c
                )
        set(LIBTENSOR_STANDALONE 1)
        set(BUILD_SHARED_LIBS ON)
        FetchContent_MakeAvailable(libtensor)
        target_compile_definitions(tensor PUBLIC USE_PTHREADS)
        add_library(libtensor::tensor ALIAS tensor)
    else ()
        message(STATUS "not building libtensor")
        target_compile_definitions(${PROJECT_NAME} PUBLIC GMB_DUMMY)
    endif ()
endif ()

target_link_libraries(${PROJECT_NAME} PUBLIC libtensor::tensor molpro::utilities molpro::Profiler molpro::FCIdump molpro::LinearAlgebra)

set(src 
    main.cpp
    init.cpp
    libtensor_utils.cpp
    get_integral.cpp)
    
add_executable(gmb ${src})

DependencyManager_Populate(Profiler)
DependencyManager_Populate(utilities)
DependencyManager_Populate(LinearAlgebra)
DependencyManager_Populate(FCIdump)

include(FetchContent)
FetchContent_Declare(libtensor GIT_REPOSITORY https://github.com/epifanovsky/libtensor
        GIT_TAG 759566bfe96a4dea4f352be6e75cf126a66acc29)

set(LIBTENSOR_STANDALONE 1)
add_definitions(-DUSE_PTHREADS)
FetchContent_MakeAvailable(libtensor)


target_link_libraries(gmb tensor FCIdump molpro::LinearAlgebra)

set(tests He-VDZ H2O-VDZ Li-VDZ-UHF)

foreach(test ${tests})
file(COPY ../../../test/${test}.fcidump DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
endforeach()


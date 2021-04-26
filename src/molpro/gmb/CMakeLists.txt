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

target_link_libraries(gmb tensor FCIdump molpro::LinearAlgebra)

set(tests He-VDZ H2O-VDZ Li-VDZ-UHF)

foreach(test ${tests})
file(COPY ../../../test/${test}.fcidump DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
endforeach()


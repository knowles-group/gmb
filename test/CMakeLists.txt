include(FetchContent)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG "release-1.10.0"
)
FetchContent_MakeAvailable(googletest)

include(GoogleTest)

file(COPY He-VDZ H2O-VDZ Li-VDZ-UHF DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

set(gmb-tests 
    eom-ccsd-test)
    
foreach (test ${gmb-tests})
    add_executable(${test} ${test}.cpp)
    target_link_libraries(${test} PUBLIC gmock_main molpro::LinearAlgebra tensor FCIdump )
    add_test(NAME ${test} COMMAND ${test})
endforeach ()


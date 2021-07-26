include(FetchContent)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG "release-1.11.0"
)
FetchContent_MakeAvailable(googletest)

include(GoogleTest)

file(COPY He-VDZ H2O-VDZ Li-VDZ-UHF hubbard DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

# foreach (test eom-ccsd-test)
#     add_executable(${test} ${test}.cpp)
#     target_link_libraries(${test} PRIVATE gmock_main molpro::gmb)
#     add_test(NAME ${test} COMMAND ${test})
# endforeach ()

add_test(NAME hubbard-electronic-ccsd COMMAND main results=-1.438005955234,1.03130920079263 dump=hubbard/hubbard.fcidump)
add_test(NAME hubbard-weak-polaritonic-ccsd COMMAND main results=-1.437954983515,1.01859157154317 dump=hubbard/hubbard.fcidump polariton_modes=1 polariton_nmax=1 polariton_gamma=.01 polariton_omega=1.028])


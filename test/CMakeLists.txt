find_package(GTest QUIET)

if(NOT GTEST_FOUND)
    configure_file(
        ${CMAKE_SOURCE_DIR}/cmake/GoogleTest-CMakeLists.txt.in
        ${CMAKE_BINARY_DIR}/googletest-download/CMakeLists.txt
    )
    execute_process(
        COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download
    )

    if(result)
        message(FATAL_ERROR "CMake step for googletest failed: ${result}")
    endif()

    execute_process(
        COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download
    )

    if(result)
        message(FATAL_ERROR "Build step for googletest failed: ${result}")
    endif()

    add_subdirectory(
        ${CMAKE_BINARY_DIR}/googletest-src
        ${CMAKE_BINARY_DIR}/googletest-build
        EXCLUDE_FROM_ALL
    )
endif()

add_executable(hmm_test
    src/hmmintrin-test.cpp
    src/fft-test.cpp
    src/main.cpp)

target_link_libraries(hmm_test gtest hmm::hmm)

add_test(
    NAME hmm_test
    COMMAND hmm_test
)
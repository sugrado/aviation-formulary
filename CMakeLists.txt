cmake_minimum_required(VERSION 3.14)
project(aviation_formulary VERSION 1.47.0 LANGUAGES CXX)

# Header-only library target
add_library(aviation_formulary INTERFACE)
add_library(avform::aviation_formulary ALIAS aviation_formulary)

target_include_directories(aviation_formulary
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_compile_features(aviation_formulary INTERFACE cxx_std_17)

option(AVFORM_BUILD_TESTS "Build Aviation Formulary unit tests" ON)

if(AVFORM_BUILD_TESTS)
    enable_testing()
    add_executable(test_worked_examples tests/test_worked_examples.cpp)
    target_link_libraries(test_worked_examples PRIVATE aviation_formulary)
    if(MSVC)
        target_compile_options(test_worked_examples PRIVATE /W4)
    else()
        target_compile_options(test_worked_examples PRIVATE -Wall -Wextra -Wpedantic)
    endif()
    add_test(NAME WorkedExamplesTest COMMAND test_worked_examples)
endif()


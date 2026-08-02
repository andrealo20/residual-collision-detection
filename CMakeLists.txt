cmake_minimum_required(VERSION 3.16)
project(qcollide VERSION 0.1.0 LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

option(QCOLLIDE_BUILD_TESTS "Build unit tests" ON)
option(QCOLLIDE_WERROR      "Treat warnings as errors" ON)
option(QCOLLIDE_BENCH       "Build the model-error sweep" ON)

add_library(qcollide
    src/dynamics.c
    src/observer.c
    src/detector.c
    src/friction.c
    src/sensors.c
)

target_include_directories(qcollide PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(qcollide PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wstrict-prototypes)
    if(QCOLLIDE_WERROR)
        target_compile_options(qcollide PRIVATE -Werror)
    endif()
    target_link_libraries(qcollide PUBLIC m)
endif()

if(QCOLLIDE_BENCH)
    add_subdirectory(bench)
endif()

if(QCOLLIDE_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

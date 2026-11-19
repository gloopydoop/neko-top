# CMakeLists.txt for the tests directory.
# This file is used to define the tests for the Neko-TOP library.

# ROM integration diagnostics are explicit, artifact-driven probes rather than
# unconditional CTest cases. Keep their sources and targets under tests even
# when the optional pFUnit suite is unavailable.
option(BUILD_ROM_RUGBY_BALL_DIAGNOSTICS
    "Build explicit ROM_rugby_ball numerical diagnostics" ON)
if (BUILD_ROM_RUGBY_BALL_DIAGNOSTICS)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/rom_rugby_ball)
endif()

# This dense ROM derivative smoke test has no pFUnit dependency.  Keep it
# available with the ordinary library build so the reduced model can be tested
# independently of the optional pFUnit suite.
add_executable(rom_model_smoke EXCLUDE_FROM_ALL rom_model_smoke.f90)
target_link_libraries(rom_model_smoke PRIVATE Neko-TOP::neko-top)

if (NOT BUILD_TESTING)
    return()
endif()

add_custom_target(Tests
    COMMENT "Build all tests."
    DEPENDS neko-top
)
set_target_properties(Tests PROPERTIES
    EXCLUDE_FROM_ALL $<IF:$<BOOL:BUILD_TESTING>,FALSE,TRUE>
    FOLDER "Tests"
)
if (TARGET ROM_rugby_ball_diagnostics)
    add_dependencies(Tests ROM_rugby_ball_diagnostics)
endif()

# Set the shared test sources directory. This is where modules and other shared
# test files are located.
set(SHARED_TEST_SOURCES ${CMAKE_SOURCE_DIR}/tests/shared)

# ============================================================================ #
# Add subdirectories for unit and regression tests.

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/unit)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/regression)

add_test(NAME rom_model_smoke COMMAND rom_model_smoke)
set_tests_properties(rom_model_smoke PROPERTIES LABELS "unit;rom" TIMEOUT 30)
add_dependencies(Tests rom_model_smoke)

# ============================================================================ #
# Lint the registered tests: every DEPENDS entry must name an existing test.
# CTest silently ignores any other entry, e.g. a build target. The lint is
# also labelled "unit" so that it runs in CI's `ctest -L unit` selection; it
# takes well under a second and needs neither MPI nor a build.

add_test(NAME lint_test_depends
    COMMAND ${CMAKE_COMMAND}
        -DCTEST_COMMAND=${CMAKE_CTEST_COMMAND}
        -DTEST_DIR=${CMAKE_BINARY_DIR}
        -P ${CMAKE_CURRENT_SOURCE_DIR}/check_test_depends.cmake
)
set_tests_properties(lint_test_depends PROPERTIES
    LABELS "lint;unit"
    TIMEOUT 60
)

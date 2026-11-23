# ROM_rugby_ball integration diagnostics.
#
# These executables are intentionally kept out of the example directory and
# out of the default build.  They are explicit numerical probes that require
# a prepared case and/or streamed POD artifacts; they are not unconditional
# CTest smoke tests.

set(CMAKE_FOLDER "Tests/ROM/ROM_rugby_ball")

set(ROM_TEST_CORE_LIBS
    ${PROJECT_NAME}::neko-top
    PkgConfig::neko
    MPI::MPI_Fortran)

add_library(ROM_rugby_ball_test_helpers STATIC rom_array_test_adapters.f90)
target_link_libraries(ROM_rugby_ball_test_helpers PRIVATE ${ROM_TEST_CORE_LIBS})

function(add_rom_rugby_ball_diagnostic target source)
    add_executable(${target} EXCLUDE_FROM_ALL ${source})
    target_link_libraries(${target} PRIVATE ${ROM_TEST_CORE_LIBS}
        ROM_rugby_ball_test_helpers ${ARGN})
    set_target_properties(${target} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
endfunction()

add_rom_rugby_ball_diagnostic(
    ROM_rugby_ball_projection_diagnostic
    projection_diagnostic.f90
    PkgConfig::json-fortran
    $<$<BOOL:${BLAS_FOUND}>:BLAS::BLAS>
    $<$<BOOL:${LAPACK_FOUND}>:LAPACK::LAPACK>
    $<$<BOOL:${CUDAToolkit_FOUND}>:CUDA::cusolver>
    $<$<BOOL:${CUDAToolkit_FOUND}>:CUDA::cudart>
    $<$<BOOL:${hipblas_FOUND}>:roc::hipblas>
    $<$<BOOL:${hipsolver_FOUND}>:roc::hipsolver>)
target_sources(ROM_rugby_ball_projection_diagnostic PRIVATE
    ${PROJECT_SOURCE_DIR}/examples/ROM_low_Re_mini/user.f90)
if(TARGET neko_cxx_support)
    target_link_libraries(ROM_rugby_ball_projection_diagnostic PRIVATE
        neko_cxx_support)
endif()

add_rom_rugby_ball_diagnostic(
    ROM_rugby_ball_forward_continuity_diagnostic
    forward_continuity_diagnostic.f90
    PkgConfig::json-fortran
    $<$<BOOL:${BLAS_FOUND}>:BLAS::BLAS>
    $<$<BOOL:${LAPACK_FOUND}>:LAPACK::LAPACK>)

add_rom_rugby_ball_diagnostic(
    ROM_rugby_ball_projection_algebra_diagnostic
    projection_algebra_diagnostic.f90
    PkgConfig::json-fortran
    $<$<BOOL:${BLAS_FOUND}>:BLAS::BLAS>
    $<$<BOOL:${LAPACK_FOUND}>:LAPACK::LAPACK>
    $<$<BOOL:${CUDAToolkit_FOUND}>:CUDA::cusolver>
    $<$<BOOL:${CUDAToolkit_FOUND}>:CUDA::cudart>
    $<$<BOOL:${hipblas_FOUND}>:roc::hipblas>
    $<$<BOOL:${hipsolver_FOUND}>:roc::hipsolver>)
if(TARGET neko_cxx_support)
    target_link_libraries(ROM_rugby_ball_projection_algebra_diagnostic PRIVATE
        neko_cxx_support)
endif()

add_rom_rugby_ball_diagnostic(
    ROM_rugby_ball_objective_diagnostic
    objective_algebra_diagnostic.f90)

add_rom_rugby_ball_diagnostic(
    ROM_rugby_ball_scalar_transport_diagnostic
    scalar_transport_diagnostic.f90)

add_rom_rugby_ball_diagnostic(
    ROM_rugby_ball_scalar_objective_equivalence_diagnostic
    scalar_objective_equivalence_diagnostic.f90
    PkgConfig::json-fortran)
target_sources(ROM_rugby_ball_scalar_objective_equivalence_diagnostic PRIVATE
    ${PROJECT_SOURCE_DIR}/examples/ROM_low_Re_mini/user.f90)

add_rom_rugby_ball_diagnostic(
    ROM_rugby_ball_residual_identity_diagnostic
    residual_identity_diagnostic.f90
    PkgConfig::json-fortran)

add_rom_rugby_ball_diagnostic(
    ROM_rugby_ball_objective_equivalence_diagnostic
    objective_equivalence_diagnostic.f90
    PkgConfig::json-fortran)

add_rom_rugby_ball_diagnostic(
    ROM_rugby_ball_production_rom_gradient_diagnostic
    production_rom_gradient_diagnostic.f90
    PkgConfig::json-fortran)

add_rom_rugby_ball_diagnostic(
    ROM_rugby_ball_validation_diagnostic
    rom_validation_diagnostic.f90)

add_rom_rugby_ball_diagnostic(
    ROM_rugby_ball_production_rom_step_diagnostic
    production_rom_step_diagnostic.f90
    PkgConfig::json-fortran)

add_custom_target(ROM_rugby_ball_diagnostics
    DEPENDS
        ROM_rugby_ball_projection_diagnostic
        ROM_rugby_ball_forward_continuity_diagnostic
        ROM_rugby_ball_projection_algebra_diagnostic
        ROM_rugby_ball_objective_diagnostic
        ROM_rugby_ball_scalar_transport_diagnostic
        ROM_rugby_ball_scalar_objective_equivalence_diagnostic
        ROM_rugby_ball_residual_identity_diagnostic
        ROM_rugby_ball_objective_equivalence_diagnostic
        ROM_rugby_ball_production_rom_gradient_diagnostic
        ROM_rugby_ball_validation_diagnostic
        ROM_rugby_ball_production_rom_step_diagnostic)
set_target_properties(ROM_rugby_ball_diagnostics PROPERTIES
    FOLDER "Tests/ROM/ROM_rugby_ball")

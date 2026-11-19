# Dense reduced-order operators.  The module is independent of the Neko field
# layout so it can be used by the production projection layer and by small
# manufactured Fortran tests.
target_sources(neko-top
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/rom_model.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/rom_assembly.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/rom_basis_orthogonalization.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/rom_objectives.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/rom_projection_operator.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/rom_basis_handoff.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/rom_problem.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/rom_validation.f90
)

# Compile the debug example with the same user-defined setup as low_Re.
set(DRIVER_TYPE "topopt-user")
set(EXTRA_SOURCES
    ${EXAMPLES_DIR}/low_Re/user.f90
)

build_example(${DRIVER_TYPE})

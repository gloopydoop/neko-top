# Compile the debug example with the same user-defined setup as POD_mixer.
set(DRIVER_TYPE "topopt-user")
set(EXTRA_SOURCES
    ${EXAMPLES_DIR}/POD_mixer/user.f90
)

build_example(${DRIVER_TYPE})

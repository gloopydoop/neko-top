# DMD uses the standard Neko-TOP topology-optimization driver; the state
# recovery implementation is selected by the case file.
set(DRIVER_TYPE "topopt")

build_example(${DRIVER_TYPE})

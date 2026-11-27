# Minimal coupled velocity/passive-scalar ROM experiment.
set(DRIVER_TYPE "topopt-user")
set(EXTRA_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/user.f90
)

build_example(${DRIVER_TYPE})

# Some Neko installations declare the SST timeout by reference in Fortran
# but by value in C++. Adapt that specific static-library ABI at the link
# boundary without modifying the installed Neko sources or library.
if(HAVE_ADIOS2 AND CMAKE_SYSTEM_NAME STREQUAL "Linux" AND
   EXISTS "${NEKO_DIR}/lib/libneko.a" AND
   EXISTS "${NEKO_DIR}/src/io/nek_adios2.cpp" AND
   EXISTS "${NEKO_DIR}/src/io/data_streamer.F90")
    file(READ "${NEKO_DIR}/src/io/nek_adios2.cpp" neko_adios_source)
    file(READ "${NEKO_DIR}/src/io/data_streamer.F90" neko_stream_source)
    if(neko_adios_source MATCHES "const int timeout_seconds" AND
       neko_stream_source MATCHES "integer\\(kind=C_INT\\) :: timeout")
        target_sources(ROM_low_Re_mini PRIVATE adios_timeout_compat.f90)
        target_link_options(ROM_low_Re_mini PRIVATE "LINKER:--wrap=adios2_initialize_")
        message(STATUS "ROM_low_Re_mini: adapting Neko ADIOS timeout ABI")
    endif()
endif()

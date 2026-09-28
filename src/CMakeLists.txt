add_library(
    amplitude
    STATIC
    amplitudelib.cpp
    datafile.cpp
    tools.cpp
    interpolation.cpp
)

add_executable(
    D0
    integrand.cpp
    int_exclusive.cpp
    int_diffractive.cpp
    alphas_running.cpp
    photon_flux.cpp
    fragmentation.cpp
    hymnd_grid.cpp
    bcfy_grid.cpp
    kk_grid.cpp
    main.cpp
    fourier.c
)

target_link_libraries(
    D0
    PRIVATE
    amplitude
    GSL::gsl
    GSL::gslcblas
    OpenMP::OpenMP_CXX
)

add_executable(
    D0_xpom
    integrand.cpp
    int_exclusive.cpp
    int_diffractive.cpp
    alphas_running.cpp
    photon_flux.cpp
    fragmentation.cpp
    hymnd_grid.cpp
    bcfy_grid.cpp
    kk_grid.cpp
    main_xpom.cpp
    fourier.c
)

target_link_libraries(
    D0_xpom
    PRIVATE
    amplitude
    GSL::gsl
    GSL::gslcblas
    OpenMP::OpenMP_CXX
)

add_executable(
    D0_fixedW
    integrand.cpp
    int_exclusive.cpp
    int_diffractive.cpp
    alphas_running.cpp
    photon_flux.cpp
    fragmentation.cpp
    hymnd_grid.cpp
    main_fixedW.cpp
    fourier.c
)

target_link_libraries(
    D0_fixedW
    PRIVATE
    amplitude
    GSL::gsl
    GSL::gslcblas
    OpenMP::OpenMP_CXX
)

add_executable(
    D0_fixedW_xpom
    integrand.cpp
    int_exclusive.cpp
    int_diffractive.cpp
    alphas_running.cpp
    photon_flux.cpp
    fragmentation.cpp
    hymnd_grid.cpp
    main_fixedW_xpom.cpp
    fourier.c
)

target_link_libraries(
    D0_fixedW_xpom
    PRIVATE
    amplitude
    GSL::gsl
    GSL::gslcblas
    OpenMP::OpenMP_CXX
)

# Diagnostic: b-integrated photon flux vs. the tabulated fluxes of
# arXiv:2404.09731 (python/flux_comparison.py). Not part of the physics pipeline.
add_executable(
	scan_flux
	scan_flux.cpp
	photon_flux.cpp
)

target_link_libraries(
	scan_flux
	PRIVATE
	GSL::gsl
	GSL::gslcblas
)

set_target_properties(scan_flux PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Diagnostic: tests the "effective flux" geometric convolution (Eq. 4 of
# arXiv:2404.09731) against Paakkinen's tabulated flux. Not part of the
# physics pipeline.
add_executable(
	scan_flux_eff
	scan_flux_eff.cpp
	photon_flux.cpp
)

target_link_libraries(
	scan_flux_eff
	PRIVATE
	GSL::gsl
	GSL::gslcblas
)

set_target_properties(scan_flux_eff PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

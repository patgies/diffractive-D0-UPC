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
    charm_fixed_qp
    integrand.cpp
    int_exclusive.cpp
    int_diffractive.cpp
    alphas_running.cpp
    photon_flux.cpp
    main_fixed_qp.cpp
    fourier.c
)

target_link_libraries(
    charm_fixed_qp
    PRIVATE
    amplitude
    GSL::gsl
    GSL::gslcblas
    OpenMP::OpenMP_CXX
)

# Test tool: photon flux integrated over b, to compare with the tables of
# arXiv:2404.09731 (python/flux_comparison.py). Not used for the physics results.
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

# Test tool: effective flux (Eq. 4 of arXiv:2404.09731) for
# python/flux_comparison.py. Not used for the physics results.
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

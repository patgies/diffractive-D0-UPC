add_library(
    amplitude
    STATIC
    amplitudelib.cpp
    datafile.cpp
    tools.cpp
    interpolation.cpp
)

# Cross sections and photon flux
set(CROSS_SECTION_SOURCES
    integrands.cpp
    cross_section_exclusive.cpp
    cross_section_diffractive.cpp
    vegas.cpp
    photon_flux.cpp
    fourier.c
)

# D0 settings and fragmentation
set(D0_SOURCES
    ${CROSS_SECTION_SOURCES}
    setup.cpp
    lhapdf_grid.cpp
    bcfy_grid.cpp
    kk_grid.cpp
)

set(CROSS_SECTION_LIBRARIES amplitude GSL::gsl GSL::gslcblas OpenMP::OpenMP_CXX)

# D0 at one (pD0, y)
add_executable(D0 ${D0_SOURCES} main.cpp)
target_link_libraries(D0 PRIVATE ${CROSS_SECTION_LIBRARIES})

# Same at fixed x_po
add_executable(D0_xpom ${D0_SOURCES} main_xpom.cpp)
target_link_libraries(D0_xpom PRIVATE ${CROSS_SECTION_LIBRARIES})

# Charm quark at fixed q+, vs pT
add_executable(charm_fixed_qp ${CROSS_SECTION_SOURCES} main_fixed_qp.cpp)
target_link_libraries(charm_fixed_qp PRIVATE ${CROSS_SECTION_LIBRARIES})

# Same vs x_po
add_executable(charm_fixed_qp_xpom ${CROSS_SECTION_SOURCES} main_fixed_qp_xpom.cpp)
target_link_libraries(charm_fixed_qp_xpom PRIVATE ${CROSS_SECTION_LIBRARIES})

# Photon flux alone (flux_comparison.py)
add_executable(scan_flux scan_flux.cpp photon_flux.cpp)
target_link_libraries(scan_flux PRIVATE GSL::gsl GSL::gslcblas)

# Saturation scale vs x (correlation_limit.py)
add_executable(saturation_scale saturation_scale.cpp fourier.c)
target_link_libraries(saturation_scale PRIVATE amplitude GSL::gsl GSL::gslcblas)

set_target_properties(D0 D0_xpom charm_fixed_qp charm_fixed_qp_xpom scan_flux saturation_scale PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

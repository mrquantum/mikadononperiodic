project(mikado)

SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wreturn-type -Wunused-variable -std=c++11")

add_executable(mikado cutsprings.cpp randomnetwork.cpp makeSpringNetworks.cpp cgwoagony.cpp shearleader.cpp clusters.cpp writefunctions.cpp nrutil.cpp BendingGrad.cpp BendingEnergy.cpp prestress.cpp minimizers.cpp EnergyandGradients.cpp makemikadonetwork.cpp main.cpp)

install(TARGETS mikado RUNTIME DESTINATION bin)

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <armadillo>

#include "cardiac_sat.hpp"
#include "numpy_wrappers.hpp"

#define ASSERT_THROW(a,msg) if (!(a)) throw std::runtime_error(msg);

namespace py = pybind11;

/* module function */
py::array BlochCardiacSat_wrap(const double sliceThickness,
			  const unsigned int nValues,
			  const unsigned int nExcitations,
			  const double T1,
			  const double T2,
			  const double heartRate,
			  const unsigned int returnType)
    
{
    arma::cx_vec signal = BlochCardiacSat(sliceThickness, nValues, nExcitations,
    					  T1, T2, heartRate, returnType);

    py::array result = CArma2CNumpy(signal);

    return result;
    
}

PYBIND11_MODULE(bloch_cardiac, m)
{
    m.doc() = "Cardiac bloch simulator for the DiBella cardiac perfusion pulse sequence";

    m.def("BlochCardiacSat",&BlochCardiacSat_wrap,"The cardiac function...");
}

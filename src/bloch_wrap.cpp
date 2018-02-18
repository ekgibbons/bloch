#include <complex>
#include <stdexcept>
#include <iostream>
#include <vector> 
#include <algorithm>

#include <pybind11/pybind11.h>

#include <armadillo>

#include "cardiac_sat.hpp"
#include "numpy_wrappers.hpp"

#define ASSERT_THROW(a,msg) if (!(a)) throw std::runtime_error(msg);

#if _MSC_VER
using boost::uint8_t;
#endif

namespace py = pybind11;

/* module function */
np::ndarray BlochCardiacSat_wrap(const double sliceThickness,
				 const unsigned int nValues,
				 const unsigned int nExcitations,
				 const double T1,
				 const double T2, const double heartRate,
				 const unsigned int returnType)

{

    arma::cx_vec signal = BlochCardiacSat(sliceThickness, nValues, nExcitations,
					  T1, T2, heartRate, returnType);

    np::ndarray result = Arma2Numpy(signal);

    return result;
    
}

PYBIND11_MODULE(bloch_cardiac, m)
{
    m.doc() = "Cardiac bloch simulator for the DiBella cardiac perfusion pulse sequence";

    m.def("BlochCardiacSat",&BlochCardiacSat_wrap,"The cardiac function...");
}

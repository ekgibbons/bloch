#include <complex>
#include <stdexcept>
#include <iostream>
#include <vector> 
#include <algorithm>

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

#include <armadillo>

#include "cardiac_sat.hpp"
#include "numpy_wrappers.hpp"

#define ASSERT_THROW(a,msg) if (!(a)) throw std::runtime_error(msg);

#if _MSC_VER
using boost::uint8_t;
#endif


namespace py = boost::python;
namespace np = boost::python::numpy;

/* module function */
np::ndarray BlochCardiacSat_wrap(const double sliceThickness, const unsigned int nValues,
				 const unsigned int nExcitations, const double T1,
				 const double T2, const unsigned int returnType)
{

    arma::cx_vec signal = BlochCardiacSat(sliceThickness, nValues, nExcitations,
					  T1, T2, returnType);

    np::ndarray result = Arma2Numpy(signal);

    return result;
    
}

BOOST_PYTHON_MODULE(bloch_cardiac)
{
    Py_Initialize();
    np::initialize();
    py::def("BlochCardiacSat",BlochCardiacSat_wrap);
}

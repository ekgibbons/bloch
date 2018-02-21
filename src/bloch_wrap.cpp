#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <armadillo>

#include "echotrain_cpmg.hpp"
#include "echotrain_alsop.hpp"
#include "echotrain_gibbons.hpp"
#include "echotrain_lrx.hpp"
#include "cardiac_sat.hpp"
#include "numpy_wrappers.hpp"


#define ASSERT_THROW(a,msg) if (!(a)) throw std::runtime_error(msg);

namespace py = pybind11;


/* module function */
py::array BlochFSEGibbons_wrap(py::array_t<double> flipsNumpy, unsigned int nValues,
			       unsigned int etl, double initialPhase, unsigned int returnType,
			       double T1, double T2, double esp)
{

    arma::vec flips = NumpyVec2ArmaVec(flipsNumpy);

    arma::cx_vec signal = BlochFSEGibbons(flips, nValues, etl, initialPhase, returnType, T1, T2, esp);

    py::array result = CArma2CNumpy(signal);

    return result;
}

/* module function */
py::array BlochFSEAlsop_wrap(py::array_t<double> flipsNumpy, unsigned int nValues,
			     unsigned int etl, double initialPhase, unsigned int returnType,
			     double T1, double T2, double esp)
{

    arma::vec flips = NumpyVec2ArmaVec(flipsNumpy);

    arma::cx_vec signal = BlochFSEAlsop(flips, nValues, etl, initialPhase, returnType, T1, T2, esp);

    py::array result = CArma2CNumpy(signal);

    return result;
}

/* module function */
py::array BlochFSELRX_wrap(unsigned int nValues, unsigned int etl, double initialPhase,
			   unsigned int returnType, double T1, double T2, double esp)
{

    arma::cx_vec signal = BlochFSELRX(nValues, etl, 2, initialPhase, returnType, T1, T2, esp);

    py::array result = CArma2CNumpy(signal);

    return result;
}


/* module function */
py::array BlochFSECPMG_wrap(py::array_t<double> flipsNumpy, unsigned int nValues,
			    unsigned int etl, unsigned int returnType,
			    double T1, double T2, double esp)
{

    arma::vec flips = NumpyVec2ArmaVec(flipsNumpy);

    arma::cx_vec signal = BlochFSECPMG(flips, nValues, etl, returnType, T1, T2, esp);

    py::array result = CArma2CNumpy(signal);

    return result;
}

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

PYBIND11_MODULE(blochsequence, m)
{
    m.doc() = R"pbdoc(
This is a library of various sequnces that are simulated through Bloch numerical methods.
				     )pbdoc";

    m.def("CardiacSat",&BlochCardiacSat_wrap,R"pbdoc(
Cardiac bloch simulator for the DiBella cardiac perfusion pulse sequence.  This sequence is essentially is a saturation pulse followed by a small flip angle readout during recovery.  This sequence assumes a beating heart (with the user-defined heart rate).  The output is the magnetization profiles--particularly the slice profile.

Arguments
---------
sliceThickness: float
    The slice thickness to use in meters
nValues: int
    The number of isochromats to use
nExcitations: int
    The number of excitations in the readout train
T1: float
    The T1 recovery constant in seconds
T2: float
    The T2 relaxation constant in seconds
heartRate: float
    Heart rate in beats per minute
returnType: int
    The type of magentization to return (0: echotrain, 1: Mxy, 2: Mz)

Returns
-------
output:  array_like
    The magentization signal based on what the returnType is
)pbdoc",
	  py::arg("sliceThickness"),py::arg("nValues"),py::arg("nExcitations"),
	  py::arg("T1"),py::arg("T2"),py::arg("heartRate"),py::arg("returnType")
	);

    m.def("FSEGibbons",&BlochFSEGibbons_wrap,"Gibbons");
    m.def("FSEAlsop",&BlochFSEAlsop_wrap,"Alsop");
    m.def("FSELRX",&BlochFSEGibbons_wrap,"LRX");
    m.def("FSECPMG",&BlochFSECPMG_wrap,"CPMG");
    
}

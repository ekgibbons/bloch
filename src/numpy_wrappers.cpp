/* numpy_wrappers.cpp
 * 
 * Just some basic wrapping functions to move numpy
 * formatted array to standard C/C++ arrays
 *
 */


#include <complex>
#include <stdexcept>
#include <iostream>
#include <vector> 
#include <algorithm>

#include <math.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <armadillo>

#include "numpy_wrappers.hpp"

namespace py = pybind11;

py::array CArma2CNumpy(const arma::cx_vec &inputVec)
{

    std::vector<std::complex<double>> tempVector(inputVec.n_elem);
    std::copy(inputVec.begin(),inputVec.end(),tempVector.begin());

    std::vector<ssize_t> shape = {(long int)inputVec.n_elem};
    std::vector<ssize_t> strides = {sizeof(std::complex<double>)};

    return py::array(py::buffer_info(
			 tempVector.data(),                                     /* data as contiguous array  */
    			 inputVec.n_elem*sizeof(std::complex<double>),          /* size of one scalar        */
    			 py::format_descriptor<std::complex<double>>::format(), /* data type                 */
    			 1,                                                     /* number of dimensions      */
    			 shape,                                                 /* shape of the matrix       */
    			 strides                                                /* strides for each axis     */
    			 )
	);
    
}


arma::vec NumpyVec2ArmaVec(py::array_t<double> input)
{
    auto r = input.unchecked<1>();

    arma::vec inputArma(input.shape(0));

    for (unsigned int ii = 0; ii < input.shape(0); ii++)
    {
	inputArma(ii) = r(ii);

    }

    return inputArma;
}




// np::ndarray Arma2Numpy(const arma::cx_mat &inputMat)
// {

//     py::tuple shape = py::make_tuple(inputMat.n_rows,inputMat.n_cols);
//     unsigned int vectorLength = ((unsigned int)(inputMat.n_rows))*(unsigned int)(inputMat.n_cols);
//     std::vector<std::complex<double> > tempVector(vectorLength);
//     std::copy(inputMat.begin(),inputMat.end(),tempVector.begin());

//     np::dtype dt1 = np::dtype::get_builtin<std::complex<double> >(); 
//     np::ndarray result = np::zeros(shape, dt1);

//     unsigned int n = 0;
//     for (unsigned int jj = 0; jj < (unsigned int)inputMat.n_cols; jj++)
//     {
// 	for (unsigned int ii = 0; ii < (unsigned int)inputMat.n_rows; ii++)
// 	{
//     	    result[ii][jj] = (std::complex<double>)tempVector[n];
// 	    n++;
// 	}
//     }
//     return result;
// }


// arma::vec Numpy2ArmaVec(const np::ndarray &input)
// {
//     unsigned int lengthInput = (unsigned int)input.shape(0);
//     std::vector<double> tempVector(lengthInput);

//     std::copy(reinterpret_cast<double*>(input.get_data()), 
// 	      reinterpret_cast<double*>(input.get_data())+lengthInput,
// 	      tempVector.begin());

//     arma::vec output(lengthInput);

//     for (unsigned int ii = 0; ii < lengthInput; ii++)
//     {
// 	output(ii) = (double)tempVector[ii];
//     }

//     return output;
// }

// np::ndarray Array2Numpy(double *array, int length)
// {
//     py::tuple shape = py::make_tuple(length);
//     np::dtype dt1 = np::dtype::get_builtin<double>(); 
//     np::ndarray result = np::zeros(shape, dt1);

//     for (unsigned int ii = 0; ii < (unsigned int)length; ii++)
//     {
// 	result[ii] = array[ii];
//     }
    
//     return result;
// }


// np::ndarray CArray2CNumpy(double *array, int length)
// {
//     std::complex<double> j = std::complex<double>(0,1.);
    
//     py::tuple shape = py::make_tuple(length);
//     np::dtype dt1 = np::dtype::get_builtin<std::complex<double> >(); 
//     np::ndarray result = np::zeros(shape, dt1);

//     for (unsigned int ii = 0; ii < (unsigned int)length; ii++)
//     {
// 	result[ii] = array[2*ii] + j*array[2*ii + 1];
//     }
    
//     return result;
// }

// void Numpy2Array(const np::ndarray &arrayNumpy,
// 		 double *arrayOut,
// 		 int length)
// {
//     std::vector<double> tempVector(length);
    
//     std::copy(reinterpret_cast<double*>(arrayNumpy.get_data()), 
// 	      reinterpret_cast<double*>(arrayNumpy.get_data())+length,
// 	      tempVector.begin());

//     for (unsigned int ii = 0; ii < (unsigned int) length; ii++)
//     {
// 	arrayOut[ii] = (double)tempVector[ii];
//     }
// }

// void Numpy2CArray(const np::ndarray &arrayNumpy,
// 		  double *arrayOut,
// 		  int length)
// {
//     std::vector<double> tempVector(length);
    
//     std::copy(reinterpret_cast<double*>(arrayNumpy.get_data()), 
// 	      reinterpret_cast<double*>(arrayNumpy.get_data())+length,
// 	      tempVector.begin());

//     for (unsigned int ii = 0; ii < (unsigned int) length; ii++)
//     {
// 	arrayOut[2*ii] = (double)tempVector[ii];
//     }
// }


// void CNumpy2CArray(const np::ndarray &arrayNumpy,
// 		  double *arrayOut,
// 		  int length)
// {
//     std::vector<std::complex<double> > tempVector(length);
    
//     std::copy(reinterpret_cast<std::complex<double>* >(arrayNumpy.get_data()), 
// 	      reinterpret_cast<std::complex<double>* >(arrayNumpy.get_data())+length,
// 	      tempVector.begin());

//     for (unsigned int ii = 0; ii < (unsigned int) length; ii++)
//     {
// 	arrayOut[2*ii] = (double)tempVector[ii].real();
// 	arrayOut[2*ii+1] = (double)tempVector[ii].imag();
//     }
// }

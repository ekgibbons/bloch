/* numpy_wrappers.hpp
 * 
 * Just some basic wrapping functions to move numpy
 * formatted array to standard C/C++ arrays
 *
 */


#ifndef __NUMPY_WRAPPERS_HPP__
#define  __NUMPY_WRAPPERS_HPP__

namespace py = pybind11;

py::array CArma2CNumpy(const arma::cx_vec &inputVec);
arma::vec NumpyVec2ArmaVec(py::array_t<double> input);
// np::ndarray Arma2Numpy(const arma::cx_mat &inputMat);
// np::ndarray Arma2Numpy(const arma::cx_vec &inputVec);
// arma::vec Numpy2ArmaVec(const np::ndarray &input);
// np::ndarray Array2Numpy(double *array, int length);
// np::ndarray CArray2CNumpy(double *array, int length);
// void Numpy2Array(const np::ndarray &arrayNumpy,
// 		 double *arrayOut,
// 		 int length);
// void Numpy2CArray(const np::ndarray &arrayNumpy,
// 		  double *arrayOut,
// 		      int length);
// void CNumpy2CArray(const np::ndarray &arrayNumpy,
// 		  double *arrayOut,
// 		  int length);

#endif

#ifndef cardiac_sat_hpp
#define cardiac_sat_hpp

#define ECHOTRAIN 0
#define MXY 1
#define MZ 2

arma::cx_vec BlochCardiacSat(const double sliceThickness,
			     const unsigned int nValues,
			     const unsigned int nExcitations, const double T1,
			     const double T2, double heartRate,
			     const unsigned int returnType);

# endif

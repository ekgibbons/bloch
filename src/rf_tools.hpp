#ifndef rf_tools_hpp
#define rf_tools_hpp

void abrx(const arma::vec &rf, const arma::vec &gx, const arma::vec &x, double phi,
	  arma::cx_vec &ALF, arma::cx_vec &BET);

void abrot(double (&a)[2], double (&b)[2], double x, double phiOriginal, 
	   const arma::vec &rf, 
	   const arma::vec &gx);

arma::mat abfm2m(const arma::cx_vec &a, const arma::cx_vec &b, arma::vec &m);

void abfm2m(std::complex<double> a, std::complex<double> b, arma::mat &m);

#endif

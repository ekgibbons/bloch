#ifndef rf_pulse_hpp
#define rf_pulse_hpp

arma::vec GetMSINC(unsigned int nPoints, double TBW,
		   double flipAngle);
arma::vec GetSLR(void);
arma::vec GetDIVERSE(void);
arma::vec GetDIVERSEGrad(void);
arma::vec GetSLRExcite(void);
arma::vec GetFourier(void);
arma::vec GetSLRLow(void);
arma::vec GetSLRHigh(void);
arma::vec GetSLRSE(void);

#endif

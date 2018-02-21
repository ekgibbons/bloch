#ifndef echotrain_lrx_hpp
#define echotrain_lrx_hpp

#define ECHOTRAIN 0
#define MXY 1
#define MZ 2

arma::cx_vec BlochFSELRX(unsigned int nValues, unsigned int etl,  unsigned int pulseType,
			 double initialPhase, unsigned int returnType, double T1, double T2,
			 double esp);

#endif

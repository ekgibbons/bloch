#ifndef echotrain_alsop_hpp
#define echotrain_alsop_hpp

#define ECHOTRAIN 0
#define MXY 1
#define MZ 2

arma::cx_vec BlochFSEAlsop(const arma::vec &flips, unsigned int nValues, unsigned int etl,
			   double initialPhase, unsigned int returnType, double T1, double T2,
			   double esp);
# endif

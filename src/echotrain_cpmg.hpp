#ifndef echotrain_cpmg_hpp
#define echotrain_cpmg_hpp

#define ECHOTRAIN 0
#define MXY 1
#define MZ 2

arma::cx_vec BlochFSECPMG(const arma::vec &flips, unsigned int nValues, unsigned int etl,
			  unsigned int returnType, double T1, double T2, double esp);

# endif

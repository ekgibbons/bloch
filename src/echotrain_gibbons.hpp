#ifndef echotrain_gibbons_hpp
#define echotrain_gibbons_hpp

#define ECHOTRAIN 0
#define MXY 1
#define MZ 2

arma::cx_vec BlochFSEGibbons(const arma::vec &flips, unsigned int nValues, unsigned int etl,
			     double initialPhase, unsigned int returnType, double T1, double T2,
			     double esp);
# endif

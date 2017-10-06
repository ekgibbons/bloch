#ifndef bloch_utilities_hpp
#define bloch_utilities_hpp

const arma::mat RotationX(double alpha);
const arma::mat RotationY(double alpha);
const arma::mat RotationZ(double alpha);
const arma::mat RotationArbitrary(double alpha, double theta);
void HardPrecession(arma::mat &M, const arma::vec &z, 
		    const double anglePrecession);
void HardPrecessionY(arma::mat &M, const arma::vec &z, 
		     const double anglePrecession);
void HardPrecessionX(arma::mat &M, const arma::vec &z, 
		     const double anglePrecession);
void Relaxation(arma::mat & A, arma::vec & b,
		double t,double T1, double T2);
void SLRRotation(arma::mat &m, const arma::vec x,
		 const arma::vec rf, const unsigned int rfType, const double phase,
		 const double dur, const double tb, const double thk);
void SLRRotation(arma::mat &m, const double x,
		 const arma::vec rf, const unsigned int rfType, const double phase,
		 const double dur, const double tb, const double thk);
void Diffusion(arma::mat &M, arma::vec &z, arma::vec &G,
	       double D, double dT);
const arma::mat Spoiling(void);

# endif

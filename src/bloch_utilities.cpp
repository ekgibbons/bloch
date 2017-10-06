#include <cmath>

#include <armadillo>

#include "rf_tools.hpp"
#include "rf_pulses.hpp"

#include "bloch_utilities.hpp"

const arma::mat RotationX(double alpha)
{
    arma::mat A = arma::zeros<arma::mat>(3,3);
    A(0,0) = 1;
    A(1,1) = cos(alpha);
    A(1,2) = sin(alpha);
    A(2,1) = -sin(alpha);
    A(2,2) = cos(alpha);
    
    return A;
}

const arma::mat RotationY(double alpha)
{
    arma::mat A = arma::zeros<arma::mat>(3,3);
    A(0,0) = cos(alpha);
    A(0,2) = sin(alpha);
    A(1,1) = 1;
    A(2,0) = -sin(alpha);
    A(2,2) = cos(alpha);
    
    return A;
}

const arma::mat RotationZ(double alpha)
{
    arma::mat A = arma::zeros<arma::mat>(3,3);
    A(0,0) = cos(alpha);
    A(0,1) = sin(alpha);
    A(1,0) = -sin(alpha);
    A(1,1) = cos(alpha);
    A(2,2) = 1;
    
    return A;
}

const arma::mat RotationArbitrary(double alpha, double theta)
{
    double u = cos(theta);
    double v = sin(theta);
    
    arma::mat A = arma::zeros<arma::mat>(3,3);
    A(0,0) = pow(u,2)*(1-cos(alpha)) + cos(alpha);
    A(0,1) = u*v*(1 - cos(alpha));
    A(0,2) = v*sin(alpha);
    A(1,0) = u*v*(1 - cos(alpha));
    A(1,1) = cos(alpha) + pow(v,2)*(1 - cos(alpha));
    A(1,2) = -u*sin(alpha);
    A(2,0) = -v*sin(alpha);
    A(2,1) = u*sin(alpha);
    A(2,2) = cos(alpha);

    return A;		  

}

void Relaxation(arma::mat & A, arma::vec & b,
		double t,double T1, double T2)
{
    double E1 = exp(-t/T1);
    double E2 = exp(-t/T2);

    A = arma::zeros<arma::mat>(3,3);
    A(0,0) = E2;
    A(1,1) = E2;
    A(2,2) = E1;
    
    b = arma::zeros<arma::vec>(3);
    b(2) = 1 - E1;
}

void HardPrecession(arma::mat &M, const arma::vec &z, 
		    const double anglePrecession)
{

    unsigned int nValues = (unsigned int)z.n_elem;
    double gamma2pi = 42.58; // kHz/mT

    for (unsigned int zIndex = 0; zIndex < nValues; zIndex++)
    {
	arma::mat crusher = RotationZ(gamma2pi*anglePrecession*z(zIndex));
	M.col(zIndex) = crusher*M.col(zIndex);
    }
    
}

void HardPrecessionY(arma::mat &M, const arma::vec &z, 
		    const double anglePrecession)
{

    unsigned int nValues = (unsigned int)z.n_elem;
    double gamma2pi = 42.58; // kHz/mT

    for (unsigned int zIndex = 0; zIndex < nValues; zIndex++)
    {
	arma::mat crusher = RotationX(gamma2pi*anglePrecession*z(zIndex));
	M.col(zIndex) = crusher*M.col(zIndex);
    }
    
}

void HardPrecessionX(arma::mat &M, const arma::vec &z, 
		    const double anglePrecession)
{

    unsigned int nValues = (unsigned int)z.n_elem;
    double gamma2pi = 42.58; // kHz/mT

    for (unsigned int zIndex = 0; zIndex < nValues; zIndex++)
    {
	arma::mat crusher = RotationX(gamma2pi*anglePrecession*z(zIndex));
	M.col(zIndex) = crusher*M.col(zIndex);
    }
    
}

void SLRRotation(arma::mat &m, const arma::vec x,
		 const arma::vec rf, const unsigned int rfType, const double phase,
		 const double dur, const double tb, const double thk)
{
    std::complex<double> j(0,1);

    double gamma = 4258;
    double bw = tb/dur;
    double a_g = bw / (thk * gamma);
    unsigned int nRF = (unsigned int)rf.n_elem;


    arma::vec gRF;
    if (rfType == 2)
    {
	double nominalThickness = 0.5; // cm
	gRF = GetDIVERSEGrad();
	gRF = nominalThickness*gRF/thk;
    }
    else 
    {
	gRF = a_g*arma::ones<arma::vec>(nRF);
    }    

    gRF = 2*M_PI*gamma*gRF*dur/nRF;

    for (unsigned int xLocation = 0 ; xLocation < x.n_elem; xLocation++)
    {
	std::complex<double> alpha;
	std::complex<double> beta;

	arma::mat mTemp = m.col(xLocation);

	double alf[2], bet[2];
	alf[0] = 1.0;
	alf[1] = 0.0; 
	bet[0] = 0.0; 
	bet[1] = 0.0;
	abrot(alf, bet, x(xLocation), phase, rf, gRF);
	alpha = alf[0] + j*alf[1];
	beta = bet[0] + j*bet[1];

	abfm2m(alpha, beta, mTemp);

	m.col(xLocation) = mTemp;

    }
}

void SLRRotation(arma::mat &m, const double x,
		 const arma::vec rf, const unsigned int rfType,  const double phase,
		 const double dur, const double tb, const double thk)
{

    std::complex<double> j(0,1);

    double gamma = 4258;
    double bw = tb/dur;
    double a_g = bw / (thk * gamma);
    unsigned int nRF = (unsigned int)rf.n_elem;

    arma::vec gRF;
    if (rfType == 2)
    {
	gRF = GetDIVERSEGrad();
    }
    else 
    {
	gRF = a_g*arma::ones<arma::vec>(nRF);
    }    
    
    gRF = 2*M_PI*gamma*gRF*dur/nRF;

    std::complex<double> alpha;
    std::complex<double> beta;
   
    double alf[2], bet[2];
    alf[0] = 1.0;
    alf[1] = 0.0; 
    bet[0] = 0.0; 
    bet[1] = 0.0;
    abrot(alf, bet, x, phase, rf, gRF);
    alpha = alf[0] + j*alf[1];
    beta = bet[0] + j*bet[1];
    
    abfm2m(alpha, beta, m);    
}

void Diffusion(arma::mat &M, arma::vec &z, arma::vec &G,
	       double D, double dT)
{
    double standardDeviation = sqrt(2*D*dT);
    arma::vec zAdd(z.n_elem);

    // rotation from the gradient+diffusion
    for (unsigned int g = 0; g < G.n_elem; g++)
    {
	z = z + standardDeviation*zAdd.randn();
	
	double gradArea = G(g)*dT;
	HardPrecession(M, z, 2*M_PI*gradArea);
    }
}

const arma::mat Spoiling(void)
{
    arma::mat A = arma::zeros<arma::mat>(3,3);
    A(2,2) = 1;

    return A;
}

#include <complex>
#include <cmath>
#include <string>

#include <armadillo>

#include "bloch_utilities.hpp"
#include "rf_pulses.hpp"

#include "echotrain_alsop.hpp"

arma::cx_vec BlochFSECPMG(const arma::vec &flips, unsigned int nValues, unsigned int etl,
			  unsigned int returnType, double T1, double T2, double esp)
{

    printf("SIMULATION STARTING: cpmg\n");

    std::complex<double> j(0,1);
    
    double gamma2pi = 42.58; // kHz/mT
    double sliceThickness = 0.005; // m

    double initialPhase = M_PI/2; // phase in rads
    
    // double TBWExcite = 3.55;
    // double pulseWidthExcite = 1.200e-3; // s
    // double BWExcite = TBWExcite/pulseWidthExcite;
    // arma::vec rfExcite = GetSLRExcite();
    // unsigned int pulseTypeExcite = 1;

    unsigned int nPointsRF = 320;

    double TBWHide = 2.;
    double pulseWidthHide = 1.288e-3; // s
    double BWHide = TBWHide/pulseWidthHide;
    double flipAngleMSINC = 90.0;
    arma::vec rfHide = GetMSINC(nPointsRF,TBWHide,flipAngleMSINC);
    unsigned int pulseTypeHide = 1;
    
    arma::mat M = arma::zeros<arma::mat>(3,nValues);
    arma::cx_vec Mxy(nValues);
    arma::cx_vec Mz(nValues);
    arma::cx_vec s(etl);
    arma::cx_vec out;

    for (unsigned int ii = 0; ii < nValues ; ii++)
    {
	M(2,ii) = 1;
    }

    arma::vec z = arma::linspace<arma::vec>(-sliceThickness*2.0,
					    sliceThickness*2.0,
					    nValues);

    double areaCrusher = 2/(gamma2pi*sliceThickness);
    // double areaRefocusingExcite = -(BWExcite*pulseWidthExcite)/(gamma2pi*sliceThickness);
    double areaRefocusingHide = -(BWHide*pulseWidthHide)/(gamma2pi*sliceThickness);
    
    arma::mat RFExcite = RotationArbitrary(M_PI/2,initialPhase);
    arma::mat RFHide = RotationX(M_PI/2);

    arma::mat Relax(3,3);
    arma::vec Recover(3);

    Relaxation(Relax,Recover,esp/2,T1,T2);

    arma::mat RecoverMat = arma::repmat(Recover,1,nValues);
    
    // Excite
    // SLRRotation(M, 100*z, rfExcite, pulseTypeExcite, 
    // 		initialPhase, pulseWidthExcite, TBWExcite, 
    // 		100*sliceThickness);
    // HardPrecession(M, z, 2*M_PI*areaRefocusingExcite*0.53));
	
    SLRRotation(M, 100*z, rfHide, pulseTypeHide, 
    		initialPhase, pulseWidthHide, TBWHide,
    		100*sliceThickness);
    HardPrecession(M, z, 2*M_PI*areaRefocusingHide*0.53);

    if (etl == 0)
    {
	for (unsigned int zIndex = 0; zIndex < nValues; zIndex++)
	{
	    Mxy(zIndex) = M(0,zIndex) + j*M(1,zIndex);
	    Mz(zIndex) = M(2,zIndex);
	}
	
	switch (returnType)
	{
	case MXY:
	    out = Mxy;
	    printf("\treturning: Mxy\n");
	    break;
	case MZ:
	    out = Mz;
	    printf("\treturning: Mz\n");
	    break;
	}

	return out;
    }
    
    // recover
    M = Relax*M + RecoverMat;

    // Start the echo train
    for (unsigned int echoNumber = 0; echoNumber < etl; echoNumber++)
    {
	arma::mat RFRefocus = RotationX(flips(echoNumber)/180.0*M_PI);
	double TBWRefocus = 1.5353;
	double pulseWidthRefocus = 1.288e-3; // s
	unsigned int pulseTypeRefocus = 1;
	arma::vec rfRefocus = GetMSINC(nPointsRF,TBWRefocus,flips(echoNumber));

	// Refocus
	HardPrecession(M, z, 2*M_PI*areaCrusher);
	SLRRotation(M, 100*z, rfRefocus, pulseTypeRefocus, 
		    0.0, pulseWidthRefocus, TBWRefocus, 
		    100*sliceThickness);
	HardPrecession(M, z, 2*M_PI*areaCrusher);

	
	M = Relax*M + RecoverMat;	

	for (unsigned int zIndex = 0; zIndex < nValues; zIndex++)
    	{
    	    Mxy(zIndex) = M(0,zIndex) + j*M(1,zIndex);
    	    Mz(zIndex) = M(2,zIndex);
	}

	s(echoNumber) = arma::mean(Mxy);
	M = Relax*M + RecoverMat;

    }

    printf("\tsimulation complete\n");

    switch (returnType)
    {
    case ECHOTRAIN:
	out = s;
	printf("\treturning: echotrain\n");
	break;
    case MXY:
	out = Mxy;
	printf("\treturning: Mxy\n");
	break;
    case MZ:
	out = Mz;
	printf("\treturning: Mz\n");
	break;
    }

    return out;
			  
}

#include <complex>
#include <cmath>
#include <string>

#include <armadillo>

#include "bloch_utilities.hpp"
#include "leroux.hpp"
#include "rf_pulses.hpp"

#include "echotrain_lrx.hpp"

arma::cx_vec BlochFSELRX(unsigned int nValues, unsigned int etl,  unsigned int pulseType,
			 double initialPhase, unsigned int returnType, double T1, double T2,
			 double esp)
{

    printf("SIMULATION STARTING: leroux\n");
	
    std::complex<double> j(0,1);

    double gamma2pi = 42.58; // kHz/mT
    double sliceThickness = 0.005; // m

    arma::mat M = arma::zeros<arma::mat>(3,nValues);
    arma::cx_vec Mxy(nValues);
    arma::cx_vec Mz(nValues);
    arma::cx_vec s(etl);
    arma::cx_vec out;
    
    double sliceScaling = 0.4;
    
    for (unsigned int ii = 0; ii < nValues ; ii++)
    {
	M(2,ii) = 1;
    }

    arma::vec z = arma::linspace<arma::vec>(-sliceThickness*4.0,
					    sliceThickness*4.0,
					    nValues);

    arma::vec transmitPhases(etl);
    arma::vec receivePhases(etl);
    CreatePhases(receivePhases,transmitPhases,etl);
    
    double areaCrusher = 2 / (gamma2pi*sliceThickness); // mT*ms/m
    
    double TBWExcite, TBWRefocus, pulseWidthExcite, pulseWidthRefocus, BWExcite, flipAngleMSINC;
    unsigned int pulseTypeExcite, pulseTypeRefocus, nPointsRF;
    arma::vec rfExcite;
    arma::vec rfRefocus;
    
    switch (pulseType)
    {
    case 1:
	TBWExcite = 2.56;
	pulseWidthExcite = 1.288e-3; // s
	BWExcite = TBWExcite/pulseWidthExcite;
	nPointsRF = 320;
	flipAngleMSINC = 90.0;
	rfExcite = GetMSINC(nPointsRF,TBWExcite,flipAngleMSINC);
	pulseTypeExcite = 1;

	TBWRefocus = 1.5353;
    	pulseWidthRefocus = 1.288e-3; // s
	pulseTypeRefocus = 1;
	flipAngleMSINC = 160.0;
    	rfRefocus = GetMSINC(nPointsRF,TBWRefocus,flipAngleMSINC);
	sliceScaling = 1/sliceScaling;
	break;

    case 2:

	TBWExcite = 3.55;
	pulseWidthExcite = 1.200e-3; // s
	BWExcite = TBWExcite/pulseWidthExcite;
	pulseTypeExcite = 1;
	rfExcite = GetSLRExcite();

	TBWRefocus = 3.55;
    	pulseWidthRefocus = 1.064e-3; // s
    	pulseTypeRefocus = 2;
    	rfRefocus = GetDIVERSE();
	sliceScaling = 1/sliceScaling;
	break;
        
    case 3:

	TBWExcite = 2.56;
	pulseWidthExcite = 1.288e-3; // s
	BWExcite = TBWExcite/pulseWidthExcite;
	nPointsRF = 320;
	flipAngleMSINC = 90.0;
	rfExcite = GetMSINC(nPointsRF,TBWExcite,flipAngleMSINC);
	pulseTypeExcite = 1;

	TBWRefocus = 1.55;
    	pulseWidthRefocus = 1.1e-3; // s
    	pulseTypeRefocus = 1;
    	rfRefocus = GetFourier();
	sliceScaling = 1/sliceScaling;
	break;

    case 4:

	TBWExcite = 3.55;
	pulseWidthExcite = 1.200e-3; // s
	BWExcite = TBWExcite/pulseWidthExcite;
	pulseTypeExcite = 1;
	rfExcite = GetSLRExcite();

	TBWRefocus = 1.55;
    	pulseWidthRefocus = 1.1e-3; // s
    	pulseTypeRefocus = 1;
    	rfRefocus = GetSLRLow();
	sliceScaling = 1/sliceScaling;
	break;

    case 5:

	TBWExcite = 3.55;
	pulseWidthExcite = 1.200e-3; // s
	BWExcite = TBWExcite/pulseWidthExcite;
	pulseTypeExcite = 1;
	rfExcite = GetSLRExcite();

	TBWRefocus = 3.55;
    	pulseWidthRefocus = 2.8-3; // s
    	pulseTypeRefocus = 1;
    	rfRefocus = GetSLRHigh();
	sliceScaling = 1/sliceScaling;

	esp += 0.0018; // penalize the ESP for the longer RF pulse
	break;

	
    }

    double TBWSpinEcho = 3.55;
    double pulseWidthSpinEcho = 3.500e-3; // s
    arma::vec rfSpinEcho = GetSLRSE();
    unsigned int pulseTypeSpinEcho = 1;
    
    arma::mat Relax(3,3);
    arma::vec Recover(3);
    double areaRefocusingExcite = -(BWExcite*pulseWidthExcite)/(gamma2pi*sliceThickness);
    

    // Excitation
    SLRRotation(M, 100*z, rfExcite, pulseTypeExcite, 
		initialPhase, pulseWidthExcite, TBWExcite, 
		100*sliceThickness);
    HardPrecession(M, z, 2*M_PI*areaRefocusingExcite*0.53);

    double spinEchoTE = (0.04-0.02) + 2*0.02; // s
    Relaxation(Relax,Recover,spinEchoTE/2,T1,T2);
    arma::mat RecoverMat = arma::repmat(Recover,1,nValues);

    // Relaxation
    M = Relax*M + RecoverMat;

    // Refocusing
    HardPrecession(M, z, 2*M_PI*areaCrusher);
    SLRRotation(M, 100*z, rfSpinEcho, pulseTypeSpinEcho, 
		0.25*2*M_PI, pulseWidthSpinEcho, TBWSpinEcho,
		100*sliceThickness);
    HardPrecession(M, z, 2*M_PI*areaCrusher);

    // Relaxation
    M = Relax*M + RecoverMat;

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
 
    // start FSE stuff
    Relaxation(Relax,Recover,esp/2,T1,T2);
    RecoverMat = arma::repmat(Recover,1,nValues);
    M = Relax*M + RecoverMat;
    
    for (unsigned int echoNumber = 0; echoNumber < etl; echoNumber++)
    {

	M = Relax*M + RecoverMat;

	HardPrecession(M, z, 2*M_PI*areaCrusher);
	SLRRotation(M, 100*z, rfRefocus, pulseTypeRefocus, 
		    transmitPhases(echoNumber), pulseWidthRefocus, TBWRefocus, 
		    100*sliceThickness*sliceScaling);
	HardPrecession(M, z, 2*M_PI*areaCrusher);

	for (unsigned int zIndex = 0; zIndex < nValues; zIndex++)
	{
	    Mxy(zIndex) = M(0,zIndex) + j*M(1,zIndex);
	    Mz(zIndex) = M(2,zIndex);
	}

	s(echoNumber) = arma::mean(Mxy)*(cos(receivePhases(echoNumber)) -
					 j*sin(receivePhases(echoNumber)));	

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

	Relaxation(Relax,Recover,0.3,T1,T2);
	RecoverMat = arma::repmat(Recover,1,nValues);
	
	M = Relax*M + RecoverMat;

	for (unsigned int zIndex = 0; zIndex < nValues; zIndex++)
	{
	    Mz(zIndex) = M(2,zIndex);
	}

	out = Mz;
	printf("\treturning: Mz\n");
	break;
    }

    return out;
}


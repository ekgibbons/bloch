#include <complex>
#include <cmath>
#include <string>

#include <armadillo>

#include "bloch_utilities.hpp"
#include "rf_pulses.hpp"

#include "echotrain_alsop.hpp"

arma::cx_vec BlochFSEAlsop(const arma::vec &flips, unsigned int nValues, unsigned int etl,
			   double initialPhase, unsigned int returnType, double T1,
			   double T2, double esp)
{

    printf("SIMULATION STARTING: alsop\n");

    std::complex<double> j(0,1);
    
    double gamma2pi = 42.58; // kHz/mT
    double sliceThickness = 0.005; // m

    double TBWExcite = 3.55;
    double pulseWidthExcite = 1.200e-3; // s
    double BWExcite = TBWExcite/pulseWidthExcite;
    arma::vec rfExcite = GetSLRExcite();
    unsigned int pulseTypeExcite = 1;

    double TBWSpinEcho = 3.55;
    double pulseWidthSpinEcho = 3.500e-3; // s
    arma::vec rfSpinEcho = GetSLRSE();
    unsigned int pulseTypeSpinEcho = 1;
    
    double TBWHide = 1.5;
    double pulseWidthHide = 1.288e-3; // s
    double BWHide = TBWHide/pulseWidthHide;
    unsigned int nPointsRF = 320;
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

    arma::vec z = arma::linspace<arma::vec>(-sliceThickness*3.0,
					    sliceThickness*3.0,
					    nValues);

    // need 8 cycles of dephase and about 2 cycles of crushing to make this stable.
    // i think.  
    double areaDephase = 4/(gamma2pi*sliceThickness);
    double areaCrusher = 2/(gamma2pi*sliceThickness) + areaDephase;
    double areaRefocusingExcite =	-(BWExcite*pulseWidthExcite)/(gamma2pi*sliceThickness/0.8);
    double areaRefocusingHide = -(BWHide*pulseWidthHide)/(gamma2pi*sliceThickness);

    arma::mat RFExcite = RotationArbitrary(M_PI/2,initialPhase);
    arma::mat RFHide = RotationX(M_PI/2);

    arma::mat Relax(3,3);
    arma::vec Recover(3);

    
    // Excite
    // SLRRotation(M, 100*z, rfExcite, pulseTypeExcite, 
    // 		initialPhase, pulseWidthExcite, TBWExcite, 
    // 		100*sliceThickness);
    // HardPrecession(M, z, 2*M_PI*areaRefocusingExcite*0.53);

    SLRRotation(M, 100*z, rfHide, pulseTypeHide, 
    		initialPhase, pulseWidthHide, TBWHide,
    		100*sliceThickness);
    HardPrecession(M, z, 2*M_PI*areaRefocusingHide*0.56);

    // M = RFExcite*M;
    
    // Diffusion
    // double D = 0.002; // mm^2/s
    // double T = 0.02; // gradient duration (s)
    // double dT = 0.00005; // grad step
    double spinEchoTE = (0.04-0.02) + 2*0.02; // s
    // double Tbetween = (spinEchoTE - 2*T)/2;
    Relaxation(Relax,Recover,spinEchoTE/2,T1,T2);
    arma::mat RecoverMat = arma::repmat(Recover,1,nValues);
    
    // arma::vec G((unsigned int)(T/dT));
    // G = 50*G.ones(); // gradient mT/m
    
    // // rotation from the gradient+diffusion
    // Diffusion(M, z, G, 0, dT);

    // // Move during the deadtime
    // double standardDeviationBetween = sqrt(2*D*Tbetween);
    // arma::vec zAdd(nValues);
    // z = z + standardDeviationBetween*zAdd.randn();

    // Relaxation
    M = Relax*M + RecoverMat;


    
    // // arma::mat RF180 = RotationX(M_PI);
    // // M = RF180*M;
    // might need to bust out the slice widening trick again...
    HardPrecession(M, z, 2*M_PI*areaCrusher);
    SLRRotation(M, 100*z, rfSpinEcho, pulseTypeSpinEcho, 
    		0.25*2*M_PI, pulseWidthSpinEcho, TBWSpinEcho,
    		100*sliceThickness/0.5);
    HardPrecession(M, z, 2*M_PI*areaCrusher);

    
    // // diffusion again
    // z = z + standardDeviationBetween*zAdd.randn();
    // Diffusion(M, z, G, D, dT);
    
    // Relaxation
    M = Relax*M + RecoverMat;

 
    
    // Dephase and Hide
    // HardPrecession(M, z, 2*M_PI*areaRefocusingExcite*0.53);
    // SLRRotation(M, 100*z, rfExcite, pulseTypeExcite, 
    // 		0.0, pulseWidthExcite, TBWExcite, 
    // 		100*sliceThickness);
    // HardPrecession(M, z, 2*M_PI*areaRefocusingExcite*0.53);
    
    // M = RFHide*M;

    

    HardPrecession(M, z, 2*M_PI*(areaRefocusingExcite*0.56-areaDephase));
    SLRRotation(M, 100*z, rfExcite, pulseTypeExcite, 
    		0.25*2*M_PI, pulseWidthExcite, TBWExcite, 
    		100*sliceThickness/0.8);
    HardPrecession(M, z, 2*M_PI*(areaRefocusingExcite*0.56));





    
    // HardPrecession(M, z, 2*M_PI*(areaRefocusingHide*0.58-areaDephase));
    // SLRRotation(M, 100*z, rfHide, pulseTypeHide, 
    // 		0.25*2*M_PI, pulseWidthHide, TBWHide,
    // 		100*sliceThickness);
    // HardPrecession(M, z, 2*M_PI*(areaRefocusingHide*0.58));    
    // HardPrecession(M, z, 2*M_PI*(areaRefocusingHide*0.53-areaDephase));

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
    Relaxation(Relax,Recover,esp/2,T1,T2);
    RecoverMat = arma::repmat(Recover,1,nValues);
    M = Relax*M + RecoverMat;

    // Start the echo train
    for (unsigned int echoNumber = 0; echoNumber < etl; echoNumber++)
    {
	// printf("echo: %i, status: computing\n", echoNumber);
	arma::mat RFRefocus = RotationX(flips(echoNumber)/180.0*M_PI);

	double TBWRefocus = 1.5353;
	double pulseWidthRefocus = 1.288e-3; // s
	unsigned int pulseTypeRefocus = 1;
	arma::vec rfRefocus = GetMSINC(nPointsRF,TBWRefocus,flips(echoNumber));

	
	if (echoNumber == 0)
	{
	    HardPrecession(M, z, 2*M_PI*(areaCrusher+areaRefocusingExcite*0.56));
	    // HardPrecession(M, z, 2*M_PI*(areaCrusher+areaRefocusingHide*0.53));
	    // HardPrecession(M, z, 2*M_PI*(areaCrusher));
	}
	else
	{
	    HardPrecession(M, z, 2*M_PI*(areaCrusher-areaDephase));
	}

	// Refocus
	SLRRotation(M, 100*z, rfRefocus, pulseTypeRefocus, 
		    0.25*2*M_PI, pulseWidthRefocus, TBWRefocus, 
		    100*sliceThickness);

	// M = RFRefocus*M;
	
	HardPrecession(M, z, 2*M_PI*(areaCrusher+areaDephase));

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

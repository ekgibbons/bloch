#include <complex>
#include <cmath>
#include <string>

#include <armadillo>

#include "bloch_utilities.hpp"
#include "rf_pulses.hpp"

#include "cardiac_sat.hpp"

arma::cx_vec BlochCardiacSat(const double sliceThickness, const unsigned int nValues,
			     const unsigned int nExcitations, const double T1,
			     const double T2, const unsigned int returnType)
{

    printf("SIMULATION STARTING: cardiac saturation\n");

    std::complex<double> j(0,1);
    
    double gamma2pi = 42.58; // kHz/mT

    unsigned int nPointsRF = 100;
    
    double initialPhase = M_PI/2; // phase in rads

    double flipAngleMSINC = 12.0;
    double TBWExcite = 1.6;
    double pulseWidthExcite = 0.400e-3; // s
    double BWExcite = TBWExcite/pulseWidthExcite;
    unsigned int pulseTypeExcite = 1; // for msinc, use 1
    arma::vec rfExcite = GetMSINC(nPointsRF,TBWExcite,flipAngleMSINC);

    double areaRefocusingExcite = -(BWExcite*pulseWidthExcite)/(gamma2pi*sliceThickness);

    arma::mat spoiler = Spoiling();
    
    arma::mat Relax(3,3);
    arma::vec Recover(3);

    Relaxation(Relax,Recover,0.0025/2,T1,T2);

    arma::mat RecoverMat = arma::repmat(Recover,1,nValues);

    arma::mat RelaxSat(3,3);
    arma::vec RecoverSat(3);

    Relaxation(RelaxSat,RecoverSat,0.02,T1,T2);

    arma::mat RecoverSatMat = arma::repmat(RecoverSat,1,nValues);

    
    arma::vec z = arma::linspace<arma::vec>(-sliceThickness*2.0,
    					    sliceThickness*2.0,
    					    nValues);
    
    arma::mat M = arma::zeros<arma::mat>(3,nValues);
    arma::cx_vec Mxy(nValues);
    arma::cx_vec Mz(nValues);
    arma::cx_vec s(nExcitations);

    for (unsigned int ii = 0; ii < nValues ; ii++)
    {
    	M(2,ii) = 1;
    }    

    arma::mat RFExcite = RotationArbitrary(M_PI/2,initialPhase);
    M = RFExcite*M;
    M = spoiler*M;

    M = RelaxSat*M + RecoverSatMat;	
    
    for (unsigned int excitationIndex = 0;
	 excitationIndex < nExcitations;
	 excitationIndex++)
    {
	printf("\texcitation number: %i\n",excitationIndex);
	
	// Excite
	SLRRotation(M, 100*z, rfExcite, pulseTypeExcite, 
		    initialPhase, pulseWidthExcite, TBWExcite, 
		    100*sliceThickness);
	HardPrecession(M, z, 2*M_PI*areaRefocusingExcite*0.53);

    	for (unsigned int zIndex = 0; zIndex < nValues; zIndex++)
    	{
    	    Mxy(zIndex) = M(0,zIndex) + j*M(1,zIndex);
    	    Mz(zIndex) = M(2,zIndex);
    	}
	
	s(excitationIndex) = arma::mean(Mxy);

	M = spoiler*M;
	
	M = Relax*M + RecoverMat;	
    } 
    

    printf("simulation complete\n");
    
    // output
    arma::cx_vec out;

    switch (returnType)
    {
    case ECHOTRAIN:
    	out = s;
    	printf("returning: echotrain\n");
    	break;
    case MXY:
    	out = Mxy;
    	printf("returning: Mxy\n");
    	break;
    case MZ:
    	out = Mz;
    	printf("returning: Mz\n");
    	break;
    }

    return out;
			  
}


    // double areaCrusher = 2/(gamma2pi*sliceThickness);
    // // double areaRefocusingExcite = -(BWExcite*pulseWidthExcite)/(gamma2pi*sliceThickness);
    // double areaRefocusingHide = -(BWHide*pulseWidthHide)/(gamma2pi*sliceThickness);
    
    // arma::mat RFExcite = RotationArbitrary(M_PI/2,initialPhase);
    // arma::mat RFHide = RotationX(M_PI/2);


    // arma::mat RecoverMat = arma::repmat(Recover,1,nValues);
    
    // // Excite
    // // SLRRotation(M, 100*z, rfExcite, pulseTypeExcite, 
    // // 		initialPhase, pulseWidthExcite, TBWExcite, 
    // // 		100*sliceThickness);
    // // HardPrecession(M, z, 2*M_PI*areaRefocusingExcite*0.53));
	
    // SLRRotation(M, 100*z, rfHide, pulseTypeHide, 
    // 		initialPhase, pulseWidthHide, TBWHide,
    // 		100*sliceThickness);
    // HardPrecession(M, z, 2*M_PI*areaRefocusingHide*0.53);

    // // recover
    // M = Relax*M + RecoverMat;

    // // Start the echo train
    // for (unsigned int echoNumber = 0; echoNumber < etl; echoNumber++)
    // {
    // 	arma::mat RFRefocus = RotationX(flips(echoNumber)/180.0*M_PI);
    // 	double TBWRefocus = 1.5353;
    // 	double pulseWidthRefocus = 1.288e-3; // s
    // 	unsigned int pulseTypeRefocus = 1;
    // 	arma::vec rfRefocus = GetMSINC(nPointsRF,TBWRefocus,flips(echoNumber));

    // 	// Refocus
    // 	HardPrecession(M, z, 2*M_PI*areaCrusher);
    // 	SLRRotation(M, 100*z, rfRefocus, pulseTypeRefocus, 
    // 		    0.0, pulseWidthRefocus, TBWRefocus, 
    // 		    100*sliceThickness);
    // 	HardPrecession(M, z, 2*M_PI*areaCrusher);

	
    // 	M = Relax*M + RecoverMat;	

    	// for (unsigned int zIndex = 0; zIndex < nValues; zIndex++)
    	// {
    	//     Mxy(zIndex) = M(0,zIndex) + j*M(1,zIndex);
    	//     Mz(zIndex) = M(2,zIndex);
    	// }

    // 	s(echoNumber) = arma::mean(Mxy);
    // 	M = Relax*M + RecoverMat;

    // }

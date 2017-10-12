#include <complex>
#include <cmath>
#include <string>

#include <armadillo>

#include "bloch_utilities.hpp"
#include "rf_pulses.hpp"

#include "cardiac_sat.hpp"

double FloatMod(double firstTerm, double secondTerm)
{
    if (firstTerm >= secondTerm)
    {
	firstTerm -= secondTerm;
	firstTerm = FloatMod(firstTerm,secondTerm);
    }

    return firstTerm;
}

double CardiacMovement(double t, double heartRate)
{

    double beatLength, lengthSystole, tUse;
    double zSystoleEnd, z, amplitude;

    
    double decayCoefficient = 0.06;
    double recoveryCoefficient = 0.12;

   
    heartRate /= 60; // bps
    beatLength = 1/(heartRate+1e-15);; // s
    lengthSystole = 0.3; // s
    amplitude = 0.01; // m

    tUse = FloatMod(t, beatLength);
    
    z = 0;
    if (t < lengthSystole)
    {
	z = amplitude*exp(-t/decayCoefficient);
    }
    if (t >= lengthSystole)
    {
	zSystoleEnd = amplitude*exp(-lengthSystole/decayCoefficient);
	z = amplitude - (amplitude - zSystoleEnd)*
	    exp(-(t - lengthSystole)/recoveryCoefficient);
    }
    
    // printf("z: %f, tOrig: %f, tUse: %f\n",z,t,tUse);

    if (heartRate == 0)
    {
	return amplitude;
    }

    
    return z;
}

arma::cx_vec BlochCardiacSat(const double sliceThickness,
			     const unsigned int nValues,
			     const unsigned int nExcitations, const double T1,
			     const double T2, double heartRate,
			     const unsigned int returnType)
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

    double runningTime = 0.0; // s
    double deltaTimeSat = 0.02; // s
    double deltaTimeTR = 0.0025; // s
    
    arma::mat Relax(3,3);
    arma::vec Recover(3);

    Relaxation(Relax,Recover,0.0025,T1,T2);

    arma::mat RecoverMat = arma::repmat(Recover,1,nValues);
    
    arma::mat RelaxSat(3,3);
    arma::vec RecoverSat(3);

    Relaxation(RelaxSat,RecoverSat,0.02,T1,T2);

    arma::mat RecoverSatMat = arma::repmat(RecoverSat,1,nValues);

    
    arma::vec z = arma::linspace<arma::vec>(-sliceThickness*4.0,
    					    sliceThickness*4.0,
    					    nValues);
    
    arma::mat M = arma::zeros<arma::mat>(3,nValues);
    arma::cx_vec Mxy(nValues);
    arma::cx_vec Mz(nValues);
    arma::cx_vec s(nExcitations);

    for (unsigned int ii = 0; ii < nValues ; ii++)
    {
    	M(2,ii) = 1;
    }    

    // movement setup
    double deltaZ = 0.0; // cm
    double zShift = 0.0; // cm
    double zShiftOld = 0.0; // cm
    
    arma::mat RFExcite = RotationArbitrary(M_PI/2,initialPhase);
    M = RFExcite*M;
    M = spoiler*M;

    M = RelaxSat*M + RecoverSatMat;
    runningTime += deltaTimeSat; 

    zShift = CardiacMovement(runningTime,heartRate);

    z += deltaZ;
    zShiftOld = zShift;

    for (unsigned int excitationIndex = 0;
	 excitationIndex < nExcitations;
	 excitationIndex++)
    {
	if (excitationIndex % 20 == 0)
	{
	    printf("\texcitation number: %i\n",excitationIndex);
	}
	
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

	runningTime += deltaTimeTR;
	zShift = CardiacMovement(runningTime,heartRate);

	deltaZ = zShift - zShiftOld;
	z += deltaZ;
	zShiftOld = zShift;
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

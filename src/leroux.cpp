#include "leroux.hpp"

void CreatePhases(arma::vec &ROut, arma::vec &XOut, unsigned int etl)
{

    double D, d, r;

    unsigned int echoNumber;

    double R[etl + 1];
    double X[etl + 1];
    
    double dX;
    double dR[] = {-0.2376, -0.2515, -0.6971, -0.2034, 0.4780, -0.4694, -1.0170, -0.3567, 0.2410, 0.4447, 
    		    0.8831, 0.4056, 0.3200, -0.1339, -0.1842, -0.5396, -0.7111, -0.5064, -0.5032, -0.3439,  
    		    -0.3182, -0.0461, 0.1491, 0.0974, 0.3490, 0.3166, 0.4383, 0.1516, 0.1695, 0.0960, 
    		    -0.0437, -0.1279, -0.1294, -0.0418, -0.0780, -0.1405, 0.0213, -0.0675, -0.0285, -0.2021, 
    		    -0.0396, -0.1423, -0.1819, -0.1700, -0.0920, -0.0318, -0.0782, 0.0477, 0.1699, 0.1377, 
    		    0.1885, 0.1458, 0.2836, 0.0535, 0.0849, -0.0242, 0.0145, -0.1608, -0.1692, -0.1521, 
    		    -0.1479, -0.2158, -0.1574, -0.1012, -0.0503, -0.1184, 0.0042, 0.0270, 0.0595, 0.0, 0.0};

    unsigned int ncor = (unsigned int)(sizeof(dR)/sizeof(dR[0]));

    D = 161.0/843.0;
    d = 0.5*D + 0.25;
    r = 0.0;


    
    for (echoNumber = 0; echoNumber < etl; echoNumber++)
    {
	r = 2*d + r;
	r = r - floor(r);
	R[echoNumber] = r;
	X[echoNumber] = R[echoNumber] - d;
	d = d + D;
	d = d - floor(d);

    }

    if (ncor > etl)
    {
	ncor = etl;
    }
    
    for (echoNumber = 0; echoNumber < (ncor - 1); echoNumber++)
    {
	dX = 0.5*0.1*(dR[echoNumber] + dR[echoNumber + 1]);
	X[echoNumber] += dX;
	X[echoNumber] -= floor(X[echoNumber]);
	R[echoNumber] += 0.1*dR[echoNumber + 1];
	R[echoNumber] -= floor(R[echoNumber]);

    }

    for (echoNumber = 0; echoNumber < etl; echoNumber++)
    {
	X[echoNumber] -= 0.1*dR[0];
	X[echoNumber] -= floor(X[echoNumber]);
	
	if (X[echoNumber] > 0.5)
	{
	    X[echoNumber] -= 1;
	}

	XOut(echoNumber) = 2*M_PI*X[echoNumber];
	
	R[echoNumber] -= 0.1*dR[0];
	R[echoNumber] -= floor(R[echoNumber]);


	
	if (R[echoNumber] > 0.5)
	{
	    R[echoNumber] -= 1;
	}
	ROut(echoNumber) = 2*M_PI*R[echoNumber];

    }


}

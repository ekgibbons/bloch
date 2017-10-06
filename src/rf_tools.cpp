#include <complex>

#include <armadillo>

#include "rf_tools.hpp"

void abrx(const arma::vec &rf, const arma::vec &gx, const arma::vec &x, 
	  double phi, arma::cx_vec &ALF, arma::cx_vec &BET) 
{
    std::complex<double> j(0,1);
    double alf[2], bet[2];
    
    unsigned int nx = (unsigned int)x.n_elem;

    for (unsigned int ix = 0; ix < nx; ix++) 
    {
	alf[0] = 1.0; 
	alf[1] = 0.0; 
	bet[0] = 0.0; 
	bet[1] = 0.0;
	double xLocation = x(ix);
	abrot(alf, bet, xLocation, phi, rf, gx);
	ALF(ix) = alf[0] + j*alf[1];
	BET(ix) = bet[0] + j*bet[1];
    }
}

void abrot(double (&a)[2], double (&b)[2], double x, double phiOriginal, 
	   const arma::vec &rf, 
	   const arma::vec &gx)
{
    
    double nx, ny, nz, snp, csp, cg, cpr, cpi;
    double al[2], be[2], ap[2], bp[2];
    unsigned int ns = (unsigned int)rf.n_elem;

    for (unsigned int k = 0; k < ns; k++)
    {
	cg = x*gx(k);
	cpr = cos(phiOriginal)*rf(k);
	cpi = sin(phiOriginal)*rf(k);

	double phi = sqrt(cg*cg+cpr*cpr+cpi*cpi);

	if (phi > 0.0)
	{
	    nx = cpr/phi; 
	    ny = cpi/phi; 
	    nz = cg/phi;
	} 
	else {
	    nx = 0.0; 
	    ny = 0.0; 
	    nz = 1.0;   /* doesn't matter, phi=0*/
	}

	csp = cos(phi/2); 
	snp = sin(phi/2);
	al[0] = csp; 
	al[1] = nz*snp;
	be[0] = ny*snp; 
	be[1] = nx*snp;
	
	bp[0] = al[0]*b[0]-al[1]*b[1]+be[0]*a[0]-be[1]*(-a[1]);
	bp[1] = al[0]*b[1]+al[1]*b[0]+be[1]*a[0]+be[0]*(-a[1]);
	    
	ap[0] =   -(  be[0] *b[0]-(-be[1])*b[1]) 
	    +   al[0] *a[0]-(-al[1])*(-a[1]);
	ap[1] = -(-(-(be[1])*b[0]+  be[0] *b[1]) 
		  + (-al[1])*a[0]+  al[0] *(-a[1]));
	
	a[0] = ap[0]; 
	a[1] = ap[1]; 
	b[0] = bp[0]; 
	b[1] = bp[1];
    }
    
    return;
}

arma::mat abfm2m(const arma::cx_vec &a, const arma::cx_vec &b, arma::vec &m)
{
    std::complex<double> j(0,1);
    
    std::complex<double> mxy = m(0) + j*m(1);
    arma::cx_vec mck = {{mxy, 
			 std::conj(mxy), 
			 m(2)}};
    
    arma::mat mout(3,a.n_elem);

    for (unsigned int k = 0; k < a.n_elem; k++)
    {
	std::complex<double> ak = a(k);
	std::complex<double> bk = b(k);
	std::complex<double> akConj = std::conj(ak);
	std::complex<double> bkConj = std::conj(bk);

	arma::cx_mat Rck(2,3);
	Rck(0,0) = pow(akConj,2);
	Rck(0,1) = -pow(bk,2);
	Rck(0,2) = (std::complex<double>(2,0))*akConj*bk;
	Rck(1,0) = -akConj*bkConj;
	Rck(1,1) = -ak*bk;
	Rck(1,2) = ak*akConj - bk*bkConj;

	arma::cx_vec mtemp = Rck*mck; 

	mout(0,k) = mtemp(0).real();
	mout(1,k) = mtemp(0).imag();
	mout(2,k) = mtemp(1).real();
    }
    return mout;
}

void abfm2m(std::complex<double> a, std::complex<double> b, arma::mat &m)
{
    std::complex<double> j(0,1);
    
    std::complex<double> mxy = m(0) + j*m(1);
    arma::cx_vec mck = {{mxy, 
			 std::conj(mxy), 
			 m(2)}};
    
    std::complex<double> aConj = std::conj(a);
    std::complex<double> bConj = std::conj(b);

    arma::cx_mat Rck(2,3);
    Rck(0,0) = pow(aConj,2);
    Rck(0,1) = -pow(b,2);
    Rck(0,2) = (std::complex<double>(2,0))*aConj*b;
    Rck(1,0) = -aConj*bConj;
    Rck(1,1) = -a*b;
    Rck(1,2) = a*aConj - b*bConj;

    arma::cx_vec mtemp = Rck*mck; 


    m(0) = mtemp(0).real();
    m(1) = mtemp(0).imag();
    m(2) = mtemp(1).real();

}

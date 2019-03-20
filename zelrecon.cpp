#include	<cstdlib>
#include	<fstream>
#include	<sstream>
#include	<iostream>
#include	<iomanip>
#include	<cmath>
#include	<vector>
#include	<string>
#include	<exception>


// These can be modified later if we convert to MPI.
void	myexit(const int flag) {
  exit(flag);
}
void	myexception(const std::exception& e) {
  std::cout<<"Exception: "<<e.what()<<std::endl;
  myexit(1);
}






class	GaussLegendre {
// Class to hold the abscissae and weights for G-L integration.
// Form:  \int_{-1}^{+1} dx f(x) = \sum w_j f(x_j)
public:
  std::vector<double>	x,w;
  int			N;
  void set(const int NN) {
    // Set the weights and abscissae, which are public variables.
    N = NN;
    try {
      x.resize(N);
      w.resize(N);
    } catch(std::exception& e) {myexception(e);}
    // First find the abscissae...this involves finding the roots
    // of the Legendre polynomials.  We use Newton-Raphson starting
    // from an analytic guess.
    int N2 = (N+1)/2;	// Weights/roots are symmetric, do half ...
    for (int i=0; i<N2; ++i) {
      // Find root by starting "close" and using N-R.
      double zz,dpdz,z=cos(M_PI*(i+0.75)/(N+0.5));
      do {
        double p1=1,p2=0;
        for (int j=0; j<N; ++j) {
          double p3 = p2;
          p2 = p1;
          p1 = ((2*j+1.)*z*p2-j*p3)/(j+1.);
        }
        // Now p1 is P_n and p2 is P_{n-1}.  Compute dP/dz
        dpdz = N*(z*p1-p2)/(z*z-1);
        // and use N-R to update:
        zz = z;
        z  = zz-p1/dpdz;
      } while(fabs(z-zz)>1e-15);
      x[  i  ] = -z;   w[  i  ] = 2.0/(1-z*z)/dpdz/dpdz;
      x[N-i-1] =  z;   w[N-i-1] = 2.0/(1-z*z)/dpdz/dpdz;
    }
  }
  GaussLegendre() {}
  GaussLegendre(const int N) {set(N);}
};






class	Spline {
private:
  std::vector<double> xa,ya,y2;
  int	nn;
public:
  Spline(const std::vector<double>& x, const std::vector<double>& y) {
    double p,qn,sig,un;
    if (x.size() != y.size()) {
      std::cout << "x and y must have the same dimensions." << std::endl;
      myexit(1);
    }
    nn = x.size();
    if (x[0]>=x[nn-1]) {
      std::cout << "x must be monotonically increasing in spline." << std::endl;
      myexit(1);
    }
    xa = x;
    ya = y;
    try {y2.resize(nn);} catch(std::exception& e) {myexception(e);}
    std::vector<double> u(nn);
    y2[0]=u[0]=0.0;
    for (int i=1;i<nn-1;++i) {
      sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
      p=sig*y2[i-1]+2.0;
      y2[i]=(sig-1.0)/p;
      u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
      u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    }
    qn=un=0.0;
    y2[nn-1]=(un-qn*u[nn-2])/(qn*y2[nn-2]+1.0);
    for (int k=nn-2;k>=0;--k)
      y2[k]=y2[k]*y2[k+1]+u[k];
  }
  double operator() (const double x) {
    int	k,klo,khi;
    if (x<xa[0] || x>xa[nn-1]) {
      std::cout << "x out of range [" << xa[0] << ","
                << xa[nn-1] << "] in Spline." << std::endl;
      myexit(1);
    }
    klo=0;
    khi=nn-1;
    while (khi-klo > 1) {	// Bisection search.
      k=(khi+klo) >> 1;
      if (xa[k] > x) khi=k; else klo=k;
    }
    double h=xa[khi]-xa[klo];
    if (h<=0.0) {std::cout << "h==0 in bisection." << std::endl; myexit(1); }
    double a=(xa[khi]-x)/h;
    double b=(x-xa[klo])/h;
    return(a*ya[klo]+b*ya[khi]+
      ((a*a*a-a)*y2[klo]+(b*b*b-b)*y2[khi])*(h*h)/6.0);
  }
  double xmin() {
    return(xa[0]);
  }
  double xmax() {
    return(xa[nn-1]);
  }
};












class	Zeldovich {
// Computes the correlation function within the Zel'dovich approximation.
// Based upon Carlson et al. (2013; MNRAS, 429, 1674): CLPT.
// This version has reconstruction capabilities.
protected:
  // Table lengths for storing functions -- should be even.
  const static int	NqTable=2048,NkTable=2048;
  GaussLegendre		gl;
  std::vector<double>	kLin,pLin,etaPer,etaPar,uVal,xiLin;
  std::vector<double>	J2Lin,J3Lin,J4Lin;
  double		qmin,qmax,Rf,delta,dkinv,sigma2;
  std::vector<double> sphBess(const double x) {
    // Returns j0(x) and j1(x)/x.
    std::vector<double> jl(2);
    if (fabs(x)<1e-3) {
      double x2=x*x;
      jl[0] = 1.0 + x2*(-1/6.  +x2/120.);
      jl[1] = 1/3.+ x2*(-1./30.+x2/840.);
    }
    else {
      jl[0] = sin(x)/x;
      jl[1] = (jl[0]-cos(x))/x/x;
    }
    return(jl);
  }
  void readPowerSpectrum(const char fname[]) {
    // Load the linear power spectrum from a file, and expand it out
    // if necessary.  This is stored log-log.
    kLin.clear();  pLin.clear();
    std::ifstream fs(fname);
    if (!fs) {
      std::cout<<"Unable to open "<<fname<<" for reading."<<std::endl;
      myexit(1);
    }
    std::string buf;
    do {	// Skip any preamble comment lines.
      getline(fs,buf);
    } while(!fs.eof() && buf[0]=='#');
    while (!fs.eof()) {
      double kval,pval;
      std::istringstream ss(buf);
      ss >> kval >> pval;
      try {
        kLin.push_back(log(kval));
        pLin.push_back(log(pval));
      } catch(std::exception& e) {myexception(e);}
      getline(fs,buf);
    }
    fs.close();
    // Now resample to an "even" log-k spacing.
    Spline ss(kLin,pLin);
    kLin.clear(); pLin.clear();
    for (int i=0; i<NkTable; ++i) {
      double x = ss.xmin() + (i+0.5)*(ss.xmax()-ss.xmin())/NkTable;
      try {
        kLin.push_back(x);
        pLin.push_back(ss(x));
      } catch(std::exception& e) {myexception(e);}
    }
    // and set up dkinv so we can interpolate quickly.
    dkinv = NkTable/(kLin[NkTable-1]-kLin[0]);
  }
  double calcSigma2(const int itype) {
    // Computes sigma^2_eta, the (square of the) 1D dispersion in the
    // displacement field.  Eq. (29) of CLPT paper.  Does the integral
    // in ln(k), assuming kLin is equally spaced in ln(k) and that there
    // are enough points in the array for such an integral.
    // Fixed dumb mistake in zero-lag for cross-term!
    int wt=4; double sum=0;
    for (int i=1; i<kLin.size(); ++i) {
      double sk1,sk2;
      double kk = exp(kLin[i]);
      switch (itype) {
        case 0: sk1=sk2=1.0;                     break;
        case 1: sk1=sk2=1.0-exp(-kk*kk*Rf*Rf/2); break;
        case 2: sk1=sk2=   -exp(-kk*kk*Rf*Rf/2); break;
        case 3: sk1=    1.0-exp(-kk*kk*Rf*Rf/2);
                sk2=       -exp(-kk*kk*Rf*Rf/2); break;
        default:
          std::cerr<<"Unknown itype="<<itype<<" in calcSigma2."<<std::endl;
          myexit(1);
      }
      if (itype==3)
        sum += exp(kLin[i]+pLin[i])*0.5*(sk1*sk1+sk2*sk2)*wt;
      else
        sum += exp(kLin[i]+pLin[i])*sk1*sk2*wt;
      wt   = 8/wt;
    }
    sum *= (kLin[2]-kLin[0])/6;
    sum /= 6*M_PI*M_PI;
    return(sum);
  }
  std::vector<double> calcQfuncs(const double q, const int itype) {
    // Returns etaPerp and etaPar, Eqs. (30-31) of CLPT paper.
    // These, up to a factor of f, are the \Psi of Reid&White, Eqs. (9-10).
    // Also, since it is useful, returns U(q) of Eq. (32) as qf[2]
    // and the linear xi as qf[3].
    int wt=4; std::vector<double> sum(4);
    int Nk=kLin.size();
    int Nint=(int)(8*exp(kLin[Nk-1])*q+1024);
    if (Nint>=20000) Nint=20000;
    const double xmax=100*M_PI;
    double lkmax=(log(xmax/(q+0.01))>kLin[Nk-1])?kLin[Nk-1]:log(xmax/(q+0.01));
    double hh=(lkmax-kLin[0])/Nint;
    for (int i=1; i<Nint; ++i) {
      double xx = kLin[0]+i*hh;
      double kk = exp(xx);
      double ap = cos(M_PI/2.*exp(xx-lkmax));
      double sk1,sk2;
      switch (itype) {
        case 0: sk1=sk2=1.0;                     break;
        case 1: sk1=sk2=1.0-exp(-kk*kk*Rf*Rf/2); break;
        case 2: sk1=sk2=   -exp(-kk*kk*Rf*Rf/2); break;
        case 3: sk1=    1.0-exp(-kk*kk*Rf*Rf/2);
                sk2=       -exp(-kk*kk*Rf*Rf/2); break;
        default:
          std::cerr<<"Unknown itype="<<itype<<" in calcQfuncs."<<std::endl;
          myexit(1);
      }
      int    jj = (int)(i*hh*dkinv);
      if (jj>=pLin.size()-2) jj=pLin.size()-2;
      double pk = exp(pLin[jj]+(xx-kLin[jj])*
                     (pLin[jj+1]-pLin[jj])/(kLin[jj+1]-kLin[jj]));
      std::vector<double> jl=sphBess(kk*q);
      sum[0] += kk*pk*sk1*sk2*(        jl[1])*wt;	// eta_per, Eq. (30)
      sum[1] += kk*pk*sk1*sk2*(jl[0]-2*jl[1])*wt;	// eta_par, Eq. (31)
      sum[2] +=-kk*pk*  1*sk2*(kk*kk*q*jl[1])*wt;	// U,       Eq. (32)
      sum[3] += kk*pk* ap*  1*(kk*kk*  jl[0])*wt;	// xi_lin
      wt      = 8/wt;
    }
    sum[0] *= hh/3.0/(2*M_PI*M_PI);
    sum[1] *= hh/3.0/(2*M_PI*M_PI);
    sum[2] *= hh/3.0/(2*M_PI*M_PI);
    sum[3] *= hh/3.0/(2*M_PI*M_PI);
    return(sum);
  }
  std::vector<double> calc_Jn(const double q) {
    // Computes the \mathcal{J}_n integrals, which are used in the shear terms.
    const int Nk=kLin.size();
    const int Nint=25000;
    const double xmax=100*M_PI;
    double lkmax=(log(xmax/(q+0.01))>kLin[Nk-1])?kLin[Nk-1]:log(xmax/(q+0.01));
    const double hh=(lkmax-kLin[0])/Nint;
    double sum1=0,sum2=0,sum3=0,sum4=0;
    for (int i=1; i<Nint; ++i) {
      double xx = kLin[0]+i*hh;
      double ap = cos(M_PI/2.*exp(xx-lkmax));
      double kk = exp(xx);
      double k2 = kk*kk;
      double kq = kk*q;
      int    jj = (int)(i*hh*dkinv);
      if (jj>=pLin.size()-2) jj=pLin.size()-2;
      double pk = exp(pLin[jj]+(xx-kLin[jj])*
                     (pLin[jj+1]-pLin[jj])/(kLin[jj+1]-kLin[jj]));
      std::vector<double> jl=sphBess(kq);
      double j0,j1,j2,j3,j4;
      j0=jl[0];
      if (kq<0.9) {
        double kq2 = kq*kq;
        j1 = kq *(1./3.+ kq2*(-1./30.+kq2*(1./840.-kq2/45360.)));
        j2 = kq2*(1./15.+kq2*(-1./210.+kq2*(1./7560-kq2/498960.)));
        j3 = kq *kq2*(1./105.+kq2*(-1./1890+kq2*(1./83160-kq2/6486480.)));
        j4 = kq2*kq2*(1./945.+kq2*(-1./20790.+kq2/1081080.));
      }
      else {
        j1 =    jl[1]*kq;
        j2 = 3.*jl[1]  -jl[0];
        j3 = 5.*j2/(kq)-j1;
        j4 = 7.*j3/(kq)-j2;
      }
      int wt= 2+2*(i%2);
      sum1 += k2*pk*kk*(j2)*wt;
      sum2 += k2*pk*(2./15.*j1-1./5.*j3)*wt * ap;
      sum3 += k2*pk*(-1./5.*j1-1./5.*j3)*wt;
      sum4 += k2*pk*(j3)*wt;
    }
    std::vector<double> sum(5);
    sum[1] = sum1 * hh/3.0/(2*M_PI*M_PI);       // mathcal{J}_1
    sum[2] = sum2 * hh/3.0/(2*M_PI*M_PI);       // mathcal{J}_2
    sum[3] = sum3 * hh/3.0/(2*M_PI*M_PI);       // mathcal{J}_3
    sum[4] = sum4 * hh/3.0/(2*M_PI*M_PI);       // mathcal{J}_4
    return(sum);
  }
  void tabulateQfuncs(const int itype=0) {
    // Finally, tabulate sigma2, etaPer, etaPar, etc.
    qmin = 1.0/exp(kLin[NkTable-1]);  if (qmin<0.2) qmin=0.2;
    qmax = 1.0/exp(kLin[    0    ]);  if (qmax>250) qmax=250;
    // First compute it on a coarse grid.
    std::vector<double> qvals;
    const int Nsample=128;
    try {
       qvals.resize(Nsample);
      etaPer.resize(Nsample);
      etaPar.resize(Nsample);
        uVal.resize(Nsample);
       xiLin.resize(Nsample);
       J2Lin.resize(Nsample);
       J3Lin.resize(Nsample);
       J4Lin.resize(Nsample);
    } catch(std::exception& e) {myexception(e);}
    delta=(qmax-qmin)/(Nsample-1);
    for (int i=0; i<Nsample; ++i) {
      double qq = qmin+i*delta;
      std::vector<double> qf=calcQfuncs(qq,itype);
      std::vector<double> Jn=calc_Jn(   qq);
       qvals[i] = qq;
      etaPer[i] = qf[0];
      etaPar[i] = qf[1];
      uVal[  i] = qf[2];
      xiLin[ i] = qf[3] * qq*qq;
      J2Lin[ i] = Jn[2];
      J3Lin[ i] = Jn[3];
      J4Lin[ i] = Jn[4];
    }
    // then fit splines and retabulate it onto a finer grid.
    Spline etaPerSpline(qvals,etaPer);
    Spline etaParSpline(qvals,etaPar);
    Spline   uValSpline(qvals,uVal);
    Spline  xiLinSpline(qvals,xiLin);
    Spline  J2LinSpline(qvals,J2Lin);
    Spline  J3LinSpline(qvals,J3Lin);
    Spline  J4LinSpline(qvals,J4Lin);
    try {
      etaPer.resize(NqTable);
      etaPar.resize(NqTable);
        uVal.resize(NqTable);
       xiLin.resize(NqTable);
       J2Lin.resize(NqTable);
       J3Lin.resize(NqTable);
       J4Lin.resize(NqTable);
    } catch(std::exception& e) {myexception(e);}
    sigma2 = calcSigma2(itype);
    delta=(qmax-qmin)/(NqTable-1);
    for (int i=0; i<NqTable; ++i) {
      double qq = qmin+i*delta;
      etaPer[i] = etaPerSpline(qq);
      etaPar[i] = etaParSpline(qq);
      uVal[  i] =   uValSpline(qq);
      xiLin[ i] =  xiLinSpline(qq)/qq/qq;
      J2Lin[ i] =  J2LinSpline(qq);
      J3Lin[ i] =  J3LinSpline(qq);
      J4Lin[ i] =  J4LinSpline(qq);
    }
  }
  std::vector<double> interpQfuncs(const double q) {
    // Does a linear interpolation to return etaPer and etaPar.
    // Also returns U(q) and xi_lin.
    std::vector<double> qf(9);
    int k=(NqTable-1)*(q-qmin)/(qmax-qmin);
    if (q>qmin && q<qmax) {
      double dq = (q-(qmin+k*delta))/delta;
      qf[0]=etaPer[k]+dq*(etaPer[k+1]-etaPer[k]);
      qf[1]=etaPar[k]+dq*(etaPar[k+1]-etaPar[k]);
      qf[2]=  uVal[k]+dq*(  uVal[k+1]-  uVal[k]);
      qf[3]= xiLin[k]+dq*( xiLin[k+1]- xiLin[k]);
      qf[4]= 0;
      qf[5]= 0;
      qf[6]= J2Lin[k]+dq*( J2Lin[k+1]- J2Lin[k]);
      qf[7]= J3Lin[k]+dq*( J3Lin[k+1]- J3Lin[k]);
      qf[8]= J4Lin[k]+dq*( J4Lin[k+1]- J4Lin[k]);
    }
    else {
      const double TINY=1e-10;
      if (q<qmin) {
        qf[0]=sigma2 - TINY;
        qf[1]=sigma2 - TINY;
        qf[2]=qf[3]=qf[4]=qf[5]=qf[6]=qf[7]=qf[8]=0;
      }
      if (q>qmax) {
        qf[0]=etaPer[NqTable-1];
        qf[1]=etaPar[NqTable-1];
        qf[2]=uVal[  NqTable-1];
        qf[3]=xiLin[ NqTable-1];
        qf[4]=qf[5]=0;
        qf[6]=J2Lin[ NqTable-1];
        qf[7]=J3Lin[ NqTable-1];
        qf[8]=J4Lin[ NqTable-1];
      }
    }
    return(qf);
  }
  std::vector<double> calcAmat(const double q[]) {
    // Returns the 3x3 matrix A (Eq. 28 of CLPT).
    double qhat[3],qq=0;
    for (int i=0; i<3; ++i) qq += q[i]*q[i];  qq=sqrt(qq);
    for (int i=0; i<3; ++i) qhat[i] = q[i]/qq;
    std::vector<double> qf=interpQfuncs(qq);
    double F =2*(sigma2-qf[0]);	// sigma_perp^2
    double G =2*(qf[0] -qf[1]);	// sigma_par ^2 - sigma_perp^2
    std::vector<double> Amat(9);
    for (int i=0; i<3; ++i)
      for (int j=i; j<3; ++j) {
        Amat[3*i+j] = Amat[3*j+i]  = G*qhat[i]*qhat[j];
        if (i==j)     Amat[3*i+i] += F;
      }
    return(Amat);
  }
  std::vector<double> calcAinv(const double q[]) {
    // Returns the inverse of the 3x3 matrix A (Eq. 28 of CLPT).
    // Also returns the determinant of Ainv as the last (extra) element.
    // The Sherman-Morrison formula says that
    //   (A+b.cT)^{-1}=Ainv - Ainv.b.cT.Ainv/[1+cT.Ainv.b]
    // so our (F.1+G.q.q)^{-1}=1/F-G.q.q/F/[F+G]
    // For real-space (r-q).Ainv.(r-q) depends only on r, q and rq.mu.
    // Moving into redshift space simply requires us to
    // divide the zz element by (1+f)^2 [and change det], but now we
    // need to do the phi integral numerically as well.
    double qhat[3],qq=0;
    for (int i=0; i<3; ++i) qq += q[i]*q[i];  qq=sqrt(qq);
    for (int i=0; i<3; ++i) qhat[i] = q[i]/qq;
    std::vector<double> qf=interpQfuncs(qq);
    double  F =2*(sigma2-qf[0]);	// sigma_perp^2
    double  G =2*(qf[0] -qf[1]);	// sigma_par ^2 - sigma_perp^2
    double FpG=2*(sigma2-qf[1]);	// sigma_par ^2
    std::vector<double> Ainv(10);
    for (int i=0; i<3; ++i)
      for (int j=i; j<3; ++j) {
        Ainv[3*i+j] = Ainv[3*j+i]  = -G*qhat[i]*qhat[j]/F/FpG;
        if (i==j)     Ainv[3*i+i] += 1.0/F;
      }
    // Now set detA.  Use det(cM)=c^n det(M) and det(I+u.vT)=1+uT.v so that
    // det(F.delta_ij+G.qhat_i.qhat_j) = F^3(1+G/F) = F^2(F+G).
    // Also note that F+G is 2[sigma^2-etaPar] for our case, which is the
    // same thing as sigma_{||}, while F is \sigma_\perp, thus
    // detA=(sigma_perp^2.sigma_par)^2
    Ainv[9] = 1.0/(F*F*FpG);
    return(Ainv);
  }
  double zeldovichIntegrand(const double r[], const double q[],
                            const double f1, const double f2){
    // Returns the integrand for the 1+\xi integral (Eq. 34 of CLPT).
    const double twoPi3=248.05021344239853;
    std::vector<double> Ainv=calcAinv(q);	// Also returns detAinv.
    // If we are in redshift space, need to correct for A->RAR.
    if (f1>0 || f2>0) {
      for (int i=0; i<3; ++i) {
        Ainv[3*i+2] /= (1+f2);
        Ainv[3*2+i] /= (1+f1);
      }
      Ainv[9] /= (1+f1)*(1+f2);
    }
    double res=0;
    for (int i=0; i<3; ++i)
      for (int j=0; j<3; ++j)
        res += (r[i]-q[i])*Ainv[3*i+j]*(r[j]-q[j]);
    res = exp(-0.5*res)*sqrt(Ainv[9]/twoPi3);
    return(res);
  }
public:
  void init(const char fname[], const double Rfilter, const int itype=0) {
    Rf = Rfilter;
    // Initialize the G-L integration points and weights.
    gl.set(128);
    // Load the linear power spectrum from a file, and expand it out
    // if necessary.
    readPowerSpectrum(fname);
    // Finally pre-tabulate the q-dependent functions.
    tabulateQfuncs(itype);
  } // End init.
  Zeldovich() {}
  Zeldovich(const char fname[], const double Rfilter, const int itype=0) {
    init(fname,Rfilter,itype);
  }
  ~Zeldovich() {}
  void print_eta() {
    // A convenience/debug feature to print sigma^2, etaPerp and etaPar.
    // Also prints useful combinations of these, and prints U(q).
    std::cout<<"sigma^2="<<std::scientific<<sigma2<<std::endl;
    std::cout<<"# "<<std::setw(10)<<"q(Mpc/h)"
                   <<std::setw(12)<<"Eta_perp"
                   <<std::setw(12)<<"Eta_par"
                   <<std::setw(12)<<"Sig2_perp"
                   <<std::setw(12)<<"Sig2_par"
                   <<std::setw(12)<<"U"
                   << std::endl;
    for (int i=0; i<NqTable; ++i) {
      double qq     = qmin + i*delta;
      double sigper = 2*(sigma2-etaPer[i]);
      double sigpar = 2*(sigma2-etaPar[i]);
      std::cout
        <<std::scientific<<std::setw(12)<<std::setprecision(4)<<qq
        <<std::scientific<<std::setw(12)<<std::setprecision(4)<<etaPer[i]
        <<std::scientific<<std::setw(12)<<std::setprecision(4)<<etaPar[i]
        <<std::scientific<<std::setw(12)<<std::setprecision(4)<<sigper
        <<std::scientific<<std::setw(12)<<std::setprecision(4)<<sigpar
        <<std::scientific<<std::setw(12)<<std::setprecision(4)<<uVal[i]
        <<std::endl;
    }
    std::cout.flush();
  }
  double xiL(const double r) {
    // The real-space, linear correlation function at r.
    // This is not tested for very large or small values of r.
    std::vector<double> qf=interpQfuncs(r);
    return(qf[3]);
  }
  double xiZ(const double rval) {
    // Returns the real-space, matter, Zel'dovich correlation function at r.
    // This is not tested for very large or small values of r.
    // The integration is over x=q-r, in length and angle with the
    // azimuthal integral being trivial.
    const double xmin=0;
    const double xmax=10*sqrt(sigma2);
    const double rr[3]={0,0,rval};
    const double r2   =rval*rval;
    const int    Nx=2000;
    const double dx=(xmax-xmin)/Nx;
    double xi=0;
    for (int ixx=0; ixx<Nx; ++ixx) {
      double xx=xmin+(ixx+0.5)*dx;
      double x2=xx*xx;
      for (int imu=0; imu<gl.N; ++imu) {
        double mu = gl.x[imu];
        // Compute vec{q}=vec{r}+vec{x} with vec{r}=r.zhat,
        // so r_z=r, x_z=x*mu, cos_rq=(r_z+x_z)/q.
        double qlen = sqrt(r2+x2+2*rval*xx*mu);
        double qcos = (rval+xx*mu)/qlen;
        double qq[3] ={qlen*sqrt(1-qcos*qcos),0,qlen*qcos};
        if (qlen>qmin && qlen<qmax)
          xi += x2 * zeldovichIntegrand(rr,qq,0,0) * gl.w[imu];
      }
    }
    xi *= dx;		// Convert sum to integral.
    xi *= 2*M_PI;	// The azimuthal integral.
    xi -= 1.0;		// Calculated 1+xi, subtract 1 for xi.
    return(xi);
  }
  std::vector<double> xiContributions(const double rval) {
    // Returns the different contributions to the real-space Zel'dovich
    // correlation function for locally biased tracers.
    // This is not tested for very large or small values of r.
    // The integration is over x=q-r, in length and angle with the
    // azimuthal integral being trivial.
    const double xmin=0;
    const double xmax=10*sqrt(sigma2);
    const double rr[3]={0,0,rval};
    const double r2   =rval*rval;
    const int    Nx=500;
    const double dx=(xmax-xmin)/Nx;
    std::vector<double> xi(7);
    for (int ixx=0; ixx<Nx; ++ixx) {
      double xx=xmin+(ixx+0.5)*dx;
      double x2=xx*xx;
      for (int imu=0; imu<gl.N; ++imu) {
        double mu = gl.x[imu];
        // Compute vec{q}=vec{r}+vec{x} with vec{r}=r.zhat,
        // so r_z=r, x_z=x*mu, cos_rq=(r_z+x_z)/q.
        double qlen = sqrt(r2+x2+2*rval*xx*mu);
        double qcos = (rval+xx*mu)/qlen;
        double qsin = sqrt(1-qcos*qcos);
        double qq[3] ={qlen*qsin,0,qlen*qcos};
        double qh[3] ={     qsin,0,     qcos};
        if (qlen>qmin && qlen<qmax) {
          // For the unbiased tracers we only need this--all other terms
          // are multiplied by this anyway.
          double pref = x2 * zeldovichIntegrand(rr,qq,0,0) * gl.w[imu];
          // For the bias terms, compute U,xi and Ainv (even though in above).
          std::vector<double> qf  =interpQfuncs(qlen);
          std::vector<double> Ainv=calcAinv(qq);
          // Construct the auxilliary matrix/vectors g, G of CLPT Eq. (45)
          double g[3],G[9];
          for (int i=0; i<3; ++i) {
            g[i]=0;
            for (int j=0; j<3; ++j)
              g[i] += Ainv[3*i+j]*(qq[j]-rr[j]);
          }
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              G[3*i+j]=Ainv[3*i+j]-g[i]*g[j];
          double Ug,UUG;  Ug=UUG=0;
          for (int i=0; i<3; ++i) {
            Ug += (qf[2]*qh[i])*g[i];
            for (int j=0; j<3; ++j)
              UUG += qf[2]*qf[2]*qh[i]*qh[j]*G[3*i+j];
          }
          // The <s^2 Delta Delta> term:
          double shear=0;
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j) {
              double upsilon= qh[i]*qh[j]*(3*qf[6]*qf[6]+4*qf[6]*qf[7]+
                              2*qf[6]*qf[8]+2*qf[7]*qf[7]+4*qf[7]*qf[8]+
                              qf[8]*qf[8]) + (i==j)*2*qf[7]*qf[7];
              shear += G[3*i+j]*upsilon;
            }
          shear *= 2;
          // Now do the 1, Fp, Fpp, Fp^2, Fp.Fpp, Fpp^2 terms.
          xi[0] +=    pref;
          xi[1] += -2*pref*Ug;
          xi[2] +=   -pref*UUG;
          xi[3] +=    pref*(qf[3]-UUG);
          xi[4] += -2*pref*qf[3]*Ug;
          xi[5] +=0.5*pref*qf[3]*qf[3];
          xi[6] +=   -pref*shear;
        }
      }
    }
    for (int j=0; j<xi.size(); ++j) {
      xi[j] *= dx;	// Convert sum to integral.
      xi[j] *= 2*M_PI;	// The azimuthal integral.
    }
    xi[0] -= 1.0;	// Calculated 1+xi, subtract 1 for xi.
    return(xi);
  }
  std::vector<double> xiContributions(const double rval, const double mu,
                                      const double f1, const double f2) {
    // Returns the different contributions to the redshift-space Zel'dovich
    // correlation function for locally biased tracers.
    // This is not tested for very large or small values of r.
    // The integration is over x=q-r, in length and angle with the
    // azimuthal integral being done explicitly.
    const double xmin=0;
    const double xmax=10*sqrt(sigma2);
    const double rr[3]={rval*sqrt(1-mu*mu),0,rval*mu};
    const double r2   =rval*rval;
    const int    Nx=256,Nphi=32;
    const double dx=(xmax-xmin)/Nx;
    const double dphi=2*M_PI/Nphi;
    std::vector<double> xi(7);
    for (int ixx=0; ixx<Nx; ++ixx) {
      double xx=xmin+(ixx+0.5)*dx;
      double x2=xx*xx;
      for (int imu=0; imu<gl.N; ++imu) {
        double mu = gl.x[imu];
        double st = sqrt(1-mu*mu);
        for (int iphi=0; iphi<Nphi; ++iphi) {
          double phi = (iphi+0.5)*dphi;
          double qq[3],qh[3],xv[3];
          xv[0]=xx*st*cos(phi);xv[1]=xx*st*sin(phi);xv[2]=xx*mu;
          double qlen=0;
          for (int i=0; i<3; ++i) {
            qq[i] = xv[i]+rr[i];
            qlen += qq[i]*qq[i];
          }
          qlen = sqrt(qlen);
          for (int i=0; i<3; ++i) qh[i] = qq[i]/qlen;
          if (qlen>qmin && qlen<qmax) {
            // For the unbiased tracers we only need this--all other terms
            // are multiplied by this anyway.
            double pref = x2 * zeldovichIntegrand(rr,qq,f1,f2) * gl.w[imu];
            // Compute U, xi and Ainv (even though in above).  Here we
            // need to multiply A by R_ij=(1+fzz).  We correct U below.
            std::vector<double> qf  =interpQfuncs(qlen);
            std::vector<double> Ainv=calcAinv(qq);
            for (int i=0; i<3; ++i) {
              Ainv[3*i+2] /= (1+f2);
              Ainv[3*2+i] /= (1+f1);
            }
            Ainv[9] /= (1+f1)*(1+f2);
            // Construct the auxilliary matrix/vectors g, G of CLPT Eq. (45)
            double g[3],U[3],G[9];
            for (int i=0; i<3; ++i) {
              g[i]=0;
              for (int j=0; j<3; ++j)
                g[i] += Ainv[3*i+j]*(qq[j]-rr[j]);
              U[i]=qf[2]*qh[i];
            }
            U[2] *= 1+f2;	// Correct U as U->RU.
            for (int i=0; i<3; ++i)
              for (int j=0; j<3; ++j)
                G[3*i+j]=Ainv[3*i+j]-g[i]*g[j];
            double Ug,UUG;  Ug=UUG=0;
            for (int i=0; i<3; ++i) {
              Ug += U[i]*g[i];
              for (int j=0; j<3; ++j)
                UUG += U[i]*U[j]*G[3*i+j];
            }
            // The <s^2 Delta Delta> term:
            double shear=0;
            for (int i=0; i<3; ++i)
              for (int j=0; j<3; ++j) {
                double upsilon= qh[i]*qh[j]*(3*qf[6]*qf[6]+4*qf[6]*qf[7]+
                                2*qf[6]*qf[8]+2*qf[7]*qf[7]+4*qf[7]*qf[8]+
                                qf[8]*qf[8]) + (i==j)*2*qf[7]*qf[7];
                shear += G[3*i+j]*upsilon;
              }
            shear *= 2;
            // Now do the 1, Fp, Fpp, Fp^2, Fp.Fpp, Fpp^2 & shear terms.
            xi[0] +=    pref;
            xi[1] += -2*pref*Ug;
            xi[2] +=   -pref*UUG;
            xi[3] +=    pref*(qf[3]-UUG);
            xi[4] += -2*pref*qf[3]*Ug;
            xi[5] +=0.5*pref*qf[3]*qf[3];
            xi[6] +=   -pref*shear;
          }
        }
      }
    }
    for (int j=0; j<xi.size(); ++j) {
      xi[j] *= dx*dphi;	// Convert sum to integral.
    }
    xi[0] -= 1.0;	// Calculated 1+xi, subtract 1 for xi.
    return(xi);
  }
  std::vector<double> xiContributions(const double rval,
                                      const double f1, const double f2) {
    // Returns the contributions to the multipoles of the redshift-space
    // correlation function for locally biased tracers.
    // This is not tested for very large or small values of r.
    const int Nmu=4;
    GaussLegendre gg = GaussLegendre(2*Nmu);	// Must be even.
    // For even lengths, can sum over half of the points.
    std::vector<double> xiell;
    try{xiell.resize(14);}catch(std::exception& e) {myexception(e);}
    for (int i=0; i<Nmu; ++i) {
      std::vector<double> ximu = xiContributions(rval,gg.x[i],f1,f2);
      double p0=1.0;
      double p2=0.5*(3*gg.x[i]*gg.x[i]-1);
      for (int j=0; j<ximu.size(); ++j) {
        xiell[0*ximu.size()+j] += ximu[j]*gg.w[i] * p0 * 1;
        xiell[1*ximu.size()+j] += ximu[j]*gg.w[i] * p2 * 5;
      }
    }
    return(xiell);
  }
};






int	main(int argc, char **argv)
{
  if (argc!=7) {
    std::cout<<"Usage: zeldovich <pkfile> <f> <b1> <b2> <bs> <R>"
             <<std::endl;
    exit(1);
  }
  const double ff = atof(argv[2]);
  const double Rf = atof(argv[6]);

  // Set up the values of (b,f) for the different "types" of
  // correlation function (0=NoRecon, 1=DD, 2=SS, 3=DS) so we
  // can just refer to them in the loops below.
  const int Ntype=4;
  double b1[Ntype],b2[Ntype],b1b1[Ntype],b1b2[Ntype],b2b2[Ntype];
  double bs[Ntype];
  double f1[Ntype],f2[Ntype];
  b1[0] = atof(argv[3]); b1[1]=  b1[0];   b1[2]=0;   b1[3]=0.5*b1[0];
  b2[0] = atof(argv[4]); b2[1]=  b2[0];   b2[2]=0;   b2[3]=0.5*b2[0];
  bs[0] = atof(argv[5]); bs[1]=  bs[0];   bs[2]=0;   bs[3]=0.5*bs[0];
  b1b1[0]=b1[0]*b1[0]; b1b1[1]=b1b1[0]; b1b1[2]=0; b1b1[3]=0;
  b1b2[0]=b1[0]*b2[0]; b1b2[1]=b1b2[0]; b1b2[2]=0; b1b2[3]=0;
  b2b2[0]=b2[0]*b2[0]; b2b2[1]=b2b2[0]; b2b2[2]=0; b2b2[3]=0;
#ifdef	SHIFTRANDOM
  f1[0]=ff;  f1[1]=ff; f1[2]=ff; f1[3]=ff;
  f2[0]=ff;  f2[1]=ff; f2[2]=ff; f2[3]=ff;
#else
  f1[0]=ff;  f1[1]=ff; f1[2]= 0; f1[3]=ff;
  f2[0]=ff;  f2[1]=ff; f2[2]= 0; f2[3]= 0;
#endif
  // Create instances for each type of correlation function.
  std::vector<Zeldovich> zel(Ntype);
  for (int itype=0; itype<Ntype; ++itype)
    zel[itype].init(argv[1],Rf,itype);

  const int Nr=40;
  const double rmin=50.0,rmax=125;
#ifdef	SIMPLE
  for (int i=0; i<Nr; ++i) {
    double rr = rmin + i*(rmax-rmin)/(Nr-1);
    std::cout<<std::fixed<<std::setw(10)<<std::setprecision(2)<<rr
             <<std::fixed<<std::setw(15)<<std::setprecision(5)
             <<zel[0].xiL(rr)*rr*rr
             <<std::fixed<<std::setw(15)<<std::setprecision(5)
             <<zel[0].xiZ(rr)*rr*rr
             <<std::endl;
  }
  return(0);
#endif

#ifdef	PIECES
  for (int i=0; i<Nr; ++i) {
    double rr = rmin + i*(rmax-rmin)/(Nr-1);
    std::vector<double> xir=zel[0].xiContributions(rr,ff,ff);
    double xitot= xir[0]+b1[0]*xir[1]+b2[0]*xir[2]
                + b1b1[0]*xir[3]+b1b2[0]*xir[4]+b2b2[0]*xir[5]+bs[0]*xir[6];
    std::cout<<std::fixed<<std::setw(10)<<std::setprecision(2)<<rr;
    for (int j=0; j<xir.size(); ++j)
      std::cout<<std::fixed<<std::setw(15)<<std::setprecision(8)
               <<xir[j]/xitot;
    std::cout<<std::endl;
  }
  return(0);
#endif

  for (int i=0; i<Nr; ++i) {
    double rr = rmin + i*(rmax-rmin)/(Nr-1);
    std::cout<<std::fixed<<std::setw(8)<<std::setprecision(2)<<rr
             <<std::fixed<<std::setw(9)<<std::setprecision(4)
             <<(1+b1[0])*(1+b1[0])*zel[0].xiL(rr)*rr*rr;
    for (int it=0; it<Ntype; ++it) {
      double xi;
      std::vector<double> xir=zel[it].xiContributions(rr);
      xi=xir[0]+b1[it]*xir[1]+b2[it]*xir[2]
        +b1b1[it]*xir[3]+b1b2[it]*xir[4]+b2b2[it]*xir[5]
        +bs[it]*xir[6];
      std::cout<<std::fixed<<std::setw(9)<<std::setprecision(4)<<xi*rr*rr;
    }
    for (int it=0; it<Ntype; ++it) {
      double xi0,xi2;
      std::vector<double> xis=zel[it].xiContributions(rr,f1[it],f2[it]);
      xi0=xis[0]+b1[it]*xis[1]+b2[it]*xis[2]
         +b1b1[it]*xis[3]+b1b2[it]*xis[4]+b2b2[it]*xis[5]
         +bs[it]*xis[6];
      xi2=xis[7]+b1[it]*xis[8]+b2[it]*xis[9]
         +b1b1[it]*xis[10]+b1b2[it]*xis[11]+b2b2[it]*xis[12]
         +bs[it]*xis[13];
      std::cout<<std::fixed<<std::setw(9)<<std::setprecision(4)<<xi0*rr*rr;
      std::cout<<std::fixed<<std::setw(9)<<std::setprecision(4)<<xi2*rr*rr;
    }
    std::cout<<std::endl;
  }

  return(0);
}

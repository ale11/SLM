#include "TurbModel_WVM.hpp"

TurbModel_WVM::TurbModel_WVM(const char *name)
{
  nshells = 32;
  nmodes = 2000;
  npar = nshells*nmodes;

  Cn = 2.2;
  Cv = 1.0;
  Ceps1 = 1.5;
  Ceps2 = 11.0/6.0;
  Ceps3 = 0.01;

  e = new double[npar][3];
  a = new double[npar][3];
  u = new double[npar][3];
}

TurbModel_WVM::~TurbModel_WVM()
{
  delete[] e;
  delete[] a;
  delete[] u;
}

void TurbModel_WVM::initialHookScalarRansTurbModel(double *struc,  
						   double &eps_init)
{
  srand(time(NULL));
  for (int i = 0; i < 6; i++)
  {
    rey[i] = 0.0;
    dim[i] = 0.0;
  }

  tke = 0.25;

  // Energy spectrum
  double mag_k[nshells], dk[nshells], Enukdk[nshells];
  double fkt = pow(4.0,0.2), Cnuk = 0.4843;

  mag_k[0] = 0.01;
  for (int i = 0; i < nshells; i++)
    mag_k[i] = pow(fkt,i)*mag_k[0];

  for (int i = 1; i < nshells; i++)
    dk[i] = mag_k[i] - mag_k[i-1];
  dk[0] = mag_k[0];
  
  double Enuk_a = 0.0, Enuk_b;
  for (int i = 0; i < nshells; i++)
  {
    Enuk_b = 2.0*tke*Cnuk*pow(mag_k[i],4.0)/
             pow(1.0 + mag_k[i]*mag_k[i],17.0/6.0);
    
    Enukdk[i] = (Enuk_a + Enuk_b)/2.0*dk[i];
    
    Enuk_a = Enuk_b;
  }

  // Initialize particle properties
  for (int ishells = 0; ishells < nshells; ishells++)
    for (int imodes = 0; imodes < nmodes; imodes++)
    {
      int ipar = ishells*nmodes + imodes;
      
      double theta;
      double r[3], s[3];

      // wave vector
      double u1 = (double)rand()/RAND_MAX;
      double u2 = (double)rand()/RAND_MAX;
      theta = acos(2.0*u2 - 1.0);
      double phi = 2.0*M_PI*u1;
      
      e[ipar][0] = sin(theta)*cos(phi);
      e[ipar][1] = sin(theta)*sin(phi);
      e[ipar][2] = cos(theta);

      // eddy-axis vector
      double eye[3] = {1.0, 0.0, 0.0};
      vecCrossVec3d(r, eye, e[ipar]);
      vecCrossVec3d(s, r, e[ipar]);
      normVec3d(r);
      normVec3d(s);
      theta = ( (double)rand()/RAND_MAX )*2.0*M_PI;

      for (int i = 0; i < 3; i++)
	a[ipar][i] = cos(theta)*r[i] + sin(theta)*s[i];

      // velocity vector
      theta = ( (double)rand()/RAND_MAX )*2.0*M_PI;

      double umag2 = 2.0*nshells*Enukdk[ishells];
      for (int i = 0; i < 3; i++)
	u[ipar][i] = sqrt(umag2)*( cos(theta)*r[i] + sin(theta)*s[i] );

      // initial Reynolds stresses and dimensionality
      rey[0] += u[ipar][0]*u[ipar][0];
      rey[1] += u[ipar][1]*u[ipar][1];
      rey[2] += u[ipar][2]*u[ipar][2];
      rey[3] += u[ipar][0]*u[ipar][1];
      rey[4] += u[ipar][0]*u[ipar][2];
      rey[5] += u[ipar][1]*u[ipar][2];

      dim[0] += umag2*e[ipar][0]*e[ipar][0];
      dim[1] += umag2*e[ipar][1]*e[ipar][1];
      dim[2] += umag2*e[ipar][2]*e[ipar][2];
      dim[3] += umag2*e[ipar][0]*e[ipar][1];
      dim[4] += umag2*e[ipar][0]*e[ipar][2];
      dim[5] += umag2*e[ipar][1]*e[ipar][2];
    }

  // structure tensors, TKE, dissipation
  for (int i = 0; i < 6; i++)
  {
    rey[i] /= (nshells*nmodes);
    dim[i] /= (nshells*nmodes);

    struc[i] = rey[i];
    struc[i+6] = dim[i];
  }
  
  double Stau = 2.36;
  tke = 0.5*(rey[0] + rey[1] + rey[2]);
  eps = tke/Stau;
  
  eps_init = eps;

  for (int i = 0; i < 3; i++) cir[i] = 2.0*tke - rey[i] - dim[i];
  for (int i = 3; i < 6; i++) cir[i] = -rey[i] - dim[i];

  // time scale
  calcTurbTimeScale();

  // output to screen
  cout << "Simulated TKE: " << tke << endl;
  cout << "Exact TKE    : " << 0.25 << endl;
  cout << "Eps          : " << eps << endl;
  cout << "Tau          : " << tau << endl;
  cout << "Sk/eps       : " << tke/eps << endl;
}

void TurbModel_WVM::bkeuler(double dt, double (*Gn)[3])
{
  double erhs[3], arhs[3], urhs[3];
  double *eipar, *aipar, *uipar;

  for (int ipar = 0; ipar < npar; ipar++)
  {
    eipar = e[ipar];
    aipar = a[ipar];
    uipar = u[ipar];
    
    calcRhs(erhs, arhs, urhs, Gn, eipar, aipar, uipar,dt);
    
    for (int i = 0; i < 3; i++)
    {
      eipar[i] += dt*erhs[i];
      aipar[i] += dt*arhs[i];
      uipar[i] += dt*urhs[i];
    }
  }
}

void TurbModel_WVM::rk4(double dt, double (*Gn)[3], double (*Gnph)[3], 
                        double (*Gnp1)[3])
{
  double erhs_1[3], erhs_2[3], erhs_3[3], erhs_4[3];
  double arhs_1[3], arhs_2[3], arhs_3[3], arhs_4[3];
  double urhs_1[3], urhs_2[3], urhs_3[3], urhs_4[3];

  for (int ipar = 0; ipar < npar; ipar++)
  {
    double *eipar, *aipar, *uipar;
    double etemp[3], atemp[3], utemp[3];

    eipar = e[ipar];
    aipar = a[ipar];
    uipar = u[ipar];

    calcRhs(erhs_1, arhs_1, urhs_1, Gn, eipar, aipar, uipar,dt);

    for (int i = 0; i < 3; i++)
    {
      etemp[i] = eipar[i] + 0.5*dt*erhs_1[i];
      atemp[i] = aipar[i] + 0.5*dt*arhs_1[i];
      utemp[i] = uipar[i] + 0.5*dt*urhs_1[i];
    }

    calcRhs(erhs_2, arhs_2, urhs_2, Gnph, etemp, atemp, utemp,dt);

    for (int i = 0; i < 3; i++)
    { 
      etemp[i] = eipar[i] + 0.5*dt*erhs_2[i];
      atemp[i] = aipar[i] + 0.5*dt*arhs_2[i];
      utemp[i] = uipar[i] + 0.5*dt*urhs_2[i];
    }

    calcRhs(erhs_3, arhs_3, urhs_3, Gnph, etemp, atemp, utemp,dt);
    
    for (int i = 0; i < 3; i++)
    {
      etemp[i] = eipar[i] + dt*erhs_3[i];
      atemp[i] = aipar[i] + dt*arhs_3[i];
      utemp[i] = uipar[i] + dt*urhs_3[i];
    }
    
    calcRhs(erhs_4, arhs_4, urhs_4, Gnp1, etemp, atemp, utemp,dt);
    
    for (int i = 0; i < 3; i++)
    {
      eipar[i] += 1.0/6.0*dt*
	          (erhs_1[i] + 2.0*erhs_2[i] + 2.0*erhs_3[i] + erhs_4[i]);
      aipar[i] += 1.0/6.0*dt*
	          (arhs_1[i] + 2.0*arhs_2[i] + 2.0*arhs_3[i] + arhs_4[i]);
      uipar[i] += 1.0/6.0*dt*
	          (urhs_1[i] + 2.0*urhs_2[i] + 2.0*urhs_3[i] + urhs_4[i]);
    }
  }  
}

void TurbModel_WVM::calcRhs(double *erhs, double *arhs, double *urhs,
                            double (*G)[3], double *eref, double *aref,
                            double *uref, double dt)
{
  // rapid part
  double Gvec[3];
  double vecGvec;

  matTransDotVec3d(Gvec, eref, G);
  vecGvec = vecDotMatDotVec3d(eref, eref, G);
  for (int i = 0; i < 3; i++)
    erhs[i] = -Gvec[i] + vecGvec*eref[i];

  matDotVec3d(Gvec, aref, G);
  vecGvec = vecDotMatDotVec3d(aref, aref, G);
  for (int i = 0; i < 3; i++)
    arhs[i] = Gvec[i] - vecGvec*aref[i];

  matDotVec3d(Gvec, uref, G);
  vecGvec = vecDotMatDotVec3d(eref, uref, G);
  for (int i = 0; i < 3; i++)
    urhs[i] = -Gvec[i] + 2.0*vecGvec*eref[i];

  // expanded gradient
  double r[3][3], d[3][3], rd[3][3];
  double q2 = 2.0*tke;

  r[0][0] = rey[0]/q2;   r[0][1] = rey[3]/q2;   r[0][2] = rey[4]/q2;
  r[1][0] = rey[3]/q2;   r[1][1] = rey[1]/q2;   r[1][2] = rey[5]/q2;
  r[2][0] = rey[4]/q2;   r[2][1] = rey[5]/q2;   r[2][2] = rey[2]/q2;

  d[0][0] = dim[0]/q2;   d[0][1] = dim[3]/q2;   d[0][2] = dim[4]/q2;
  d[1][0] = dim[3]/q2;   d[1][1] = dim[1]/q2;   d[1][2] = dim[5]/q2;
  d[2][0] = dim[4]/q2;   d[2][1] = dim[5]/q2;   d[2][2] = dim[2]/q2;

  matTimesMat3d(rd, r, d);

  matTransDotVec3d(Gvec, eref, rd);
  vecGvec = vecDotMatDotVec3d(eref, eref, rd);
  for (int i = 0; i < 3; i++)
    erhs[i] += -Cn/tau*(Gvec[i] - vecGvec*eref[i]);

  matDotVec3d(Gvec, uref, rd);
  vecGvec = vecDotMatDotVec3d(eref, uref, rd);
  for (int i = 0; i < 3; i++)
    urhs[i] += -Cv/tau*Gvec[i] + (Cv+Cn)/tau*vecGvec*eref[i];

  // slow rotational randomization
  double f[3][3];
  f[0][0] = cir[0]/q2;   f[0][1] = cir[3]/q2;   f[0][2] = cir[4]/q2;
  f[1][0] = cir[3]/q2;   f[1][1] = cir[1]/q2;   f[1][2] = cir[5]/q2;
  f[2][0] = cir[4]/q2;   f[2][1] = cir[5]/q2;   f[2][2] = cir[2]/q2;

  double fnn = vecDotMatDotVec3d(eref, eref, f);
  
  double Ovec[3], Omag, C1, C2;
  Ovec[0] = rd[2][1] - rd[1][2];
  Ovec[1] = rd[0][2] - rd[2][0];
  Ovec[2] = rd[1][0] - rd[0][1];

  Omag = sqrt(vecDotVec3d(Ovec,Ovec));

  C1 = 8.5/tau*Omag*fnn;
  C2 = sqrt(C1);

  double dW[3];
  for (int i = 0; i < 3; i++)
    dW[i] = distribution(generator)/sqrt(dt);

  double umag = sqrt(vecDotVec3d(uref,uref));

  urhs[0] -= C1*uref[0] + C2*umag*(dW[2]*eref[1] - dW[1]*eref[2]);
  urhs[1] -= C1*uref[1] + C2*umag*(dW[0]*eref[2] - dW[2]*eref[0]);
  urhs[2] -= C1*uref[2] + C2*umag*(dW[1]*eref[0] - dW[0]*eref[1]);

}

double TurbModel_WVM::updateDissipation(double (*Gn)[3], double dt)
{
  double prod = -rey[0]*Gn[0][0] - rey[3]*Gn[0][1] - rey[4]*Gn[0][2]
                -rey[3]*Gn[1][0] - rey[1]*Gn[1][1] - rey[5]*Gn[1][2]
                -rey[4]*Gn[2][0] - rey[5]*Gn[2][1] - rey[2]*Gn[2][2];
  
  double vort[3];
  vort[0] = Gn[2][1] - Gn[1][2];
  vort[1] = Gn[0][2] - Gn[2][0];
  vort[2] = Gn[1][0] - Gn[0][1];

  double d[3][3];
  double q2 = 2.0*tke;
  d[0][0] = dim[0]/q2;   d[0][1] = dim[3]/q2;   d[0][2] = dim[4]/q2;
  d[1][0] = dim[3]/q2;   d[1][1] = dim[1]/q2;   d[1][2] = dim[5]/q2;
  d[2][0] = dim[4]/q2;   d[2][1] = dim[5]/q2;   d[2][2] = dim[2]/q2;

  double OdO = vecDotMatDotVec3d(vort, vort, d);
  
  double rhs = Ceps1*eps/tke*prod - Ceps2*eps*eps/tke - Ceps3*eps*sqrt(OdO);

  eps += dt*rhs;

  return eps;
}

void TurbModel_WVM::calcTurbTimeScale()
{
  double r[3][3], d[3][3], rd[3][3];
  double q2 = 2.0*tke;

  r[0][0] = rey[0]/q2;   r[0][1] = rey[3]/q2;   r[0][2] = rey[4]/q2;
  r[1][0] = rey[3]/q2;   r[1][1] = rey[1]/q2;   r[1][2] = rey[5]/q2;
  r[2][0] = rey[4]/q2;   r[2][1] = rey[5]/q2;   r[2][2] = rey[2]/q2;

  d[0][0] = dim[0]/q2;   d[0][1] = dim[3]/q2;   d[0][2] = dim[4]/q2;
  d[1][0] = dim[3]/q2;   d[1][1] = dim[1]/q2;   d[1][2] = dim[5]/q2;
  d[2][0] = dim[4]/q2;   d[2][1] = dim[5]/q2;   d[2][2] = dim[2]/q2;
 
  matTimesMat3d(rd, r, d);
  tau = 0.0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      tau += rd[i][j]*r[j][i];
  tau *= 2.0*Cv*tke/eps;
}

void TurbModel_WVM::calcReStress(double *struc, double (*Gn)[3], 
                                 double (*Gnph)[3], double (*Gnp1)[3], 
                                 double dt)
{  
  // time integration
  bkeuler(dt, Gn);
  //rk4(dt, Gn, Gnph, Gnp1);

  // structure tensors and tke
  for (int i = 0; i < 6; i++)
  {
    rey[i] = 0.0;
    dim[i] = 0.0;
  }

  for (int ipar = 0; ipar < npar; ipar ++)
  {
    rey[0] += u[ipar][0]*u[ipar][0];
    rey[1] += u[ipar][1]*u[ipar][1];
    rey[2] += u[ipar][2]*u[ipar][2];
    rey[3] += u[ipar][0]*u[ipar][1];
    rey[4] += u[ipar][0]*u[ipar][2];
    rey[5] += u[ipar][1]*u[ipar][2];

    double umag2 = vecDotVec3d(u[ipar],u[ipar]);
    dim[0] += umag2*e[ipar][0]*e[ipar][0];
    dim[1] += umag2*e[ipar][1]*e[ipar][1];
    dim[2] += umag2*e[ipar][2]*e[ipar][2];
    dim[3] += umag2*e[ipar][0]*e[ipar][1];
    dim[4] += umag2*e[ipar][0]*e[ipar][2];
    dim[5] += umag2*e[ipar][1]*e[ipar][2];

    /*double emag = e[ipar][0]*e[ipar][0] + 
                  e[ipar][1]*e[ipar][1] +
                  e[ipar][2]*e[ipar][2];

    if ((emag > 1.0001) || (emag < 0.9999))
      cout << "Error, |e|^2 = " << emag << endl;*/
  }

  for (int i = 0; i < 6; i++)
  {
    rey[i] /= (nshells*nmodes);
    dim[i] /= (nshells*nmodes);

    struc[i] = rey[i];
    struc[i+6] = dim[i];
  }
  
  tke = 0.5*(rey[0] + rey[1] + rey[2]);

  //cout << "tke rij: " << tke << endl;
  //cout << "tke dij: " << 0.5*(dim[0] + dim[1] + dim[2]) << endl;

  for (int i = 0; i < 3; i++) cir[i] = 2.0*tke - rey[i] - dim[i];
  for (int i = 3; i < 6; i++) cir[i] = -rey[i] - dim[i];

  calcTurbTimeScale();
}

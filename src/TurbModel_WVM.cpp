#include "TurbModel_WVM.hpp"

void TurbModel_WVM::bkeuler(int ipar, double dt, double (*Gn)[3])
{
  double erhs[3], arhs[3], urhs[3];
  double *eipar, *aipar, *uipar;

  eipar = e[ipar];
  aipar = a[ipar];
  uipar = u[ipar];

  calcRhs(erhs, arhs, urhs, Gn, eipar, aipar, uipar);

  for (int i = 0; i < 3; i++)
  {
    eipar[i] += dt*erhs[i];
    aipar[i] += dt*arhs[i];
    uipar[i] += dt*urhs[i];
  }
}

void TurbModel_WVM::rk4(int ipar, double dt, double (*Gn)[3], 
                        double (*Gnph)[3], double (*Gnp1)[3])
{
  double erhs_1[3], erhs_2[3], erhs_3[3], erhs_4[3];
  double arhs_1[3], arhs_2[3], arhs_3[3], arhs_4[3];
  double urhs_1[3], urhs_2[3], urhs_3[3], urhs_4[3];

  double *eipar, *aipar, *uipar;
  double etemp[3], atemp[3], utemp[3];

  eipar = e[ipar];
  aipar = a[ipar];
  uipar = u[ipar];

  calcRhs(erhs_1, arhs_1, urhs_1, Gn, eipar, aipar, uipar);

  for (int i = 0; i < 3; i++)
  {
    etemp[i] = eipar[i] + 0.5*dt*erhs_1[i];
    atemp[i] = aipar[i] + 0.5*dt*arhs_1[i];
    utemp[i] = uipar[i] + 0.5*dt*urhs_1[i];
  }

  calcRhs(erhs_2, arhs_2, urhs_2, Gnph, etemp, atemp, utemp);

  for (int i = 0; i < 3; i++)
  { 
    etemp[i] = eipar[i] + 0.5*dt*erhs_2[i];
    atemp[i] = aipar[i] + 0.5*dt*arhs_2[i];
    utemp[i] = uipar[i] + 0.5*dt*urhs_2[i];
  }

  calcRhs(erhs_3, arhs_3, urhs_3, Gnph, etemp, atemp, utemp);

  for (int i = 0; i < 3; i++)
  {
    etemp[i] = eipar[i] + dt*erhs_3[i];
    atemp[i] = aipar[i] + dt*arhs_3[i];
    utemp[i] = uipar[i] + dt*urhs_3[i];
  }

  calcRhs(erhs_4, arhs_4, urhs_4, Gnp1, etemp, atemp, utemp);

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

void TurbModel_WVM::calcRhs(double *erhs, double *arhs, double *urhs, 
                            double (*G)[3], double *eref, double *aref, 
                            double *uref)
{
  double tempA[3];
  double tempB;
  
  matTransDotVec3d(tempA, eref, G);
  tempB = vecDotMatDotVec3d(eref, eref, G);  
  for (int i = 0; i < 3; i++)
    erhs[i] = -tempA[i] + tempB*eref[i];

  matDotVec3d(tempA, aref, G);
  tempB = vecDotMatDotVec3d(aref, aref, G);
  for (int i = 0; i < 3; i++)
    arhs[i] = tempA[i] - tempB*aref[i];

  matDotVec3d(tempA, uref, G);
  tempB = vecDotMatDotVec3d(eref, uref, G);
  for (int i = 0; i < 3; i++)
    urhs[i] = -tempA[i] + 2.0*tempB*eref[i];
}

TurbModel_WVM::TurbModel_WVM(const char *name)
{
  nshells = 32;
  nmodes = 2000;
  npar = nshells*nmodes;

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

void TurbModel_WVM::initialHookScalarRansTurbModel(double *rey)
{
  srand(time(NULL));

  // Energy spectrum
  double mag_k[nshells], dk[nshells], Enukdk[nshells];
  double fkt = pow(4.0,0.2), Cnuk = 0.4843, qsqd = 0.5;

  mag_k[0] = 0.01;
  for (int i = 0; i < nshells; i++)
    mag_k[i] = pow(fkt,i)*mag_k[0];

  for (int i = 1; i < (nshells - 1); i++)
    dk[i] = (mag_k[i+1] - mag_k[i-1])/2.0;
  dk[0] = 2.0*dk[1] - dk[2];
  dk[nshells-1] = 2.0*dk[nshells-2] - dk[nshells-3];
  
  for (int i = 0; i < nshells; i++)
  {
    double Enuk = qsqd*Cnuk*
                  pow(mag_k[i],4.0)/pow(1.0 + mag_k[i]*mag_k[i], 17.0/6.0);
    Enukdk[i] = Enuk*dk[i];  
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

      double umag = sqrt(2.0*nshells*Enukdk[ishells]);
      for (int i = 0; i < 3; i++)
	u[ipar][i] = umag*( cos(theta)*r[i] + sin(theta)*s[i] );

      // initial Reynolds stresses
      rey[0] += u[ipar][0]*u[ipar][0];
      rey[1] += u[ipar][1]*u[ipar][1];
      rey[2] += u[ipar][2]*u[ipar][2];
      rey[3] += u[ipar][0]*u[ipar][1];
      rey[4] += u[ipar][0]*u[ipar][2];
      rey[5] += u[ipar][1]*u[ipar][2];
    }

  // normalize
  for (int i = 0; i < 6; i++)
    rey[i] /= (nshells*nmodes);

  cout << "Simulated TKE: " << 0.5*(rey[0] + rey[1] + rey[2]) << endl;
  cout << "Exact TKE    : " << 0.5*qsqd << endl;
}

void TurbModel_WVM::calcReStress(double *rey, double (*Gn)[3], 
                                 double (*Gnph)[3], double (*Gnp1)[3], 
                                 double dt)
{
  // loop through all the particles
  for (int ipar = 0; ipar < npar; ipar++)
  {
    rk4(ipar, dt, Gn, Gnph, Gnp1);

    rey[0] += u[ipar][0]*u[ipar][0];
    rey[1] += u[ipar][1]*u[ipar][1];
    rey[2] += u[ipar][2]*u[ipar][2];
    rey[3] += u[ipar][0]*u[ipar][1];
    rey[4] += u[ipar][0]*u[ipar][2];
    rey[5] += u[ipar][1]*u[ipar][2];
  }

  // normalize
  for (int i = 0; i < 6; i++)
    rey[i] /= (nshells*nmodes);
}

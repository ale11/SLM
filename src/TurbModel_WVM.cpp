#include "TurbModel_WVM.hpp"

/*##############################################################################
#
# Stochastic Lagrangian Models
#
##############################################################################*/
TurbModel_SLM::TurbModel_SLM(void) {}

TurbModel_SLM::~TurbModel_SLM() {}

void TurbModel_SLM::initialHookScalarRansTurbModel(double Stau, double *struc, 
						   double &eps_init, double *prod, 
						   double *rRedi, double (*Gn)[3]) {}

void TurbModel_SLM::calcReStress(double *struc, double &eps_main, double *prod, double *rRedi,
				 double *sRediEps, double (*Gn)[3], double (*Gnph)[3],
				 double (*Gnp1)[3], double dt, char* tIntName) {}

void TurbModel_SLM::writeData(double St) {};

/*############################################################################################
#
# Particle Simulation of SLM
#
############################################################################################*/

TurbModel_PAR::TurbModel_PAR(void) : TurbModel_SLM()
{
  nshells = 32;
  nmodes = 2000;
  npar = nshells*nmodes;

  e = new double[npar][3];
  a = new double[npar][3];
  u = new double[npar][3];
}

TurbModel_PAR::~TurbModel_PAR()
{
  delete[] e;
  delete[] a;
  delete[] u;
}

void TurbModel_PAR::initialHookScalarRansTurbModel(double Stau, double *struc,
						   double &eps_init, 
                                                   double *prod, 
                                                   double *rRedi, 
                                                   double (*Gn)[3])
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
    }
  
  calcTurbStatistics();
  eps = tke/Stau;
  calcTurbTimeScale();

  // output to main
  for (int i = 0; i < 6; i++)
  {
    struc[i] = rey[i];
    struc[i+6] = dim[i];
  }
  eps_init = eps;

  int index[3][3] = { {0, 3, 4}, {3, 1, 5}, {4, 5, 2} };
  // production
  for (int i = 0; i < 3; i++)
    for (int j = i; j < 3; j++)
      for (int k = 0; k < 3; k++)
      {
	prod[index[i][j]] = -rey[index[i][k]]*Gn[j][k]
	                    -rey[index[j][k]]*Gn[i][k];
      }

  // rapid redistribution
  for (int i = 0; i < 3; i++)
    for (int j = i; j < 3; j++)
    {
      rRedi[index[i][j]] = 0.0;
      for (int k = 0; k < 3; k++)
	for (int l = 0; l < 3; l++)
          rRedi[index[i][j]] += 2.0*Gn[k][l]*(M[i][l][j][k] + M[j][l][i][k]);
    }

  // output to screen
  cout << "Simulated TKE : " << tke << endl;
  cout << "Exact TKE     : " << 0.25 << endl;
  cout << "Eps           : " << eps << endl;
  cout << "Tau           : " << tau << endl;
  cout << "Sk/eps        : " << tke/eps << endl;
  cout << "--------------------------------" << endl;
}

void TurbModel_PAR::fwEuler(double dt, double (*Gn)[3])
{
  for (int ipar = 0; ipar < npar; ipar++)
  {
    double *eipar = e[ipar];
    double *uipar = u[ipar];

    double drift[6], diff[6][6];
    driftCoeff(drift, Gn, eipar, uipar);
    diffCoeff(diff, eipar, uipar);

    double dW, bdW[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    for (int i = 0; i < 6; i++)
    {
      dW = distribution(generator)*sqrt(dt);
      for (int j = 0; j < 6; j++)
        bdW[j] += diff[j][i]*dW;
    }

    for (int i = 0; i < 3; i++)
    {
      eipar[i] += drift[i]*dt + bdW[i];
      uipar[i] += drift[i+3]*dt + bdW[i+3];
    }

    normVec3d(eipar);
    VecOrthoVec3d(uipar,eipar);
    //double small = 1.0e-1;
    //double e2 = vecDotVec3d(eipar, eipar);
    //if ((e2 < (1.0 - small)) || (e2 > (1.0 + small)))
    //  cout << "Wrong e2: " << e2 << endl;
  }

  double eps_rhs = rhsDissipation(Gn, dt);
  eps += dt*eps_rhs;

  calcTurbStatistics();
  calcTurbTimeScale();
}

void TurbModel_PAR::semiHeun(double dt, double (*Gn)[3])
{
  for (int ipar = 0; ipar < npar; ipar++)
  {    
    double *evec = e[ipar];
    double *uvec = u[ipar];

    double drift_0[6], diff_0[6][6]; 
    driftCoeff(drift_0, Gn, evec, uvec);
    diffCoeff(diff_0, evec, uvec);

    double dW_0, bdW_0[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    for (int i = 0; i < 6; i++)
    {
      dW_0 = distribution(generator)*sqrt(dt);
      for (int j = 0; j < 6; j++)
        bdW_0[j] += diff_0[j][i]*dW_0;
    }

    double ebar[3], ubar[3];
    for (int i = 0; i < 3; i++)
    {
      ebar[i] = evec[i] + drift_0[i]*dt + bdW_0[i];
      ubar[i] = uvec[i] + drift_0[i+3]*dt + bdW_0[i+3];
    }
    //normVec3d(ebar); // forced renormalization
    double drift_bar[6];
    driftCoeff(drift_bar, Gn, ebar, ubar);

    for (int i = 0; i < 3; i++)
    {
      evec[i] += 0.5*(drift_bar[i] + drift_0[i])*dt + bdW_0[i];
      uvec[i] += 0.5*(drift_bar[i+3] + drift_0[i+3])*dt + bdW_0[i+3];
    }

    normVec3d(evec);
    VecOrthoVec3d(uvec,evec);
    //double small = 1.0e-1;
    //double e2 = vecDotVec3d(evec, evec);
    //if ((e2 < (1.0 - small)) || (e2 > (1.0 + small)))
    //  cout << "Wrong e2: " << e2 << endl;
    //double ue = vecDotVec3d(uvec, evec);
    //if ( fabs(ue) > small)
    //	cout << "No ortho: " << ue << endl;
  }

  double eps_rhs = rhsDissipation(Gn, dt);
  eps += dt*eps_rhs;

  calcTurbStatistics();
  calcTurbTimeScale();
}

void TurbModel_PAR::Heun(double dt, double (*Gn)[3])
{
  for (int ipar = 0; ipar < npar; ipar++)
  {    
    double *evec = e[ipar];
    double *uvec = u[ipar];

    // -------------------------------------------------------------------------
    // Deterministic component
    // -------------------------------------------------------------------------
    double drift_0[6], diff_0[6][6]; 
    driftCoeff(drift_0, Gn, evec, uvec);
    diffCoeff(diff_0, evec, uvec);

    double dW_0, bdW_0[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    for (int i = 0; i < 6; i++)
    {
      dW_0 = distribution(generator)*sqrt(dt);
      for (int j = 0; j < 6; j++)
        bdW_0[j] += diff_0[j][i]*dW_0;
    }

    double ebar[3], ubar[3];
    for (int i = 0; i < 3; i++)
    {
      ebar[i] = evec[i] + drift_0[i]*dt + bdW_0[i];
      ubar[i] = uvec[i] + drift_0[i+3]*dt + bdW_0[i+3];
    }

    double drift_bar[6];
    driftCoeff(drift_bar, Gn, ebar, ubar);

    // -------------------------------------------------------------------------
    // Stochastic component
    // -------------------------------------------------------------------------
    double dW[3];
    for (int i = 0; i < 3; i++)
      dW[i] = distribution(generator)*sqrt(dt);

    double Vij[3][3];
    Vij[0][0] = -dt;   
    Vij[1][1] = -dt;
    Vij[2][2] = -dt;
    
    double uniform;
    uniform = (double)rand()/RAND_MAX;
    if (uniform < 0.5)   Vij[1][0] = -dt;
    else                 Vij[1][0] = dt;
    uniform = (double)rand()/RAND_MAX;
    if (uniform < 0.5)   Vij[2][0] = -dt;
    else                 Vij[2][0] = dt;
    uniform = (double)rand()/RAND_MAX;
    if (uniform < 0.5)   Vij[2][1] = -dt;
    else                 Vij[2][1] = dt;

    Vij[0][1] = -Vij[1][0];
    Vij[0][2] = -Vij[2][0];
    Vij[1][2] = -Vij[2][1];

    double bdW[3] = {0.0, 0.0, 0.0};
    /*double diff_j[3][3], rje[3], rju[3], ure[3], uru[3];
    
    // first sum
    for (int j = 0; j < 3; j++)
    {
      for (int i = 0; i < 3; i++)
      {
	rje[i] = evec[i] + drift_0[i]*dt;
	rju[i] = uvec[i] + drift_0[i+3]*dt + diff_0[i][j]*sqrt(dt);
      }
      
      diffCoeff(diff_j, rje, rju);
      
      for (int i = 0; i < 3; i++)
	bdW[i] += diff_j[i][j]*dW[j];

      for (int i = 0; i < 3; i++)
      {
	rje[i] = evec[i] + drift_0[i]*dt;
	rju[i] = uvec[i] + drift_0[i+3]*dt - diff_0[i][j]*sqrt(dt);
      }

      diffCoeff(diff_j, rje, rju);

      for (int i = 0; i < 3; i++)
	bdW[i] += diff_j[i][j]*dW[j];

      for (int i = 0; i < 3; i++)
	bdW[i] += 2.0*diff_0[i][j]*dW[j];

      for (int r = 0; r < 3; r++)
      {
	if (r != j)
	{
	  for (int i = 0; i < 3; i++)
	  {
	    ure[i] = evec[i];
	    uru[i] = uvec[i] + diff_0[i][r]*sqrt(dt);
	  }

	  diffCoeff(diff_j, ure, uru);

	  for (int i = 0; i < 3; i++)
	    bdW[i] += diff_j[i][j]*dW[j]/sqrt(dt);

	  for (int i = 0; i< 3; i++)
	  {
	    ure[i] = evec[i];
	    uru[i] = uvec[i] - diff_0[i][r]*sqrt(dt);
	  }

	  diffCoeff(diff_j, ure, uru);

	  for (int i = 0; i< 3; i++)
	    bdW[i] += diff_j[i][j]*dW[j]/sqrt(dt);

	  for (int i = 0; i< 3; i++)
	    bdW[i] -= 2.0*diff_0[i][j]*dW[j]/sqrt(dt);
	} // end: if(r != j)
      } // end: for (r = 0; r < 3; r++)
    } // end: for (j = 0; j < 3; j++)

    // second sum
    for (int j = 0; j < 3; j++)
    {
      for (int i = 0; i < 3; i++)
      {
	rje[i] = evec[i] + drift_0[i]*dt;
	rju[i] = uvec[i] + drift_0[i+3]*dt + diff_0[i][j]*sqrt(dt);
      }
      
      diffCoeff(diff_j, rje, rju);
      
      for (int i = 0; i < 3; i++)
	bdW[i] += diff_j[i][j]*(dW[j]*dW[j] - dt)/sqrt(dt);

      for (int i = 0; i < 3; i++)
      {
	rje[i] = evec[i] + drift_0[i]*dt;
	rju[i] = uvec[i] + drift_0[i+3]*dt - diff_0[i][j]*sqrt(dt);
      }

      diffCoeff(diff_j, rje, rju);

      for (int i = 0; i < 3; i++)
	bdW[i] -= diff_j[i][j]*(dW[j]*dW[j] - dt)/sqrt(dt);

      for (int r = 0; r < 3; r++)
      {
	if (r != j)
	{
	  for (int i = 0; i < 3; i++)
	  {
	    ure[i] = evec[i];
	    uru[i] = uvec[i] + diff_0[i][r]*sqrt(dt);
	  }

	  diffCoeff(diff_j, ure, uru);

	  for (int i = 0; i < 3; i++)
	    bdW[i] += diff_j[i][j]*(dW[j]*dW[r] + Vij[r][j])/sqrt(dt);

	  for (int i = 0; i< 3; i++)
	  {
	    ure[i] = evec[i];
	    uru[i] = uvec[i] - diff_0[i][r]*sqrt(dt);
	  }

	  diffCoeff(diff_j, ure, uru);

	  for (int i = 0; i< 3; i++)
	    bdW[i] -= diff_j[i][j]*(dW[j]*dW[r] + Vij[r][j])/sqrt(dt);

	} // end: if(r != j)
      } // end: for (r = 0; r < 3; r++)
    } // end: for (j = 0; j < 3; j++)
    */
    // -------------------------------------------------------------------------
    // Update Solution
    // -------------------------------------------------------------------------
    for (int i = 0; i < 3; i++)
    {
      evec[i] += 0.5*(drift_bar[i] + drift_0[i])*dt;
      uvec[i] += 0.5*(drift_bar[i+3] + drift_0[i+3])*dt + 0.25*bdW[i];
    }
  }

  double eps_rhs = rhsDissipation(Gn, dt);
  eps += dt*eps_rhs;

  calcTurbStatistics();
  calcTurbTimeScale();
}

void TurbModel_PAR::driftCoeff(double *drift, double (*G)[3], double *eref,
                               double *uref)
{
  for (int i = 0; i < 6; i++)
    drift[i] = 0.0;
}

void TurbModel_PAR::diffCoeff(double (*diff)[6], double *eref, double *uref)
{
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      diff[i][j] = 0.0;
}

double TurbModel_PAR::rhsDissipation(double(*Gn)[3], double dt)
{
  return 0.0;
}

void TurbModel_PAR::calcTurbTimeScale()
{
  tau = 0.0;
}

void TurbModel_PAR:: calcTurbStatistics()
{
  // initialize to zero
  for (int i = 0; i < 6; i++)
  {
    rey[i] = 0.0;
    dim[i] = 0.0;
  }
  for (int i = 0; i < 3; i++)
    for (int j = i; j < 3; j++)
    {
      e2[i][j] = 0.0;
      for (int k = 0; k < 3; k++)
      	for (int l = k; l < 3; l++)
      	{
      		M[i][j][k][l] = 0.0;
      		e4[i][j][k][l] = 0.0;
      	}
    }

  // reynolds stresses and dimensionality
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

    for (int i = 0; i < 3; i++)
      for(int j = i; j <3; j++)
      {
      	e2[i][j] += e[ipar][i]*e[ipar][j];
      	for (int k = 0; k< 3; k++)
      		for (int l = k; l < 3; l++)
      		{
      			M[i][j][k][l] += u[ipar][i]*u[ipar][j]*e[ipar][k]*e[ipar][l];
      			e4[i][j][k][l] += e[ipar][i]*e[ipar][j]*e[ipar][k]*e[ipar][l];
      		}
      }
  }

  // normalize
  for (int i = 0; i < 6; i++)
  {
    rey[i] /= (nshells*nmodes);
    dim[i] /= (nshells*nmodes);
  }

  tke = 0.5*(rey[0] + rey[1] + rey[2]);

  for (int i = 0; i < 3; i++) cir[i] = 2.0*tke - rey[i] - dim[i];
  for (int i = 3; i < 6; i++) cir[i] = -rey[i] - dim[i];

  for (int i = 0; i < 3; i++)
    for(int j = i; j <3; j++)
    {
      e2[i][j] /= (nshells*nmodes);
      for (int k = 0; k< 3; k++)
      	for (int l = k;l < 3; l++)
      	{
      		M[i][j][k][l] /= (nshells*nmodes);
      		e4[i][j][k][l] /= (nshells*nmodes);
      	}
    }

  // make symmetric
  for (int i = 0; i < 3; i++)
    for (int j = i; j < 3; j++)
    {
      e2[i][j] = e2[j][i];

      M[i][j][1][0] = M[i][j][0][1];
      M[i][j][2][0] = M[i][j][0][2];
      M[i][j][2][1] = M[i][j][1][2];

      e4[i][j][1][0] = e4[i][j][0][1];
      e4[i][j][2][0] = e4[i][j][0][2];
      e4[i][j][2][1] = e4[i][j][1][2];
    }

  for (int k = 0; k < 3; k++)
    for (int l = 0; l < 3; l++)
    {
      M[1][0][k][l] = M[0][1][k][l];
      M[2][0][k][l] = M[0][2][k][l];
      M[2][1][k][l] = M[1][2][k][l];

      e4[1][0][k][l] = e4[0][1][k][l];
      e4[2][0][k][l] = e4[0][2][k][l];
      e4[2][1][k][l] = e4[1][2][k][l];
    }
}

void TurbModel_PAR::writeData(double St)
{
	FILE *fid;
	string fname = "pData." + to_string(St) + ".dat";

  fid = fopen(fname.c_str(), "w");
  if (fid != NULL)
  {
  	fprintf(fid, "e[0]\te[1]\te[2]\tu[0]\tu[1]\tu[2]\n");
    for (int ipar = 0; ipar < npar; ipar++)
    {
      double *evec = e[ipar];
      double *uvec = u[ipar];

      fprintf(fid, "%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n",
    	  evec[0], evec[1], evec[2], uvec[0], uvec[1], uvec[2]);
    }
  }
  else
  	cout << "Could not open file " << fname << endl;

  fclose(fid);

}

void TurbModel_PAR::calcReStress(double *struc, double &eps_main, double *prod,
                                 double *rRedi, double *sRediEps, 
                                 double (*Gn)[3], double (*Gnph)[3], 
                                 double (*Gnp1)[3], double dt, char* tIntName)
{  
  // time integration
  if      (strcmp(tIntName, "fwEuler") == 0)   fwEuler(dt, Gn);
  else if (strcmp(tIntName, "Heun") == 0)      Heun(dt, Gn);
  else if (strcmp(tIntName, "semiHeun") == 0)  semiHeun(dt, Gn);
  else    cout << "Not a valid time integration scheme" << endl;
  
  // structure tensors and eps to main
  for (int i = 0; i < 6; i++)
  {
    struc[i] = rey[i];
    struc[i+6] = dim[i];
  }
  eps_main = eps;

  int index[3][3] = { {0, 3, 4}, {3, 1, 5}, {4, 5, 2} };
  // production
  for (int i = 0; i < 3; i++)
    for (int j = i; j < 3; j++)
      for (int k = 0; k < 3; k++)
      {
	prod[index[i][j]] = -rey[index[i][k]]*Gn[j][k] 
                            -rey[index[j][k]]*Gn[i][k];
      }

  // rapid redistribution
  for (int i = 0; i < 3; i++)
    for (int j = i; j < 3; j++)
    {
      rRedi[index[i][j]] = 0.0;
      for (int k = 0; k < 3; k++)
	for (int l = 0; l < 3; l++)
	  rRedi[index[i][j]] += 2.0*Gn[k][l]*(M[i][l][j][k] + M[j][l][i][k]);
    }

  // slow redistribution + dissipation
  double RD[3][3] = { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0} };
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
	RD[i][j] += rey[index[i][k]]*dim[index[k][j]];

  for (int i = 0; i < 3; i++)
    for (int j = i; j < 3; j++)
    {
      sRediEps[index[i][j]] = 0.0;      
    }

}

/*##############################################################################
#        
# Interacting Particle Representation Model
#
##############################################################################*/

TurbModel_IPRM::TurbModel_IPRM() : TurbModel_PAR()
{
	cout << "Model         : IPRM" << endl;

  Cn = 2.2;
  Cv = 1.0; //1.5;
  Ceps1 = 1.5;
  Ceps2 = 11.0/6.0;
  Ceps3 = 0.01;
}

void TurbModel_IPRM::driftCoeff(double *drift, double (*G)[3], double *eref,
                               double *uref)
{
  // rapid part
  double Gvec[3];
  double vecGvec;

  matTransDotVec3d(Gvec, eref, G);
  vecGvec = vecDotMatDotVec3d(eref, eref, G);
  for (int i = 0; i < 3; i++)
    drift[i] = -Gvec[i] + vecGvec*eref[i];

  matDotVec3d(Gvec, uref, G);
  vecGvec = vecDotMatDotVec3d(eref, uref, G);
  for (int i = 0; i < 3; i++)
    drift[i+3] = -Gvec[i] + 2.0*vecGvec*eref[i];

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
    drift[i] += -Cn/tau*(Gvec[i] - vecGvec*eref[i]);

  matDotVec3d(Gvec, uref, rd);
  vecGvec = vecDotMatDotVec3d(eref, uref, rd);
  for (int i = 0; i < 3; i++)
    drift[i+3] += -Cv/tau*Gvec[i] + (Cv+Cn)/tau*vecGvec*eref[i];

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

  drift[3] -= C1*uref[0];
  drift[4] -= C1*uref[1];
  drift[5] -= C1*uref[2];
}

void TurbModel_IPRM::diffCoeff(double (*diff)[6], double *eref, double *uref)
{
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++)
      diff[i][j] = 0.0;

  // slow rotational randomization
  double r[3][3], d[3][3], rd[3][3];
  double q2 = 2.0*tke;

  r[0][0] = rey[0]/q2;   r[0][1] = rey[3]/q2;   r[0][2] = rey[4]/q2;
  r[1][0] = rey[3]/q2;   r[1][1] = rey[1]/q2;   r[1][2] = rey[5]/q2;
  r[2][0] = rey[4]/q2;   r[2][1] = rey[5]/q2;   r[2][2] = rey[2]/q2;

  d[0][0] = dim[0]/q2;   d[0][1] = dim[3]/q2;   d[0][2] = dim[4]/q2;
  d[1][0] = dim[3]/q2;   d[1][1] = dim[1]/q2;   d[1][2] = dim[5]/q2;
  d[2][0] = dim[4]/q2;   d[2][1] = dim[5]/q2;   d[2][2] = dim[2]/q2;

  matTimesMat3d(rd, r, d);

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

  double umag = sqrt(vecDotVec3d(uref,uref));
  diff[3][0] = C2*umag*(1.0 - eref[0]*eref[0]);
  diff[3][1] = C2*umag*(-eref[0]*eref[1]);
  diff[3][2] = C2*umag*(-eref[0]*eref[2]);

  diff[4][0] = C2*umag*(-eref[1]*eref[0]);
  diff[4][1] = C2*umag*(1.0 - eref[1]*eref[1]);
  diff[4][2] = C2*umag*(-eref[1]*eref[2]);

  diff[5][0] = C2*umag*(-eref[2]*eref[0]);
  diff[5][1] = C2*umag*(-eref[2]*eref[1]);
  diff[5][2] = C2*umag*(1.0 - eref[2]*eref[2]);

  /*diff[3][0] = C2*sqrt(2.0*tke)*(1.0 - eref[0]*eref[0]);
  diff[3][1] = C2*sqrt(2.0*tke)*(-eref[0]*eref[1]);
  diff[3][2] = C2*sqrt(2.0*tke)*(-eref[0]*eref[2]);

  diff[4][0] = C2*sqrt(2.0*tke)*(-eref[1]*eref[0]);
  diff[4][1] = C2*sqrt(2.0*tke)*(1.0 - eref[1]*eref[1]);
  diff[4][2] = C2*sqrt(2.0*tke)*(-eref[1]*eref[2]);

  diff[5][0] = C2*sqrt(2.0*tke)*(-eref[2]*eref[0]);
  diff[5][1] = C2*sqrt(2.0*tke)*(-eref[2]*eref[1]);
  diff[5][2] = C2*sqrt(2.0*tke)*(1.0 - eref[2]*eref[2]);*/

}

double TurbModel_IPRM::rhsDissipation(double(*Gn)[3], double dt)
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

  return rhs;
}

void TurbModel_IPRM::calcTurbTimeScale()
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

/*##############################################################################
#        
# Langevin Velocity Model
#
##############################################################################*/

TurbModel_LANG::TurbModel_LANG() : TurbModel_PAR()
{
	cout << "Model         : LANG" << endl;

  ae = 0.03;
  au = 2.1;
  gamma = 2.0;
  Ceps1 = 1.5;
  Ceps2 = 11.0/6.0;
  Ceps3 = 0.01;
}

void TurbModel_LANG::driftCoeff(double *drift, double (*G)[3], double *eref,
                                double *uref)
{
  // rapid part
  double Gvec[3];
  double vecGvec;

  matTransDotVec3d(Gvec, eref, G);
  vecGvec = vecDotMatDotVec3d(eref, eref, G);
  for (int i = 0; i < 3; i++)
    drift[i] = -Gvec[i] + vecGvec*eref[i];

  matDotVec3d(Gvec, uref, G);
  vecGvec = vecDotMatDotVec3d(eref, uref, G);
  for (int i = 0; i < 3; i++)
    drift[i+3] = -Gvec[i] + 2.0*vecGvec*eref[i];

  // slow part
  double b[3][3];
  double q2 = 2.0*tke;

  b[0][0] = rey[0]/q2 - 1.0/3.0;   
  b[0][1] = rey[3]/q2;   
  b[0][2] = rey[4]/q2;
  
  b[1][0] = rey[3]/q2;   
  b[1][1] = rey[1]/q2 - 1.0/3.0;   
  b[1][2] = rey[5]/q2;
  
  b[2][0] = rey[4]/q2;   
  b[2][1] = rey[5]/q2;   
  b[2][2] = rey[2]/q2 - 1.0/3.0;

  double u2 = vecDotVec3d(uref, uref);
  double ebe = vecDotMatDotVec3d(eref, eref, b);
  double be[3]; matDotVec3d(be, eref, b);
  double bu[3]; matDotVec3d(bu, uref, b);
  double b2 = traceMatTimesMat3d(b, b);

  for (int i = 0; i < 3; i++)
  {
    drift[i] += -0.5*eps/tke*(ae + au*tke/u2)*eref[i]
                -gamma*eps/tke*(be[i] - eref[i]*ebe);

    drift[i+3] += -0.5*eps/tke*(1.0 + 3.0/2.0*au)*uref[i]
                  +gamma*eps/tke*(bu[i] - b2*uref[i]);
  }
}

void TurbModel_LANG::diffCoeff(double (*diff)[6], double *eref, double *uref)
{
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++)
      diff[i][j] = 0.0;

  double u2 = vecDotVec3d(uref,uref);
  double c1 = sqrt(au*eps);
  double c2 = sqrt(ae*eps/tke);

  // wave vector
  diff[0][0] = -c1*uref[0]*eref[0]/u2;
  diff[0][1] = -c1*uref[0]*eref[1]/u2;
  diff[0][2] = -c1*uref[0]*eref[2]/u2;

  diff[1][0] = -c1*uref[1]*eref[0]/u2;
  diff[1][1] = -c1*uref[1]*eref[1]/u2;
  diff[1][2] = -c1*uref[1]*eref[2]/u2;

  diff[2][0] = -c1*uref[2]*eref[0]/u2;
  diff[2][1] = -c1*uref[2]*eref[1]/u2;
  diff[2][2] = -c1*uref[2]*eref[2]/u2;
  
  diff[0][3] = c2*(1.0 - eref[0]*eref[0] - uref[0]*uref[0]/u2);
  diff[0][4] = c2*(    - eref[0]*eref[1] - uref[0]*uref[1]/u2);
  diff[0][5] = c2*(    - eref[0]*eref[2] - uref[0]*uref[2]/u2);
  
  diff[1][3] = c2*(    - eref[1]*eref[0] - uref[1]*uref[0]/u2);
  diff[1][4] = c2*(1.0 - eref[1]*eref[1] - uref[1]*uref[1]/u2);
  diff[1][5] = c2*(    - eref[1]*eref[2] - uref[1]*uref[2]/u2);
  
  diff[2][3] = c2*(    - eref[2]*eref[0] - uref[2]*uref[0]/u2);
  diff[2][4] = c2*(    - eref[2]*eref[1] - uref[2]*uref[1]/u2);
  diff[2][5] = c2*(1.0 - eref[2]*eref[2] - uref[2]*uref[2]/u2);
  
  double check1[3] = {0.0, 0.0, 0.0};
  double check2[3] = {0.0, 0.0, 0.0};
  for (int j = 0; j < 3; j++)
    for (int i = 0; i < 3; i++)
    {
      check1[j] += eref[i]*diff[i][j];
      check2[j] += eref[i]*diff[i][j+3];
    }
  if ( vecDotVec3d(check1,check1) > 1.0e-06)
    cout << "failed1: " << check1[0] << "," << check1[1] << "," << check1[2] << endl;
  if ( vecDotVec3d(check2,check2) > 1.0e-06)
    cout << "failed2: " << check2[0] << "," << check2[1] << "," << check2[2] << endl; 

  double GG = 0.0;
  double HH = 0.0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      {
	GG += diff[i][j]*diff[i][j];
	HH += diff[i][j+3]*diff[i][j+3];
      }
  double check3 = -0.5*eps/tke*(ae + au*tke/u2)*vecDotVec3d(eref,eref) + 0.5*(GG + HH);
  if (fabs(check3) > 1.0e-06)
    cout << "failed3: " << check3 << endl;

  double check4 = 0.0;
  double GH[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
    {
      GH[i][j] = 0.0;
      for (int k = 0; k < 3; k++)
	GH[i][j] += diff[k][i]*diff[k][j+3];
      check4 += GH[i][j];
    }
  if (fabs(check4) > 1.0e-06)
    cout << "failed4: " << check4 << endl;

	 
  // velocity
  diff[3][0] = c1;
  diff[4][1] = c1;
  diff[5][2] = c1;

}

double TurbModel_LANG::rhsDissipation(double(*Gn)[3], double dt)
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

  return rhs;
}

/*############################################################################################
#
# Cluster Simulation of SLM
#
############################################################################################*/

TurbModel_CLS::TurbModel_CLS(void) : TurbModel_SLM()
{
	cout << "Model         : CLS" << endl;

  nshells = 32;
  nmodes = 1000; //200;
  ncls = nshells*nmodes;

  e = new double[ncls][3];
  c = new double[ncls][6];

  Cn = 2.2;
  Cv = 1.0;
  Ceps1 = 1.5;
  Ceps2 = 11.0/6.0;
  Ceps3 = 0.01;
}

TurbModel_CLS::~TurbModel_CLS()
{
  delete[] e;
  delete[] c;
}

void TurbModel_CLS::initialHookScalarRansTurbModel(double Stau, double *struc,
			                                       double &eps_init, double *prod,
                                                   double *rRedi, double (*Gn)[3])
{
  srand(time(NULL));
  tke = 0.25;
  var.set_size(ncls,9);

  // Energy spectrum
  double mag_k[nshells], dk[nshells], Enukdk[nshells];
  double fkt = pow(6.0,0.2), Cnuk = 0.4843;

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

  int icls = 0;
  while (icls < ncls)
  {
    double u1 = (double)rand()/RAND_MAX*2.0 - 1.0;
    double u2 = (double)rand()/RAND_MAX*2.0 - 1.0;
  	if (u1*u1 + u2*u2 < 1.0)
  	{
  		var(icls,0) = 2.0*u1*sqrt(1.0 - u1*u1 - u2*u2);
  		var(icls,1) = 2.0*u2*sqrt(1.0 - u1*u1 - u2*u2);
  		var(icls,2) = 1.0 - 2.0*(u1*u1 + u2*u2);
  		icls += 1;
  	}
  }

  // Initialize cluster properties
  for (int ishells = 0; ishells < nshells; ishells++)
    for (int imodes = 0; imodes < nmodes; imodes++)
    {
      int icls = ishells*nmodes + imodes;

      // wave vector
      /*double theta;
      double r[3], s[3];

      double u1 = (double)rand()/RAND_MAX;
      double u2 = (double)rand()/RAND_MAX;
      theta = acos(2.0*u2 - 1.0);
      double phi = 2.0*M_PI*u1;

      e[icls][0] = sin(theta)*cos(phi);
      e[icls][1] = sin(theta)*sin(phi);
      e[icls][2] = cos(theta);*/

      /*double u = (double)rand()/RAND_MAX*2.0 - 1.0;
      double phi = 2.0*M_PI*(double)rand()/RAND_MAX;

      e[icls][0] = sqrt(1.0 - u*u)*cos(phi);
      e[icls][1] = sqrt(1.0 - u*u)*sin(phi);
      e[icls][2] = u; */


      // cluster tensor

      double umag2 = 2.0*nshells*Enukdk[ishells];
      var(icls,3) = 0.5*umag2*( 1.0 - var(icls,0)*var(icls,0) );
      var(icls,4) = 0.5*umag2*( 1.0 - var(icls,1)*var(icls,1) );
      var(icls,5) = 0.5*umag2*( 1.0 - var(icls,2)*var(icls,2) );
      var(icls,6) = 0.5*umag2*( -var(icls,0)*var(icls,1) );
      var(icls,7) = 0.5*umag2*( -var(icls,0)*var(icls,2) );
      var(icls,8) = 0.5*umag2*( -var(icls,1)*var(icls,2) );
    }

  outputs();
  eps = tke/Stau;

  // output to main
  for (int i = 0; i < 6; i++)
  {
    struc[i] = rey[i];
    struc[i+6] = dim[i];
  }
  eps_init = eps;

  int index[3][3] = { {0, 3, 4}, {3, 1, 5}, {4, 5, 2} };
  // production
  for (int i = 0; i < 3; i++)
    for (int j = i; j < 3; j++)
      for (int k = 0; k < 3; k++)
      {
      	prod[index[i][j]] = -rey[index[i][k]]*Gn[j][k] - rey[index[j][k]]*Gn[i][k];
      }

  // rapid redistribution
  for (int i = 0; i < 3; i++)
    for (int j = i; j < 3; j++)
    {
      rRedi[index[i][j]] = 0.0;
      for (int k = 0; k < 3; k++)
      	for (int l = 0; l < 3; l++)
      		rRedi[index[i][j]] += 2.0*Gn[k][l]*(M[i][l][j][k] + M[j][l][i][k]);
    }

  // output to screen
  cout << "Simulated TKE : " << tke << endl;
  cout << "Exact TKE     : " << 0.25 << endl;
  cout << "Eps           : " << eps << endl;
  cout << "Sk/eps        : " << tke/eps << endl;
  cout << "--------------------------------" << endl;
}

void TurbModel_CLS::fwEuler(double dt, double (*G)[3])
{
	mat rhs1(ncls,9);

  inputs(G);

	calcRhs(rhs1, var, dt);
	var += rhs1;

  double eps_rhs = calcRhsEps(G, dt);
  eps += dt*eps_rhs;

  correction();
	outputs();
}

void TurbModel_CLS::Heun(double dt, double (*G)[3])
{
	mat rhs1(ncls,9), rhs2(ncls,9);
	mat ivar(ncls,9);

  inputs(G);

	ivar = var;
	calcRhs(rhs1, ivar, dt);

	ivar = var + rhs1;
	calcRhs(rhs2, ivar, dt);

	var += 0.5*(rhs1 + rhs2);

  double eps_rhs = calcRhsEps(G, dt);
  eps += dt*eps_rhs;

  correction();
	outputs();
}

void TurbModel_CLS::CrankN(double dt, double (*Gn)[3])
{/*
  for (int icls = 0; icls < ncls; icls++)
  {
    double *eicls = e[icls];
    double *cicls = c[icls];

    double de[3], dc[6], Ae[3][3], Ac[6][6];

    // right and left hand side
    calcRhs(de, dc, eicls, cicls, Gn);
    calcRhsJacob(Ae, Ac, eicls, cicls, Gn, dt);

    // solve linear systems
    GaussSolver(3, Ae[0], de);
    GaussSolver(6, Ac[0], dc);

    // update solutions
    for (int i = 0; i < 3; i++) eicls[i] += de[i];
    for (int i = 0; i < 6; i++) cicls[i] += dc[i];

    normVec3d(eicls);
    SymMatOrthoVec3d(cicls, eicls);
  }

  double eps_rhs = calcRhsEps(Gn, dt);
  eps += dt*eps_rhs;

  calcTurbStatistics();
  calcTurbTimeScale();*/
}

void TurbModel_CLS::calcRhs(mat &rhs, mat &var, double dt)
{
  int iIndex[6] = {0, 1, 2, 0, 0, 1};
  int jIndex[6] = {0, 1, 2, 1, 2, 2};

	for (int icls = 0; icls < ncls; icls++)
	{
		for (int i = 0; i < 9; i++) rhs(icls,i) = 0.0;

		double eref[3], cref[3][3];

		eref[0] = var(icls,0); eref[1] = var(icls,1); eref[2] = var(icls,2);

		cref[0][0] = var(icls,3);   cref[0][1] = var(icls,6);   cref[0][2] = var(icls,7);
		cref[1][0] = var(icls,6);   cref[1][1] = var(icls,4);   cref[1][2] = var(icls,8);
		cref[2][0] = var(icls,7);   cref[2][1] = var(icls,8);   cref[2][2] = var(icls,5);

		double Ge[3];
		double eGe;
		double Gc[3][3], eGce1[3][3], eGce2[3][3];

		matTransDotVec3d(Ge, eref, Gn);
		eGe = vecDotMatDotVec3d(eref, eref, Gn);
		matTimesMat3d(Gc, Gv, cref);
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
			{
				eGce1[i][j] = 0.0;
				eGce2[i][j] = 0.0;
				for (int q = 0; q < 3; q++)
					for (int s = 0; s < 3; s++)
					{
						eGce1[i][j] += eref[q]*Gn[q][s]*cref[s][i]*eref[j];
						eGce2[i][j] += eref[q]*Gv[q][s]*cref[s][i]*eref[j];
					}
			}

		for (int i = 0; i < 3; i++)
			rhs(icls,i) += -Ge[i] + eGe*eref[i];

		for (int n = 0; n < 6; n++)
		{
			int i = iIndex[n];
			int j = jIndex[n];
			rhs(icls,n+3) += -Gc[i][j] - Gc[j][i] + eGce1[i][j] + eGce1[j][i]
                                            + eGce2[i][j] + eGce2[j][i];
		}

		// slow rotational randomization
		double C1_icls = 0.0;
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				C1_icls += C1[i][j]*eref[i]*eref[j];

		double psikk = cref[0][0] + cref[1][1] + cref[2][2];
		rhs(icls,3) -= 2.0*C1_icls*cref[0][0] - C1_icls*psikk + C1_icls*psikk*eref[0]*eref[0];
		rhs(icls,4) -= 2.0*C1_icls*cref[1][1] - C1_icls*psikk + C1_icls*psikk*eref[1]*eref[1];
		rhs(icls,5) -= 2.0*C1_icls*cref[2][2] - C1_icls*psikk + C1_icls*psikk*eref[2]*eref[2];
		rhs(icls,6) -= 2.0*C1_icls*cref[0][1] + C1_icls*psikk*eref[0]*eref[1];
		rhs(icls,7) -= 2.0*C1_icls*cref[0][2] + C1_icls*psikk*eref[0]*eref[2];
		rhs(icls,8) -= 2.0*C1_icls*cref[1][2] + C1_icls*psikk*eref[1]*eref[2];

		for (int i = 0; i < 9; i++) rhs(icls,i) *= dt;
	}
}

void TurbModel_CLS::calcRhsJacob(double (*eDrift)[3], double (*cDrift)[6], double *eref,
		                           double *cref, double (*G)[3], double dt)
{/*
  int iIndex[6] = {0, 1, 2, 0, 0, 1};
  int jIndex[6] = {0, 1, 2, 1, 2, 2};
	double delta[3][3] = { {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0} };

	// -----------------------------------------------------------------------------------------
	// expanded gradients
	// -----------------------------------------------------------------------------------------
  double r[3][3], d[3][3], rd[3][3];
  double q2 = 2.0*tke;
  double Gn[3][3], Gv[3][3];

  r[0][0] = rey[0]/q2;   r[0][1] = rey[3]/q2;   r[0][2] = rey[4]/q2;
  r[1][0] = rey[3]/q2;   r[1][1] = rey[1]/q2;   r[1][2] = rey[5]/q2;
  r[2][0] = rey[4]/q2;   r[2][1] = rey[5]/q2;   r[2][2] = rey[2]/q2;

  d[0][0] = dim[0]/q2;   d[0][1] = dim[3]/q2;   d[0][2] = dim[4]/q2;
  d[1][0] = dim[3]/q2;   d[1][1] = dim[1]/q2;   d[1][2] = dim[5]/q2;
  d[2][0] = dim[4]/q2;   d[2][1] = dim[5]/q2;   d[2][2] = dim[2]/q2;

  matTimesMat3d(rd, r, d);

  for (int i = 0; i < 3; i++)
  	for (int j = 0; j < 3; j++)
  	{
  		Gn[i][j] = G[i][j] + Cn/tau*rd[i][j];
  		Gv[i][j] = G[i][j] + Cv/tau*rd[i][j];
  	}

	// -----------------------------------------------------------------------------------------
  // slow rotational randomization
	// -----------------------------------------------------------------------------------------
  double f[3][3];
  f[0][0] = cir[0]/q2;   f[0][1] = cir[3]/q2;   f[0][2] = cir[4]/q2;
  f[1][0] = cir[3]/q2;   f[1][1] = cir[1]/q2;   f[1][2] = cir[5]/q2;
  f[2][0] = cir[4]/q2;   f[2][1] = cir[5]/q2;   f[2][2] = cir[2]/q2;

  double fnn = vecDotMatDotVec3d(eref, eref, f);

  double Ovec[3], Omag, C1;
  Ovec[0] = rd[2][1] - rd[1][2];
  Ovec[1] = rd[0][2] - rd[2][0];
  Ovec[2] = rd[1][0] - rd[0][1];

  Omag = sqrt(vecDotVec3d(Ovec,Ovec));

  C1 = 8.5/tau*Omag*fnn;

	// -----------------------------------------------------------------------------------------
	// Form Jacobian for wave vector drift
	// -----------------------------------------------------------------------------------------
	double Gn_e[3], e_Gn[3], e_Gv[3], e_Gn_e = 0.0;
	for (int i = 0; i < 3; i++)
	{
		Gn_e[i] = 0.0;
		e_Gn[i] = 0.0;
		e_Gv[i] = 0.0;
		for (int j = 0; j < 3; j++)
		{
			Gn_e[i] += Gn[i][j]*eref[j];
			e_Gn[i] += eref[j]*Gn[j][i];
			e_Gv[i] += eref[j]*Gv[j][i];
		}
		e_Gn_e += eref[i]*Gn_e[i];
	}

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			eDrift[i][j] = delta[i][j]/dt
			             - 0.5*(-Gn[j][i] + Gn_e[j]*eref[i] + e_Gn[j]*eref[i] + delta[i][j]*e_Gn_e);

	// -----------------------------------------------------------------------------------------
	// Form Jacobian for cluster drift
	// -----------------------------------------------------------------------------------------
	// expanded gradient
	double b[3];
	for (int i = 0; i < 3; i++)
		b[i] = e_Gn[i] + e_Gv[i];

  cDrift[0][0] = 2.0*Gv[0][0] - 2.0*b[0]*eref[0];
  cDrift[0][1] = 0.0;
  cDrift[0][2] = 0.0;
  cDrift[0][3] = 2.0*Gv[0][1] - 2.0*b[1]*eref[0];
  cDrift[0][4] = 2.0*Gv[0][2] - 2.0*b[2]*eref[0];
  cDrift[0][5] = 0.0;

  cDrift[1][0] = 0.0;
  cDrift[1][1] = 2.0*Gv[1][1] - 2.0*b[1]*eref[1];
  cDrift[1][2] = 0.0;
  cDrift[1][3] = 2.0*Gv[1][0] - 2.0*b[0]*eref[1];
  cDrift[1][4] = 0.0;
  cDrift[1][5] = 2.0*Gv[1][2] - 2.0*b[2]*eref[1];

  cDrift[2][0] = 0.0;
  cDrift[2][1] = 0.0;
  cDrift[2][2] = 2.0*Gv[2][2] - 2.0*b[2]*eref[2];
  cDrift[2][3] = 0.0;
  cDrift[2][4] = 2.0*Gv[2][0] - 2.0*b[0]*eref[2];
  cDrift[2][5] = 2.0*Gv[2][1] - 2.0*b[1]*eref[2];

  cDrift[3][0] = Gv[1][0] - b[0]*eref[1];
  cDrift[3][1] = Gv[0][1] - b[1]*eref[0];
  cDrift[3][2] = 0.0;
  cDrift[3][3] = Gv[0][0] - b[0]*eref[0] + Gv[1][1] - b[1]*eref[1];
  cDrift[3][4] = Gv[1][2] - b[2]*eref[1];
  cDrift[3][5] = Gv[0][2] - b[2]*eref[0];

  cDrift[4][0] = Gv[2][0] - b[0]*eref[2];
  cDrift[4][1] = 0.0;
  cDrift[4][2] = Gv[0][2] - b[2]*eref[0];
  cDrift[4][3] = Gv[2][1] - b[1]*eref[2];
  cDrift[4][4] = Gv[0][0] - b[0]*eref[0] + Gv[2][2] - b[2]*eref[2];
  cDrift[4][5] = Gv[0][1] - b[1]*eref[0];

  cDrift[5][0] = 0.0;
  cDrift[5][1] = Gv[2][1] - b[1]*eref[2];
  cDrift[5][2] = Gv[1][2] - b[2]*eref[1];
  cDrift[5][3] = Gv[2][0] - b[0]*eref[2];
  cDrift[5][4] = Gv[1][0] - b[0]*eref[1];
  cDrift[5][5] = Gv[1][1] - b[1]*eref[1] + Gv[2][2] - b[2]*eref[2];

  // slow rotational randomization model
  for (int n = 0; n < 6; n++)
  {
  	cDrift[n][n] -= 2.0*C1;
  	for (int m = 0; m < 3; m++)
  	{
  		int i = iIndex[n];
  		int j = jIndex[n];
  		cDrift[n][m] += C1*(delta[i][j] - eref[i]*eref[j]);
  	}
  }

  // form jacobian
  for (int n = 0; n < 6; n++)
  	for (int m = 0; m < 6; m++)
  		cDrift[n][m] *= -0.5;

  for (int n = 0; n < 6; n++)
  	cDrift[n][n] += 1.0/dt; */

}

double TurbModel_CLS::calcRhsEps(double(*Gn)[3], double dt)
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

  return rhs;
}

void TurbModel_CLS::correction()
{
	for (int icls = 0; icls < ncls; icls++)
	{
		double eicls[3], cicls[6];

		eicls[0] = var(icls,0); eicls[1] = var(icls,1); eicls[2] = var(icls,2);

		cicls[0] = var(icls,3); cicls[1] = var(icls,4); cicls[2] = var(icls,5);
    cicls[3] = var(icls,6); cicls[4] = var(icls,7); cicls[5] = var(icls,8);

		normVec3d(eicls);
		SymMatOrthoVec3d(cicls, eicls);

		var(icls,0) = eicls[0]; var(icls,1) = eicls[1]; var(icls,2) = eicls[2];

		var(icls,3) = cicls[0]; var(icls,4) = cicls[1]; var(icls,5) = cicls[2];
		var(icls,6) = cicls[3]; var(icls,7) = cicls[4]; var(icls,8) = cicls[5];
	}
}

void TurbModel_CLS::inputs(double(*G)[3])
{
  double r[3][3], d[3][3], f[3][3], rd[3][3];
  double q2 = 2.0*tke;

  r[0][0] = rey[0]/q2;   r[0][1] = rey[3]/q2;   r[0][2] = rey[4]/q2;
  r[1][0] = rey[3]/q2;   r[1][1] = rey[1]/q2;   r[1][2] = rey[5]/q2;
  r[2][0] = rey[4]/q2;   r[2][1] = rey[5]/q2;   r[2][2] = rey[2]/q2;

  d[0][0] = dim[0]/q2;   d[0][1] = dim[3]/q2;   d[0][2] = dim[4]/q2;
  d[1][0] = dim[3]/q2;   d[1][1] = dim[1]/q2;   d[1][2] = dim[5]/q2;
  d[2][0] = dim[4]/q2;   d[2][1] = dim[5]/q2;   d[2][2] = dim[2]/q2;

  f[0][0] = cir[0]/q2;   f[0][1] = cir[3]/q2;   f[0][2] = cir[4]/q2;
  f[1][0] = cir[3]/q2;   f[1][1] = cir[1]/q2;   f[1][2] = cir[5]/q2;
  f[2][0] = cir[4]/q2;   f[2][1] = cir[5]/q2;   f[2][2] = cir[2]/q2;

  matTimesMat3d(rd, r, d);
  double tau = 0.0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      tau += rd[i][j]*r[j][i];
  tau *= 2.0*Cv*tke/eps;

  for (int i = 0; i < 3; i++)
  	for (int j = 0; j < 3; j++)
  	{
  		Gn[i][j] = G[i][j] + Cn/tau*rd[i][j];
  		Gv[i][j] = G[i][j] + Cv/tau*rd[i][j];
  	}

  double Ovec[3], Omag;
  Ovec[0] = rd[2][1] - rd[1][2];
  Ovec[1] = rd[0][2] - rd[2][0];
  Ovec[2] = rd[1][0] - rd[0][1];

  Omag = sqrt(vecDotVec3d(Ovec,Ovec));

  for (int i = 0; i < 3; i++)
  	for (int j = 0; j < 3; j++)
  		C1[i][j] = 8.5/tau*Omag*f[i][j];
}

void TurbModel_CLS:: outputs()
{
	int index[3][3] = { {0, 3, 4}, {3, 1, 5}, {4, 5, 2} };

	// initialize to zero
  for (int i = 0; i < 6; i++)
    rey[i] = 0.0;

  for (int i = 0; i < 3; i++)
    for (int j = i; j < 3; j++)
    {
      e2[i][j] = 0.0;
      for (int k = 0; k < 3; k++)
      	for (int l = k; l < 3; l++)
      	{
      		M[i][j][k][l] = 0.0;
      		e4[i][j][k][l] = 0.0;
      	}
    }

  // structure tensors
  for (int icls = 0; icls < ncls; icls ++)
  {
  	for (int i = 0; i < 6; i++)
  		rey[i] += var(icls,i+3);

    for (int i = 0; i < 3; i++)
      for(int j = i; j <3; j++)
      {
      	e2[i][j] += var(icls,i)*var(icls,j);
      	for (int k = 0; k< 3; k++)
      		for (int l = k; l < 3; l++)
      		{
      		  M[i][j][k][l] += var(icls,index[i][j]+3)*var(icls,k)*var(icls,l);
      			e4[i][j][k][l] += var(icls,i)*var(icls,j)*var(icls,k)*var(icls,l);
      		}
      }
  }

  // normalize
  for (int i = 0; i < 6; i++)
    rey[i] /= (ncls);

  for (int i = 0; i < 3; i++)
    for(int j = i; j <3; j++)
    {
      e2[i][j] /= (ncls);
      for (int k = 0; k< 3; k++)
      	for (int l = k;l < 3; l++)
      	{
      		M[i][j][k][l] /= (ncls);
      		e4[i][j][k][l] /= (ncls);
      	}
    }

  // make symmetric
  for (int i = 0; i < 3; i++)
    for (int j = i; j < 3; j++)
    {
      e2[i][j] = e2[j][i];

      M[i][j][1][0] = M[i][j][0][1];
      M[i][j][2][0] = M[i][j][0][2];
      M[i][j][2][1] = M[i][j][1][2];

      e4[i][j][1][0] = e4[i][j][0][1];
      e4[i][j][2][0] = e4[i][j][0][2];
      e4[i][j][2][1] = e4[i][j][1][2];
    }

  for (int k = 0; k < 3; k++)
    for (int l = 0; l < 3; l++)
    {
      M[1][0][k][l] = M[0][1][k][l];
      M[2][0][k][l] = M[0][2][k][l];
      M[2][1][k][l] = M[1][2][k][l];

      e4[1][0][k][l] = e4[0][1][k][l];
      e4[2][0][k][l] = e4[0][2][k][l];
      e4[2][1][k][l] = e4[1][2][k][l];
    }

  // tke, dimensionality and circulicity
  tke = 0.5*(rey[0] + rey[1] + rey[2]);

  for (int i = 0; i < 3; i++)
  	for (int j = i; j < 3; j++)
  		dim[index[i][j]] = M[0][0][i][j] + M[1][1][i][j] + M[2][2][i][j];

  for (int i = 0; i < 3; i++) cir[i] = 2.0*tke - rey[i] - dim[i];
  for (int i = 3; i < 6; i++) cir[i] = -rey[i] - dim[i];
}

void TurbModel_CLS::calcReStress(double *struc, double &eps_main, double *prod,
                                 double *rRedi, double *sRediEps,
                                 double (*Gn)[3], double (*Gnph)[3],
                                 double (*Gnp1)[3], double dt, char* tIntName)
{
  // time integration
  if      (strcmp(tIntName, "fwEuler") == 0) fwEuler(dt, Gn);
  else if (strcmp(tIntName, "Heun") == 0)    Heun(dt, Gn);
  else if (strcmp(tIntName, "CrankN") == 0)  CrankN(dt, Gn);
  else    cout << "Not a valid time integration scheme" << endl;

  // structure tensors and eps to main
  for (int i = 0; i < 6; i++)
  {
    struc[i] = rey[i];
    struc[i+6] = dim[i];
  }
  eps_main = eps;

  int index[3][3] = { {0, 3, 4}, {3, 1, 5}, {4, 5, 2} };
  // production
  for (int i = 0; i < 3; i++)
    for (int j = i; j < 3; j++)
      for (int k = 0; k < 3; k++)
      {
      	prod[index[i][j]] = -rey[index[i][k]]*Gn[j][k] - rey[index[j][k]]*Gn[i][k];
      }

  // rapid redistribution
  for (int i = 0; i < 3; i++)
    for (int j = i; j < 3; j++)
    {
    	rRedi[index[i][j]] = 0.0;
    	for (int k = 0; k < 3; k++)
    		for (int l = 0; l < 3; l++)
    			rRedi[index[i][j]] += 2.0*Gn[k][l]*(M[i][l][j][k] + M[j][l][i][k]);
    }

  // slow redistribution + dissipation
  double RD[3][3] = { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0} };
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
      	RD[i][j] += rey[index[i][k]]*dim[index[k][j]];

  for (int i = 0; i < 3; i++)
    for (int j = i; j < 3; j++)
    {
      sRediEps[index[i][j]] = 0.0;
    }

}

/*############################################################################################
#
# Eulerian Simulation of SLM
#
############################################################################################*/

TurbModel_EUL::TurbModel_EUL(void) : TurbModel_SLM()
{
	cout << "Model         : EUL" << endl;

	epsilon = 6.0;

  Cn = 2.2;
  Cv = 1.0;
  Ceps1 = 1.5;
  Ceps2 = 11.0/6.0;
  Ceps3 = 0.01;
}

TurbModel_EUL::~TurbModel_EUL()
{
	delete [] grid;
	delete [] elem;
	delete [] lambda;
	delete [] theta;
}

void TurbModel_EUL::initialHookScalarRansTurbModel(double Stau, double *struc,
			                                             double &eps_init, double *prod,
                                                   double *rRedi, double (*G)[3])
{
  tke = 0.25;

  // Open the mesh
	string fileGrid = "grid.dat";
	FILE *fp;
	if ((fp = fopen(fileGrid.c_str(), "rt")) == NULL)
		cout << "Could not open the file "<< fileGrid << endl;

	int count = 0;
	count += fscanf(fp,"variables = \"x\", \"y\", \"z\"\n");
	count += fscanf(fp,"zone T=\"grid\"\n");
	count += fscanf(fp,"Nodes=%d, Elements=%d, ZONETYPE=FETriangle\n", &nnodes, &nelems);
	count += fscanf(fp,"DATAPACKING=POINT\n");
	count += fscanf(fp,"DT=(SINGLE SINGLE SINGLE)\n");
	if (count != 2) cout << "Error: nnodes or/and nelems were not read" << endl;
	count = 0;

	cout << "Nodes         : " << nnodes << endl;
	cout << "Elems         : " << nelems << endl;

	// Initialize variables
  grid   = new double [nnodes][3];
  elem   = new int [nelems][3];
  lambda = new double [nnodes];
  theta  = new double [nnodes];

	Ainv.set_size(nnodes,nnodes);
	U.set_size(nnodes,nnodes);
	V.set_size(nnodes,nnodes);
	W.set_size(nnodes,nnodes,6);
	J.set_size(nnodes,3,3);
	K.set_size(nnodes,3,3);
	H.set_size(nnodes,3,3);

	var.set_size(nnodes,6);

	for (int i = 0; i < nnodes; i++)
	{
		// Read mesh
		for (int j = 0; j < 3; j++)
			count += fscanf(fp, "%lf", &grid[i][j]);

		// Compute angles
  	lambda[i] = atan2(grid[i][1],grid[i][0]);
  	theta[i]  = atan2(grid[i][2],sqrt(grid[i][0]*grid[i][0] + grid[i][1]*grid[i][1]));

  	double sinLa = sin(lambda[i]);
    double cosLa = cos(lambda[i]);
    double sinTh = sin(theta[i]);
    double cosTh = cos(theta[i]);

    // Compute J, K, H
    J(i,0,0) = -sinLa*cosLa;
    J(i,0,1) = cosLa*cosLa;
    J(i,0,2) = 0.0;

    J(i,1,0) = -sinLa*sinLa;
    J(i,1,1) = sinLa*cosLa;
    J(i,1,2) = 0.0;

    J(i,2,0) = -sinLa*sinTh/cosTh;
    J(i,2,1) = cosLa*sinTh/cosTh;
    J(i,2,2) = 0.0;

    K(i,0,0) = -cosLa*cosLa*sinTh*cosTh;
    K(i,0,1) = -sinLa*cosLa*sinTh*cosTh;
    K(i,0,2) = cosLa*cosTh*cosTh;

    K(i,1,0) = -sinLa*cosLa*sinTh*cosTh;
    K(i,1,1) = -sinLa*sinLa*sinTh*cosTh;
    K(i,1,2) = sinLa*cosTh*cosTh;

    K(i,2,0) = -cosLa*sinTh*sinTh;
    K(i,2,1) = -sinLa*sinTh*sinTh;
    K(i,2,2) = sinTh*cosTh;

    H(i,0,0) = cosLa*cosLa*cosTh*cosTh;
    H(i,0,1) = sinLa*cosLa*cosTh*cosTh;
    H(i,0,2) = cosLa*sinTh*cosTh;

    H(i,1,0) = sinLa*cosLa*cosTh*cosTh;
    H(i,1,1) = sinLa*sinLa*cosTh*cosTh;
    H(i,1,2) = sinLa*sinTh*cosTh;

    H(i,2,0) = cosLa*sinTh*cosTh;
    H(i,2,1) = sinLa*sinTh*cosTh;
    H(i,2,2) = sinTh*sinTh;

    // Initial condition
    var(i,0) = tke/(4.0*M_PI)*(1.0 - H(i,0,0));
    var(i,1) = tke/(4.0*M_PI)*(1.0 - H(i,1,1));
    var(i,2) = tke/(4.0*M_PI)*(1.0 - H(i,2,2));
    var(i,3) = tke/(4.0*M_PI)*(-H(i,0,1));
    var(i,4) = tke/(4.0*M_PI)*(-H(i,0,2));
    var(i,5) = tke/(4.0*M_PI)*(-H(i,1,2));
	}

	for (int i = 0; i < nelems; i++)
		for (int j = 0; j < 3; j++)
			count += fscanf(fp, "%d", &elem[i][j]);

	if ( count != 3*(nnodes + nelems) )
		cout << "Error: incorrect number of nodes or elems read" << endl;

	fclose(fp);
	cout << "Read grid file: " << fileGrid << endl;

  // Compute the metrics
	for (int i = 0; i < nnodes; i++)
		for (int j = 0; j < nnodes; j++)
		{
			double sinLa_i = sin(lambda[i]);   double sinLa_j = sin(lambda[j]);
			double cosLa_i = cos(lambda[i]);   double cosLa_j = cos(lambda[j]);
			double sinTh_i = sin(theta[i]);    double sinTh_j = sin(theta[j]);
			double cosTh_i = cos(theta[i]);    double cosTh_j = cos(theta[j]);

			double sinLa_ij = sinLa_i*cosLa_j - cosLa_i*sinLa_j;
			double cosLa_ij = sinLa_i*sinLa_j + cosLa_i*cosLa_j;

			double r2 = 2.0 - 2.0*( cosTh_i*cosTh_j*cosLa_ij + sinTh_i*sinTh_j );

			Ainv(i,j) = exp( -epsilon*epsilon*r2 );

			U(i,j) = cosTh_i*cosTh_j*sinLa_ij;
			V(i,j) = sinTh_i*cosTh_j*cosLa_ij - cosTh_i*sinTh_j;

			U(i,j) *= -2.0*epsilon*epsilon*exp( -epsilon*epsilon*r2 );
			V(i,j) *= -2.0*epsilon*epsilon*exp( -epsilon*epsilon*r2 );
		}

	Ainv = inv(Ainv);

	U = U*Ainv;
	V = V*Ainv;
  //coeff.print("coeff:");

	for (int inode = 0; inode < nnodes; inode++)
		for (int jnode = 0; jnode < nnodes; jnode++)
		{
				W(inode,jnode,0) = Ainv(inode,jnode)*H(jnode,0,0);
				W(inode,jnode,1) = Ainv(inode,jnode)*H(jnode,1,1);
				W(inode,jnode,2) = Ainv(inode,jnode)*H(jnode,2,2);
				W(inode,jnode,3) = Ainv(inode,jnode)*H(jnode,0,1);
				W(inode,jnode,4) = Ainv(inode,jnode)*H(jnode,0,2);
				W(inode,jnode,5) = Ainv(inode,jnode)*H(jnode,1,2);
		}

  outputs();
  eps = tke/Stau;

  // output to main
  for (int i = 0; i < 6; i++)
  {
    struc[i] = rey[i];
    struc[i+6] = dim[i];
  }
  eps_init = eps;

  int index[3][3] = { {0, 3, 4}, {3, 1, 5}, {4, 5, 2} };
  // production
  for (int i = 0; i < 3; i++)
    for (int j = i; j < 3; j++)
      for (int k = 0; k < 3; k++)
      {
      	prod[index[i][j]] = -rey[index[i][k]]*G[j][k] - rey[index[j][k]]*G[i][k];
      }

  // rapid redistribution
  for (int i = 0; i < 3; i++)
    for (int j = i; j < 3; j++)
    {
      rRedi[index[i][j]] = 0.0;
      for (int k = 0; k < 3; k++)
      	for (int l = 0; l < 3; l++)
      		rRedi[index[i][j]] += 2.0*G[k][l]*(M[i][l][j][k] + M[j][l][i][k]);
    }

  // output to screen
  cout << "Simulated TKE : " << tke << endl;
  cout << "Exact TKE     : " << 0.25 << endl;
  cout << "Eps           : " << eps << endl;
  cout << "Sk/eps        : " << tke/eps << endl;
  cout << "--------------------------------" << endl;
}

void TurbModel_EUL::fwEuler(double dt, double (*G)[3])
{
	mat rhs(nnodes,6);

  inputs(G);

	calcRhs(rhs, var, dt);
	var += rhs;

  double eps_rhs = calcRhsEps(G, dt);
  eps += dt*eps_rhs;

	outputs();
}

void TurbModel_EUL::Heun(double dt, double (*G)[3])
{

}

void TurbModel_EUL::CrankN(double dt, double (*G)[3])
{

}

void TurbModel_EUL::RK4(double dt, double (*G)[3])
{
	mat rhs1(nnodes,6), rhs2(nnodes,6), rhs3(nnodes,6), rhs4(nnodes,6);
	mat ivar(nnodes,6);

  inputs(G);

	ivar = var;
	calcRhs(rhs1, ivar, dt);

	ivar = var + 0.5*rhs1;
	calcRhs(rhs2, ivar, dt);

	ivar = var + 0.5*rhs2;
	calcRhs(rhs3, ivar, dt);

	ivar = var + rhs3;
	calcRhs(rhs4, ivar, dt);

	var += 1.0/6.0*(rhs1 + 2.0*rhs2 + 2.0*rhs3 + rhs4);

  double eps_rhs = calcRhsEps(G, dt);
  eps += dt*eps_rhs;

	outputs();
}

void TurbModel_EUL::calcRhs(mat &rhs, mat &var, double dt)
{
	mat conv1;  conv1 = U*var;
	mat conv2;  conv2 = V*var;

	int index[3][3] = { {0, 3, 4}, {3, 1, 5}, {4, 5, 2} };
  int iIndex[6] = {0, 1, 2, 0, 0, 1};
  int jIndex[6] = {0, 1, 2, 1, 2, 2};
  double delta[3][3] = { {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0} };

	// gradients
  double Gnkk = Gn[0][0] + Gn[1][1] + Gn[2][2];
  double Gt[3][3];
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			Gt[i][j] = Gn[i][j] + Gv[i][j];

	// loop through nodes
	for (int inode = 0; inode < nnodes; inode++)
	{
		// coefficients
		double GJ = 0.0, GK = 0.0, GH = 0.0, CH = 0.0;
		double GtH[3][3];

		for (int q = 0; q < 3; q++)
			for (int r = 0; r < 3; r++)
			{
				GJ += Gn[q][r]*J(inode,q,r);
				GK += Gn[q][r]*K(inode,q,r);
				GH += Gn[q][r]*H(inode,q,r);
				CH += C1[q][r]*H(inode,q,r);

				GtH[q][r] = 0.0;
				for (int p = 0; p < 3; p++)
					GtH[q][r] += Gt[p][q]*H(inode,p,r);
			}

		double GGtH[3][3];
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				GGtH[i][j] = -Gv[i][j] + GtH[j][i];

		double I1[6], I2[6], I3[6];
		for (int n = 0; n < 6; n++)
		{
			int i = iIndex[n];
			int j = jIndex[n];

			// convection
			conv1(inode,n) *= GJ;
			conv2(inode,n) *= GK;

			// sources (L1)
			I1[n] = (Gnkk - 3.0*GH - 2.0*CH)*var(inode,n);

			// sources (L2)
			I2[n] = 0.0;
			for (int s = 0; s < 3; s++)
				I2[n] += GGtH[j][s]*var(inode,index[i][s]) + GGtH[i][s]*var(inode,index[j][s]);

			// sources (L3)
			I3[n] = CH*(delta[i][j] - H(inode,i,j))*(var(inode,0) + var(inode,1) + var(inode,2));

			// right hand side
			rhs(inode,n) = (conv1(inode,n) + conv2(inode,n) + I1[n] + I2[n] + I3[n])*dt;
		}
	}
}

double TurbModel_EUL::calcRhsEps(double(*G)[3], double dt)
{
  double prod = -rey[0]*G[0][0] - rey[3]*G[0][1] - rey[4]*G[0][2]
                -rey[3]*G[1][0] - rey[1]*G[1][1] - rey[5]*G[1][2]
                -rey[4]*G[2][0] - rey[5]*G[2][1] - rey[2]*G[2][2];

  double vort[3];
  vort[0] = G[2][1] - G[1][2];
  vort[1] = G[0][2] - G[2][0];
  vort[2] = G[1][0] - G[0][1];

  double d[3][3];
  double q2 = 2.0*tke;
  d[0][0] = dim[0]/q2;   d[0][1] = dim[3]/q2;   d[0][2] = dim[4]/q2;
  d[1][0] = dim[3]/q2;   d[1][1] = dim[1]/q2;   d[1][2] = dim[5]/q2;
  d[2][0] = dim[4]/q2;   d[2][1] = dim[5]/q2;   d[2][2] = dim[2]/q2;

  double OdO = vecDotMatDotVec3d(vort, vort, d);

  double rhs = Ceps1*eps/tke*prod - Ceps2*eps*eps/tke - Ceps3*eps*sqrt(OdO);

  return rhs;
}

void TurbModel_EUL::inputs(double(*G)[3])
{
  double r[3][3], d[3][3], f[3][3], rd[3][3];
  double q2 = 2.0*tke;

  r[0][0] = rey[0]/q2;   r[0][1] = rey[3]/q2;   r[0][2] = rey[4]/q2;
  r[1][0] = rey[3]/q2;   r[1][1] = rey[1]/q2;   r[1][2] = rey[5]/q2;
  r[2][0] = rey[4]/q2;   r[2][1] = rey[5]/q2;   r[2][2] = rey[2]/q2;

  d[0][0] = dim[0]/q2;   d[0][1] = dim[3]/q2;   d[0][2] = dim[4]/q2;
  d[1][0] = dim[3]/q2;   d[1][1] = dim[1]/q2;   d[1][2] = dim[5]/q2;
  d[2][0] = dim[4]/q2;   d[2][1] = dim[5]/q2;   d[2][2] = dim[2]/q2;

  f[0][0] = cir[0]/q2;   f[0][1] = cir[3]/q2;   f[0][2] = cir[4]/q2;
  f[1][0] = cir[3]/q2;   f[1][1] = cir[1]/q2;   f[1][2] = cir[5]/q2;
  f[2][0] = cir[4]/q2;   f[2][1] = cir[5]/q2;   f[2][2] = cir[2]/q2;

  matTimesMat3d(rd, r, d);
  double tau = 0.0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      tau += rd[i][j]*r[j][i];
  tau *= 2.0*Cv*tke/eps;

  for (int i = 0; i < 3; i++)
  	for (int j = 0; j < 3; j++)
  	{
  		Gn[i][j] = G[i][j] + Cn/tau*rd[i][j];
  		Gv[i][j] = G[i][j] + Cv/tau*rd[i][j];
  	}

  double Ovec[3], Omag;
  Ovec[0] = rd[2][1] - rd[1][2];
  Ovec[1] = rd[0][2] - rd[2][0];
  Ovec[2] = rd[1][0] - rd[0][1];

  Omag = sqrt(vecDotVec3d(Ovec,Ovec));

  for (int i = 0; i < 3; i++)
  	for (int j = 0; j < 3; j++)
  		C1[i][j] = 8.5/tau*Omag*f[i][j];
}

void TurbModel_EUL:: outputs()
{
	int index[3][3] = { {0, 3, 4}, {3, 1, 5}, {4, 5, 2} };

	// initialize to zero
  for (int i = 0; i < 6; i++)
  {
    rey[i] = 0.0;
    dim[i] = 0.0;
  }

  // structure tensors
  mat c(nnodes,6), d(nnodes,6);

  c = Ainv*var;
  for (int i = 0; i < 6; i++)
  	d.col(i) = W.slice(i)*(var.col(0) + var.col(1) + var.col(2));

  for (int inode = 0; inode < nnodes; inode++)
  {
  	for (int i = 0; i < 6; i++)
  	{
  		rey[i] += c(inode,i);
  		dim[i] += d(inode,i);
  	}
  }

  // normalize
  double I1 = M_PI/(epsilon*epsilon)*(1.0 - exp(-4.0*epsilon*epsilon));
  for (int i = 0; i < 6; i++)
  {
    rey[i] *= I1;
    dim[i] *= I1;
  }

  // tke, dimensionality and circulicity
  tke = 0.5*(rey[0] + rey[1] + rey[2]);

  for (int i = 0; i < 3; i++) cir[i] = 2.0*tke - rey[i] - dim[i];
  for (int i = 3; i < 6; i++) cir[i] = -rey[i] - dim[i];
}

void TurbModel_EUL::writeData(double St)
{
	FILE *fsph;
	string fname;

	// sphere data
	fname = "sphere." + to_string(St) + ".dat";
  fsph = fopen(fname.c_str(), "w");
  if (fsph != NULL)
  {
  	fprintf(fsph,"variables = \"x\", \"y\", \"z\" \"coeff_1\" \"coeff_2\" \"coeff_3\"");
  	fprintf(fsph,"\"coeff_4\" \"coeff_5\" \"coeff_6\"\n");
  	fprintf(fsph,"zone T=\"sphere\"\n");
  	fprintf(fsph,"Nodes=%d, Elements=%d, ZONETYPE=FETriangle\n", nnodes, nelems);
  	fprintf(fsph,"DATAPACKING=POINT\n");
  	fprintf(fsph,"DT=(SINGLE SINGLE SINGLE SINGLE)\n");

    for (int i = 0; i < nnodes; i++)
    {
      fprintf(fsph, "%.8e\t%.8e\t%.8e\t", grid[i][0], grid[i][1], grid[i][2]);
      fprintf(fsph, "%.8e\t%.8e\t%.8e\t", var(i,0), var(i,1), var(i,2));
      fprintf(fsph, "%.8e\t%.8e\t%.8e\n", var(i,3), var(i,4), var(i,5));
    }
    for (int i = 0; i < nelems; i++)
    	fprintf(fsph, "%d\t%d\t%d\n", elem[i][0], elem[i][1], elem[i][2]);
  }
  else
  	cout << "Could not open file " << fname << endl;

  fclose(fsph);

  // plane data
	/*fname = "points." + to_string(St) + ".dat";
  fid = fopen(fname.c_str(), "w");
  if (fid != NULL)
  {
    for (int i = 0; i < nnodes; i++)
    	fprintf(fid, "%lf\t%lf\t%.8e\n", lambda[i], theta[i], var(i,0));
  }
  else
  	cout << "Could not open file " << fname << endl;

  fclose(fid);*/
}

void TurbModel_EUL::calcReStress(double *struc, double &eps_main, double *prod,
                                 double *rRedi, double *sRediEps,
                                 double (*G)[3], double (*Gph)[3],
                                 double (*Gp1)[3], double dt, char* tIntName)
{
  // time integration
  if      (strcmp(tIntName, "fwEuler") == 0) fwEuler(dt, G);
  else if (strcmp(tIntName, "Heun") == 0)    Heun(dt, G);
  else if (strcmp(tIntName, "CrankN") == 0)  CrankN(dt, G);
  else if (strcmp(tIntName, "RK4") == 0)     RK4(dt, G);
  else    cout << "Not a valid time integration scheme" << endl;

  // structure tensors and eps to main
  for (int i = 0; i < 6; i++)
  {
    struc[i] = rey[i];
    struc[i+6] = dim[i];
  }
  eps_main = eps;

  int index[3][3] = { {0, 3, 4}, {3, 1, 5}, {4, 5, 2} };
  // production
  for (int i = 0; i < 3; i++)
    for (int j = i; j < 3; j++)
      for (int k = 0; k < 3; k++)
      {
      	prod[index[i][j]] = -rey[index[i][k]]*G[j][k] - rey[index[j][k]]*G[i][k];
      }

  // rapid redistribution
  for (int i = 0; i < 3; i++)
    for (int j = i; j < 3; j++)
    {
    	rRedi[index[i][j]] = 0.0;
    	for (int k = 0; k < 3; k++)
    		for (int l = 0; l < 3; l++)
    			rRedi[index[i][j]] += 2.0*G[k][l]*(M[i][l][j][k] + M[j][l][i][k]);
    }

  // slow redistribution + dissipation
  double RD[3][3] = { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0} };
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
      	RD[i][j] += rey[index[i][k]]*dim[index[k][j]];

  for (int i = 0; i < 3; i++)
    for (int j = i; j < 3; j++)
    {
      sRediEps[index[i][j]] = 0.0;
    }

}

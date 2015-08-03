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
  cout << "Simulated TKE: " << tke << endl;
  cout << "Exact TKE    : " << 0.25 << endl;
  cout << "Eps          : " << eps << endl;
  cout << "Tau          : " << tau << endl;
  cout << "Sk/eps       : " << tke/eps << endl;
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
    double small = 1.0e-1;
    double e2 = vecDotVec3d(eipar, eipar);
    if ((e2 < (1.0 - small)) || (e2 > (1.0 + small)))
      cout << "Wrong e2: " << e2 << endl;
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
    //normVec3d(evec); // forced renormalization
    double small = 1.0e-1;
    double e2 = vecDotVec3d(evec, evec);
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
  nshells = 32;
  nmodes = 200;
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

  // Initialize cluster properties
  for (int ishells = 0; ishells < nshells; ishells++)
    for (int imodes = 0; imodes < nmodes; imodes++)
    {
      int icls = ishells*nmodes + imodes;

      double theta;
      double r[3], s[3];

      // wave vector
      double u1 = (double)rand()/RAND_MAX;
      double u2 = (double)rand()/RAND_MAX;
      theta = acos(2.0*u2 - 1.0);
      double phi = 2.0*M_PI*u1;

      e[icls][0] = sin(theta)*cos(phi);
      e[icls][1] = sin(theta)*sin(phi);
      e[icls][2] = cos(theta);

      // cluster tensor
      theta = ( (double)rand()/RAND_MAX )*2.0*M_PI;

      double umag2 = 2.0*nshells*Enukdk[ishells];
      c[icls][0] = 0.5*umag2*( 1.0 - e[icls][0]*e[icls][0] );
      c[icls][1] = 0.5*umag2*( 1.0 - e[icls][1]*e[icls][1] );
      c[icls][2] = 0.5*umag2*( 1.0 - e[icls][2]*e[icls][2] );
      c[icls][3] = 0.5*umag2*( -e[icls][0]*e[icls][1] );
      c[icls][4] = 0.5*umag2*( -e[icls][0]*e[icls][2] );
      c[icls][5] = 0.5*umag2*( -e[icls][1]*e[icls][2] );
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
  cout << "Simulated TKE: " << tke << endl;
  cout << "Exact TKE    : " << 0.25 << endl;
  cout << "Eps          : " << eps << endl;
  cout << "Tau          : " << tau << endl;
  cout << "Sk/eps       : " << tke/eps << endl;
  cout << "--------------------------------" << endl;
}

void TurbModel_CLS::fwEuler(double dt, double (*Gn)[3])
{
  for (int icls = 0; icls < ncls; icls++)
  {
    double *eicls = e[icls];
    double *cicls = c[icls];

    double drift[9];
    driftCoeff(drift, Gn, eicls, cicls);

    for (int i = 0; i < 3; i++) eicls[i] += drift[i]*dt;
    for (int i = 0; i < 6; i++) cicls[i] += drift[i+3]*dt;

    double small = 1.0e-1;
    double e2 = vecDotVec3d(eicls, eicls);
    if ((e2 < (1.0 - small)) || (e2 > (1.0 + small)))
      cout << "Wrong e2: " << e2 << endl;
  }

  double eps_rhs = rhsDissipation(Gn, dt);
  eps += dt*eps_rhs;

  calcTurbStatistics();
  calcTurbTimeScale();
}

void TurbModel_CLS::Heun(double dt, double (*Gn)[3])
{
  for (int icls = 0; icls < ncls; icls++)
  {
    double *eicls = e[icls];
    double *cicls = c[icls];

    double drift_0[9];
    driftCoeff(drift_0, Gn, eicls, cicls);

    double ebar[3],cbar[6];
    for (int i = 0; i < 3; i++) ebar[i] = eicls[i] + drift_0[i]*dt;
    for (int i = 0; i < 6; i++) cbar[i] = cicls[i] + drift_0[i+3]*dt;

    double drift_bar[9];
    driftCoeff(drift_bar, Gn, ebar, cbar);

    for (int i = 0; i < 3; i++) eicls[i] += 0.5*(drift_bar[i] + drift_0[i])*dt;
    for (int i = 0; i < 6; i++) cicls[i] += 0.5*(drift_bar[i+3] + drift_0[i+3])*dt;

    double small = 1.0e-1;
    double e2 = vecDotVec3d(eicls, eicls);
    //if ((e2 < (1.0 - small)) || (e2 > (1.0 + small)))
    //  cout << "Wrong e2: " << e2 << endl;
  }

  double eps_rhs = rhsDissipation(Gn, dt);
  eps += dt*eps_rhs;

  calcTurbStatistics();
  calcTurbTimeScale();
}

void TurbModel_CLS::driftCoeff(double *drift, double (*G)[3], double *eref,
                               double *cref)
{
  for (int i = 0; i < 9; i++)
    drift[i] = 0.0;

  double psi[3][3];
  psi[0][0] = cref[0];   psi[0][1] = cref[3];   psi[0][2] = cref[4];
  psi[1][0] = cref[3];   psi[1][1] = cref[1];   psi[1][2] = cref[5];
  psi[2][0] = cref[4];   psi[2][1] = cref[5];   psi[2][2] = cref[2];

  int iIndex[6] = {0, 1, 2, 0, 0, 1};
  int jIndex[6] = {0, 1, 2, 1, 2, 2};

  double Ge[3];
  double eGe;
  double Gpsi[3][3], eGpsie1[3][3], eGpsie2[3][3];

  // rapid part
  matTransDotVec3d(Ge, eref, G);
  eGe = vecDotMatDotVec3d(eref, eref, G);
  matTimesMat3d(Gpsi, G, psi);
  for (int i = 0; i < 3; i++)
  	for (int j = 0; j < 3; j++)
  	{
  		eGpsie1[i][j] = 0.0;
  		for (int q = 0; q < 3; q++)
  			for (int s = 0; s < 3; s++)
  				eGpsie1[i][j] += eref[q]*G[q][s]*psi[s][i]*eref[j];
  	}

  for (int i = 0; i < 3; i++)
    drift[i] = -Ge[i] + eGe*eref[i];

  for (int n = 0; n < 6; n++)
	{
	  int i = iIndex[n];
	  int j = jIndex[n];
    drift[n+3] = -Gpsi[i][j] - Gpsi[j][i] + 2.0*eGpsie1[i][j] + 2.0*eGpsie1[j][i];
	}

  // expanded gradient
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
  		Gn[i][j] = Cn/tau*rd[i][j];
  		Gv[i][j] = Cv/tau*rd[i][j];
  	}

  matTransDotVec3d(Ge, eref, Gn);
  eGe = vecDotMatDotVec3d(eref, eref, Gn);
  matTimesMat3d(Gpsi, Gv, psi);
  for (int i = 0; i < 3; i++)
  	for (int j = 0; j < 3; j++)
  	{
  		eGpsie1[i][j] = 0.0;
  		eGpsie2[i][j] = 0.0;
  		for (int q = 0; q < 3; q++)
  			for (int s = 0; s < 3; s++)
  			{
  				eGpsie1[i][j] += eref[q]*Gn[q][s]*psi[s][i]*eref[j];
  				eGpsie2[i][j] += eref[q]*Gv[q][s]*psi[s][i]*eref[j];
  			}
  	}

  for (int i = 0; i < 3; i++)
    drift[i] += -Ge[i] + eGe*eref[i];

  for (int n = 0; n < 6; n++)
	{
	  int i = iIndex[n];
	  int j = jIndex[n];
    drift[n+3] += -Gpsi[i][j] - Gpsi[j][i] + eGpsie1[i][j] + eGpsie1[j][i]
                                           + eGpsie2[i][j] + eGpsie2[j][i];
	}

  // slow rotational randomization
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

  double psikk = cref[0] + cref[1] + cref[2];

  drift[3] -= 2.0*C1*cref[0] - C1*psikk + C1*psikk*eref[0]*eref[0];
  drift[4] -= 2.0*C1*cref[1] - C1*psikk + C1*psikk*eref[1]*eref[1];
  drift[5] -= 2.0*C1*cref[2] - C1*psikk + C1*psikk*eref[2]*eref[2];
  drift[6] -= 2.0*C1*cref[3] + C1*psikk*eref[0]*eref[1];
  drift[7] -= 2.0*C1*cref[4] + C1*psikk*eref[0]*eref[2];
  drift[8] -= 2.0*C1*cref[5] + C1*psikk*eref[1]*eref[2];
}

double TurbModel_CLS::rhsDissipation(double(*Gn)[3], double dt)
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

void TurbModel_CLS::calcTurbTimeScale()
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

void TurbModel_CLS:: calcTurbStatistics()
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
  		rey[i] += c[icls][i];

    for (int i = 0; i < 3; i++)
      for(int j = i; j <3; j++)
      {
      	e2[i][j] += e[icls][i]*e[icls][j];
      	for (int k = 0; k< 3; k++)
      		for (int l = k; l < 3; l++)
      		{
      		  M[i][j][k][l] += c[icls][index[i][j]]*e[icls][k]*e[icls][l];
      			e4[i][j][k][l] += e[icls][i]*e[icls][j]*e[icls][k]*e[icls][l];
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
  if      (strcmp(tIntName, "fwEuler") == 0)   fwEuler(dt, Gn);
  else if (strcmp(tIntName, "Heun") == 0)      Heun(dt, Gn);
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

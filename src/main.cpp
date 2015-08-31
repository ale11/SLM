#include "TurbModel_WVM.hpp"

/* Syntax for calling main:
   ./main input_file
   where input_file is the name of the input file.
*/

void WifromWij(double *Wi, double (*Wij)[3]);
void WijfromWi(double *Wi, double (*Wij)[3]);
void evolutionWi(double *Wi, double *Wihat, double St);

int main(int argc, char *argv[])
{

  // Parameters & variables
  int nt;               // number of time steps
  int def;              // type of deformation
  int ierr;             // error flag for reading inputs

  bool seq_axe;                   // used to switch from axe to axc def
  FILE *fid;                      // file ID
  char history[] = "history.dat"; // history file
  char buffer[20];
  char model[10];                 // particle model
  char def_name[10];              // name for the type of deformation
  char tIntName[10];              // time integration name

  double dt;            // increment for the time steps
  double t;             // time
  double ststar = 0.0;  // for sequential simulations
  double S;             // rate-of-strain parameter
  double W;             // rate-of-rotation parameter
  double WoverS;        // W/S
  double Sij[3][3] = {{0.0, 0.0, 0.0},  // rate-of-strain tensor 
                      {0.0, 0.0, 0.0}, 
                      {0.0, 0.0, 0.0}};
  double Wij[3][3] = {{0.0, 0.0, 0.0},  // rate-of-rotation tensor 
                      {0.0, 0.0, 0.0}, 
		      {0.0, 0.0, 0.0}};
  double Wi[3] = {0.0, 0.0, 0.0};    // vorticity vec
  double Wihat[3] = {0.0, 0.0, 0.0}; // vorticity vec used for diagonalization
  double *st;         // normalized time at each time step
  
  double Stau_init;   // initial Sk/eps 
  double struc[12];   // comp and dim tensors
  double tke_ndim;    // turbulent kinetic energy (non-dimensionalized)
  double tke_init;    // initial turbulent kinetic energy
  double eps_ndim;    // dissipation of TKE (non-dimensionalized)
  double eps_init;    // initial dissipation of TKE
  double Sstar;       // Shear parameter: Sq2/eps
  double P_eps;       // Production over dissipation
  double prod[6];     // Production of rij
  double rRedi[6];    // Rapid redistribution of rij
  double sRediEps[6]; // Slow redistribution + dissipation of rij
  
  double (*G)[3][3];  // velocity deformation tensor
  double (*it_Wi)[3]; // vorticity vec at each time step

  TurbModel_SLM *slm;

  enum def_types {axc, axe, ps, rs, seq, thd};

  // Read input file
  fid = fopen(argv[1],"r");
  if (fid == NULL)
  {
    cout << "Could not open file input.dat" << endl;
    return 0;
  }

  ierr = 0;
  ierr += fscanf(fid,"%s = %s", buffer, model);
  ierr += fscanf(fid,"%s = %s", buffer, def_name);
  ierr += fscanf(fid,"%s = %d", buffer, &nt);
  ierr += fscanf(fid,"%s = %lf", buffer, &dt);
  ierr += fscanf(fid,"%s = %lf", buffer, &WoverS);
  ierr += fscanf(fid,"%s = %lf", buffer, &Stau_init);
  ierr += fscanf(fid,"%s = %s", buffer, tIntName);

  if (ierr !=14)
  {
    cout << "Could not read all of the inputs" << endl;
    return 0;
  };

  fclose(fid);

  cout << "--------------------------------" << endl;
  cout << "Deformation   : " << def_name << endl;
  cout << "Time integ.   : " << tIntName << endl;
  cout << "nt            : " << nt << endl;
  cout << "dt            : " << dt << endl;
  cout << "W/S           : " << WoverS << endl;

  if (strcmp(def_name, "axc") == 0)
    def = axc;
  else if (strcmp(def_name, "axe") == 0)
    def = axe;
  else if (strcmp(def_name, "ps") == 0)
    def = ps;
  else if (strcmp(def_name, "rs") == 0)
    def = rs;
  else if (strcmp(def_name, "seq") == 0)
    def = seq;
  else if (strcmp(def_name, "thd") == 0)
    def = thd;
  else
    cout << "Deformation specified not available" << endl;

  switch (def)
  {
  case axc:
    S = 1.0;
    W = 0.0;
    Sij[0][0] = S;
    Sij[1][1] = -0.5*S;
    Sij[2][2] = -0.5*S;
    break;
  case axe:
    S = 1.0;
    W = 0.0;
    Sij[0][0] = -S;
    Sij[1][1] = 0.5*S;
    Sij[2][2] = 0.5*S;
    break;
  case ps:
    S = 1.0;
    W = 0.0;
    Sij[0][0] = S;
    Sij[1][1] = -S;
    break;
  case rs:
    S = 1.0;
    W = WoverS*S;
    Sij[0][1] = 0.5*S;
    Sij[1][0] = 0.5*S;
    Wij[0][1] = 0.5*W;
    Wij[1][0] = -0.5*W;
    break;
  case seq:
    S = -1.0;
    W = WoverS*S;
    Sij[0][0] = S;
    Sij[1][1] = -0.5*S;
    Sij[2][2] = -0.5*S;
    Wij[1][2] = -0.5*W;
    Wij[2][1] = 0.5*W;
    WifromWij(Wi,Wij);
    seq_axe = true;
    break;
  case thd:
    S = 1.0;
    W = WoverS*S;
    Sij[0][1] = 0.5*S;
    Sij[1][0] = 0.5*S;
    Sij[1][2] = 0.5*S;
    Sij[2][1] = 0.5*S;
    Wij[0][1] = 0.5*W;
    Wij[1][0] = -0.5*W;
    Wij[1][2] = -0.125*W;
    Wij[2][1] = 0.125*W;
    WifromWij(Wi,Wij);
    Wihat[0] = -0.5*Wi[0] + 0.0*Wi[1] + 0.5*Wi[2];
    Wihat[1] = 0.25*Wi[0] - 0.25*sqrt(2.0)*Wi[1] + 0.25*Wi[2];
    Wihat[2] = 0.25*Wi[0] + 0.25*sqrt(2.0)*Wi[1] + 0.25*Wi[2];
    break;
  default:
    cout << "Not a valid deformation" << endl;
    break;
  }

  // Data allocation
  st    = new double [nt+1];
  G     = new double [2*nt+1][3][3]; 
  it_Wi = new double [nt+1][3];
  if (strcmp(model, "iprm") == 0)      slm = new TurbModel_IPRM;
  else if (strcmp(model, "lang") == 0) slm = new TurbModel_LANG;
  else if (strcmp(model, "cls") == 0)  slm = new TurbModel_CLS;
  else cout << "Model " << model << " not available." << endl;

  // Initial condition
  t = 0.0;
  st[0] = fabs(S)*t;
  for (int i = 0; i < 3; i++)
  {
    it_Wi[0][i] = Wi[i];
    for (int j = 0; j < 3; j++)
      G[0][i][j] = Sij[i][j] + Wij[i][j];
  }

  slm->initialHookScalarRansTurbModel(Stau_init,struc,eps_init, prod, rRedi,
                                      G[0]);
  tke_init = 0.5*(struc[0] + struc[1] + struc[2]);
  double q2 = 2.0*tke_init;
  Sstar = S*q2/eps_init;
  P_eps = -S*struc[3]/eps_init;

  fid = fopen(history, "w");
  if (fid == NULL)
  {
    cout << "Could not open file " << history << endl;
    return 0;
  }
  fprintf(fid,"Variables = \"st\" \"e_st\" \"rey0_2k\" \"rey1_2k\" ");
  fprintf(fid,"\"rey2_2k\" \"rey3_2k\" \"rey4_2k\" \"rey5_2k\" ");
  fprintf(fid,"\"dim0_2k\" \"dim1_2k\" \"dim2_2k\" \"dim3_2k\" \"dim4_2k\" ");
  fprintf(fid,"\"dim5_2k\" \"tke\" \"eps\" \"S*\" \"P/eps\" ");
  fprintf(fid,"\"prod_0\" \"prod_1\" \"prod_2\" \"prod_3\" \"prod_4\" ");
  fprintf(fid,"\"prod_5\" \"rredi_0\" \"rredi_1\" \"rredi_2\" \"rredi_3\" ");
  fprintf(fid,"\"rredi_4\" \"rredi_5\" \n");
  fprintf(fid,"Zone T = \"%s\" F = point\n", model);
  fprintf(fid,"%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t",
	  st[0], exp(st[0]), struc[0]/q2, struc[1]/q2, struc[2]/q2, 
          struc[3]/q2, struc[4]/q2, struc[5]/q2);
  fprintf(fid,"%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t",
          struc[6]/q2, struc[7]/q2, struc[8]/q2, struc[9]/q2, struc[10]/q2, 
          struc[11]/q2, 1.0, 1.0);
  fprintf(fid,"%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t",
	  Sstar, P_eps, prod[0]/(S*q2), prod[1]/(S*q2), prod[2]/(S*q2),
	  prod[3]/(S*q2), prod[4]/(S*q2), prod[5]/(S*q2));
  fprintf(fid, "%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n",
	  rRedi[0]/(S*q2), rRedi[1]/(S*q2), rRedi[2]/(S*q2),
	  rRedi[3]/(S*q2), rRedi[4]/(S*q2), rRedi[5]/(S*q2));

  //slm->writeData(st[0]);
  // -----------------------------------------------------------------------------------------
  // Loop through time
  // -----------------------------------------------------------------------------------------
  for (int n = 0; n < nt; n++) // 1 is replacing nt
  {
    // compute time
    t += dt;
    st[n+1] = fabs(S)*t;

    // compute deformations
    switch (def)
    {
    case seq:
    	double W_temp;
      if ( (exp(st[n+1]) > 2.72) && (seq_axe) )
      {
      	if ( exp(st[n] + fabs(S)*0.5*dt) > 2.72 )
      	{
      		ststar = S*(t - 0.5*dt);
      		S = 1.0;
      		Sij[0][0] = S;
      		Sij[1][1] = -0.5*S;
      		Sij[2][2] = -0.5*S;

      		W_temp = W*exp(S*(t - 0.5*dt) + 2.0*ststar);
      		Wij[1][2] = -0.5*W_temp;
      		Wij[2][1] = 0.5*W_temp;
      		matPlusMat3d(G[2*n+1], Sij, Wij);

      		W_temp = W*exp(S*t + 2.0*ststar);
      		Wij[1][2] = -0.5*W_temp;
      		Wij[2][1] = 0.5*W_temp;
      		WifromWij(Wi, Wij);
      		matPlusMat3d(G[2*n+2], Sij, Wij);
      	}
      	else
      	{
      		ststar = S*t;
      		S = 1.0;
          Sij[0][0] = S;
          Sij[1][1] = -0.5*S;
          Sij[2][2] = -0.5*S;

          W_temp = W*exp(-S*(t - 0.5*dt));
          Wij[1][2] = -0.5*W_temp;
          Wij[2][1] = 0.5*W_temp;
          matPlusMat3d(G[2*n+1], Sij, Wij);

          W_temp = W*exp(S*t + 2.0*ststar);
          Wij[1][2] = -0.5*W_temp;
          Wij[2][1] = 0.5*W_temp;
          WifromWij(Wi, Wij);
          matPlusMat3d(G[2*n+2], Sij, Wij);
      	}
      	seq_axe = false;
      }
      else 
      {
      	W_temp = W*exp(S*(t - 0.5*dt) + 2.0*ststar);
      	Wij[1][2] = -0.5*W_temp;
      	Wij[2][1] = 0.5*W_temp;
      	matPlusMat3d(G[2*n+1], Sij, Wij);

        W_temp = W*exp(S*t + 2.0*ststar);
        Wij[1][2] = -0.5*W_temp;
        Wij[2][1] = 0.5*W_temp;
        WifromWij(Wi, Wij);
        matPlusMat3d(G[2*n+2], Sij, Wij);
      }
      break;
    case thd:
      evolutionWi(Wi, Wihat, S*(t - 0.5*dt));
      WijfromWi(Wi, Wij);
      matPlusMat3d(G[2*n+1], Sij, Wij);

      evolutionWi(Wi, Wihat, S*t);
      WijfromWi(Wi, Wij);
      matPlusMat3d(G[2*n+2], Sij, Wij);
      break;
    case axc:
    case axe:
    case ps:
    case rs:
      for (int i = 0; i < 3; i++)
	for (int j = 0; j < 3; j++)
	{
	  G[2*n+1][i][j] = Sij[i][j] + Wij[i][j];
	  G[2*n+2][i][j] = Sij[i][j] + Wij[i][j];
	}
      break;
    }
    
    // compute Reynolds stresses and TKE
    slm->calcReStress(struc, eps_ndim, prod, rRedi, sRediEps, G[2*n], G[2*n+1], 
                      G[2*n+2], dt, tIntName);
    tke_ndim = 0.5*(struc[0] + struc[1] + struc[2]);
    double q2 = 2.0*tke_ndim;
    Sstar = S*q2/eps_ndim;
    P_eps = -S*struc[3]/eps_ndim;
    tke_ndim /= tke_init;
    eps_ndim /= eps_init;

    // write data
    fprintf(fid,"%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t",
	    st[n+1], exp(st[n+1]), struc[0]/q2, struc[1]/q2, struc[2]/q2,
	    struc[3]/q2, struc[4]/q2, struc[5]/q2);
    fprintf(fid,"%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t",
	    struc[6]/q2, struc[7]/q2, struc[8]/q2, struc[9]/q2, struc[10]/q2,
	    struc[11]/q2, tke_ndim, eps_ndim);
    fprintf(fid,"%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t", 
            Sstar, P_eps, prod[0]/(S*q2), prod[1]/(S*q2), prod[2]/(S*q2),
            prod[3]/(S*q2), prod[4]/(S*q2), prod[5]/(S*q2));
    fprintf(fid, "%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n",
	    rRedi[0]/(S*q2), rRedi[1]/(S*q2), rRedi[2]/(S*q2),
            rRedi[3]/(S*q2), rRedi[4]/(S*q2), rRedi[5]/(S*q2));

    // Output iteration number
    cout << "iter: " << setw(6) << st[n+1] << setw(12) << tke_ndim << endl;

    // Write data
    //if ( (n+1) % 200 == 0) slm->writeData(st[n+1]);
  }

  fclose(fid);

  // Cleanup
  delete [] st;
  delete [] G;
  delete [] it_Wi;
  delete slm;

  return (0);
}

void WifromWij(double *Wi, double (*Wij)[3])
{
  Wi[0] = Wij[2][1] - Wij[1][2];
  Wi[1] = Wij[0][2] - Wij[2][0];
  Wi[2] = Wij[1][0] - Wij[0][1];
}

void WijfromWi(double *Wi, double (*Wij)[3])
{
  Wij[0][1] = -0.5*Wi[2];
  Wij[0][2] = 0.5*Wi[1];
  Wij[1][2] = -0.5*Wi[0];

  Wij[1][0] = -Wij[0][1];
  Wij[2][0] = -Wij[0][2];
  Wij[2][1] = -Wij[1][2];
}

void evolutionWi(double *Wi, double * const Wihat, double St)
{
  double Wihat_t[3];

  Wihat_t[0] = Wihat[0];
  Wihat_t[1] = Wihat[1]*exp(-1.0/sqrt(2.0)*St);
  Wihat_t[2] = Wihat[2]*exp(1.0/sqrt(2.0)*St);

  Wi[0] = -Wihat_t[0] + Wihat_t[1] + Wihat_t[2];
  Wi[1] = -sqrt(2.0)*Wihat_t[1] + sqrt(2.0)*Wihat_t[2];
  Wi[2] = Wihat_t[0] + Wihat_t[1] + Wihat_t[2];
}

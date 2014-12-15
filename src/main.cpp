#include "TurbModel_WVM.hpp"

/* Syntax for calling main:
   ./main model def ratio
   where model: particle model to be used
	 def  : "axc", "axe", etc.
         ratio: is the ratio of rotation to strain, which needs to be specified
                for some deformations only.
*/

void WifromWij(double *Wi, double (*Wij)[3]);
void WijfromWi(double *Wi, double (*Wij)[3]);
void evolutionWi(double *Wi, double *Wihat, double St);

int main(int argc, char *argv[])
{

  // Parameters & variables
  int nt;               // number of time steps
  int def;              // type of deformation

  bool seq_axe;                   // used to switch from axe to axc def
  FILE *fid;                      // file ID
  char history[] = "history.dat"; // history file

  double dt;            // increment for the time steps
  double t;             // time
  double ststar = 0.0;  // for sequential simulations
  double S;             // rate-of-strain parameter
  double W;             // rate-of-rotation parameter
  double Sij[3][3] = {{0.0, 0.0, 0.0},  // rate-of-strain tensor 
                      {0.0, 0.0, 0.0}, 
                      {0.0, 0.0, 0.0}};
  double Wij[3][3] = {{0.0, 0.0, 0.0},  // rate-of-rotation tensor 
                      {0.0, 0.0, 0.0}, 
		      {0.0, 0.0, 0.0}};
  double Wi[3] = {0.0, 0.0, 0.0};    // vorticity vec
  double Wihat[3] = {0.0, 0.0, 0.0}; // vorticity vec used for diagonalization
  double *st;         // normalized time at each time step
  double struc[12];   // comp and dim tensors
  double tke_ndim;    // turbulent kinetic energy (non-dimensionalized)
  double tke_init;    // initial turbulent kinetic energy
  double eps_ndim;    // dissipation of TKE (non-dimensionalized)
  double eps_init;    // initial dissipation of TKE
  double (*G)[3][3];  // velocity deformation tensor
  double (*it_Wi)[3]; // vorticity vec at each time step

  TurbModel_WVM *wvm;

  enum def_types {axc, axe, ps, rs, seq, thd};

  if (argc < 3)
  {
    cout << "Not enough input arguments." << endl;
    return 0;
  }

  if (strcmp(argv[2], "axc") == 0)
    def = axc;
  else if (strcmp(argv[2], "axe") == 0)
    def = axe;
  else if (strcmp(argv[2], "ps") == 0)
    def = ps;
  else if (strcmp(argv[2], "rs") == 0)
    def = rs;
  else if (strcmp(argv[2], "seq") == 0)
    def = seq;
  else if (strcmp(argv[2], "thd") == 0)
    def = thd;
  else
    cout << "Deformation specified not available" << endl;

  switch (def)
  {
  case axc:
    if (argc != 3) cout << "Wrong number of inputs specified." << endl;
    cout << "Deformation: axc" << endl;
    S = 1.0;
    W = 0.0;
    Sij[0][0] = S;
    Sij[1][1] = -0.5*S;
    Sij[2][2] = -0.5*S;
    nt = 100;
    dt = 0.025;
    break;
  case axe:
    if (argc != 3) cout << "Wrong number of inputs specified." << endl;
    cout << "Deformation: axe" << endl;
    S = 1.0;
    W = 0.0;
    Sij[0][0] = -S;
    Sij[1][1] = 0.5*S;
    Sij[2][2] = 0.5*S;
    nt = 100;
    dt = 0.025;
    break;
  case ps:
    if (argc != 3) cout << "Wrong number of inputs specified." << endl;
    cout << "Deformation: ps" << endl;
    S = 1.0;
    W = 0.0;
    Sij[0][0] = S;
    Sij[1][1] = -S;
    nt = 500;
    dt = 0.005;
    break;
  case rs:
    if (argc != 4) cout << "Wrong number of inputs specified." << endl;
    cout << "Deformation: rs " << argv[3] << endl;
    S = 1.0;
    W = atof(argv[3])*S;
    Sij[0][1] = 0.5*S;
    Sij[1][0] = 0.5*S;
    Wij[0][1] = 0.5*W;
    Wij[1][0] = -0.5*W;
    nt = 1000;
    if (atof(argv[3]) >= 4.0) dt = 0.025;
    else dt = 0.005; 
    break;
  case seq:
    if (argc != 3) cout << "Wrong number of inputs specified." << endl;
    cout << "Deformation: seq " << endl;
    S = -1.0;
    W = 1.5*S;
    Sij[0][0] = S;
    Sij[1][1] = -0.5*S;
    Sij[2][2] = -0.5*S;
    Wij[1][2] = -0.5*W;
    Wij[2][1] = 0.5*W;
    WifromWij(Wi,Wij);
    seq_axe = true;
    nt = 100;
    dt = 0.05;
    break;
  case thd:
    if (argc != 4) cout << "Wrong number of inputs specified." << endl;
    cout << "Deformation: thd "<< argv[3] << endl;
    S = 1.0;
    W = atof(argv[3])*S;
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
    nt = 200;
    dt = 0.025;
    break;
  default:
    cout << "Not a valid deformation" << endl;
    break;
  }

  // Data allocation
  st    = new double [nt+1];
  G     = new double [2*nt+1][3][3]; 
  it_Wi = new double [nt+1][3];    
  wvm   = new TurbModel_WVM(argv[1]);

  // Initial condition
  t = 0.0;
  st[0] = fabs(S)*t;
  for (int i = 0; i < 3; i++)
  {
    it_Wi[0][i] = Wi[i];
    for (int j = 0; j < 3; j++)
      G[0][i][j] = Sij[i][j] + Wij[i][j];
  }

  wvm->initialHookScalarRansTurbModel(struc,eps_init);
  tke_init = 0.5*(struc[0] + struc[1] + struc[2]);
  double q2 = 2.0*tke_init;

  fid = fopen(history, "w");
  if (fid == NULL)
  {
    cout << "Could not open file " << history << endl;
    return 0;
  }
  fprintf(fid,"Variables = \"st\" \"e_st\" \"rey0_2k\" \"rey1_2k\" ");
  fprintf(fid,"\"rey2_2k\" \"rey3_2k\" \"rey4_2k\" \"rey5_2k\" ");
  fprintf(fid,"\"dim0_2k\" \"dim1_2k\" \"dim2_2k\" \"dim3_2k\" \"dim4_2k\" ");
  fprintf(fid,"\"dim5_2k\" \"tke\" \"eps\" \n");
  fprintf(fid,"Zone T = \"%s\" F = point\n",argv[1]);
  fprintf(fid,"%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t",
	  st[0], exp(st[0]), struc[0]/q2, struc[1]/q2, struc[2]/q2, 
          struc[3]/q2, struc[4]/q2, struc[5]/q2);
  fprintf(fid,"%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n",
          struc[6]/q2, struc[7]/q2, struc[8]/q2, struc[9]/q2, struc[10]/q2, 
          struc[11]/q2, 1.0, 1.0);

  // Loop through time
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

    // update dissipation
    eps_ndim = wvm->updateDissipation(G[2*n], dt);
    eps_ndim /= eps_init;
    
    // compute Reynolds stresses and TKE
    wvm->calcReStress(struc, G[2*n], G[2*n+1], G[2*n+2], dt);
    tke_ndim = 0.5*(struc[0] + struc[1] + struc[2]);
    double q2 = 2.0*tke_ndim;
    tke_ndim /= tke_init;

    // write data
    fprintf(fid,"%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t",
	    st[n+1], exp(st[n+1]), struc[0]/q2, struc[1]/q2, struc[2]/q2,
	    struc[3]/q2, struc[4]/q2, struc[5]/q2);
    fprintf(fid,"%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n",
	    struc[6]/q2, struc[7]/q2, struc[8]/q2, struc[9]/q2, struc[10]/q2,
	    struc[11]/q2, tke_ndim, eps_ndim);

    // Output iteration number
    cout << "iter: " << st[n+1] << endl;
  }

  fclose(fid);

  // Cleanup
  delete [] st;
  delete [] G;
  delete [] it_Wi;
  delete wvm;

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

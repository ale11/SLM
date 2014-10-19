#include <iostream>
#include <cstdlib>
#include <math.h>

using namespace std;

double normVec3d(const double *v, double *vNorm);
double vecDotVec3d(const double *v1, const double *v2);
void vecCrossVec3d(const double *va, const double *vb, double *vc);

int main(int argc, char*argv[])
{
  int nshells;           // number of shells
  int nmodes;            // number of modes per shell
  int nt;                // number of time steps

  double *st;            // S times timescale
  double *e_st;          // exp(st)
  double (*it_e)[3];     // container for the unit wave vector
  double (*it_a)[3];     // container for the eddy-axis vector
  double (*it_u)[3];     // container for the velocity vector
  double (*it_Wi)[3];    // container for the mean-vorticity vector
  double (*rij_2k)[6];   // reynolds stresses divided by 2k
  double *tke;           // turbulent kinetic energy
  double (*vxa)[6];      // < V^2 * xi_star * a_i a_i > / < V^2 >
  double *jet;           // < V^2 * xi_star > / < V^2 >
  double *hel;           // < V^2 * gamma_star > / < V^2 >
  double (*aij)[6];      // < V^2 *a_i a_j > / <V^2>
  double (*sbm_rey)[6];  // SBM predicted reynolds stresses, divided by 2k
  double normWi;         // norm of the mean-vorticity vector
  double unitWi[3];      // unitary mean-vorticity vector
  double unitUi[3];      // unitary velocity vector
  double (*debug)[3];    // debugging holder

  char outputfile[] = "outputfile.dat";
  FILE *fid;

  fid = fopen(argv[2], "r");
  if (fid == NULL)
  {
    cout << "could not open file " << argv[2] << endl;
    return 0;
  }
  
  fscanf(fid,"nt=%d\tnshells=%d\tnmodes=%d", &nt, &nshells, &nmodes);
  cout << "Opened file: " << argv[2] << endl;  
  cout << "nshells: " << nshells << " nmodes: " << nmodes << endl;  
  // allocate memory for data
  st      = new double[nt];
  e_st    = new double[nt];
  it_e    = new double[nt*nshells*nmodes][3];
  it_a    = new double[nt*nshells*nmodes][3];
  it_u    = new double[nt*nshells*nmodes][3];
  it_Wi   = new double[nt][3];
  rij_2k  = new double[nt][6];
  tke     = new double[nt];
  vxa     = new double[nt][6];
  jet     = new double[nt];
  hel     = new double[nt];
  aij     = new double[nt][6];
  sbm_rey = new double[nt][6];
  debug   = new double[nt][3];

  // initialize some data
  for (int it = 0; it < nt; it++)
  {
    vxa[it][0] = 0.0;
    vxa[it][1] = 0.0;
    vxa[it][2] = 0.0;
    vxa[it][3] = 0.0;
    vxa[it][4] = 0.0;
    vxa[it][5] = 0.0;

    jet[it] = 0.0;
    hel[it] = 0.0;

    aij[it][0] = 0.0;
    aij[it][1] = 0.0;
    aij[it][2] = 0.0;
    aij[it][3] = 0.0;
    aij[it][4] = 0.0;
    aij[it][5] = 0.0;

    debug[it][0] = 0.0;
    debug[it][1] = 0.0;
    debug[it][2] = 0.0;
  }

  // read data
  for (int it = 0; it < nt; it++)
  {
    fscanf(fid,"%lf",&st[it]);

    for (int ishells = 0; ishells < nshells; ishells++)
      for (int imodes = 0; imodes < nmodes; imodes++)
      {
	int index = it*nshells*nmodes + ishells*nmodes + imodes;
	
	fscanf(fid,"%lf",&it_e[index][0]);
	fscanf(fid,"%lf",&it_e[index][1]);
	fscanf(fid,"%lf",&it_e[index][2]);

	fscanf(fid,"%lf",&it_a[index][0]);
	fscanf(fid,"%lf",&it_a[index][1]);
	fscanf(fid,"%lf",&it_a[index][2]);

	fscanf(fid,"%lf",&it_u[index][0]);
	fscanf(fid,"%lf",&it_u[index][1]);
	fscanf(fid,"%lf",&it_u[index][2]);
      }

    fscanf(fid,"%lf",&it_Wi[it][0]);
    fscanf(fid,"%lf",&it_Wi[it][1]);
    fscanf(fid,"%lf",&it_Wi[it][2]);

    fscanf(fid,"%lf",&rij_2k[it][0]);
    fscanf(fid,"%lf",&rij_2k[it][1]);
    fscanf(fid,"%lf",&rij_2k[it][2]);
    fscanf(fid,"%lf",&rij_2k[it][3]);
    fscanf(fid,"%lf",&rij_2k[it][4]);
    fscanf(fid,"%lf",&rij_2k[it][5]);

    fscanf(fid,"%lf",&tke[it]);  
  } 

  fclose(fid);

  // iterate through time
  for (int it = 0; it < nt; it++)
  {
    e_st[it] = exp(st[it]);

    double unitWi[3];
    double magWi = normVec3d(it_Wi[it], unitWi);

    //debug[it][0] = unitWi[0];
    //debug[it][1] = unitWi[1];
    //debug[it][2] = unitWi[2];
    double counter = 0.0;

    for (int ishells = 0; ishells < nshells; ishells++)
      for (int imodes = 0; imodes < nmodes; imodes++)
      {
	int index = it*nshells*nmodes + ishells*nmodes + imodes;

	double unitUi[3];
	double magUi = normVec3d(it_u[index], unitUi);
	double magUi2 = magUi*magUi;

	double *a = it_a[index];
        double eCrossa[3];
	vecCrossVec3d(it_e[index], a, eCrossa);

	double u1 = vecDotVec3d(unitUi, eCrossa);
	double u2 = vecDotVec3d(unitUi, a);

	double xi = u1*u1;
	
	double unitWini = vecDotVec3d(unitWi, it_e[index]);
	double ga = u1*u2/unitWini;

	vxa[it][0] += magUi2*xi*a[0]*a[0];
	vxa[it][1] += magUi2*xi*a[1]*a[1];
	vxa[it][2] += magUi2*xi*a[2]*a[2];
	vxa[it][3] += magUi2*xi*a[0]*a[1];
	vxa[it][4] += magUi2*xi*a[0]*a[2];
	vxa[it][5] += magUi2*xi*a[1]*a[2];

	jet[it] += magUi2*xi;
	if (unitWini > 1.0e-2)
	  hel[it] += magUi2*ga;
	else
	  counter += 1.0;
	

	aij[it][0] += magUi2*a[0]*a[0];
	aij[it][1] += magUi2*a[1]*a[1];
	aij[it][2] += magUi2*a[2]*a[2];
	aij[it][3] += magUi2*a[0]*a[1];
	aij[it][4] += magUi2*a[0]*a[2];
	aij[it][5] += magUi2*a[1]*a[2];
      }

    debug[it][0] = counter;

    double twiceTke = 2.0*tke[it];
    vxa[it][0] /= (nshells*nmodes*twiceTke);
    vxa[it][1] /= (nshells*nmodes*twiceTke);
    vxa[it][2] /= (nshells*nmodes*twiceTke);
    vxa[it][3] /= (nshells*nmodes*twiceTke);
    vxa[it][4] /= (nshells*nmodes*twiceTke);
    vxa[it][5] /= (nshells*nmodes*twiceTke);

    jet[it] /= (nshells*nmodes*twiceTke);
    hel[it] /= (nshells*nmodes*twiceTke);

    aij[it][0] /= (nshells*nmodes*twiceTke);
    aij[it][1] /= (nshells*nmodes*twiceTke);
    aij[it][2] /= (nshells*nmodes*twiceTke);
    aij[it][3] /= (nshells*nmodes*twiceTke);
    aij[it][4] /= (nshells*nmodes*twiceTke);
    aij[it][5] /= (nshells*nmodes*twiceTke);

    sbm_rey[it][0] = 0.5*jet[it]*(1.0 - aij[it][0]) 
                   + (1.0 - jet[it])*aij[it][0]
                   + 0.5*hel[it]*(unitWi[1]*aij[it][4] - unitWi[2]*aij[it][3]
				+ unitWi[1]*aij[it][4] - unitWi[2]*aij[it][3]);

    sbm_rey[it][1] = 0.5*jet[it]*(1.0 - aij[it][1])
                   + (1.0 - jet[it])*aij[it][1]
                   + 0.5*hel[it]*(unitWi[2]*aij[it][3] - unitWi[0]*aij[it][5]
		                + unitWi[2]*aij[it][3] - unitWi[0]*aij[it][5]);

    sbm_rey[it][2] = 0.5*jet[it]*(1.0 - aij[it][2])
                   + (1.0 - jet[it])*aij[it][2]
                   + 0.5*hel[it]*(unitWi[0]*aij[it][5] - unitWi[1]*aij[it][4]
		                + unitWi[0]*aij[it][5] - unitWi[1]*aij[it][4]);

    sbm_rey[it][3] = -0.5*jet[it]*aij[it][3]
                   + (1.0 - jet[it])*aij[it][3]
                   + 0.5*hel[it]*(unitWi[1]*aij[it][5] - unitWi[2]*aij[it][1]
		                + unitWi[2]*aij[it][0] - unitWi[0]*aij[it][4]);

    sbm_rey[it][4] = -0.5*jet[it]*aij[it][4]
                   + (1.0 - jet[it])*aij[it][4]
                   + 0.5*hel[it]*(unitWi[1]*aij[it][2] - unitWi[2]*aij[it][5]
		                + unitWi[0]*aij[it][3] - unitWi[1]*aij[it][0]);

    sbm_rey[it][5] = -0.5*jet[it]*aij[it][5]
                   + (1.0 - jet[it])*aij[it][5]
                   + 0.5*hel[it]*(unitWi[2]*aij[it][4] - unitWi[0]*aij[it][2]
		                + unitWi[0]*aij[it][1] - unitWi[1]*aij[it][3]);

    cout << "st: " << st[it] << endl;
  }

  // write data
  fid = fopen(outputfile,"w");

  if (fid == NULL)
  { 
    cout << "could not open file " << outputfile << endl;
    return 0;
  }
  fprintf(fid,"Variables = \"st\" \"e_st\" \"r11_2k\" \"r22_2k\" \"r33_2k\" ");
  fprintf(fid,"\"r12_2k\" \"r13_2k\" \"r23_2k\" \"sbm_r11\" \"sbm_r22\" ");
  fprintf(fid," \"sbm_r33\" \"sbm_r12\" \"sbm_r13\" \"sbm_r23\" ");
  fprintf(fid," \"v2_xi_aa11_2k\" \"v2_xi_aa22_2k\" \"v2_xi_aa33_2k\" ");
  fprintf(fid," \"v2_xi_aa12_2k\" \"v2_xi_aa13_2k\" \"v2_xi_aa23_2k\" ");
  fprintf(fid," \"jet\" \"hel\" \"a11\" \"a22\" \"a33\" \"a12\" \"a13 \" ");
  fprintf(fid," \"a23\" \"debug1\" \"debug2\" \"debug3\"\n");

  fprintf(fid,"Zone T = \"%s\" F = point\n",argv[1]);
  for (int it = 0; it < nt; it++)
  {
    fprintf(fid,"%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t",
	    st[it], e_st[it], rij_2k[it][0], rij_2k[it][1], rij_2k[it][2],
            rij_2k[it][3], rij_2k[it][4], rij_2k[it][5], sbm_rey[it][0],
            sbm_rey[it][1]);
    fprintf(fid,"%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t",
	    sbm_rey[it][2], sbm_rey[it][3], sbm_rey[it][4], sbm_rey[it][5],
            vxa[it][0], vxa[it][1], vxa[it][2], vxa[it][3], vxa[it][4],
            vxa[it][5]);
    fprintf(fid,"%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t",
	    jet[it], hel[it], aij[it][0], aij[it][1], aij[it][2], aij[it][3],
            aij[it][4], aij[it][5], debug[it][0], debug[it][1]);
    fprintf(fid,"%.8e\n",
	    debug[it][2]);
  } 
  fclose(fid);

  // clean up data
  delete [] st;
  delete [] e_st;
  delete [] it_e;
  delete [] it_a;
  delete [] it_u;
  delete [] it_Wi;
  delete [] rij_2k;
  delete [] tke;
  delete [] vxa;
  delete [] jet;
  delete [] hel;
  delete [] aij;
  delete [] sbm_rey;
  delete [] debug;

  return 0;
}

double normVec3d(const double *v, double *vNorm)
{
  double mag = vecDotVec3d(v,v);
  mag = sqrt(mag);

  if ( mag > 1.0e-12)
  {
    vNorm[0] = v[0]/mag;
    vNorm[1] = v[1]/mag;
    vNorm[2] = v[2]/mag;
  }
  else
  {
    vNorm[0] = 0.0;
    vNorm[1] = 0.0;
    vNorm[2] = 0.0;
  }
  return mag;
}

double vecDotVec3d(const double *va, const double *vb)
{
  double dotproduct = va[0]*vb[0] + va[1]*vb[1] + va[2]*vb[2];
  return dotproduct;
}

void vecCrossVec3d(const double *va, const double *vb, double *vc)
{
  vc[0] = va[1]*vb[2] - va[2]*vb[1];
  vc[1] = va[2]*vb[0] - va[0]*vb[2];
  vc[2] = va[0]*vb[1] - va[1]*vb[0];
}

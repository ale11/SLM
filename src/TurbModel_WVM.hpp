#ifndef TURBMODEL_WVM_H
#define TURBMODEL_WVM_H

#include "UgpWithCvCompFlow.hpp"
#include <random>

class TurbModel_SLM: public UgpWithCvCompFlow
{ 
// Constructors & Destructor
public:
  TurbModel_SLM();
  virtual ~TurbModel_SLM();

  // Functions
public:
  virtual void initialHookScalarRansTurbModel(double Stau, double *struc, double &eps_init, 
                                              double *prod, double *rapRedi, double (*Gn)[3]);
  
  virtual void calcReStress(double *struc, double &eps_main, double *prod, double *rRedi, 
                            double *sRediEps, double (*Gn)[3], double (*Gnph)[3], 
                            double (*Gnp1)[3], double dt, char* tIntName);

  virtual void writeData(double St);
};

class TurbModel_PAR: public TurbModel_SLM
{
  // Constructors & Destructor
public:
  TurbModel_PAR();
  virtual ~TurbModel_PAR();

  // Data memebers
public:
  int nshells; // number of shells
  int nmodes;  // number of modes per shell
  int npar;    // number of particles

  double (*e)[3];
  double (*a)[3];
  double (*u)[3];

  double eps, tau;
  double M[3][3][3][3], e2[3][3], e4[3][3][3][3];

  default_random_engine generator;
  normal_distribution<double> distribution;

  // Memeber functions
public:
  void initialHookScalarRansTurbModel(double Stau, double *struc, 
                                      double &eps_init, double *prod,
                                      double *rapRedi, double (*Gn)[3]);

  void calcReStress(double *struc, double &eps_main, double *prod,
                    double *rRedi, double *sRediEps, double (*Gn)[3],
                    double (*Gnph)[3], double (*Gnp1)[3], double dt,
                    char* tIntName);

  void fwEuler(double dt, double (*Gn)[3]);

  void Heun(double dt, double (*Gn)[3]);

  void semiHeun(double dt, double (*Gn)[3]);

  void calcTurbStatistics();

  void writeData(double St);

  virtual void driftCoeff(double *drift, double (*G)[3], double *eref,
                          double *uref);

  virtual void diffCoeff(double (*diff)[6], double *eref, double *uref);

  virtual double rhsDissipation(double (*Gn)[3], double dt);

  virtual void calcTurbTimeScale();

};

class TurbModel_IPRM: public TurbModel_PAR
{
  // Constructors & Destructor
public:
  TurbModel_IPRM();
  ~TurbModel_IPRM() {};

  // Data members
public:
  double Cn, Cv, Ceps1, Ceps2, Ceps3;
  
  // Memeber functions
public:
  void driftCoeff(double *drift, double (*G)[3], double *eref, double *uref);
  void diffCoeff(double (*diff)[6], double *eref, double *uref);
  double rhsDissipation(double (*Gn)[3], double dt);
  void calcTurbTimeScale();
};

class TurbModel_LANG: public TurbModel_PAR
{
  // Constructors & Destructor
public:
  TurbModel_LANG();
  ~TurbModel_LANG() {};

  // Data members
public:
  double ae, au, gamma, Ceps1, Ceps2, Ceps3;

  // Memeber functions
public:
  void driftCoeff(double *drift, double (*G)[3], double *eref, double *uref);
  void diffCoeff(double (*diff)[6], double *eref, double *uref);
  double rhsDissipation(double (*Gn)[3], double dt);
};

class TurbModel_CLS: public TurbModel_SLM
{
  // Constructors & Destructor
public:
  TurbModel_CLS();
  virtual ~TurbModel_CLS();

  // Data memebers
public:
  int nshells; // number of shells
  int nmodes;  // number of modes per shell
  int ncls;    // number of clusters

  double (*e)[3];
  double (*c)[6];
  mat var;

  double Gn[3][3], Gv[3][3], C1[3][3];
  double eps;
  double M[3][3][3][3], e2[3][3], e4[3][3][3][3];

  double Cn, Cv, Ceps1, Ceps2, Ceps3;

  // Memeber functions
public:
  void initialHookScalarRansTurbModel(double Stau, double *struc, double &eps_init,
  		                                double *prod, double *rapRedi, double (*Gn)[3]);

  void calcReStress(double *struc, double &eps_main, double *prod, double *rRedi,
  		              double *sRediEps, double (*Gn)[3], double (*Gnph)[3], double (*Gnp1)[3],
  		              double dt, char* tIntName);

  void fwEuler(double dt, double (*Gn)[3]);

  void Heun(double dt, double (*Gn)[3]);

  void CrankN(double dt, double (*Gn)[3]);

  virtual void calcRhs(mat &rhs, mat &var, double dt);

  virtual void calcRhsJacob(double (*eDrift)[3], double (*cDrift)[6], double *eref,
  		                  double *cref, double (*G)[3], double dt);

  virtual void correction();

  virtual void inputs(double (*G)[3]);

  virtual void outputs();

  virtual double calcRhsEps(double (*Gn)[3], double dt);

};

class TurbModel_EUL: public TurbModel_SLM
{
  // Constructors & Destructor
public:
  TurbModel_EUL();
  virtual ~TurbModel_EUL();

  // Data memebers
public:
  int nnodes;         // number of nodes for RBFs
  int nelems;         // number of triangulated elements
  int (*elem)[3];     // contains connectivity for triangulated elements

  double epsilon;
  double (*grid)[3];  // contains grid points for RBFs
  double *lambda;
  double *theta;

  mat Ainv;
  mat U, V;
  cube J, K, H, W;
  mat var;

  double Gn[3][3], Gv[3][3], C1[3][3];
  double eps;

  double M[3][3][3][3], e2[3][3], e4[3][3][3][3];

  double Cn, Cv, Ceps1, Ceps2, Ceps3;

  // Memeber functions
public:
  void initialHookScalarRansTurbModel(double Stau, double *struc, double &eps_init,
  		                                double *prod, double *rapRedi, double (*G)[3]);

  void calcReStress(double *struc, double &eps_main, double *prod, double *rRedi,
  		              double *sRediEps, double (*G)[3], double (*Gph)[3], double (*Gp1)[3],
  		              double dt, char* tIntName);

  void fwEuler(double dt, double (*G)[3]);

  void Heun(double dt, double (*G)[3]);

  void CrankN(double dt, double (*G)[3]);

  void RK4(double dt, double (*G)[3]);

  virtual void calcRhs(mat &rhs, mat &var, double dt);

  virtual void inputs(double (*G)[3]);

  virtual void outputs();

  virtual double calcRhsEps(double (*G)[3], double dt);

  void writeData(double St);

};

#endif

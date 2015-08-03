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

  double eps, tau;
  double M[3][3][3][3], e2[3][3], e4[3][3][3][3];

  double Cn, Cv, Ceps1, Ceps2, Ceps3;

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

  void calcTurbStatistics();

  virtual void driftCoeff(double *drift, double (*G)[3], double *eref,
                          double *cref);

  virtual double rhsDissipation(double (*Gn)[3], double dt);

  virtual void calcTurbTimeScale();

};

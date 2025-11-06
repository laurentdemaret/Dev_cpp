#ifndef _NL_OPT_H_
#define _NL_OPT_H_

#include "matrix.h"
#include "symbolic.h"
#include "numeric.h"
#include "utils.h"
#include <chrono>
#include <fftw3.h>
//#include "tiffio.h"

typedef std::complex<double> Complex;
typedef std::valarray<Complex> ComplexArray;

using namespace std;
using namespace GiNaC;


// Class to solve Least Square problems
class NL_LeastSquares
{
public:
  //Default constructor and destructor
  NL_LeastSquares();
  virtual ~NL_LeastSquares();
  NL_LeastSquares(vector<string> equations,vector<string> var_names,vector<string> param_names);

  void Initialise(vector<string> equations,vector<string> var_names,vector<string> param_names);

  void SetToHomography();
  //void SetToHomography_3D();

  void DumpSolution(bool verbose=false);
  void LM_Minimisation();
  void SetGroundTruth(M3Matrix Xv, vector<double> p_gt_val);
  void SetGroundTruthParameters(vector<double> p_gt_val);
  void SetInput(M3Matrix Xv);
  void SetObservations(vector<double> y_obs_v);
  void SetInitialGuess(vector<double> p0_val);

  void SimulateOutput(double sigma=0, int seed=1); //computes Y from f, p_gt
  void SimulateOutput(M3Matrix& X_val, vector<double>& p_gt_val); //set ground truth and compute the output

  double ReprojectionError();

  //members:
  Function_Multivariate f; //for a least square multivariate function
  M3Matrix X; //input
  //M3Matrix Y; //measurements from X
  vector<double> p_gt; //ground truth parameters
  vector<double> p_sol; //computed (solution) parameters output by the minimization algorithm
  vector<double> p0; //initial values

  vector<double> y_obs; //observed data
  vector<double> y_gt; //ground truth if available
  vector<double> y_sol; //computed solution by using the computed p_sol

  bool flag_gt_known = false;
  int NbIt = 10;
};


vector<double> LM_Solve(Function_Multivariate f, M3Matrix x_val,  vector<double> p0, vector<double> y_obs, int NbIt = 10, bool verbose = false);


//**************************
//**** Tests ***************
//**************************
void test_homography_minimization();
void test_LM();
void test_LM2();
void test_LM3();
void test_LM4();
#endif // of _NL_OPT_H_

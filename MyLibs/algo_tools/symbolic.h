#ifndef _SYMBOLIC_H_
#define _SYMBOLIC_H_

#include <stdlib.h>
#include <cstring>
#include <ginac/ginac.h>
#include <ginac/excompiler.h>
#include <cstdio>
#include <cstdlib>
#include <math.h>
#include "matrix.h"
#include "utils.h"

using namespace std;
using namespace GiNaC;

//This file gathers wrapper classes for symbolic calculus based on the ginac library

class Function_Univ
{
public:
  Function_Univ(); //Default constructor
  Function_Univ(string eq_str, string v_str); //Constructor from equation and variable names

  ~Function_Univ(){} //Default destructor

  //Initialisation
  void init(string eq_str, string v_str);

  //Modifiers
  void set(string eq_str,string v_str);
  void set_equation(string eq_str)
  {
    equation_str = eq_str;
  }
  void set_varnames(string v_str)
  {
    var_str = v_str;
  }
  void set_ex(ex expr)
  {
    expression = expr;
  }

  //Compilation of FUNC_1P
  void compile();
  //Evaluation tool
  double eval( double x){return equation_p(x);}

  Function_Univ diff(); //partial derivative according to the first variable

  //These are the symbols and the corresponding variable names
  symbol symbol1;
  string var_str;

  string equation_str;
  ex expression;
  FUNCP_1P equation_p;
};


//Wrapper for a bivarate function based on ginac
//It includes ex and FUNC_2P ginac objects:
//and allows for facilitated use
class Function_Biv
{
public:
  Function_Biv(); //Default constructor
  Function_Biv(string eq_str, string v1_str,string v2_str); //Constructor from equation and variable names

  ~Function_Biv(){} //Default destructor

  //Initialisation
  void init(string eq_str, string v1_str,string v2_str);

  //Modifiers
  void set(string eq_str,string v1_str,string v2_str);
  void set_equation(string eq_str)
  {
    equation_str = eq_str;
  }
  void set_varnames(string v1_str,string v2_str)
  {
    var1_str = v1_str;
    var2_str = v2_str;
  }

  void set_ex(ex expr)
  {
    expression = expr;
  }

  //I/O functionalitiess
  void WriteMatrixValuesForMATLAB(string& filenameOut);

  //Dump functionalities
  void dump_equation(){std::cout << equation_str << std::endl;}
  void dump_variables(){cout << var1_str << " " <<  var2_str<< endl;}

  Function_Biv diff(string var_str);
  Function_Biv diff1(); //partial derivative according to the first variable
  Function_Biv diff2(); //partial derivative according to the first variable

  vector<Function_Biv> gradient();
  vector<vector<Function_Biv>> hessian();

  //Compilation of FUNC_2P
  void compile();
  //Evaluation tool
  double eval( double x, double y){return equation_p(x,y);}
  //These are the symbols and the corresponding variable names
  symbol symbol1, symbol2;
  string var1_str, var2_str;

  string equation_str;
  ex expression;
  FUNCP_2P equation_p;
};


//******************************************************
//** class for multivariate functions ******************
//** basically a wrapper of expression and parameters **
//******************************************************
class Function_Multivariate
{
public:
  Function_Multivariate();
  Function_Multivariate(vector<string> equations,vector<string> var_names,vector<string> param_names);
  virtual ~Function_Multivariate();

  void setVariables(vector<string> var_names);
  void setParameters(vector<string> param_names);
  void setEquations(vector<string> equations);
  void Initialise(vector<string> equations,vector<string> var_names,vector<string> param_names);

  M3Matrix JacobiAnalytic(M3Matrix x_val, vector<double> p_val);//Jacobi
  vector<double> evalf(M3Matrix x_val, vector<double> p_val={});
  vector<double> operator()(M3Matrix x_val, vector<double> p_val={});

  int GetNbParams();
  int GetNbVars();
  int GetNbOutputSize();

  //Members
  vector<ex> f; // f_1,f_2,...
  vector<symbol> vars; //variables: x_1,...,x_m
  vector<symbol> params; //parameters: p_1, ... , p_n
};


// TODO: write a small wrapper class for multivariate functions of the type
// f(x_1, ..., x_n, p1, ..., p_n)
// p_1, ..., p_n are parameters to be optimized

//Tests for symbolic computations : Function_Biv class
void test_Function_Biv();
void test_ginac();
void test_ginac2();
void test_ginac_composition();
void test_function_multivariate();

#endif // of _SYMBOLIC_H_

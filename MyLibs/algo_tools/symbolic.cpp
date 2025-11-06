#include "symbolic.h"

//*************************************
//*************************************
// Function_Univ class *****************
//*************************************
//*************************************
//Constructors
//Default constructor
Function_Univ::Function_Univ(){}

// Constructor with dimensions
Function_Univ::Function_Univ(string eq_str, string v_str)
{
  var_str = v_str;
  equation_str = eq_str;
  compile();
}


void Function_Univ::init(string eq_str, string v_str)
{
  set(eq_str,v_str);
  compile();
}

void Function_Univ::set(string eq_str,string v_str)
{
  set_equation(eq_str);
  set_varnames(v_str);
}


void Function_Univ::compile()
{
  symbol1.set_name(var_str);

  //ex expression_tmp = ex(equation_str, lst(symbol1));
  ex expression_tmp = ex(equation_str, symbol1);
  expression = expression_tmp;
  compile_ex(expression,symbol1,equation_p);
}


Function_Univ Function_Univ ::diff() //derivative
{
  Function_Univ deriv;
  deriv.set_varnames(var_str);
  deriv.set_ex(expression.diff(symbol1));

  compile_ex(deriv.expression,symbol1,deriv.equation_p);
  return deriv;
}



//*************************************
//*************************************
// Function_Biv class *****************
//*************************************
//*************************************
//Constructors
//Default constructor
Function_Biv::Function_Biv() {}
// Constructor with dimensions
Function_Biv::Function_Biv(string eq_str, string v1_str,string v2_str)
{
  var1_str = v1_str;
  var2_str = v2_str;
  equation_str = eq_str;
  compile();
}


void Function_Biv::init(string eq_str, string v1_str,string v2_str)
{
  set(eq_str,v1_str,v2_str);
  compile();
}


void Function_Biv::set(string eq_str,string v1_str,string v2_str)
{
  set_equation(eq_str);
  set_varnames(v1_str,v2_str);
}

vector<Function_Biv> Function_Biv::gradient()
{
    vector<Function_Biv> gradient_tmp{this->diff1(),this->diff2()};
    return gradient_tmp;
}

vector<vector<Function_Biv>> Function_Biv::hessian()
{
    vector<vector<Function_Biv>> hessian_tmp{
        {this->diff1().diff1(),this->diff1().diff2()},
        {this->diff2().diff1(),this->diff2().diff2()}
    };
    return hessian_tmp;
}

Function_Biv Function_Biv::diff(string var_str)
{
  if(var_str == var1_str)
  {
    return diff1();
  }
  //else if(var_str == var2_str)
  else
  {
    return diff2();
  }
}


Function_Biv Function_Biv::diff1() //partial derivative according to the first variable
{
  Function_Biv partial_deriv;
  std::stringstream eq_diffstream;
  eq_diffstream << expression.diff(symbol1);
  //std::string eq_diff_str = eq_diffstream.str();
  partial_deriv.init(eq_diffstream.str(),var1_str,var2_str);
  /*partial_deriv.set_varnames(var1_str,var2_str);
  partial_deriv.symbol1.set_name(var1_str);
  partial_deriv.symbol2.set_name(var2_str);
  partial_deriv.set_ex(expression.diff(symbol1));
  //string ex_diff(expression.diff(symbol1));
   //ex expression_tmp = ex(equation_str, lst(symbol1,symbol2));
  //compile_ex(expression,symbol1,symbol2,equation_p);
  FUNCP_2P equation;
  compile_ex(partial_deriv.expression,symbol1,symbol2,equation);
  partial_deriv.equation_p =equation;*/

  return partial_deriv;
}


Function_Biv Function_Biv::diff2() //partial derivative according to the second variable
{
  Function_Biv partial_deriv;
  std::stringstream eq_diffstream;
  eq_diffstream << expression.diff(symbol2);
  std::string eq_diff_str = eq_diffstream.str();
  partial_deriv.init(eq_diff_str,var1_str,var2_str);

  //compile_ex(partial_deriv.expression,symbol1,symbol2,partial_deriv.equation_p);
  return partial_deriv;
}




void Function_Biv::compile()
{
  symbol1.set_name(var1_str);
  symbol2.set_name(var2_str);

  lst l = {symbol1,symbol2};
  //ex expression_tmp = ex(equation_str, lst(symbol1,symbol2));
  ex expression_tmp = ex(equation_str, l);
  expression = expression_tmp;
  compile_ex(expression,symbol1,symbol2,equation_p);
}

// ****************************************
// * Write the matrices on external files
// ****************************************
void Function_Biv::WriteMatrixValuesForMATLAB(string& filenameOut)
{
   //Default
    ofstream OutStr(filenameOut.c_str());
    double x, y;
    double max1 = 2.0;
    double max2 = 2.0;
    double NbX = 40;
    double NbY = 40;

    //10 als Parameter eingeben
   for(int i = 0; i < NbX; i++)
   {
       x = double(i)/ max1;
       for(int j = 0; j < NbY; j++)
       {
          //stdout
           y = double(j)/ max2;
           OutStr <<  (*this).eval(x,y) <<",";
       }
       OutStr << std::endl;
   }
}

//***************************************
//* FUNCTION_MULTIVARIATE
//***************************************
Function_Multivariate::~Function_Multivariate(){}
Function_Multivariate::Function_Multivariate(){}
Function_Multivariate::Function_Multivariate(vector<string> equations, vector<string> var_names, vector<string> param_names)
{
  Initialise(equations, var_names, param_names);
}

//Accessors
int Function_Multivariate::GetNbParams()
{
  return params.size();
}


int Function_Multivariate::GetNbVars()
{
  return vars.size();
}


int Function_Multivariate::GetNbOutputSize()
{
  return  f.size();
}


//Modifiers
void Function_Multivariate::setVariables(vector<string> var_names)
{
  if(vars.size() != var_names.size())
    vars.resize(var_names.size());

  for(int i = 0;i<int(var_names.size());i++)
  {
    vars[i].set_name(var_names[i]);
  }
}


void Function_Multivariate::setParameters(vector<string> param_names)
{
  if(params.size() != param_names.size())
      params.resize(param_names.size());

    for(int i = 0;i<int(param_names.size());i++)
    {
      params[i].set_name(param_names[i]);
    }
}


void Function_Multivariate::setEquations(vector<string> equations)
{
  lst symbol_list;
  //Append the symbols of the variables
  for(int i = 0; i < int(vars.size());i++)
    symbol_list.append(vars[i]);

  //Append the symbols of the parameters
  for(int i = 0; i < int(params.size());i++)
    symbol_list.append(params[i]);

  f.resize(equations.size());
  for(int i=0; i < int(equations.size());i++)
    f[i] = ex(equations[i],symbol_list);
}


void Function_Multivariate::Initialise(vector<string> equations,vector<string> var_names,vector<string> param_names)
{
  setVariables(var_names);
  setParameters(param_names);
  setEquations(equations);
}



vector<double> Function_Multivariate::operator()(M3Matrix x_val, vector<double> p_val)
{
  vector<double> y;
  y = evalf(x_val,p_val);
  return y;
}


//Function evaluation
vector<double> Function_Multivariate::evalf(M3Matrix x_val, vector<double> p_val)
{
  vector<double> y_val;

  cout << "evalf :: " << "x_val.GetNbRows() : " << x_val.GetNbRows() << endl;
  cout << "evalf :: " << "vars.size() : " << vars.size() << endl;
  cout << "evalf :: " << "p_val.size() : " << p_val.size() << endl;
  cout << "evalf :: " << "params.size() : " << params.size() << endl;

  if(x_val.GetNbRows() != int(vars.size()) || p_val.size() != params.size() )
  {
    cout << "JacobiAnalytic: something went wrong !!" << endl;
    exit(-1);
  }

  // This is the number of measured points : Each column is a measurement input
  int NbObs = x_val.GetNbCols();
  y_val.resize(NbObs*f.size());

  for(int i = 0;i<NbObs;i++) // along the observations
  {
    //compute value at the different measured points (get this out of the loop ? at least in part ? )
    lst symbol_value_list;
    for(int k = 0; k < int(vars.size()); k++)
      symbol_value_list.append(vars[k]==x_val[k][i]);
    for(int k = 0; k < int(params.size()); k++)
      symbol_value_list.append(params[k]==p_val[k]);

    for(int e=0; e < int(f.size());e++) //along the dimension of the output vector
    {
      ex tmp = f[e].subs(symbol_value_list);
      y_val[i*f.size()+e] = ex_to<numeric>(tmp).to_double();
   }
  }
  return y_val;
}


M3Matrix Function_Multivariate::JacobiAnalytic(M3Matrix x_val, vector<double> p_val)//Jacobi
{
  if(x_val.GetNbRows() != int(vars.size()) || p_val.size() != params.size() )
  {
    cout << "JacobiAnalytic: something went wrong !!" << endl;
  }

  // This is the number of measured points : Each column is a measurement input
  int NbObs = x_val.GetNbCols();

  //cout << "Jacobi Analytic" << endl;
  M3Matrix Jf(NbObs*int(f.size()),int(params.size()));
  ex dfdp; // temporary variable to store the derivatives

  for(int j = 0;j<Jf.GetNbCols();j++) //Along the parameters
  {
    for(int e=0; e < int(f.size());e++) //along the dimension of the output vector
    {
      //Symbolic derivation
      dfdp = f[e].diff(params[j]);
      //cout << "dfdp_" << j << "_" << e <<" : " << dfdp << endl;

      for(int i = 0;i<NbObs;i++) // along the observations
      {
        //compute value at the different measured points (get this out of the loop ? at least in part ? )
        lst symbol_value_list;
        for(int k = 0; k < int(vars.size()); k++)
          symbol_value_list.append(vars[k]==x_val[k][i]);
        for(int k = 0; k < int(params.size()); k++)
          symbol_value_list.append(params[k]==p_val[k]);

        ex tmp = dfdp.subs(symbol_value_list);

        Jf[i*f.size()+e][j] = ex_to<numeric>(tmp).to_double();
     }
    }
  }

  return Jf;
}



//*************************************
//*************************************
// Tests ******************************
//*************************************
//*************************************
double MySqr(double x)
{
  return sqrt(fabs(x));
}



void test_ginac_composition()
{
  cout<< "GinaC minimal test for function composition" << endl;
  string myeqsq("z*z"), myeqsum("x*x+y");
  symbol x("x"), y("y"), z("z");

  lst lz ={z};
  lst l_xy ={x,y};
  ex f = ex(myeqsq,lz);
  ex g = ex(myeqsum,l_xy);
  cout << "f : " << f << endl;
  cout << "g : " << g << endl;

  lst list;
  list.append(z==g);
  ex h = ex(f.subs(list));

  symbol a("a");
  cout << "f.subs(list) : " << f.subs(list) << endl;
  cout << "h : " << h << endl;
  cout << "h : " << h.diff(x) << endl;

  //Image formation model
  symbol h00("h00"), h01("h01"), h02("h02"), tx("tx"),
          h10("h10"), h11("h11"), h12("h12"), ty("ty"),
          h20("h20"), h21("h21"), h22("h22"), tz("tz");
  //matrix M(3,3);
  matrix M = {{h00, h01, h02},
                  {h10, h11, h12},
                  {h20, h21, h22}};

  cout << M << endl;
  symbol x1("x1"), x2("x2"),x3("x3");

  matrix X = {{x1},{x2},{x3}};

  matrix Y = M.mul(X);
  cout << "Y : " << Y << endl;

  //cout << " Y.diff(x1) : "<<Y.diff(h00) << endl;

  //Calibration matrix
  symbol fd("f"), u("u"), v("v");

  //matrix K(3,3);
  matrix K = {{fd,0,u},{0,fd,v},{0,0,1}};

  //XS: X on the sensor
  matrix XS = K.mul(Y);
  cout << "XS: " << XS << endl;

}


void test_ginac()
{
  //Minimal Ginac test program to test the parsing, symbolic calculus and numeric evaluation capacities of ginac
  cout<< "GinaC TEST" << endl;

  string myeq("x+y*y"), myeqsph("cos(theta) + phi *phi *phi");
  symbol x("x"), y("y"), theta("theta"), phi("phi");

  //test pour une formule développée soi-même
  /*string myeq2("MySqrt(x)");
  cout <<"string"<< endl;
  ex f_test = ex(myeq2,x);
  cout <<"ex"<< endl;
  FUNCP_1P fp_test;
  cout <<"FUNCT_1P"<< endl;
  compile_ex(f_test, x, fp_test); */

  /*ex e = 2*MySqr(x) +3;
  FUNCP_1P fp_test;
  cout <<"FUNCT_1P"<< endl;
  compile_ex(e, x, fp_test);
  cout << fp_test(2.5) << endl; */

  //fin du test

  lst l_xy={x,y};
  lst l_theta_phi={theta,phi};
  ex f1 = ex(myeq,l_xy);
  ex g1 = f1.diff(x);
  ex f2 = ex(myeqsph,l_theta_phi);
  ex g2 = (f2.diff(theta) + 3*theta*theta).diff(theta);

  cout << "f1 : " << f1 << endl;
  cout << "f2 : " << f2 << endl;
  cout << "g2 : " << g2 << endl;

  //Univariate function
  string expr_str("x*x");

  lst lx ={x};
  ex f = ex(expr_str,lx);

  FUNCP_1P fp;
  compile_ex(f, x, fp);
  cout << fp(2.5) << endl;

  FUNCP_2P gp;
  compile_ex(f1,x,y,gp);
  double xv = 3.5,yv = 1.2;
  cout <<"For x =" << xv << " , y = " << yv << ", "<< f1 << ":" << gp(xv,yv) << endl;

  compile_ex(g1,x,y,gp);
  cout <<"For x =" << xv << " , y = " << yv << ", "<< g1 << ":" << gp(xv,yv) << endl;
}

void test_Function_Biv()
{
  Function_Biv  f("x*x*y + 0.5*x*x + x*y+x*y*y +y^3","x","y");
  std::cout << "this is the function : " << endl;
  std::cout<< f.expression << endl;
  std::cout << "f(2,3): " << f.eval(2,3) << endl;
  Function_Biv f_x;
  f_x = f.diff1();
  std::cout << "this is the partial derivative f_x : " <<  std::endl << f_x.expression << std::endl;
  std::cout<< "... and its evaluation at the point (2,3)" << std::endl;
  std::cout << "f(2,3): " << f_x.eval(2,3) << std::endl;
  Function_Biv f_xx;
  f_xx = f_x.diff1();
  std::cout<< "f_xx : " << std::endl;
  std::cout<< f_xx.expression << endl;
  Function_Biv f_yy;
  f_yy = f.diff2().diff2();
  std::cout<< "f_yy : " << std::endl;
  std::cout<< f_yy.expression << endl;

  string filenameOut("f.txt");
  f.WriteMatrixValuesForMATLAB(filenameOut);
  //std::stringstream a;
  //a << f_x.expression  << endl;
  //std::string a_str = a.str();
  //std::cout << "a_str : "  << a_str << std::endl;
  /*Function_Biv myfunc("sin(theta)+3*theta*phi*phi","theta","phi");
  cout << "this is the function : " << endl;
  cout<< myfunc.expression << endl;
  cout << "f(1,2): " << myfunc.eval(1,2) << endl;*/

  /* Function_Biv myfunc_theta = myfunc.diff("theta");
  cout << "its partial derivative along theta : " << endl;
  cout<< myfunc_theta.expression << endl;

  Function_Biv myfunc_phi = myfunc.diff("phi");
  cout << "its partial derivative along phi : " << endl;
  cout << myfunc_phi.expression << endl;

  std::cout  << std::endl;
  std::cout  << std::endl;
  std::cout  << std::endl;*/

  //*****************************+
  Function_Biv myfunc2("x^2+x*y+y^3","x","y");
  std::cout << "this is the second function : " << endl;
  std::cout<< myfunc2.expression << endl;

  vector<Function_Biv> gradient2;
  gradient2 = myfunc2.gradient();
  std::cout << "this is the gradient" << std::endl;
  std::cout << "(" << gradient2.at(0).expression << ", " << gradient2.at(1).expression << ")" << std::endl;

  Function_Biv f2_x = myfunc2.diff1();
  /*Function_Biv f2_xx = f2_x.diff1();
  std::cout << "this is f2_xx" << std::endl;
  std::cout << f2_xx.expression << std::endl;*/

  vector<vector<Function_Biv>> hessian2;
  hessian2 = myfunc2.hessian();

  //Write the hessian matrix components, each in a file
  string filename_fx("f_x.txt");
  string filename_fy("f_y.txt");
  myfunc2.diff1().WriteMatrixValuesForMATLAB(filename_fx);
  myfunc2.diff2().WriteMatrixValuesForMATLAB(filename_fy);

   //Write the hessian matrix components, each in a file
  string filename_fxx("f_xx.txt");
  string filename_fxy("f_xy.txt");
  string filename_fyy("f_yy.txt");
  hessian2.at(0).at(0).WriteMatrixValuesForMATLAB(filename_fxx);
  hessian2.at(0).at(1).WriteMatrixValuesForMATLAB(filename_fxy);
  hessian2.at(1).at(1).WriteMatrixValuesForMATLAB(filename_fyy);

  std::cout << "this is the hessian" << std::endl;
  //std::cout << (myfunc2.diff1()).diff1().expression << std::endl;
  std::cout << "("  << std::endl << hessian2.at(0).at(0).expression << "   " << hessian2.at(0).at(1).expression << " " << std::endl;
  std::cout << "  "  << hessian2.at(1).at(0).expression << "  " << hessian2.at(1).at(1).expression << " " << std::endl;
  std::cout  <<  ")" << std::endl;
}



void test_function_multivariate()
{
  /*M3Matrix A(2,3);
  A.SetValues({1,2,3,
               4,5,6}
  );
  A.Dump();*/

  //Test JacobiAnalytic
  int NbParams = 2; // NbVariables = 2;
  int NbObs = 3; // number of observations
  //Values
  vector<double> p_val(NbParams);
  p_val[0] = 1.;
  p_val[1] = 2.;

  M3Matrix x_val(2,NbObs);
  x_val[0][0] = 1;
  x_val[1][0] = 1;

  x_val[0][1] = 2;
  x_val[1][1] = 3;

  x_val[0][2] = 3;
  x_val[1][2] = 2;


  //Test the class Function_Multivariate
  Function_Multivariate f({"a1*x1^2 + 0.25*a2*a2*x2"},{"x1","x2"},{"a1","a2"});

  M3Matrix Jf2;
  Jf2 = f.JacobiAnalytic(x_val, p_val);
  cout << "Jf2.Dump()" << endl;
  Jf2.Dump();

  //Test of evalf
  Function_Multivariate f2({"x+y","x-y"},{"x","y"},{});
  M3Matrix X(2,3,{1,2,3,4,5,6});

  vector<double> f_val;
  f_val = f2.evalf(X);
  dump_vec(f_val,"f2.eval(X) : ");
  //dump_vec(c,"c");

  //Test block
  /* symbol x("x"), y("y"), z("z");
  //string linear_eq2("x+y+z");
  //ex f2 = ex(linear_eq2,lst(x, y, z));
  ex f2 = x+y+2*z;

  double x_val = 1;
  double y_val = 2.1;
  double z_val = 3.2;

  lst list;
  list.append(x==x_val);
  list.append(y==y_val);
  list.append(z==z_val);

  //variables.push_back(x1_val);

  cout << "f2.subs(lst(x == 1, y == 2.1, z == 3.2)) : "
         << f2.subs(list) << endl;
  //cout << "f2.subs(lst(x == 1, y == 2.1, z == 3.2)) : "
       //<< f2.subs(lst(x == 1, y == 2.1, z == 3.2)) << endl;
  //x = numeric(0.2);


  //cout << f2.subs(x==1.0).subs(y==2.0).subs(z==3.275)<< endl;
  cout << "f2 : " << f2 << endl; */
  //End of Test block
}


void test_ginac2()
{
    int example_indx = 1;
    if(example_indx == 0)
    {
      //Below minimal working example
      symbol x("x"), y("y");
      ex poly;

      for (int i=0; i<3; ++i)
          poly += factorial(i+16)*pow(x,i)*pow(y,2-i);

      cout << poly << endl;
    }
    else if(example_indx == 1)
    {
      symbol x("x");
      symbol y("y");
      //string expr_str = "3*pow(x,3)+2*x*y+2*sin(x) + (x+y)/(x-y)";
      string expr_str = "1/(x-y)";
      //ex e = 3*pow(x,3)+2*y*x+sin(x);
      lst l_xy ={x,y};
      ex e = ex(expr_str,l_xy);
      cout << "The partial derivative along x of the function " << expr_str << " is : " << endl;
      cout << e.diff(x) << endl;
      cout << "The partial derivative along y of the function " <<  expr_str << "is : " << endl;
      cout << e.diff(y) << endl;
    }
}

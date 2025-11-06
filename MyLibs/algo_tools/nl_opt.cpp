#include "nl_opt.h"

using namespace std::chrono;


// Levenberg-Marquardt non-linear optimization
// Input parameters:
// Function_Multivariate f
// initial parameter p0
// This is the standard LM algo.
// Here we follow the lecture script
// Einfuehrung in die Optimierung von Harrach
// Uni. Stuttgart, WS 13/14 - www.mathematik.uni-stuttgart.de/oip
vector<double> LM_Solve(Function_Multivariate f, M3Matrix x_val,  vector<double> p0, vector<double> y_obs, int NbIt, bool verbose)
{
    //Variables
    double beta_0 = 0.3, beta_1 = 0.9, mu = 0.01, tol = 0.0001;
    double err_current, err_new, epsilon = 0.;
    vector<double> p_current, p_new, r, r_new, y_val, JtF_vec;
    M3Matrix J, Jt, H, Jt_r, d;

    p_current = p0;
    int it;
    for(it = 0; it < NbIt; it++)
    {
        cout << "LM_Solve:: iter : " << it << endl;
        epsilon = 0.;
        //while(epsilon < beta_0)
        //{
        //J = \nabla_p f the Jacobian of the function f relatively to the parameter variables p
        J = f.JacobiAnalytic(x_val, p_current);

        if(verbose)
        {
            cout << "Jacobi matrix" << endl;
            J.Dump();
        }

        Jt = Transpose(J);
        if(verbose)
        {
            cout << "its transpose" << endl;
            Jt.Dump();
        }

        H = Jt * J; //Approximate hessian

        if(verbose)
        {
            cout << "... and the approximate hessian" << endl;
            H.Dump();
        }

        //reweighting of the diagonal term of the approximate Hessian
        for(int i = 0;i<H.GetNbRows();i++)
            H[i][i] *= (1.+mu);

        if(verbose)
        {
            cout << "after reweighting" << endl;
            H.Dump();
        }

        // Computes the residual r_i
        y_val = f.evalf(x_val, p_current);
        if(verbose)
            dump_vec(y_val,"y_val\n");

        r = y_obs - y_val;
        if(verbose)
            dump_vec(r,"r\n");

        JtF_vec = Transpose(J) * r; //r = y_obs -y_val
        if(verbose)
            dump_vec(JtF_vec,"JtF_vec \n");

        M3Matrix JtF(JtF_vec);
        if(verbose)
        {
            cout << "JtF.Dump()" << endl;
            JtF.Dump();
        }

        CholeskySolve(H, d, JtF);
        if(verbose)
        {
            cout << "d.Dump()" << endl;
            d.Dump();
        }

        //update the parameters
        p_new = p_current + d;
        if(verbose)
        {
            dump_vec(p_current,"p_current \n");
            dump_vec(p_new,"p_new \n");
        }

        //vector<double> fxk = f.evalf(x_val,p_current);
        vector<double> y_val_new = f.evalf(x_val,p_new);
        r_new = y_obs - y_val_new;

        err_current = Norm(r);
        err_new = Norm(r_new);
        if(err_new < tol)
        {
            cout << "err_new : " << err_new << " lower than " << tol << endl;
            break;
        }

        if(verbose)
        {
            dump_vec(r_new,"r_new \n");
            dump_vec(y_val_new,"y_val_new \n");
        }

        //vector<double > r_model = r - J*(p_new-p_current); //is this really correct ??? cf. p.58 in von Harrach script !!
        //dump_vec(r_model,"r_model \n");
        if(verbose)
        {
            cout << "error : " << err_current << endl;
            cout << "new error : " << err_new << endl;
        }

        //cout << "err_model : " << Norm(r_model) << endl;
        epsilon = (err_current*err_current - err_new*err_new)/(err_current*err_current);
        //epsilon /= (Norm(r)*Norm(r) - Norm(r_model)*Norm(r_model));
        if(verbose)
        {
            cout << "epsilon : " << epsilon << endl;
        }

        /*if(epsilon < beta_0)
      {
        mu/=1.75; //?????? this works with the inverse
      }*/
        cout << "error at iteration "<< it << " : " << err_new << endl;
        //}

        //Update parameters
        p_current = p_new;
        if(err_new < tol)
        {
            cout << "err_new : " << err_new << " lower than " << tol << endl;
            break;
        }

        if(epsilon > beta_1)
        {
            //mu/=2.;
            mu*=2; //?????? this works with the inverse
        }
    }

    cout << endl;
    cout << "LM_Solve used " << it << " iterations" << endl;
    return p_current;
}


//***********************************************
//* Class Non Linear Least Squares Minimization
//***********************************************
NL_LeastSquares::NL_LeastSquares(){}
NL_LeastSquares::~NL_LeastSquares(){}

NL_LeastSquares::NL_LeastSquares(vector<string> equations,vector<string> var_names,vector<string> param_names)
{
    Initialise(equations, var_names, param_names);
}


void NL_LeastSquares::Initialise(vector<string> equations,vector<string> var_names,vector<string> param_names)
{
    f.Initialise(equations, var_names, param_names);
}


// Classical notations - see for instance :
// http://www.corrmap.com/features/homography_transformation.php
void NL_LeastSquares::SetToHomography()
{
    vector<string> var_names({"x","y"});
    vector<string> param_names({"a","b","tx","d","e","ty","g","h","tz"});
    vector<string> equations({"(a*x+b*y+tx)/(g*x+h*y+tz)","(d*x+e*y+ty)/(g*x+h*y+tz)"});
    Initialise(equations,var_names,param_names);
}


//Notations like in paper
// Zhang's camera calibration algorithm:
// In-Depth Tutorial and Implementation
// WARNING !!!!
// This is the situation where all the input points have Z coordinate Z=0 (plane target)
/*void NL_LeastSquares::SetToHomography_3D_Plane()
{
  vector<string> var_names({"X","Y","Z"});
  //vector<string> param_names({"rxx","rxy","rxz","ryx","ryy","ryz","rzx","rzy","rzz","tx","ty","tz"});
  //vector<string> equations({"(rxx*X+rxy*Y+rxz*Z+tx)/(rzx*X+rzy*Y+rzz*Z+tz)","(ryx*X+ryy*Y+ryz*Z+ty)/(rzx*X+rzy*Y+rzz*Z+tz)"});
  vector<string> param_names({"rxx","rxy","rxz","ryx","ryy","ryz","tx","ty","tz"});
  vector<string> equations({"(rxx*X+rxy*Y+tx)/(rzx*X+rzy*Y+tz)","(ryx*X+ryy*Y+ty)/(rzx*X+rzy*Y+tz)"});
  Initialise(equations,var_names,param_names);
}*/



//Use of the Levenberg-Marquardt algorithm to solve
void NL_LeastSquares::LM_Minimisation()
{
    cout << "LM_Minimisation : X:" << endl;
    X.Dump();

    cout << "LM_Minimisation : y_obs:" << endl;
    dump_vec(y_obs);

    p_sol = LM_Solve(f, X, p0, y_obs, NbIt);
    y_sol = f(X,p_sol);
}


void NL_LeastSquares::DumpSolution(bool verbose)
{
    dump_vec(p_gt,"ground truth parameters: \n");
    dump_vec(p_sol,"computed parameters: \n");
    dump_vec(y_obs,"observed function values: \n");
    dump_vec(y_sol,"computed function values: \n");

    if(verbose)
    {
        cout << "Parameter norm of the diff:" << Norm(p_gt - p_sol) << endl;
        cout << "reprojection error:" << Norm(y_gt - y_sol) << endl;
    }
}


//computes (ground truth) Y from f, p_gt and the input points
void NL_LeastSquares::SimulateOutput(M3Matrix& X_val, vector<double>& p_gt_val)
{
    SetGroundTruth(X_val, p_gt_val);
    y_gt = f(X,p_gt);
}


//computes (ground truth) Y from f, p_gt and the input points
void NL_LeastSquares::SimulateOutput(double sigma, int seed)
{
    y_gt = f(X,p_gt);

    //Here possibly add simulated noise
    y_obs = y_gt + sigma *randn(y_gt.size(),seed);
}


//Reprojection Error
double NL_LeastSquares::ReprojectionError()
{
    double reproj_err;
    reproj_err = Norm(y_sol-y_obs);
    return reproj_err;
}



void NL_LeastSquares::SetInitialGuess(vector<double> p0_v)
{
    p0 = p0_v;
}


void NL_LeastSquares::SetGroundTruth(M3Matrix Xv, vector<double> p_gt_val)
{
    flag_gt_known = true;
    X = Xv;
    p_gt = p_gt_val;
}



void NL_LeastSquares::SetGroundTruthParameters(vector<double> p_gt_val)
{
    flag_gt_known = true;
    p_gt = p_gt_val;
}


void NL_LeastSquares::SetInput(M3Matrix Xv)
{
    X = Xv;
}

void NL_LeastSquares::SetObservations(vector<double> y_obs_v)
{
    y_obs = y_obs_v;
}


//*************************
// Test routines
//*************************

void test_homography_minimization()
{
    NL_LeastSquares NL;
    NL.SetToHomography();

    int NbObs = 6;
    M3Matrix X(NbObs,2);

    X={0,0,
       0,1,
       0,2,
       1,0,
       1,1,
       1,2
      };
    X.Transpose();

    cout << "x after transpose" << endl;
    X.Dump();

    vector<double> p_gt({1,0,1,0,1,2,0.8,1.2});
    dump_vec(p_gt,"ground truth parameters\n");

    vector<double> p0({0.5,-1.1,1.5,1,0.5,1.1,1.5,1.2});
    dump_vec(p0,"initial guess\n");

    NL.SetInitialGuess(p0);
    NL.SetGroundTruth(X,p_gt);

    //Generate (simulated) output taken as observations
    NL.SimulateOutput();

    //Compute the minimisation
    NL.NbIt=50;
    NL.LM_Minimisation();
    NL.DumpSolution(true);

    //TODO add some simulated noise to see the robustness ?
}


void test_LM() //This was verified and tested (9. november 2018): do not touch !!
{
    cout << "minimal test of the Levenberg-Marquardt algorithm" << endl;
    //Test the class Function_Multivariate
    Function_Multivariate f({"a1*x1 + a2*x2"},{"x1","x2"},{"a1","a2"});
    int NbObs = 3;
    M3Matrix X(f.GetNbVars(),NbObs);
    X={0, 1, -1,
       2, 0, 1
      };
    X.Dump();

    vector<double> p_gt({0.17,0.17}); //ground truth

    dump_vec(f(X,p_gt), "resultat \nr");

    //vector<double> p0({0.4,0.7}); //initial guess for the parameters
    vector<double> p0({0.34,0.75}); //initial guess for the parameters
    vector<double> y_obs;
    y_obs = f(X,p_gt); //note that in reality such measurements will be noisy
    dump_vec(y_obs,"observed data, generated from the gt:");

    //Now invoke LM algo
    vector<double> p_sol;
    p_sol = LM_Solve(f, X,  p0, y_obs, 20);
    dump_vec(p_gt, "GT Parameters \n");
    dump_vec(p_sol, "Parameters computed by LM_Solve\n");
}


void test_LM2() //This was verified and tested (9. november 2018): do not touch !!
{
    cout << "minimal test of the Levenberg-Marquardt algorithm" << endl;
    //Test the class Function_Multivariate
    Function_Multivariate f({"a1*sin(x1) + a2*cos(x2)*x2+a3*x1*x3"},{"x1","x2","x3"},{"a1","a2","a3"});
    int NbObs = 5;
    M3Matrix X(f.GetNbVars(),NbObs);
    X={0, 1, -1, 0,1,
       2, 0, 1, 3,4,
       0,0,0,1,1
      };
    X.Dump();

    vector<double> p_gt({0.3,0.5,0.7}); //ground truth

    dump_vec(f(X,p_gt), "resultat \nr");

    //vector<double> p0({0.4,0.7}); //initial guess for the parameters
    vector<double> p0({0.2,0.6,0.8}); //initial guess for the parameters
    vector<double> y_obs;
    y_obs = f(X,p_gt); //note that in reality such measurements will be noisy
    dump_vec(y_obs,"observed data, generated from the gt:");

    //Now invoke Levenberg-Marwuardt algorithm
    vector<double> p_sol;
    p_sol = LM_Solve(f, X,  p0, y_obs, 25);
    dump_vec(p_gt, "GT Parameters \n");
    dump_vec(p_sol, "Parameters computed by LM_Solve\n");
}


void test_LM3() //This was verified and tested (9. november 2018): do not touch !!
{
    cout << "minimal test of the Levenberg-Marquardt algorithm" << endl;
    //Test the class Function_Multivariate
    Function_Multivariate f({"a1*x1+a2*x2","a1*x1-a2*x2"},{"x1","x2"},{"a1","a2"});
    //Function_Multivariate f({"a1*sin(x1) + a2*cos(x2)*x2+a3*x1*x3","a1+a2*+a3+x1+x2+x3"},{"x1","x2","x3"},{"a1","a2","a3"});

    int NbObs = 5;
    M3Matrix X(f.GetNbVars(),NbObs);
    X={0, 1, -1, 0,1,
       2, 0, 1, 3,4
      };
    X.Dump();

    vector<double> p_gt({0.1,0.6}); //ground truth
    //vector<double> y;
    //y = f(X,p_gt);
    //dump_vec(f(X,p_gt), "ground truth values\n");

    //vector<double> p0({0.4,0.7}); //initial guess for the parameters
    vector<double> p0({0.3,0.4}); //initial guess for the parameters
    vector<double> y_obs;
    y_obs = f(X,p_gt); //note that in reality such measurements will be noisy
    dump_vec(y_obs,"output generated from the gt:");

    //Now invoke LM algo
    vector<double> p_sol;
    p_sol = LM_Solve(f, X,  p0, y_obs, 25);
    dump_vec(p_gt, "GT Parameters \n");
    dump_vec(p_sol, "Parameters computed by LM_Solve\n");
}


void test_LM4() //This is the homography example !!
// TODOs: lundi 12 novembre 2018:
// comprendre pourquoi les criteres pour adapter le mu ne correspondent pas a ceux donnes par von Harrach
{
    cout << "minimal test of the Levenberg-Marquardt algorithm - for a homography" << endl;

    vector<string> param_names({"h00","h01","h02","h10","h11","h12","h20","h21","h22"});
    vector<string> var_names({"X","Y"});
    vector<string> equations({"(h00*X+h01*Y+h02)/(h20*X+h21*Y+h22)","(h10*X+h11*Y+h12)/(h20*X+h21*Y+h22)","h00+h01+h02+h10+h11+h12+h20+h21+h22"});

    Function_Multivariate f(equations,var_names, param_names);
    int NbObs = 9;
    M3Matrix X(f.GetNbVars(),NbObs);
    X={0, 0, 0, 1, 1, 1, 2, 2, 2,
       0, 1, 2, 0, 1, 2, 0, 1, 2
      };
    X.Dump();

    vector<double> p_gt({1, 0.1,0,0,1,0,0,0,1}); //ground truth

    vector<double> p0({0.7, 0.05,-0.04,0.32,1.3,0,0,0.1,0.6}); //initial guess for the parameters
    vector<double> y_obs;
    y_obs = f(X,p_gt); //note that in reality such measurements will be noisy
    dump_vec(y_obs,"output generated from the gt:");

    //Now invoke LM algo
    vector<double> p_sol;
    p_sol = LM_Solve(f, X,  p0, y_obs, 1000);
    dump_vec(p_gt, "GT Parameters \n");
    dump_vec(p_sol, "Parameters computed by LM_Solve\n");
    vector<double> y_sol;
    y_sol = f(X,p_sol); //note that in reality such measurements will be noisy
    dump_vec(y_sol,"output generated from the computed solution:");
}


/*void test_LM5() //Homography example : to test
// TODOs: 20 april 2020
// TBD : PAS ENCORE COMMENCE !!!!!
{
    cout << "minimal test of the Levenberg-Marquardt algorithm - for a homography" << endl;

    vector<string> param_names({"h00","h01","h02","h10","h11","h12","h20","h21","h22"});
    vector<string> var_names({"X","Y"});
    vector<string> equations({"(h00*X+h01*Y+h02)/(h20*X+h21*Y+h22)","(h10*X+h11*Y+h12)/(h20*X+h21*Y+h22)","h00+h01+h02+h10+h11+h12+h20+h21+h22"});

    Function_Multivariate f(equations,var_names, param_names);
    int NbObs = 9;
    M3Matrix X(f.GetNbVars(),NbObs);
    X={0, 0, 0, 1, 1, 1, 2, 2, 2,
       0, 1, 2, 0, 1, 2, 0, 1, 2
      };
    X.Dump();

    vector<double> p_gt({1, 0.1,0,0,1,0,0,0,1}); //ground truth

    vector<double> p0({0.7, 0.05,-0.04,0.32,1.3,0,0,0.1,0.6}); //initial guess for the parameters
    vector<double> y_obs;
    y_obs = f(X,p_gt); //note that in reality such measurements will be noisy
    dump_vec(y_obs,"output generated from the gt:");

    //Now invoke LM algo
    vector<double> p_sol;
    p_sol = LM_Solve(f, X,  p0, y_obs, 1000);
    dump_vec(p_gt, "GT Parameters \n");
    dump_vec(p_sol, "Parameters computed by LM_Solve\n");
    vector<double> y_sol;
    y_sol = f(X,p_sol); //note that in reality such measurements will be noisy
    dump_vec(y_sol,"output generated from the computed solution:");
}*/




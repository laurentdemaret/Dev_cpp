

// Main inspired from the ceres documentation
// Author: sameeragarwal@google.com (Sameer Agarwal)
// Original notice
// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2015 Google Inc. All rights reserved.
// http://ceres-solver.org/
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * Neither the name of Google Inc. nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

// Modification : Author Laurent Demaret 2022-2023


#include <vector>

#include "ceres/ceres.h"
#include "glog/logging.h"

using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solve;
using ceres::Solver;

#define TWO_PI 6.2831853071795864769252866
#define PI     3.14159265358979323846264338

//TODO: see arc4random : better random generator ?
double randn(int seed)
{
    if(seed != -1)
        srand(seed);

    double rand1, rand2;
    rand1 = rand() / ((double) RAND_MAX);
    if(rand1 < 1e-100) rand1 = 1e-100;

    rand1 = -2 * log(rand1);
    rand2 = (rand() / ((double) RAND_MAX)) * TWO_PI;

    return sqrt(rand1) * cos(rand2);
}

//***********************
// Functor
//***********************
// this is a functor
struct affine {
  affine(double val) : x(val) {}  // Constructor
  double operator()(const double& a, const double& b, double& result) const {
      result = a*x+b;
      return x + a; }

private:
  double x;
};

void test_functor()
{
    // Now you can use it like this:
    affine affine_a(42.); // create an instance of the functor class
    double a = 2, b= -3;
    double result;
    affine_a(a,b,result); // and "call" i
    std::cout << "result of functor call affine_a(8) : " << result << std:: endl;
}
//assert(i == 50); // and it added 42 to its argument


//***********************
// Function
//***********************
struct LinRegResidual
{    
	LinRegResidual(double x, double y) : x_(x), y_(y) {}
    /*bool operator()( const double& a0, const double& a1,  double &  residual) const {
        double model_eval;
        model_eval = a1 * x_ - a0;
        residual = y_ - model_eval;
        return true;
      }*/
    //Alte Version ...
    template <typename T>
    bool operator()( const T* const& a0, const T* const& a1, T* residual) const {
      T model_eval;
      model_eval = a1[0] * x_ + a0[0];
      residual[0] = y_ - model_eval;
      return true;
    }

   private:
    double x_;
    double y_;
};


template<typename T>
T quadratic_model(double x, T a0, T a1, T a2)
{
    T y = a2* x * x - a1 * x - a0;
    return y;
}

struct QuadRegResidual
{    
	QuadRegResidual(double x, double y) : x_(x), y_(y) {}
    template <typename T>
    bool operator()(const T* const a0, const T* const a1, const T* const a2, T* residual) const
    {
      //double a0_ = a0[0]*1.0;
      //double a1_ = a1[0]*1.0;
      //double a2_ = a2[0]*1.0;

      //residual[0] = y_ - quadratic_model(x_,a0_,a1_,a2_);
        T y_mod =  quadratic_model(x_,a0[0],a1[0],a2[0]);
      //T y_mod =   a2[0] * x_ * x_ - a1[0] * x_ - a0[0];
      residual[0] = y_ - y_mod;

      //residual[0] = y_ - a2[0] * x_ * x_ - a1[0] * x_ - a0[0];
      return true;
    }

   private:
    const double x_;
    const double y_;
};


struct ExponentialResidual
{
  ExponentialResidual(double x, double y) : x_(x), y_(y) {}

  template <typename T>
  bool operator()(const T* const m, const T* const c, T* residual) const {
    residual[0] = y_ - exp(m[0] * x_ + c[0]);
    return true;
  }

 private:
  const double x_;
  const double y_;
};


struct LinBivariateResidual
{
  LinBivariateResidual(double x1, double x2, double y) : x1_(x1), x2_(x2),y_(y) {}

  template <typename T>
  bool operator()(const T* const a1, const T* const a2, T* residual) const {
    residual[0] = y_ - a1[0] * x1_ - a2[0] * x2_ ;
    return true;
  }

 private:
   const double x1_;
   const double x2_;
   const double y_;
};

struct HomographicResidual
{
  HomographicResidual(double x1, double x2, double x3, double y) : x1_(x1), x2_(x2), x3_(x3),y_(y) {}

  template <typename T>
  bool operator()(const T* const a1, const T* const a2, const T* const a3,
                  const T* const a4, const T* const a5, const T* const a6,
                                    T* residual) const {
      residual[0] = y_ - (a1[0] * x1_ + a2[0] * x2_ + a3[0] * x3_)/(1.+a6[0]*x3_); // /(a4[0] * x1_ + a5[0] * x2_ + a6[0] * x3_);
    return true;
  }

  private:
    const double x1_;
    const double x2_;
    const double x3_;
    const double y_;
};

void curve_fitting_examples()
{
    std::cout << "Debut de curve_fitting " << std::endl;	

    //Initial values for the parameters of the linear model
    double a0_init = 0.;
    double a1_init = 0.;
    double a0 = a0_init;
    double a1 = a1_init;

    //Initial values for the parameters of the linear model
    double b0_init = 0.;
    double b1_init = 0.;
    double b2_init = 0.;
  
    double b0 = b0_init;
    double b1 = b1_init;
    double b2 = b2_init;


    // Data 
    int nb_obs = 20;
    std::vector<double> xv(nb_obs),yv(nb_obs); 
  
    //double 
    for(int i = 0;i<nb_obs;i++)
    {
  	  xv[i] = (double)i;
  	  yv[i] = 0.3*xv[i]*xv[i] + 3.2*xv[i] + 1.3;
    }

    //**************************
    // * TEST
    //**************************

    //LinRegResidual lbr(xv[0],yv[0]);
    
    Problem problem1,problem2;
    for (int i = 0; i < nb_obs; i++) 
    {
  	  //Linear fitting model
  	  problem1.AddResidualBlock(new AutoDiffCostFunction<LinRegResidual, 1, 1, 1>
  			    (new LinRegResidual(xv[i], yv[i])),nullptr,  &a0, &a1);

  	  //Quadratic fitting model			  
  	  problem2.AddResidualBlock(new AutoDiffCostFunction<QuadRegResidual, 1, 1, 1, 1>
  			    (new QuadRegResidual(xv[i], yv[i])),nullptr,  &b0, &b1, &b2);
			  
    }

    Solver::Options options;
    options.max_num_iterations = 125;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = true;

    Solver::Summary summary1,summary2;
    Solve(options, &problem1, &summary1);
    Solve(options, &problem2, &summary2);
    std::cout << summary1.BriefReport() << "\n";
    std::cout << "Initial model : " << 0.0 << "x + " << 0.0 << "\n";
    if(a0 >0) 
  	  std::cout << "Estimated model(linear): " << a1 << "x+" << a0 << "\n";
    else
  	  std::cout << "Estimated model(linear): " << a1 << "x" << a0 << "\n";
    std::cout << "Estimated model(quadratic): " << b2 << "x^2+"<< b1 << "x+" << b0 << "\n";
  
    std::cout << "curve fitting termine" << std::endl;	
}


void linear_bivariate_example()
{
	//This example simulates the situation of a function R^2 to R y = a1*x1+a2*x2
	
	//Data generation 
    int SideSize = 5;
    int nb_obs =  SideSize * SideSize;
	std::vector<double> x1(nb_obs), x2(nb_obs);
	std::vector<double> y(nb_obs);
	
	
	double a1_gt= 1.,a2_gt = 2.4;
	for(int i = 0;i<SideSize;i++)
		for(int j = 0;j<SideSize;j++)
     	{
			std::cout << "i*SideSize+j : " << i*SideSize+j << std::endl;
			x1[i*SideSize+j] = (double)i;
			x2[i*SideSize+j] = (double)j;
			y[i*SideSize+j] = a1_gt*x1[i*SideSize+j] + a2_gt*x2[i*SideSize+j];
	    }
	 
    std::cout << "Lnear_bivariate_example" <<std::endl;

  
    Problem problem1;
	
	double a1 = 0;
	double a2 = 0;
	
	
    for (int i = 0; i < nb_obs; i++) 
    {
	problem1.AddResidualBlock(new AutoDiffCostFunction<LinBivariateResidual, 1, 1, 1>
		        (new LinBivariateResidual(x1[i],x2[i], y[i])),nullptr,  &a1, &a2);
	}
	  
    Solver::Options options;
    Solver::Summary summary1;
    
	options.max_num_iterations = 125;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = true;
    
	Solve(options, &problem1, &summary1);

    std::cout << "Estimated factors of the linear bivariate. a1 : " << a1 << " - a2 : " << a2 << "\n";		
}


void linear_transformation_estimate()
{
	//This example simulates the situation of a function R^2 to R^2
	// y1 = a11*x1+a12*x2
	// y2 = a21*x1+a22*x2
	
	//Data generation 
    int SideSize = 5;
    int nb_obs =  SideSize * SideSize;
	std::vector<double> x1(nb_obs), x2(nb_obs);
	std::vector<double> y1(nb_obs), y2(nb_obs);
	
	double a11_gt= 1.,  a12_gt = 2.4;
	double a21_gt= 0.9, a22_gt = -4.4;
	for(int i = 0;i<SideSize;i++)
		for(int j = 0;j<SideSize;j++)
     	{
			std::cout << "i*SideSize+j : " << i*SideSize+j << std::endl;
			x1[i*SideSize+j] = (double)i;
			x2[i*SideSize+j] = (double)j;
			y1[i*SideSize+j] = a11_gt*x1[i*SideSize+j] + a12_gt*x2[i*SideSize+j];
			y2[i*SideSize+j] = a21_gt*x1[i*SideSize+j] + a22_gt*x2[i*SideSize+j];
	    }
	 
    std::cout << "Lnear_transformation_estimate" <<std::endl;
  
    Problem problem1;
	
	double a11 = 0, a12 = 0, a21 = 0, a22 = 0;
		
    for (int i = 0; i < nb_obs; i++) 
    {
	    problem1.AddResidualBlock(new AutoDiffCostFunction<LinBivariateResidual, 1, 1, 1>
		        (new LinBivariateResidual(x1[i],x2[i], y1[i])),nullptr,  &a11, &a12);
    	problem1.AddResidualBlock(new AutoDiffCostFunction<LinBivariateResidual, 1, 1, 1>
		        (new LinBivariateResidual(x1[i],x2[i], y2[i])),nullptr,  &a21, &a22);
	}
	  
    Solver::Options options;
    Solver::Summary summary1;
    
	options.max_num_iterations = 125;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = true;
    
	Solve(options, &problem1, &summary1);


    std::cout << "Estimated matrix : " << std::endl;
	std::cout <<	"a11 : " << a11 << " - a12 : " << a12 << "\n";		
	std::cout <<	"a21 : " << a21 << " - a22 : " << a22 << "\n";		
}


void non_linear_transformation_estimate()
{
   std::cout << "non linear transformation" << std::endl;
   //This example simulates the situation of a function R^2 to R^2
   // y1 = (a11*x1+a12*x2+a13)/(a31*x1+a32+x2+a33)
   // y2 = (a21*x1+a22*x2+a23 )/(a31*x1+a32+x2+a33)

   //TODO: ecrire cette partie
   //Data generation
   int SideSize = 5;
   int nb_obs =  SideSize * SideSize * SideSize;
   std::vector<double> x1(nb_obs), x2(nb_obs), x3(nb_obs);
   std::vector<double> y1(nb_obs), y2(nb_obs);

   double r11_gt= 1.,  r12_gt=0, r13_gt = 0;
   double r21_gt= 0.,  r22_gt=1.5, r23_gt = 0.0;
   double r31_gt= 0.,  r32_gt=0, r33_gt = 1.2;
   int indx;
   for(int i = 0;i<SideSize;i++)
       for(int j = 0;j<SideSize;j++)
       {
           for(int k = 0;k<SideSize;k++)
           {
             indx = SideSize*SideSize*i+SideSize*j+k;
             x1[indx] = (double)i;
             x2[indx] = (double)j;
             x3[indx] = (double)k;
             y1[indx] = (r11_gt*x1[indx] + r12_gt*x2[indx] + r13_gt*x3[indx])/(1.0+r33_gt*x3[indx]); // /(r31_gt*x1[indx]+r32_gt*x2[indx]+r33_gt*x3[indx]);
             y2[indx] = (r21_gt*x1[indx] + r22_gt*x2[indx] + r23_gt*x3[indx])/(1.0+r33_gt*x3[indx]); // /(r31_gt*x1[indx]+r32_gt*x2[indx]+r33_gt*x3[indx]);
           }
       }

   std::cout << "non Linear_transformation_estimate" <<std::endl;
   Problem problem;

   double r11 = 1, r12 = 0, r13 = 0;
   double r21 = 0, r22 = 1, r23 = 0;
   double r31 = 0, r32 = 0, r33 = 1;

   for (int i = 0; i < nb_obs; i++)
   {
       problem.AddResidualBlock(new AutoDiffCostFunction<HomographicResidual, 1, 1, 1, 1, 1, 1, 1>
               (new HomographicResidual(x1[i],x2[i], x3[i],  y1[i])),nullptr,  &r11, &r12, & r13, &r31, &r32, &r33);
       problem.AddResidualBlock(new AutoDiffCostFunction<HomographicResidual, 1, 1, 1, 1, 1, 1, 1>
               (new HomographicResidual(x1[i],x2[i],x3[i], y2[i])),nullptr,  &r21, &r22, & r23, &r31, &r32, &r33);
   }

   Solver::Options options;
   Solver::Summary summary;

   options.max_num_iterations = 125;
   options.linear_solver_type = ceres::DENSE_QR;
   options.minimizer_progress_to_stdout = true;

   Solve(options, &problem, &summary);
   double eps = 0.00001;
   if(fabs(r11) < eps)  r11= 0.;
   if(fabs(r12) < eps)  r12= 0.;
   if(fabs(r13) < eps)  r13= 0.;
   if(fabs(r21) < eps)  r21= 0.;
   if(fabs(r22) < eps)  r22= 0.;
   if(fabs(r23) < eps)  r23= 0.;
   if(fabs(r31) < eps)  r31= 0.;
   if(fabs(r32) < eps)  r32= 0.;
   if(fabs(r33) < eps)  r33= 0.;

   std::cout << "Estimated non linear model: " << std::endl;
   std::cout <<	"r11 : " << r11 << ",  r12 : " << r12 << ",  r13 : " << r13 << "\n";
   std::cout <<	"r21 : " << r21 << ",  r22 : " << r22 << ",  r23 : " << r23 << "\n";
   std::cout <<	"r31 : " << r31 << ",  r32 : " << r32 << ",  r33 : " << r33 << "\n";
}

// 1D homographic function
// Notation are those used in the document Stereo3D
// Parameters: tx ty alpha and f
struct HomographyResidual1D
{
  HomographyResidual1D(double X, double Z, double x) : X_(X),  Z_(Z),x_(x) {}

  template <typename T>
  bool operator()(const T* tx, const T* tz, const T* f, const T* alpha, T* residual
          ) const {
    residual[0] = x_ - f[0]*(cos(alpha[0])*X_ - sin(alpha[0])*Z_ + tx[0])/(sin(alpha[0])*X_ + cos(alpha[0])*Z_ + tz[0]);
    return true;
  }

  private:
    const double X_;
    const double Z_;
    const double x_;
};


void homography1D_estimate()
{
    std::cout << std::endl;
    std::cout << "*******************************" << std::endl;
    std::cout << "Homography - 1D Version" << std::endl;
    std::cout << "*******************************" << std::endl;
   //This example simulates the situation of a function R^2 to R^2
   // y1 = (a11*x1+a12*x2+a13)/(a31*x1+a32+x2+a33)
   // y2 = (a21*x1+a22*x2+a23 )/(a31*x1+a32+x2+a33)

    //Data generation
    int SideSize = 5;
    int nb_obs =  SideSize * SideSize;
    std::vector<double> X(nb_obs), Z(nb_obs); // Coordinate of the real points
    std::vector<double> x(nb_obs); // Image coordinates
    int indx;

    //ground truth parameters
    double alpha_gt = -0.03;
    double f_gt = 90.;
    double tx_gt = 0.;
    double tz_gt = 380.;

    double XOffset=-2., XStep=1.;
    double ZOffset=360, ZStep=2;

    //Initial values
    double alpha_init = 0.02;
    double f_init = 70.;
    double tx_init = 0.1;
    double tz_init = 310.;

    //Initial values
    double alpha = alpha_init;
    double f = f_init;
    double tx = tx_init;
    double tz = tz_init;


    for(int i = 0;i<SideSize;i++)
       for(int j = 0;j<SideSize;j++)
       {
           indx = SideSize*i + j;

           X[indx] = XOffset + (double)i*XStep;
           Z[indx] = ZOffset + (double)j*ZStep;

           x[indx] = f_gt * (cos(alpha_gt) *X[indx] - sin(alpha_gt)*Z[indx] + tx_gt)/(sin(alpha_gt) *X[indx] + cos(alpha_gt)*Z[indx] + tz_gt);
       }

   Problem problem;

   for (int indx = 0; indx < nb_obs; indx++)
   {
       problem.AddResidualBlock(new AutoDiffCostFunction<HomographyResidual1D, 1, 1, 1, 1, 1>
               (new HomographyResidual1D(X[indx],Z[indx],x[indx])),nullptr, &tx,&tz,&f,&alpha);
   }

   Solver::Options options;
   Solver::Summary summary;

   options.max_num_iterations = 125;
   options.linear_solver_type = ceres::DENSE_QR;
   options.minimizer_progress_to_stdout = true;

   Solve(options, &problem, &summary);
   std::cout << "Ground truth non linear model: " << std::endl;
   std::cout <<	"tx    : "  << tx_gt << std::endl;
   std::cout <<	"tz    : "  << tz_gt<< std::endl;
   std::cout <<	"alpha : "  << alpha_gt << std::endl;
   std::cout <<	"f     : "  << f_gt << std::endl;

   std::cout << "Parameter initialisation : " << std::endl;
   std::cout <<	"tx    : "  << tx_init << std::endl;
   std::cout <<	"tz    : "  << tz_init<< std::endl;
   std::cout <<	"alpha : "  << alpha_init << std::endl;
   std::cout <<	"f     : "  << f_init << std::endl;

   std::cout << "Estimated non linear model: " << std::endl;
   std::cout <<	"tx    : "  << tx << std::endl;
   std::cout <<	"tz    : "  << tz<< std::endl;
   std::cout <<	"alpha : "  << alpha << std::endl;
   std::cout <<	"f     : "  << f << std::endl;

   double x_gt,x_est;
   for(int indx=0;indx<(int)x.size();indx++)
   {
       x_gt = f_gt * (cos(alpha_gt) *X[indx] - sin(alpha_gt)*Z[indx] + tx_gt)/(sin(alpha_gt) *X[indx] + cos(alpha_gt)*Z[indx] + tz_gt);
       x_est = f * (cos(alpha) *X[indx] - sin(alpha)*Z[indx] + tx)/(sin(alpha) *X[indx] + cos(alpha)*Z[indx] + tz);

       std::cout << "(X,Z) : "<< X[indx] << ", "<< Z[indx]  << " -> x_gt : " << x_gt << ", x_est : " << x_est << std::endl;
   }

   //Comparison of the values obtained with the ground truth model and with the
}


// *************************************************
// Reduced 2D homographic function
// *************************************************
// Notation are those used in the document Stereo3D
// Parameters: tx ty alpha and f
struct Homography2DReducedResidual_x
{
  Homography2DReducedResidual_x(double X, double Y, double Z, double x) : X_(X),  Y_(Y), Z_(Z),x_(x) {}
  template <typename T>
  bool operator()(const T* tx, const T* ty, const T* tz, const T* f, const T* alpha, T* residual
          ) const {
    residual[0] = x_ - f[0]*(cos(alpha[0])*X_ - sin(alpha[0])*Z_ + tx[0])/(sin(alpha[0])*X_ + cos(alpha[0])*Z_ + tz[0]);
    return true;
  }

  private:
    const double X_;
    const double Y_;
    const double Z_;
    const double x_;
};

struct Homography2DReducedResidual_y
{
  Homography2DReducedResidual_y(double X, double Y, double Z, double y) : X_(X),  Y_(Y), Z_(Z),y_(y) {}
  template <typename T>
  bool operator()(const T* tx, const T* ty, const T* tz, const T* f, const T* alpha, T* residual
          ) const {
    residual[0] = y_ - f[0]*(Y_ + ty[0])/(sin(alpha[0])*X_ + cos(alpha[0])*Z_ + tz[0]);
    return true;
  }

  private:
    const double X_;
    const double Y_;
    const double Z_;
    const double y_;
};




void homography2D_reduced_estimate()
{
    std::cout << std::endl;
    std::cout << "*******************************" << std::endl;
    std::cout << "Homography - 2D Version" << std::endl;
    std::cout << "*******************************" << std::endl;

    //Data generation
    int SizeX = 5;
    int SizeY = 5;
    int SizeZ = 5;
    int nb_obs =  SizeX * SizeY * SizeZ;
    std::vector<double> X(nb_obs), Y(nb_obs), Z(nb_obs); // Coordinate of the real points
    std::vector<double> x(nb_obs),y(nb_obs); // Image coordinates
    int indx;

    //ground truth parameters
    double tx_gt = 0.32;
    double ty_gt = 0.05;
    double tz_gt = 380.;
    double alpha_gt = -0.03;
    double beta_gt = 0.0;
    double gamma_gt = 0.0;
    double f_gt = 90.;

    double XOffset=-2., XStep=1.;
    double YOffset=-2., YStep=1.;
    double ZOffset=360, ZStep=2;

    //Initial values
    double tx_init = -0.3;
    double ty_init = 0.5;
    double tz_init = 310.;
    double alpha_init = 0.02;
    double beta_init = 0.0;
    double gamma_init = 0.0;
    double f_init = 70.;

    //Initial values
    double tx = tx_init;
    double ty = ty_init;
    double tz = tz_init;
    double alpha = alpha_init;
    double beta = beta_init;
    double gamma  = gamma_init;
    double f = f_init;


    for(int i = 0;i<SizeX;i++)
       for(int j = 0;j<SizeY;j++)
       {
           for(int k = 0;k<SizeZ;k++)
           {
             indx = SizeX*SizeY*i + j*SizeZ +k;

             X[indx] = XOffset + (double)i*XStep;
             Y[indx] = YOffset + (double)j*YStep;
             Z[indx] = ZOffset + (double)k*ZStep;

             x[indx] = f_gt * (cos(alpha_gt) *X[indx] - sin(alpha_gt)*Z[indx] + tx_gt)/(sin(alpha_gt) *X[indx] + cos(alpha_gt)*Z[indx] + tz_gt);
             y[indx] = f_gt * (Y[indx] + ty_gt)/(sin(alpha_gt) *X[indx] + cos(alpha_gt)*Z[indx] + tz_gt);
           }
       }

   Problem problem;

   //Reduced model
   for (int indx = 0; indx < nb_obs; indx++)
   {
       problem.AddResidualBlock(new AutoDiffCostFunction<Homography2DReducedResidual_x, 1, 1, 1, 1, 1, 1>
               (new Homography2DReducedResidual_x(X[indx],Y[indx],Z[indx],x[indx])),nullptr, &tx,&ty,&tz,&f,&alpha);
       problem.AddResidualBlock(new AutoDiffCostFunction<Homography2DReducedResidual_y, 1, 1, 1, 1, 1, 1>
               (new Homography2DReducedResidual_y(X[indx],Y[indx],Z[indx],y[indx])),nullptr, &tx,&ty,&tz,&f,&alpha);
   }

   Solver::Options options;
   Solver::Summary summary;

   options.max_num_iterations = 125;
   options.linear_solver_type = ceres::DENSE_QR;
   options.minimizer_progress_to_stdout = true;

   Solve(options, &problem, &summary);
   std::cout << std::endl;
   std::cout << "*******************************" << std::endl;
   std::cout << "*** HOMOGRAPHY 2D MODEL ***" << std::endl;
   std::cout << "*******************************" << std::endl;

   std::cout << "Ground truth 2d homography model: " << std::endl;
   std::cout <<	"tx    : "  << tx_gt << std::endl;
   std::cout <<	"ty    : "  << ty_gt << std::endl;
   std::cout <<	"tz    : "  << tz_gt<< std::endl;
   std::cout <<	"alpha : "  << alpha_gt << std::endl;
   std::cout <<	"f     : "  << f_gt << std::endl;

   std::cout << "Parameter initialisation : " << std::endl;
   std::cout <<	"tx    : "  << tx_init << std::endl;
   std::cout <<	"ty    : "  << ty_init << std::endl;
   std::cout <<	"tz    : "  << tz_init<< std::endl;
   std::cout <<	"alpha : "  << alpha_init << std::endl;
   std::cout <<	"f     : "  << f_init << std::endl;

   std::cout << "Estimated non linear model: " << std::endl;
   std::cout <<	"tx    : "  << tx << std::endl;
   std::cout <<	"ty    : "  << ty << std::endl;
   std::cout <<	"tz    : "  << tz<< std::endl;
   std::cout <<	"alpha : "  << alpha << std::endl;
   std::cout <<	"f     : "  << f << std::endl;

   double x_gt,x_est,y_gt,y_est;
   for(int indx=0;indx<(int)x.size();indx++)
   {
       x_gt = f_gt * (cos(alpha_gt) *X[indx] - sin(alpha_gt)*Z[indx] + tx_gt)/(sin(alpha_gt) *X[indx] + cos(alpha_gt)*Z[indx] + tz_gt);
       x_est = f * (cos(alpha) *X[indx] - sin(alpha)*Z[indx] + tx)/(sin(alpha) *X[indx] + cos(alpha)*Z[indx] + tz);
       y_gt = f_gt * (Y[indx] +ty_gt)/(sin(alpha_gt) *X[indx] + cos(alpha_gt)*Z[indx] + tz_gt);
       y_est = f * (Y[indx]+ ty)/(sin(alpha) *X[indx] + cos(alpha)*Z[indx] + tz);

       std::cout << "(X,Y,Z) : "<< X[indx] << ", "<< Y[indx] << ", "<< Z[indx]  << " -> x_gt : " << x_gt << ", x_est : " << x_est <<" ; y_gt : " << y_gt << ", y_est : " << y_est << std::endl;
   }

   //Comparison of the values obtained with the ground truth model and with the
}



// *************************************************
// Complete 2D homographic function
// *************************************************
// Notation are those used in the document Stereo3D
// Parameters: tx ty alpha and f
struct Homography2DResidual_x
{
  Homography2DResidual_x(double X, double Y, double Z, double x) : X_(X),  Y_(Y), Z_(Z),x_(x) {}
  template <typename T>
  bool operator()(const T* tx, const T* ty, const T* tz, const T* f, const T* alpha, const T* beta, const T* gamma, T* residual
          ) const {
      residual[0] = (cos(alpha[0])*cos(gamma[0]) + sin(gamma[0])*sin(alpha[0])*sin(beta[0]) ) *X_;
      residual[0] += -cos(beta[0])*sin(gamma[0])*Y_;
      residual[0] +=  (sin(gamma[0])*cos(alpha[0])*sin(beta[0]) -sin(alpha[0])*cos(gamma[0])) *Z_;
      residual[0] += tx[0];
      residual[0] /=  ( cos(beta[0]) * sin(alpha[0]) * X_ -sin(beta[0])*Y_ + cos(alpha[0]) *cos(beta[0]) *Z_ + tz[0]);
      residual[0] *=f[0];
      residual[0] = x_ - residual[0];
      return true;
  }

  private:
    const double X_;
    const double Y_;
    const double Z_;
    const double x_;
};

struct Homography2DResidual_y
{
  Homography2DResidual_y(double X, double Y, double Z, double y) : X_(X),  Y_(Y), Z_(Z),y_(y) {}
  template <typename T>
  bool operator()(const T* tx, const T* ty, const T* tz, const T* f, const T* alpha, const T* beta, const T* gamma,T* residual
          ) const {

      residual[0] = (sin(gamma[0])*cos(alpha[0]) - cos(gamma[0]) *sin(alpha[0])* sin(beta[0])) *X_;
      residual[0] += cos(gamma[0])*cos(beta[0])*Y_;
      residual[0] +=  (-sin(alpha[0])*sin(gamma[0]) -cos(gamma[0])*sin(beta[0])*cos(alpha[0]) ) *Z_;
      residual[0] += ty[0];
      residual[0] /=  (cos(beta[0]) * sin(alpha[0]) * X_ -sin(beta[0])*Y_ + cos(alpha[0]) *cos(beta[0]) *Z_ + tz[0]);
      residual[0] *= f[0];
      residual[0] = y_ - residual[0];
      return true;
  }

  private:
    const double X_;
    const double Y_;
    const double Z_;
    const double y_;
};


void homography2D_estimate()
{
    std::cout << std::endl;
    std::cout << "*******************************" << std::endl;
    std::cout << "Homography - 2D Version" << std::endl;
    std::cout << "*******************************" << std::endl;

    //Data generation
    int SizeX = 21;
    int SizeY = 21;
    int SizeZ = 11;
    int nb_obs =  SizeX * SizeY * SizeZ;
    std::vector<double> X(nb_obs), Y(nb_obs), Z(nb_obs); // Coordinate of the real points
    std::vector<double> x(nb_obs),y(nb_obs); // Image coordinates
    int indx;

    //ground truth parameters
    double tx_gt = 20.32;
    double ty_gt = 0.05;
    double tz_gt = 380.;
    double alpha_gt = -0.03;
    double beta_gt = -0.2;
    double gamma_gt = 0.1;
    double f_gt = 90.;

    //Initial values
    double tx_init = -0.3;
    double ty_init = 1.5;
    double tz_init = 310.;
    double alpha_init = 0.02;
    double beta_init = 0.65;
    double gamma_init = -0.2;
    double f_init = 60.;

    //Initial values
    double tx = tx_init;
    double ty = ty_init;
    double tz = tz_init;
    double alpha = alpha_init;
    double beta = beta_init;
    double gamma  = gamma_init;
    double f = f_init;

    //Data generation (simulation)
    double XOffset=-4., XStep=0.5;
    double YOffset=-4., YStep=0.5;
    double ZOffset=370, ZStep=1;

    //Perturbation level
    double sigma = 0.02;

    for(int i = 0;i<SizeX;i++)
       for(int j = 0;j<SizeY;j++)
       {
           for(int k = 0;k<SizeZ;k++)
           {
             indx = SizeY*SizeZ*i + j*SizeZ +k;

             X[indx] = XOffset + (double)i*XStep;
             Y[indx] = YOffset + (double)j*YStep;
             Z[indx] = ZOffset + (double)k*ZStep;

             //first image coordinate computation
             x[indx]= (cos(alpha_gt)*cos(gamma_gt) + sin(gamma_gt)*sin(alpha_gt)*sin(beta_gt) ) *X[indx];
             x[indx] += -cos(beta_gt)*sin(gamma_gt)*Y[indx];
             x[indx] += (sin(gamma_gt)*cos(alpha_gt)*sin(beta_gt) -sin(alpha_gt)*cos(gamma_gt)) *Z[indx];
             x[indx] += tx_gt;
             x[indx] /=  (cos(beta_gt)*sin(alpha_gt)*X[indx] -sin(beta_gt)*Y[indx] + cos(alpha_gt)*cos(beta_gt) *Z[indx] + tz_gt);
             x[indx] = f_gt *x[indx];

             //second image coordinate computation
             //y[indx] = f_gt * (Y[indx] + ty_gt)/(sin(alpha_gt) *X[indx] + cos(alpha_gt)*Z[indx] + tz_gt);
             y[indx] = (sin(gamma_gt)*cos(alpha_gt) - cos(gamma_gt) *sin(alpha_gt)* sin(beta_gt)) *X[indx];
             y[indx] += cos(gamma_gt)*cos(beta_gt)*Y[indx];
             y[indx] +=  (-sin(alpha_gt)*sin(gamma_gt) -cos(gamma_gt)*sin(beta_gt)*cos(alpha_gt) ) *Z[indx];
             y[indx] += ty_gt;
             y[indx] /=  (cos(beta_gt)*sin(alpha_gt)*X[indx] -sin(beta_gt)*Y[indx] + cos(alpha_gt)*cos(beta_gt) *Z[indx] + tz_gt);
             y[indx] = f_gt *y[indx];

             //Add some noise
             std::cout << "x[indx] : " << x[indx] << std::endl;
             std::cout << "y[indx] : " << y[indx] << std::endl;
             x[indx] += sigma * randn(2*indx);
             y[indx] += sigma * randn(2*indx+1);
           }
       }

   Problem problem;

   //Reduced model
   for (int indx = 0; indx < nb_obs; indx++)
   {
       problem.AddResidualBlock(new AutoDiffCostFunction<Homography2DResidual_x, 1,1,1, 1, 1, 1, 1, 1>
               (new Homography2DResidual_x(X[indx],Y[indx],Z[indx],x[indx])),nullptr, &tx,&ty,&tz,&f,&alpha,&beta,&gamma);
       problem.AddResidualBlock(new AutoDiffCostFunction<Homography2DResidual_y, 1,1,1, 1, 1, 1, 1, 1>
               (new Homography2DResidual_y(X[indx],Y[indx],Z[indx],y[indx])),nullptr, &tx,&ty,&tz,&f,&alpha,&beta,&gamma);
   }

   Solver::Options options;
   Solver::Summary summary;

   options.max_num_iterations = 125;
   options.linear_solver_type = ceres::DENSE_QR;
   options.minimizer_progress_to_stdout = true;

   Solve(options, &problem, &summary);
   std::cout << std::endl;
   std::cout << "*******************************" << std::endl;
   std::cout << "*** HOMOGRAPHY 2D MODEL ***" << std::endl;
   std::cout << "*******************************" << std::endl;

   std::cout << "Ground truth 2d homography model: " << std::endl;
   std::cout <<	"tx    : "  << tx_gt << std::endl;
   std::cout <<	"ty    : "  << ty_gt << std::endl;
   std::cout <<	"tz    : "  << tz_gt<< std::endl;
   std::cout <<	"alpha : "  << alpha_gt << std::endl;
   std::cout <<	"beta  : "  << beta_gt << std::endl;
   std::cout <<	"gamma : "  << gamma_gt << std::endl;
   std::cout <<	"f     : "  << f_gt << std::endl;

   std::cout << "Parameter initialisation : " << std::endl;
   std::cout <<	"tx    : "  << tx_init << std::endl;
   std::cout <<	"ty    : "  << ty_init << std::endl;
   std::cout <<	"tz    : "  << tz_init<< std::endl;
   std::cout <<	"alpha : "  << alpha_init << std::endl;
   std::cout <<	"beta  : "  << beta_init << std::endl;
   std::cout <<	"gamma : "  << gamma_init << std::endl;
   std::cout <<	"f     : "  << f_init << std::endl;

   std::cout << "Estimated non linear model: " << std::endl;
   std::cout <<	"tx    : "  << tx << std::endl;
   std::cout <<	"ty    : "  << ty << std::endl;
   std::cout <<	"tz    : "  << tz<< std::endl;
   std::cout <<	"alpha : "  << alpha << std::endl;
   std::cout <<	"beta  : "  << beta << std::endl;
   std::cout <<	"gamma : "  << gamma << std::endl;
   std::cout <<	"f     : "  << f << std::endl;

   double x_gt,x_est,y_gt,y_est;
   //TODO: write the complete function here
   /*for(int indx=0;indx<(int)x.size();indx++)
   {
       //x_gt = f_gt * (cos(alpha_gt) *X[indx] - sin(alpha_gt)*Z[indx] + tx_gt)/(sin(alpha_gt) *X[indx] + cos(alpha_gt)*Z[indx] + tz_gt);
       //x_est = f * (cos(alpha) *X[indx] - sin(alpha)*Z[indx] + tx)/(sin(alpha) *X[indx] + cos(alpha)*Z[indx] + tz);
       //y_gt = f_gt * (Y[indx] +ty_gt)/(sin(alpha_gt) *X[indx] + cos(alpha_gt)*Z[indx] + tz_gt);
       //y_est = f * (Y[indx]+ ty)/(sin(alpha) *X[indx] + cos(alpha)*Z[indx] + tz);

       std::cout << "(X,Y,Z) : "<< X[indx] << ", "<< Y[indx] << ", "<< Z[indx]  << " -> x_gt : " << x_gt << ", x_est : " << x_est <<" ; y_gt : " << y_gt << ", y_est : " << y_est << std::endl;
   }*/

   //Comparison of the values obtained with the ground truth model and with the
}




double my_fabs(double x,double y)
{
    return fabs(x-y);
}



// ****************************************************************
// TODO (JCR-LD): Here check the angle convention
// Warning: any change here affects also the forward model
//          (in simulation situation) and the convention for the
//          parametrization of the camera poses (s. Euler angles)
// ****************************************************************
template <typename T>
T ComputeRay(double xL, double yL,
             T f,
             T alpha, T beta, T gamma, //extrinsic left camera parameters
             T& ray_x, T& ray_y, T& ray_z
             )
{
    ray_x =  (cos(alpha)*cos(gamma) +sin(gamma)*sin(alpha)*sin(beta)) * xL;
    ray_x += (sin(gamma)*cos(alpha)-cos(gamma)*sin(alpha)*sin(beta)) * yL;
    ray_x += cos(beta) * sin(alpha)* f;

    ray_y =  -cos(beta)*sin(gamma) * xL;
    ray_y += cos(gamma)*cos(beta)* yL;
    ray_y +=  sin(beta) * f;

    ray_z =  (-sin(alpha)*cos(gamma) + sin(gamma)*cos(alpha)*sin(beta)) * xL;
    ray_z += (-sin(alpha)*sin(gamma) -cos(gamma)*sin(beta)*cos(alpha)) * yL;
    ray_z += cos(alpha) * cos(beta)* f;
}


//adapted from IntersectLines_PointDirParametrisation
template <typename T>
bool IntersectLines
(T txL, T tyL,T tzL, T RayLx,T RayLy,T RayLz,
 T txR, T tyR,T tzR, T RayRx,T RayRy,T RayRz,
 T& Ox, T &Oy, T& Oz)
{
    T p13x = txL -txR;
    T p13y = tyL -tyR;
    T p13z = tzL -tzR;

    T p21x = RayLx;
    T p21y = RayLy;
    T p21z = RayLz;

    T p43x = RayRx;
    T p43y = RayRy;
    T p43z = RayRz;

    T d1343,d4321, d1321,d4343, d2121;

    d1343 = p13x * p43x + p13y * p43y + p13z * p43z;
    d4321 = p43x * p21x + p43y * p21y + p43z * p21z;
    d1321 = p13x * p21x + p13y * p21y + p13z * p21z;
    d4343 = p43x * p43x + p43y * p43y + p43z * p43z;
    d2121 = p21x * p21x + p21y * p21y + p21z * p21z;

    T denom = d2121 * d4343 - d4321 * d4321;

    //TODO: this is not clean. to be cleaned !!
    double  eps = 0.00000001;
    if (abs(denom) < eps) {  std::cout << "Intersect Lines small determinant !" << std::endl; return(false);}

    T numer = d1343 * d4321 - d1321 * d4343;

    T mua = numer / denom;
    T mub = (d1343 + d4321 * mua) / d4343;

    Ox = 0.5*(txL+ mua * p21x + txR + mub *p43x);
    Oy = 0.5*(tyL+ mua * p21y + tyR + mub *p43y);
    Oz = 0.5*(tzL+ mua * p21z + tzR + mub *p43z);

    return true;
}


// ************************************************************************************
// Reconstruction function of a 3d point from image projections and camera parameters
// Input:  Bilddaten (L + R)
//         left and right camera parameters
// Output: reconstructed 3D point
// ************************************************************************************
template <typename T>
void triangulate(double xL, double yL, double xR, double yR, //Koordinaten auf dem Bild
               T fL, //intrinsic (focal distance) left camera parameter
               T txL, T tyL, T tzL, T alphaL, T betaL, T gammaL, //extrinsic left camera parameters
                T fR, //intrinsic (focal distance) camera parameter
                T txR, T tyR, T tzR, T alphaR, T betaR, T gammaR, //extrinsic right camera parameters
                T& ORec_x, T& ORec_y, T& ORec_z //reconstructed 3D point
              )
{
    //T x_rec;
    //see void StereoCamera::SensorsToObject(Point2D& P1, Point2D& P2, Point3D& O)
    //    Point3D P1_virtual(-P1.x,-P1.y,m_Cam1.m_f), P2_virtual(-P2.x,-P2.y,m_Cam2.m_f);
    //T PLz = -yL;
    //    Point3D Ray1, Ray2; //virtual points in camera coordinates
    //    Ray1 = m_Cam1.Orientation.Rt * P1_virtual;
    //    Ray2 = m_Cam2.Orientation.Rt * P2_virtual;
    T RayLx,RayLy, RayLz;
    T RayRx,RayRy, RayRz;

    ComputeRay(xL, yL, fL, /*txL, tyL, tzL, */ alphaL, betaL, gammaL,RayLx,RayLy,RayLz); //extrinsic left camera parameters
    ComputeRay(xR, yR, fR, /* txR, tyR, tzR,*/ alphaR, betaR, gammaR,RayRx,RayRy,RayRz); //extrinsic left camera parameters

    //std::cout << "Left Ray : " << RayLx <<", "  << RayLy <<", "  << RayLz   << std::endl;
    //std::cout << "Right Ray : " << RayRx <<", "  << RayRy <<", "  << RayRz   << std::endl;

    //    IntersectLines_PointDirParametrisation(m_Cam1.Center,Ray1,m_Cam2.Center,Ray2,O);
    IntersectLines(txL, tyL, tzL, RayLx, RayLy, RayLz,
                   txR, tyR, tzR, RayRx, RayRy, RayRz,
                   ORec_x, ORec_y, ORec_z);
}



template <typename T>
T triangulate_x(double xL, double yL, double xR, double yR,
               T fL, //intrinsic (focal distance) left camera parameter
               T txL, T tyL, T tzL, T alphaL, T betaL, T gammaL, //extrinsic left camera parameters
                T fR, //intrinsic (focal distance) camera parameter
                T txR, T tyR, T tzR, T alphaR, T betaR, T gammaR //extrinsic right camera parameters
              )
{
    // T x_rec;
    // see void StereoCamera::SensorsToObject(Point2D& P1, Point2D& P2, Point3D& O)
    //    Point3D P1_virtual(-P1.x,-P1.y,m_Cam1.m_f), P2_virtual(-P2.x,-P2.y,m_Cam2.m_f);
    //T PLz = -yL;
    //    Point3D Ray1, Ray2; //virtual points in camera coordinates
    //    Ray1 = m_Cam1.Orientation.Rt * P1_virtual;
    //    Ray2 = m_Cam2.Orientation.Rt * P2_virtual;
    T RayLx,RayLy, RayLz;
    T RayRx,RayRy, RayRz;

    ComputeRay(xL, yL, fL, txL, tyL, tzL, alphaL, betaL, gammaL,RayLx,RayLy,RayLz); //extrinsic left camera parameters
    ComputeRay(xR, yR, fR, txR, tyR, tzR, alphaR, betaR, gammaR,RayRx,RayRy,RayRz); //extrinsic left camera parameters

    T ORec_x,ORec_y,ORec_z;
    //    IntersectLines_PointDirParametrisation(m_Cam1.Center,Ray1,m_Cam2.Center,Ray2,O);
    IntersectLines(txL, tyL, tzL, RayLx, RayLy, RayLz,
                   txR, tyR, tzR, RayRx, RayRy, RayRz,
                   ORec_x, ORec_y, ORec_z);
    return ORec_x; //reconstructed x
}


struct Triangulation_x
{
  Triangulation_x(double xL, double yL, double xR, double yR, double X) : xL_(xL), yL_(yL), xR_(xR), yR_(yR), X_(X) {}
  template <typename T>
  bool operator()(
          const T* fL, const T* txL, const T* tyL, const T* tzL, const T* alphaL, const T* betaL, const T* gammaL,
          const T* fR, const T* txR, const T* tyR, const T* tzR, const T* alphaR, const T* betaR, const T* gammaR,
          T* residual
          ) const {
      T Y_Rec,Z_Rec, X_Rec;
      triangulate(xL_,yL_,xR_,yR_,
                          fL[0], txL[0], tyL[0], tzL[0],  alphaL[0], betaL[0], gammaL[0],
                          fR[0], txR[0], tyR[0], tzR[0],  alphaR[0], betaR[0], gammaR[0],
                           X_Rec, Y_Rec, Z_Rec
                         );

    residual[0] = X_Rec-X_;
      return true;
  }

  private:
     const double xL_;
     const double yL_;
     const double xR_;
     const double yR_;
     const double X_;
};


struct Triangulation_y
{
  Triangulation_y(double xL, double yL, double xR, double yR, double Y) : xL_(xL), yL_(yL), xR_(xR), yR_(yR), Y_(Y) {}
  template <typename T>
  bool operator()(
          const T* fL, const T* txL, const T* tyL, const T* tzL, const T* alphaL, const T* betaL, const T* gammaL, //parameter "left" camera
          const T* fR, const T* txR, const T* tyR, const T* tzR, const T* alphaR, const T* betaR, const T* gammaR, //parameter "right" camera
          T* residual
          ) const {
      T X_Rec, Y_Rec, Z_Rec;
      triangulate(xL_,yL_,xR_,yR_,
                          fL[0], txL[0], tyL[0], tzL[0],  alphaL[0], betaL[0], gammaL[0], //extrinsic left camera parameters
                          fR[0], txR[0], tyR[0], tzR[0],  alphaR[0], betaR[0], gammaR[0], //extrinsic left camera parameters
                           X_Rec, Y_Rec, Z_Rec
                         );

      residual[0] = Y_Rec-Y_;
      return true;
  }

  private:
     const double xL_;
     const double yL_;
     const double xR_;
     const double yR_;
     const double Y_;
};


struct Triangulation_z
{
  Triangulation_z(double xL, double yL, double xR, double yR, double Z) : xL_(xL), yL_(yL), xR_(xR), yR_(yR), Z_(Z) {}
  template <typename T>
  bool operator()(
          const T* fL, const T* txL, const T* tyL, const T* tzL, const T* alphaL, const T* betaL, const T* gammaL, //parameter "left" camera
          const T* fR, const T* txR, const T* tyR, const T* tzR, const T* alphaR, const T* betaR, const T* gammaR, //parameter "right" camera
          T* residual
          ) const {
      T X_Rec, Y_Rec, Z_Rec;
      triangulate(xL_,yL_,xR_,yR_,
                          fL[0], txL[0], tyL[0], tzL[0],  alphaL[0], betaL[0], gammaL[0], //extrinsic left camera parameters
                          fR[0], txR[0], tyR[0], tzR[0],  alphaR[0], betaR[0], gammaR[0], //extrinsic left camera parameters
                           X_Rec, Y_Rec, Z_Rec
                         );

      residual[0] = Z_Rec-Z_;
      return true;
  }

  private:
     const double xL_;
     const double yL_;
     const double xR_;
     const double yR_;
     const double Z_;
};


void  stereo_test()
{
    std::cout << "********************" << std::endl;
    std::cout << "Stereo test - 2 cameras" << std::endl;
    std::cout << "********************" << std::endl;


    //TODO: generate two cameras.

    // Approximation
    //Generation of the triangulated 3D reconstruction

}

void triangulate_test()
{
    std::cout << std::endl;
    std::cout << "*******************************" << std::endl;
    std::cout << "*******  Triangulate  **********" << std::endl;
    std::cout << "*******************************" << std::endl;

    // **************************
    // PARAMETERS
    // **************************

    // GROUND TRUTH PARAMETERS
    //ground truth parameters: left camera
    double txL_gt = -20.32;
    double tyL_gt  = 0.0;
    double tzL_gt = 380.;
    double alphaL_gt = 0.03;
    double betaL_gt = 0.;
    double gammaL_gt = 0.;
    double fL_gt = 90.;

    /*double txL_gt = -10.;//-10.;
    double tyL_gt = 0.0;
    double tzL_gt = 0.;
    double alphaL_gt = 0.;
    //double alphaL_gt = -atan(txL_gt/(400.));
    //double alphaL_gt = M_PI*0.25;
    double betaL_gt = -0.4;
    double gammaL_gt = 0.;
    double fL_gt = sqrt(2.)/2;//100.;*/

    //ground truth parameters: right camera
    double txR_gt = 20.32;
    double tyR_gt  = 0.0;
    double tzR_gt = 380.;
    double alphaR_gt = -0.03;
    double betaR_gt = 0.;
    double gammaR_gt = 0.;
    double fR_gt = 90.;
    /*double txR_gt = 0;
    double tyR_gt = 0.0;
    double tzR_gt = 0.;
    double alphaR_gt = 0.0;
    double betaR_gt = 0.;
    double gammaR_gt = 0.;
    double fR_gt = 100.0;  */ //100.;

    // INITIAL VALUES
    //initial parameters: left camera
    double txL_init = -20.32;
    double tyL_init = 0.0;
    double tzL_init = 380.;
    double alphaL_init = 0.02;
    double betaL_init = 0.;
    double gammaL_init = 0.;
    double fL_init = 100; //gt: 90

    //initial parameters: right camera
    double txR_init = 20.32;
    double tyR_init = 0.0;
    double tzR_init = 380.;
    double alphaR_init = -0.04;
    double betaR_init = 0.;
    double gammaR_init = 0.;
    double fR_init = 80.; //gt:90

    //initial parameters: left camera
    double txL = txL_init;
    double tyL = tyL_init;
    double tzL = tzL_init;
    double alphaL = alphaL_init;
    double betaL = betaL_init;
    double gammaL = gammaL_init;
    double fL = fL_init;

    //initial parameters: right camera
    double txR = txR_init;
    double tyR = tyR_init;
    double tzR = tzR_init;
    double alphaR = alphaR_init;
    double betaR = betaR_init;
    double gammaR = gammaR_init;
    double fR = fR_init;

    //Data generation
    //Data generation (simulation)
    /*double XOffset=-4., XStep=0.5;
    double YOffset=-4., YStep=0.5;
    double ZOffset=390, ZStep=1; */

    int SizeX = 21;
    int SizeY = 21;
    int SizeZ = 11;
    /*int SizeX = 1;
     int SizeY = 1;
     int SizeZ = 1;*/
     double XOffset=1., XStep=1;
     double YOffset=0., YStep=1;
     //double ZOffset=2, ZStep=1;
     double ZOffset=400, ZStep=1;

    int nb_obs =  SizeX * SizeY * SizeZ;
    std::vector<double> X(nb_obs), Y(nb_obs), Z(nb_obs); // Coordinate of the real points
    std::vector<double> xL(nb_obs),yL(nb_obs); // (Simulated) Image coordinates
    std::vector<double> xR(nb_obs),yR(nb_obs); // (Simulated) Image coordinates
    int indx;

    //Perturbation level - here noise can be added to the simulated measurements
    double sigma = 0.;

    for(int i = 0;i<SizeX;i++)
       for(int j = 0;j<SizeY;j++)
       {
           for(int k = 0;k<SizeZ;k++)
           {
             indx = SizeY*SizeZ*i + j*SizeZ +k;

             //grid generation for object points (ground truth)
             X[indx] = XOffset + (double)i*XStep;
             Y[indx] = YOffset + (double)j*YStep;
             Z[indx] = ZOffset + (double)k*ZStep;

             // Coordinates of the vector camera center -> object (in the real world coordinate system)
             // Warning: in literature ambiguously documented: the translation vector must be transformed in the camera coordinates
             double XL_d = X[indx]-txL_gt;
             double YL_d = Y[indx]-tyL_gt;
             double ZL_d = Z[indx]-tzL_gt;
             double XR_d = X[indx]-txR_gt;
             double YR_d = Y[indx]-tyR_gt;
             double ZR_d = Z[indx]-tzR_gt;

             //first left image coordinate computation
             xL[indx]= (cos(alphaL_gt)*cos(gammaL_gt) + sin(gammaL_gt)*sin(alphaL_gt)*sin(betaL_gt) ) * XL_d;
             xL[indx] += -cos(betaL_gt)*sin(gammaL_gt)* YL_d;
             xL[indx] += (sin(gammaL_gt)*cos(alphaL_gt)*sin(betaL_gt) -sin(alphaL_gt)*cos(gammaL_gt)) * ZL_d;
             //xL[indx] -= txL_gt;

             //std::cout << "lbl 1. xL[indx] : " << xL[indx] << std::endl;
             //xL[indx] /=  (cos(betaL_gt)*sin(alphaL_gt)*X_d -sin(betaL_gt)*Y_d + cos(alphaL_gt)*cos(betaL_gt) * Z_d); //- tzL_gt);
             xL[indx] /=  (cos(betaL_gt)*sin(alphaL_gt)*XL_d -sin(betaL_gt)*YL_d + cos(alphaL_gt)*cos(betaL_gt) * ZL_d);// - tzL_gt); //- tzL_gt);
             xL[indx] = fL_gt *xL[indx];
             //std::cout << "lbl 2. xL[indx] : " << xL[indx] << std::endl;

             //second left image coordinate computation
             //y[indx] = f_gt * (Y[indx] + ty_gt)/(sin(alpha_gt) *X[indx] + cos(alpha_gt)*Z[indx] + tz_gt);
             yL[indx] = (sin(gammaL_gt)*cos(alphaL_gt) - cos(gammaL_gt) *sin(alphaL_gt)* sin(betaL_gt)) * XL_d;
             yL[indx] += cos(gammaL_gt)*cos(betaL_gt) * YL_d;
             yL[indx] +=  (-sin(alphaL_gt)*sin(gammaL_gt) -cos(gammaL_gt)*sin(betaL_gt)*cos(alphaL_gt) ) * ZL_d;
             //std::cout << "lbl 3. yL[indx] : " << yL[indx] << std::endl;

             //yL[indx] -= tyL_gt;
             yL[indx] /=  (cos(betaL_gt)*sin(alphaL_gt)*XL_d -sin(betaL_gt)*YL_d + cos(alphaL_gt)*cos(betaL_gt) *ZL_d);// - tzL_gt);
             yL[indx] = fL_gt *yL[indx];
             //std::cout << "lbl 4. yL[indx] : " << yL[indx] << std::endl;

             //first right image coordinate computation
             xR[indx] = (cos(alphaR_gt)*cos(gammaR_gt) + sin(gammaR_gt)*sin(alphaR_gt)*sin(betaR_gt) ) *XR_d;//X[indx];
             xR[indx] += -cos(betaR_gt)*sin(gammaR_gt)*YR_d;//Y[indx];
             xR[indx] += (sin(gammaR_gt)*cos(alphaR_gt)*sin(betaR_gt) -sin(alphaR_gt)*cos(gammaR_gt)) * ZR_d;//Z[indx];
             //std::cout << "lbl 5. xR[indx] : " << xR[indx] << std::endl;

             //xR[indx] -= txR_gt;
             xR[indx] /=  (cos(betaR_gt)*sin(alphaR_gt)*XR_d -sin(betaR_gt)*YR_d + cos(alphaR_gt)*cos(betaR_gt) *ZR_d); //- tzR_gt);
             xR[indx] = fR_gt *xR[indx];
             //std::cout << "lbl 6. xR[indx] : " << xR[indx] << std::endl;

             //second image coordinate computation
             //y[indx] = f_gt * (Y[indx] + ty_gt)/(sin(alpha_gt) *X[indx] + cos(alpha_gt)*Z[indx] + tz_gt);
             yR[indx] = (sin(gammaR_gt)*cos(alphaR_gt) - cos(gammaR_gt) *sin(alphaR_gt)* sin(betaR_gt)) *XR_d;//X[indx];
             yR[indx] += cos(gammaR_gt)*cos(betaR_gt)*YR_d;//Y[indx];
             yR[indx] +=  (-sin(alphaR_gt)*sin(gammaR_gt) -cos(gammaR_gt)*sin(betaR_gt)*cos(alphaR_gt) ) *ZR_d;//Z[indx];
             //yR[indx] -= tyR_gt;
             //std::cout << "lbl 7. yR[indx] : " << yR[indx] << std::endl;

             yR[indx] /=  (cos(betaR_gt)*sin(alphaR_gt)*XR_d -sin(betaR_gt)*YR_d + cos(alphaR_gt)*cos(betaR_gt) *ZR_d);// - tzR_gt);
             yR[indx] = fR_gt * yR[indx];
             //std::cout << "lbl 8. yR[indx] : " << yR[indx] << std::endl;

             //Add some noise
             //std::cout << "x[indx] : " << x[indx] << std::endl;
             //std::cout << "y[indx] : " << y[indx] << std::endl;
             xL[indx] += sigma * randn(4*indx);
             yL[indx] += sigma * randn(4*indx+1);
             xR[indx] += sigma * randn(4*indx+2);
             yR[indx] += sigma * randn(4*indx+3);
           }
       }

  //Test that the triangulation is correct - debug purpose
  double XRec, YRec, ZRec;
  for(int indx = 0; indx < nb_obs; indx++)
  {
     std::cout << "indx :" << indx << std::endl;
     std::cout << "(xL,yL,xR,yR): "  << xL[indx] <<  ", "<< yL[indx] << ", "<< xR[indx] << ", " << yR[indx] << std::endl;
     std::cout << "(X,Y,Z)         : "  << X[indx] <<  ", "<< Y[indx] << ", "<< Z[indx]  << std::endl;

     triangulate(xL[indx],yL[indx],xR[indx],yR[indx],
                         fL_gt, txL_gt, tyL_gt, tzL_gt,  alphaL_gt, betaL_gt, gammaL_gt, //extrinsic left camera parameters
                         fR_gt, txR_gt, tyR_gt, tzR_gt,  alphaR_gt, betaR_gt, gammaR_gt, //extrinsic left camera parameters
                          XRec, YRec, ZRec
                        );

     std::cout << "(XRec,YRec,ZRec): "  << XRec <<  ", "<< YRec << ", "<< ZRec  << std::endl;

     std::cout << std::endl;
  }


  std::cout << "Ground truth stereo 3d reconstruction (triangulation) model: " << std::endl;
  std::cout <<	"txL    : "  << txL_gt << std::endl;
  std::cout <<	"tyL    : "  << tyL_gt << std::endl;
  std::cout <<	"tzL    : "  << tzL_gt<< std::endl;
  std::cout <<	"alphaL : "  << alphaL_gt << std::endl;
  std::cout <<	"betaL  : "  << betaL_gt << std::endl;
  std::cout <<	"gammaL : "  << gammaL_gt << std::endl;
  std::cout <<	"fL    : "  << fL_gt << std::endl;
  std::cout <<	"txR    : "  << txR_gt << std::endl;
  std::cout <<	"tyR    : "  << tyR_gt << std::endl;
  std::cout <<	"tzR    : "  << tzR_gt<< std::endl;
  std::cout <<	"alphaR : "  << alphaR_gt << std::endl;
  std::cout <<	"betaR  : "  << betaR_gt << std::endl;
  std::cout <<	"gammaR : "  << gammaR_gt << std::endl;
  std::cout <<	"fR    : "  << fR_gt << std::endl;
  std::cout << "just before testing" << std::endl;
  std::cout << "just before testing" << std::endl;

   //Now problem solver (CERES)
   Problem problem;
   for(int indx=0;indx < nb_obs;indx++)
   {
       problem.AddResidualBlock(new AutoDiffCostFunction<Triangulation_x, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1>
               (new Triangulation_x(xL[indx],yL[indx],xR[indx],yR[indx],X[indx])),
                nullptr, &fL, &txL, &tyL, &tzL, &alphaL, &betaL, &gammaL,
                &fR, &txR, &tyR, &tzR, &alphaR, &betaR, &gammaR ) ;
       problem.AddResidualBlock(new AutoDiffCostFunction<Triangulation_y, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1>
               (new Triangulation_y(xL[indx],yL[indx],xR[indx],yR[indx],Y[indx])),
                nullptr, &fL, &txL, &tyL, &tzL, &alphaL, &betaL, &gammaL,
                &fR, &txR, &tyR, &tzR, &alphaR, &betaR, &gammaR ) ;
       problem.AddResidualBlock(new AutoDiffCostFunction<Triangulation_z, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1>
               (new Triangulation_z(xL[indx],yL[indx],xR[indx],yR[indx],Z[indx])),
                nullptr, &fL, &txL, &tyL, &tzL, &alphaL, &betaL, &gammaL,
                &fR, &txR, &tyR, &tzR, &alphaR, &betaR, &gammaR ) ;
   }

   Solver::Options options;
   Solver::Summary summary;

   options.max_num_iterations = 125;
   options.linear_solver_type = ceres::DENSE_QR;
   options.minimizer_progress_to_stdout = true;

   Solve(options, &problem, &summary);
   std::cout << std::endl;
   std::cout << "*******************************" << std::endl;
   std::cout << "*** TRIANGULATE TEST***" << std::endl;
   std::cout << "*******************************" << std::endl;


   std::cout << "Parameter initialisation : " << std::endl;
   std::cout <<	"txL    : "  << txL_init << std::endl;
   std::cout <<	"tyL    : "  << tyL_init << std::endl;
   std::cout <<	"tzL    : "  << tzL_init<< std::endl;
   std::cout <<	"alphaL : "  << alphaL_init << std::endl;
   std::cout <<	"betaL  : "  << betaL_init << std::endl;
   std::cout <<	"gammaL : "  << gammaL_init << std::endl;
   std::cout <<	"fL     : "  << fL_init << std::endl;
   std::cout <<	"txR    : "  << txR_init << std::endl;
   std::cout <<	"tyR    : "  << tyR_init << std::endl;
   std::cout <<	"tzR    : "  << tzR_init<< std::endl;
   std::cout <<	"alphaR : "  << alphaR_init << std::endl;
   std::cout <<	"betaR  : "  << betaR_init << std::endl;
   std::cout <<	"gammaR : "  << gammaR_init << std::endl;
   std::cout <<	"fR     : "  << fR_init << std::endl;

   std::cout << "Estimated non linear model: " << std::endl;
   std::cout <<	"txL    : "  << txL << std::endl;
   std::cout <<	"tyL    : "  << tyL << std::endl;
   std::cout <<	"tzL    : "  << tzL << std::endl;
   std::cout <<	"alphaL : "  << alphaL << std::endl;
   std::cout <<	"betaL  : "  << betaL << std::endl;
   std::cout <<	"gammaL : "  << gammaL << std::endl;
   std::cout <<	"fL     : "  << fL << std::endl;
   std::cout <<	"txR    : "  << txR << std::endl;
   std::cout <<	"tyR    : "  << tyR << std::endl;
   std::cout <<	"tzR    : "  << tzR << std::endl;
   std::cout <<	"alphaR : "  << alphaR << std::endl;
   std::cout <<	"betaR  : "  << betaR << std::endl;
   std::cout <<	"gammaR : "  << gammaR << std::endl;
   std::cout <<	"fR     : "  << fR << std::endl;
}

int main(int argc, char** argv) 
{	
  //test_functor();

  //First examples - adapted from the original examples in the ceres examples (come with installation)
  //curve_fitting_examples();
  
  // Example with multivariate optimisation - inspired from the powell example
  // Linear model - estimate of a linear form y = a1*x1+a2*x2
  //linear_bivariate_example();
  
  // Linear model - estimate of a 2x2 matrix transformation
  //linear_transformation_estimate();
  
  // Non linear model - first very elementary example
  //non_linear_transformation_estimate();

  // Homography : elementary 1D variant
  //homography1D_estimate();

  // Homography : 2D variant only alpha as a parameter (beta,gamma not optimized)
  //homography2D_reduced_estimate();

  // Homography : complete 2D variant
  //homography2D_estimate();

  //Triangulate
  triangulate_test();

  return 0;
}






//*************************
// From the original file
//*************************

// Data generated using the following octave code.
//   randn('seed', 23497);
//   m = 0.3;
//   c = 0.1;
//   x=[0:0.075:5];
//   y = exp(m * x + c);
//   noise = randn(size(x)) * 0.2;
//   y_observed = y + noise;
//   data = [x', y_observed'];

/*const int kNumObservations = 67;
// clang-format off
double data[] = {
  0.000000e+00, 1.133898e+00,
  7.500000e-02, 1.334902e+00,
  1.500000e-01, 1.213546e+00,
  2.250000e-01, 1.252016e+00,
  3.000000e-01, 1.392265e+00,
  3.750000e-01, 1.314458e+00,
  4.500000e-01, 1.472541e+00,
  5.250000e-01, 1.536218e+00,
  6.000000e-01, 1.355679e+00,
  6.750000e-01, 1.463566e+00,
  7.500000e-01, 1.490201e+00,
  8.250000e-01, 1.658699e+00,
  9.000000e-01, 1.067574e+00,
  9.750000e-01, 1.464629e+00,
  1.050000e+00, 1.402653e+00,
  1.125000e+00, 1.713141e+00,
  1.200000e+00, 1.527021e+00,
  1.275000e+00, 1.702632e+00,
  1.350000e+00, 1.423899e+00,
  1.425000e+00, 1.543078e+00,
  1.500000e+00, 1.664015e+00,
  1.575000e+00, 1.732484e+00,
  1.650000e+00, 1.543296e+00,
  1.725000e+00, 1.959523e+00,
  1.800000e+00, 1.685132e+00,
  1.875000e+00, 1.951791e+00,
  1.950000e+00, 2.095346e+00,
  2.025000e+00, 2.361460e+00,
  2.100000e+00, 2.169119e+00,
  2.175000e+00, 2.061745e+00,
  2.250000e+00, 2.178641e+00,
  2.325000e+00, 2.104346e+00,
  2.400000e+00, 2.584470e+00,
  2.475000e+00, 1.914158e+00,
  2.550000e+00, 2.368375e+00,
  2.625000e+00, 2.686125e+00,
  2.700000e+00, 2.712395e+00,
  2.775000e+00, 2.499511e+00,
  2.850000e+00, 2.558897e+00,
  2.925000e+00, 2.309154e+00,
  3.000000e+00, 2.869503e+00,
  3.075000e+00, 3.116645e+00,
  3.150000e+00, 3.094907e+00,
  3.225000e+00, 2.471759e+00,
  3.300000e+00, 3.017131e+00,
  3.375000e+00, 3.232381e+00,
  3.450000e+00, 2.944596e+00,
  3.525000e+00, 3.385343e+00,
  3.600000e+00, 3.199826e+00,
  3.675000e+00, 3.423039e+00,
  3.750000e+00, 3.621552e+00,
  3.825000e+00, 3.559255e+00,
  3.900000e+00, 3.530713e+00,
  3.975000e+00, 3.561766e+00,
  4.050000e+00, 3.544574e+00,
  4.125000e+00, 3.867945e+00,
  4.200000e+00, 4.049776e+00,
  4.275000e+00, 3.885601e+00,
  4.350000e+00, 4.110505e+00,
  4.425000e+00, 4.345320e+00,
  4.500000e+00, 4.161241e+00,
  4.575000e+00, 4.363407e+00,
  4.650000e+00, 4.161576e+00,
  4.725000e+00, 4.619728e+00,
  4.800000e+00, 4.737410e+00,
  4.875000e+00, 4.727863e+00,
  4.950000e+00, 4.669206e+00,
};*/
// clang-format on


// *******************************
// LD: July 2022 - add my data
// *******************************
 //int kMyNumObservations = 10;

/*double mydata[] = 
{
	0, 2,
	1, 3,
	2, 4,	
	3, 5,
	4, 2,
	5, 3,
	6, 4,	
	7, 5,
	7, 5,
	7, 5,
};*/	
//

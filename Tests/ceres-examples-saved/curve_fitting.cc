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
//
// Author: sameeragarwal@google.com (Sameer Agarwal)

#include <vector>

#include "ceres/ceres.h"
#include "glog/logging.h"

using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solve;
using ceres::Solver;


struct LinRegResidual
{    
	LinRegResidual(double x, double y) : x_(x), y_(y) {}
    template <typename T>
    bool operator()(const T* const a0, const T* const a1, T* residual) const {
      residual[0] = y_ - a1[0] * x_ - a0[0];
      return true;
    }

   private:
    const double x_;
    const double y_;
};


struct QuadRegResidual
{    
	QuadRegResidual(double x, double y) : x_(x), y_(y) {}
    template <typename T>
    bool operator()(const T* const a0, const T* const a1, const T* const a2, T* residual) const {
      residual[0] = y_ - a2[0] * x_ * x_ - a1[0] * x_ - a0[0];
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



void curve_fitting_examples()
{
    std::cout << "Debut de curve_fitting " << std::endl;	
    //google::InitGoogleLogging(argv[0]);

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




int main(int argc, char** argv) 
{
  //First examples - adapted from the original examples in the ceres examples (come with installation)
  curve_fitting_examples();
  
  // TODO: example with multivariate optimisation - inspired from the powell example
  //Linear model - estimate of a linear form y = a1*x1+a2*x2
  linear_bivariate_example();
  
  //Linear model - estimate of a 2x2 matrix transformation
  linear_transformation_estimate();
  
  return 0;
}



//**********%%%%%%%%%%%%%%%
// From the original file
//**********%%%%%%%%%%%%%%%

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
// LD: juillet 2022 - add my data
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

#ifndef _NUMERIC_H_
#define _NUMERIC_H_

#include <fftw3.h>
#include <chrono>
#include <complex>
//#include "tiffio.h"

#include "../algo_tools/matrix.h"
#include "../algo_tools/utils.h"

using namespace std;

double linear_interpolate(double x, double xl, double xr,
                          double yl, double yr);

//1D piecewise linear approximation
// maybe TODO: pack this in a class
void PiecewiseLinearInterpolation(vector<double>&y, vector<double>& x_part,
                                  vector<double>&y_approx, double a=0, double b=1);

//Solution of the Poisson problem in 2D
void PoissonSolver2D(M3Matrix& b,M3Matrix& U,vector<double>& left,vector<double>& right,
                     vector<double>& top,vector<double>& bottom, double h = 0.1);

// Generate the (separable, 1D) DCT matrix: this implements the DCT-I
void MakeDCT5(M3Matrix& DCTMatrix);

// DCT
void dct5_2d_sep(M3Matrix& ImBlock,M3Matrix& DCTBlock);
void idct5_2d_sep(M3Matrix& DCTBlock,M3Matrix& ImBlock);

void dct5_1d(M3Matrix& InSegment,M3Matrix& DCTSegment);
void idct5_1d(M3Matrix& DCTSegment,M3Matrix& InSegment);


//*************************************************
// Wrapper class for fftw functionalities
//*************************************************
class FFTW_2D
{
public:
  //Default constructor
  FFTW_2D();
  //Default destructor
  virtual ~FFTW_2D();

  FFTW_2D(int Nx, int Ny);

  void SetSizes(int Nx, int Ny);
  void GetData(vector<vector<complex<double>>>& input_vec);
  void GetData(M3Matrix& Im);

  vector<vector<complex<double>>> fft2();
  void fft2(vector<vector<complex<double>>>& fft_mat);
  vector<vector<complex<double>>> ifft2();
  void ifft2(vector<vector<complex<double>>>& ifft_mat);

  void SetForward();
  void SetBackward();

  // members
  int m_Nx = 0;
  int m_Ny = 0;

  fftw_plan m_plan = nullptr;
  fftw_complex* m_data = nullptr;
};


//*********************************************************************
// dft / fft
// adapted from:
// https://www.nayuki.io/page/free-small-fft-in-multiple-languages
// free small fft in multiple languages
//*********************************************************************
vector<complex<double> > fft(vector<complex<double> > &vec);
vector<complex<double>> ifft(vector<complex<double> >& x);
vector<complex<double> > inverseTransform(vector<complex<double> > &vec);
vector<complex<double> > transformRadix2(vector<complex<double> > &vec);
vector<complex<double> > transformBluestein(vector<complex<double> > &vec);
void convolve( const vector<complex<double> > &xvec,
        const vector<complex<double> > &yvec,
        vector<complex<double> > &outvec);
//static size_t reverseBits(size_t x, int n);

vector<vector<complex<double>>> fft2(M3Matrix& Im);
vector<vector<complex<double>>> fft2(vector<vector<complex<double>>>&  input_mat);
vector<vector<complex<double>>> fft2(vector<vector<complex<double>>>&  input_mat, vector<vector<complex<double>>>&  tmp_mat, vector<vector<complex<double>>>&  fft_mat);
vector<vector<complex<double>>> ifft2(vector<vector<complex<double>>>&  input_mat);

vector<vector<complex<double>>> fftshift(vector<vector<complex<double>>>&  input_mat);
void fftshift(vector<vector<complex<double>>>&  input_mat,vector<vector<complex<double>>>&  shift_fft_mat);

//void LoadTiff(string filename,vector<M3Matrix>& RGB);


//************************
//* Windows
//************************
vector<complex<double>> HannWindow(int m);
vector<vector<complex<double>>> gaussianWindow(int m, int n, double sigma);


//**************************
//**** Tests ***************
//**************************
void test_cholesky();
void test_fft();
void test_fft_time();
void test_poisson_solver();
void test_fftw();

void test_fftw_vs_fft();
void test_ifft2();

//void test_tiff();
//void test_HannPadded();
void test_subpixel_translation();

#endif // of _NUMERIC_H_

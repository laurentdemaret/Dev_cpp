#include "numeric.h"

using namespace std::chrono;

static size_t reverseBits(size_t x, int n);


double linear_interpolate(double x, double xl, double xr,
                          double yl, double yr)
{
    double y = 0.;
    //TODO: insert the tests
    if(xl > xr)
    {
        return y;
    }

    y = yr * (x-xl)/(xr-xl) + yl*(xr-x)/(xr-xl);
    return y;
}


// ************************************
//1D piecewise linear approximation
// * a,b: intervall bounds
// * y: (n+1) values, for a spacing of
// 1/(n+1)
// * Default: a=0,b=1,
// *
// ************************************
void PiecewiseLinearInterpolation(vector<double>&y, vector<double>& x_part,
                                  vector<double>&y_approx, double a, double b)
{
   //Check the consistency of the partition
    if(x_part[0] != a)
    {
       std::cout << "Error: ... tbd" << std::endl;
               return;
    }
    if(x_part[x_part.size()-1] != b)
    {
       std::cout << "Error: ... tbd" << std::endl;
               return;
    }

   for(int i = 0;i<x_part.size()-1;i++)
   {
       if(x_part[i+1]<=x_part[i])
       {
         std::cout << "Error: strict monotony of the partition not satisfied" << std::endl;
       }
   }
   //End of check

   //Initialise y_approx
   y_approx.resize(y.size());

   y_approx[0]= a;
   y_approx[y_approx.size()-1]= b;

   //Generate the indices of the vertex for each data point
   int index = 0;
   double step = (b-a)/((double) (y_approx.size()-1));
   std::cout << "step : " << step << std::endl;

   double x_current = a;

   vector<int> indices(y.size());
   vector<double> part_values(x_part.size());
   part_values[0] = y[0];
   part_values[x_part.size()-1] = y[y.size()-1];

   for(int i = 0;i < (int)y.size();i++)
   {
     // here should be a while-loop and not a single if
     // in order to prevent two or more  partition points between two consecutive data points
     x_current = a + i*step;
     if(index <  (int)(x_part.size()-1))
     {
         if(x_current > x_part[index+1])
         {
             index++;
             //Compute the value taken at the coordinate point to interpolate
             double t = (x_current - x_part[index])/ step;
             std::cout << "index : " << index  << ", t : " << t << std::endl;
             part_values[index] = (1-t)*y[i] + t*y[i+1];
         }
     }
     indices[i]= index;
   }

   index = 0;

   //TODO: a first loop to decide the values at the points.
   double left_val,right_val;

   for(int i=1;i<(int)y_approx.size();i++)
   {
       x_current = a + i*step;
       index = indices[i];
       //Now interpolate
       double intervall_length = x_part[index+1] - x_part[index];
       double t = (x_current - x_part[index])/ (intervall_length);
       left_val = part_values[index];
       right_val = part_values[index+1];
       y_approx[i] = t * right_val + (1-t) * left_val;

       std::cout  << "i : " << i << ", x_current : " <<  x_current <<", interval_length"  <<
                  intervall_length << ",  t : " << t
                  << ", left_val : " << left_val << ", right_val : " <<
                  right_val << std::endl;
       std::cout  << "y_approx[" << i << "] : " << y_approx[i] << std::endl;
   }

   return;
}


/*double ReprojectionError(Function_Multivariate f, M3Matrix x_val,  vector<double> p0, vector<double> y_obs)
{
  double error;
  //Generate vector y : function output for x_val
  vector<double> y;

  //compute y
  return error;
}*/


// Quickly written Poisson solver !
// Solves Poisson equation \Delta u = b, for Dirichlet conditions u_{\partial \Omega = 0}
// Input arguments: a matrix b containing the data (as a matrix !!!)
// left, right,top, bottom
void PoissonSolver2D(M3Matrix& b,M3Matrix& U,vector<double>& left,vector<double>& right,vector<double>& top,vector<double>& bottom, double h)
{
    int NbRows = b.GetNbRows();
    int NbCols = b.GetNbCols();

    U.Reshape(NbRows,NbCols);

    M3Matrix uold(NbRows,NbCols), unew(NbRows,NbCols);
    uold = b;
    int NbIt = 500; // Maximal number of iterations
    //Jacobi Iterations
    for(int it = 0 ;it<NbIt;it++)
    {
        for(int i =  0;i<NbRows;i++)
            for(int j =  0;j<NbCols;j++)
            {
                //Jacobi iteration
                double val = 0.;

                if(i>0)
                    val -= uold[i-1][j];
                else
                    val -= top[j];

                if(j>0)
                    val -= uold[i][j-1];
                else
                    val -= left[i];

                if(i<NbRows-1)
                    val -= uold[i+1][j];
                else
                    val -= bottom[j];

                if(j<NbCols-1)
                    val -= uold[i][j+1];
                else
                    val -= right[i];


                val += h*h*b[i][j];

                val /= -4.;
                unew[i][j] = val;
            }
        cout << "Error at iteration " << it << " : "<< (uold-unew).LpNorm(2) << endl;
        uold = unew;
    }

    U = unew;
}


//************************************************************
//* Minimalistic (non-optimized) Fourier routines (DCT, FFT)
//************************************************************
void MakeDCT5(M3Matrix& DCTMatrix)
{
    DCTMatrix.Reshape(5,5);
    double DCT[5][5] =
    {
        {0.5, 1., 1., 1., 0.5},
        {0.5, 1./sqrt(2.), 0., -1./sqrt(2.), -0.5},
        {0.5, 0, -1, 0, 0.5},
        {0.5, -1/sqrt(2), 0, 1/sqrt(2), -0.5},
        { 0.5, -1, 1, -1, 0.5}
    };

    for(int i = 0;i<5;i++)
        for(int j = 0;j<5;j++)
        {
            DCTMatrix[i][j] = DCT[i][j];
        }
}



void dct5_2d_sep(M3Matrix& ImBlock,M3Matrix& DCTBlock)
{
    M3Matrix D,Dt;
    MakeDCT5(D);
    D.Transpose(Dt);

    DCTBlock = D * ImBlock * Dt;
}


void idct5_2d_sep(M3Matrix& DCTBlock,M3Matrix& ImBlock)
{
    M3Matrix D,Dt;
    MakeDCT5(D);
    D.Transpose(Dt);

    ImBlock = D * DCTBlock * Dt;
}


void dct5_1d(M3Matrix& InSegment,M3Matrix& DCTSegment)
{
    M3Matrix D, Dt;
    MakeDCT5(D);
    D.Transpose(Dt);

    DCTSegment = InSegment * Dt;
}


void idct5_1d(M3Matrix& DCTSegment,M3Matrix& OutSegment)
{
    M3Matrix D,Dt;
    MakeDCT5(D);
    D.Transpose(Dt);

    OutSegment = DCTSegment * Dt;
}

//****************************************
//* Windowin
//****************************************
//Gaussian window: see gaussianWindow.m
vector<vector<complex<double>>> gaussianWindow(int m, int n, double sigma)
{
    vector<vector<complex<double>>> g;
    g.resize(m);
    for(int k = 0;k<m;k++)
        g[k].resize(n);

    double centerX = ((n+1.0)/2);
    double centerY = ((m+1.0)/2);
    double a = -2 *PI*PI*sigma*sigma;

    double b, val;

    vector<double> gx, gy;
    for(int j = 0;j<n;j++)
    {
        b = ((j+1)-centerX)/n;
        val = exp(a*b*b);
        gx.push_back(val);
    }

    for(int k = 0;k<m;k++)
    {
        b = ((k+1)-centerY)/m;
        val = exp(a*b*b);
        gy.push_back(val);
    }

    // g = gy'*gx;
    for(int k = 0;k<m;k++)
    {
        for(int j = 0;j<n;j++)
        {
            complex<double> cval = gx[j]*gy[k];
            g[k][j] = cval;
        }
    }

    return g;
}



vector<complex<double>> HannWindow(int m)
{
    vector<complex<double>> hw(m);
    //hw.resize(m);

    for(int k = 0;k < m;k++)
    {
      complex<double> cval = 0.5 - 0.5 * cos(2*PI*k/m);
      hw[k] = cval;
    }

    return hw;
}


//*******************************************************************
// Codes below are from:
// https://www.nayuki.io/page/free-small-fft-in-multiple-languages
// free small fft in multiple languages
//*******************************************************************
vector<complex<double> > fft(vector<complex<double> > &vec)
{
    vector<complex<double> > y;
    size_t n = vec.size();
    if (n == 0)
        return vec;
    else if ((n & (n - 1)) == 0)  // Is power of 2
        y = transformRadix2(vec);
    else  // More complicated algorithm for arbitrary sizes
        y = transformBluestein(vec);
    return y;
}


vector<complex<double> > inverseTransform(vector<complex<double> > &vec)
{
    vector<complex<double> > y;
    y = conj(vec);
    y = fft(y);
    y = conj(y);
    return y;
}


vector<complex<double> > transformRadix2(vector<complex<double> > &vec)
{
    vector<complex<double> >  y = vec;
    // Length variables
    size_t n = vec.size();
    int levels = 0;  // Compute levels = floor(log2(n))
    for (size_t temp = n; temp > 1U; temp >>= 1)
        levels++;
    if (static_cast<size_t>(1U) << levels != n)
        throw std::domain_error("Length is not a power of 2");

    // Trignometric table
    vector<complex<double> > expTable(n / 2);
    for (size_t i = 0; i < n / 2; i++)
        expTable[i] = std::exp(complex<double>(0, -2 * M_PI * i / n));

    // Bit-reversed addressing permutation
    for (size_t i = 0; i < n; i++) {
        size_t j = reverseBits(i, levels);
        if (j > i)
            std::swap(y[i], y[j]);
    }

    // Cooley-Tukey decimation-in-time radix-2 FFT
    for (size_t size = 2; size <= n; size *= 2) {
        size_t halfsize = size / 2;
        size_t tablestep = n / size;
        for (size_t i = 0; i < n; i += size) {
            for (size_t j = i, k = 0; j < i + halfsize; j++, k += tablestep) {
                complex<double> temp = y[j + halfsize] * expTable[k];
                y[j + halfsize] = y[j] - temp;
                y[j] += temp;
            }
        }
        if (size == n)  // Prevent overflow in 'size *= 2'
            break;
    }
    return y;
}


vector<complex<double> > transformBluestein(vector<complex<double> > &vec)
{
    vector<complex<double> >y = vec;

    // Find a power-of-2 convolution length m such that m >= n * 2 + 1
    size_t n = vec.size();
    size_t m = 1;
    while (m / 2 <= n) {
        if (m > SIZE_MAX / 2)
            throw std::length_error("Vector too large");
        m *= 2;
    }

    // Trignometric table
    vector<complex<double> > expTable(n);
    for (size_t i = 0; i < n; i++) {
        unsigned long long temp = static_cast<unsigned long long>(i) * i;
        temp %= static_cast<unsigned long long>(n) * 2;
        double angle = M_PI * temp / n;
        // Less accurate alternative if long long is unavailable: double angle = M_PI * i * i / n;
        expTable[i] = std::exp(complex<double>(0, -angle));
    }

    // Temporary vectors and preprocessing
    vector<complex<double> > av(m);
    for (size_t i = 0; i < n; i++)
        av[i] = vec[i] * expTable[i];
    vector<complex<double> > bv(m);
    bv[0] = expTable[0];
    for (size_t i = 1; i < n; i++)
        bv[i] = bv[m - i] = std::conj(expTable[i]);

    // Convolution
    vector<complex<double> > cv(m);
    convolve(av, bv, cv);

    // Postprocessing
    for (size_t i = 0; i < n; i++)
        y[i] = cv[i] * expTable[i];
    return y;
}


void convolve(const vector<complex<double> > &xvec,
              const vector<complex<double> > &yvec,
              vector<complex<double> > &outvec)
{
    size_t n = xvec.size();
    if (n != yvec.size() || n != outvec.size())
        throw std::domain_error("Mismatched lengths");
    vector<complex<double> > xv = xvec;
    vector<complex<double> > yv = yvec;
    xv = fft(xv);
    yv = fft(yv);
    for (size_t i = 0; i < n; i++)
        xv[i] *= yv[i];
    xv = inverseTransform(xv);
    for (size_t i = 0; i < n; i++)  // Scaling (because this FFT implementation omits it)
        outvec[i] = xv[i] / static_cast<double>(n);
}


static size_t reverseBits(size_t x, int n)
{
    size_t result = 0;
    for (int i = 0; i < n; i++, x >>= 1)
        result = (result << 1) | (x & 1U);
    return result;
}


// inverse fft (in-place)
vector<complex<double>> ifft(vector<complex<double> >& x)
{
    vector<complex<double>> y;
    y = inverseTransform(x);
    y = y / (double(x.size()));
    return y;
}


void fftshift(vector<vector<complex<double>>>&  input_mat,vector<vector<complex<double>>>&  shift_fft_mat)
{
    int height = input_mat.size();
    int width = input_mat[0].size();

    //shifted indices
    int shift_h = (height+1)/2; //vertical shift
    int shift_w = (width+1)/2; //horizontal shift
    int is,js;
    for(int i = 0;i<height;i++)
    {
        is = i + shift_h;
        if(is >= height)
            is -= height;

        for(int j = 0;j<width;j++)
        {
            js = j + shift_w;
            if(js >= width)
                js -= width;
            shift_fft_mat[i][j] = input_mat[is][js];
        }
    }
}


vector<vector<complex<double>>> fftshift(vector<vector<complex<double>>>&  input_mat)
{
    vector<vector<complex<double>>> shift_fft_mat;
    int height = input_mat.size();
    int width = input_mat[0].size();

    //Allocation
    shift_fft_mat.resize(height);
    for(int i = 0;i<(int)input_mat.size();i++)
    {
        shift_fft_mat[i].resize(width);
    }

    //shifted indices
    int shift_h = (height+1)/2; //vertical shift
    int shift_w = (width+1)/2; //horizontal shift
    int is,js;
    for(int i = 0;i<height;i++)
    {
        is = i + shift_h;
        if(is >= height)
            is -= height;

        for(int j = 0;j<width;j++)
        {
            js = j + shift_w;
            if(js >= width)
                js -= width;
            shift_fft_mat[i][j] = input_mat[is][js];
        }
    }
    return shift_fft_mat;
}



vector<vector<complex<double>>> fft2(vector<vector<complex<double>>>&  input_mat, vector<vector<complex<double>>>&  tmp_mat, vector<vector<complex<double>>>&  fft_mat)
{
    //1D fft
    for(int i = 0;i<(int)input_mat.size();i++)
    {
        tmp_mat[i] = fft(input_mat[i]);
    }

    //Transpose
    complex<double> tmp;
    for(int i = 0;i<(int)tmp_mat.size();i++)
    {
        for(int j = i+1;j<(int)tmp_mat[i].size();j++)
        {
            tmp = tmp_mat[i][j];
            tmp_mat[i][j] = tmp_mat[j][i];
            tmp_mat[j][i] = tmp;
        }
    }

    //1D fft of the transpose
    for(int i = 0;i<(int)tmp_mat.size();i++)
    {
        fft_mat[i] = fft(tmp_mat[i]);
    }

    return fft_mat;
}



vector<vector<complex<double>>> fft2(vector<vector<complex<double>>>&  input_mat)
{
    vector<vector<complex<double>>>  tmp_mat, fft_mat;
    tmp_mat.resize(input_mat.size()), fft_mat.resize(input_mat.size());
    //1D fft
    for(int i = 0;i<(int)input_mat.size();i++)
    {
        tmp_mat[i] = fft(input_mat[i]);
    }

    //transpose
    complex<double> tmp;
    for(int i = 0;i<(int)tmp_mat.size();i++)
    {
        for(int j = i+1;j<(int)tmp_mat[i].size();j++)
        {
            tmp = tmp_mat[i][j];
            tmp_mat[i][j] = tmp_mat[j][i];
            tmp_mat[j][i] = tmp;
        }
    }

    //1D fft of the transpose
    for(int i = 0;i<(int)tmp_mat.size();i++)
    {
        fft_mat[i] = fft(tmp_mat[i]);
    }

    return fft_mat;
}


vector<vector<complex<double>>> fft2(M3Matrix& Im)
{
    vector<vector<complex<double>>> input_mat, fft_mat;
    input_mat.resize(Im.GetNbRows()), fft_mat.resize(Im.GetNbRows());

    //This data copy could be optimised
    for(int i = 0;i<Im.GetNbRows();i++)
        for(int j = 0;j<Im.GetNbCols();j++)
        {
            complex<double> z(Im[i][j]);
            input_mat[i].push_back(z);
        }

    return fft_mat;
}



//******************************************************
//* Class FFTW_Wrapper is a wrapper connected to the
//* fftw library
//******************************************************
FFTW_2D::FFTW_2D()
{
  m_Nx = -1;
  m_Ny = -1;
}

FFTW_2D::~FFTW_2D()
{
    fftw_destroy_plan(m_plan);
    fftw_free(m_data);
}

FFTW_2D::FFTW_2D(int Nx, int Ny)
{
    m_Nx = Nx;
    m_Ny = Ny;
    m_data = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * m_Nx * m_Ny);
    //Default : forward is computed
    m_plan = fftw_plan_dft_2d(m_Nx, m_Ny, m_data, m_data, FFTW_FORWARD, FFTW_ESTIMATE);
}


void FFTW_2D::SetSizes(int Nx, int Ny)
{
    if(Nx != m_Nx || Ny != m_Ny)
    {
        // LD: condition added to fix a segmentation fault occurring only in the make version !?!
        if(m_Nx != -1 && m_Ny != -1)
        {
            fftw_free(m_data);
        }
        m_Nx = Nx;
        m_Ny = Ny;

        m_data = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * m_Nx * m_Ny);
    }

    m_plan = fftw_plan_dft_2d(m_Nx, m_Ny, m_data, m_data, FFTW_FORWARD, FFTW_ESTIMATE);
}

void FFTW_2D::SetForward()
{
    m_plan = fftw_plan_dft_2d(m_Nx, m_Ny, m_data, m_data, FFTW_FORWARD, FFTW_ESTIMATE);
}

void FFTW_2D::SetBackward()
{
    m_plan = fftw_plan_dft_2d(m_Nx, m_Ny, m_data, m_data, FFTW_BACKWARD, FFTW_ESTIMATE);
}


void FFTW_2D::GetData(M3Matrix& Im)
{
    if(Im.GetNbRows() * Im.GetNbCols() > 0)
    {
        SetSizes(Im.GetNbRows(),Im.GetNbCols());
        //Naive, non optimized data transfer
        for(int i = 0; i < m_Nx; i++)
        {
            for(int j = 0; j < m_Ny; j++)
            {
                m_data[i*m_Ny+j][0] = Im[i][j];
                m_data[i*m_Ny+j][1] = 0.;
            }
        }

    }
}


void FFTW_2D::GetData(vector<vector<complex<double>>>& input_vec)
{
    if(input_vec.size() > 0)
    {
        SetSizes(input_vec.size(),input_vec[0].size());

        //Naive, non optimized data transfer
        for(int i = 0; i < m_Nx; i++)
        {
            for(int j = 0; j < m_Ny; j++)
            {
                m_data[i*m_Ny+j][0] = input_vec[i][j].real();
                m_data[i*m_Ny+j][1] = input_vec[i][j].imag();
            }
        }
    }
}


vector<vector<complex<double>>> FFTW_2D::fft2()
{
    SetForward();
    vector<vector<complex<double>>> fft;
    fftw_execute(m_plan);

    fft.resize(m_Nx);
    for(int i = 0; i < m_Nx; i++)
        fft[i].resize(m_Ny);

    //Naive, non optimized data transfer
    //double norm_factor = m_Nx *m_Ny;
    for(int i = 0; i < m_Nx; i++)
    {
        for(int j = 0; j < m_Ny; j++)
        {
            fft[i][j] = {m_data[i*m_Ny+j][0] ,m_data[i*m_Ny+j][1] } ;
        }
    }
    return fft;
}


void FFTW_2D::fft2(vector<vector<complex<double>>>& fft_mat)
{
    SetForward();
    fftw_execute(m_plan);

    //Naive, non optimized data transfer
    for(int i = 0; i < m_Nx; i++)
    {
        for(int j = 0; j < m_Ny; j++)
        {
            fft_mat[i][j] = {m_data[i*m_Ny+j][0] ,m_data[i*m_Ny+j][1] } ;
        }
    }
}


void FFTW_2D::ifft2(vector<vector<complex<double>>>& ifft_mat)
{
    //The following line is required to set fftw in backward modus
    SetBackward();

    //vector<vector<complex<double>>> ifft_mat;
    fftw_execute(m_plan);

    double norm_factor = 1.0/(m_Nx*m_Ny);

    //Naive, non optimized data transfer : TODO: improve performance
    for(int i = 0; i < m_Nx; i++)
    {
        for(int j = 0; j < m_Ny; j++)
        {
            ifft_mat[i][j] = {m_data[i*m_Ny+j][0]* norm_factor,m_data[i*m_Ny+j][1]* norm_factor};
        }
    }
}


vector<vector<complex<double>>> FFTW_2D::ifft2()
{
    //The following line is required to set fftw in backward modus
    SetBackward();

    vector<vector<complex<double>>> ifft_mat;
    fftw_execute(m_plan);

    ifft_mat.resize(m_Nx);
    for(int i = 0; i < m_Nx; i++)
        ifft_mat[i].resize(m_Ny);

    //Naive, non optimized data transfer
    double norm_factor = 1.0/(m_Nx*m_Ny);
    for(int i = 0; i < m_Nx; i++)
    {
        for(int j = 0; j < m_Ny; j++)
        {
            ifft_mat[i][j] = {m_data[i*m_Ny+j][0]* norm_factor,m_data[i*m_Ny+j][1]* norm_factor};
        }
    }
    return ifft_mat;
}


//**************************************
//* tests
//**************************************
void test_cholesky()
{
    M3Matrix A(4,4);
    //First block : test of CholeskyDecomp
    A={18, 22, 54, 42,
       22, 70, 86, 62,
       54, 86, 174, 134,
       42, 62, 134, 106,
      };

    //The expected result should be the following matrix:
    //  4.24264    0.00000    0.00000    0.00000
    // 5.18545    6.56591    0.00000    0.00000
    // 12.72792    3.04604    1.64974    0.00000
    // 9.89949    1.62455    1.84971    1.39262

    M3Matrix L;
    CholeskyDecomp(A,L);
    //End of first block

    cout << "L.Dump() : " << endl;
    L.Dump();

    M3Matrix b(4,1),x(4,1);
    b[0][0] = 1;
    b[1][0] = 1.201;
    b[2][0] = 1.411;
    b[3][0] = 1.7;

    CholeskySolve(A,x,b);

    cout << "b : " <<  endl;
    b.Dump();

    M3Matrix b_rec(4,1);
    b_rec = A * x;
    cout << "b_rec" << endl;
    b_rec.Dump();

    cout << endl << "**************** test for a non symmetric matrix ******************" << endl;
    A.Reshape(7,7,0);
    b.Reshape(7,1,2.0);
    double sigma = 1.0;
    A.AddGaussianWhiteNoise(sigma);
    CholeskySolve_GenericMatrix(A,x,b);

    cout << "b : " <<  endl;
    b.Dump();

    b_rec = A * x;
    cout << "b_rec" << endl;
    b_rec.Dump();
}






//Computes cumulative sums
//Note: this is the _accumulateTranslation_variable_prec in MATLAB
/*vector<double> accumulateTranslation(vector<double> tr)
{
    vector<double> tr_cum_sum(tr.size());
    for(int i = 1;i<tr.size(); i++)
    {
        tr_cum_sum[i] = tr_cum_sum[i-1] + tr[i];
    }
    return tr_cum_sum;
}*/





//*********************************************************
//* From now on tests for the fft / inverse fft and so on
//*********************************************************
void test_fft()
{
    // vector<complex<double>> x({  0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0 });
    //vector<complex<double>> x({  0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.3 });
    vector<complex<double>> x({  0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.3,-3.4 });

    vector<complex<double>> y, x_inv;
    // forward fft
    y = fft(x);
    // inverse fft
    x_inv = ifft(y);

    //Dumps
    dump_vec(x,"original\n");
    dump_vec(y,"fft\n");
    dump_vec(x_inv,"ifft\n");

    //2D FFT Tests
    cout << endl;
    cout << "*****************************" << endl;
    cout << "* 2D-FFT ********************" << endl;
    cout << "*****************************" << endl;
    M3Matrix InputIm(8,8,0.);
    InputIm.DrawRectangle(0.2,0.2,0.65,0.65,128.);

    //InputIm[0][1] = 3;
    //InputIm[1][0] = 5;
    InputIm.Dump();

    vector<vector<complex<double>>>  fft_mat = fft2(InputIm);

    //cout << "InputIm.GetVal(1) : "  << InputIm.GetVal(1) << endl;
    //cout << "InputIm.GetVal(8) : "  << InputIm.GetVal(8) << endl;

    //vector<vector<complex<double>>> input_vec;
    //input_vec

    //vector<vector<complex<double>>> input({  0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.3,-3.4 });
}


void test_fft_time()
{
    int NbRows = 128;
    int NbCols = 100;

    M3Matrix X(NbRows,NbCols,0.);
    double sigma = 1.0;
    X.AddGaussianWhiteNoise(sigma);

    vector<vector<complex<double>>> x_mat;
    x_mat.resize(X.GetNbRows());

    //This data copy could be optimised
    for(int i = 0;i<X.GetNbRows();i++)
        for(int j = 0;j<X.GetNbCols();j++)
        {
            complex<double> z(X[i][j]);
            x_mat[i].push_back(z);
        }


    int NbFFT = 21000;
    //int size = 128;

    vector<vector<complex<double>>> fft_mat, x_rec;
    vector<vector<complex<double>>> tmp_mat, tmp2_mat;
    tmp_mat = x_mat;
    tmp2_mat = x_mat;

    struct timespec Start; //, Stop;
    get_start_time(Start);

    for(int i = 0;i<NbFFT;i++)
    {
        fft_mat = fft2(x_mat,tmp_mat,tmp2_mat);
        if(i % 1000 == 0)
            cout << "i" << i  << endl;
        //TODO: implement ifft2
        //x_rec = ifft2(fft_mat);
    }

    get_elapsed_time(Start);

    //seconds = difftime(timer,mktime(&y2k));
    //printf ("%.f seconds since January 1, 2000 in the current timezone", seconds);
    //printf("computing time for %d 2D-FFT of size (%d,%d) : %f sec", NbFFT, NbRows, NbCols, seconds);
}


void test_fftw_vs_fft()
{
    cout << "test fftw_Vs_fft" << endl;
    const ptrdiff_t N0 = 9, N1 = 11;
    fftw_plan plan;
    fftw_complex *data;
    data = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N0 * N1);
    // create plan for forward DFT
    plan = fftw_plan_dft_2d(N0, N1, data, data, FFTW_FORWARD, FFTW_ESTIMATE);

    cout << endl << "original data " << endl;
    // initialize data to some function my_function(x,y)
    int i, j;
    //double pdata=0;
    for (i = 0; i < N0; ++i)
    {
        for (j = 0; j < N1; ++j)
        {
            //data[i*N1 + j][0]=i;
            //data[i*N1 + j][1]=0;

            if(i>=2 && j>=2 && i<=5 && j<=5)
                data[i*N1 + j][0]=1.0;
            else
                data[i*N1 + j][0]=0.;

            data[i*N1 + j][1]=0;

            //pdata+=data[i*N1 + j][0]*data[i*N1 + j][0]+data[i*N1 + j][1]*data[i*N1 + j][1];
            cout << "(" << data[i*N1+j][0] << " , " << data[i*N1+j][1] << "), ";
        }
        cout << endl;
    }

    // compute transforms, in-place, as many times as desired
    fftw_execute(plan);

    double normalization=sqrt((double)N0*N1);
    cout << endl << " fft" << endl;

    // normalize data and calculate power of transform
    for (i = 0; i < N0; ++i)
    {
        for (j = 0; j < N1; ++j)
        {
            data[i*N1+j][0]/=normalization;
            data[i*N1+j][1]/=normalization;

            cout << "(" << data[i*N1+j][0] << " , " << data[i*N1+j][1] << "), ";
        }
        cout << endl;
    }

    fftw_destroy_plan(plan);
    fftw_free(data);
}


void test_fftw()
{
    cout << "test fftw" << endl;
    const ptrdiff_t N0 = 120, N1 = 120;
    int NbFFT = 7200;

    fftw_plan plan;
    fftw_complex *data;

    struct timespec Start; //, Stop;

    data = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N0 * N1);
    // create plan for forward DFT
    plan = fftw_plan_dft_2d(N0, N1, data, data, FFTW_FORWARD, FFTW_ESTIMATE);

    get_start_time(Start);

    // initialize data to some function my_function(x,y)
    int i, j;
    double pdata=0;
    for (i = 0; i < N0; ++i){
        for (j = 0; j < N1; ++j){
            data[i*N1 + j][0]=i;
            data[i*N1 + j][1]=0;
            pdata+=data[i*N1 + j][0]*data[i*N1 + j][0]+data[i*N1 + j][1]*data[i*N1 + j][1];
        }
    }
    printf("power of original data is %f\n", pdata);

    // compute transforms, in-place, as many times as desired
    for(int i=0;i<NbFFT;i++)
        fftw_execute(plan);

    double elapsed_time =  get_elapsed_time(Start);
    cout << NbFFT << " fftw (2D)  of size (" <<  N0 <<","<< N1 << ") : " << elapsed_time << " sec "<< endl;

    double normalization=sqrt((double)N0*N1);
    double ptransform = 0;

    // normalize data and calculate power of transform
    for (i = 0; i < N0; ++i){
        for (j = 0; j < N1; ++j){
            data[i*N1+j][0]/=normalization;
            data[i*N1+j][1]/=normalization;
            ptransform+=data[i*N1 + j][0]*data[i*N1 + j][0]+data[i*N1 + j][1]*data[i*N1 + j][1];
        }
    }

    printf("power of transform is %f\n", pdata);

    fftw_destroy_plan(plan);
    fftw_free(data);


    //*********************************************
    //*********************************************
    //* Test FFTW_Wrapper class
    //*********************************************
    //*********************************************
    M3Matrix X(32,32,0.);
    double sigma = 1.0;
    X.AddGaussianWhiteNoise(sigma);

    vector<vector<complex<double>>> x_mat;
    x_mat.resize(X.GetNbRows());

    //This data copy could be optimised
    for(int i = 0;i<X.GetNbRows();i++)
        for(int j = 0;j<X.GetNbCols();j++)
        {
            complex<double> z(X[i][j]);
            x_mat[i].push_back(z);
        }


    FFTW_2D fftww(N0,N1); //must contain
    fftww.GetData(x_mat);

    //M3Matrix InIm;
    vector<vector<complex<double>>> fft = fftww.fft2();
}


//TODO: COPIER DANS UN AUTRE FICHIER ET TESTER
/*void LoadTiff(string filename,vector<M3Matrix>& RGB)
{
    cout << "LoadTiff" << endl;
    TIFF* tif = TIFFOpen(filename.c_str(),"r");
    cout << "TIFFOpen worked" << endl;
    uint32 w, h;
    tdata_t buf;

    if(tif)
    {

        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
        cout << "w : " << w << " h : " << h << endl;
        //int linesize = TIFFScanlineSize(tif);
        int nsamples = 3;
        TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &nsamples);
        cout << "nsamples : " << nsamples << endl;

        RGB.resize(nsamples);
        for(int c=0;c<nsamples;c++)
        {
            RGB[c].Reshape(h,w);
        }

        int bitspersample = 0;
        if (TIFFGetFieldDefaulted(tif, TIFFTAG_BITSPERSAMPLE, &bitspersample) == 0)
        {
            // assume 8 bit
            bitspersample = 8;
        }
        
        if(bitspersample == 16)
        {
            for(int c=0;c<nsamples;c++)
            {
                RGB[c].InitialiseBitDepth(bitspersample);
            }
        }



        uint16 config;
        TIFFGetField(tif, TIFFTAG_PLANARCONFIG, &config);
        //if(config == PLANARCONFIG_CONTIG)
       //{
       //   cout <<  "PLANARCONFIG_CONTIG recognized" << endl;
       //}


        buf = _TIFFmalloc(TIFFScanlineSize(tif));

        for(int i = 0; i < (int)h; i++) //image rows
        {
            TIFFReadScanline(tif, buf, i);
            for(int col = 0;col < (int)w;col++)
            {
                if(nsamples == 3)
                {
                    uint16 red  = static_cast<uint16*>(buf)[col*nsamples+0];
                    uint16 green  = static_cast<uint16*>(buf)[col*nsamples+1];
                    uint16 blue = static_cast<uint16*>(buf)[col*nsamples+2];

                    RGB[0][i][col] = red;
                    RGB[1][i][col] = green;
                    RGB[2][i][col] = blue;
                }
                else if (nsamples == 1)
                {
                    uint16 val = static_cast<uint16*>(buf)[col*nsamples+0];
                    RGB[0][i][col] = val;
                }
            }
        }


        _TIFFfree(buf);
        //TODO (04 july 2019): see why this crashed
        TIFFClose(tif);
    }
}*/


/*void test_tiff()
{
    cout << "Test of the tiff5 library " << endl;
    string filename("/home/ash/Desktop/Test/C++/TIFF/22000.tif");
    vector<M3Matrix> ColorIm;
    LoadTiff(filename,ColorIm);

    string filenameOut("/home/ash/Desktop/Test/C++/TIFF/22000_saved.ppm");
    SavePPM(filenameOut,ColorIm[0],ColorIm[1],ColorIm[2]);
}*/


void test_ifft2()
{
    int  NbRows=8,NbCols=8;

    M3Matrix X(NbRows,NbCols);

    for(int i = 0;i<NbRows;i++)
        for(int j = 0;j<NbRows;j++)
        {
            X[i][j] = i*NbCols+j+1;
        }

    vector<vector<complex<double>>> x,f,x_rec;
    FFTW_2D fftww(NbRows,NbCols);
    cout << "FFTW_2D created " << endl;
    fftww.GetData(X);
    cout << "Data in fftw " << endl;

    cout << "X.Dump()" << endl;
    X.Dump();

    f = fftww.fft2();
    x_rec = fftww.ifft2();

    cout << "x_rec : " << endl;
    dump_vec(x_rec);
    cout << endl;
}



#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <stdlib.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdint>
#include <vector>
#include <iomanip>

using namespace std;

#define DBL2INT(x)       ((int)(x+0.000000000001))
#define DBL_ROUND(x)     ((double)(int)(x+0.000000000001))

#define TOL 0.00000001

//color pixel  transformation
void RGBtoYUV( double r, double g, double b, double &y, double &u, double &v);
void YUVtoRGB( double y, double u, double v, double &r, double &g, double &b);
void RGBtoYCoCg( double r, double g, double b, double &y, double &co, double &cg);

double compute_median(double* x, int length);

enum TXTDIRECTION{vertical,horizontal};
enum LAPLACIANTYPE {isotropic,anisotropic,lothar};

class M3Matrix3D;

class LUTParams
{
public:
  LUTParams();
  virtual ~LUTParams();
  double LUT_f(double x);

  int getLUTType() const {return m_LUTType; }
  void setLUTType(const int LUTType) {  m_LUTType = LUTType;}

  double gamma;
  double linear_slope;
  double linear_shift;

private:
  int m_LUTType; //0:gamma, 1: linear multiplicative
};


class BilParams //This is a simple container class which contains the parameters used in bilateral filtering
{
public:
  BilParams();
  virtual ~BilParams ();

  BilParams& operator=(const BilParams& Src);

  int getAWType() const {return m_AWType; }
  void setAWType(const int AWType) {  m_AWType = AWType;}

  enum class NormType { off, constant, exact, approx, fpga };

  //These are the parameters for the adaptive weight function phi = a + b/x
  // if not LUT, this corresponds to the fpga codes psi = ax +b (currently used together with option)
  double Threshold = 1.0;
  double a = 0.0;
  double b = 0.0;
  double Dynamic = 255.0;
  NormType norm_type = NormType::exact;
  int Border_Type = 0; //Option to indicate how the boundary is extended; 0: mirroring (default), TODO: implement it !
  int nbBits = 8;

  bool GreenOnly = false;
  int BPattern = 0; //Only used when green only is on

  //These are the thresholds used when considering the bilateral ratio
  double t1_br = 0.0;
  double t2_br = 0.0;

  double AW_f(double x);

  int add_bits = 0; //this is used to add bit precision

private:
  int m_AWType; //0:
};

class Image;

class M3Matrix
{
public:
  // Constructors
  M3Matrix();
  M3Matrix(int NbRows, int NbCols);
  M3Matrix(int NbRows, int NbCols, const double& r_Value);
  M3Matrix(int NbRows, int NbCols, const vector<double>& vals);
  M3Matrix(const vector<double>& vals);
  M3Matrix(const Image &im);
  M3Matrix(string image_name);

  // Copy
  M3Matrix(const M3Matrix& r_Source);

  // Destructor
  virtual ~M3Matrix();

  // Accessors
  int GetNbRows() const;
  int GetNbCols() const;
  int GetNbBits() const {return t_NbBits;}
  int GetRangeMax() const  {return (int)t_Max;}
  int GetRangeMin() const {return (int)t_Min;}
  double* GetData() {return Data;}
  void SetData(int NbRows, int NbCols, double* x);
  void SetData(double* x);
  double* operator[](int Row);
  double* operator[](int Row) const;

  double GetVal(int indx);

  // Extraction of rows and columns
  void SetLine(int j, M3Matrix& Row);
  void SetCol(int j, M3Matrix& Col);

  // Transformers
  void Set2Id(int n);
  void SetValues(double const& Value);
  void SetValues(vector<double> vals);

  void Reshape(int NbRows, int NbCols, double Val = 0.);
  void CatMat(M3Matrix& A,M3Matrix& B,int Dim = 2);

  void Initialise(const M3Matrix &r_Source);
  void InitialiseDynamic(const M3Matrix& r_Source);

  void Transpose();
  void Transpose(M3Matrix& AT);
  void GetRevertRows(M3Matrix B);
  void GetRevertCols(M3Matrix B);

  // Operators
  void operator=(const vector<double>& flat);
  M3Matrix& operator=(const M3Matrix& Src);

  M3Matrix operator+(const M3Matrix& B);
  M3Matrix operator*(const M3Matrix& B);
  M3Matrix operator-(const M3Matrix& B);

  //operator right multiplication by a vector of double
  vector<double> operator*(const vector<double>& b);
  vector<double> operator+(const vector<double>& v1);
  vector<double> operator-(const vector<double>& v1);

  //operations with scalars
  M3Matrix operator+(const double& t);
  M3Matrix operator*(const double& t);
  M3Matrix operator/(const double& t);
  M3Matrix operator-(const double& t);

  void operator+=(const double& t);
  void operator*=(const double& t);
  void operator/=(const double& t);
  void operator-=(const double& t);

  void operator+=(M3Matrix& B);
  void operator-=(M3Matrix& B);
  void operator*=(M3Matrix& B);

  void PointMult(M3Matrix& A,M3Matrix& B);

  M3Matrix GradientH();
  M3Matrix GradientV();

  M3Matrix GradientH_Centered();
  M3Matrix GradientV_Centered();

  // Display - save
  void DumpSize();


  void DumpBlock(int rStart, int cStart, int rEnd, int cEnd, ostream& Out = cout, int setw_val = 4, double tol = 0.0);
  void Dump(ostream& Out = cout, int setw_val = 4, double tol = 0.0);
  void DumpRGB(ostream& Out = cout, int setw_val = 4);
  void Dump_Binary(ostream& Out = cout);
  void DumpRGB_Binary(ostream& Out = cout);
  void SavePPM(const string& filenameOut);
  void SavePPM(ostream& Out);
  void SavePGM(const string& filenameOut);
  void SavePGM(ostream& Out);
  void SaveOff(const string& filenameOut, double factor = 50.);
  void SaveOff(ostream& Out, double factor = 50.);

  void SavePGM_Binary(const string& filenameOut);
  void SavePPM_Binary(const string& filenameOut);
  void SavePGM_Binary(ostream& Out);
  void SavePPM_Binary(ostream& Out);

  void Dump_From_Data();
  
  void ExchangeRows(int i1= 0, int i2=0);
  void ExchangeCols(int i1, int i2);
  
  // Sorting functions
  void SortByCol(int c);
  void SortByCols(int c1, int c2);
  void SortByBlocs(int c1, int c2);

  // Gets the minimal value of the array
  double GetMin();
  double GetMax();
  int GetDynamic();

  void Clip();

  // Quantization of the coefficients
  void Quantize(double& QStep, M3Matrix& Symbol);
  //Round the values to int and cast between 0 and the dynamic
  void Digitalize();

  void SetBitDepth(unsigned int NbBitsOut);
  void InitialiseBitDepth(int NbBitsOut);

  void SetNbBits(int NbBits) {t_NbBits = NbBits; return;}
  void SetRangeMax(double RangeMax) {t_Max = RangeMax; return;}
  void SetRangeMin(double RangeMin) {t_Min = RangeMin; return;}
  void Rescale(double factor);

  void Negative(M3Matrix&  Negative);

  // Standard computation tools
  double Median();
  double MedMed();

  double GetSum();
  void Normalise();
  double LpNorm(double p);
  double L0Norm();
  double Kurtosis(int Margin = 0);
  double MeanKurtosisOfGradients(int Margin = 0);
  double ScaledMeanKurtosis(int Margin = 0);

  double MeanKurtosisOfHGradient(int Margin = 0);
  double MeanKurtosisOfVGradient(int Margin = 0);

  //Stastistical tools
  void Autocorrelation(M3Matrix& A,double mux = 0., int MaxShift=3);
  void Correlation(M3Matrix& A,M3Matrix& B, double mux = 0., double muy = 0., int MaxShift=3);
  
  double MSE(M3Matrix& M3Matrix2, int Margin = 0);
  double PSNR(M3Matrix& M3Matrix2, int Margin = 0);
  double ColorPSNR(M3Matrix& RGB, int Margin);

  double SSIM(M3Matrix& B,int WindowSize, double& lum, double& contrast, double& structure, double K1=0.1,double K2=0.3);
  double SobolevNorm(int p = 1);
  double Sobolev2Norm(int p = 1);
  double MeanOfMaxGradient(double p = 0.01);

  double Histogram(M3Matrix& Histo);

  // Loading the M3Matrix from a file
  void LoadFromTXT(const string& filenameIn);
  void WriteTXT(const string& filenameOut, TXTDIRECTION DataDirection=vertical);

  void ParsePXMMetaData(string line_str);

  void LoadFromXRGB(const char* infile);
  void LoadFromXRGB(const string&  infile);

  void LoadRawFromXRGB(const char* infile);
  void LoadRawFromXRGB(const string& infile);

  void LoadFromPGM(const char* infile);
  void LoadFromPGM(const string& infile);
  void LoadFromPPM(const char* infile);
  void LoadFromPPM(const string& infile);

  void ExtractRGBChannels(M3Matrix& R,M3Matrix& G,M3Matrix& B);
  void CombineRGBChannels(M3Matrix& R,M3Matrix& G,M3Matrix& B);
  void CombineInterlacedRGBChannels(M3Matrix& R,M3Matrix& G,M3Matrix& B);

  //Image Processing toolbox
  M3Matrix SimulateRawBayer(int BayerPattern = 0);
  void SimulateRawBayer(M3Matrix& R,M3Matrix& G,M3Matrix& B, int BayerPattern = 0);

  void HAJapanDebayering(M3Matrix& R,M3Matrix& G,M3Matrix& B, int BayerPattern = 0, int ThresholdJapan = 0, int epsilon = 5); //Akishima version of Hamilton-Adams
  void HAPrefDebayering(M3Matrix& R,M3Matrix& G,M3Matrix& B, int BayerPattern = 0); //Prefiltered Hamilton-Adams
  void HADebayering(M3Matrix& R,M3Matrix& G,M3Matrix& B, int BayerPattern = 0);
  M3Matrix HAJapanDebayering(int BayerPattern = 0, int Threshold = 0, int ThresholdJapan = 5); //Akishima version of Hamilton-Adams
  M3Matrix HADebayering(int BayerPattern = 0);
  M3Matrix HAPrefDebayering(int BayerPattern = 0);


  void ZhangWuDebayering(M3Matrix& R,M3Matrix& G,M3Matrix& B, int BayerPattern = 0);
  void ExperimentalDebayering(M3Matrix& R,M3Matrix& G,M3Matrix& B, int BayerPattern = 0);
  void Debayering(M3Matrix& Image_Debayered);
  void DPC(M3Matrix& Corrected, M3Matrix& Mask, double threshold, int BayerPattern = 0, int opt = 0);

  void PixelTransfo(M3Matrix& TransformedImage, LUTParams& Parameters);
  void GammaCorrection(M3Matrix& GammaCorrectedImage, double gamma, int Max = 255);

  void Shift(M3Matrix& ShiftedImage, int tx, int ty);
  void Extend(M3Matrix& MirroredImage, int Height, int Width, int ExtensionType=0);
  void Mirror(M3Matrix& MirroredImage, int Height, int Width); //Warning: this is only valid for grayscale images

  void AddBlock(M3Matrix& BlockMatrix, int rStart, int cStart);
  void SetBlock(M3Matrix& BlockMatrix, int rStart, int cStart);
  void ExtractBlock(M3Matrix& BlockMatrix, int rStart, int cStart, int NbRows,int NbCols);
  M3Matrix ExtractBlock(int rStart, int cStart, int NbRows,int NbCols);


  void LinearTemporalFilter(M3Matrix& FilteredImage, M3Matrix3D& Mask3D, M3Matrix3D& CircPast, int CircIndx=0/*, int NbPastFrames=0*/);
  void LinearFilter(M3Matrix& FilteredImage, M3Matrix& Mask);
  void LinearFilter(M3Matrix& Mask);
  M3Matrix LinearFiltering(M3Matrix& Mask);
  M3Matrix MaxFilter(int rSize, int cSize);
  void FilterMax(int rSize, int cSize);

  void LinearFilter_Ver_DS(M3Matrix& FilteredImage, M3Matrix& Mask);
  void LinearFilter_Hor_DS(M3Matrix& FilteredImage, M3Matrix& Mask);

  void DoG(M3Matrix& B,double sigmaLow=1,double sigmaHigh=2);
  void GaussianBlur(double sigma=1);
  M3Matrix GaussianBlurring(double sigma=1);

  void DifferenceImage(M3Matrix& A, M3Matrix&B, double Factor=1.0, double Dynamic=255);

  //Noise simulation
  void AddGaussianWhiteNoise(double sigma, int seed = 0);
  void AddGaussianShotNoise(double sigma=1, int seed = 0);
  void AddImpulsiveNoise(double sigma=1, double prob = 0.02, int seed = 0);

  void AddGaussianWhiteNoise(M3Matrix& NoisyImage, double sigma=1, int seed = 0);
  void AddGaussianShotNoise(M3Matrix& NoisyImage, double sigma=1, int seed = 0);
  void AddImpulsiveNoise(M3Matrix& NoisyImage, double sigma=1, double prob = 0.02, int seed = 0);

  void AddImpulsiveNoise(M3Matrix& NoisyImage, M3Matrix& ImpulsiveMask, double sigma=1, double prob = 0.02, int seed = 0);

  //Local Fractal Dimension
  M3Matrix LocalFractalDimension();

  //Forward shading operator
  M3Matrix ForwardShading();

  //Some classical filters
  //Edge detectors ("differential" operators)
  void MakeGradientC();
  void MakeGradientR();

  void MakePrewittC();
  void MakePrewittR();
  void MakeSobelC();
  void MakeSobelR();

  void MakeKirschFilter(int KirschIndx=1);
  void MakeFreiChenFilter(int FreiChenIndx=1);

  void MakeLaws1D(int LawsIndx);
  void MakeLawsFilter(int VertLawsIndx, int HorLawsIndx);

  void MinMedianMaxFilter(int S=5, double threshold = 10., double factor = 1.0);

  void DirectionalMedianFilter(M3Matrix& Filtered, int Size);

  void MedianFilter(M3Matrix& MedianFilteredImage, int Height=3, int Width=3);
  void MedMedFilter(M3Matrix& MedianFilteredImage, int Height=3, int Width=3);

  //Smoothing (averaging) kernels
  void MakeRampFilter(int Height,int Width);

  void MakeBoxFilter(int Height,int Width);
  void MakeBoxFilter(double Value, int Height,int Width);
  void MakeSparseBoxFilter(double Value, int Height,int Width, int SamplingStep=2);
  
  void MakeCircBoxFilter(double Value, int Height,int Width);
  void MakeLaplaceFilter(LAPLACIANTYPE laplaciantype=isotropic);
  void MakeGaussianFilter3(); //discrete version of a gaussian filter of size 3x3
  void MakeGaussianFilter5(); //discrete version of a gaussian filter of size 5x5
  void MakeGaussianFilter7(); //discrete version of a gaussian filter of size 7x7
  
  void MakeGaussianFilter7_1D(); //in order to test on signals

  void MakeGaussianFilter(double sigma, int Size = 11);
  void MakeGaborRealFilter(double sigma, int Size = 11, double omega=1, double theta=0, double gamma = 1);
  void MakeGaborImagFilter(double sigma, int Size = 11, double omega=1, double theta=0, double gamma = 1);

  void MakeGaussianFilter_1D(double sigma, int Size = 7);
  
  //Filter manipulation
  void ChessBoard();
  void Sparsify(int SamplingStep);
  void Hollow();
  
  void Rescale(M3Matrix& RescaledMatrix, double MinVal=0, double MaxVal=255);
  void LaplaceEnhancement(M3Matrix& FilteredImage, double alpha = 1.0);
  void UnsharpMasking_Box(M3Matrix& FilteredImage, int Height=3,int Width=3);

  // Wedge filters
  void MakeWedgeLeft(int Size);
  void MakeWedgeRight(int Size);
  void MakeWedgeTop(int Size);
  void MakeWedgeBottom(int Size);
  void MakeWedgeNE(int Size);
  void MakeWedgeNW(int Size);
  void MakeWedgeSW(int Size);
  void MakeWedgeSE(int Size);
  
  //Computes the variance of the "active" wedge
  double WedgeAbsDev(M3Matrix& WedgeFilter,int i,int j);
  double WedgeVariance(M3Matrix& WedgeFilter,int i,int j);
  double WedgeDeviation(M3Matrix& WedgeFilter,int i, int j, int Criterion = 0);
  void SelectiveWedgeFilter(M3Matrix& WedgeFilteredImage, int Size, int FilterType = 0, double Tau = 0.5, int Criterion = 0, double Sigma= 1.0);

  //Reference bilateral filter
  void AdaptiveBilateralFilter(M3Matrix& FilteredImage, M3Matrix& Mask, BilParams Parameters);

  //Experimental bilateral filters
  void AdaptiveBilateralFilter_Median(M3Matrix& FilteredImage, int Size, double Threshold, double Dynamic = 255); //Bilateral filter with median as underlying filter
  void AdaptiveBilateralFilter_MedMed(M3Matrix& FilteredImage, int Size, double Threshold, double Dynamic = 255); //Bilateral filter with median of median as underlying filter
  void AdaptiveBilateralFilter_EdgeMasked(M3Matrix& FilteredImage, M3Matrix& Mask, M3Matrix& EdgeMask, BilParams Parameters); //
  void AdaptiveBilateralFilter_BilRatio(M3Matrix& FilteredImage, M3Matrix& Mask, BilParams Parameters);
  void BilateralTemporalFilter(M3Matrix& FilteredImage, M3Matrix3D& Mask3D, M3Matrix3D& Buffer, int CircIndx, BilParams Parameters);

  //Hessian kernel
  void Hessian(M3Matrix3D& H, double sigma);

  //Wavelet tools
  void UndecimatedWaveletFilterRow(M3Matrix& SubbandMatrix, M3Matrix& LowFilter, M3Matrix& HighFilter);
  void UndecimatedWaveletFilterCol(M3Matrix& SubbandMatrix,M3Matrix& LowFilter,M3Matrix& HighFilter);
  void WaveletFilterRow(M3Matrix& SubbandMatrix,M3Matrix& LowFilter, M3Matrix& HighFilter);
  void WaveletFilterCol(M3Matrix& SubbandMatrix,M3Matrix& LowFilter,M3Matrix& HighFilter);
  void InverseWaveletFilterRow(M3Matrix& SubbandMatrix, M3Matrix& LowFilter,M3Matrix& HighFilter);
  void InverseWaveletFilterCol(M3Matrix& SubbandMatrix,M3Matrix& LowFilter,M3Matrix& HighFilter);

  void InverseUndecimatedWaveletFilterRow(M3Matrix& SubbandMatrix,M3Matrix& LowFilter,M3Matrix& HighFilter);
  void InverseUndecimatedWaveletFilterCol(M3Matrix& SubbandMatrix,M3Matrix& LowFilter,M3Matrix& HighFilter);

  //decimated wavelet transform
  void WaveletDecomposition(M3Matrix& SubbandMatrix,M3Matrix& LowFilter,M3Matrix& HighFilter,int NbDecomp);
  void InverseWaveletDecomposition(M3Matrix& SubbandMatrix,M3Matrix& LowFilter,M3Matrix& HighFilter,int NbDecomp);
  void GenerateWaveletVisualisation(M3Matrix& Visualisation, int NbDec);

  //undecimated wavelet transform
  void UndecimatedWaveletDecomposition(M3Matrix3D& USubband,M3Matrix& LowFilter,M3Matrix& HighFilter,int NbDecomp);
  void InverseUndecimatedWaveletDecomposition(M3Matrix3D& USubband,M3Matrix& LowFilter,M3Matrix& HighFilter,int NbDecomp);


  void WaveletEnhancement(double AmplificationFactor,int NbDec);
  int WaveletHardThreshold(double Threshold, int NbDecomp);
  int WaveletSoftThreshold(double Threshold, int NbDecomp);
  int WaveletThreshold(double Threshold,int NbDecomp, int ThresholdType);
  int WaveletGarroteThreshold(double Threshold,int NbDecomp);

  int HardThreshold(double Threshold);
  int SoftThreshold(double Threshold);
  int Threshold(double Threshold, int ThresholdType=1); //default hard thresholding
  int GarroteThreshold(double Threshold);

  //Geometric features
  void AdditiveShadowing(double a);
  void AdditiveShadowing(M3Matrix& Shadowed, double a = 45);
  void MultiplicativeShadowing(M3Matrix& Shadowed,double a = 45);
  void MultiplicativeShadowing(double a = 45);

  void SlantedEdge(double theta = 0.4);

  void DrawRandomCrackTree(double xStart, double yStart, double Value);
  void DrawLine(double xStart,double yStart,double Value, double theta, double Length,double Width);
  void DrawEllipticCrown(double r1, double r2, double cx, double cy, double a, double b, double theta, double Value);
  void DrawCrown(double r1, double r2, double cx, double cy, double Value);
  void DrawCircle(double r=0.25, double cx=0.5, double cy=0.5, double Value=128.);
  void DrawEllipse(double r, double cx=0.5, double cy=0.5, double a=1.0, double b=1.0, double theta=0., double Value=128.);
  void DrawRectangle(double ulcx, double ulcy, double lrcx, double lrcy, double Value=128.);
  void DrawDiagonalBand(double xStart, double Width, double Value=128., double Slope = 1);
  void DrawCapillary(double Value);
  void AddCapillary(double Factor,double xShift=0,double yShift=0);
  void AddBlob(double r, double cx=0.5, double cy=0.5, double Value=100);
  void AddGaussianBlob(double cx=0.5, double cy=0.5, double sigma=1, double MaxValue=100);

  void Eig2(double& lambda1, double& lambda2, M3Matrix& e1, M3Matrix& e2);

  string getBayerPattern() const { return m_BayerPattern; }
  void setBayerPattern(const string BayerPattern) { m_BayerPattern = BayerPattern; }
  void dumpBayerPattern(){cout << "m_BayerPattern : " << m_BayerPattern << endl; }

  double MeshGrid(int NbXSamples=11, int NbYSamples=11, double xMin=-1., double xMax=1., double yMin=-1., double yMax=1.);

  void SolveLinear(M3Matrix& A, M3Matrix&  b, M3Matrix& Ainv);

protected:
  unsigned int t_NbBits;
  double t_Min;
  double t_Max;
  int t_NbRows;
  int t_NbCols;
  double *Data = nullptr;
  string m_BayerPattern="Undefined";

private:
  // Memory deallocation
  void Alloc();
  void Free();

  // Sorting tools
  int vCompareFirstCol(const void *a1, const void *a2);
};


class M3Matrix3D
{
public:
  // Constructors
  M3Matrix3D();
  M3Matrix3D(int NbFrames);
  M3Matrix3D(int NbFrames, int NbRows, int NbCols);
  M3Matrix3D(int NbFrames, int NbRows, int NbCols,  const double& Value);

  // Destructor
  virtual ~M3Matrix3D();

  // Accessors
  int GetNbRows();
  int GetNbCols();
  int GetNbFrames();

  double  GetSum();
  void Normalise();

  //Operators
  M3Matrix& operator[](int Frame);

  void operator+=(const double& t);
  void operator*=(const double& t);
  void operator/=(const double& t);
  void operator-=(const double& t);

  M3Matrix3D operator/(const double& t);

  // Transformers
  void SetValues(const double &Value);
  void Reshape(int NbFrames);
  void Reshape(int NbRows, int NbCols);
  void Reshape(int NbFrames, int NbRows, int NbCols);

  // Display - save
  void Dump(ostream& Out = cout);

  //Visualisation tool for the undecimated wavelet transform
  void GenerateWaveletVisualisation(M3Matrix& Visualisation, int NbDec);

  void WaveletEnhancement(double AmplificationFactor, int NbDec);

protected:
  M3Matrix* m_ImageStack;
  int t_NbFrames;
};


int CholeskyDecomp(M3Matrix& A, M3Matrix& L);
void CholeskySolve(M3Matrix&  A,  M3Matrix& x, M3Matrix& b);
void CholeskySolve_GenericMatrix(M3Matrix&  A,  M3Matrix& x, M3Matrix& b);

double ColorPSNR(M3Matrix& R1, M3Matrix& G1, M3Matrix& B1, M3Matrix& R2, M3Matrix& G2, M3Matrix& B2, int Margin);
double PSNR(M3Matrix& M3Matrix1, M3Matrix& M3Matrix2, int Margin);

M3Matrix HannPadded(int n, int m,  double factor = 0.1);
M3Matrix HadamardProduct(M3Matrix& A, M3Matrix& B);
void HadamardProduct_NoAlloc(M3Matrix& A, M3Matrix& B,M3Matrix& C);

M3Matrix Transpose(M3Matrix& A);
M3Matrix Diag(vector<double> d);

//Left multiplication
M3Matrix operator*(const double& scalar, M3Matrix A);
M3Matrix operator+(const double& scalar, M3Matrix A);
M3Matrix operator-(const double& scalar, M3Matrix A);

//vector / Matrix operators
vector<double> operator+(  const vector<double>& v1, const M3Matrix& V2);
vector<double> operator-(  const vector<double>& v1, const M3Matrix& V2);

//Synthetic generators of signals and images
void generate_synthetic_image(int option = 0);
void generate_color_image();
void generate_and_test_line_image(int option = 0);

//Wavelet filter generation
void MakeWavelet(M3Matrix& LowFilter,M3Matrix&  HighFilter,M3Matrix&  LowInvFilter,M3Matrix&  HighInvFilter,M3Matrix&  LowFilterL1, M3Matrix& HighFilterL1, int WaveletType = 0);
void MakeHaar(M3Matrix& LowFilter,M3Matrix&  HighFilter,M3Matrix&  LowInvFilter,M3Matrix&  HighInvFilter,M3Matrix&  LowFilterL1, M3Matrix& HighFilterL1);
void MakeDaub4(M3Matrix& LowFilter,M3Matrix&  HighFilter,M3Matrix&  LowInvFilter,M3Matrix&  HighInvFilter,M3Matrix&  LowFilterL1, M3Matrix& HighFilterL1);

//Color space transformations
void RGBtoYUV(M3Matrix& R,M3Matrix& G,M3Matrix& B,M3Matrix& Y,M3Matrix& U,M3Matrix& V);
void YUVtoRGB(M3Matrix& Y, M3Matrix& U, M3Matrix& V, M3Matrix& R,M3Matrix& G,M3Matrix& B);
void RGB2YUVDCT(const M3Matrix &R, const M3Matrix &G, const M3Matrix &B, M3Matrix& Y, M3Matrix& U, M3Matrix& V);
void YUVDCT2RGB(const M3Matrix &Y, const M3Matrix &U, const M3Matrix &V, M3Matrix& R, M3Matrix& G, M3Matrix& B);

void SavePPM(string& filenameOut, M3Matrix& R, M3Matrix& G, M3Matrix& B, int NbBits = -1);
void SavePPM_Binary(string& filenameOut, M3Matrix& R, M3Matrix& G, M3Matrix& B, int NbBits = -1);


//Useful functions for bilateral filtering
double ComputeAdaptiveWeight(double AbsDiff, int Option = 0, double Threshold = 0.5, double Dynamic=255.0);
double ComputeAdaptiveWeight(double AbsDiff, BilParams& Parameters);
double ComputeAdaptiveWeight2(double AbsDiff, BilParams& Parameters);
double ComputePsi(double Diff, BilParams& Parameters);

//Tests
void test_for_median();
void test_digitalize();
void test_rgb();
void test_yuv();
void test_matrix_operators();
void test_formats();
void test_randn();
void test_exact_rgb2hsv();

void test_line_image(M3Matrix& Signal,string& filename);
void system_engineering_test();
void test_RGB2XYZ_matrices();
void test_mat_vec_mult();

#endif /* of _MATRIX_H_ */

#include "utils.h"
#include "matrix.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <assert.h>

#ifndef M_PI
#define M_PI		3.14159265358979323846
#endif

//***************************************
// Author: Laurent Demaret, 2010-2025 ***
//***************************************
using namespace std;

//Default constructor / destructor
LUTParams::LUTParams() {}
LUTParams::~LUTParams() {}

double LUTParams::LUT_f(double x)
{
    double val = linear_slope * pow(fabs(x),gamma) + linear_shift;
    val = myclip(val,0,1);
    return val;
}


//Default constructor / destructor
BilParams::BilParams() {}
BilParams::~BilParams() {}

BilParams& BilParams::operator=(const BilParams& Src)
{
    (*this).a = Src.a;
    (*this).b = Src.b;
    (*this).Dynamic = Src.Dynamic;
    (*this).setAWType(Src.getAWType());
    (*this).nbBits = Src.nbBits;

    (*this).GreenOnly = Src.GreenOnly;
    (*this).BPattern = Src.BPattern;
    (*this).add_bits = Src.add_bits;

    return *this;
}


double BilParams::AW_f(double x)
{
    double val; //= x;
    if (x == 0)
    {
        val = 1;
    }
    else
    {
        val = a + b / x;
        if (val > 1)
            val = 1.;
    }

    val = myclip(val, 0, 1);

    return val;
}


//*****************************
//* M3Matrix class
//*****************************
//Default constructor
M3Matrix::M3Matrix()
{
    t_NbRows = 0;
    t_NbCols = 0;
    //Default values for bit depths and range
    t_NbBits = 8;
    t_Min = 0;
    t_Max = 255;

    Data = NULL;
}


// Constructor with dimensions
M3Matrix::M3Matrix(int NbRows, int NbCols)
{
    t_NbRows = NbRows;
    t_NbCols = NbCols;

    //Default values for bit depths and range
    t_NbBits = 8;
    t_Min = 0;
    t_Max = 255;

    Alloc();
    SetValues(0.);
}


//Constructor with dimensions and vector
M3Matrix::M3Matrix(const vector<double>& vals)
{
  Reshape(vals.size(),1);
  SetValues(vals);
}


//Constructor with dimensions and vector
M3Matrix::M3Matrix(int NbRows, int NbCols, const vector<double>& vals)
{
  Reshape(NbRows,NbCols);
  SetValues(vals);
}


// Constructor with dimensions and initial value
M3Matrix::M3Matrix(int NbRows, int NbCols, const double& r_Value)
{
    t_NbRows = NbRows;
    t_NbCols = NbCols;

    //Default values for bit depths and range
    t_NbBits = 8;
    t_Min = 0;
    t_Max = 255;

    Alloc();
    SetValues(r_Value);
}


M3Matrix::M3Matrix(string image_name)
{
    if(extract_ext(image_name) == "pgm")
        (*this).LoadFromPGM(image_name);
    else if(extract_ext(image_name) == "ppm")
    {
        (*this).LoadFromPPM(image_name);
    }
    else
    {
        t_NbRows = 0;
        t_NbCols = 0;
        //Default values for bit depths and range
        t_NbBits = 8;
        t_Min = 0;
        t_Max = 255;

        Data = NULL;
    }
}


// Copy constructor : r_Source copied into this
M3Matrix::M3Matrix(const M3Matrix& r_Source)
{
    t_NbRows = r_Source.t_NbRows;
    t_NbCols = r_Source.t_NbCols;

    InitialiseDynamic(r_Source);

    Alloc();
    if(Data)
    {
        for(int i=0;i<t_NbRows;i++)
            for(int j=0;j<t_NbCols;j++)
                (*this)[i][j] = r_Source.Data[i*t_NbCols+j];
    }
}


// Destructor
M3Matrix::~M3Matrix()
{
    Free();
    Data = NULL; //No clean memory
}


// Deallocation
void M3Matrix::Free()
{
    free(Data);
    Data = 0;
}


// Copy constructor : r_Source copied into this
void M3Matrix::Initialise(const M3Matrix& Src)
{
  t_NbRows = Src.t_NbRows;
  t_NbCols = Src.t_NbCols;

  InitialiseDynamic(Src);
  Alloc();
}


void M3Matrix::InitialiseDynamic(const M3Matrix& Src)
{
    t_Max = Src.t_Max;
    t_Min = Src.t_Min;
    t_NbBits = Src.t_NbBits;
}


//Concatenation of two matrices
//dim = 1: along the first dimension (i.e. vertical)
//dim = 2: along the second dimension (i.e. horizontal), default
void M3Matrix::CatMat(M3Matrix& A,M3Matrix& B,int dim)
{
    int NbRows,NbCols;
    if(dim == 1) //vertical concatenation
    {
        NbRows = A.GetNbRows() + B.GetNbRows();
        NbCols = A.GetNbCols();
        if(B.GetNbCols()<NbCols)
            NbCols = B.GetNbCols();
    }
    else if(dim == 2)
    {
        NbCols = A.GetNbCols() + B.GetNbCols();
        NbRows = A.GetNbRows();
        if(B.GetNbRows()<NbRows)
            NbRows = B.GetNbRows();
    }
    else
    {
        return;
    }

    (*this).Reshape(NbRows,NbCols);
    if(dim == 1) //vertical concatenation
    {
        for(int j=0;j<NbCols;j++)
        {
            for(int i=0;i<A.GetNbRows();i++)
            {
                (*this)[i][j] = A[i][j];
            }
            for(int i=0;i<B.GetNbRows();i++)
            {
                (*this)[i+A.GetNbRows()][j] = B[i][j];
            }
        }
    }
    else
    {
        for(int i=0;i<NbRows;i++)
        {
            for(int j=0;j<A.GetNbCols();j++)
            {
                (*this)[i][j] = A[i][j];
            }
            for(int j=0;j<B.GetNbCols();j++)
            {
                (*this)[i][j+A.GetNbCols()] =   B[i][j];
            }
        }
    }
}


// Destroys all the elements of the array and reshape the dimensions of the array
// Values are initialized to Val (default value = 0)
void M3Matrix::Reshape(int NbRows, int NbCols, double Val)
{
  if( t_NbRows == NbRows && t_NbCols == NbCols )
    return;

  if(t_NbRows * t_NbCols == NbRows * NbCols) //non allocation required
  {
    t_NbRows = NbRows;
    t_NbCols = NbCols;
    //Default: values from the already existing matrix are used
    if(Val != 0.)
        SetValues(Val);
  }
  else
  {
    t_NbRows = NbRows;
    t_NbCols = NbCols;
    // Memory allocation
    Alloc();
    SetValues(Val);
  }
}


void M3Matrix::SetData(int NbRows, int NbCols, double* x)
{
    Reshape(NbRows,NbCols);
    memcpy(Data,x,sizeof(double)*t_NbRows*t_NbCols);
}


void M3Matrix::SetData(double* x)
{
    memcpy(Data,x,sizeof(double)*t_NbRows*t_NbCols);
}

// Memory allocation
void M3Matrix::Alloc()
{
    // Always free before alloc
    Free();
    int size = sizeof(double) * t_NbRows *t_NbCols;
    Data = (double*) malloc(size);
    if( size > 0 )
        assert( Data != nullptr );
}


// Accessors
double M3Matrix::GetVal(int indx)
{
    return Data[indx];
}

int M3Matrix::GetNbRows() const
{
    return t_NbRows;
}

int M3Matrix::GetNbCols() const
{
    return t_NbCols;
}


double* M3Matrix::operator[] (int Row)
{
    return Data + Row*t_NbCols;
}

double* M3Matrix::operator[] (int Row) const
{
    return Data + Row*t_NbCols;
}


// Sets the line i in the M3Matrix Line
void M3Matrix::SetLine(int i, M3Matrix& Line)
{
    Line.SetData(1,t_NbCols,Data+i*t_NbCols);
}


// Sets the col j in the M3Matrix Col
void M3Matrix::SetCol(int j, M3Matrix& Col)
{
    Col.Reshape(t_NbRows,1);
    for(int i=0;i<t_NbRows;i++)
        Col[i][0] = (*this)[i][j];
}


// Operators
void M3Matrix::operator=(const vector<double>& flat)
{
  SetValues(flat);
}


M3Matrix& M3Matrix::operator=(const M3Matrix& Src)
{
  if(t_NbRows != Src.GetNbRows() || t_NbCols != Src.GetNbCols())
    (*this).Reshape(Src.GetNbRows(),Src.GetNbCols());

  for(int i=0;i<t_NbRows;i++)
    for(int j=0;j<t_NbCols;j++)
      (*this)[i][j] = Src.Data[i*t_NbCols+j];

  (*this).InitialiseDynamic(Src);
  return *this;
}


M3Matrix M3Matrix::operator+(const M3Matrix& B)
{
    M3Matrix C(B.t_NbRows,B.t_NbCols);

    if((t_NbRows)&&( t_NbCols))
    {
        for(int i=0;i<t_NbRows;i++)
            for(int j=0;j<t_NbCols;j++)
                C[i][j] = (*this)[i][j] + B[i][j];
    }

    C.InitialiseDynamic(*this);
    return C;
}


vector<double> M3Matrix::operator+( const vector<double>& v1)
{
  vector<double> v2;
  if(GetNbRows() == int(v1.size()))
  {
    v2.resize(v1.size());
    for(int i = 0; i< int(v2.size()); i++)
        v2[i] =  (*this)[i][0] + v1[i];
  }
  else
  {
    cout << "operator addition of vector<double> and M3Matrix used for incompatible sizes" << endl;
    exit(-1);
  }

  return v2;
}



vector<double> M3Matrix::operator-( const vector<double>& v1)
{
  vector<double> v2;
  if(GetNbRows() == int(v1.size()))
  {
    v2.resize(v1.size());
    for(int i = 0; i< int(v2.size()); i++)
        v2[i] =  (*this)[i][0] - v1[i];
  }
  else
  {
    cout << "operator substraction of vector<double> and M3Matrix used for incompatible sizes" << endl;
    exit(-1);
  }

  return v2;
}


//Multiplication of two matrices: use C = A*B
M3Matrix M3Matrix::operator*(const M3Matrix& B)
{
    M3Matrix C(t_NbRows,B.t_NbCols);
    double tmp;
    //Warning: user should ensure that B.t_NbRows = t_NbCols
    if((t_NbRows)&&( t_NbCols)&&(B.t_NbCols))
    {
        for(int i=0;i<t_NbRows;i++)
            for(int j=0;j<B.t_NbCols;j++)
            {
                tmp = 0.;
                for(int k=0;k<t_NbCols;k++)
                    tmp += (*this)[i][k] * B[k][j];

                C[i][j] = tmp;
            }
    }
    return C;
}


//Difference between two matrices: use C = A-B
M3Matrix M3Matrix::operator-(const M3Matrix& B)
{
    M3Matrix C(B.t_NbRows,B.t_NbCols);
    if((t_NbRows)&&( t_NbCols))
    {
        for(int i=0;i<t_NbRows;i++)
            for(int j=0;j<t_NbCols;j++)
                C[i][j] = (*this)[i][j] - B[i][j];
    }

    C.SetNbBits((*this).GetNbBits()+1);
    C.SetRangeMax((*this).GetRangeMax());
    if((*this).GetRangeMin() == 0)
    {
        C.SetRangeMin(-(*this).GetRangeMax()-1);
    }

    return C;
}


vector<double> M3Matrix::operator*(const vector<double>& b)
{
  vector<double> c;
  double sum;
  if(t_NbCols == int(b.size()))
  {
    c.resize(t_NbRows);

    for(int i = 0;i<t_NbRows;i++)
    {
       sum = 0.;
       for(int k = 0; k<t_NbCols; k++)
       {
         sum += (*this)[i][k] * b[k];
       }
       c[i] = sum;
    }
  }
  else
  {
   cout << "M3Matrix::operator*(const double<vector>& b) : you tried to multiply Matrix and vector with incompatible size" << endl;
   exit(1);
  }
  return c;
}


//Addition of a matrix and a scalar: use C = A+t
M3Matrix M3Matrix::operator+(const double& t)
{
    M3Matrix C(t_NbRows,t_NbCols);

    for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols;j++)
            C[i][j] = (*this)[i][j] + t;

    return C;
}

//Addition of a matrix and a scalar: use C = A*t
M3Matrix M3Matrix::operator*(const double& t)
{
    M3Matrix C(t_NbRows,t_NbCols);

    for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols;j++)
            C[i][j] = (*this)[i][j] * t;

    C.InitialiseDynamic(*this);
    return C;
}

//Addition of a matrix and a scalar: use C = A/t
M3Matrix M3Matrix::operator/(const double& t)
{
    M3Matrix C(t_NbRows,t_NbCols);

    for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols;j++)
            C[i][j] = (*this)[i][j] / t;

    C.InitialiseDynamic(*this);
    return C;
}


//Substraction of a scalar from a matrix: use C = A-t
M3Matrix M3Matrix::operator-(const double& t)
{
    M3Matrix C(t_NbRows,t_NbCols);

    for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols;j++)
            C[i][j] = (*this)[i][j] - t;

    C.InitialiseDynamic(*this);
    return C;
}

//Addition of a matrix and a scalar: use C = A+t
void M3Matrix::operator+=(const double& t)
{
    for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols;j++)
        {
            (*this)[i][j] += t;
        }
}


//Division of a matrix by a scalar: use C = A/t
void M3Matrix::operator/=(const double& t)
{
    for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols;j++)
        {
            (*this)[i][j] /= t;
        }
}


//Multiplication of a matrix by a scalar: use C = A*t
void M3Matrix::operator*=(const double& t)
{
    for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols;j++)
        {
            (*this)[i][j] *= t;
        }
}


//Addition of a matrix and a scalar: use C = A+t
void M3Matrix::operator-=(const double& t)
{
    for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols;j++)
        {
            (*this)[i][j] -= t;
        }
}


//Addition of two matrices: use A += B
void M3Matrix::operator+=(M3Matrix& B)
{
    if((t_NbRows)&&( t_NbCols))
    {
        for(int i=0;i<t_NbRows;i++)
            for(int j=0;j<t_NbCols;j++)
                (*this)[i][j] +=   B[i][j];
    }
    return;
}


//Substraction of two matrices: use A -= B
void M3Matrix::operator-=(M3Matrix& B)
{
    if((t_NbRows)&&( t_NbCols))
    {
        for(int i=0;i<t_NbRows;i++)
            for(int j=0;j<t_NbCols;j++)
                (*this)[i][j] -=   B[i][j];
    }
    return;
}


//Multiplication of two matrices: use A *= B
void M3Matrix::operator*=(M3Matrix& B)
{
    M3Matrix C(t_NbRows,B.t_NbCols);
    C=(*this)*B;
    (*this) = C;
    return;
}


M3Matrix M3Matrix::GradientH_Centered()
{
    M3Matrix Extended;
    (*this).Mirror(Extended,1,1);

    M3Matrix Grad((*this).t_NbRows,(*this).t_NbCols);
    for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols-1;j++)
            Grad[i][j] = Extended[i+1][j+2] - Extended[i+1][j];

    return Grad;
}


M3Matrix M3Matrix::GradientV_Centered()
{
    M3Matrix Extended;
    (*this).Mirror(Extended,1,1);

    M3Matrix Grad((*this).t_NbRows,(*this).t_NbCols);
    for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols-1;j++)
            Grad[i][j] = Extended[i+2][j+1] -Extended[i][j+1];

    return Grad;
}


M3Matrix M3Matrix::GradientH()
{
    M3Matrix Grad((*this).t_NbRows,(*this).t_NbCols-1);
    for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols-1;j++)
            Grad[i][j] = (*this)[i][j+1] -(*this)[i][j];

    return Grad;
}


M3Matrix M3Matrix::GradientV()
{
    M3Matrix Grad((*this).t_NbRows-1,(*this).t_NbCols);
    for(int i=0;i<t_NbRows-1;i++)
        for(int j=0;j<t_NbCols;j++)
            Grad[i][j] = (*this)[i+1][j] -(*this)[i][j];

    return Grad;
}


//Pointwise multiplication of two matrices (i.e. C = A.*B in MATLAB)
//Use: C.PointMult(A,B) instead
void M3Matrix::PointMult(M3Matrix& A,M3Matrix& B)
{
    if(this != &A)
        Reshape(B.t_NbRows,B.t_NbCols);

    if((t_NbRows)&&( t_NbCols))
    {
        for(int i=0;i<t_NbRows;i++)
            for(int j=0;j<t_NbCols;j++)
                (*this)[i][j] = A[i][j] * B[i][j];
    }
    return;
}


void M3Matrix::DumpSize()
{
    cout << "DumpSize(): " << (*this).GetNbRows() << " , " << (*this).GetNbCols() << endl;
}


// Dumps the coefficients
void M3Matrix::DumpRGB(ostream& Out, int setw_val)
{
    int NbRows = t_NbRows;
    int NbCols = t_NbCols/3;

    for (int i=0;i<NbRows;i++)
    {
        for (int j=0;j<NbCols;j++)
        {
            Out << setw(setw_val) << (*this)[i][j] << " ";
            Out << setw(setw_val) << (*this)[i][j+NbCols] << " ";
            Out << setw(setw_val) << (*this)[i][j+2*NbCols] << " ";
        }
        Out << endl;
    }
}


// Dumps the matrix coefficients on the output Out
void M3Matrix::Dump(ostream& Out, int setw_val, double tol)
{
    if(tol == 0)
    {
        for (int i=0;i<GetNbRows();i++)
        {
            for (int j=0;j<GetNbCols();j++)
            {
                Out << setw(setw_val) << (*this)[i][j] << " ";
            }
            Out << endl;
        }
    }
    else if (tol > 0)
    {
        for (int i=0;i<GetNbRows();i++)
        {
            for (int j=0;j<GetNbCols();j++)
            {
                double val = (*this)[i][j];
                if(fabs(val)< tol)
                    val = 0.;

                Out << setw(setw_val) <<  val << " ";
            }
            Out << endl;
        }
    }
}


// Dumps the matrix coefficients on the output Out
void M3Matrix::DumpBlock(int rStart, int cStart, int rEnd, int cEnd, ostream& Out, int setw_val, double tol)
{
    if(rStart < 0) rStart = 0;
    if(cStart < 0) cStart = 0;

    if(rEnd > GetNbRows()-1) rStart = GetNbRows()-1;
    if(cEnd > GetNbCols()-1) cStart = GetNbCols()-1;

    if(tol == 0)
    {
        for (int i=rStart;i<=rEnd;i++)
        {
            for (int j=cStart;j<=cEnd;j++)
            {
                Out << setw(setw_val) << (*this)[i][j] << " ";
            }
            Out << endl;
        }
    }
    else if (tol > 0)
    {
        for (int i=rStart;i<=rEnd;i++)
        {
            for (int j=cStart;j<=cEnd;j++)
            {
                double val = (*this)[i][j];
                if(fabs(val)< tol)
                    val = 0.;

                Out << setw(setw_val) <<  val << " ";
            }
            Out << endl;
        }
    }
}



void M3Matrix::DumpRGB_Binary(ostream& Out)
{
    int NbRows = t_NbRows;
    int NbCols = t_NbCols/3;

    if(t_NbBits<=8)
    {
        for (int i=0;i<NbRows;i++)
        {
            for (int j=0;j<NbCols;j++)
            {
                int indx1 = (int)((*this)[i][j]);
                int indx2 = (int)((*this)[i][j+NbCols]);
                int indx3 = (int)((*this)[i][j+2*NbCols]);
                char c1 = (char) indx1;
                char c2 = (char) indx2;
                char c3 = (char) indx3;
                Out << c1;
                Out << c2;
                Out << c3;
            }
        }
    }
    else if(t_NbBits<=16)
    {
        for (int i=0;i<NbRows;i++)
        {
            for (int j=0;j<NbCols;j++)
            {
                int indx1 = (int)((*this)[i][j]);
                int indx2 = (int)((*this)[i][j+NbCols]);
                int indx3 = (int)((*this)[i][j+2*NbCols]);
                int left1 = indx1 >> 8;
                int right1 = indx1 & 255;
                int left2 = indx2 >> 8;
                int right2 = indx2 & 255;
                int left3 = indx3 >> 8;
                int right3 = indx3 & 255;

                char cl1 = (char) left1;
                char cr1 = (char) right1;
                char cl2 = (char) left2;
                char cr2 = (char) right2;
                char cl3 = (char) left3;
                char cr3 = (char) right3;
                Out << cl1 << cr1;
                Out << cl2 << cr2;
                Out << cl3 << cr3;
            }
        }
    }
    else
    {
        cout << "More than 16 bits. This conversion is not supported." << endl;
    }
}


// Dumps the matrix coefficients on the output Out
void M3Matrix::Dump_Binary(ostream& Out)
{
    if(t_NbBits<=8)
    {
        for (int i=0;i<GetNbRows();i++)
        {
            for (int j=0;j<GetNbCols();j++)
            {
                int indx = (int)((*this)[i][j]);
                char c = (char) indx;
                Out << c;
            }
        }
    }
    else if(t_NbBits<=16)
    {
        for (int i=0;i<GetNbRows();i++)
        {
            for (int j=0;j<GetNbCols();j++)
            {
                int indx = (int)((*this)[i][j]);
                int left = indx >> 8;
                int right = indx & 255;

                char c1 = (char) left;
                char c2= (char) right;
                Out << c1 << c2;
            }
        }
    }
    else
    {
        cout << "More than 16 bits. This conversion is not supported." << endl;
    }
}


//This is an out of class function
void SavePPM(string& filenameOut, M3Matrix& R, M3Matrix& G, M3Matrix& B, int NbBits)
{
    M3Matrix RGB;
    if(NbBits == -1) //Default: the number of bits
        NbBits = R.GetNbBits();

    RGB.SetNbBits(NbBits);
    RGB.SetRangeMax(R.GetRangeMax());
    RGB.SetRangeMin(R.GetRangeMin());
    RGB.CombineRGBChannels(R,G,B);
    RGB.SavePPM(filenameOut);
}


//This is an out of class function
void SavePPM_Binary(string& filenameOut, M3Matrix& R, M3Matrix& G, M3Matrix& B, int NbBits)
{
    M3Matrix RGB;
    if(NbBits == -1) //Default: the number of bits
        NbBits = R.GetNbBits();

    RGB.SetNbBits(NbBits);
    RGB.SetRangeMax(R.GetRangeMax());
    RGB.SetRangeMin(R.GetRangeMin());
    RGB.CombineRGBChannels(R,G,B);
    RGB.SavePPM_Binary(filenameOut);
}

// Saves the M3Matrix color values into a PPM file
void M3Matrix::SavePPM_Binary(const string& filenameOut)
{
    M3Matrix A;
    A.Initialise(*this);
    A = (*this);

    ofstream OutStr(filenameOut.c_str());
    A.SavePPM_Binary(OutStr);
    OutStr.close();
}



// Saves the M3Matrix color values into a PPM file
void M3Matrix::SavePPM(const string& filenameOut)
{
    M3Matrix A;
    A.Initialise(*this);
    A = (*this);

    ofstream OutStr(filenameOut.c_str());
    A.SavePPM(OutStr);
    OutStr.close();
}


// Saves the M3Matrix color values into an (ASCII) PPM file
void M3Matrix::SavePPM(ostream& Out)
{
    (*this).Digitalize();

    //Header
    Out << "P3 "<< endl;
    Out << "#CREATOR : SavePPM" << endl;
    Out << t_NbCols/3 << " " << t_NbRows << endl;
    Out << (int)(t_Max+0.5) << endl;

    //Data
    DumpRGB(Out);
}


void M3Matrix::SaveOff(const string& filenameOut, double AmplificationFactor)
{
    ofstream OutStr(filenameOut.c_str());
    (*this).SaveOff(OutStr,AmplificationFactor);
    OutStr.close();
}


void M3Matrix::SaveOff(ostream& Out, double AmplificationFactor)
{
    //Header
    Out << "OFF" << endl;
    int NbTriangles = (t_NbRows-1)*(t_NbCols-1)*2;
    int NbVertices = t_NbRows*t_NbCols;
    int NbEdges = 3*t_NbRows*t_NbCols- 2*t_NbCols -2*t_NbRows+1;

    Out << NbVertices << " ";
    Out << NbTriangles << " ";
    Out << NbEdges << endl;

    //Vertex information
    for(int i = 0;i<t_NbRows;i++)
        for(int j = 0;j<t_NbCols;j++)
        {
            Out << i << " " << j << " " << AmplificationFactor*(*this)[i][j]/t_Max << endl;
        }

    //Face information (here simply triangles)
    for(int i = 0;i<t_NbRows-1;i++)
        for(int j = 0;j<t_NbCols-1;j++)
        {
            //triangle above left to the diagonal
            Out << 3 << " "<< j*t_NbCols + i  << " "  << (j+1)*t_NbCols + i << " " << j*t_NbCols + i+1  << endl;

            //triangle below right to the diagonal
            Out << 3 << " "<< (j+1)*t_NbCols + i+1  << " "  << (j+1)*t_NbCols + i << " " << j*t_NbCols + i+1  << endl;
        }
}


void M3Matrix::SavePGM(const string& filenameOut)
{
    M3Matrix A;
    A.Initialise(*this);
    A = (*this);

    ofstream OutStr(filenameOut.c_str());
    (*this).SavePGM(OutStr);
    OutStr.close();
}


// Saves the M3Matrix values into a PGM file
void M3Matrix::SavePGM(ostream& Out)
{
    (*this).Digitalize();

    //Header
    Out << "P2" << endl;
    Out << "#CREATOR : SavePGM" << endl;
    if(strcmp(m_BayerPattern.c_str(),""))
        Out << "#BP " <<  m_BayerPattern << endl;

    Out << t_NbCols << " " << t_NbRows << endl;
    Out << (int)(t_Max+0.5) << endl;

    //Data
    Dump(Out,3);
}


void M3Matrix::SavePGM_Binary(const string& filenameOut)
{
    M3Matrix A;
    A.Initialise(*this);
    A = (*this);

    ofstream OutStr(filenameOut.c_str());
    (*this).SavePGM_Binary(OutStr);
    OutStr.close();
}


void M3Matrix::SavePGM_Binary(ostream& Out)
{
    Out << "P5" << endl;
    Out << "#CREATOR : SavePGM_Binary" << endl;

    Out << t_NbCols << " " << t_NbRows << endl;
    Out << (int)(t_Max+0.5) << endl;

    Dump_Binary(Out);
}


void M3Matrix::SavePPM_Binary(ostream& Out)
{
    Out << "P6 "; //<< endl;

    Out << t_NbCols/3 << " " << t_NbRows << " "; //<< endl;
    Out << (int)(t_Max+0.5) << endl;

    DumpRGB_Binary(Out);
}


//*****************
// Transformers
//*****************
// Initializes all the M3Matrix coefficients with r_Value
void M3Matrix::Set2Id(int n)
{
    Reshape(n,n,0.);
    for(int i = 0;i<n;i++)
        (*this)[i][i] = 1.0;
}


void M3Matrix::SetValues(double const& Value)
{
    for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols;j++)
            (*this)[i][j] = Value;
}


void M3Matrix::SetValues(vector<double> vals)
{
  if(t_NbRows * t_NbCols == int(vals.size()))
  {
    for(int i=0;i<t_NbRows;i++)
      for(int j=0;j<t_NbCols;j++)
        (*this)[i][j] = vals[i*t_NbCols+j];
  }
  else
  {
    Reshape(vals.size(),1);
    for(int i=0;i<t_NbRows;i++)
        (*this)[i][0] = vals[i];
  }
}


// Exchange columns i and j of the M3Matrix
void M3Matrix::ExchangeCols(int j1, int j2)
{
    double tmp;
    for(int i=0;i<t_NbRows;i++)
    {
        tmp = (*this)[i][j1];
        (*this)[i][j1] = (*this)[i][j2];
        (*this)[i][j2] = tmp;
    }
}


// Exchange rows i and j of the M3Matrix
void M3Matrix::ExchangeRows(int i1, int i2)
{
    double tmp;
    for(int j=0;j<t_NbCols;j++)
    {
        tmp = (*this)[i1][j];
        (*this)[i1][j] = (*this)[i2][j];
        (*this)[i2][j] = tmp;
    }
}


void M3Matrix::Clip()
{
    for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols;j++)
        {
            (*this)[i][j] = myclip((*this)[i][j],t_Min,t_Max);
        }
}

// Computes the minimal coefficient of the array
double M3Matrix::GetMin()
{
    double min = (*this)[0][0];
    for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols;j++)
        {
            if((*this)[i][j] < min)
                min = (*this)[i][j];
        }
    return min;
}


// Computes the maximal coefficient of the array
double M3Matrix::GetMax()
{
    double max = (*this)[0][0];
    for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols;j++)
        {
            if((*this)[i][j] > max)
                max = (*this)[i][j];
        }
    return max;
}


int M3Matrix::GetDynamic()
{
    return  (1<<(GetNbBits()))-1;
}


// Quantization of the coefficients of the M3Matrix according to the QStep
// The Symbol M3Matrix has a size of the original M3Matrix and contains the symbols
// used to be coded
void M3Matrix::Quantize(double& QStep, M3Matrix& Symbol)
{
    Symbol.Reshape(t_NbRows,t_NbCols);

    double tmp,qvalue,symbol;
    double HalfStep = QStep/2;

    for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols;j++)
        {
            tmp = (*this)[i][j] + HalfStep;
            symbol = DBL_ROUND(tmp/QStep);
            if (symbol < 0.) symbol -= 1;

            qvalue  = symbol * QStep;

            (*this)[i][j] = qvalue;
            Symbol[i][j] = symbol;
        }
}


void M3Matrix::Digitalize()
{
    //The following is added for the bit depth
    //First offset the data
    if(t_Min != 0)
    {
        (*this) = (*this) - t_Min;
        t_Max = t_Max - t_Min;
        t_Min = 0;
    }

    //Then factor
    if(static_cast<unsigned long long>(t_Max+0.5) != (1ULL<<static_cast<unsigned long long>(t_NbBits))-1ULL)
    {
        double Step = t_Max/((1ULL<<static_cast<unsigned long long>(t_NbBits))-1);
        double factor = ((double)(1ULL<<static_cast<unsigned long long>(t_NbBits)))/(t_Max+Step);
        Rescale(factor);
    }

    //From now on the original digitalize function
    double value, MaxValue = t_Max;
    for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols;j++)
        {
            value = (int)((*this)[i][j] + 0.5);
            value = myclip(value,0,MaxValue);

            (*this)[i][j] = value;
        }
}


void M3Matrix::InitialiseBitDepth(int NbBitsOut)
{
    SetNbBits(NbBitsOut);
    SetRangeMax((1ULL<<static_cast<unsigned long long>(NbBitsOut))-1ULL);
    SetRangeMin(0);
}


void M3Matrix::SetBitDepth( unsigned int NbBitsOut)
{
    unsigned int NbBits = t_NbBits;
    M3Matrix A = (*this);

    if(NbBitsOut >0 && NbBitsOut <= 64)
    {
        double factor = 0.;
        if(NbBitsOut < NbBits) // cut off
        {
            factor = 1ULL<< static_cast<unsigned long long>( NbBits-NbBitsOut);
            factor *= (t_Max-t_Min)/((1ULL<<static_cast<unsigned long long>(NbBits)) - 1);

            A = (A+t_Min);
            A = A/factor;

            for(int i=0;i<A.GetNbRows();i++)
            {
                for(int j=0;j<A.GetNbCols();j++)
                {
                    A[i][j] = (double)(int)(A[i][j]);
                }
            }

            A = A*factor;
            A = A - t_Min;
        }

        t_NbBits = NbBitsOut;

        double Max = t_Max-t_Min;
        double f = ((1ULL<<static_cast<unsigned long long>(NbBitsOut))-1.0) / ((1ULL<<static_cast<unsigned long long>(NbBits))-1.0);
        Max = Max * f;
        if(NbBitsOut > NbBits)
        {
            t_Max = Max / (1ULL<<static_cast<unsigned long long>(NbBitsOut - NbBits));
        }
        else if ((NbBits > NbBitsOut) )
        {
            t_Max =Max*(1ULL<<static_cast<unsigned long long>(NbBits-NbBitsOut));
        }
        else
        {
            t_Max = Max;
        }

        t_Max += t_Min;

        for(int i=0;i<A.GetNbRows();i++)
            for(int j=0;j<A.GetNbCols();j++)
            {
                (*this)[i][j] = A[i][j];
            }
    }
}


void M3Matrix::Rescale(double factor)
{
    (*this) *= factor;
    t_Max = t_Max*factor;
    t_Min = t_Min*factor;
}


// Function to sort the coefficients of the M3Matrix by the i column
void M3Matrix::SortByCol(int c)
{
    // The column i is set as first col
    if(c != 0)
        ExchangeCols(c,0);

    qsort(Data,t_NbRows,sizeof(double *),fncompare);

    // Column newly exchanged to provide a good
    if(c != 0)
        ExchangeCols(c,0);
}


// Function to sort the coefficients of the M3Matrix by the i
// and then the j columns
//Warning : c1 MUST be different from c2
void M3Matrix::SortByCols(int c1, int c2)
{
    // The column i is set as first col
    if(c1 != 0)
        ExchangeCols(c1,0);
    if(c2 != 1)
        ExchangeCols(c2,1);

    qsort(Data,t_NbRows,sizeof(double *),fncompare2);

    // Column newly exchanged to provide a good
    if(c2 != 1)
        ExchangeCols(c2,1);
    if(c1 != 0)
        ExchangeCols(c1,0);
}


// Function to sort the coefficients of the M3Matrix according to blocs
void M3Matrix::SortByBlocs(int c1, int c2)
{
    // The column i is set as first col
    if(c1 != 0)
        ExchangeCols(c1,0);
    if(c2 != 1)
        ExchangeCols(c2,1);

    qsort(Data,t_NbRows,sizeof(double *),fncomparebloc);

    if(c2 != 1)
        ExchangeCols(c2,1);
    if(c1 != 0)
        ExchangeCols(c1,0);
}


// Private function used for sorting purpose
int M3Matrix::vCompareFirstCol(const void *a1, const void *a2)
{
    int res;
    int i1 = *(int *) a1;
    int i2 = *(int *) a2;

    if((*this)[i1][0] > (*this)[i2][0])
        res = -1;
    else
    {
        if((*this)[i1][0] < (*this)[i2][0])
            res = 1;
        else
            res = 0;
    }
    return res;
}


double M3Matrix::GetSum()
{
    double Sum_val = 0.;

    for(int i = 0;i<t_NbRows;i++)
        for(int j = 0;j<t_NbCols;j++)
        {
            Sum_val += (*this)[i][j];
        }

    return Sum_val;
}


// Computes the Lp norm or quasi-norm of a matrix
// typically makes sense for  wavelet coefficients or similar
// p must be strictly positive !!
double M3Matrix::LpNorm(double p)
{
    double tmp = 0.;
    double threshold = 0.0000000001;

    for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols;j++)
        {
            if(fabs((*this)[i][j]) > threshold)
                tmp += pow((fabs((*this)[i][j])),p);
        }

    tmp = pow(tmp,1.0/p);

    return tmp;
}


// Computes the number of non zero values: useful for measuring the sparsity
double M3Matrix::L0Norm()
{
    double tmp = 0.;
    for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols;j++)
        {
            if(fabs((*this)[i][j]) > 0.0)
                tmp += 1.0;
        }
    return tmp;
}



double M3Matrix::ScaledMeanKurtosis(int Margin)
{
    double GradHKurtosis = (*this).GradientH().Kurtosis(Margin);
    double GradVKurtosis = (*this).GradientV().Kurtosis(Margin);

    double ScaledMeanKurtosis = (GradHKurtosis + GradVKurtosis)*0.5;
    ScaledMeanKurtosis *= (*this).GetMax() - (*this).GetMin();

    return ScaledMeanKurtosis;
}


double M3Matrix::MeanKurtosisOfHGradient(int Margin)
{
    double GradHKurtosis = (*this).GradientH().Kurtosis(Margin);
    return GradHKurtosis;
}


double M3Matrix::MeanKurtosisOfVGradient(int Margin)
{
    double GradVKurtosis = (*this).GradientV().Kurtosis(Margin);
    return GradVKurtosis;
}


double M3Matrix::MeanKurtosisOfGradients(int Margin)
{
    double GradHKurtosis = (*this).GradientH().Kurtosis(Margin);
    double GradVKurtosis = (*this).GradientV().Kurtosis(Margin);

    double MeanKurtosis = (GradHKurtosis + GradVKurtosis)*0.5;

    return MeanKurtosis;
}


//Computes the statistical kurtosis of the matrix values (without the pixels in a margin of size Margin)
double M3Matrix::Kurtosis(int Margin)
{
    double kurtosis = 0, mu4 = 0., sigma = 0.;

    double mean = GetSum()/(t_NbRows*t_NbCols);
    M3Matrix Centered;
    Centered = (*this) - mean;
    for(int i = Margin;i<t_NbRows-Margin;i++)
        for(int j = Margin;j<t_NbCols-Margin;j++)
        {
            mu4  += pow(fabs(Centered[i][j]),4);
            sigma += Centered[i][j]*Centered[i][j];
        }

    mu4 /= ((t_NbRows-2*Margin)*(t_NbCols-2*Margin));
    sigma /= ((t_NbRows-2*Margin)*(t_NbCols-2*Margin));

    if(sigma != 0)
        kurtosis = mu4/(sigma*sigma);
    else
        kurtosis = 0;

    return kurtosis;
}


//Warning: assumes positive entries
void M3Matrix::Normalise()
{
    double Sum_val = (*this).GetSum();
    if(Sum_val > 0)
    {
        (*this) /= Sum_val;
    }
}


//Compute (*this) = A-B, and rescale the image to be visualized between 0 and Dynamic
//WARNING: this works for images whose values are in the range [0,Dynamic]
void M3Matrix::DifferenceImage(M3Matrix& A, M3Matrix& B, double Factor, double Dynamic)
{
    (*this) = (A - B) * Factor + Dynamic;
    (*this) = (*this) * 0.5;
    (*this).Digitalize();
}


//*******************************************************
//Bilateral filter with median (instead of box filter)
//*******************************************************
void M3Matrix::AdaptiveBilateralFilter_MedMed(M3Matrix& FilteredImage, int Size, double Threshold, double Dynamic)
{
    Dynamic = t_Max; // Check this

    M3Matrix MirroredImage,  X;
    double /*FilteredValue,*/Diff;
    int HalfHeight = Size/2,HalfWidth, Count;
    if(t_NbRows > 1)
        HalfWidth = Size/2;
    else
        HalfWidth = 0;

    (*this).Mirror(MirroredImage,HalfHeight,HalfWidth);
    FilteredImage.Reshape(t_NbRows,t_NbCols);

    for(int i = 0;i<t_NbRows;i++)
    {
        for(int j = 0;j<t_NbCols;j++)
        {
            Count = 0;
            X.Reshape(Size,Size);

            for(int n=0;n<Size;n++)
            {
                Count = 0;
                for(int m=0;m<Size;m++)
                {
                    Diff = MirroredImage[i+m][j+n] - (*this)[i][j];
                    if(fabs(Diff) < Threshold*Dynamic)
                    {
                        X[m][n] = MirroredImage[i+m][j+n];
                    }
                    else
                    {
                        Count++;
                        if((Count+n)%2==0)
                        {
                            X[m][n] = Dynamic;
                        }
                        else
                        {
                            X[m][n] = 0;
                        }
                    }
                }
            }

            FilteredImage[i][j] = X.MedMed();
        }
    }
}


//**********************************+
// Max Filter
//**********************************+
M3Matrix M3Matrix::MaxFilter(int rSize, int cSize)
{
  M3Matrix MirroredImage, FilteredImage;
  int HalfHeight = rSize/2, HalfWidth = cSize/2;

  FilteredImage.Reshape(t_NbRows,t_NbCols);
  FilteredImage.SetValues(0.);
  (*this).Mirror(MirroredImage,HalfHeight,HalfWidth);

  for(int i = 0;i<(*this).GetNbRows();i++)
    for(int j = 0;j<(*this).GetNbCols();j++)
    {
      double MaxVal = MirroredImage[i][j];
      for(int m=0;m<rSize;m++)
      {
        for(int n=0;n<cSize;n++)
        {
          if(MirroredImage[i+m][j+n] > MaxVal)
            MaxVal  = MirroredImage[i+m][j+n];
        }
      }
      FilteredImage[i][j] = MaxVal;
  }
  return FilteredImage;
}


void M3Matrix::FilterMax(int rSize, int cSize)
{
    M3Matrix MirroredImage, FilteredImage;
    int HalfHeight = rSize/2, HalfWidth = cSize/2;

    FilteredImage.Reshape(t_NbRows,t_NbCols);
    FilteredImage.SetValues(0.);
    (*this).Mirror(MirroredImage,HalfHeight,HalfWidth);

    for(int i = 0;i<(*this).GetNbRows();i++)
        for(int j = 0;j<(*this).GetNbCols();j++)
        {
            double MaxVal = MirroredImage[i][j];
            for(int m=0;m<rSize;m++)
            {
                for(int n=0;n<cSize;n++)
                {
                    if(MirroredImage[i+m][j+n] > MaxVal)
                        MaxVal  = MirroredImage[i+m][j+n];
                }
            }
            FilteredImage[i][j] = MaxVal;
        }

    for(int i = 0;i<(*this).GetNbRows();i++)
        for(int j = 0;j<(*this).GetNbCols();j++)
            (*this)[i][j] = FilteredImage[i][j];
}


void M3Matrix::AdaptiveBilateralFilter_Median(M3Matrix& FilteredImage, int Size, double Threshold, double Dynamic)
{
    Dynamic = t_Max;

    M3Matrix MirroredImage, AdaptiveMask, X;
    double Diff;
    int HalfHeight = Size/2,HalfWidth, Count;
    if(t_NbRows > 1)
        HalfWidth = Size/2;
    else
        HalfWidth = 0;

    (*this).Mirror(MirroredImage,HalfHeight,HalfWidth);
    FilteredImage.Reshape(t_NbRows,t_NbCols);

    if(GetNbRows() > 1)
        AdaptiveMask.Reshape(Size,Size);
    else
        AdaptiveMask.Reshape(1,Size);

    for(int i = 0;i<t_NbRows;i++)
    {
        for(int j = 0;j<t_NbCols;j++)
        {
            AdaptiveMask.SetValues(0.);
            Count = 0;

            for(int m=0;m<AdaptiveMask.GetNbRows();m++)
            {
                for(int n=0;n<AdaptiveMask.GetNbCols();n++)
                {
                    Diff = MirroredImage[i+m][j+n] - (*this)[i][j];
                    if(fabs(Diff) < Threshold*Dynamic)
                    {
                        AdaptiveMask[m][n] = 1.;
                        Count++;
                    }
                }
            }

            X.Reshape(Count,1);
            Count = 0;

            for(int m=0;m<AdaptiveMask.GetNbRows();m++)
            {
                for(int n=0;n<AdaptiveMask.GetNbCols();n++)
                {
                    if(AdaptiveMask[m][n] == 1)
                    {
                        X[Count][0] = MirroredImage[i+m][j+n];
                        Count++;
                    }
                }
            }

            FilteredImage[i][j] = X.Median();
        }
    }
}


// 3D bilateral filter
void M3Matrix::BilateralTemporalFilter(M3Matrix& FilteredImage, M3Matrix3D& Mask3D, M3Matrix3D& Buffer, int CircIndx,  BilParams Parameters)
{
    FilteredImage.Reshape(t_NbRows,t_NbCols);
    FilteredImage.InitialiseDynamic(*this);

    int MaxPastFrames = Buffer.GetNbFrames();
    int LocalIndx;

    //Copy from 2D-bilateral
    BilParams ParametersApp;
    ParametersApp = Parameters; //copies a,b, and Dynamic into ParametersApp
    ParametersApp.setAWType(4);

    int HalfHeight = Mask3D[0].GetNbRows()/2;
    int HalfWidth = Mask3D[0].GetNbCols()/2;

    M3Matrix3D MirroredImageStack(MaxPastFrames+1);
    (*this).Mirror(MirroredImageStack[0],HalfHeight,HalfWidth); // mirroring of the current image
    //mirroring of the past images
    for(int f = 0;f<MaxPastFrames;f++)
    {
        LocalIndx = (CircIndx-f)%MaxPastFrames;
        Buffer[LocalIndx].Mirror(MirroredImageStack[f+1],HalfHeight,HalfWidth);
    }

    M3Matrix3D AdaptiveMask3D(MaxPastFrames+1,Mask3D.GetNbRows(),Mask3D.GetNbCols());

    double FilteredDiff,Diff,AdaptiveWeightFactor,NormalisationFactorApprox;
    for(int i = 0;i<t_NbRows;i++)
    {
        for(int j = 0;j<t_NbCols;j++)
        {
            //Compute first the adaptive weights
            FilteredDiff = 0.;
            NormalisationFactorApprox = 0.;

            double CentralPixelValue = MirroredImageStack[0][i+HalfHeight][j+HalfWidth];

            for(int f= 0;f<=MaxPastFrames;f++)
            {
                for(int m=0;m<Mask3D.GetNbRows();m++)
                {
                    for(int n=0;n<Mask3D.GetNbCols();n++)
                    {
                        Diff = MirroredImageStack[f][i+m][j+n] - CentralPixelValue;
                        AdaptiveWeightFactor = ComputeAdaptiveWeight(fabs(Diff),Parameters); //phi(\delta_{m,n})
                        NormalisationFactorApprox +=  ComputeAdaptiveWeight(fabs(Diff),ParametersApp) * Mask3D[f][m][n]; //phi(\delta_{m,n}) w_{m,n}
                        AdaptiveMask3D[f][m][n] = AdaptiveWeightFactor * Mask3D[f][m][n]; //
                    }
                }
            }


            AdaptiveMask3D.Normalise();

            //Perform the actual filtering
            for(int f= 0;f<=MaxPastFrames;f++)
            {
                for(int m=0;m<Mask3D.GetNbRows();m++)
                {
                    for(int n=0;n<Mask3D.GetNbCols();n++)
                    {
                        Diff = MirroredImageStack[f][i+m][j+n] - CentralPixelValue;
                        FilteredDiff += AdaptiveMask3D[f][m][n]*Diff;
                    }
                }
            }

            double Value = CentralPixelValue + FilteredDiff;
            FilteredImage[i][j]  = myclip(Value,t_Min,t_Max);
        }
    }
}


// Normalised version of the bilateral filter: version with bilateral ratio (Added July 2015)
void M3Matrix::AdaptiveBilateralFilter_BilRatio(M3Matrix& FilteredImage, M3Matrix& Mask, BilParams Parameters)
{
    cout << "AdaptiveBilateralFilter_BilRatio : " << endl;

    BilParams ParametersApp;
    ParametersApp = Parameters; //copies a,b, and Dynamic into ParametersApp
    ParametersApp.setAWType(4);

    M3Matrix MirroredImage;
    int HalfHeight = Mask.GetNbRows()/2;
    int HalfWidth = Mask.GetNbCols()/2;

    (*this).Mirror(MirroredImage,HalfHeight,HalfWidth);
    FilteredImage.Reshape(t_NbRows,t_NbCols);
    FilteredImage.InitialiseDynamic(*this);

    M3Matrix AdaptiveMask(Mask.GetNbRows(),Mask.GetNbCols());

    double FilteredValue,Diff,AdaptiveWeightFactor,NormalisationFactorApprox;
    for(int i = 0;i<t_NbRows;i++)
        for(int j = 0;j<t_NbCols;j++)
        {
            //Compute first the adaptive weights
            FilteredValue = 0.;
            NormalisationFactorApprox = 0.;

            double SumPhi = 0.;

            for(int m=0;m<Mask.GetNbRows();m++)
            {
                for(int n=0;n<Mask.GetNbCols();n++)
                {
                    Diff = MirroredImage[i+m][j+n] - MirroredImage[i+HalfHeight][j+HalfWidth];
                    AdaptiveWeightFactor = ComputeAdaptiveWeight(fabs(Diff),Parameters); //phi(\delta_{m,n})
                    SumPhi += AdaptiveWeightFactor;
                    NormalisationFactorApprox +=  ComputeAdaptiveWeight(fabs(Diff),ParametersApp) * Mask[m][n]; //phi(\delta_{m,n}) w_{m,n}
                    AdaptiveMask[m][n] = AdaptiveWeightFactor * Mask[m][n]; //
                }
            }

            //Here insert the computation of the bilateral ratio
            double r = 1 - SumPhi / ((double)(Mask.GetNbRows()*Mask.GetNbCols())); //This is the "Friedberg" formula

            //Different types of normalisations
            if(Parameters.norm_type == BilParams::NormType::off) // off (because "Mask" is normalized already, undo it)
            {
                AdaptiveMask *= Mask.GetNbRows() * Mask.GetNbCols();
            }
            else if(Parameters.norm_type == BilParams::NormType::constant) //divide by sum w_i - usually this should be already equal to 1
            {
                AdaptiveMask /= Mask.GetSum();
            }
            else if(Parameters.norm_type == BilParams::NormType::exact) //corresponds to norm_exact, the exact continuous normalisation
            {
                AdaptiveMask.Normalise();
            }
            else if(Parameters.norm_type == BilParams::NormType::approx) //corresponds to norm_approx
            {
                double dummy = 1./NormalisationFactorApprox;
                AdaptiveMask = AdaptiveMask*dummy;
            }
            //Perform the actual filtering
            double t1,t2;
            t1 = Parameters.t1_br;
            t2 = Parameters.t2_br;

            double LinearFilteredValue = 0.;
            double BilateralFilteredValue = 0.;

            for(int m=0;m<Mask.GetNbRows();m++)
            {
                for(int n=0;n<Mask.GetNbCols();n++)
                {
                    Diff = MirroredImage[i+m][j+n] - MirroredImage[i+HalfHeight][j+HalfWidth];
                    BilateralFilteredValue += AdaptiveMask[m][n]*Diff;
                    LinearFilteredValue += Mask[m][n]*Diff;
                }
            }

            //Now FilteredValue is computed, according to the bilateral ratio
            if(r>=t2)  //This is a texture
            {
                FilteredValue = LinearFilteredValue;
            }
            else if(t1 <r && r < t2)
            {
                //Hybrid modus
                FilteredValue = (r-t1)*LinearFilteredValue/(t2-t1) + (t2-r)*BilateralFilteredValue/(t2-t1);
            }
            else //more or less edges and flat areas ?
            {
                FilteredValue = BilateralFilteredValue;
            }

            double Value = MirroredImage[i+HalfHeight][j+HalfWidth] + FilteredValue;
            FilteredImage[i][j]  = myclip(Value,t_Min,t_Max);
        }
}


// Normalised version of the bilateral filter
void M3Matrix::AdaptiveBilateralFilter_EdgeMasked(M3Matrix& FilteredImage, M3Matrix& Mask, M3Matrix& EdgeMask, BilParams Parameters)
{
    BilParams ParametersApp;
    ParametersApp = Parameters; //copies a,b, and Dynamic into ParametersApp
    ParametersApp.setAWType(4);

    M3Matrix MirroredImage, Mirrored_EdgeMask;
    int HalfHeight = Mask.GetNbRows()/2;
    int HalfWidth = Mask.GetNbCols()/2;

    (*this).Mirror(MirroredImage,HalfHeight,HalfWidth);
    EdgeMask.Mirror(Mirrored_EdgeMask,HalfHeight,HalfWidth);

    FilteredImage.Reshape(t_NbRows,t_NbCols);
    FilteredImage.InitialiseDynamic(*this);

    M3Matrix AdaptiveMask(Mask.GetNbRows(),Mask.GetNbCols());

    double FilteredValue,Diff,AdaptiveWeightFactor,NormalisationFactorApprox;
    for(int i = 0;i<t_NbRows;i++)
        for(int j = 0;j<t_NbCols;j++)
        {
            //Compute first the adaptive weights
            FilteredValue = 0.;
            NormalisationFactorApprox = 0.;
            AdaptiveMask.SetValues(0.);

            if(EdgeMask[i][j] == 0)
            {
                double SumPhi = 0.;

                for(int m=0;m<Mask.GetNbRows();m++)
                {
                    for(int n=0;n<Mask.GetNbCols();n++)
                    {
                        if(Mirrored_EdgeMask[i+m][j+n] == 0) //the filter is performed only on non edge neighbour pixels
                        {
                            Diff = MirroredImage[i+m][j+n] - MirroredImage[i+HalfHeight][j+HalfWidth];
                            AdaptiveWeightFactor = ComputeAdaptiveWeight(fabs(Diff),Parameters); //phi(\delta_{m,n})
                            SumPhi += AdaptiveWeightFactor;
                            NormalisationFactorApprox +=  ComputeAdaptiveWeight(fabs(Diff),ParametersApp) * Mask[m][n]; //phi(\delta_{m,n}) w_{m,n}
                            AdaptiveMask[m][n] = AdaptiveWeightFactor * Mask[m][n]; //
                        }
                    }
                }

                //Here insert the computation of the bilateral ratio
                double r = 1 - SumPhi / ((double)(Mask.GetNbRows()*Mask.GetNbCols())); //This is the "Friedberg" formula
                r = r;

                //Different types of normalisations
                if(Parameters.norm_type == BilParams::NormType::off) // off (because "Mask" is normalized already, undo it): should not be used in the usual cases
                {
                    AdaptiveMask *= Mask.GetNbRows() * Mask.GetNbCols();
                }
                else if(Parameters.norm_type == BilParams::NormType::constant) //divide by sum w_i - usually this should be already equal to 1
                {
                    AdaptiveMask /= Mask.GetSum();
                }
                else if(Parameters.norm_type == BilParams::NormType::exact) //corresponds to norm_exact, the exact continuous normalisation
                {
                    AdaptiveMask.Normalise();
                }
                else if(Parameters.norm_type == BilParams::NormType::approx) //corresponds to norm_approx
                {
                    double dummy = 1./NormalisationFactorApprox;
                    AdaptiveMask = AdaptiveMask*dummy;
                }

                //Perform the actual filtering
                for(int m=0;m<Mask.GetNbRows();m++)
                {
                    for(int n=0;n<Mask.GetNbCols();n++)
                    {
                        Diff = MirroredImage[i+m][j+n] - MirroredImage[i+HalfHeight][j+HalfWidth];
                        FilteredValue += AdaptiveMask[m][n]*Diff;
                    }
                }
            }
            double Value = MirroredImage[i+HalfHeight][j+HalfWidth] + FilteredValue;
            FilteredImage[i][j]  = myclip(Value,t_Min,t_Max);

            if(EdgeMask[i][j] != 0 )//edge
            {
                cout << "edge" << endl;
                if(FilteredImage[i][j] != (*this)[i][j])
                {
                    cout << "strange " << endl;
                }
            }
        }
}


// Normalised version of the bilateral filter
void M3Matrix::AdaptiveBilateralFilter(M3Matrix& FilteredImage, M3Matrix& Mask, BilParams Parameters)
{
    if(Parameters.norm_type == BilParams::NormType::fpga) //this should be like in vhdl !
    {
        cout << "AdaptiveBilateral: vhdl  division modus recognized" << endl;
    }

    int ri=0,rj=0,bi=0,bj=0;
    if(Parameters.GreenOnly == true)
    {
        cout << "AdaptiveBilateralFilter :: GreenOnly recognized" << endl;
        int BayerPattern = Parameters.BPattern;
        if(BayerPattern == 0)
        {
            ri = 0; rj = 0;
        }
        else if(BayerPattern == 1)
        {
            ri = 0; rj = 1;
        }
        else if(BayerPattern == 2)
        {
            ri = 1; rj = 0;
        }
        else if(BayerPattern == 3)
        {
            ri = 1;
            rj = 1;
        }

        bi = 1-ri;
        bj = 1-rj;
    }


    BilParams ParametersApp;
    ParametersApp = Parameters; //copies a,b, and Dynamic into ParametersApp
    ParametersApp.setAWType(4);

    M3Matrix MirroredImage;
    int HalfHeight = Mask.GetNbRows()/2;
    int HalfWidth = Mask.GetNbCols()/2;

    (*this).Mirror(MirroredImage,HalfHeight,HalfWidth);
    FilteredImage.Reshape(t_NbRows,t_NbCols);
    FilteredImage.InitialiseDynamic(*this);

    M3Matrix AdaptiveMask(Mask.GetNbRows(),Mask.GetNbCols());
    int add_bits = Parameters.add_bits;

    double FilteredValue,Diff,AdaptiveWeightFactor,NormalisationFactorApprox;
    for(int i = 0;i<t_NbRows;i++)
        for(int j = 0;j<t_NbCols;j++)
        {
            if( Parameters.GreenOnly == true && (  ( i%2 == ri && j%2 == rj) || (i%2 == bi && j%2 == bj) )  )
            {
                FilteredImage[i][j] = (*this)[i][j];
                continue;
            }

            //Compute first the adaptive weights
            FilteredValue = 0.;
            NormalisationFactorApprox = 0.;

            double SumPhi = 0.;

            for(int m=0;m<Mask.GetNbRows();m++)
            {
                for(int n=0;n<Mask.GetNbCols();n++)
                {
                    Diff = MirroredImage[i+m][j+n] - MirroredImage[i+HalfHeight][j+HalfWidth];
                    AdaptiveWeightFactor = ComputeAdaptiveWeight(fabs(Diff),Parameters); //phi(\delta_{m,n})
                    SumPhi += AdaptiveWeightFactor;
                    NormalisationFactorApprox +=  ComputeAdaptiveWeight(fabs(Diff),ParametersApp) * Mask[m][n]; //phi(\delta_{m,n}) w_{m,n}
                    AdaptiveMask[m][n] = AdaptiveWeightFactor * Mask[m][n]; //
                }
            }

            //Here insert the computation of the bilateral ratio
            double r = 1 - SumPhi / ((double)(Mask.GetNbRows()*Mask.GetNbCols())); //This is the "Friedberg" formula

            //Different types of normalisations
            if(Parameters.norm_type == BilParams::NormType::off) // off (because "Mask" is normalized already, undo it): should not be used in the usual cases
            {
                AdaptiveMask *= Mask.GetNbRows() * Mask.GetNbCols();
            }
            else if(Parameters.norm_type == BilParams::NormType::constant) //divide by sum w_i - usually this should be already equal to 1
            {
                AdaptiveMask /= Mask.GetSum();
            }
            else if(Parameters.norm_type == BilParams::NormType::exact) //corresponds to norm_exact, the exact continuous normalisation
            {
                AdaptiveMask.Normalise();
            }
            else if(Parameters.norm_type == BilParams::NormType::approx) //corresponds to norm_approx
            {
                double dummy = 1./NormalisationFactorApprox;
                AdaptiveMask = AdaptiveMask*dummy;
            }

            //Perform the actual filtering
            if(Parameters.norm_type != BilParams::NormType::fpga)
            {
                for(int m=0;m<Mask.GetNbRows();m++)
                {
                    for(int n=0;n<Mask.GetNbCols();n++)
                    {
                        Diff = MirroredImage[i+m][j+n] - MirroredImage[i+HalfHeight][j+HalfWidth];
                        FilteredValue += AdaptiveMask[m][n]*Diff; // \psi(\delta_{m,n})
                    }
                }
            }
            else if(Parameters.norm_type == BilParams::NormType::fpga) //this should be like in vhdl !
            {
                int NumberOfPixels = Mask.GetNbRows() * Mask.GetNbCols();
                for(int m=0;m<Mask.GetNbRows();m++)
                {
                    for(int n=0;n<Mask.GetNbCols();n++)
                    {
                        Diff = MirroredImage[i+m][j+n] - MirroredImage[i+HalfHeight][j+HalfWidth];

                        double psi_val = ComputePsi(Diff, Parameters); //Direct psi computation, like in vhdl
                        FilteredValue += psi_val;
                    }
                }

                //Normalisation: precision to be fixed !!!
                int div_bits = 20;

                uint64_t Normalisation_inv = (uint64_t) round((1.0 / (double) NumberOfPixels) * (1 << div_bits));

                int64_t IntFilteredValue = FilteredValue;

                IntFilteredValue  *= Normalisation_inv;
                IntFilteredValue >>= div_bits - add_bits - 1;
                IntFilteredValue  += 1;
                IntFilteredValue >>= 1;

                FilteredValue = IntFilteredValue;
            }

            int tmp_val = (int) MirroredImage[i+HalfHeight][j+HalfWidth];
            double Value = (tmp_val << add_bits) + FilteredValue;
            FilteredImage[i][j]  = myclip(Value,(int)t_Min << add_bits, (int)t_Max << add_bits);
        }
}


double ComputeAdaptiveWeight(double AbsDiff,BilParams& Parameters)
{
    double Dynamic = Parameters.Dynamic;
    double AdaptiveWeightFactor;
    double x = AbsDiff/Dynamic;

    int Option = Parameters.getAWType();
    double Threshold = Parameters.Threshold;
    Option = 3;

    if(Option == 0)
    {
        //Linear Modus
        AdaptiveWeightFactor = 1-x/Threshold;

        if(AdaptiveWeightFactor<0)
            AdaptiveWeightFactor = 0;
    }
    else if(Option == 1) //
    {
        if(Threshold > 0)
            x=x/Threshold;

        if(x<0.25)
            AdaptiveWeightFactor = 1;
        else if(x < 0.75)
            AdaptiveWeightFactor =1/(4*x);
        else if(x<1)
            AdaptiveWeightFactor =1/x-1;
        else
            AdaptiveWeightFactor = 0;
    }
    else if(Option == 2) //step function (cf. original FPGA version)
    {
        if(x < Threshold/2.0)
            AdaptiveWeightFactor = 1;
        else if(x < Threshold)
            AdaptiveWeightFactor = 0.5;
        else
            AdaptiveWeightFactor = 0.;
    }
    else if(Option == 3) //This is  a continuous version (feb. 2015) thought for a suitable fpga implementation but still with full accuracy
        //in particular the according phi functions have larger tails. cf. test_fpga_approximative_weights.m
    {
        //Implements the function phi = min(1,a+b/x);
        if(x == 0)
        {
            AdaptiveWeightFactor = 1;
        }
        else
        {
            double a = Parameters.a;  //standard value:  a = -8./255.;
            double b = Parameters.b;  //standard value : double b = 12./255.;

            AdaptiveWeightFactor = a + b / x;

            if(AdaptiveWeightFactor >= 1)
                AdaptiveWeightFactor  = 1.;
            if(AdaptiveWeightFactor < 0)
                AdaptiveWeightFactor  = 0.;
        }
    }
    else if(Option == 4) //fpga approx. (april 2015)
    {
        //Implements the piewise linear approximation  of the function phi = min(1,a+b/x);
        if(x == 0)
        {
            AdaptiveWeightFactor = 1;
        }
        else
        {
            double a = Parameters.a;
            double b = Parameters.b;
            double t = 2./3.;
            double x0 = -b/(a-1); //BIG WARNING: no test is done so far. To be corrected if this bugs !!
            double x2 = -b/(a); //BIG WARNING: no test is done so far. To be corrected if this bugs !!
            double y0 = 1;
            double y2 = 0;
            double x1 = t *x0 +(1-t)*x2;
            double y1 = a+b/x1;

            if(x == 0)
                AdaptiveWeightFactor = 1;
            else
            {
                if(x<x0)
                    AdaptiveWeightFactor = 1;
                else if(x>=x2)
                    AdaptiveWeightFactor = 0;
                else if (x<x1)
                {
                    AdaptiveWeightFactor=y0+(y1-y0)*(x-x0)/(x1-x0);
                }
                else
                {
                    AdaptiveWeightFactor=y1+(y2-y1)*(x-x1)/(x2-x1);
                }
            }
        }
    }
    else
    {
        AdaptiveWeightFactor = 1;
    }

    return AdaptiveWeightFactor;
}


// Diff: 0 .. 2^bits - 1
double ComputePsi(double Diff, BilParams& Parameters)
{
    int64_t PsiValue;
    int64_t IntDiff = Diff;
    uint64_t AbsDiff = abs(IntDiff);

    //Implements the function psi = min(x,ax+b);
    int nbBits = Parameters.nbBits;
    int64_t a = round(Parameters.a * (1ULL << static_cast<unsigned long long>(nbBits))); // 0 .. -2^bits - 1
    int64_t b = round(Parameters.b * (1ULL << static_cast<unsigned long long>((2 * nbBits)))); // 0 .. 2^(2 * bits) - 1 (normalized to multiplication result bit count)

    uint64_t x0 = round((double)b / (double)((1ULL << static_cast<unsigned long long>(nbBits)) - a));
    uint64_t x1 = round((double)b / (double)(-a));

    if(a>=0 )
    {
        exit(1);   //todo exception better ?
        //throw string e("a parameter is positive or zero: not allowed"));
    }

    if (AbsDiff < x0)
    {
        PsiValue = Diff;
    }
    else if (AbsDiff >= x1)
    {
        PsiValue = 0;
    }
    else
    {
        // VV shouldn't this be  (1ULL << static_cast<unsigned long long>(nbBits))- 1)
        PsiValue   = a * AbsDiff + (b + (1ULL << static_cast<unsigned long long>((nbBits - 1))));             // calculate psi(delta) with rounding preparation
        PsiValue   = (IntDiff >= 0) ? PsiValue : -PsiValue;               // differenciate between positive and negative difference
        PsiValue >>= nbBits ;                                             // normalize and round
        PsiValue   = myclip(PsiValue, -(1ULL << static_cast<unsigned long long>(nbBits)), (1ULL << static_cast<unsigned long long>(nbBits)) - 1); // clipping
    }
    return PsiValue; // 0 .. 2^bit - 1
}


double ComputeAdaptiveWeight(double AbsDiff, int Option, double Threshold, double Dynamic)
{
    double AdaptiveWeightFactor;
    double x = AbsDiff/Dynamic;

    if(Option == 0)
    {
        //Linear Modus
        AdaptiveWeightFactor = 1-x/Threshold;

        if(AdaptiveWeightFactor<0)
            AdaptiveWeightFactor = 0;
    }
    else if(Option == 1) //
    {
        if(Threshold > 0)
            x=x/Threshold;

        if(x<0.25)
            AdaptiveWeightFactor = 1;
        else if(x < 0.75)
            AdaptiveWeightFactor =1/(4*x);
        else if(x<1)
            AdaptiveWeightFactor =1/x-1;
        else
            AdaptiveWeightFactor = 0;
    }
    else if(Option == 2) //step function (cf. original FPGA version)
    {
        if(x < Threshold/2.0)
            AdaptiveWeightFactor = 1;
        else if(x < Threshold)
            AdaptiveWeightFactor = 0.5;
        else
            AdaptiveWeightFactor = 0.;
    }
    else if(Option == 3) //FPGA version: February 2015. Contains a larger tail. cf. test_fpga_approximative_weights.m
    {
        //Implements the function phi = min(1,a+b/x);
        if(x == 0)
        {
            AdaptiveWeightFactor = 1;
        }
        else
        {
            double a = -8./255.; //a multiplication will occur.
            double b = 12./255.;

            AdaptiveWeightFactor = a + b / x;
            if(AdaptiveWeightFactor > 1)
                AdaptiveWeightFactor  = 1.;
            if(AdaptiveWeightFactor < 0)
                AdaptiveWeightFactor  = 0.;
        }
    }
    else
    {
        AdaptiveWeightFactor = 1.;
    }

    return AdaptiveWeightFactor;
}


//Median of median
double M3Matrix::MedMed()
{
    M3Matrix Col;
    M3Matrix MedianLine(1,GetNbCols());

    for(int j=0;j<GetNbCols();j++)
    {
        (*this).SetCol(j,Col);
        MedianLine[0][j] = Col.Median();
    }

    double medmed = MedianLine.Median();
    return medmed;
}


//Computes the median value of a matrix
double M3Matrix::Median()
{
    int Length = t_NbRows * t_NbCols;
    double* y = new double[Length];
    for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols;j++)
        {
            y[i*t_NbCols+j] = (*this)[i][j];
        }

    qsort (y, Length, sizeof(double), my_cmp);

    double median;
    if(Length%2==0) //even number
        median = (y[Length/2]+y[Length/2-1])/2.;
    else
        median = y[Length/2];

    delete[] y;
    return median;
}


// Computes the mean square error between two matrices
// (used for measuring the error between images)
double M3Matrix::MSE(M3Matrix& B, int Margin)
{
    double tmp = 0.;

    for(int i=Margin;i<t_NbRows-Margin;i++)
        for(int j=Margin;j<t_NbCols-Margin;j++)
        {
            tmp += ((*this)[i][j] - B[i][j])*((*this)[i][j] - B[i][j]);
        }

    tmp /= (double)((t_NbRows-2*Margin)*(t_NbCols-2*Margin));

    return tmp;
}


//Computes the SSIM between two matrices
//This implementation directly follows the paper
//Image Quality Assessment: From Error Visibility to Structural Similarity
// Z. Wang, A. Bovik, H. Sheikh, E.Simoncelli
//It actually computes the mean SSIM (MSSSIM) between the images (*this) and Image2
double M3Matrix::SSIM(M3Matrix& B,int WindowSize, double& lum, double& contrast, double& structure, double K1,double K2)
{
    double L = (1ULL << static_cast<unsigned long long>(t_NbBits))-1;
    double C1 = K1*K1*(L*L);
    double C2 = K2*K2*(L*L);

    lum = 0;
    contrast = 0;
    structure = 0;

    double SSIM = 0.;
    double C3 = C2/2; //This is what they do in the paper !

    int S = WindowSize, HS = S/2; //Half size of the windows use to compute the local
    M3Matrix X,Y;
    (*this).Mirror(X,HS,HS);
    B.Mirror(Y,HS,HS);

    double mu_X, mu_Y, sigma_X,sigma_Y,sigma_XY,l,c,s;

    for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols;j++)
        {
            mu_X = 0.;
            mu_Y = 0.;
            sigma_X = 0.;
            sigma_Y = 0.;
            sigma_XY = 0.;

            for(int m=0;m<S;m++)
                for(int n=0;n<S;n++)
                {
                    mu_X +=  X[i+m][j+n];
                    mu_Y +=  Y[i+m][j+n];

                    sigma_X +=  X[i+m][j+n]*X[i+m][j+n];
                    sigma_Y +=  Y[i+m][j+n]*Y[i+m][j+n];
                }

            mu_X /= (S*S);
            mu_Y /= (S*S);
            sigma_X /= (S*S);
            sigma_Y /= (S*S);

            sigma_X -= mu_X*mu_X;
            sigma_Y -= mu_Y*mu_Y;

            for(int m=0;m<S;m++)
                for(int n=0;n<S;n++)
                {
                    sigma_XY +=  (X[i+m][j+n] - mu_X)*(Y[i+m][j+n]-mu_Y);
                }

            sigma_XY /= (S*S); //covariance

            //standard deviations
            sigma_X = sqrt(sigma_X);
            sigma_Y = sqrt(sigma_Y);

            //now compute the three factors
            l = (2 * mu_X *mu_Y +C1)/(mu_X*mu_X + mu_Y*mu_Y + C1);
            c = (2 * sigma_X *sigma_Y +C2)/(sigma_X*sigma_X + sigma_Y*sigma_Y + C2);
            s = (sigma_XY + C3)/(sigma_X*sigma_Y + C3);

            lum +=l;
            contrast += c;
            structure += s;

            SSIM += l*c*s; //the local SSIM is added
        }


    lum /= (t_NbRows*t_NbCols);
    contrast /= (t_NbRows*t_NbCols);
    structure /= (t_NbRows*t_NbCols);
    SSIM /= (t_NbRows*t_NbCols);

    return SSIM;
}


// Computes the PSNR between two matrices
double M3Matrix::ColorPSNR(M3Matrix& RGB, int Margin)
{
    double PSNR;
    M3Matrix R,G,B,Ro,Go,Bo;
    RGB.ExtractRGBChannels(R,G,B);
    (*this).ExtractRGBChannels(Ro,Go,Bo);

    double MSE = ( Ro.MSE(R, Margin) + Go.MSE(G, Margin)  + Bo.MSE(B, Margin) )/3.0;
    PSNR = 10. * log10(t_Max*t_Max/MSE);

    return PSNR;
}


// Computes the PSNR between two matrices
double M3Matrix::PSNR(M3Matrix& M3Matrix2, int Margin)
{
    double PSNR;
    double MSE = (*this).MSE(M3Matrix2, Margin);
    PSNR = 10. * log10(t_Max*t_Max/MSE);

    return PSNR;
}


double PSNR(M3Matrix& Im1, M3Matrix& Im2, int Margin)
{
    double PSNR;
    int Max = Im1.GetRangeMax();
    double MSE = Im1.MSE( Im2, Margin);
    PSNR = 10. * log10(Max*Max/MSE);

    return PSNR;
}


double ColorPSNR(M3Matrix& R1, M3Matrix& G1, M3Matrix& B1, M3Matrix& R2, M3Matrix& G2, M3Matrix& B2, int Margin)
{
    int Max = R1.GetRangeMax();
    double PSNR;
    double MSE_R = R1.MSE( R2, Margin);
    double MSE_G = G1.MSE( G2, Margin);
    double MSE_B = B1.MSE( B2, Margin);

    double MSE = (MSE_R + MSE_G + MSE_B)/3.0;
    PSNR = 10. * log10(Max*Max/MSE);

    return PSNR;
}

// Computes the histogram of a signal and  returns the entropy of the signal
// WARNING (VERY IMPORTANT): only useful for symbol matrices
double M3Matrix::Histogram(M3Matrix& Histo)
{
    int Min = DBL2INT(GetMin());
    int Max = DBL2INT(GetMax());
    double entropy = 0.;
    Histo.Reshape(Max-Min+1,1);

    for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols;j++)
            Histo[DBL2INT((*this)[i][j])-Min][0] += 1;

    double tmp = 0.;
    double Proba = 0.;
    double NbCoeff = (double)(t_NbRows * t_NbCols);

    for(int i=0;i<Histo.GetNbRows();i++)
    {
        Proba = Histo[i][0]/NbCoeff;
        if(Proba != 0.)
            tmp += -Proba * log(Proba)/log(2.);
    }

    entropy = tmp;
    //double cost = entropy * NbCoeff;
    return entropy;
}


void M3Matrix::Autocorrelation(M3Matrix& A,double mux, int MaxShift)
{
    M3Matrix Tmp;
    Tmp.Correlation(A,A,mux,mux,MaxShift);

    (*this).Reshape(MaxShift+1,MaxShift+1);
    for(int i = 0;i<=MaxShift;i++)
        for(int j = 0;j<=MaxShift;j++)
        {
            (*this)[i][j] = Tmp[i+MaxShift][j+MaxShift];
        }
}


// Computes the autocorrelation matrix between the matrices A and B
void M3Matrix::Correlation(M3Matrix& A,M3Matrix& B, double mux, double muy, int MaxShift)
{
    (*this).Reshape(2*MaxShift+1,2*MaxShift+1);
    double tmp;
    int mMin,mMax,nMin,nMax;

    for(int m=-MaxShift;m<=MaxShift;m++)
        for(int n=-MaxShift;n<=MaxShift;n++)
        {
            tmp = 0.;

            mMin = max(0,-m);
            nMin = max(0,-n);
            mMax = A.GetNbRows()-max(0,m);
            nMax = A.GetNbCols()-max(0,n);

            for(int i = mMin;i<mMax;i++)
                for(int j = nMin;j<nMax;j++)
                {
                    tmp += (A[i][j] -mux) * (B[i+m][j+n]-muy);
                }

            (*this)[m+MaxShift][n+MaxShift] = tmp;
        }
}



double M3Matrix::MeanOfMaxGradient(double p)
{
    double tmpV=0,tmpH = 0,tmp;
    double* HorizontalGradients = new double[t_NbRows*(t_NbCols-1)];
    double* VerticalGradients = new double[(t_NbRows-1)*t_NbCols];


    for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols-1;j++)
            HorizontalGradients[i*(t_NbCols-1)+j] = fabs ( (*this)[i][j+1] - (*this)[i][j]);

    for(int i=0;i<t_NbRows-1;i++)
        for(int j=0;j<t_NbCols;j++)
            VerticalGradients[i*t_NbCols+j] = fabs ((*this)[i+1][j] - (*this)[i][j]);


    qsort(HorizontalGradients,t_NbRows*(t_NbCols-1),sizeof(double),my_cmp_descent);
    qsort(VerticalGradients,(t_NbRows-1)*t_NbCols,sizeof(double),my_cmp_descent);

    if(p<0) p=0;
    if(p<1) p=1;

    int  NbMaxValuesH = int(p * t_NbRows* (t_NbCols-1)); //default 1% of the largest values is taken
    int  NbMaxValuesV = int(p * (t_NbRows-1)* t_NbCols); //default 1% of the largest values is taken

    for(int i=0;i<NbMaxValuesH;i++)
    {
        tmpH += HorizontalGradients[i];
    }
    for(int i=0;i<NbMaxValuesH;i++)
    {
        tmpV += VerticalGradients[i];
    }

    tmpH /= NbMaxValuesH;
    tmpV /= NbMaxValuesV;

    tmp = 0.5*(tmpH+tmpV);

    delete[] HorizontalGradients;
    delete[] VerticalGradients;
    return tmp;
}


double M3Matrix::Sobolev2Norm(int p)
{
    double tmp = 0., Diff;

    //Variation along the vertical directions
    for(int i = 0;i<t_NbRows-2;i++)
    {
        for(int j = 0;j<t_NbCols;j++)
        {
            Diff = fabs((*this)[i][j] - 2* (*this)[i+1][j]  +  (*this)[i+2][j]);
            tmp += pow(Diff,(double)p);
        }
    }

    for(int i = 0;i<t_NbRows;i++)
    {
        for(int j = 0;j<t_NbCols-2;j++)
        {
            Diff = fabs((*this)[i][j] + (*this)[i][j+2]  - 2 * (*this)[i][j+1]);
            tmp += pow(Diff,(double)p);
        }
    }

    tmp = pow(tmp,1./((double)p));
    return tmp;
}


double M3Matrix::SobolevNorm(int p)
{
    double tmp = 0., Diff;

    //Finite difference along the vertical directions
    for(int i = 0;i<t_NbRows-1;i++)
    {
        for(int j = 0;j<t_NbCols;j++)
        {
            Diff = fabs((*this)[i+1][j]  - (*this)[i][j] );
            tmp += pow(Diff,(double)p);
        }
    }

    //Finite difference along the horizontal directions
    for(int i = 0;i<t_NbRows;i++)
    {
        for(int j = 0;j<t_NbCols-1;j++)
        {
            Diff = fabs((*this)[i][j+1] - (*this)[i][j]);
            tmp += pow(Diff,(double)p);
        }
    }

    tmp = pow(tmp,1./((double)p));
    return tmp;
}


// Fills up the matrix with the values
// from a TXT container with following convention
// rows cols
// x_{1,1} x_{1,2} ... _x_{1,cols}
// ...
// x_{rows,1} x_{rows,2} ... _x_{rows,cols}
void M3Matrix::LoadFromTXT(const string& filenameIn)
{
    ifstream InStr(filenameIn.c_str());
    int NbRows,NbCols;

    InStr >> NbRows;
    InStr >> NbCols;

    Reshape(NbRows,NbCols);

    double val;
    for(int i = 0;i<NbRows;i++)
        for(int j = 0;j<NbCols;j++)
        {
            InStr >> val;
            (*this)[i][j] = val;
        }
}


void M3Matrix::WriteTXT(const string& filenameOut, TXTDIRECTION DataDirection)
{
    ofstream out_str(filenameOut);

    if(DataDirection == vertical)
    {
        out_str << t_NbRows << endl;
        out_str << t_NbCols << endl;

        for(int i=0;i<t_NbRows;i++)
        {
            for(int j=0;j<t_NbCols;j++)
            {
                out_str << (*this)[i][j] << " ";
            }
            out_str << endl;
        }
    }
    else if(DataDirection == horizontal)
    {
        out_str <<  t_NbCols << endl;
        out_str <<  t_NbRows << endl;

        for(int j=0;j<t_NbCols;j++)
        {
            for(int i=0;i<t_NbRows;i++)
            {
                out_str << (*this)[i][j] << " ";
            }
            out_str << endl;
        }
    }
}


void M3Matrix::ParsePXMMetaData(string line_str)
{
    stringstream ss_line(line_str);
    string word1;
    string word2;

    ss_line >> word1;
    ss_line >> word2;

    if(!strcmp(word1.c_str(),"BP"))
    {
        if(strcmp(word2.c_str(),"RGGB") && strcmp(word2.c_str(),"GRBG")
                && strcmp(word2.c_str(),"GBRG") && strcmp(word2.c_str(),"BGGR")
                && strcmp(word2.c_str(),"None") && strcmp(word2.c_str(),"Undefined") )
            m_BayerPattern.assign("Undefined");
        else
            m_BayerPattern.assign(word2);
    }
}


void M3Matrix::LoadRawFromXRGB(const char* infile)
{
    ifstream pFile(infile,std::ios_base::binary);
    if (!pFile)
    {
        cout << "Can't open the file "<< infile << " (does the file exist ?). Abort ... " << endl;
        exit(1);
    }

    int NbRows = 968, NbCols = 1288; //Currently hard coded, this corresponds to the size of the sensor
    InitialiseBitDepth(10);
    Reshape(NbRows,NbCols);

    char val1,val2,val3,val4;

    for(int i = 0;i < NbRows; i++)
        for(int j = 0;j < NbCols; j++)
        {
            pFile.get(val1);
            pFile.get(val2);
            pFile.get(val3);
            pFile.get(val4);

            int msb = (int) ((unsigned char) val1);
            int lsb = (int) ((unsigned char) (val4 & 0x3));

            int pixel = (msb << 2) + lsb;
            (*this)[i][j] = pixel;
        }
}


// At each block of 4 bytes XRGB stores the information as follows
// byte 1: MSB_B (blue most significant bits)
// byte 2: MSB_G (green most significant bits)
// byte 3: MSB_R (red most significant bits)
// byte 4: 2 dummy bits, 2 least significant bits for r,g, and b respectively
void M3Matrix::LoadFromXRGB(const char* infile)
{
    ifstream pFile(infile,std::ios_base::binary);
    if (!pFile)
    {
        cout << "Can't open the file "<< infile << " (does the file exist ?). Abort ... " << endl;
        exit(1);
    }

    int NbRows = 968, NbCols = 1288; //Currently hard coded, this corresponds to the size of the sensor
    M3Matrix R(NbRows,NbCols),G(NbRows,NbCols),B(NbRows,NbCols);

    R.InitialiseBitDepth(10);
    G.InitialiseBitDepth(10);
    B.InitialiseBitDepth(10);

    char val_b, val_g, val_r, val_lsb;

    for(int i = 0;i < NbRows; i++)
        for(int j = 0;j < NbCols; j++)
        {
            pFile.get(val_b);
            pFile.get(val_g);
            pFile.get(val_r);
            pFile.get(val_lsb);

            int msbr = (int) ((unsigned char) val_r);
            int msbg = (int) ((unsigned char) val_g);
            int msbb = (int) ((unsigned char) val_b);
            int lsbr = ((int) ((unsigned char) (val_lsb)) >> 4) &0x3;
            int lsbg = ((int) ((unsigned char) (val_lsb)) >> 2) &0x3;
            int lsbb = ((int) ((unsigned char) (val_lsb)) >> 0) &0x3;

            R[i][j] = (msbr << 2) + lsbr;
            G[i][j] = (msbg << 2) + lsbg;
            B[i][j] = (msbb << 2) + lsbb;
        }

    CombineRGBChannels(R,G,B);
}


// Fills up the matrix with the luminance values
// read in a PGM file
void M3Matrix::LoadFromPGM(const char* infile)
{
    ifstream pFile(infile,std::ios_base::binary);
    if (!pFile)
    {
        cout << "Can't open the file "<< infile << " (does the file exist ?). Abort ... " << endl;
        exit(1);
    }

    char Buffer;
    string line_str;

    pFile.seekg(0);
    pFile.get(Buffer); //First character
    while(Buffer == '#') //Skip comment lines
    {
        getline(pFile,line_str);
        pFile.get(Buffer);
    }

    //Buffer now contains the first non commented character
    if(Buffer != 'P') //not correct pgm file
    {
        cout << "Error: this is not a  pgm file, abort ..." << endl;
        exit(1);
    }

    int binary;
    pFile.get(Buffer);
    if(Buffer !='2' && Buffer !='5')
    {
        cout 	<< "Error in LoadFromPGM: file header is not compatible with pgm specifications, abort ... " << endl;
        exit(1);
    }
    else
    {
        if(Buffer =='2') //classical, ASCII case
        {
            binary = 0;
        }
        else //binary case
        {
            binary = 1;
        }
    }

    Buffer = pFile.peek();
    while(isspace(Buffer))
    {
        Buffer = pFile.peek();
        if(isspace(Buffer))
        {
            pFile.get(Buffer);
        }
    }

    Buffer = pFile.peek();
    while(Buffer == '#')
    {
        pFile.get(Buffer);
        getline(pFile,line_str);
        ParsePXMMetaData(line_str);
        Buffer = pFile.peek();
    }

    int NbRows,NbCols;
    double MaxValue;

    pFile >> NbCols;

    Buffer = pFile.peek();
    while(isspace(Buffer))
    {
        Buffer = pFile.peek();
        if(isspace(Buffer))
        {
            pFile.get(Buffer);
        }
    }

    Buffer = pFile.peek();
    while(Buffer == '#')
    {
        pFile.get(Buffer);
        getline(pFile,line_str);
        ParsePXMMetaData(line_str);
        Buffer = pFile.peek();
    }

    pFile >> NbRows;

    Buffer = pFile.peek();
    while(isspace(Buffer))
    {
        Buffer = pFile.peek();
        if(isspace(Buffer))
        {
            pFile.get(Buffer);
        }
    }

    Buffer = pFile.peek();
    while(Buffer == '#')
    {
        pFile.get(Buffer);
        getline(pFile,line_str);
        ParsePXMMetaData(line_str);
        Buffer = pFile.peek();
    }

    pFile >> MaxValue;

    t_Min = 0;
    t_Max = MaxValue;
    t_NbBits =(int) (log(MaxValue + 1)/log(2.0) + 0.5);

    (*this).Reshape(NbRows,NbCols);
    if(binary==0) //ASCII file, any it depth is allowed
    {
        for(int i=0;i<t_NbRows;i++)
            for(int j=0;j<t_NbCols;j++)
                pFile >> (*this)[i][j];
    }
    else //binary
    {
        if(t_NbBits <=8 )
        {
            char valc, dummy;
            pFile.get(dummy); //This should be the last line

            for(int i=0;i<t_NbRows;i++)
                for(int j=0;j<t_NbCols;j++)
                {
                    pFile.get(valc);
                    (*this)[i][j] = (int) ((unsigned char) valc);
                }
        }
        else if(t_NbBits <= 16)
        {
            char valcl,valcr, dummy;
            pFile.get(dummy); //This should be the last line

            for(int i=0;i<t_NbRows;i++)
                for(int j=0;j<t_NbCols;j++)
                {
                    pFile.get(valcl);
                    pFile.get(valcr);
                    (*this)[i][j] = (int) ((unsigned char) valcr);
                    int vall = (int) ((unsigned char) valcl);
                    (*this)[i][j] += (vall << 8);
                }
        }
    }
}



void M3Matrix::LoadFromXRGB(const string& infile)
{
    LoadFromXRGB(infile.c_str());
}

void M3Matrix::LoadRawFromXRGB(const string& infile)
{
    LoadRawFromXRGB(infile.c_str());
}


void M3Matrix::LoadFromPGM(const string& infile)
{
    LoadFromPGM(infile.c_str());
}


void M3Matrix::LoadFromPPM(const string& infile)
{
    LoadFromPPM(infile.c_str());
}


// Fills up the matrix with the color values
// read in a PPM file
void M3Matrix::LoadFromPPM(const char* infile)
{    
    ifstream pFile(infile,std::ios_base::binary);
    if (!pFile)
    {
        cout << "Can't open the file "<< infile << " (does the file exist ?). Abort ... " << endl;
        exit(1);
    }

    char Buffer;
    string line_str;

    pFile.seekg(0);
    pFile.get(Buffer); //First character
    while(Buffer == '#') //Skip comment lines
    {
        while(Buffer != '\n')
        {
            pFile.get(Buffer);
        }
        pFile.get(Buffer);
    }

    //Buffer now contains the first non commented character
    if(Buffer != 'P') //file not compatible with the ppm specification
    {
        cout << "Error: this is not a  ppm file, abort ..." << endl;
        exit(1);
    }

    int binary;
    pFile.get(Buffer);
    if(Buffer !='3' && Buffer !='6')
    {
        cout 	<< "Error in LoadFromPPM: file header is not compatible with ppm specifications, abort ... " << endl;
        exit(1);
    }
    else
    {
        if(Buffer =='3') //classical, ASCII case
        {
            binary = 0;
        }
        else //binary case
        {
            binary = 1;
        }
    }

    Buffer = pFile.peek();
    while(isspace(Buffer))
    {
        Buffer = pFile.peek();
        if(isspace(Buffer))
        {
            pFile.get(Buffer);
        }
    }

    Buffer = pFile.peek();
    while(Buffer == '#')
    {
        pFile.get(Buffer);
        getline(pFile,line_str);
        ParsePXMMetaData(line_str);
        Buffer = pFile.peek();
    }


    int NbRows,NbCols;
    double MaxValue;
    pFile >> NbCols;

    Buffer = pFile.peek();
    while(isspace(Buffer))
    {
        Buffer = pFile.peek();
        if(isspace(Buffer))
        {
            pFile.get(Buffer);
        }
    }

    Buffer = pFile.peek();
    while(Buffer == '#')
    {
        pFile.get(Buffer);
        getline(pFile,line_str);
        ParsePXMMetaData(line_str);
        Buffer = pFile.peek();
    }

    pFile >> NbRows;

    Buffer = pFile.peek();
    while(isspace(Buffer))
    {
        Buffer = pFile.peek();
        if(isspace(Buffer))
        {
            pFile.get(Buffer);
        }
    }

    Buffer = pFile.peek();
    while(Buffer == '#')
    {
        pFile.get(Buffer);
        getline(pFile,line_str);
        ParsePXMMetaData(line_str);
        Buffer = pFile.peek();
    }

    pFile >> MaxValue;


    t_Min = 0;
    t_Max = MaxValue;
    t_NbBits =(int) (log(MaxValue + 1)/log(2.0) + 0.5);

    (*this).Reshape(NbRows,3*NbCols); //NbCols is the number of columns of one single channel
    if(binary==0) //ASCII modus
    {
        for(int i=0;i<t_NbRows;i++)
            for(int j=0;j<NbCols;j++)
            {
                pFile >> (*this)[i][j];
                pFile >> (*this)[i][j+NbCols];
                pFile >> (*this)[i][j+2*NbCols];
            }
    }
    else //binary modus
    {
        if(t_NbBits <= 8)
        {
            char valcR,valcG,valcB;
            pFile.get(valcR);

            for(int i=0;i<t_NbRows;i++)
                for(int j=0;j<NbCols;j++)
                {
                    //Red Value
                    pFile.get(valcR);
                    (*this)[i][j] = (int) ((unsigned char) valcR);

                    //Green Value
                    pFile.get(valcG);
                    (*this)[i][j+NbCols] = (int) ((unsigned char) valcG);

                    //Blue Value
                    pFile.get(valcB);
                    (*this)[i][j+2*NbCols] = (int) ((unsigned char) valcB);
                }
        }
        else if(t_NbBits <= 16 )
        {
            char valcRl,valcRr,valcGl,valcGr,valcBl,valcBr;
            char Dummy;
            pFile.get(Dummy);

            for(int i=0;i<t_NbRows;i++)
                for(int j=0;j<NbCols;j++)
                {
                    //Red Value
                    pFile.get(valcRl);
                    pFile.get(valcRr);
                    (*this)[i][j] = (int) ((unsigned char) valcRr);
                    int valRl = (int) ((unsigned char) valcRl);
                    (*this)[i][j] += (valRl << 8);

                    //Green Value
                    pFile.get(valcGl);
                    pFile.get(valcGr);
                    (*this)[i][j+NbCols] = (int) ((unsigned char) valcGr);
                    int valGl = (int) ((unsigned char) valcGl);
                    (*this)[i][j+NbCols] += (valGl << 8);

                    //Blue Value
                    pFile.get(valcBl);
                    pFile.get(valcBr);
                    (*this)[i][j+2*NbCols] = (int) ((unsigned char) valcBr);
                    int valBl = (int) ((unsigned char) valcBl);
                    (*this)[i][j+2*NbCols] += (valBl << 8);
                }
        }
    }
}


void M3Matrix::ExtractRGBChannels(M3Matrix& R,M3Matrix& G,M3Matrix& B)
{
    //Warning: this must be a matrix
    int NbRows = t_NbRows;
    int NbCols = t_NbCols/3;

    R.Reshape(NbRows,NbCols);
    G.Reshape(NbRows,NbCols);
    B.Reshape(NbRows,NbCols);

    R.InitialiseDynamic(*this);
    G.InitialiseDynamic(*this);
    B.InitialiseDynamic(*this);

    for(int i=0;i<NbRows;i++)
        for(int j=0;j<NbCols;j++)
        {
            R[i][j] = (*this)[i][j];
            G[i][j] = (*this)[i][j+NbCols];
            B[i][j] = (*this)[i][j+2*NbCols];
        }
}


void M3Matrix::CombineInterlacedRGBChannels(M3Matrix& R,M3Matrix& G,M3Matrix& B)
{
    (*this).InitialiseDynamic(R);

    int NbRows = R.GetNbRows();
    int NbCols = R.GetNbCols();

    (*this).Reshape(NbRows,3*NbCols);

    for(int i=0;i<NbRows;i++)
        for(int j=0;j<NbCols;j++)
        {
            (*this)[i][3*j] = R[i][j];
            (*this)[i][3*j+1] = G[i][j] ;
            (*this)[i][3*j+2] = B[i][j];
        }
}


void M3Matrix::CombineRGBChannels(M3Matrix& R,M3Matrix& G,M3Matrix& B)
{
    (*this).InitialiseDynamic(R);

    int NbRows = R.GetNbRows();
    int NbCols = R.GetNbCols();

    (*this).Reshape(NbRows,3*NbCols);

    for(int i=0;i<NbRows;i++)
        for(int j=0;j<NbCols;j++)
        {
            (*this)[i][j] = R[i][j];
            (*this)[i][j+NbCols] = G[i][j] ;
            (*this)[i][j+2*NbCols] = B[i][j];
        }
}


M3Matrix M3Matrix::HAJapanDebayering(int BayerPattern, int Threshold, int epsilon)
{
    M3Matrix Debayered,R,G,B;
    Debayered.InitialiseDynamic(*this);
    (*this).HAJapanDebayering(R,G,B,BayerPattern, Threshold, epsilon);
    Debayered.CombineRGBChannels(R,G,B);
    return Debayered;
}


M3Matrix M3Matrix::HAPrefDebayering(int BayerPattern)
{
    M3Matrix Debayered,R,G,B;
    Debayered.InitialiseDynamic(*this);
    (*this).HAPrefDebayering(R,G,B,BayerPattern);
    Debayered.CombineRGBChannels(R,G,B);
    return Debayered;
}


M3Matrix M3Matrix::HADebayering(int BayerPattern)
{
    M3Matrix Debayered,R,G,B;
    Debayered.InitialiseDynamic(*this);
    (*this).HADebayering(R,G,B,BayerPattern);
    Debayered.CombineRGBChannels(R,G,B);
    return Debayered;
}


// Implementation taken from Liu, Shen, Yu
// An Improved Demosaicing Algorithm by Adopting Color Correlation Aided Gradients
void M3Matrix::HADebayering(M3Matrix& R,M3Matrix& G,M3Matrix& B, int BayerPattern)
{
    R.Reshape(t_NbRows,t_NbCols);
    G.Reshape(t_NbRows,t_NbCols);
    B.Reshape(t_NbRows,t_NbCols);

    int ExtSize = 2; //Must be even !!!
    int NbRowsExt = t_NbRows+2*ExtSize;
    int NbColsExt = t_NbCols+2*ExtSize;

    M3Matrix RawMirrored;
    (*this).Mirror(RawMirrored,ExtSize,ExtSize);

    M3Matrix RE,GE,BE; // extended color matrices
    RE.Reshape(NbRowsExt,NbColsExt);
    GE.Reshape(NbRowsExt,NbColsExt);
    BE.Reshape(NbRowsExt,NbColsExt);

    int ri=0,rj=0,bi=0,bj=0;
    if(BayerPattern == 0)
    {
        ri = 0;  rj = 0;
    }
    else if(BayerPattern == 1)
    {
        ri = 0;  rj = 1;
    }
    else if(BayerPattern == 2)
    {
        ri = 1;  rj = 0;
    }
    else if(BayerPattern == 3)
    {
        ri = 1;  rj = 1;
    }

    bi = 1-ri;  bj = 1-rj;

    //Fill the colors where available
    for(int i=0;i<NbRowsExt;i++)
        for(int j=0;j<NbColsExt;j++)
        {
            if(i%2 == ri && j%2 == rj) //red pixels
            {
                RE[i][j] = RawMirrored[i][j];
            }
            else if(i%2 == bi && j%2 == bj) //blue pixels
            {
                BE[i][j] = RawMirrored[i][j];
            }
            else //green pixels
            {
                GE[i][j] = RawMirrored[i][j];
            }
        }

    //Now begin the interpolation algorithm for the green channel
    double Lh,Lv;

    //double alphaR = 1.0, alphaG = 1.0, alphaB = 1.0;
    for(int i=ExtSize;i<t_NbRows+ExtSize;i++)
        for(int j=ExtSize;j<t_NbCols+ExtSize;j++)
        {
            if(i%2 == ri && j%2 == rj) //Red pixels
            {
                Lh = fabs(GE[i][j+1]-GE[i][j-1]) + fabs(2*RE[i][j] - RE[i][j+2] -RE[i][j-2]);
                Lv = fabs(GE[i+1][j]-GE[i-1][j]) + fabs(2*RE[i][j] - RE[i+2][j] -RE[i-2][j]);

                if(Lh <Lv)
                    GE[i][j] = 0.5*(GE[i][j+1] + GE[i][j-1]) +0.25*(2*RE[i][j] - RE[i][j+2] -RE[i][j-2]);
                else if(Lv < Lh)
                    GE[i][j] = 0.5*(GE[i+1][j] + GE[i-1][j]) +0.25*(2*RE[i][j] - RE[i+2][j] -RE[i-2][j]);
                else
                    GE[i][j] = 0.25*(GE[i][j+1] + GE[i][j-1] + GE[i+1][j] + GE[i-1][j]) + 0.125*(4*RE[i][j] - RE[i][j+2] -RE[i][j-2] - RE[i+2][j] -RE[i-2][j]);
            }
            else if(i%2 == bi && j%2 == bj) //Blue pixels
            {
                Lh = fabs(GE[i][j+1]-GE[i][j-1]) + fabs(2*BE[i][j] - BE[i][j+2] -BE[i][j-2]);
                Lv = fabs(GE[i+1][j]-GE[i-1][j]) + fabs(2*BE[i][j] - BE[i+2][j] -BE[i-2][j]);

                if(Lh <Lv)
                    GE[i][j] = 0.5*(GE[i][j+1] + GE[i][j-1]) + 0.25*(2*BE[i][j] - BE[i][j+2] - BE[i][j-2]);
                else if(Lv < Lh)
                    GE[i][j] = 0.5*(GE[i+1][j] + GE[i-1][j]) + 0.25*(2*BE[i][j] - BE[i+2][j] -BE[i-2][j]);
                else
                    GE[i][j] = 0.25*(GE[i][j+1] + GE[i][j-1] + GE[i+1][j] + GE[i-1][j]) + 0.125*(4*BE[i][j] - BE[i][j+2] -BE[i][j-2] - BE[i+2][j] -BE[i-2][j]);
            }

            if(GE[i][j]<0)
                GE[i][j] = 0.;
            if(GE[i][j]>t_Max)
                GE[i][j] = t_Max;
        }


    //Now interpolation of red and blue value at green pixels
    double delta_n,delta_p; // this is for the diagonal gradients
    for(int i=ExtSize;i<t_NbRows+ExtSize;i++)
        for(int j=ExtSize;j<t_NbCols+ExtSize;j++)
        {
            if(i%2 == ri && j%2 == rj) //red pixel
            {
                //Estimates diagonal gradients
                delta_n = fabs(BE[i-1][j-1]-BE[i+1][j+1]) + fabs(2*GE[i][j] - GE[i-1][j-1] - GE[i+1][j+1]);
                delta_p = fabs(BE[i+1][j-1]-BE[i-1][j+1]) + fabs(2*GE[i][j] - GE[i+1][j-1] - GE[i-1][j+1]);

                //Estimates blue value
                if(delta_n < delta_p)
                    BE[i][j] = 0.5*(BE[i-1][j-1] + BE[i+1][j+1]) + 0.5*(2*GE[i][j] - GE[i-1][j-1] - GE[i+1][j+1]);
                else if(delta_p < delta_n)
                    BE[i][j] = 0.5*(BE[i+1][j-1] + BE[i-1][j+1]) + 0.5*(2*GE[i][j] - GE[i+1][j-1] - GE[i-1][j+1]);
                else
                    BE[i][j] = 0.25*(BE[i-1][j-1] + BE[i+1][j+1] + BE[i+1][j-1] + BE[i-1][j+1]) + 0.25*(4*GE[i][j] - GE[i+1][j-1] - GE[i-1][j+1] - GE[i-1][j-1] - GE[i+1][j+1]);
            }
            else if(i%2 == bi && j%2 == bj) //blue pixel
            {
                //Estimates diagonal gradients
                delta_n = fabs(RE[i-1][j-1]-RE[i+1][j+1]) + fabs(2*GE[i][j] - GE[i-1][j-1] - GE[i+1][j+1]);
                delta_p = fabs(RE[i+1][j-1]-RE[i-1][j+1]) + fabs(2*GE[i][j] - GE[i+1][j-1] - GE[i-1][j+1]);

                //Estimates red value
                if(delta_n < delta_p)
                    RE[i][j] = 0.5*(RE[i-1][j-1] + RE[i+1][j+1]) + 0.5*(2*GE[i][j] - GE[i-1][j-1] - GE[i+1][j+1]);
                else if(delta_p < delta_n)
                    RE[i][j] = 0.5*(RE[i+1][j-1] + RE[i-1][j+1]) + 0.5*(2*GE[i][j] - GE[i+1][j-1] - GE[i-1][j+1]);
                else
                    RE[i][j] = 0.25*(RE[i-1][j-1] + RE[i+1][j+1] + RE[i+1][j-1] + RE[i-1][j+1]) + 0.25*(4*GE[i][j] - GE[i+1][j-1] - GE[i-1][j+1] - GE[i-1][j-1] - GE[i+1][j+1]);
            }
            else if(i%2 == ri && j%2 == bj) //green pixel, red line, blue column
            {
                //Estimates red value
                RE[i][j] = 0.5*(RE[i][j+1] + RE[i][j-1]) + 0.5*(2*GE[i][j] - GE[i][j+1] - GE[i][j-1]);

                //Estimates blue value
                BE[i][j] = 0.5*(BE[i+1][j] + BE[i-1][j]) + 0.5*(2*GE[i][j] - GE[i+1][j] -  GE[i-1][j]);
            }
            else if(i%2 == bi && j%2 == rj) //green pixel, blue line, red column
            {
                //Estimates red value
                RE[i][j] = 0.5*(RE[i+1][j] + RE[i-1][j]) + 0.5*(2*GE[i][j] - GE[i+1][j] -  GE[i-1][j]);

                //Estimates blue value
                BE[i][j] = 0.5*(BE[i][j+1] + BE[i][j-1]) + 0.5*(2*GE[i][j] - GE[i][j+1] -  GE[i][j-1]);
            }

            if(BE[i][j]<0)
                BE[i][j] = 0.;
            if(BE[i][j]>t_Max)
                BE[i][j] = t_Max;

            if(RE[i][j]<0)
                RE[i][j] = 0.;
            if(RE[i][j]>t_Max)
                RE[i][j] = t_Max;
        }

    RE.ExtractBlock(R,ExtSize,ExtSize,t_NbRows,t_NbCols);
    GE.ExtractBlock(G,ExtSize,ExtSize,t_NbRows,t_NbCols);
    BE.ExtractBlock(B,ExtSize,ExtSize,t_NbRows,t_NbCols);
}


void M3Matrix::ExperimentalDebayering(M3Matrix& R,M3Matrix& G,M3Matrix& B, int BayerPattern)
{  
    R.Reshape(t_NbRows,t_NbCols);
    G.Reshape(t_NbRows,t_NbCols);
    B.Reshape(t_NbRows,t_NbCols);

    int ExtSize = 2; //Must be even !!!
    int NbRowsExt = t_NbRows+2*ExtSize;
    int NbColsExt = t_NbCols+2*ExtSize;

    M3Matrix RawMirrored;
    (*this).Mirror(RawMirrored,ExtSize,ExtSize);

    M3Matrix RE,GE,BE; // Extended color matrices
    RE.Reshape(NbRowsExt,NbColsExt);
    GE.Reshape(NbRowsExt,NbColsExt);
    BE.Reshape(NbRowsExt,NbColsExt);

    int ri=0,rj=0,bi=0,bj=0;
    if(BayerPattern == 0)
    {
        ri = 0;  rj = 0;
    }
    else if(BayerPattern == 1)
    {
        ri = 0;  rj = 1;
    }
    else if(BayerPattern == 2)
    {
        ri = 1;  rj = 0;
    }
    else if(BayerPattern == 3)
    {
        ri = 1;  rj = 1;
    }

    bi = 1-ri;  bj = 1-rj;

    //Fill the colors where available
    for(int i=0;i<NbRowsExt;i++)
        for(int j=0;j<NbColsExt;j++)
        {
            if(i%2 == ri && j%2 == rj) //red pixels
            {
                RE[i][j] = RawMirrored[i][j];
            }
            else if(i%2 == bi && j%2 == bj) //blue pixels
            {
                BE[i][j] = RawMirrored[i][j];
            }
            else //green pixels
            {
                GE[i][j] = RawMirrored[i][j];
            }
        }

    M3Matrix Rh(NbRowsExt,NbColsExt),Rv(NbRowsExt,NbColsExt);
    M3Matrix Gh(NbRowsExt,NbColsExt),Gv(NbRowsExt,NbColsExt);
    M3Matrix Bh(NbRowsExt,NbColsExt),Bv(NbRowsExt,NbColsExt);

    double grad_v,lap_v,grad_h,lap_h;
    double grad_d1,lap_d1,grad_d2,lap_d2;

    //Now begin the interpolation algorithm for the green channel
    for(int i=ExtSize;i<t_NbRows+ExtSize;i++)
        for(int j=ExtSize;j<t_NbCols+ExtSize;j++)
        {
            if(i%2 == ri && j%2 == rj) //Red original pixel
            {
                Gh[i][j] = 0.5 * (GE[i][j-1] + GE[i][j+1]) + (2*RE[i][j] - RE[i][j-2] - RE[i][j+2]) * 0.25; //Horizontal estimate
                Gv[i][j] = 0.5 * (GE[i-1][j] + GE[i+1][j]) + (2*RE[i][j] - RE[i-2][j] - RE[i+2][j]) * 0.25; //Vertical estimate

                Rh[i][j] = RE[i][j];
                Rv[i][j] = RE[i][j];
            }
            else if(i%2 == bi && j%2 == bj) //Blue original pixel
            {
                Gv[i][j] = 0.5*(GE[i-1][j] + GE[i+1][j]) + 0.25 *  (2*BE[i][j] - BE[i-2][j] - BE[i+2][j]); //Vertical estimate
                Gh[i][j] = 0.5*(GE[i][j-1] + GE[i][j+1]) + 0.25 *  (2*BE[i][j] - BE[i][j-2] - BE[i][j+2]); //Horizontal estimate

                Bh[i][j] = BE[i][j];
                Bv[i][j] = BE[i][j];
            }
            else if(i%2 == ri && j%2 == bj) //Green original pixel, red line, blue column
            {
                // Green original value available
                Gv[i][j] = GE[i][j];
                Gh[i][j] = GE[i][j];

                Rh[i][j] = 0.5*(RE[i][j-1] + RE[i][j+1]) + 0.25 * (2*GE[i][j] - GE[i][j-2] - GE[i][j+2]); //Horizontal red estimate
                Bv[i][j] = 0.5*(BE[i-1][j] + BE[i+1][j]) + 0.25 * (2*GE[i][j] - GE[i-2][j] - GE[i+2][j]); //Vertical blue estimate
            }
            else if(i%2 == bi && j%2 == rj) //Green original pixel, blue line, red column
            {
                // Green original value available
                Gv[i][j] = GE[i][j];
                Gh[i][j] = GE[i][j];

                Rv[i][j] = 0.5*(RE[i-1][j] + RE[i+1][j]) + 0.25 * (2*GE[i][j] - GE[i-2][j] - GE[i+2][j]); //Vertical red estimate
                Bh[i][j] = 0.5*(BE[i][j-1] + BE[i][j+1]) + 0.25 * (2*GE[i][j] - GE[i][j-2] - GE[i][j+2]); //Horizontal estimate
            }
        }

    M3Matrix Delta_GR_h,Delta_GR_v,Delta_GB_h,Delta_GB_v;
    Delta_GR_h = Gh - Rh;
    Delta_GR_v = Gv - Rv;
    Delta_GB_h = Gh - Bh;
    Delta_GB_v = Gv - Bv;

    //Now smooth the differences
    M3Matrix HF, VF;
    //double sigma = 3.6;

    int BilMaskSize = 11;
    HF.Reshape(1,BilMaskSize); HF.SetValues(1.0/((double) BilMaskSize ));
    VF.Reshape(BilMaskSize,1); VF.SetValues(1.0/((double) BilMaskSize ));

    M3Matrix Tmp;
    BilParams Params;
    Params.a = -0.3;
    Params.b = 0.1;
    Params.Dynamic = R.GetDynamic();
    Params.setAWType(3);

    Params.norm_type = BilParams::NormType::exact;

    Delta_GR_h.AdaptiveBilateralFilter(Tmp,HF,Params);
    Delta_GR_h = Tmp;
    Delta_GB_h.AdaptiveBilateralFilter(Tmp,HF,Params);
    Delta_GB_h = Tmp;

    //Vertical smoothing of the vertical estimate for the differences
    //HF.MakeGaussianFilter_1D(sigma,5);
    //HF.Transpose(VF);

    Delta_GR_v.AdaptiveBilateralFilter(Tmp,VF,Params);
    Delta_GR_v = Tmp;
    Delta_GB_v.AdaptiveBilateralFilter(Tmp,VF,Params);
    Delta_GB_v = Tmp;

    for(int i=ExtSize;i<t_NbRows+ExtSize;i++)
        for(int j=ExtSize;j<t_NbCols+ExtSize;j++)
        {
            if( ( i%2 == ri && j%2 == rj) || (i%2 == bi && j%2 == bj)  ) //Red or blue pixel
            {
                grad_v = GE[i+1][j] - GE[i-1][j]; //Vertical gradient
                grad_h = GE[i][j+1] - GE[i][j-1]; //Horizontal gradient
                if(i%2 == ri && j%2 == rj) //Red pixels
                {
                    lap_v = 2*RE[i][j] - RE[i-2][j] - RE[i+2][j]; //Vertical Laplace
                    lap_h = 2*RE[i][j] - RE[i][j-2] - RE[i][j+2]; //Horizontal Laplace

                    //This is similar to the Hamilton-Adams test
                    if(fabs(grad_v) +  fabs(lap_v)  <   (fabs(grad_h) +  fabs(lap_h)) )
                    {
                        GE[i][j] = Rv[i][j] + Delta_GR_v[i][j];
                    }
                    else if( fabs(grad_h) + fabs(lap_h)   < ( fabs(grad_v) +  fabs(lap_v) ) )
                    {
                        GE[i][j] = Rh[i][j] + Delta_GR_h[i][j];
                    }
                    else
                    {
                        GE[i][j] = 0.5*(Gv[i][j] + Gh[i][j]);
                    }
                }
                else if(i%2 == bi && j%2 == bj) //Blue pixels
                {
                    lap_v = 2*BE[i][j] - BE[i-2][j] - BE[i+2][j]; //Vertical Laplace
                    lap_h = 2*BE[i][j] - BE[i][j-2] - BE[i][j+2]; //Horizontal Laplace

                    //This is similar to the Hamilton-Adams test
                    if(fabs(grad_v) + fabs(lap_v)  < fabs(grad_h) + fabs(lap_h))
                    {
                        GE[i][j] = Bv[i][j] + Delta_GB_v[i][j];
                    }
                    else if(fabs(grad_h) + fabs(lap_h) < fabs(grad_v) +  fabs(lap_v)   )
                    {
                        GE[i][j] = Bh[i][j] + Delta_GB_h[i][j];
                    }
                    else
                    {
                        GE[i][j] = 0.5*(Gv[i][j] + Gh[i][j]);
                    }
                }
            }

            if(GE[i][j] < 0)
                GE[i][j] = 0.;

            if(GE[i][j] > t_Max)
                GE[i][j] = t_Max;
        }


    //Now begin the interpolation algorithm for the blue and red channel
    for(int i=ExtSize;i<t_NbRows+ExtSize;i++)
        for(int j=ExtSize;j<t_NbCols+ExtSize;j++)
        {
            if(i%2 == ri && j%2 == rj) //Red pixels
            {
                grad_d1 = BE[i+1][j+1]-BE[i-1][j-1];
                lap_d1 = 2*GE[i][j] - GE[i+1][j+1] -GE[i-1][j-1];
                grad_d2 = BE[i-1][j+1]-BE[i+1][j-1];
                lap_d2 = 2*GE[i][j] - GE[i-1][j+1] -GE[i+1][j-1];

                //Interpolation of blue
                if(fabs(grad_d1)+ fabs(lap_d1) < fabs(grad_d2)+ fabs(lap_d2) ) //along the first diagonal
                {
                    BE[i][j] = 0.5*(BE[i+1][j+1]+BE[i-1][j-1]) + 0.5*lap_d1;
                }
                else if(fabs(grad_d1)+ fabs(lap_d1) > fabs(grad_d2)+ fabs(lap_d2) ) //along the second diagonal
                {
                    BE[i][j] = 0.5*(BE[i-1][j+1]+BE[i+1][j-1]) + 0.5*lap_d2;
                }
                else
                {
                    BE[i][j] = 0.25*(BE[i+1][j+1]+BE[i-1][j+1]+BE[i+1][j-1]+BE[i-1][j-1]);

                    //Adds the laplacian of green
                    BE[i][j] += 0.25*(4*GE[i][j] - (GE[i+1][j+1]+GE[i-1][j+1]+GE[i+1][j-1]+GE[i-1][j-1]) );
                }
            }
            else if(i%2 == bi && j%2 == bj) //Blue pixels
            {
                grad_d1 = RE[i+1][j+1]-RE[i-1][j-1];
                lap_d1 = 2*GE[i][j] - GE[i+1][j+1] -GE[i-1][j-1];
                grad_d2 = RE[i-1][j+1]-RE[i+1][j-1];
                lap_d2 = 2*GE[i][j] - GE[i-1][j+1] -GE[i+1][j-1];

                //Interpolation of red
                if(fabs(grad_d1)+ fabs(lap_d1) < fabs(grad_d2)+ fabs(lap_d2) ) //along the first diagonal
                {
                    RE[i][j] = 0.5*(RE[i+1][j+1]+RE[i-1][j-1]) + 0.5*lap_d1;
                }
                else if(fabs(grad_d1)+ fabs(lap_d1) > fabs(grad_d2)+ fabs(lap_d2) ) //along the second diagonal
                {
                    RE[i][j] = 0.5*(RE[i-1][j+1]+RE[i+1][j-1]) + 0.5*lap_d2;
                }
                else
                {
                    RE[i][j] = 0.25*(RE[i+1][j+1]+RE[i-1][j+1]+RE[i+1][j-1]+RE[i-1][j-1]);

                    //Adds the laplacian of green
                    RE[i][j] += 0.25*(4*GE[i][j] - (GE[i+1][j+1]+GE[i-1][j+1]+GE[i+1][j-1]+GE[i-1][j-1]) );
                }
            }
            else //Green pixels
            {
                if(i%2 == ri) // red line
                {
                    //red interpolation along the line
                    RE[i][j] = 0.5*(RE[i][j+1]+RE[i][j-1]);
                    //adds the laplacian of green
                    RE[i][j] += 0.5*(2*GE[i][j] - GE[i][j+1] -GE[i][j-1]);

                    //blue interpolation along the column
                    BE[i][j] = 0.5*(BE[i+1][j]+BE[i-1][j]);
                    //adds the laplacian of green
                    BE[i][j] += 0.5*(2*GE[i][j] - GE[i+1][j] -GE[i-1][j]);
                }
                else // blue line
                {
                    //red interpolation along the column
                    RE[i][j] = 0.5*(RE[i+1][j]+RE[i-1][j]);
                    RE[i][j] += 0.5*(2*GE[i][j] - GE[i+1][j] -GE[i-1][j]);

                    //blue interpolation along the line
                    BE[i][j] = 0.5*(BE[i][j+1]+BE[i][j-1]);
                    BE[i][j] += 0.5*(2*GE[i][j] - GE[i][j+1] -GE[i][j-1]);
                }
            }

            if(RE[i][j] < 0)
                RE[i][j] = 0.;

            if(RE[i][j] > t_Max)
                RE[i][j] = t_Max;

            if(BE[i][j] < 0)
                BE[i][j] = 0.;

            if(BE[i][j] > t_Max)
                BE[i][j] = t_Max;
        }

    RE.ExtractBlock(R,ExtSize,ExtSize,t_NbRows,t_NbCols);
    GE.ExtractBlock(G,ExtSize,ExtSize,t_NbRows,t_NbCols);
    BE.ExtractBlock(B,ExtSize,ExtSize,t_NbRows,t_NbCols);
}


void M3Matrix::HAPrefDebayering(M3Matrix& R,M3Matrix& G,M3Matrix& B, int BayerPattern)
{
    R.Reshape(t_NbRows,t_NbCols);
    G.Reshape(t_NbRows,t_NbCols);
    B.Reshape(t_NbRows,t_NbCols);

    int ExtSize = 2; //Must be even !!!
    int NbRowsExt = t_NbRows+2*ExtSize;
    int NbColsExt = t_NbCols+2*ExtSize;

    M3Matrix RawMirrored;
    (*this).Mirror(RawMirrored,ExtSize,ExtSize);

    M3Matrix RE(NbRowsExt,NbColsExt),GE(NbRowsExt,NbColsExt),BE(NbRowsExt,NbColsExt); // extended color matrices

    int ri=0,rj=0,bi=0,bj=0;
    if(BayerPattern == 0)
    {
        ri = 0;  rj = 0;
    }
    else if(BayerPattern == 1)
    {
        ri = 0;  rj = 1;
    }
    else if(BayerPattern == 2)
    {
        ri = 1;  rj = 0;
    }
    else if(BayerPattern == 3)
    {
        ri = 1;  rj = 1;
    }

    bi = 1-ri;  bj = 1-rj;

    //Fill the colors where available
    for(int i=0;i<NbRowsExt;i++)
        for(int j=0;j<NbColsExt;j++)
        {
            if(i%2 == ri && j%2 == rj) //red pixels
            {
                RE[i][j] = RawMirrored[i][j];
            }
            else if(i%2 == bi && j%2 == bj) //blue pixels
            {
                BE[i][j] = RawMirrored[i][j];
            }
            else //green pixels
            {
                GE[i][j] = RawMirrored[i][j];
            }
        }

    M3Matrix Rh(NbRowsExt,NbColsExt),Rv(NbRowsExt,NbColsExt);
    M3Matrix Gh(NbRowsExt,NbColsExt),Gv(NbRowsExt,NbColsExt);
    M3Matrix Bh(NbRowsExt,NbColsExt),Bv(NbRowsExt,NbColsExt);

    double grad_v,lap_v,grad_h,lap_h,grad_d1,lap_d1,grad_d2,lap_d2;

    //Now begin the interpolation algorithm for the green channel
    for(int i=ExtSize;i<t_NbRows+ExtSize;i++)
        for(int j=ExtSize;j<t_NbCols+ExtSize;j++)
        {
            if(i%2 == ri && j%2 == rj) //Red pixels
            {
                Gv[i][j] = 0.5 * (GE[i-1][j] + GE[i+1][j]) + (2*RE[i][j] - RE[i-2][j] - RE[i+2][j]) * 0.25; //Vertical estimate
                Gh[i][j] = 0.5 * (GE[i][j-1] + GE[i][j+1]) + (2*RE[i][j] - RE[i][j-2] - RE[i][j+2]) * 0.25; //Horizontal estimate

                Rh[i][j] = RE[i][j];
                Rv[i][j] = RE[i][j];
            }
            else if(i%2 == bi && j%2 == bj) //Blue pixels
            {
                Bh[i][j] = BE[i][j];
                Bv[i][j] = BE[i][j];

                Gv[i][j] = 0.5*(GE[i-1][j] + GE[i+1][j]) + 0.25 *  (2*BE[i][j] - BE[i-2][j] - BE[i+2][j]); //Vertical estimate
                Gh[i][j] = 0.5*(GE[i][j-1] + GE[i][j+1]) + 0.25 *  (2*BE[i][j] - BE[i][j-2] - BE[i][j+2]); //Horizontal estimate
            }
            else if(i%2 == ri && j%2 == bj) //Green original pixel, red line, blue column
            {
                // Green original value available
                Gv[i][j] = GE[i][j];
                Gh[i][j] = GE[i][j];

                Rh[i][j] = 0.5*(RE[i][j-1] + RE[i][j+1]) + 0.25 * (2*GE[i][j] - GE[i][j-2] - GE[i][j+2]); //Horizontal red estimate
                Bv[i][j] = 0.5*(BE[i-1][j] + BE[i+1][j]) + 0.25 * (2*GE[i][j] - GE[i-2][j] - GE[i+2][j]); //Vertical blue estimate
            }
            else if(i%2 == bi && j%2 == rj) //Green original pixel, blue line, red column
            {
                // Green original value available
                Gv[i][j] = GE[i][j];
                Gh[i][j] = GE[i][j];

                Rv[i][j] = 0.5*(RE[i-1][j] + RE[i+1][j]) + 0.25 * (2*GE[i][j] - GE[i-2][j] - GE[i+2][j]); //Vertical red estimate
                Bh[i][j] = 0.5*(BE[i][j-1] + BE[i][j+1]) + 0.25 * (2*GE[i][j] - GE[i][j-2] - GE[i][j+2]); //Horizontal estimate
            }
        }

    //Step 2: Compute the differences
    M3Matrix Delta_GR_h,Delta_GR_v,Delta_GB_h,Delta_GB_v;
    Delta_GR_h = Gh - Rh;
    Delta_GR_v = Gv - Rv;
    Delta_GB_h = Gh - Bh;
    Delta_GB_v = Gv - Bv;

    //Step 3: Now smooth the differences
    M3Matrix HF,VF;
    double sigma = 3.6;
    //Horizontal smoothing of the horizontal ones for the differences
    HF.MakeGaussianFilter_1D(sigma,5);
    Delta_GR_h.LinearFilter(HF);
    Delta_GB_h.LinearFilter(HF);

    //Vertical smoothing of the vertical estimate for the differences
    HF.MakeGaussianFilter_1D(sigma,7);
    HF.Transpose(VF);
    Delta_GR_v.LinearFilter(VF);
    Delta_GB_v.LinearFilter(VF);

    //Step 4: Final estimates for initially missing green pixels
    for(int i=ExtSize;i<t_NbRows+ExtSize;i++)
        for(int j=ExtSize;j<t_NbCols+ExtSize;j++)
        {
            if( ( i%2 == ri && j%2 == rj) || (i%2 == bi && j%2 == bj)  ) //Red and blue pixels
            {
                grad_v = GE[i+1][j] - GE[i-1][j]; //Vertical gradient
                grad_h = GE[i][j+1] - GE[i][j-1]; //Horizontal gradient
                if(i%2 == ri && j%2 == rj) //Red pixels
                {
                    lap_v = 2*RE[i][j] - RE[i-2][j] - RE[i+2][j]; //Vertical Laplace
                    grad_h = GE[i][j+1] - GE[i][j-1]; //Horizontal gradient
                    lap_h = 2*RE[i][j] - RE[i][j-2] - RE[i][j+2]; //Horizontal Laplace

                    //This is similar to the Hamilton-Adams test
                    if(fabs(grad_v) +  fabs(lap_v)  <   (fabs(grad_h) +  fabs(lap_h)) )
                    {
                        GE[i][j] = Rv[i][j] + Delta_GR_v[i][j];
                    }
                    else if( fabs(grad_h) + fabs(lap_h)   < ( fabs(grad_v) +  fabs(lap_v) ) )
                    {
                        GE[i][j] = Rh[i][j] + Delta_GR_h[i][j];
                    }
                    else
                    {
                        GE[i][j] = 0.5*(Gv[i][j] + Gh[i][j]);
                    }
                }
                else if(i%2 == bi && j%2 == bj) //Blue pixels
                {
                    lap_v = 2*BE[i][j] - BE[i-2][j] - BE[i+2][j]; //Vertical Laplace
                    lap_h = 2*BE[i][j] - BE[i][j-2] - BE[i][j+2]; //Horizontal Laplace

                    //This is similar to the Hamilton-Adams test
                    if(fabs(grad_v) + fabs(lap_v)  < fabs(grad_h) + fabs(lap_h))
                    {
                        GE[i][j] = Bv[i][j] + Delta_GB_v[i][j];
                    }
                    else if(fabs(grad_h) + fabs(lap_h) < fabs(grad_v) +  fabs(lap_v)   )
                    {
                        GE[i][j] = Bh[i][j] + Delta_GB_h[i][j];
                    }
                    else
                    {
                        GE[i][j] = 0.5*(Gv[i][j] + Gh[i][j]);
                    }
                }
            }

            if(GE[i][j] < 0)
                GE[i][j] = 0.;

            if(GE[i][j] > t_Max)
                GE[i][j] = t_Max;
        }


    //Now begin the interpolation algorithm for the blue and red channel
    for(int i=ExtSize;i<t_NbRows+ExtSize;i++)
        for(int j=ExtSize;j<t_NbCols+ExtSize;j++)
        {
            if(i%2 == ri && j%2 == rj) //Red pixels
            {
                grad_d1 = BE[i+1][j+1]-BE[i-1][j-1];
                lap_d1 = 2*GE[i][j] - GE[i+1][j+1] -GE[i-1][j-1];
                grad_d2 = BE[i-1][j+1]-BE[i+1][j-1];
                lap_d2 = 2*GE[i][j] - GE[i-1][j+1] -GE[i+1][j-1];

                //Interpolation of blue
                if(fabs(grad_d1)+ fabs(lap_d1) < fabs(grad_d2)+ fabs(lap_d2) ) //along the first diagonal
                {
                    BE[i][j] = 0.5*(BE[i+1][j+1]+BE[i-1][j-1]) + 0.5*lap_d1;
                }
                else if(fabs(grad_d1)+ fabs(lap_d1) > fabs(grad_d2)+ fabs(lap_d2) ) //along the second diagonal
                {
                    BE[i][j] = 0.5*(BE[i-1][j+1]+BE[i+1][j-1]) + 0.5*lap_d2;
                }
                else
                {
                    BE[i][j] = 0.25*(BE[i+1][j+1]+BE[i-1][j+1]+BE[i+1][j-1]+BE[i-1][j-1]);

                    //Adds the laplacian of green
                    BE[i][j] += 0.25*(4*GE[i][j] - (GE[i+1][j+1]+GE[i-1][j+1]+GE[i+1][j-1]+GE[i-1][j-1]) );
                }
            }
            else if(i%2 == bi && j%2 == bj) //Blue pixels
            {
                grad_d1 = RE[i+1][j+1]-RE[i-1][j-1];
                lap_d1 = 2*GE[i][j] - GE[i+1][j+1] -GE[i-1][j-1];
                grad_d2 = RE[i-1][j+1]-RE[i+1][j-1];
                lap_d2 = 2*GE[i][j] - GE[i-1][j+1] -GE[i+1][j-1];

                //Interpolation of red
                if(fabs(grad_d1)+ fabs(lap_d1) < fabs(grad_d2)+ fabs(lap_d2) ) //along the first diagonal
                {
                    RE[i][j] = 0.5*(RE[i+1][j+1]+RE[i-1][j-1]) + 0.5*lap_d1;
                }
                else if(fabs(grad_d1)+ fabs(lap_d1) > fabs(grad_d2)+ fabs(lap_d2) ) //along the second diagonal
                {
                    RE[i][j] = 0.5*(RE[i-1][j+1]+RE[i+1][j-1]) + 0.5*lap_d2;
                }
                else
                {
                    RE[i][j] = 0.25*(RE[i+1][j+1]+RE[i-1][j+1]+RE[i+1][j-1]+RE[i-1][j-1]);

                    //Adds the laplacian of green
                    RE[i][j] += 0.25*(4*GE[i][j] - (GE[i+1][j+1]+GE[i-1][j+1]+GE[i+1][j-1]+GE[i-1][j-1]) );
                }
            }
            else //Green pixels
            {
                if(i%2 == ri) // red line
                {
                    //red interpolation along the line
                    RE[i][j] = 0.5*(RE[i][j+1]+RE[i][j-1]);
                    //adds the laplacian of green
                    RE[i][j] += 0.5*(2*GE[i][j] - GE[i][j+1] -GE[i][j-1]);

                    //blue interpolation along the column
                    BE[i][j] = 0.5*(BE[i+1][j]+BE[i-1][j]);
                    //adds the laplacian of green
                    BE[i][j] += 0.5*(2*GE[i][j] - GE[i+1][j] -GE[i-1][j]);
                }
                else // blue line
                {
                    //red interpolation along the column
                    RE[i][j] = 0.5*(RE[i+1][j]+RE[i-1][j]);
                    RE[i][j] += 0.5*(2*GE[i][j] - GE[i+1][j] -GE[i-1][j]);

                    //blue interpolation along the line
                    BE[i][j] = 0.5*(BE[i][j+1]+BE[i][j-1]);
                    BE[i][j] += 0.5*(2*GE[i][j] - GE[i][j+1] -GE[i][j-1]);
                }
            }

            if(RE[i][j] < 0)
                RE[i][j] = 0.;

            if(RE[i][j] > t_Max)
                RE[i][j] = t_Max;

            if(BE[i][j] < 0)
                BE[i][j] = 0.;

            if(BE[i][j] > t_Max)
                BE[i][j] = t_Max;
        }

    RE.ExtractBlock(R,ExtSize,ExtSize,t_NbRows,t_NbCols);
    GE.ExtractBlock(G,ExtSize,ExtSize,t_NbRows,t_NbCols);
    BE.ExtractBlock(B,ExtSize,ExtSize,t_NbRows,t_NbCols);
}


void M3Matrix::HAJapanDebayering(M3Matrix& R,M3Matrix& G,M3Matrix& B, int BayerPattern, int Threshold, int epsilon)
{
    cout << "Threshold : " << Threshold << endl;
    int ThresholdR = 0;
    int ThresholdB = 0;

    R.Reshape(t_NbRows,t_NbCols);
    G.Reshape(t_NbRows,t_NbCols);
    B.Reshape(t_NbRows,t_NbCols);

    int ExtSize = 2; //Must be even !!!
    int NbRowsExt = t_NbRows+2*ExtSize;
    int NbColsExt = t_NbCols+2*ExtSize;

    M3Matrix RawMirrored;
    (*this).Mirror(RawMirrored,ExtSize,ExtSize);

    M3Matrix GEJ, RE,GE,BE; // extended color matrices
    RE.Reshape(NbRowsExt,NbColsExt);
    GE.Reshape(NbRowsExt,NbColsExt);
    GEJ.Reshape(NbRowsExt,NbColsExt);
    BE.Reshape(NbRowsExt,NbColsExt);

    int ri=0,rj=0,bi=0,bj=0;
    if(BayerPattern == 0)
    {
        ri = 0;  rj = 0;
    }
    else if(BayerPattern == 1)
    {
        ri = 0;  rj = 1;
    }
    else if(BayerPattern == 2)
    {
        ri = 1;  rj = 0;
    }
    else if(BayerPattern == 3)
    {
        ri = 1;  rj = 1;
    }

    bi = 1-ri;  bj = 1-rj;

    //Fill the colors where available
    for(int i=0;i<NbRowsExt;i++)
        for(int j=0;j<NbColsExt;j++)
        {
            if(i%2 == ri && j%2 == rj) //red pixels
            {
                RE[i][j] = RawMirrored[i][j];
            }
            else if(i%2 == bi && j%2 == bj) //blue pixels
            {
                BE[i][j] = RawMirrored[i][j];
            }
            else //green pixels
            {
                GE[i][j] = RawMirrored[i][j];
            }
        }

    //Japanese addition
    for(int i=ExtSize;i<t_NbRows+ExtSize;i++)
        for(int j=ExtSize;j<t_NbCols+ExtSize;j++)
        {
            if( (i%2 == ri && j%2 == bj) || (i%2 == bi && j%2 == rj)) //green pixels
            {
                int sum_of_weights = 4;
                double sum_val  = 4*RawMirrored[i][j];

                if(fabs(RawMirrored[i][j] - RawMirrored[i-1][j-1] ) < epsilon)
                {
                    sum_val += RawMirrored[i-1][j-1];
                    sum_of_weights++;
                }

                if(fabs(RawMirrored[i][j] - RawMirrored[i-1][j+1] ) < epsilon)
                {
                    sum_val += RawMirrored[i-1][j+1];
                    sum_of_weights++;
                }

                if(fabs(RawMirrored[i][j] - RawMirrored[i+1][j+1] ) < epsilon)
                {
                    sum_val += RawMirrored[i+1][j+1];
                    sum_of_weights++;
                }

                if(fabs(RawMirrored[i][j] - RawMirrored[i+1][j-1] ) < epsilon)
                {
                    sum_val += RawMirrored[i+1][j-1];
                    sum_of_weights++;
                }

                GEJ[i][j] = sum_val/sum_of_weights;
            }
        }

    GE = GEJ; //too slow ?

    //Now begin the interpolation algorithm for the green channel
    double Lh,Lv;

    //double alphaR = 1.0, alphaG = 1.0, alphaB = 1.0;
    for(int i=ExtSize;i<t_NbRows+ExtSize;i++)
        for(int j=ExtSize;j<t_NbCols+ExtSize;j++)
        {
            if(i%2 == ri && j%2 == rj) //Red pixels
            {
                Lh = fabs(GE[i][j+1]-GE[i][j-1]) + fabs(2*RE[i][j] - RE[i][j+2] -RE[i][j-2]);
                Lv = fabs(GE[i+1][j]-GE[i-1][j]) + fabs(2*RE[i][j] - RE[i+2][j] -RE[i-2][j]);

                if(Lh < Lv - Threshold)
                    GE[i][j] = 0.5*(GE[i][j+1] + GE[i][j-1]) + 0.25*(2*RE[i][j] - RE[i][j+2] -RE[i][j-2]);
                else if(Lv < Lh - Threshold)
                    GE[i][j] = 0.5*(GE[i+1][j] + GE[i-1][j]) + 0.25*(2*RE[i][j] - RE[i+2][j]  -RE[i-2][j]);
                else
                    GE[i][j] = 0.25*(GE[i][j+1] + GE[i][j-1] + GE[i+1][j] + GE[i-1][j]) + 0.125*(4*RE[i][j] - RE[i][j+2] -RE[i][j-2] - RE[i+2][j] -RE[i-2][j]);
            }
            else if(i%2 == bi && j%2 == bj) //Blue pixels
            {
                Lh = fabs(GE[i][j+1]-GE[i][j-1]) + fabs(2*BE[i][j] - BE[i][j+2] -BE[i][j-2]);
                Lv = fabs(GE[i+1][j]-GE[i-1][j]) + fabs(2*BE[i][j] - BE[i+2][j] -BE[i-2][j]);

                if(Lh < Lv - Threshold)
                    GE[i][j] = 0.5*(GE[i][j+1] + GE[i][j-1]) + 0.25*(2*BE[i][j] - BE[i][j+2] - BE[i][j-2]);
                else if(Lv < Lh - Threshold)
                    GE[i][j] = 0.5*(GE[i+1][j] + GE[i-1][j]) + 0.25*(2*BE[i][j] - BE[i+2][j] -BE[i-2][j]);
                else
                    GE[i][j] = 0.25*(GE[i][j+1] + GE[i][j-1] + GE[i+1][j] + GE[i-1][j]) + 0.125*(4*BE[i][j] - BE[i][j+2] -BE[i][j-2] - BE[i+2][j] -BE[i-2][j]);
            }

            if(GE[i][j]<0)
                GE[i][j] = 0.;
            if(GE[i][j]>t_Max)
                GE[i][j] = t_Max;
        }


    //Now interpolation of red and blue value at green pixels
    double delta_n,delta_p; // this is for the diagonal gradients
    for(int i=ExtSize;i<t_NbRows+ExtSize;i++)
        for(int j=ExtSize;j<t_NbCols+ExtSize;j++)
        {
            if(i%2 == ri && j%2 == rj) //red pixel
            {
                //Estimates diagonal gradients
                delta_n = fabs(BE[i-1][j-1]-BE[i+1][j+1]) + fabs(2*GE[i][j] - GE[i-1][j-1] - GE[i+1][j+1]);
                delta_p = fabs(BE[i+1][j-1]-BE[i-1][j+1]) + fabs(2*GE[i][j] - GE[i+1][j-1] - GE[i-1][j+1]);

                //Estimates blue value
                if(delta_n < delta_p -ThresholdB)
                    BE[i][j] = 0.5*(BE[i-1][j-1] + BE[i+1][j+1]) + 0.5*(2*GE[i][j] - GE[i-1][j-1] - GE[i+1][j+1]);
                else if(delta_p < delta_n - ThresholdB)
                    BE[i][j] = 0.5*(BE[i+1][j-1] + BE[i-1][j+1]) + 0.5*(2*GE[i][j] - GE[i+1][j-1] - GE[i-1][j+1]);
                else
                    BE[i][j] = 0.25*(BE[i-1][j-1] + BE[i+1][j+1] + BE[i+1][j-1] + BE[i-1][j+1]) + 0.25*(4*GE[i][j] - GE[i+1][j-1] - GE[i-1][j+1] - GE[i-1][j-1] - GE[i+1][j+1]);
            }
            else if(i%2 == bi && j%2 == bj) //blue pixel
            {
                //Estimates diagonal gradients
                delta_n = fabs(RE[i-1][j-1]-RE[i+1][j+1]) + fabs(2*GE[i][j] - GE[i-1][j-1] - GE[i+1][j+1]);
                delta_p = fabs(RE[i+1][j-1]-RE[i-1][j+1]) + fabs(2*GE[i][j] - GE[i+1][j-1] - GE[i-1][j+1]);

                //Estimates red value
                if(delta_n < delta_p -ThresholdR)
                    RE[i][j] = 0.5*(RE[i-1][j-1] + RE[i+1][j+1]) + 0.25*(2*GE[i][j] - GE[i-1][j-1] - GE[i+1][j+1]);
                else if(delta_p < delta_n -ThresholdR)
                    RE[i][j] = 0.5*(RE[i+1][j-1] + RE[i-1][j+1]) + 0.25*(2*GE[i][j] - GE[i+1][j-1] - GE[i-1][j+1]);
                else
                    RE[i][j] = 0.25*(RE[i-1][j-1] + RE[i+1][j+1] + RE[i+1][j-1] + RE[i-1][j+1]) + 0.25*(4*GE[i][j] - GE[i+1][j-1] - GE[i-1][j+1] - GE[i-1][j-1] - GE[i+1][j+1]);
            }
            else if(i%2 == ri && j%2 == bj) //green pixel, red line, blue column
            {
                //Estimates red value
                RE[i][j] = 0.5*(RE[i][j+1] + RE[i][j-1]) + 0.5*(2*GE[i][j] - GE[i][j+1] - GE[i][j-1]);

                //Estimates blue value
                BE[i][j] = 0.5*(BE[i+1][j] + BE[i-1][j]) + 0.5*(2*GE[i][j] - GE[i+1][j] -  GE[i-1][j]);
            }
            else if(i%2 == bi && j%2 == rj) //green pixel, blue line, red column
            {
                //Estimates red value
                RE[i][j] = 0.5*(RE[i+1][j] + RE[i-1][j]) + 0.5*(2*GE[i][j] - GE[i+1][j] -  GE[i-1][j]);

                //Estimates blue value
                BE[i][j] = 0.5*(BE[i][j+1] + BE[i][j-1]) + 0.5*(2*GE[i][j] - GE[i][j+1] -  GE[i][j-1]);
            }

            if(BE[i][j]<0)
                BE[i][j] = 0.;
            if(BE[i][j]>t_Max)
                BE[i][j] = t_Max;

            if(RE[i][j]<0)
                RE[i][j] = 0.;
            if(RE[i][j]>t_Max)
                RE[i][j] = t_Max;
        }

    RE.ExtractBlock(R,ExtSize,ExtSize,t_NbRows,t_NbCols);
    GE.ExtractBlock(G,ExtSize,ExtSize,t_NbRows,t_NbCols);
    BE.ExtractBlock(B,ExtSize,ExtSize,t_NbRows,t_NbCols);
}


void M3Matrix::ZhangWuDebayering(M3Matrix& R,M3Matrix& G,M3Matrix& B, int BayerPattern)
{
    R.Reshape(t_NbRows,t_NbCols);
    G.Reshape(t_NbRows,t_NbCols);
    B.Reshape(t_NbRows,t_NbCols);

    int ExtSize = 2; //Must be even !!!
    M3Matrix RawMirrored;
    (*this).Mirror(RawMirrored,ExtSize,ExtSize);

    M3Matrix RE,GE,BE; // extended color matrices
    RE.Reshape(t_NbRows+2*ExtSize,t_NbCols+2*ExtSize);
    GE.Reshape(t_NbRows+2*ExtSize,t_NbCols+2*ExtSize);
    BE.Reshape(t_NbRows+2*ExtSize,t_NbCols+2*ExtSize);

    int ri=0,rj=0,bi,bj;
    if(BayerPattern == 0)
    {
        ri = 0;    rj = 0;
    }
    else if(BayerPattern == 1)
    {
        ri = 0;    rj = 1;
    }
    else if(BayerPattern == 2)
    {
        ri = 1;    rj = 0;
    }
    else if(BayerPattern == 3)
    {
        ri = 1;    rj = 1;
    }

    bi = 1-ri;  bj = 1-rj;


    //Fill the colors where available
    for(int i=0;i<t_NbRows+2*ExtSize;i++)
        for(int j=0;j<t_NbCols+2*ExtSize;j++)
        {
            if(i%2 == ri && j%2 == rj) //red pixels
            {
                RE[i][j] = RawMirrored[i][j];
            }
            else if(i%2 == bi && j%2 == bj) //blue pixels
            {
                BE[i][j] = RawMirrored[i][j];
            }
            else //green pixels
            {
                GE[i][j] = RawMirrored[i][j];
            }
        }

    M3Matrix Gh(t_NbRows+2*ExtSize,t_NbCols+2*ExtSize),Gv(t_NbRows+2*ExtSize,t_NbCols+2*ExtSize);
    M3Matrix Rh(t_NbRows+2*ExtSize,t_NbCols+2*ExtSize),Rv(t_NbRows+2*ExtSize,t_NbCols+2*ExtSize);
    M3Matrix Bh(t_NbRows+2*ExtSize,t_NbCols+2*ExtSize),Bv(t_NbRows+2*ExtSize,t_NbCols+2*ExtSize);

    double grad_v, lap_v, grad_h, lap_h, grad_d1, lap_d1, grad_d2, lap_d2;
    double fac = 1.0;

    //Now begin the interpolation algorithm for the green channel
    for(int i=ExtSize;i<t_NbRows+ExtSize;i++)
        for(int j=ExtSize;j<t_NbCols+ExtSize;j++)
        {
            if(i%2 == ri && j%2 == rj) //Red pixels
            {
                Gv[i][j] = 0.5 * (GE[i-1][j] + GE[i+1][j]) + (2*RE[i][j] - RE[i-2][j] - RE[i+2][j]) * 0.25; //Vertical estimate
                Gh[i][j] = 0.5 * (GE[i][j-1] + GE[i][j+1]) + (2*RE[i][j] - RE[i][j-2] - RE[i][j+2]) * 0.25; //Horizontal estimate

                Rh[i][j] = RE[i][j];
                Rv[i][j] = RE[i][j];
            }
            else if(i%2 == bi && j%2 == bj) //Blue pixels
            {
                Bh[i][j] = BE[i][j];
                Bv[i][j] = BE[i][j];

                Gv[i][j] = 0.5*(GE[i-1][j] + GE[i+1][j]) + 0.25 *  (2*BE[i][j] - BE[i-2][j] - BE[i+2][j]); //Horizontal estimate
                Gh[i][j] = 0.5*(GE[i][j-1] + GE[i][j+1]) + 0.25 *  (2*BE[i][j] - BE[i][j-2] - BE[i][j+2]); //Horizontal estimate
            }
            else //red and blue estimates at green pixels
            {
                Gv[i][j] = GE[i][j];
                Gh[i][j] = GE[i][j];

                Rv[i][j] = 0.5*(RE[i-1][j] + RE[i+1][j]) + 0.25 * (2*GE[i][j] - GE[i-2][j] - GE[i+2][j]); //Horizontal estimate
                Rh[i][j] = 0.5*(RE[i][j-1] + RE[i][j+1]) + 0.25 * (2*GE[i][j] - GE[i][j-2] - GE[i][j+2]); //Vertical estimate

                Bv[i][j] = 0.5*(BE[i-1][j] + BE[i+1][j]) + 0.25 * (2*GE[i][j] - GE[i-2][j] - GE[i+2][j]); //Horizontal estimate
                Bh[i][j] = 0.5*(BE[i][j-1] + BE[i][j+1]) + 0.25 * (2*GE[i][j] - GE[i][j-2] - GE[i][j+2]); //Horizontal estimate
            }
        }

    M3Matrix Delta_GR_h,Delta_GR_v,Delta_GB_h,Delta_GB_v;
    Delta_GR_h = Gh - Rh;
    Delta_GR_v = Gv - Rv;
    Delta_GB_h = Gh - Bh;
    Delta_GB_v = Gv - Bv;

    //Now smooth the differences
    M3Matrix HF,VF;
    double sigma = 3.5;

    //Horizontal smoothing of the horizontal ones for the differences
    HF.MakeGaussianFilter_1D(sigma,7);
    Delta_GR_h.LinearFilter(HF);
    Delta_GB_h.LinearFilter(HF);

    //Vertical smoothing of the vertical estimate for the differences
    HF.MakeGaussianFilter_1D(sigma,7); HF.Transpose(VF);
    Delta_GR_v.LinearFilter(VF);
    Delta_GB_v.LinearFilter(VF);

    //TODO: add the "fine features of Zhang-Wu here" << endl;

    for(int i=ExtSize;i<t_NbRows+ExtSize;i++)
        for(int j=ExtSize;j<t_NbCols+ExtSize;j++)
        {
            if(i%2 == ri && j%2 == rj) //Red pixels
            {
                grad_v = GE[i+1][j] - GE[i-1][j]; //Vertical gradient
                lap_v = 2*RE[i][j] - RE[i-2][j] - RE[i+2][j]; //Vertical Laplace
                grad_h = GE[i][j+1] - GE[i][j-1]; //Horizontal gradient
                lap_h = 2*RE[i][j] - RE[i][j-2] - RE[i][j+2]; //Horizontal Laplace

                //This corresponds to the test used in pure Hamilton-Adams debayering
                if(fabs(grad_v) +  fabs(lap_v)  <   (fabs(grad_h) +  fabs(lap_h)) )
                {
                    GE[i][j] = Rv[i][j] + Delta_GR_v[i][j];
                }
                else if( fabs(grad_h) + fabs(lap_h)   < ( fabs(grad_v) +  fabs(lap_v) ) )
                {
                    GE[i][j] = Rh[i][j] + Delta_GR_h[i][j];
                }
                else
                {
                    GE[i][j] = 0.5*(Gv[i][j] + Gh[i][j]);
                }
            }
            else if(i%2 == bi && j%2 == bj) //Blue pixels
            {
                grad_v = GE[i+1][j] - GE[i-1][j]; //Vertical gradient
                lap_v = 2*BE[i][j] - BE[i-2][j] - BE[i+2][j]; //Vertical Laplace
                grad_h = GE[i][j+1] - GE[i][j-1]; //Horizontal gradient
                lap_h = 2*BE[i][j] - BE[i][j-2] - BE[i][j+2]; //Horizontal Laplace

                //Temporary: this is the Hamilton-Adams test
                if(fabs(grad_v) + fabs(lap_v)  < fabs(grad_h) + fabs(lap_h))
                {
                    //GE[i][j] = Gv[i][j];
                    GE[i][j] = Bv[i][j] + Delta_GB_v[i][j];
                }
                else if(fabs(grad_h) + fabs(lap_h) < fabs(grad_v) +  fabs(lap_v)   )
                {
                    //GE[i][j] = Gh[i][j];
                    GE[i][j] = Bh[i][j] + Delta_GB_h[i][j];
                }
                else
                {
                    GE[i][j] = 0.5*(Gv[i][j] + Gh[i][j]);
                }
            }

            if(GE[i][j] < 0)
                GE[i][j] = 0.;
            if(GE[i][j] > t_Max)
                GE[i][j] = t_Max;

        }

    //Now begin the interpolation algorithm for the blue and red channel
    for(int i=ExtSize;i<t_NbRows+ExtSize;i++)
        for(int j=ExtSize;j<t_NbCols+ExtSize;j++)
        {
            if(i%2 == ri && j%2 == rj) //Red pixels
            {
                grad_d1 = BE[i+1][j+1]-BE[i-1][j-1];
                lap_d1 = 2*GE[i][j] - GE[i+1][j+1] -GE[i-1][j-1];
                grad_d2 = BE[i-1][j+1]-BE[i+1][j-1];
                lap_d2 = 2*GE[i][j] - GE[i-1][j+1] -GE[i+1][j-1];

                //Interpolation of blue
                if(fabs(grad_d1)+ fabs(lap_d1) < fabs(grad_d2)+ fabs(lap_d2) ) //along the first diagonal
                {
                    BE[i][j] = 0.5*(BE[i+1][j+1]+BE[i-1][j-1]) + 0.5*fac*lap_d1;
                }
                else if(fabs(grad_d1)+ fabs(lap_d1) > fabs(grad_d2)+ fabs(lap_d2) ) //along the second diagonal
                {
                    BE[i][j] = 0.5*(BE[i-1][j+1]+BE[i+1][j-1]) + 0.5*fac*lap_d2;
                }
                else
                {
                    BE[i][j] = 0.25*(BE[i+1][j+1]+BE[i-1][j+1]+BE[i+1][j-1]+BE[i-1][j-1]);

                    //Adds the laplacian of green
                    BE[i][j] += 0.25*fac*(4*GE[i][j] - (GE[i+1][j+1]+GE[i-1][j+1]+GE[i+1][j-1]+GE[i-1][j-1]) );
                }
            }
            else if(i%2 == bi && j%2 == bj) //Blue pixels
            {
                grad_d1 = RE[i+1][j+1]-RE[i-1][j-1];
                lap_d1 = 2*GE[i][j] - GE[i+1][j+1] -GE[i-1][j-1];
                grad_d2 = RE[i-1][j+1]-RE[i+1][j-1];
                lap_d2 = 2*GE[i][j] - GE[i-1][j+1] -GE[i+1][j-1];

                //Interpolation of red
                if(fabs(grad_d1)+ fabs(lap_d1) < fabs(grad_d2)+ fabs(lap_d2) ) //along the first diagonal
                {
                    RE[i][j] = 0.5*(RE[i+1][j+1]+RE[i-1][j-1]) + 0.5*lap_d1;
                }
                else if(fabs(grad_d1)+ fabs(lap_d1) > fabs(grad_d2)+ fabs(lap_d2) ) //along the second diagonal
                {
                    RE[i][j] = 0.5*(RE[i-1][j+1]+RE[i+1][j-1]) + 0.5*lap_d2;
                }
                else
                {
                    RE[i][j] = 0.25*(RE[i+1][j+1]+RE[i-1][j+1]+RE[i+1][j-1]+RE[i-1][j-1]);

                    //Adds the laplacian of green
                    RE[i][j] += 0.25*(4*GE[i][j] - (GE[i+1][j+1]+GE[i-1][j+1]+GE[i+1][j-1]+GE[i-1][j-1]) );
                }
            }
            else //Green pixels
            {
                if(i%2 == ri) // red line
                {
                    //red interpolation along the line
                    RE[i][j] = 0.5*(RE[i][j+1]+RE[i][j-1]);
                    //adds the laplacian of green
                    RE[i][j] += 0.5*(2*GE[i][j] - GE[i][j+1] -GE[i][j-1]);

                    //blue interpolation along the column
                    BE[i][j] = 0.5*(BE[i+1][j]+BE[i-1][j]);
                    //adds the laplacian of green
                    BE[i][j] += 0.5*fac*(2*GE[i][j] - GE[i+1][j] -GE[i-1][j]);
                }
                else // blue line
                {
                    //red interpolation along the column
                    RE[i][j] = 0.5*(RE[i+1][j]+RE[i-1][j]);
                    RE[i][j] += 0.5*(2*GE[i][j] - GE[i+1][j] -GE[i-1][j]);

                    //blue interpolation along the line
                    BE[i][j] = 0.5*(BE[i][j+1]+BE[i][j-1]);
                    BE[i][j] += 0.5*fac*(2*GE[i][j] - GE[i][j+1] -GE[i][j-1]);
                }
            }

            if(BE[i][j] < 0)
                BE[i][j] = 0.;
            if(BE[i][j] > t_Max)
                BE[i][j] = t_Max;
            if(RE[i][j] < 0)
                RE[i][j] = 0.;
            if(RE[i][j] > t_Max)
                RE[i][j] = t_Max;
        }

    RE.ExtractBlock(R,ExtSize,ExtSize,t_NbRows,t_NbCols);
    GE.ExtractBlock(G,ExtSize,ExtSize,t_NbRows,t_NbCols);
    BE.ExtractBlock(B,ExtSize,ExtSize,t_NbRows,t_NbCols);
}


//Defect Pixel Correction
//option 0 corresponds to a median correction; 1: to a mean one
void M3Matrix::DPC(M3Matrix& Corrected, M3Matrix& Mask, double threshold, int BayerPattern, int opt)
{
    Mask.Reshape((*this).GetNbRows(),(*this).GetNbCols(),0.);
    Corrected = (*this);

    if(opt == 0)
        cout << "DPC: Median Modus" << endl;
    else if(opt == 1)
        cout << "DPC: Mean Modus" << endl;
    else if(opt == 2)
        cout << "DPC: MinMax Modus" << endl;
    else if(opt == 3)
        cout << "DPC: Directed Modus" << endl;

    int NbDefectPixels = 0, ri=0,rj=0,bi=0,bj=0;
    if(BayerPattern == 0)
    {
        ri = 0;
        rj = 0;
    }
    else if(BayerPattern == 1)
    {
        ri = 0;
        rj = 1;
    }
    else if(BayerPattern == 2)
    {
        ri = 1;
        rj = 0;
    }
    else if(BayerPattern == 3)
    {
        ri = 1;
        rj = 1;
    }

    bi = 1-ri;
    bj = 1-rj;


    for(int i= 2;i<t_NbRows-2;i++)
        for(int j= 2;j<t_NbCols-2;j++)
        {
            double* y = new double[9];

            if((i%2 == ri && j%2 == rj) || (i%2 == bi && j%2 == bj)) //red pixel (first two conditions) or blue pixel (last two conditions)
            {
                for(int m=-1;m<=1;m++)
                    for(int n=-1;n<=1;n++)
                        y[3*(m+1)+n+1] = (*this)[i+2*m][j+2*n];
            }
            else //green pixel: to be implemented
            {
                //quincunx
                y[0] = (*this)[i-2][j];
                y[1] = (*this)[i-1][j-1];
                y[2] = (*this)[i-1][j+1];
                y[3] = (*this)[i][j-2];
                y[4] = (*this)[i][j];
                y[5] = (*this)[i][j+2];
                y[6] = (*this)[i+1][j-1];
                y[7] = (*this)[i+1][j+1];
                y[8] = (*this)[i+2][j];
            }


            qsort(y, 9, sizeof(double), my_cmp);

            if((*this)[i][j] == y[8]) //max value
            {
                if(fabs((*this)[i][j] - y[7]) > threshold)
                {
                    if(opt == 0)
                        Corrected[i][j] = y[5];//5: Median
                    else if(opt == 1)
                    {
                        double val;
                        if((i%2 == ri && j%2 == rj) || (i%2 == bi && j%2 == bj)) //red pixel (first two conditions) or blue pixel (last two conditions)
                        {
                            //Mean of the 4 neighbours (same colour)
                            val = ( (*this)[i-2][j] + (*this)[i+2][j] + (*this)[i][j-2] + (*this)[i][j+2]) / 4.0;
                        }
                        else
                        {
                            val = ( (*this)[i+1][j+1] + (*this)[i-1][j+1] + (*this)[i-1][j-1] + (*this)[i+1][j-1])/4.0;
                        }
                        Corrected[i][j] = val;
                    }
                    else if(opt == 2)
                        Corrected[i][j] = y[7];
                    else if(opt == 3)
                    {
                        double val;
                        if((i%2 == ri && j%2 == rj) || (i%2 == bi && j%2 == bj)) //red pixel (first two conditions) or blue pixel (last two conditions)
                        {
                            //Mean of the 4 neighbours (same colour)
                            val = ( (*this)[i-2][j] + (*this)[i+2][j] + (*this)[i][j-2] + (*this)[i][j+2]) / 4.0;
                        }
                        else
                        {
                            val = ( (*this)[i+1][j+1] + (*this)[i-1][j+1] + (*this)[i-1][j-1] + (*this)[i+1][j-1])/4.0;
                        }
                        Corrected[i][j] = val;
                    }

                    NbDefectPixels++;
                    Mask[i][j] = 255;
                }
            }

            delete[] y;
        }

    cout << NbDefectPixels << " defect pixels detected." << endl;
}


M3Matrix M3Matrix::SimulateRawBayer(int BayerPattern)
{
    M3Matrix Raw,R,G,B;
    Raw.InitialiseDynamic((*this));

    (*this).ExtractRGBChannels(R,G,B);
    Raw.SimulateRawBayer(R,G,B,BayerPattern);
    return Raw;
}


void M3Matrix::SimulateRawBayer(M3Matrix& R,M3Matrix& G,M3Matrix& B, int BayerPattern)
{
    int NbRows = R.GetNbRows();
    int NbCols = R.GetNbCols();
    Reshape(NbRows,NbCols);

    int ri=0,rj=0,bi=0,bj=0;
    if(BayerPattern == 0)
    {
        ri = 0;
        rj = 0;
    }
    else if(BayerPattern == 1)
    {
        ri = 0;
        rj = 1;
    }
    else if(BayerPattern == 2)
    {
        ri = 1;
        rj = 0;
    }
    else if(BayerPattern == 3)
    {
        ri = 1;
        rj = 1;
    }

    bi = 1-ri;
    bj = 1-rj;

    for(int i=0;i<NbRows;i++)
        for(int j=0;j<NbCols;j++)
        {
            if(i%2 == ri && j%2 == rj) //red pixels
            {
                (*this)[i][j] = R[i][j];
            }
            else if(i%2 == bi && j%2 == bj) //blue pixels
            {
                (*this)[i][j] = B[i][j];
            }
            else //green pixels
            {
                (*this)[i][j] = G[i][j];
            }
        }
}


//This is an implementation of the debayering via bilinear interpolation
void M3Matrix::Debayering(M3Matrix& Image_Debayered)
{
    M3Matrix R_Org,G_Org,B_Org;
    M3Matrix R,G,B;
    M3Matrix RM,GM,BM;

    (*this).ExtractRGBChannels(R_Org,G_Org,B_Org);
    double sumG, sumB, sumR;

    R_Org.Mirror(R,2,2);
    G_Org.Mirror(G,2,2);
    B_Org.Mirror(B,2,2);

    int NbRowsM = R.GetNbRows();
    int NbColsM = R.GetNbCols();

    //Naive interpolation
    for(int i=2;i<NbRowsM-2;i++)
        for(int j=2;j<NbColsM-2;j++)
        {
            if(i%2 == 0 && j%2 == 0) //red pixels
            {
                //Green interpolation
                sumG = 0.;
                sumB = 0.;

                sumG += G[i-1][j];
                sumG += G[i][j-1];
                sumG += G[i+1][j];
                sumG += G[i][j+1];
                G[i][j] = sumG/4.; //interpolation of the missing green values with the mean value of the available pixels

                //Blue interpolation
                sumB += B[i-1][j-1];
                sumB += B[i-1][j+1];
                sumB += B[i+1][j-1];
                sumB += B[i+1][j+1];

                B[i][j] = sumB/4.; //interpolation of the missing green values with the mean value of the
            }
            else if(i%2 == 1 && j%2 == 1) //blue pixels
            {
                //Green interpolation
                sumG = 0.;
                sumR = 0.;

                sumG += G[i-1][j];
                sumG += G[i][j-1];
                sumG += G[i+1][j];
                sumG += G[i][j+1];
                G[i][j] = sumG/4.; //interpolation of the missing green values with the mean value of the available pixels

                //Red interpolation
                sumR += R[i-1][j-1];
                sumR += R[i-1][j+1];
                sumR += R[i+1][j-1];
                sumR += R[i+1][j+1];
                R[i][j] = sumR/4.; //interpolation of the missing green values with the mean value of the
            }
            else //green pixels
            {
                if(i%2==0) //on a "red" line
                {
                    //Red interpolation
                    R[i][j] = 0.5*(R[i][j-1]+R[i][j+1]);

                    //Blue interpolation
                    B[i][j] = 0.5*(B[i-1][j]+B[i+1][j]);
                }
                else //on a "blue" line
                {
                    //Red interpolation
                    R[i][j] = 0.5*(R[i-1][j]+R[i+1][j]);

                    //Blue interpolation
                    B[i][j] = 0.5*(B[i][j-1]+B[i][j+1]);
                }
            }
        }

    M3Matrix ROut,GOut,BOut;
    R.ExtractBlock(ROut,2,2,NbRowsM-4,NbColsM-4);
    G.ExtractBlock(GOut,2,2,NbRowsM-4,NbColsM-4);
    B.ExtractBlock(BOut,2,2,NbRowsM-4,NbColsM-4);

    Image_Debayered.CombineRGBChannels(ROut,GOut,BOut);
}


void M3Matrix::GetRevertRows(M3Matrix B)
{
    (*this).Reshape(B.GetNbRows(),B.GetNbCols());

    for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols;j++)
        {
            (*this)[t_NbRows-i-1][j] = B[i][j];
        }
}


void M3Matrix::GetRevertCols(M3Matrix B)
{
    (*this).Reshape(B.GetNbRows(),B.GetNbCols());

    for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols;j++)
        {
            (*this)[i][t_NbCols-j-1] = B[i][j];
        }
}


void M3Matrix::Transpose()
{
  M3Matrix AT;
  (*this).Transpose(AT);
  (*this) = AT;
}


void M3Matrix::Transpose(M3Matrix& AT)
{    
  AT.Reshape(t_NbCols,t_NbRows);
  AT.InitialiseDynamic(*this);
  for(int i=0;i<t_NbRows;i++)
    for(int j=0;j<t_NbCols;j++)
    {
      AT[j][i] = (*this)[i][j];
    }
}


void M3Matrix::Shift(M3Matrix& ShiftedImage, int tx, int ty)
{
  ShiftedImage.Reshape((*this).GetNbRows(),(*this).GetNbCols());
  ShiftedImage.SetValues(0.);

  for(int i=0;i<t_NbRows;i++)
    for(int j=0;j<t_NbCols;j++)
    {
      if(i-tx >= 0 && i-tx <t_NbRows && j-ty >= 0 && j-ty <t_NbCols)
         ShiftedImage[i - tx][j - ty] = (*this)[i][j];
    }
}


void M3Matrix::Extend(M3Matrix& ExtendedImage, int Height, int Width, int ExtensionType)
{
    int H=Height;
    int W=Width;
    if(ExtensionType == 0) //0: mirroring
    {
        (*this).Mirror(ExtendedImage,Height,Width);
    }
    else if(ExtensionType == 1) //1: zero extension
    {
        ExtendedImage.Reshape(t_NbRows+2*Height,t_NbCols+2*Width);
        ExtendedImage.SetValues(0.);

        //Central part
        for(int i=0;i<t_NbRows;i++)
            for(int j=0;j<t_NbCols;j++)
            {
                ExtendedImage[i+H][j+W] = (*this)[i][j];
            }
    }
}


void M3Matrix::Mirror(M3Matrix& MirroredImage, int Height, int Width)
{
    int H=Height;
    int W=Width;
    int NbRows = GetNbRows();
    int NbCols = GetNbCols();
    MirroredImage.Reshape(NbRows+2*H,NbCols+2*W);

    //central part
    for(int i=0;i<NbRows;i++)
        for(int j=0;j<NbCols;j++)
        {
            MirroredImage[i+H][j+W] = (*this)[i][j];
        }

    //left and right mirrors
    for(int i=0;i<NbRows;i++)
        for(int j=0;j<W;j++)
        {
            MirroredImage[i+H][j] =(*this)[i][W-j]; //left mirror
            MirroredImage[i+H][j+NbCols+W] =(*this)[i][NbCols-2-j]; //right mirror
        }

    //top and bottom mirrors
    for(int j=0;j<NbCols;j++)
        for(int i=0;i<H;i++)
        {
            MirroredImage[i][j+W] =(*this)[H-i][j]; //top mirror
            MirroredImage[i+NbRows+H][j+W] = (*this)[NbRows-2-i][j]; //bottom mirror
        }

    //corners
    for(int i=0;i<H;i++)
        for(int j=0;j<W;j++)
        {
            MirroredImage[i][j] = (*this)[H-i][W-j]; //NW corner
            MirroredImage[i+NbRows+H][j] = (*this)[NbRows-2-i][W-j]; //SW corner
            MirroredImage[i][j+NbCols+W] = (*this)[H-i][NbCols-2-j]; //NE corner
            MirroredImage[i+NbRows+H][j+NbCols+W] = (*this)[NbRows-2-i][NbCols-2-j]; //SE corner
        }
}


//Set a block to the values
void M3Matrix::AddBlock(M3Matrix& BlockMatrix, int rStart, int cStart)
{
    for(int i=0;i<BlockMatrix.GetNbRows();i++)
        for(int j=0;j<BlockMatrix.GetNbCols();j++)
        {
            if(i+rStart < GetNbRows() && j+rStart < GetNbCols())
                (*this)[i+rStart][j+cStart] += BlockMatrix[i][j];
        }

}


//Set a block to the values
void M3Matrix::SetBlock(M3Matrix& BlockMatrix, int rStart, int cStart)
{
    for(int i=0;i<BlockMatrix.GetNbRows();i++)
        for(int j=0;j<BlockMatrix.GetNbCols();j++)
        {
            if(i+rStart < GetNbRows() && j+cStart < GetNbCols())
            {
                (*this)[i+rStart][j+cStart] = BlockMatrix[i][j];
            }
        }
}


// Extraction of  a block of size (NbRows,NbCols)
// starting at (rStart,cStart)
void M3Matrix::ExtractBlock(M3Matrix& BlockMatrix, int rStart, int cStart, int NbRows,int NbCols)
{
    if(BlockMatrix.GetNbRows() != (*this).GetNbRows() || BlockMatrix.GetNbCols() != (*this).GetNbCols())
        BlockMatrix.Reshape(NbRows,NbCols);

    for(int i=0;i<NbRows;i++)
        for(int j=0;j<NbCols;j++)
        {
            BlockMatrix[i][j] = (*this)[i+rStart][j+cStart];
        }
}


M3Matrix M3Matrix::ExtractBlock(int rStart, int cStart, int NbRows,int NbCols)
{
    M3Matrix BlockMatrix;
    if(rStart+NbRows > (*this).GetNbRows())
        NbRows =(*this).GetNbRows() -rStart;
    if(cStart+NbCols > (*this).GetNbCols())
        NbCols =(*this).GetNbCols() -cStart;

    (*this).ExtractBlock(BlockMatrix,rStart, cStart, NbRows, NbCols);
    return BlockMatrix;
}


void M3Matrix::LinearFilter_Hor_DS(M3Matrix& FilteredImage, M3Matrix& Mask)
{
    M3Matrix MirroredImage;
    int HalfHeight = Mask.GetNbRows()/2;
    int HalfWidth = Mask.GetNbCols()/2;
    (*this).Mirror(MirroredImage,HalfHeight,HalfWidth);
    FilteredImage.Reshape(t_NbRows,t_NbCols/2);

    for(int i = 0;i<t_NbRows;i++)
        for(int j = 0;j<t_NbCols/2;j++)
        {
            double FilteredValue = 0.;
            for(int m=0;m<Mask.GetNbRows();m++)
                for(int n=0;n<Mask.GetNbCols();n++)
                {
                    FilteredValue += MirroredImage[i+m][2*j+n+1] * Mask[m][n];
                }

            FilteredImage[i][j] = FilteredValue;
        }
}


void M3Matrix::LinearFilter_Ver_DS(M3Matrix& FilteredImage, M3Matrix& Mask)
{
    M3Matrix MirroredImage;
    int HalfHeight = Mask.GetNbRows()/2;
    int HalfWidth = Mask.GetNbCols()/2;

    (*this).Mirror(MirroredImage,HalfHeight,HalfWidth);
    FilteredImage.Reshape(t_NbRows/2,t_NbCols);

    for(int i = 0;i<t_NbRows/2;i++)
        for(int j = 0;j<t_NbCols;j++)
        {
            double FilteredValue = 0.;
            for(int m=0;m<Mask.GetNbRows();m++)
                for(int n=0;n<Mask.GetNbCols();n++)
                {
                    FilteredValue += MirroredImage[2*i+m+1][j+n] * Mask[m][n];
                }

            FilteredImage[i][j] = FilteredValue;
        }
}


void M3Matrix::LinearTemporalFilter(M3Matrix& FilteredImage, M3Matrix3D& Mask3D, M3Matrix3D& CircPast, int CircIndx /*, int NbPastFrames*/)
{
    int LocalIndx;
    int NbPastFrames = CircPast.GetNbFrames();
    (*this).LinearFilter(FilteredImage,Mask3D[0]);

    M3Matrix Tmp,TmpFiltered;
    for(int i = 1;i<=NbPastFrames;i++)
    {
        LocalIndx = (CircIndx-i)%NbPastFrames; //TODO: check this line !!! is the index correct (07.05.2015)

        Tmp = CircPast[LocalIndx];
        Tmp.LinearFilter(TmpFiltered,Mask3D[i]);

        FilteredImage += TmpFiltered;
    }
}


void M3Matrix::DoG(M3Matrix& B, double sigmaLow,double sigmaHigh)
{
    M3Matrix C,D;
    C = (*this).GaussianBlurring(sigmaLow);
    D =  (*this).GaussianBlurring(sigmaHigh);
    B=(C-D);
}



M3Matrix M3Matrix::GaussianBlurring(double sigma)
{
    M3Matrix BlurFilter;
    BlurFilter.MakeGaussianFilter(sigma,9);

    M3Matrix Filtered;
    Filtered = LinearFiltering(BlurFilter);
    return Filtered;
}

void M3Matrix::GaussianBlur(double sigma)
{
    M3Matrix BlurFilter;
    BlurFilter.MakeGaussianFilter(sigma,9);

    LinearFilter(BlurFilter);
}



M3Matrix M3Matrix::LinearFiltering(M3Matrix& Mask)
{
    M3Matrix MirroredImage;
    int HalfHeight = Mask.GetNbRows()/2;
    int HalfWidth = Mask.GetNbCols()/2;
    (*this).Mirror(MirroredImage,HalfHeight,HalfWidth);
    M3Matrix FilteredImage(t_NbRows,t_NbCols);
    FilteredImage.InitialiseDynamic(*this);

    for(int i = 0;i<t_NbRows;i++)
        for(int j = 0;j<t_NbCols;j++)
        {
            double FilteredValue = 0.;
            for(int m=0;m<Mask.GetNbRows();m++)
                for(int n=0;n<Mask.GetNbCols();n++)
                {
                    FilteredValue += MirroredImage[i+m][j+n] * Mask[m][n];
                }

            FilteredImage[i][j] = FilteredValue;
        }

    return FilteredImage;
}



void M3Matrix::LinearFilter(M3Matrix& Mask)
{
    M3Matrix MirroredImage;
    int HalfHeight = Mask.GetNbRows()/2;
    int HalfWidth = Mask.GetNbCols()/2;
    (*this).Mirror(MirroredImage,HalfHeight,HalfWidth);
    M3Matrix FilteredImage(t_NbRows,t_NbCols);
    FilteredImage.InitialiseDynamic(*this);

    for(int i = 0;i<t_NbRows;i++)
        for(int j = 0;j<t_NbCols;j++)
        {
            double FilteredValue = 0.;
            for(int m=0;m<Mask.GetNbRows();m++)
                for(int n=0;n<Mask.GetNbCols();n++)
                {
                    FilteredValue += MirroredImage[i+m][j+n] * Mask[m][n];
                }

            FilteredImage[i][j] = FilteredValue;
        }

    (*this) = FilteredImage;
}


void M3Matrix::LinearFilter(M3Matrix& FilteredImage, M3Matrix& Mask)
{
    cout << "Maks: NbRows: " << Mask.GetNbRows() << "  Mask.GetNbcols : " << Mask.GetNbCols() << endl;

    M3Matrix MirroredImage;
    int HalfHeight = Mask.GetNbRows()/2;
    int HalfWidth = Mask.GetNbCols()/2;
    (*this).Mirror(MirroredImage,HalfHeight,HalfWidth);  FilteredImage.Reshape(t_NbRows,t_NbCols);
    FilteredImage.InitialiseDynamic(*this);


    for(int i = 0;i<t_NbRows;i++)
        for(int j = 0;j<t_NbCols;j++)
        {
            double FilteredValue = 0.;
            for(int m=0;m<Mask.GetNbRows();m++)
                for(int n=0;n<Mask.GetNbCols();n++)
                {
                    FilteredValue += MirroredImage[i+m][j+n] * Mask[m][n];
                }

            FilteredImage[i][j] = FilteredValue;
        }
}



void M3Matrix::MakeLaws1D(int LawsIndx)
{
    double* T;

    if(LawsIndx == 1) //Level
    {
        T = new double[5] {1./16,4./16.,6./16.,4./16.,1./16.};
    }
    else if(LawsIndx == 2) //Edge
    {
        T = new double[5] {-1,-2,0,2,1};
    }
    else if(LawsIndx == 3) //Spot
    {
        T = new double[5] {-1,0,2,0,-1};
    }
    else if(LawsIndx == 4) //Ripple
    {
        T = new double[5] {1,-4,6,-4,1};
    }
    else if(LawsIndx == 5) //Wave
    {
        T = new double[5] {-1,2,0,-2,-1};
    }
    else //default set to Level
        T = new double[5] {1./16,4./16.,6./16.,4./16.,1./16.};

    SetData(1,5,T);
    delete T;
}


void M3Matrix::MakeLawsFilter(int VertLawsIndx, int HorLawsIndx)
{
    Reshape(5,5);

    M3Matrix VertLaws, HorLaws;
    VertLaws.MakeLaws1D(VertLawsIndx);
    HorLaws.MakeLaws1D(HorLawsIndx);

    for(int i = 0;i<5;i++)
    {
        for(int j = 0;j<5;j++)
        {
            (*this)[i][j] = VertLaws[0][i] * HorLaws[0][j];
        }
    }
}


//Kirsch Filters
void M3Matrix::MakeKirschFilter(int KirschIndx)
{
    double* T;
    if(KirschIndx == 1)
    {
        T = new double[9]
        {
                5, -3, -3,
                5,  0, -3,
                5, -3, -3
    };
    }
    else if(KirschIndx == 2)
    {
        T = new double[9]
        {
                5, 5, -3,
                5,  0, -3,
                -3, -3, -3
    };
    }
    else if(KirschIndx == 3)
    {
        T = new double[9]
        {
                5, 5, 5,
                -3,  0, -3,
                -3, -3, -3
    };
    }
    else if(KirschIndx == 4)
    {
        T = new double[9]
        {
                -3, 5, 5,
                -3,  0, 5,
                -3, -3, -3
    };
    }
    else if(KirschIndx == 5)
    {
        T = new double[9]
        {
                -3, -3, 5,
                -3,  0, 5,
                -3, -3, 5
    };
    }
    else if(KirschIndx == 6)
    {
        T = new double[9]
        {
                -3, -3, -3,
                -3,  0, 5,
                -3, 5, 5
    };
    }
    else if(KirschIndx == 7)
    {
        T = new double[9]
        {
                -3, -3, -3,
                -3,  0, -3,
                5, 5, 5
    };
    }
    else if(KirschIndx == 8)
    {
        T = new double[9]
        {
                -3, -3, -3,
                5,  0, -3,
                5, 5, -3
    };
    }
    else
    {
        T = new double[9]
        {
                5, -3, -3,
                5,  0, -3,
                5, -3, -3
        };
    }

    SetData(3,3,T);
    delete T;
}


void M3Matrix::MakeFreiChenFilter(int FreiChenIndx)
{
    double* T;
    if(FreiChenIndx == 0)
    {
        T = new double[9]
        {
                1, 1, 1,
                1,  1, 1,
                1, 1, 1
    };
    }
    else if(FreiChenIndx == 1)
    {
        T = new double[9]
        {
                -1, sqrt(2.0), -1,
                0,  0, 0,
                1, sqrt(2.0), 1
    };
    }
    else if(FreiChenIndx == 2)
    {
        T = new double[9]
        {
                -1, 0, 1,
                sqrt(2),  0, sqrt(2),
                -1, 0, 1
    };
    }
    else if(FreiChenIndx == 3)
    {
        T = new double[9]
        {
                0, -1, sqrt(2),
                1,  0, -1,
                -sqrt(2), 1, 0
    };
    }
    else if(FreiChenIndx == 4)
    {
        T = new double[9]
        {
                sqrt(2), -1, 0,
                -1,  0, 1,
                0, 1, sqrt(2)
    };
    }
    else if(FreiChenIndx == 5)
    {
        T = new double[9]
        {
                0, 1, 0,
                -1,  0, 1,
                0, -1, 0
    };
    }
    else if(FreiChenIndx == 6)
    {
        T = new double[9]
        {
                -1, 0, 1,
                0,  0, 0,
                1, 0, -1
    };
    }
    else if(FreiChenIndx == 7)
    {
        T = new double[9]
        {
                1, -2, 1,
                -2,  4, -2,
                1, -2, 1
    };
    }
    else if(FreiChenIndx == 8)
    {
        T = new double[9]
        {
                -2, 1, -2,
                1,  4, 1,
                -2, 1, -2
    };
    }
    else
    {
      T = new double[9]
      {
            1, 1, 1,
            1,  1, 1,
            1, 1, 1
      };
    }

    SetData(3,3,T);
    (*this).Normalise();

    delete T;
}


void M3Matrix::MakeGaussianFilter3()
{
    (*this).Reshape(3,3);

    (*this)[0][0]=1./16.;
    (*this)[0][1]=2./16.;
    (*this)[0][2]=1./16.;
    (*this)[1][0]=2./16.;
    (*this)[1][1]=4./16.;
    (*this)[1][2]=2./16.;
    (*this)[2][0]=1./16.;
    (*this)[2][1]=2./16.;
    (*this)[2][2]=1./16.;
}


void M3Matrix::MakeGaussianFilter5()
{
    (*this).Reshape(5,5);
    (*this).SetValues({1., 4., 7., 4.,1.,
                       4.,16.,26.,16.,4.,
                       7.,26.,41.,26.,7.,
                       4.,16.,26.,16.,4.,
                       1., 4., 7., 4.,1.,
                      });
    (*this) = (*this)/273;
}


void M3Matrix::MakeGaussianFilter7()
{
    (*this).Reshape(7,7);

    (*this)[0][0] = 0.;
    (*this)[0][1] = 0.;
    (*this)[0][2] = 1.;
    (*this)[0][3] = 2.;
    (*this)[0][4] = 1.;
    (*this)[0][5] = 0.;
    (*this)[0][6] = 0.;

    (*this)[1][0] = 0.;
    (*this)[1][1] = 3.;
    (*this)[1][2] = 13.;
    (*this)[1][3] = 22.;
    (*this)[1][4] = 13.;
    (*this)[1][5] = 3.;
    (*this)[1][6] = 0.;

    (*this)[2][0] = 1.;
    (*this)[2][1] = 13.;
    (*this)[2][2] = 59.;
    (*this)[2][3] = 97.;
    (*this)[2][4] = 59.;
    (*this)[2][5] = 13.;
    (*this)[2][6] = 1.;

    (*this)[3][0] = 2.;
    (*this)[3][1] = 22.;
    (*this)[3][2] = 97.;
    (*this)[3][3] = 159.;
    (*this)[3][4] = 97.;
    (*this)[3][5] = 22.;
    (*this)[3][6] = 2.;

    (*this)[4][0] = 1.;
    (*this)[4][1] = 13.;
    (*this)[4][2] = 59.;
    (*this)[4][3] = 97.;
    (*this)[4][4] = 59.;
    (*this)[4][5] = 13.;
    (*this)[4][6] = 1.;

    (*this)[5][0] = 0.;
    (*this)[5][1] = 3.;
    (*this)[5][2] = 13.;
    (*this)[5][3] = 22.;
    (*this)[5][4] = 13.;
    (*this)[5][5] = 3.;
    (*this)[5][6] = 0.;

    (*this)[6][0] = 0.;
    (*this)[6][1] = 0.;
    (*this)[6][2] = 1.;
    (*this)[6][3] = 2.;
    (*this)[6][4] = 1.;
    (*this)[6][5] = 0.;
    (*this)[6][6] = 0.;

    (*this) = (*this)/1000.;
}


void M3Matrix::MakeGaussianFilter7_1D()
{
    (*this).Reshape(1,7);

    (*this)[0][0] =  0.004;
    (*this)[0][1] =  0.054;
    (*this)[0][2] =  0.242;
    (*this)[0][3] =  0.399;
    (*this)[0][4] =  0.242;
    (*this)[0][5] =  0.054;
    (*this)[0][6]  = 0.004;

    (*this).Normalise();
}


void M3Matrix::MakeGaussianFilter(double sigma, int Size)
{
    (*this).Reshape(Size,Size);

    if(sigma > 0)
    {
        for(int i = 0;i<Size;i++)
            for(int j = 0;j<Size;j++)
            {
                double  id = double(i-Size/2);
                double  jd = double(j-Size/2);
                (*this)[i][j] = exp(-0.5*(id*id+jd*jd)/(sigma*sigma));
            }
        (*this).Normalise();
    }
    else
    {
        (*this).SetValues(0.);
        (*this)[Size/2][Size/2] = 1.0;
    }
}


void M3Matrix::MakeGaborImagFilter(double sigma, int Size, double omega, double theta, double gamma)
{
    (*this).Reshape(Size,Size);

    if(sigma > 0)
    {
        for(int i = 0;i<Size;i++)
            for(int j = 0;j<Size;j++)
            {
                double  id = double(i-Size/2);
                double  jd = double(j-Size/2);
                double x = jd;
                double y = -id;

                double xp = x * cos(theta) + y * sin(theta);
                double yp = - x * sin(theta) + y * cos(theta);

                (*this)[i][j] = exp(-0.5*(xp*xp+gamma*yp*yp)/(sigma*sigma)) * sin(omega*xp);
            }
    }
    else
    {
        (*this).SetValues(0.);
        (*this)[Size/2][Size/2] = 1.0;
    }
}


void M3Matrix::MakeGaborRealFilter(double sigma, int Size, double omega, double theta, double gamma)
{
    (*this).Reshape(Size,Size);

    if(sigma > 0)
    {
        for(int i = 0;i<Size;i++)
            for(int j = 0;j<Size;j++)
            {
                double  id = double(i-Size/2);
                double  jd = double(j-Size/2);
                double x = jd;
                double y = -id;

                double xp = x * cos(theta) + y * sin(theta);
                double yp = - x * sin(theta) + y * cos(theta);

                (*this)[i][j] = exp(-0.5*(xp*xp+gamma*yp*yp)/(sigma*sigma)) * cos(omega*xp);
            }
    }
    else
    {
        (*this).SetValues(0.);
        (*this)[Size/2][Size/2] = 1.0;
    }
}



void M3Matrix::MakeGaussianFilter_1D(double sigma, int Size)
{
    (*this).Reshape(1,Size);

    if(sigma > 0)
    {
        for(int j = 0;j<Size;j++)
        {
            double  jd = double(j-Size/2);
            (*this)[0][j] = exp(-0.5*(jd*jd)/(sigma*sigma));
        }
        (*this).Normalise();
    }
    else
    {
        (*this).SetValues(0.);
        (*this)[0][Size/2] = 1.0;
    }
}


//Local fractal dimension
//s. Haefner, Tamaki,Tanaka, Uhl, Wimmer, Yoshida
M3Matrix M3Matrix::LocalFractalDimension()
{
    M3Matrix LFD(t_NbRows,t_NbCols);//Local fractal dimension

    enum MuDef {mu1,mu2,mu3};
    MuDef option = mu1;
    int MaxSize = 9; //Maximal size of the boxes
    vector<double> LogMeasures(MaxSize);

    M3Matrix Extension;
    (*this).Mirror(Extension,MaxSize,MaxSize);

    for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols;j++)
        {
            for(int s = 1;s<MaxSize;s++)
            {
                double tmp = 0.;
                if(option == mu1)
                {
                    for(int h=-s;h<=s;h++)
                    {

                        for(int k=-s;k<=s;k++)
                            tmp += Extension[s+i+h][s+j+k]; //implement other measures !
                    }
                }
                else if(option == mu2) //This is formula (4) in the above mentioned paper
                {

                }

                LogMeasures[s-1] = log(tmp)/log(2*s+1);
            }

            double FracDim = 0.;
            for(int s = 0;s<MaxSize;s++)
            {
                FracDim += LogMeasures[s];
            }
            FracDim /= MaxSize;
            LFD[i][j] = FracDim;
        }


    LFD /= LFD.GetMax();
    LFD *= 255;

    LFD -= 160;
    LFD *=1.6;

    LFD.Digitalize();

    return LFD;
}


//This function is the forward shading operator under assumption of vertical light and
//homogeneous reflectance of the target surface
M3Matrix M3Matrix::ForwardShading()
{
    M3Matrix L(t_NbRows,t_NbCols), GradH(t_NbRows,t_NbCols), GradV(t_NbRows,t_NbCols);
    GradH = (*this).GradientH_Centered();
    GradV = (*this).GradientV_Centered();

    double factor=0.25;
    for(int i = 0 ;i<t_NbRows;i++)
        for(int j = 0 ;j<t_NbCols;j++)
        {
            double gradienL2Norm = factor*factor*(GradH[i][j] * GradH[i][j] + GradV[i][j] * GradV[i][j]);
            L[i][j] = 1./sqrt(1 + gradienL2Norm);
        }

    L /= L.GetMax();
    L*=255;

    return L;
}


//Edge detectors ("differential" operators)
void M3Matrix::MakeGradientC()
{
    (*this).Reshape(1,3);
    (*this)[0][0] = 0.;
    (*this)[0][1] = -1.;
    (*this)[0][2] = 1.;
}


void M3Matrix::MakeGradientR()
{
    (*this).Reshape(3,1);
    (*this)[0][0] = 0.;
    (*this)[1][0] = -1.;
    (*this)[2][0] = 1.;
}


void M3Matrix::MakePrewittC()
{
    (*this).Reshape(3,3);

    (*this)[0][0] = -1;
    (*this)[0][1] = 0;
    (*this)[0][1] = 1;
    (*this)[1][0] = -1;
    (*this)[1][1] = 0;
    (*this)[1][1] = 1;
    (*this)[2][0] = -1;
    (*this)[2][1] = 0;
    (*this)[2][1] = 1;
}


void M3Matrix::MakePrewittR()
{
    (*this).Reshape(3,3);

    (*this)[0][0] = -1;
    (*this)[0][1] = -1;
    (*this)[0][1] = -1;
    (*this)[1][0] = 0;
    (*this)[1][1] = 0;
    (*this)[1][1] = 0;
    (*this)[2][0] = 1;
    (*this)[2][1] = 1;
    (*this)[2][1] = 1;
}


void M3Matrix::MakeSobelC()
{
    (*this).Reshape(3,3);

    (*this)[0][0] = -1;
    (*this)[0][1] = 0;
    (*this)[0][1] = 1;
    (*this)[1][0] = -2;
    (*this)[1][1] = 0;
    (*this)[1][1] = 2;
    (*this)[2][0] = -1;
    (*this)[2][1] = 0;
    (*this)[2][1] = 1;
}


void M3Matrix::MakeSobelR()
{
    (*this).Reshape(3,3);

    (*this)[0][0] = -1;
    (*this)[0][1] = -2;
    (*this)[0][1] = -1;
    (*this)[1][0] = 0;
    (*this)[1][1] = 0;
    (*this)[1][1] = 0;
    (*this)[2][0] = 1;
    (*this)[2][1] = 2;
    (*this)[2][1] = 1;
}


void M3Matrix::MinMedianMaxFilter(int S,double threshold, double factor)
{
    M3Matrix MirroredImage, Block;
    double median,dmax,dmin;

    (*this).Mirror(MirroredImage,S/2,S/2);

    M3Matrix FilteredImage;
    FilteredImage.Reshape(t_NbRows,t_NbCols);

    for(int i = 0;i<t_NbRows;i++)
        for(int j = 0;j<t_NbCols;j++)
        {
            MirroredImage.ExtractBlock(Block,i,j,S,S);

            median = Block.Median();
            dmin = Block.GetMin();
            dmax = Block.GetMax();

            if(median - dmin > threshold && dmax - median > threshold) //test for overshoots in the block ?  look for a better one ?
            {
                FilteredImage[i][j] = median + factor*(MirroredImage[i+S/2][j+S/2] - median);
            }
            else
                FilteredImage[i][j] = MirroredImage[i+S/2][j+S/2];
        }

    for(int i = 0;i < t_NbRows;i++)
        for(int j = 0;j < t_NbCols;j++)
        {
            (*this)[i][j] = FilteredImage[i][j];
        }
}


void M3Matrix::DirectionalMedianFilter(M3Matrix& Filtered, int Size)
{
    int Height = Size;
    int Width = Size;

    M3Matrix MirroredImage, Block,Line;
    double median;
    double medianV,medianH,medianD1,medianD2;
    (*this).Mirror(MirroredImage,Height/2,Width/2);
    Filtered.Reshape(t_NbRows,t_NbCols);

    for(int i = 0;i<t_NbRows;i++)
        for(int j = 0;j<t_NbCols;j++)
        {
            MirroredImage.ExtractBlock(Block,i,j,Height,Width);

            //Compute the directional medians
            Line.Reshape(1,Size);

            //Horizontal
            for(int k=0;k<Size;k++)
                Line[0][k] = Block[Height/2][k];

            medianH = Line.Median();

            //Vertical
            for(int k=0;k<Size;k++)
                Line[0][k] = Block[k][Width/2];

            medianV = Line.Median();

            //Principal diagonal
            for(int k=0;k<Size;k++)
                Line[0][k] = Block[k][k];

            medianD1 = Line.Median();

            //Secondary diagonal
            for(int k=0;k<Size;k++)
                Line[0][k] = Block[Height-k-1][k];

            medianD2 = Line.Median();

            median = min4(medianH,medianV,medianD1,medianD2);
            Filtered[i][j] = median;
        }
}


void M3Matrix::MedianFilter(M3Matrix& MedianFilteredImage, int Height, int Width)
{
    int H = Height, W = Width;
    M3Matrix MirroredImage, Block;
    double median;

    (*this).Mirror(MirroredImage,H/2,W/2);
    MedianFilteredImage.Reshape(t_NbRows,t_NbCols);

    for(int i = 0;i<t_NbRows;i++)
        for(int j = 0;j<t_NbCols;j++)
        {
            MirroredImage.ExtractBlock(Block,i,j,H,W);
            median = Block.Median();
            MedianFilteredImage[i][j] = median;
        }
}


void M3Matrix::MedMedFilter(M3Matrix& MedianFilteredImage, int Height, int Width)
{
    int H = Height, W = Width;
    M3Matrix MirroredImage, Block;
    double medmed;

    (*this).Mirror(MirroredImage,H/2,W/2);
    MedianFilteredImage.Reshape(t_NbRows,t_NbCols);

    for(int i = 0;i<t_NbRows;i++)
        for(int j = 0;j<t_NbCols;j++)
        {
            MirroredImage.ExtractBlock(Block,i,j,H,W);
            medmed = Block.MedMed();
            MedianFilteredImage[i][j] = medmed;
        }
}


void M3Matrix::MakeRampFilter(int Height,int Width)
{
    Reshape(Height,Width);

    for(int i = 0;i <Height;i++)
        for(int j = 0;j < Width;j++)
        {
            if(i == Height/2 && j==Height/2)
            {
                (*this)[i][j] = 0;
            }
            else
            {
                double id = (double)(i-Height/2);
                double jd = (double)(j-Width/2);

                (*this)[i][j] = sqrt((fabs(id) + fabs(jd)));
            }
        }

    Normalise();
}


void M3Matrix::MakeBoxFilter(int Height,int Width)
{
    (*this).Reshape(Height,Width);
    (*this).SetValues(1.0/(Height*Width));
}


//Smoothing filters
void M3Matrix::MakeBoxFilter(double /*Value*/, int Height,int Width)
{
    (*this).Reshape(Height,Width);
    (*this).SetValues(1.0/(Height*Width));
}


//Like a box filter but with 4 zero values at the corners
void M3Matrix::MakeCircBoxFilter(double /*Value*/, int Height,int Width)
{
    (*this).Reshape(Height,Width);
    (*this).SetValues(1.0/(Height*Width-4));
    (*this)[0][0] = 0.;
    (*this)[Height-1][0] = 0.;
    (*this)[0][Width-1] = 0.;
    (*this)[Height-1][Width-1] = 0.;
}


void M3Matrix::MakeLaplaceFilter(LAPLACIANTYPE laplacian_type)
{
    (*this).Reshape(3,3);

    double* T;
    if(laplacian_type == anisotropic)
    {
        T = new double[9]
        {
                -1., -1., -1.,
                -1.,  8., -1.,
                -1., -1., -1.
    };
    }
    else if (laplacian_type == isotropic)
    {
        T = new double[9]
        {
                -1, -2., -1.,
                -2.,  12., -2.,
                -1., -2., -1.
    };
    }
    else if(laplacian_type == lothar)
    {
        T = new double[9]
        {
                -1, -4., -1.,
                -4.,  20., -4.,
                -1., -4., -1.
    };
    }
    else
    {
        T = new double[9]
        {
                -1, -2., -1.,
                -2.,  12., -2.,
                -1., -2., -1.
    };
    }

    SetData(3,3,T);
    delete T;

    if(laplacian_type == isotropic)
        (*this) /= 24;
    if(laplacian_type == lothar)
        (*this) /= 40;
}


void M3Matrix::ChessBoard()
{
    int H = GetNbRows();
    int W = GetNbCols();
    for(int i = 0; i < H; i++)
        for(int j = 0; j < W; j++)
        {
            if((i+j)%2 == 1)
                (*this)[i][j] = 0.;
        }
    (*this).Normalise();
}


void M3Matrix::Sparsify(int SamplingStep)
{
    int H = GetNbRows();
    int W = GetNbCols();
    for(int i = 0; i < H; i++)
        for(int j = 0; j < W; j++)
        {
            if(!((i-H/2)%SamplingStep==0 && (j-W/2)%SamplingStep==0))
                (*this)[i][j] = 0.;
        }
    (*this).Normalise();
}


void M3Matrix::Hollow()
{
    int H = GetNbRows();
    int W = GetNbCols();
    for(int i = H/2-1; i <= H/2+1; i++)
        for(int j = W/2-1; j <= W/2+1; j++)
        {
            if(i != H/2 || j != H/2)
                (*this)[i][j] = 0.;
        }

    (*this).Normalise();
}


void M3Matrix::Rescale(M3Matrix& RescaledMatrix, double MinVal, double MaxVal)
{
    double OrgMax = GetMax();
    double OrgMin = GetMin();
    double OrgRange = OrgMax - OrgMin;

    double val;
    double RescaledVal = 0.;
    RescaledMatrix.Reshape(GetNbRows(),GetNbCols());

    for(int i = 0;i<GetNbRows();i++)
        for(int j = 0;j<GetNbCols();j++)
        {
            val = (*this)[i][j];
            RescaledVal = MinVal + (val-OrgMin)/OrgRange*(MaxVal-MinVal);
            RescaledMatrix[i][j] = RescaledVal;
        }
}


void M3Matrix::LaplaceEnhancement(M3Matrix& FilteredImage, double alpha)
{
    M3Matrix LaplaceFilter;
    LaplaceFilter.MakeLaplaceFilter();
    (*this).LinearFilter(FilteredImage,LaplaceFilter);
    FilteredImage = (*this) + FilteredImage * alpha;
}


void M3Matrix::UnsharpMasking_Box(M3Matrix& FilteredImage, int Height,int Width)
{
    M3Matrix Filter;
    Filter.MakeBoxFilter(1./(Height*Width),Height,Width);
    (*this).LinearFilter(FilteredImage,Filter);
    FilteredImage = (*this) - FilteredImage;
}


// Wedge filters
void M3Matrix::MakeWedgeLeft(int Size)
{
    (*this).Reshape(Size,Size);
    (*this).SetValues(0.);

    int HalfSize = Size/2;
    double Area = Size*(1.+HalfSize);

    for(int i=0;i<Size;i++)
        for(int j=0;j<1+Size/2;j++)
        {
            (*this)[i][j] = 1./Area;
        }
}


void M3Matrix::MakeWedgeRight(int Size)
{
    (*this).Reshape(Size,Size);
    int HalfSize = Size/2;
    double Area = Size*(1.+HalfSize);

    for(int i=0;i<Size;i++)
        for(int j=Size/2;j<Size;j++)
        {
            (*this)[i][j] = 1./Area;
        }
}


void M3Matrix::MakeWedgeTop(int Size)
{
    (*this).Reshape(Size,Size);
    int HalfSize = Size/2;
    double Area = Size*(1.+HalfSize);

    for(int i=0;i<Size/2+1;i++)
        for(int j=0;j<Size;j++)
        {
            (*this)[i][j] = 1./Area;
        }
}


void M3Matrix::MakeWedgeBottom(int Size)
{
    (*this).Reshape(Size,Size);
    int HalfSize = Size/2;
    double Area = Size*(1.+HalfSize);

    for(int i=Size/2;i<Size;i++)
        for(int j=0;j<Size;j++)
        {
            (*this)[i][j] = 1./Area;
        }
}


void M3Matrix::MakeWedgeNE(int Size)
{
    (*this).Reshape(Size,Size);
    double Area = 0.5*Size*(Size+1);

    for(int i=0;i<Size;i++)
        for(int j=i;j<Size;j++)
        {
            (*this)[i][j] = 1./Area;
        }
}


void M3Matrix::MakeWedgeNW(int Size)
{
    (*this).Reshape(Size,Size);

    double Area = 0.5*Size*(Size+1);
    for(int i=0;i<Size;i++)
        for(int j=0;j<Size-i;j++)
        {
            (*this)[i][j] = 1./Area;
        }
}


void M3Matrix::MakeWedgeSW(int Size)
{
    (*this).Reshape(Size,Size);
    double Area = 0.5*Size*(Size+1);

    for(int i=0;i<Size;i++)
        for(int j=0;j<=i;j++)
        {
            (*this)[i][j] = 1./Area;
        }
}


void M3Matrix::MakeWedgeSE(int Size)
{
    (*this).Reshape(Size,Size);
    double Area = 0.5*Size*(Size+1);

    for(int i=0;i<Size;i++)
        for(int j=Size-1-i;j<Size;j++)
        {
            (*this)[i][j] = 1./Area;
        }
}


double M3Matrix::WedgeDeviation(M3Matrix& WedgeFilter,int i, int j, int Criterion)
{
    double Deviation = 0.;
    if(Criterion == 0) //Variance
    {
        Deviation =   WedgeVariance(WedgeFilter,i,j);
    }
    else if(Criterion == 1)
    {
        Deviation =   WedgeAbsDev(WedgeFilter,i,j);
    }

    return Deviation;
}


double M3Matrix::WedgeAbsDev(M3Matrix& WedgeFilter,int i,int j)
{
    int NbPixels = 0;
    int H = WedgeFilter.GetNbRows();
    int W = WedgeFilter.GetNbCols();
    double val = 0., Sum=0.;

    for(int m=0;m<H;m++)
    {
        for(int n=0;n<W;n++)
        {
            if(WedgeFilter[m][n] > 0)
            {
                NbPixels++;
                //WARNING: This must be a mirrored image, and i,j must be between w and NbRows+w etc...
                val = (*this)[i+m-W/2][j+n-H/2];
                Sum += fabs(val);
            }
        }
    }
    double absdev = Sum/NbPixels;
    return absdev;
}


double M3Matrix::WedgeVariance(M3Matrix& WedgeFilter,int i,int j)
{
    int NbPixels = 0;
    int H = WedgeFilter.GetNbRows();
    int W = WedgeFilter.GetNbCols();
    double val, Sum=0.,SumSq=0.;

    for(int m=0;m<H;m++)
    {
        for(int n=0;n<W;n++)
        {
            if(WedgeFilter[m][n] > 0)
            {
                NbPixels++;
                //WARNING: This must be a mirrored image, and i,j must be between w and NbRows+w etc...
                val = (*this)[i+m-W/2][j+n-H/2];
                Sum += val;
                SumSq += val*val;
            }
        }
    }
    double variance = SumSq/NbPixels-(Sum/NbPixels)*(Sum/NbPixels);
    return variance;
}


//FilterType: 0: box filter, 1 gaussian filter
//Criterion: 0: variance, 1 abs. deviation
void M3Matrix::SelectiveWedgeFilter(M3Matrix& WedgeFilteredImage, int Size, int FilterType, double tau, int Criterion, double sigma)
{
    M3Matrix MirroredImage,WedgeFilter,Block,BestWedgeFilter,GaussFilter;

    (*this).Mirror(MirroredImage,Size/2,Size/2);
    WedgeFilteredImage.Reshape(t_NbRows,t_NbCols);

    int H = Size/2;
    int W = Size/2;
    double CurrentDeviation,MinDeviation,GlobalDeviation;

    cout << "SelectiveWedgeFilter :"<< endl << "Criterion : ";
    if(Criterion == 0)
        cout << "Variance" << endl;
    else
        cout << "Absolute deviation" << endl;

    if(FilterType == 0)
        cout << " FilterType for the final smoothing : box" << endl;
    else
        cout << " FilterType for the final smoothing : gaussian" << endl;

    //Loop on the image pixels
    for(int i=0;i<t_NbRows;i++)
    {
        for(int j=0;j<t_NbCols;j++)
        {
            //Switch between Box filter and gaussian
            WedgeFilter.MakeBoxFilter(1./(Size*Size),Size,Size);
            GlobalDeviation = MirroredImage.WedgeDeviation(WedgeFilter,i+H,j+W,Criterion);

            MinDeviation = GlobalDeviation*tau;
            BestWedgeFilter = WedgeFilter;

            //Select wedge filter
            WedgeFilter.MakeWedgeLeft(Size);
            CurrentDeviation = MirroredImage.WedgeDeviation(WedgeFilter,i+H,j+W,Criterion);

            if(CurrentDeviation < MinDeviation)
            {
                MinDeviation = CurrentDeviation;
                BestWedgeFilter = WedgeFilter;
            }

            WedgeFilter.MakeWedgeRight(Size);
            CurrentDeviation = MirroredImage.WedgeDeviation(WedgeFilter,i+H,j+W,Criterion);

            if(CurrentDeviation < MinDeviation)
            {
                MinDeviation = CurrentDeviation;
                BestWedgeFilter = WedgeFilter;
            }

            WedgeFilter.MakeWedgeTop(Size);
            CurrentDeviation = MirroredImage.WedgeDeviation(WedgeFilter,i+H,j+W,Criterion);

            if(CurrentDeviation < MinDeviation)
            {
                MinDeviation = CurrentDeviation;
                BestWedgeFilter = WedgeFilter;
            }

            WedgeFilter.MakeWedgeBottom(Size);
            CurrentDeviation = MirroredImage.WedgeDeviation(WedgeFilter,i+H,j+W,Criterion);

            if(CurrentDeviation < MinDeviation)
            {
                MinDeviation = CurrentDeviation;
                BestWedgeFilter = WedgeFilter;
            }

            WedgeFilter.MakeWedgeNW(Size);
            CurrentDeviation = MirroredImage.WedgeDeviation(WedgeFilter,i+H,j+W,Criterion);

            if(CurrentDeviation < MinDeviation)
            {
                MinDeviation = CurrentDeviation;
                BestWedgeFilter = WedgeFilter;
            }

            WedgeFilter.MakeWedgeNE(Size);
            CurrentDeviation = MirroredImage.WedgeDeviation(WedgeFilter,i+H,j+W,Criterion);

            if(CurrentDeviation < MinDeviation)
            {
                MinDeviation = CurrentDeviation;
                BestWedgeFilter = WedgeFilter;
            }

            WedgeFilter.MakeWedgeSE(Size);
            CurrentDeviation = MirroredImage.WedgeDeviation(WedgeFilter,i+H,j+W,Criterion);

            if(CurrentDeviation < MinDeviation)
            {
                MinDeviation = CurrentDeviation;
                BestWedgeFilter = WedgeFilter;
            }

            WedgeFilter.MakeWedgeSW(Size);
            CurrentDeviation = MirroredImage.WedgeDeviation(WedgeFilter,i+H,j+W,Criterion);

            if(CurrentDeviation < MinDeviation)
            {
                MinDeviation = CurrentDeviation;
                BestWedgeFilter = WedgeFilter;
            }

            //Apply the wedge filter to the current pixel of the image
            MirroredImage.ExtractBlock(Block,i,j,Size,Size);

            if(FilterType == 1) // 0  : box;
            {
                if(Size == 3)
                {
                    GaussFilter.MakeGaussianFilter3();
                }
                else if(Size == 5)
                {
                    GaussFilter.MakeGaussianFilter5();
                }
                else if(Size == 7)
                {
                    GaussFilter.MakeGaussianFilter7();
                }
                else
                {
                    GaussFilter.MakeGaussianFilter(sigma,Size);
                }

                for(int m=0;m<BestWedgeFilter.GetNbRows();m++)
                    for(int n=0;n<BestWedgeFilter.GetNbCols();n++)
                    {
                        BestWedgeFilter[m][n] = GaussFilter[m][n] * BestWedgeFilter[m][n];
                    }
            }

            BestWedgeFilter.Normalise();

            double FilteredValue = 0.;
            for(int m=0;m<BestWedgeFilter.GetNbRows();m++)
                for(int n=0;n<BestWedgeFilter.GetNbCols();n++)
                {
                    FilteredValue += MirroredImage[i+m][j+n] * BestWedgeFilter[m][n];
                }

            WedgeFilteredImage[i][j] = FilteredValue;
        }
    }
}


void M3Matrix::PixelTransfo(M3Matrix& TransformedImage, LUTParams& Parameters)
{
    TransformedImage.Reshape(t_NbRows,t_NbCols);
    TransformedImage.InitialiseDynamic((*this));

    //Read parameters
    double gamma = 1.0, slope = 1.0, shift = 0.;
    string opt;

    //Check type of transformation
    // 0: bypass --  1: gamma -- 2: linear multiplication
    // 3:gamnma for fpga -- 4: combined gamma and linear
    if(Parameters.getLUTType() == 0)
    {
        opt = "bypass";
    }
    else if(Parameters.getLUTType() == 1)
    {
        opt = "gamma";
        gamma = Parameters.gamma;
    }
    else if(Parameters.getLUTType() == 2)
    {
        opt = "linear";
        slope= Parameters.linear_slope;
        shift= Parameters.linear_shift;
    }
    else if(Parameters.getLUTType() == 3)
    {
        //Approximation for a fixed gamma = 0.65
        opt = "gamma_fpga";
    }
    else if(Parameters.getLUTType() == 4)
    {
        opt = "gamma_linear";
        slope= Parameters.linear_slope;
        shift= Parameters.linear_shift;
        gamma = Parameters.gamma;
    }
    else if(Parameters.getLUTType() == 5)
    {
        opt = "gamma_srgb";
    }
    else if (Parameters.getLUTType() == 6)
    {
        opt = "gamma_HDTV";
    }
    else if (Parameters.getLUTType() == 7)
    {
        opt = "linear_fpga";
        slope= Parameters.linear_slope;
        shift= Parameters.linear_shift;
    }

    double oldval,newval; //original value and gamma corrected value

    for(int i = 0;i<t_NbRows;i++)
        for(int j = 0;j<t_NbCols;j++)
        {
            oldval = (*this)[i][j]/(t_Max); //normalized auf [0...1]

            if(opt == "gamma")
            {
                newval = pow(fabs(oldval),gamma);

                //clip the values
                if(newval <0)
                    newval = 0;
                if(newval > 1)
                    newval = 1;

                if(oldval <0)
                    TransformedImage[i][j] = -newval*t_Max;
                else
                    TransformedImage[i][j] = newval*t_Max;
            }
            else if(opt == "linear")
            {
                newval = slope * oldval + shift;
                TransformedImage[i][j] = newval*t_Max;
            }
            else if(opt == "gamma_fpga")
            {
                //see in gamma_approx.m
                if(oldval <= 1./16.)
                    newval = 3.*oldval;
                else if(oldval <= 1./2.)
                    newval = oldval + 2./16.;
                else
                    newval =(6./8)*oldval + 4./16.;

                TransformedImage[i][j] = newval*t_Max;
            }
            else if(opt == "gamma_linear")
            {
                newval = sign(oldval) * pow(fabs(oldval),gamma);
                newval = slope * newval + shift;
                newval = myclip(newval,-1,1);

                TransformedImage[i][j] = newval*t_Max;
            }
            else if(opt == "gamma_srgb")
            {
                double oldval_abs = fabs(oldval);
                if(oldval_abs <= 0.0031308)
                    newval = 12.92 * oldval_abs;
                else
                    newval = (1+0.055) * pow(oldval_abs,1./2.4) - 0.055;

                newval = sign(oldval) * newval;
                newval = myclip(newval,-1,1);

                TransformedImage[i][j] = newval*t_Max;
            }
            else if (opt == "gamma_HDTV")
            {
                double oldval_abs = fabs(oldval);
                if (oldval_abs <= 0.018)
                    newval = 4.5068 * oldval_abs;
                else
                    newval = (1 + 0.09915) * pow(oldval_abs, 1. / 2.222) - 0.09915;

                newval = sign(oldval) * newval;
                newval = myclip(newval, -1, 1);

                TransformedImage[i][j] = newval*t_Max;
            }
            else if(opt == "bypass")
            {
                TransformedImage[i][j] = (*this)[i][j];
            }
        }


    if(opt == "linear_fpga")
    {
        int nbbits_frac = 8; //number of bits after the decimal point for the slope
        int64_t oldval_int, newval_int;
        int64_t slope_int = round(slope * (1 << nbbits_frac));
        int64_t offset_int = round(shift * (t_Max - t_Min));

        for(int i = 0;i<t_NbRows;i++)
            for(int j = 0;j<t_NbCols;j++)
            {
                oldval_int = (int64_t) ((*this)[i][j]);

                newval_int   = slope_int * oldval_int;
                newval_int  += offset_int << nbbits_frac;   // plus offset as Q-number
                newval_int  += 1 << (nbbits_frac - 1);      // plus "0.5" as Q-number
                newval_int >>= nbbits_frac;                 // round and transform Q-number to signed
                newval_int   = (int64_t) (myclip(newval_int, t_Min, t_Max));

                TransformedImage[i][j] = newval_int;
            }
    }
}


//First shoot, done without any LUT
void M3Matrix::GammaCorrection(M3Matrix& GammaCorrectedImage, double gamma, int Max)
{
    Max = Max;
    GammaCorrectedImage.Reshape(t_NbRows,t_NbCols);
    double val,gc_val; //original value and gamma corrected value

    for(int i = 0;i<t_NbRows;i++)
        for(int j = 0;j<t_NbCols;j++)
        {
            val = (*this)[i][j]/t_Max; //normalized auf [0...1]
            gc_val = pow(fabs(val),gamma);
            if(val < 0)
                GammaCorrectedImage[i][j] = -gc_val*t_Max;
            else
                GammaCorrectedImage[i][j] = gc_val*t_Max;
        }
}


void M3Matrix::AddGaussianShotNoise(double sigma, int seed)
{
    srand(seed);
    double rand_val;

    for(int i = 0;i<t_NbRows;i++)
        for(int j = 0;j<t_NbCols;j++)
        {
            rand_val = sigma * sqrt(fabs((*this)[i][j])/(double)t_Max) *randn();
            (*this)[i][j]  += rand_val;
        }
}


void M3Matrix::AddGaussianShotNoise(M3Matrix& NoisyImage, double sigma, int seed)
{
    NoisyImage.Reshape(t_NbRows,t_NbCols);
    srand(seed);
    double rand_val;

    for(int i = 0;i<t_NbRows;i++)
        for(int j = 0;j<t_NbCols;j++)
        {
            rand_val = sigma * sqrt(fabs((*this)[i][j])/(double)t_Max) *randn();
            NoisyImage[i][j] = (*this)[i][j] + rand_val;
        }
}


void M3Matrix::AddGaussianWhiteNoise(double sigma, int seed)
{
    srand(seed);
    double rand_val;

    for(int i = 0;i<t_NbRows;i++)
        for(int j = 0;j<t_NbCols;j++)
        {
            rand_val = sigma*randn();
            (*this)[i][j] += rand_val;
        }
}


void M3Matrix::AddGaussianWhiteNoise(M3Matrix& NoisyImage, double sigma, int seed)
{
    NoisyImage.Reshape(t_NbRows,t_NbCols);
    NoisyImage.InitialiseDynamic(*this);

    srand(seed);
    double rand_val;

    for(int i = 0;i<t_NbRows;i++)
        for(int j = 0;j<t_NbCols;j++)
        {
            rand_val = sigma*randn();
            NoisyImage[i][j] = (*this)[i][j] + rand_val;
        }
}


void M3Matrix::AddImpulsiveNoise(M3Matrix& NoisyImage, M3Matrix& ImpulsiveMask, double sigma, double prob, int seed)
{
    seed = seed;
    ImpulsiveMask.Reshape((*this).GetNbRows(),(*this).GetNbCols(),0.);
    NoisyImage=(*this);
    int NbDefectPixels = 0;
    double p1,rand_val;
    for(int i = 0;i<t_NbRows;i++)
        for(int j = 0;j<t_NbCols;j++)
        {
            p1 = rand_unif();
            if(p1 < prob)
            {
                rand_val = fabs(sigma*randn());
                NoisyImage[i][j] = (*this)[i][j] + rand_val;
                ImpulsiveMask[i][j] = 255;
                NbDefectPixels++;
            }
        }

    cout << "number of defect pixels generated in AddImpulsiveNoise"  << NbDefectPixels << endl;
}



void M3Matrix::AddImpulsiveNoise(double sigma, double prob, int seed)
{
    seed = seed;
    double p1,rand_val;
    for(int i = 0;i<t_NbRows;i++)
        for(int j = 0;j<t_NbCols;j++)
        {
            p1 = rand_unif();
            if(p1 < prob)
            {
                rand_val = sigma*randn();
                (*this)[i][j] += rand_val;
            }
        }
}


void M3Matrix::AddImpulsiveNoise(M3Matrix& NoisyImage, double sigma, double prob, int seed)
{
    seed = seed;
    double p1,rand_val;
    for(int i = 0;i<t_NbRows;i++)
        for(int j = 0;j<t_NbCols;j++)
        {
            p1 = rand_unif();
            if(p1 < prob)
            {
                rand_val = sigma*randn();
                NoisyImage[i][j] = (*this)[i][j] + rand_val;
            }
        }
}


void M3Matrix::Dump_From_Data()
{
    for(int i=0;i<t_NbRows*t_NbCols;i++)
    {
        cout <<" i  " << i << "data : "<< (*this)[i][0] << endl;
    }
}


void M3Matrix::AdditiveShadowing(double a)
{
    double x0=0.5;
    double fac=100;
    double cst=0;
    double jd,tmpval;

    for(int i = 0;i<t_NbRows;i++)
        for(int j = 0;j<t_NbCols;j++)
        {
            jd = ((double)j)/ ((double)t_NbCols);

            tmpval = cst -fac/2 +fac/(1 + exp(-a*(jd-x0)));
            (*this)[i][j] +=  tmpval;
        }
}


// Adds a shadow according to a sigmoid function defined by
// f(x) = cst -fac/2 +fac/(1 +exp(-a*x(i)-x0))
void M3Matrix::AdditiveShadowing(M3Matrix& Shadowed, double a)
{
    Shadowed.Reshape(t_NbRows,t_NbCols);
    Shadowed.InitialiseDynamic((*this));

    double x0=0.5;
    double fac=100;
    double cst=0;
    double jd,tmpval;

    for(int i = 0;i<t_NbRows;i++)
        for(int j = 0;j<t_NbCols;j++)
        {
            jd = ((double)j)/ ((double)t_NbCols);

            tmpval = cst -fac/2 +fac/(1 + exp(-a*(jd-x0)));
            Shadowed[i][j] = (*this)[i][j] + tmpval;
        }
}


// Adds a shadow according to a sigmoid function defined by
// f(x) = cst -fac/2 +fac/(1 +exp(-a*x(i)-x0))
void M3Matrix::MultiplicativeShadowing(M3Matrix& Shadowed, double a)
{
    Shadowed.Reshape(t_NbRows,t_NbCols);
    Shadowed.InitialiseDynamic((*this));

    double x0=0.5, fac=0.5, cst=0.5;
    double jd,tmpval;

    for(int i = 0;i<t_NbRows;i++)
        for(int j = 0;j<t_NbCols;j++)
        {
            jd = ((double)j)/ ((double)t_NbCols);

            tmpval = cst +fac/(1 + exp(-a*(jd-x0)));
            Shadowed[i][j] = (*this)[i][j] * tmpval;
        }
}


// Adds a shadow according to a sigmoid function defined by
// f(x) = cst -fac/2 +fac/(1 +exp(-a*x(i)-x0))
void M3Matrix::MultiplicativeShadowing(double a)
{
    double x0=0.5, fac=0.5, cst=0.5;
    double jd,tmpval;

    for(int i = 0;i<t_NbRows;i++)
        for(int j = 0;j<t_NbCols;j++)
        {
            jd = ((double)j)/ ((double)t_NbCols);

            tmpval = cst +fac/(1 + exp(-a*(jd-x0)));
            (*this)[i][j] *=  tmpval;
        }
}


void M3Matrix::SlantedEdge(double theta)
{
    double id,jd,id_r,jd_r,id_theta,jd_theta;

    for(int i = 0;i<t_NbRows;i++)
        for(int j = 0;j<t_NbCols;j++)
        {
            id = ((double)i)/ ((double)t_NbRows);
            jd = ((double)j)/ ((double)t_NbCols);

            id_r = id-0.5;
            jd_r = jd-0.5;

            id_theta = cos(theta)*id_r + sin(theta)*jd_r;
            jd_theta = -sin(theta)*id_r + cos(theta)*jd_r;

            if((id_theta >0 && jd_theta>0) || (id_theta<0 && jd_theta<0))
                (*this)[i][j] = 192;
            else
                (*this)[i][j] = 64;
        }
}


void M3Matrix::DrawRandomCrackTree(double xStart, double yStart, double Value)
{
    double LineWidth=0.001;
    double Length=0.07;

    double theta = 3.14159/2.+0.2*randn();
    (*this).DrawLine(xStart,yStart,Value,theta,Length,LineWidth);

    xStart += Length*cos(theta);
    yStart += Length*sin(theta);

    //Two lines after the first junction
    double theta_1p = theta + 0.35;// + 0.2*randn();
    (*this).DrawLine(xStart,yStart,Value,theta_1p,Length,LineWidth);
    double theta_1n = theta - 0.55;// - 0.2*randn();
    (*this).DrawLine(xStart,yStart,Value,theta_1n,Length,LineWidth);

    //Two lines after the second junction pos
    double xStart_1p =   xStart + Length*cos(theta_1p);
    double yStart_1p =   yStart + Length*sin(theta_1p);

    theta_1p += theta_1p + 0.25;

    (*this).DrawLine(xStart_1p,yStart_1p,Value,theta_1p,Length,LineWidth);
    xStart_1p += Length*cos(theta_1p);
    yStart_1p += Length*sin(theta_1p);

    double theta_1p_2p  = theta_1p + 0.25 + 1*randn();

    (*this).DrawLine(xStart_1p,yStart_1p,Value,theta_1p_2p,Length,LineWidth);
    double theta_1p_2n = theta_1p - 0.05-1*randn();
    (*this).DrawLine(xStart_1p,yStart_1p,Value,theta_1p_2n,Length,LineWidth);

    //Two lines after the second junction neg
    double xStart_1n =   xStart + Length*cos(theta_1n);
    double yStart_1n =   yStart + Length*sin(theta_1n);
    double theta_1n_2p = theta_1n + 0.05 + 1*randn();
    (*this).DrawLine(xStart_1n,yStart_1n,Value,theta_1n_2p,Length,LineWidth);
    double theta_1n_2n = theta_1n - 0.05-1*randn();
    (*this).DrawLine(xStart_1n,yStart_1n,Value,theta_1n_2n,Length,LineWidth);
}


void M3Matrix::DrawCrown(double r1, double r2, double cx, double cy, double Value)
{
    for(int i=0;i<GetNbRows();i++)
        for(int j=0;j<GetNbCols();j++)
        {
            double id = (double)i/((double) GetNbRows());
            double jd = (double)j/((double) GetNbCols());
            if((id-cx)*(id-cx)+(jd-cy)*(jd-cy)< r2*r2  && (id-cx)*(id-cx)+(jd-cy)*(jd-cy) > r1*r1)
                (*this)[i][j] = Value;
        }
}


void M3Matrix::DrawCircle(double r, double cx, double cy, double Value)
{
    for(int i=0;i<GetNbRows();i++)
        for(int j=0;j<GetNbCols();j++)
        {
            double id = (double)i/((double) GetNbRows());
            double jd = (double)j/((double) GetNbCols());
            if((id-cx)*(id-cx)+(jd-cy)*(jd-cy)< r*r )
                (*this)[i][j] = Value;
        }
}


void M3Matrix::AddBlob(double r, double cx, double cy, double Value)
{
    for(int i=0;i<GetNbRows();i++)
        for(int j=0;j<GetNbCols();j++)
        {
            double id = (double)i/((double) (GetNbRows()-1));
            double jd = (double)j/((double) (GetNbCols()-1));
            double dist_to_center = sqrt((id-cx)*(id-cx)+(jd-cy)*(jd-cy));
            if((id-cx)*(id-cx)+(jd-cy)*(jd-cy)< r*r )
                (*this)[i][j] += Value*fabs(1-dist_to_center);
        }
}


void M3Matrix::AddGaussianBlob(double cx, double cy, double sigma, double MaxValue)
{
    for(int i=0;i<GetNbRows();i++)
        for(int j=0;j<GetNbCols();j++)
        {
            double id = (double)i/((double) (GetNbRows()-1));
            double jd = (double)j/((double) (GetNbCols()-1));
            double dist2_to_center = (id-cx)*(id-cx)+(jd-cy)*(jd-cy);

            (*this)[i][j] += exp(-dist2_to_center/(sigma*sigma))*MaxValue;
        }
}


//Draws an elliptic crown in the image
void M3Matrix::DrawEllipticCrown(double r1, double r2, double cx, double cy, double a, double b, double theta, double Value)
{

    for(int i=0;i<GetNbRows();i++)
        for(int j=0;j<GetNbCols();j++)
        {
            double id = (double)i/((double) GetNbRows());
            double jd = (double)j/((double) GetNbCols());

            double xd = cos(theta)*(id-cx) + sin(theta)*(jd-cy);
            double yd = -sin(theta)*(id-cx) + cos(theta)*(jd-cy);

            if(xd*xd/(a*a) + yd*yd/(b*b) < r2*r2 && xd*xd/(a*a) + yd*yd/(b*b) > r1*r1)
                (*this)[i][j] = Value;
        }
}


//Draws an ellipse in the image
void M3Matrix::DrawEllipse(double r, double cx, double cy, double a, double b, double theta, double Value)
{
    for(int i=0;i<GetNbRows();i++)
        for(int j=0;j<GetNbCols();j++)
        {
            double id = (double)i/((double) GetNbRows());
            double jd = (double)j/((double) GetNbCols());

            double xd = cos(theta)*(id-cx) + sin(theta)*(jd-cy);
            double yd = -sin(theta)*(id-cx) + cos(theta)*(jd-cy);

            if(xd*xd/(a*a) + yd*yd/(b*b) < r*r )
                (*this)[i][j] = Value;
        }
}


//ulc: Upper Left Corner; lrc: Lower Right Corner
void M3Matrix::DrawRectangle(double ulcx, double ulcy, double lrcx, double lrcy, double Value)
{
    for(int i=0;i<GetNbRows();i++)
        for(int j=0;j<GetNbCols();j++)
        {
            double id = (double)i/((double) t_NbRows);
            double jd = (double)j/((double) t_NbCols);
            if(id > ulcx && id < lrcx && jd > ulcy && jd < lrcy)
                (*this)[i][j] = Value;
        }
}


//Draws a line
void M3Matrix::DrawLine(double xStart,double yStart,double Value, double theta, double Length,double Width)
{
    double id,jd,scal_prod,lat_prod;

    for(int i=0;i<GetNbRows();i++)
        for(int j=0;j<GetNbCols();j++)
        {
            id = (double)i/((double) GetNbRows());
            jd = (double)j/((double) GetNbCols());

            scal_prod =  (id-xStart)*cos(theta) +(jd-yStart)*sin(theta);
            lat_prod =  -(id-xStart)*sin(theta) +(jd-yStart)*cos(theta);
            if(scal_prod >0 && scal_prod <Length )
            {
                if(fabs(lat_prod) < Width)
                    (*this)[i][j] = Value;
            }
        }
}


void M3Matrix::DrawDiagonalBand(double xStart, double Width, double Value, double Slope)
{
    //Slope = 2;
    for(int i=0;i<GetNbRows();i++)
        for(int j=0;j<GetNbCols();j++)
        {
            double id = (double)i/((double) GetNbRows());
            double jd = (double)j/((double) GetNbCols());

            if(id-jd/Slope > xStart && id-jd/Slope < xStart+Width)
                (*this)[i][j] = Value;
        }
}


void M3Matrix::AddCapillary(double Factor, double xShift,double yShift)
{
    M3Matrix Mask(t_NbRows,t_NbCols);
    Mask.SetValues(0.);

    //Original parameters
    // double fx = 7;   double fy = 8.5;   double fz = 2;
    //Some test parameters
    double fx = 7.8;
    double fy = 23;
    double fz = 0.5;

    int NbSteps = 46500;
    for(int t=0; t < NbSteps; t++)
    {
        double td = t*2*6.3/((double)NbSteps);

        double xt = 0.3*sin(fx*td)+0.5;
        double yt = 0.3*sin(fy*td)+0.5;
        double zt = (128+13*cos(fz*td));

        int xd = (int)(t_NbRows*(xt-xShift));
        int yd = (int)(t_NbCols*(yt-yShift));

        if(xd >=0 &&  xd < t_NbRows && yd >=0 && yd < t_NbCols)
        {
            Mask[xd][yd] = (int)(zt+0.5) - 128;
        }
    }

    for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols;j++)
        {
            (*this)[i][j] += Factor*Mask[i][j];
        }
}


void M3Matrix::DrawCapillary(double Value)
{
    M3Matrix Mask(t_NbRows,t_NbCols);
    Mask.SetValues(0.);

    int NbSteps = 24500;
    for(int t=0; t < NbSteps; t++)
    {
        double td = t*2*6.3/((double)NbSteps);

        double xt = 0.3*sin(7*td)+0.5;
        double yt = 0.3*sin(8.5*td)+0.5;
        double zt = 128+13*cos(2*td);

        int xd = (int)(t_NbRows*xt);
        int yd = (int)(t_NbCols*yt);

        if(xd >=0 &&  xd < t_NbRows && yd >=0 && yd < t_NbCols)
        {
            Mask[xd][yd] = (int)(zt+0.5);
        }
    }

    for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols;j++)
        {
            if(Mask[i][j]>0)
            {
                (*this)[i][j] = Value;
                (*this)[i][j] = Mask[i][j];
            }
        }
}


//*************************************
//* Negative : used mainly for testing
//*************************************
void M3Matrix::Negative(M3Matrix&  Negative)
{
    Negative.Reshape(GetNbRows(),GetNbCols());
    Negative = t_Max -(*this);
}


//Left multiplication
M3Matrix operator*(const double& scalar, M3Matrix A)
{
    M3Matrix B(A.GetNbRows(),A.GetNbCols());
    for(int i = 0;i<A.GetNbRows();i++)
        for(int j = 0;j<A.GetNbCols();j++)
            B[i][j] = A[i][j] * scalar;
    return B;
}


//Left addition
M3Matrix operator+(const double& scalar, M3Matrix A)
{
    M3Matrix B(A.GetNbRows(),A.GetNbCols());
    for(int i = 0;i<A.GetNbRows();i++)
        for(int j = 0;j<A.GetNbCols();j++)
            B[i][j] = A[i][j] + scalar;
    return B;
}


//Left substraction
M3Matrix operator-(const double& scalar, M3Matrix A)
{
    M3Matrix B(A.GetNbRows(),A.GetNbCols());
    for(int i = 0;i<A.GetNbRows();i++)
        for(int j = 0;j<A.GetNbCols();j++)
            B[i][j] = scalar - A[i][j] ;
    return B;
}


M3Matrix Diag(vector<double> d)
{
    M3Matrix D(d.size(),d.size());

    for(int i = 0;i<(int)d.size();i++)
        D[i][i] = d[i];

    return D;
}


//Taken from the MATLAB function hannPadded.m
M3Matrix HannPadded(int m, int n, double factor)
{
    M3Matrix WH(m,n);
    int add_m = floor(factor*m);
    int add_n = floor(factor*n);

    if ( add_m % 2 != 0)
        add_m -= 1;

    if ( add_n % 2 != 0)
        add_n -= 1;

    int nn = n + add_n;
    int mm = m + add_m;

    double cosx, cosy;
    int ii,jj;
    for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
        {
            ii = i + add_m/2;
            jj = j + add_n/2;
            cosx = cos(2*PI*jj/(nn-1));
            cosy = cos(2*PI*ii/(mm-1));
            WH[i][j] = 0.25*(1-cosx-cosy+cosx*cosy);
        }

    return WH;
}


//Pointwise product of two matrices
M3Matrix HadamardProduct(M3Matrix& A, M3Matrix& B)
{
    M3Matrix AB;
    int NbRows1 = A.GetNbRows(), NbCols1= A.GetNbCols();
    int NbRows2 = B.GetNbRows(), NbCols2= B.GetNbCols();

    if(NbRows1 == NbRows2 && NbCols1 == NbCols2)
    {
        AB.Reshape(NbRows1,NbCols1);
        for(int i=0;i<NbRows1;i++)
          for(int j=0;j<NbCols1;j++)
          {
              AB[i][j] = A[i][j] * B[i][j];
          }
    }

    return AB;
}


//Pointwise product of two matrices
//Warning: this is the efficient, not memory controlled version !!!
void HadamardProduct_NoAlloc(M3Matrix& A, M3Matrix& B, M3Matrix& C)
{
  for(int i=0;i<A.GetNbRows();i++)
    for(int j=0;j<A.GetNbCols();j++)
    {
      C[i][j] = A[i][j] * B[i][j];
    }
}


M3Matrix Transpose(M3Matrix& A)
{
  M3Matrix AT;
  int NbRows = A.GetNbRows(), NbCols= A.GetNbCols();

  AT.Reshape(NbCols,NbRows);
  AT.InitialiseDynamic(A);
  for(int i=0;i<NbRows;i++)
    for(int j=0;j<NbCols;j++)
    {
      AT[j][i] = A[i][j];
    }

  return AT;
}


//*******************************************
//* Useful numerical tools
//*******************************************
void M3Matrix::Hessian(M3Matrix3D& H, double sigma)
{
  //H[0]: is a matrix containing H_xx,H[1]: H_yy and H[2]: H_xy
  H.Reshape(3,t_NbRows,t_NbCols);
  int HalfSize = round(3 *sigma), KernelSize = 2*HalfSize + 1,i = 0, j = 0; //Inspired from the opencv Frangi filters

  //Creates the hessian kernels
  M3Matrix HK_xx(KernelSize,KernelSize), HK_yy(KernelSize,KernelSize), HK_xy(KernelSize,KernelSize);
  for (int x = -HalfSize; x <= HalfSize; x++)
  {
    j = 0;
    for (int y = -HalfSize; y <= HalfSize; y++)
    {
      HK_xx[i][j] = 1.0f / (2.0f*M_PI*sigma*sigma*sigma*sigma) * (x*x / (sigma*sigma) - 1) * exp(-(x*x + y*y) / (2.0f*sigma*sigma));
      HK_xy[i][j] = 1.0f / (2.0f*M_PI*sigma*sigma*sigma*sigma*sigma*sigma)*(x*y)*exp(-(x*x + y*y) / (2.0f*sigma*sigma));
      j++;
    }
    i++;
  }

  HK_xx.Transpose(HK_yy);

  //Filtering
  (*this).LinearFilter(H[0],HK_xx);
  (*this).LinearFilter(H[1],HK_yy);
  (*this).LinearFilter(H[2],HK_xy);
}


//******************************************
// WAVELETS TOOLS
//******************************************
// Horizontal undecimated wavelet filtering step  of an image (downsampling included)
void M3Matrix::UndecimatedWaveletFilterRow(M3Matrix& SubbandMatrix, M3Matrix& LowFilter, M3Matrix& HighFilter)
{
    //Mirroring
    M3Matrix EnlargedMatrix;
    (*this).Extend(EnlargedMatrix,0,LowFilter.GetNbRows()/2); //TODO: change this to some mirroring ?

    M3Matrix LowFilterRev,HighFilterRev;
    LowFilter.Transpose(LowFilterRev);
    LowFilterRev.GetRevertCols(LowFilterRev);
    HighFilter.Transpose(HighFilterRev);
    HighFilterRev.GetRevertCols(HighFilterRev);

    M3Matrix LowPassImage, HighPassImage;
    (*this).LinearFilter(LowPassImage,LowFilterRev);
    (*this).LinearFilter(HighPassImage,HighFilterRev);
    SubbandMatrix.CatMat(LowPassImage,HighPassImage,2);
}


// Wavelet filtering of an image according to the rows
void M3Matrix::UndecimatedWaveletFilterCol(M3Matrix& SubbandMatrix,M3Matrix& LowFilter,M3Matrix& HighFilter)
{
    //Mirrorisation
    M3Matrix EnlargedMatrix;
    (*this).Extend(EnlargedMatrix,LowFilter.GetNbRows()/2,0); //TODO: change this to some mirroring ?

    M3Matrix LowFilterRev,HighFilterRev;
    LowFilterRev.GetRevertRows(LowFilter);
    HighFilterRev.GetRevertRows(HighFilter);

    M3Matrix LowPassImage,HighPassImage;
    (*this).LinearFilter(LowPassImage,LowFilterRev);
    (*this).LinearFilter(HighPassImage,HighFilterRev);
    SubbandMatrix.CatMat(LowPassImage,HighPassImage,1);
}


// Horizontal wavelet filtering step  of an image (downsampling included)
void M3Matrix::WaveletFilterRow(M3Matrix& SubbandMatrix,M3Matrix& LowFilter, M3Matrix& HighFilter)
{
    //Mirroring
    M3Matrix EnlargedMatrix;
    (*this).Extend(EnlargedMatrix,0,LowFilter.GetNbRows()/2);

    M3Matrix LowFilterRev,HighFilterRev;
    LowFilter.Transpose(LowFilterRev);
    LowFilterRev.GetRevertCols(LowFilterRev);
    HighFilter.Transpose(HighFilterRev);
    HighFilterRev.GetRevertCols(HighFilterRev);

    M3Matrix LowPassImage(t_NbRows,t_NbCols/2),HighPassImage(t_NbRows,t_NbCols/2);
    (*this).LinearFilter_Hor_DS(LowPassImage,LowFilterRev);
    (*this).LinearFilter_Hor_DS(HighPassImage,HighFilterRev);

    SubbandMatrix.CatMat(LowPassImage,HighPassImage,2);
}


// Wavelet filtering of an image according to the rows
void M3Matrix::WaveletFilterCol(M3Matrix& SubbandMatrix,M3Matrix& LowFilter,M3Matrix& HighFilter)
{
    //Mirrorisation
    M3Matrix EnlargedMatrix;
    (*this).Extend(EnlargedMatrix,LowFilter.GetNbRows()/2,0);

    M3Matrix LowFilterRev,HighFilterRev;
    LowFilterRev.GetRevertRows(LowFilter);
    HighFilterRev.GetRevertRows(HighFilter);

    M3Matrix LowPassImage(t_NbRows/2,t_NbCols),HighPassImage(t_NbRows/2,t_NbCols);
    (*this).LinearFilter_Ver_DS(LowPassImage,LowFilterRev);
    (*this).LinearFilter_Ver_DS(HighPassImage,HighFilterRev);

    SubbandMatrix.CatMat(LowPassImage,HighPassImage,1);
}


void M3Matrix::InverseUndecimatedWaveletFilterRow(M3Matrix& SubbandMatrix,M3Matrix& LowFilter,M3Matrix& HighFilter)
{
    int LowFSize = LowFilter.GetNbRows();
    int HighFSize = HighFilter.GetNbRows();
    SubbandMatrix.Reshape(t_NbRows,t_NbCols/2);

    M3Matrix LMatrix(t_NbRows,t_NbCols/2), HMatrix(t_NbRows,t_NbCols/2);
    int HalfWidth = t_NbCols/2;

    //No Upsampling
    for(int Row = 0;Row < t_NbRows ;Row++)
        for(int Col = 0;Col < t_NbCols/2;Col++)
        {
            LMatrix[Row][Col] = (*this)[Row][Col];
            HMatrix[Row][Col] = (*this)[Row][Col + HalfWidth];
        }

    //Mirroring
    M3Matrix LEnlargedMatrix, HEnlargedMatrix;
    LMatrix.Extend(LEnlargedMatrix,0,LowFSize/2);
    HMatrix.Extend(HEnlargedMatrix,0,HighFSize/2);

    //Wavelet inverse filtering
    double tmp1,tmp2,tmp3,tmp4;
    for(int i=0;i < t_NbRows;i++)
        for(int j=1;j <= t_NbCols/2;j++)
        {
            tmp1 = 0.,tmp2 = 0.,tmp3 = 0.,tmp4 = 0.;

            if(LowFSize == 2 && HighFSize == 2) //this case should be only used for Haar wavelets
            {
                tmp1 =  LowFilter[0][0] * LEnlargedMatrix[i][j];
                tmp2 =  LowFilter[1][0] * LEnlargedMatrix[i][j+1];
                tmp3 = HighFilter[0][0] * HEnlargedMatrix[i][j];
                tmp4 = HighFilter[1][0] * HEnlargedMatrix[i][j+1];

                if(j<t_NbCols/2)
                    SubbandMatrix[i][j-1] = 0.5*(tmp1 + tmp2 + tmp3 +tmp4); // + tmp3 +tmp4); //the factor 0.5 is due to the oversampling while using undecimated wavelet dec.
                else
                    SubbandMatrix[i][j-1] = (tmp1 + tmp3); // + tmp3 +tmp4); //the factor 0.5 is due to the oversampling while using undecimated wavelet dec.

            }
            else
            {
                for(int k=0;k<LowFSize;k++)
                {
                    tmp1 += LowFilter[k][0] * LEnlargedMatrix[i][j+2*(LowFSize/2)-k];
                }

                for(int k=0;k<HighFSize;k++)
                {
                    tmp3 += HighFilter[k][0] * HEnlargedMatrix[i][j+2*(HighFSize/2)-k];
                }

                SubbandMatrix[i][j-1] = 0.5*(tmp1 + tmp3); // + tmp3 +tmp4); //the factor 0.5 is due to the oversampling while using undecimated wavelet dec.
            }
        }
}


void M3Matrix::InverseUndecimatedWaveletFilterCol(M3Matrix& SubbandMatrix,M3Matrix& LowFilter,M3Matrix& HighFilter)
{
    int LowFSize = LowFilter.GetNbRows();
    int HighFSize = HighFilter.GetNbRows();
    SubbandMatrix.Reshape(t_NbRows/2,t_NbCols);

    M3Matrix LMatrix(t_NbRows/2,t_NbCols), HMatrix(t_NbRows/2,t_NbCols);
    int HalfHeight = t_NbRows/2;

    //no upsampling
    for(int Row = 0;Row < t_NbRows/2;Row++)
        for(int Col = 0;Col < t_NbCols;Col++)
        {
            LMatrix[Row][Col] = (*this)[Row][Col];
            HMatrix[Row][Col] = (*this)[Row + HalfHeight][Col];
        }

    //Mirroring
    M3Matrix LEnlargedMatrix, HEnlargedMatrix;
    LMatrix.Extend(LEnlargedMatrix,LowFSize/2+1,0);
    HMatrix.Extend(HEnlargedMatrix,HighFSize/2+1,0);

    //Wavelet inverse filtering
    double tmp1,tmp2, tmp3,tmp4;
    for(int i=1; i <= t_NbRows/2;i++)
        for(int j=0;j < t_NbCols;j++)
        {
            tmp1 = 0.,tmp2 = 0.,tmp3 = 0.,tmp4 = 0.;

            // Haar wavelets
            if(LowFSize == 2 && HighFSize == 2)
            {
                tmp1 =  LowFilter[0][0] * LEnlargedMatrix[i][j];
                tmp2 =  LowFilter[1][0] * LEnlargedMatrix[i+1][j];
                tmp3 = HighFilter[0][0] * HEnlargedMatrix[i][j];
                tmp4 = HighFilter[1][0] * HEnlargedMatrix[i+1][j];
            }
            else       // other than Haar wavelets
            {
                for(int k=0;k<LowFSize;k++)
                {
                    tmp1 += LowFilter[k][0]* LEnlargedMatrix[i+1+2*(LowFSize/2)-k][j];
                }
                for(int k=0;k<HighFSize;k++)
                {
                    tmp3 += HighFilter[k][0]* HEnlargedMatrix[i+1+2*(HighFSize/2)-k][j];
                }
            }

            if(LowFSize == 2 && HighFSize == 2)
            {
                if(i<t_NbRows/2)
                    SubbandMatrix[i-1][j] = 0.5*(tmp1 + tmp2 + tmp3 +tmp4); // + tmp3 +tmp4); //the factor 0.5 is due to the oversampling while using undecimated wavelet dec.
                else
                    SubbandMatrix[i-1][j] = (tmp1 + tmp3); // + tmp3 +tmp4);
            }
            else
            {
                SubbandMatrix[i-1][j] = 0.5*(tmp1 + tmp3); // + tmp3 +tmp4);
            }
        }
}


// Wavelet filtering of an image according to the rows
void M3Matrix::InverseWaveletFilterRow(M3Matrix& SubbandMatrix,M3Matrix& LowFilter,M3Matrix& HighFilter)
{
    int LowFSize = LowFilter.GetNbRows();
    int HighFSize = HighFilter.GetNbRows();
    SubbandMatrix.Reshape(t_NbRows,t_NbCols);

    M3Matrix LMatrix(t_NbRows,t_NbCols), HMatrix(t_NbRows,t_NbCols);
    int HalfWidth = t_NbCols/2;

    //Upsampling
    for(int Row = 0;Row < t_NbRows ;Row++)
        for(int Col = 0;Col < HalfWidth ;Col++)
        {
            LMatrix[Row][2*Col] = (*this)[Row][Col];
            LMatrix[Row][2*Col+1] = 0.;

            HMatrix[Row][2*Col] = (*this)[Row][Col + HalfWidth];
            HMatrix[Row][2*Col+1] = 0.;
        }

    //Mirroring
    M3Matrix LEnlargedMatrix, HEnlargedMatrix;
    LMatrix.Extend(LEnlargedMatrix,0,LowFSize/2);
    HMatrix.Extend(HEnlargedMatrix,0,HighFSize/2);

    //Wavelet inverse filtering
    double tmp1,tmp2;
    for(int i=0;i < t_NbRows;i++)
        for(int j=0;j < t_NbCols;j++)
        {
            tmp1 = 0., tmp2 = 0.;

            if(LowFSize == 2 && HighFSize == 2) //this case should be only used for Haar wavelets
            {
                if(j%2==0)
                {
                    tmp1 =  LowFilter[0][0] * LMatrix[i][j] + LowFilter[1][0] * LMatrix[i][j+1];
                    tmp2 = HighFilter[0][0] * HMatrix[i][j+1] + HighFilter[1][0] * HMatrix[i][j];
                }
                if(j%2==1)
                {
                    tmp1 =  LowFilter[0][0] * LMatrix[i][j] +  LowFilter[1][0] * LMatrix[i][j-1];
                    tmp2 = HighFilter[0][0] * HMatrix[i][j-1] + HighFilter[1][0] * HMatrix[i][j];
                }
            }
            else
            {
                for(int k=0;k<LowFSize;k++)
                    tmp1 += LowFilter[k][0] * LEnlargedMatrix[i][j+2*(LowFSize/2)-k];

                for(int k=0;k<HighFSize;k++)
                    tmp2 += HighFilter[k][0] * HEnlargedMatrix[i][j+2*(HighFSize/2)-k];
            }
            SubbandMatrix[i][j] = tmp1 + tmp2;
        }
}


// Wavelet filtering of an image according to the cols
void M3Matrix::InverseWaveletFilterCol(M3Matrix& SubbandMatrix,M3Matrix& LowFilter,M3Matrix& HighFilter)
{
    int LowFSize = LowFilter.GetNbRows();
    int HighFSize = HighFilter.GetNbRows();

    SubbandMatrix.Reshape(t_NbRows,t_NbCols);

    M3Matrix LMatrix(t_NbRows,t_NbCols);
    M3Matrix HMatrix(t_NbRows,t_NbCols);
    int HalfHeight = t_NbRows/2;

    //Upsampling operator
    for(int Row = 0;Row < HalfHeight ;Row++)
        for(int Col = 0;Col < t_NbCols;Col++)
        {
            LMatrix[2*Row][Col] = (*this)[Row][Col];
            LMatrix[2*Row+1][Col] = 0.;

            HMatrix[2*Row][Col] = (*this)[Row + HalfHeight][Col];
            HMatrix[2*Row+1][Col] = 0.;
        }

    //Mirroring
    M3Matrix LEnlargedMatrix, HEnlargedMatrix;
    LMatrix.Extend(LEnlargedMatrix,LowFSize/2,0);
    HMatrix.Extend(HEnlargedMatrix,HighFSize/2,0);

    //Wavelet inverse filtering
    double tmp1,tmp2;
    for(int i=0;i < t_NbRows;i++)
        for(int j=0;j < t_NbCols;j++)
        {
            tmp1 = 0., tmp2 = 0.;

            if(LowFSize == 2 && HighFSize == 2) //this case should be only used for Haar wavelets
            {
                if(i%2==0)
                {
                    tmp1 =  LowFilter[0][0] * LMatrix[i][j] + LowFilter[1][0] * LMatrix[i+1][j];
                    tmp2 = HighFilter[0][0] * HMatrix[i+1][j] + HighFilter[1][0] * HMatrix[i][j];
                }
                if(i%2==1)
                {
                    tmp1 =  LowFilter[0][0] * LMatrix[i][j] + LowFilter[1][0] * LMatrix[i-1][j];
                    tmp2 = HighFilter[0][0] * HMatrix[i-1][j] + HighFilter[1][0] * HMatrix[i][j];
                }
            }
            else
            {
                for(int k=0;k<LowFSize;k++)
                {
                    tmp1 += LowFilter[k][0]* LEnlargedMatrix[i+2*(LowFSize/2)-k][j];
                }
                for(int k=0;k<HighFSize;k++)
                {
                    tmp2 += HighFilter[k][0]* HEnlargedMatrix[i+2*(HighFSize/2)-k][j];
                }
            }
            SubbandMatrix[i][j] = tmp1 + tmp2;
        }
}


void M3Matrix::InverseUndecimatedWaveletDecomposition(M3Matrix3D& USubband,M3Matrix& LowFilter,M3Matrix& HighFilter,int NbDecomp)
{
    int NbRows =  USubband[0].GetNbRows();
    int NbCols =  USubband[0].GetNbCols();

    M3Matrix TmpLL,TmpHH,TmpHL,TmpLH,TmpT,TmpB,Tmp, Tmp2;
    TmpLL = USubband[0];

    for(int step = 1;step <= NbDecomp;step++)
    {
        USubband[step].ExtractBlock(TmpLH,0,0,NbRows,NbCols);
        USubband[step].ExtractBlock(TmpHH,0,NbCols,NbRows,NbCols);
        USubband[step].ExtractBlock(TmpHL,0,2*NbCols,NbRows,NbCols);

        TmpT.CatMat(TmpLL,TmpHL); //concatenation of the two top images
        TmpB.CatMat(TmpLH,TmpHH); //concatenation of the two bottom images
        Tmp.CatMat(TmpT,TmpB,1);

        Tmp.InverseUndecimatedWaveletFilterCol(Tmp2,LowFilter,HighFilter);
        Tmp2.InverseUndecimatedWaveletFilterRow(TmpLL,LowFilter,HighFilter);
    }

    (*this) = TmpLL;
}


//Redundant wavelet decomposition (without downsampling)
void M3Matrix::UndecimatedWaveletDecomposition(M3Matrix3D& USubband,M3Matrix& LowFilter,M3Matrix& HighFilter,int NbDecomp)
{
    int NbRows =  t_NbRows;
    int NbCols =  t_NbCols;

    USubband.Reshape(NbDecomp+1);

    M3Matrix Tmp, Tmp2, Tmp3, Tmp4, Tmp5, TmpLL, TmpHL,TmpLH,TmpHH;
    Tmp = (*this);

    for(int step = 1;step <= NbDecomp;step++)
    {
        Tmp.UndecimatedWaveletFilterRow(Tmp2,LowFilter,HighFilter);
        Tmp2.UndecimatedWaveletFilterCol(Tmp3,LowFilter,HighFilter);

        Tmp3.ExtractBlock(TmpHH,NbRows,NbCols,NbRows,NbCols);
        Tmp3.ExtractBlock(TmpHL,0,NbCols,NbRows,NbCols);
        Tmp3.ExtractBlock(TmpLH,NbRows,0,NbRows,NbCols);

        Tmp4.CatMat(TmpLH,TmpHH);
        Tmp5.CatMat(Tmp4,TmpHL);
        USubband[NbDecomp-step+1].Reshape(Tmp5.GetNbRows(),Tmp5.GetNbCols());
        USubband[NbDecomp-step+1] = Tmp5;

        Tmp3.ExtractBlock(Tmp,0,0,NbRows,NbCols);
    }

    USubband[0] = Tmp;
}


// Wavelet decomposition
void M3Matrix::WaveletDecomposition(M3Matrix& SubbandMatrix,M3Matrix& LowFilter,M3Matrix& HighFilter,int NbDecomp)
{
    int NbRows =  t_NbRows;
    int NbCols =  t_NbCols;

    M3Matrix TmpMatrix,Tmp2Matrix;
    SubbandMatrix = (*this);

    for(int step = 1;step <= NbDecomp;step++)
    {
        cout << "step : " << step << " Nb Rows : " << NbRows << endl;
        SubbandMatrix.ExtractBlock(TmpMatrix, 0, 0, NbRows, NbCols);

        TmpMatrix.WaveletFilterRow(Tmp2Matrix,LowFilter,HighFilter);
        Tmp2Matrix.WaveletFilterCol(TmpMatrix,LowFilter,HighFilter);

        for(int i = 0;i < NbRows;i++)
            for(int j = 0;j < NbCols;j++)
                SubbandMatrix[i][j] = TmpMatrix[i][j];

        NbRows /= 2;
        NbCols /= 2;
    }
}


// Inverse Wavelet decomposition
// Param : ro_SubbandVisualisation, made         from the ro_subbandMatrix
void M3Matrix::InverseWaveletDecomposition(M3Matrix& SubbandMatrix,M3Matrix& LowFilter,M3Matrix& HighFilter,int NbDecomp)
{
    int NbRows =  SubbandMatrix.GetNbRows();
    int NbCols =  SubbandMatrix.GetNbCols();

    //SubbandVisualisationImage.Reshape(NbRows,NbCols);
    M3Matrix TmpMatrix,Tmp2Matrix;
    this->Reshape(NbRows,NbCols);

    NbRows =  NbRows/((int)pow(2.0,(double)(NbDecomp-1)));
    NbCols =  NbCols/((int)pow(2.0,(double)(NbDecomp-1)));

    for(int step = 1;step <= NbDecomp;step++)
    {
        TmpMatrix.Reshape(NbRows,NbCols);
        for(int i = 0;i < NbRows;i++)
            for(int j = 0;j < NbCols;j++)
                TmpMatrix[i][j] = SubbandMatrix[i][j];

        TmpMatrix.InverseWaveletFilterCol(Tmp2Matrix,LowFilter,HighFilter);
        Tmp2Matrix.InverseWaveletFilterRow(TmpMatrix,LowFilter,HighFilter);

        for(int i = 0;i < NbRows;i++)
            for(int j = 0;j < NbCols;j++)
                SubbandMatrix[i][j] = TmpMatrix[i][j];

        NbRows *= 2;
        NbCols *= 2;
    }

    NbRows /= 2;
    NbCols /= 2;

    for(int i = 0;i < NbRows;i++)
        for(int j = 0;j < NbCols;j++)
            (*this)[i][j] = TmpMatrix[i][j];
}


void M3Matrix::WaveletEnhancement(double AmplificationFactor, int NbDec)
{
    //First (06 march 2015), naive version
    int dyadic = (1<<NbDec);
    double factor = 1;

    for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols;j++)
        {
            int ij=i/(t_NbRows/dyadic);
            if(j>i) {
                ij=j/(t_NbCols/dyadic);
            }

            int ij_d = int(log(double(ij))/log(2))+1;


            if(ij == 0)
            {
                factor = 1;
            }
            else
            {
                factor = AmplificationFactor*ij_d;
            }

            (*this)[i][j] = factor*( (*this)[i][j]);
        }
}


void M3Matrix::GenerateWaveletVisualisation(M3Matrix& Visualisation, int NbDec)
{
    int dyadic = (1<<NbDec);
    double factor = 1;
    double shift = 0;
    Visualisation.Reshape(t_NbRows,t_NbCols);
    Visualisation.SetNbBits(8);
    Visualisation.SetRangeMax(255);
    Visualisation.SetRangeMin(0);

    for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols;j++)
        {
            int ij=i/(t_NbRows/dyadic);
            if(ij<j/(t_NbCols/dyadic)) {
                ij=j/(t_NbCols/dyadic);
            }

            if(ij == 0)
            {
                factor = 1/(double(1<<NbDec));
                shift = 0.;
            }
            else
            {
                factor = 2;
                shift = 128;
            }

            Visualisation[i][j] = factor*( (*this)[i][j]) + shift;
        }
}


//This function sets to zero coefficients for which absolute value is below threshold
//counts number of non-zero coefficients
int M3Matrix::WaveletHardThreshold(double Threshold,int NbDecomp)
{
    int dyadic = (1<<NbDecomp);
    cout << "dyadic hard threshold: " << dyadic << endl;

    int HeightNonDecomposedBlock = t_NbRows/dyadic;
    int WidthNonDecomposedBlock = t_NbCols/dyadic;
    int NonZero = 0;

    for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols;j++)
        {
            if(i >= HeightNonDecomposedBlock || j  >= WidthNonDecomposedBlock)
            {
                if(fabs((*this)[i][j])<Threshold)
                    (*this)[i][j] = 0.;
                else
                    NonZero++;
            }
            else //we are in the "coarse non-detail" block
            {
                NonZero++;
            }
        }

    cout << "Number of NonZero coefficients after hard thresholding: " << NonZero << endl;
    return NonZero;
}


int M3Matrix::WaveletSoftThreshold(double Threshold,int NbDecomp)
{
    int dyadic = (1<<NbDecomp);

    int HeightNonDecomposedBlock = t_NbRows/dyadic;
    int WidthNonDecomposedBlock = t_NbCols/dyadic;
    int NonZero = 0;

    for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols;j++)
        {
            if(i >= HeightNonDecomposedBlock || j  >= WidthNonDecomposedBlock)
            {
                if(fabs((*this)[i][j])<Threshold)
                    (*this)[i][j] = 0.;
                else if((*this)[i][j] <0)
                {
                    (*this)[i][j] += Threshold;
                    NonZero++;
                }
                else if((*this)[i][j] > 0)
                {
                    (*this)[i][j] -= Threshold;
                    NonZero++;
                }
            }
            else //we are in the "coarse non-detail" block
            {
                NonZero++;
            }
        }

    cout << "Number of NonZero coefficients after softthresholding: " << NonZero << endl;
    return NonZero;
}


//Interesting because of "ideal" theoretical properties
int M3Matrix::WaveletGarroteThreshold(double Threshold,int NbDecomp)
{
    int dyadic = (1<<NbDecomp);
    int HeightNonDecomposedBlock = t_NbRows/dyadic;
    int WidthNonDecomposedBlock = t_NbCols/dyadic;
    int NonZero = 0;

    for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols;j++)
        {
            if(i >= HeightNonDecomposedBlock || j  >= WidthNonDecomposedBlock)
            {
                if(fabs((*this)[i][j]) <= Threshold)
                    (*this)[i][j] = 0.;
                else
                {
                    double factor = ((*this)[i][j]*(*this)[i][j] - Threshold*Threshold)/((*this)[i][j]*(*this)[i][j]);
                    (*this)[i][j] *= factor;
                    NonZero++;
                }
            }
            else //we are in the "coarse non-detail" block
            {
                NonZero++;
            }
        }

    cout << "Number of NonZero coefficients after soft thresholding: " << NonZero << endl;
    return NonZero;
}


//This function performs wavelet thresholding
//counts number of non-zero coefficients
int M3Matrix::WaveletThreshold(double Threshold,int NbDecomp,int ThresholdType)
{
    int NonZero = 0;
    if(ThresholdType == 1)
    {
        NonZero = (*this).WaveletHardThreshold(Threshold,NbDecomp);
    }
    else if(ThresholdType == 2)
    {
        NonZero = (*this).WaveletSoftThreshold(Threshold,NbDecomp);
    }
    else if(ThresholdType == 3)
    {
        NonZero = (*this).WaveletGarroteThreshold(Threshold,NbDecomp);
    }
    else
    {
        cout << "Warning: Thresholding type not recognized. Nothing is done " << endl;
    }

    return NonZero;
}


//This function sets to zero coefficients for which absolute value is below threshold
//counts number of non-zero coefficients
int M3Matrix::HardThreshold(double Threshold)
{
    int NonZero = 0;

    for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols;j++)
        {
            if(fabs((*this)[i][j])<Threshold)
                (*this)[i][j] = 0.;
            else
                NonZero++;
        }

    return NonZero;
}


int M3Matrix::SoftThreshold(double Threshold)
{
    int NonZero = 0;

    for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols;j++)
        {
            if(fabs((*this)[i][j])<Threshold)
                (*this)[i][j] = 0.;
            else if((*this)[i][j] <0)
            {
                (*this)[i][j] += Threshold;
                NonZero++;
            }
            else if((*this)[i][j] > 0)
            {
                (*this)[i][j] -= Threshold;
                NonZero++;
            }
        }

    cout << "Number of NonZero coefficients after soft thresholding: " << NonZero << endl;
    return NonZero;
}


//Interesting because of "ideal" theoretical properties
int M3Matrix::GarroteThreshold(double Threshold)
{
    int NonZero = 0;

    for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols;j++)
        {
            if(fabs((*this)[i][j]) <= Threshold)
                (*this)[i][j] = 0.;
            else
            {
                double factor = ((*this)[i][j]*(*this)[i][j] - Threshold*Threshold)/((*this)[i][j]*(*this)[i][j]);
                (*this)[i][j] *= factor;
                NonZero++;
            }
        }

    cout << "Number of NonZero coefficients after garrote thresholding: " << NonZero << endl;
    return NonZero;
}


//This function performs wavelet thresholding
//counts number of non-zero coefficients
int M3Matrix::Threshold(double Threshold,int ThresholdType)
{
    int NonZero = 0;
    if(ThresholdType == 1)
    {
        NonZero = (*this).HardThreshold(Threshold);
    }
    else if(ThresholdType == 2)
    {
        NonZero = (*this).SoftThreshold(Threshold);
    }
    else if(ThresholdType == 3)
    {
        NonZero = (*this).GarroteThreshold(Threshold);
    }
    else
    {
        cout << "Warning: Thresholding type not recognized. Nothing is done " << endl;
    }

    return NonZero;
}


void M3Matrix::Eig2( double& lambda1, double& lambda2, M3Matrix& e1, M3Matrix& e2)
{
    if((*this).GetNbRows() != 2  || (*this).GetNbCols() != 2)
        return;

    e1.Reshape(2,1);
    e2.Reshape(2,1);

    double tr,det,a,b,c,d;
    a = (*this)[0][0];
    b = (*this)[0][1];
    c = (*this)[1][0];
    d = (*this)[1][1];

    tr= a+d;
    det = a*d -c*b;
    lambda1 = 0.5*tr + sqrt(tr*tr*0.25-det);
    lambda2 = 0.5*tr - sqrt(tr*tr*0.25-det);

    if(c!=0)
    {
        e1[0][0] = lambda1-d;
        e1[0][1] = c;

        e2[0][0] = lambda2-d;
        e2[0][1] = c;
    }
    else if(b!=0)
    {
        e1[0][0] = b;
        e1[0][1] = lambda1-a;

        e2[0][0] = b;
        e2[0][1] = lambda2-a;
    }
    else
    {
        e1[0][0] = 1;
        e1[0][1] = 0;

        e2[0][0] = 0;
        e2[0][1] = 1;
    }


    if(fabs(lambda1) > fabs(lambda2))
    {
        double tmp;
        M3Matrix e_tmp;
        tmp = lambda1;
        lambda1 = lambda2;
        lambda2 = tmp;

        e_tmp = e1;
        e1 = e2;
        e2 = e_tmp;
    }

    e1 /= e1.LpNorm(2);
    e2 /= e2.LpNorm(2);
}


double M3Matrix::MeshGrid(int NbXSamples,int NbYSamples, double xMin, double xMax, double yMin, double yMax)
{
    (*this).Reshape(NbXSamples*NbYSamples,3);

    double Dx = fabs(xMax-xMin);
    double Dy = fabs(yMax-yMin);

    for(int i=0;i<NbXSamples;i++)
        for(int j=0;j<NbYSamples;j++)
        {
            double x = xMin + Dx*((double)i)/((double)(NbXSamples-1));
            double y = yMin + Dy * ((double)j)/((double)(NbYSamples-1));

            (*this)[i*NbXSamples+j][0] = x; // x coordinate
            (*this)[i*NbXSamples+j][1] = y; // y coordinate
        }
    return Dx/((double)(NbXSamples-1));
}


// Resolution of a simple linear system with a simple Gaussian pivot, solve A*this = b
// Ainv: will contain the inverse
// Warning: this is a quick and dirty (i.e. not optimized) implementation of the Gauss-Jordan elimination
// Use it only for first versions and small matrices !!! (approximately: less than 20x20)
void M3Matrix::SolveLinear(M3Matrix& A, M3Matrix&  b, M3Matrix& Ainv)
{
    M3Matrix Temp, I;
    int NbRows = A.GetNbRows(), NbCols = A.GetNbCols();
    int NbColsb = b.GetNbCols();
    int  i, j, k, l, iMax;
    double Pivot, tmp, tmp2;

    if(NbRows != NbCols)
    {
        cout << "Warning: SolveLinear: attempt to solve a non square linear system" << endl;
        return;
    }

    Temp = A;

    //Ainv initialised to identity
    Ainv.Set2Id(NbRows);

    (*this).Reshape(NbRows,NbColsb);
    for(i=0; i<NbRows;i++)
    {
        for(j=0; j<NbColsb;j++)
        {
            (*this)[i][j] = b[i][j];
        }
    }


    //Triangulation of the matrix
    for( i = 0; i < NbRows ; i++ )
    {
        Pivot = Temp[i][i];  iMax = i;
        for(j = i + 1; j < NbRows;j++)
        {
            if (fabs( Temp[j][i] ) > fabs( Pivot ))
            { Pivot = Temp[j][i];  iMax = j; }
        }
        if(Pivot == 0) // non-invertible
        {
            return;
        }

        if( iMax != i) //Exchange of lines i and iMax
        {
            for( k = i; k <= NbCols - 1; k++)
            {
                tmp = Temp[i][k];  Temp[i][k] = Temp[iMax][k];  Temp[iMax][k] = tmp;
                tmp = Ainv[i][k];  Ainv[i][k] = Ainv[iMax][k];  Ainv[iMax][k] = tmp;
            }
            for (k = 0; k<= (*this).GetNbCols() - 1;k++)
            {
                tmp = (*this)[i][k];  (*this)[i][k] = (*this)[iMax][k];  (*this)[iMax][k] = tmp;
            }
        }

        //Substraction of the rows
        for(l = 0; l < NbRows;l++)
        {
            if(l != i)
            {
                tmp = Temp[l][i] / Temp[i][i];
                for( k = 0; k < NbCols ; k++)
                {
                    Ainv[l][k] = Ainv[l][k] - tmp * Ainv[i][k];
                    tmp2 = Temp[i][k] * tmp;  Temp[l][k] = Temp[l][k] - tmp2;
                }
                for( k = 0; k< (*this).GetNbCols();k++)
                {
                    tmp2 = (*this)[i][k] * tmp;  (*this)[l][k] = (*this)[l][k] - tmp2;
                }
            }
        }

        //Normalisation
        double tmp = Temp[i][i];
        for(k=0;k<NbCols;k++)
        {
            Ainv[i][k] /= tmp;
            Temp[i][k] /= tmp;
        }
        for(k=0;k<NbColsb;k++)
        {
            (*this)[i][k] /= tmp;
        }
    }

    cout << "Temp"  << endl;
    Temp.Dump();
    cout << endl;

    return;
}


//*****************************************************
//** M3Matrix3D class
//*****************************************************
//Default constructor
M3Matrix3D::M3Matrix3D()
{
    t_NbFrames = 0;
    m_ImageStack = 0;
}


M3Matrix3D::M3Matrix3D(int NbFrames)
{
    t_NbFrames = NbFrames;
    m_ImageStack = new M3Matrix[NbFrames];
}


M3Matrix3D::M3Matrix3D(int NbFrames, int NbRows, int NbCols)
{
    t_NbFrames = NbFrames;
    m_ImageStack = new M3Matrix[NbFrames];

    for(int f = 0;f<t_NbFrames;f++)
    {
        m_ImageStack[f].Reshape(NbRows,NbCols);
    }
}


M3Matrix3D::M3Matrix3D(int NbFrames, int NbRows, int NbCols,  const double& Value)
{
    t_NbFrames = NbFrames;
    m_ImageStack = new M3Matrix[NbFrames];

    for(int f = 0;f<t_NbFrames;f++)
    {
        m_ImageStack[f].Reshape(NbRows,NbCols);
        m_ImageStack[f].SetValues(Value);
    }
}



M3Matrix3D::~M3Matrix3D()
{
    delete[] m_ImageStack;
    m_ImageStack = NULL;
}


void M3Matrix3D::Reshape(int NbRows, int NbCols)
{
    for(int f=0;f<t_NbFrames;f++)
    {
        m_ImageStack[f].Reshape(NbRows,NbCols);
    }
}


void M3Matrix3D::Reshape(int NbFrames)
{
    if(t_NbFrames != 0) // if(m_ImageStack != 0) ?
    {
        delete [] m_ImageStack;
        m_ImageStack = 0;
    }

    t_NbFrames = NbFrames;
    m_ImageStack = new M3Matrix[NbFrames];
}


void M3Matrix3D::Reshape(int NbFrames, int NbRows, int NbCols)
{
    if(t_NbFrames != 0) // if(m_ImageStack != 0)
    {
        delete [] m_ImageStack;
        m_ImageStack = 0;
    }

    t_NbFrames = NbFrames;
    if(NbRows*NbCols*t_NbFrames != 0 )
    {
        m_ImageStack = new M3Matrix[t_NbFrames];
        for(int f=0;f<t_NbFrames;f++)
        {
            m_ImageStack[f].Reshape(NbRows,NbCols);
        }
    }
}


int M3Matrix3D::GetNbRows()
{
    if(t_NbFrames > 0)
        return m_ImageStack[0].GetNbRows();
    else
        return -1;
}


int M3Matrix3D::GetNbCols()
{
    if(t_NbFrames >0)
        return m_ImageStack[0].GetNbCols();
    else
        return -1;
}


int M3Matrix3D::GetNbFrames()
{
    return t_NbFrames;
}


double  M3Matrix3D::GetSum()
{
    double Sum_val = 0.;

    for(int f = 0;f<t_NbFrames;f++)
    {
        for(int i = 0;i<(*this)[0].GetNbRows();i++)
        {
            for(int j = 0;j<(*this)[0].GetNbCols();j++)
            {
                Sum_val += (*this)[f][i][j];
            }
        }
    }

    return Sum_val;
}


//Warning: assumes positive entries
void M3Matrix3D::Normalise()
{
    double Sum_val = (*this).GetSum();
    if(Sum_val > 0)
    {
        (*this) /= Sum_val;
    }
}


void M3Matrix3D::Dump(ostream& Out)
{
    for(int f=0;f<t_NbFrames;f++)
    {
        Out << "Matrix " << f << endl;
        m_ImageStack[f].Dump(Out);
    }
}


void M3Matrix3D::SetValues(const double& Value)
{
    for(int f=0;f<t_NbFrames;f++)
    {
        m_ImageStack[f].SetValues(Value);
    }
}


M3Matrix& M3Matrix3D::operator[](int Frame)
{
    return m_ImageStack[Frame];
}


//Addition of a matrix and a scalar: use C = A+t
void M3Matrix3D::operator+=(const double& t)
{
    for(int f=0;f<t_NbFrames;f++)
        (*this)[f] += t;

    return;
}


//Addition of a matrix and a scalar: use C = A+t
void M3Matrix3D::operator/=(const double& t)
{
    for(int f=0;f<t_NbFrames;f++)
        (*this)[f] /= t;

    return;
}


//Addition of a matrix and a scalar: use C = A+t
void M3Matrix3D::operator*=(const double& t)
{
    for(int f=0;f<t_NbFrames;f++)
        (*this)[f] *= t;

    return;
}


//Addition of a matrix and a scalar: use C = A-t
void M3Matrix3D::operator-=(const double& t)
{
    for(int f=0;f<t_NbFrames;f++)
        (*this)[f] -= t;

    return;
}


//Addition of a matrix and a scalar: use C = A/t
M3Matrix3D M3Matrix3D::operator/(const double& t)
{
    int NbRows = (*this)[0].GetNbRows();
    int NbCols = (*this)[0].GetNbCols();

    M3Matrix3D C(t_NbFrames,NbRows,NbCols);

    for(int f=0;f<t_NbFrames;f++)
        for(int i=0;i<NbRows;i++)
            for(int j=0;j<NbCols;j++)
                C[f][i][j] =(*this)[f][i][j]/t;

    return C;
}


void M3Matrix3D::WaveletEnhancement(double AmplificationFactor, int NbDec)
{
    double factor = 1;
    NbDec = NbDec;

    for(int f=1;f<t_NbFrames;f++)
    {
        int f_d = int(log(double(f))/log(2))+1;
        factor =  AmplificationFactor * f_d;

        (*this)[f] *=  factor;
    }
}


void M3Matrix3D::GenerateWaveletVisualisation(M3Matrix& Visualisation, int NbDec)
{
    //int dyadic = (1<<NbDec);
    double factor = 1;
    double shift = 0;
    int NbRows = (*this)[0].GetNbRows();
    int NbCols = (*this)[0].GetNbCols();

    Visualisation.Reshape(t_NbFrames*NbRows,NbCols*3);

    Visualisation.SetNbBits(8);
    Visualisation.SetRangeMax(255);
    Visualisation.SetRangeMin(0);

    M3Matrix Tmp,Tmp1,Tmp2,Dummy;
    Dummy.Reshape(NbRows,NbCols);

    Tmp1.CatMat((*this)[0],Dummy,2);
    Tmp.CatMat(Tmp1,Dummy,2);
    Tmp /= (double(1<<NbDec)); //why ?

    cout << "NbFrames : " << t_NbFrames << endl;

    shift = 128;
    for(int f = 1;f<t_NbFrames;f++)
    {
        factor = 2;
        Tmp1 = (*this)[f];

        Tmp1 *= factor;
        Tmp1 += shift;

        Tmp2.CatMat(Tmp,Tmp1,1);
        Tmp = Tmp2;
    }

    Visualisation = Tmp;
}
//*****************************************************
//** End of  M3Matrix3D class
//*****************************************************

//******************************************************
// Out of class utilities ******************************
//******************************************************


// Computes the Cholesky decomposition of a positive definite symmetric matrix A=LL^t
int CholeskyDecomp(M3Matrix& A, M3Matrix& L)
{
    int i,j,k;
    double sum;
    int n = A.GetNbRows();
    L.Reshape(n,n);

    for (i=0; i<n; i++)
    {
        for (j=0; j<i; j++)
        {
            sum = 0.;
            for (k=0; k<j; k++)
                sum += L[i][k] * L[j][k];

            L[i][j] = (A[i][j] - sum)/L[j][j];
        }

        sum = 0;
        for (k=0; k<i; k++)
            sum += L[i][k] * L[i][k];
        sum = A[i][i] - sum;
        if (sum<TOL) return 1; // not positive-definite
        L[i][i] = sqrt(sum);
    }

    return 0;
}


// Solve the equation Ax=b for a symmetric positive-definite matrix A,
// using the Cholesky decomposition A=LL^T.  The matrix L is passed in "l".
// Elements above the diagonal are ignored.
void CholeskySolve(M3Matrix&  A,  M3Matrix& x, M3Matrix& b)
{
    int n = A.GetNbRows();
    if(x.GetNbRows() != n || x.GetNbCols() != 1 )
        x.Reshape(n,1);

    M3Matrix L;
    CholeskyDecomp(A,L);
    double sum;

    // solve L*y = b for y (where x[] is used to store y)
    for (int i=0; i < n; i++)
    {
        sum = 0;
        for (int j=0; j<i; j++)
            sum += L[i][j] * x[j][0];

        x[i][0] = (b[i][0] - sum)/L[i][i];
    }

    // solve L^T*x = y for x (where x[] is used to store both y and x)
    for (int i=n-1; i>=0; i--)
    {
        sum = 0;
        for (int j=i+1; j<n; j++)
            sum += L[j][i] * x[j][0];

        x[i][0] = (x[i][0] - sum)/L[i][i];
    }
}


//This is for a general matrix
void CholeskySolve_GenericMatrix(M3Matrix&  A,  M3Matrix& x, M3Matrix& b)
{
    M3Matrix AAt, Atb;
    AAt = Transpose(A) * A;
    Atb = Transpose(A) *b;
    CholeskySolve(AAt,x,Atb);
}



//Wavelet filter handling
void MakeWavelet(M3Matrix& LowFilter,M3Matrix&  HighFilter,M3Matrix&  LowInvFilter,M3Matrix&  HighInvFilter,M3Matrix&  LowFilterL1, M3Matrix& HighFilterL1,int WaveletType)
{
    if(WaveletType == 0)
    {
        MakeHaar(LowFilter, HighFilter, LowInvFilter, HighInvFilter, LowFilterL1, HighFilterL1);
    }
    else if(WaveletType == 1)
    {
        MakeDaub4(LowFilter, HighFilter, LowInvFilter, HighInvFilter, LowFilterL1, HighFilterL1);
    }
}


//Wavelet filters
void MakeHaar(M3Matrix& LowFilter,M3Matrix&  HighFilter,M3Matrix&  LowInvFilter,M3Matrix&  HighInvFilter,M3Matrix&  LowFilterL1, M3Matrix& HighFilterL1)
{
    LowFilter.Reshape(2,1);
    HighFilter.Reshape(2,1);
    LowInvFilter.Reshape(2,1);
    HighInvFilter.Reshape(2,1);

    //L2 normalised Haar filter
    LowFilter[0][0] = 1.0/sqrt(2.0);
    LowFilter[1][0] = 1.0/sqrt(2.0);
    HighFilter[0][0] = 1.0/sqrt(2.0);
    HighFilter[1][0] = -1.0/sqrt(2.0);
    //and its inverse
    LowInvFilter[0][0] = 1.0/sqrt(2.0);
    LowInvFilter[1][0] = 1.0/sqrt(2.0);
    HighInvFilter[0][0] = 1.0/sqrt(2.0);
    HighInvFilter[1][0] = -1.0/sqrt(2.0);

    //L1 normalised Haar filter
    LowFilterL1.Reshape(2,1);
    HighFilterL1.Reshape(2,1);

    LowFilterL1[0][0] = 0.5;
    LowFilterL1[1][0] = 0.5;
    HighFilterL1[0][0] = 0.5;
    HighFilterL1[1][0] = -0.5;
}


void MakeDaub4(M3Matrix& LowFilter,M3Matrix&  HighFilter,M3Matrix&  LowInvFilter,M3Matrix&  HighInvFilter,M3Matrix&  LowFilterL1, M3Matrix& HighFilterL1)
{
    LowFilter.Reshape(4,1);
    HighFilter.Reshape(4,1);
    LowInvFilter.Reshape(5,1);
    HighInvFilter.Reshape(5,1);

    LowFilter[3][0] = (1.0+sqrt(3.0))/(4*sqrt(2.0));
    LowFilter[2][0] = (3.0+sqrt(3.0))/(4*sqrt(2.0));
    LowFilter[1][0] = (3.0-sqrt(3.0))/(4*sqrt(2.0));
    LowFilter[0][0]= (1.0-sqrt(3.0))/(4*sqrt(2.0));

    HighFilter[3][0] = (1.0-sqrt(3.0))/(4*sqrt(2.0));
    HighFilter[2][0] = (-3.0+sqrt(3.0))/(4*sqrt(2.0));
    HighFilter[1][0] = (3.0+sqrt(3.0))/(4*sqrt(2.0));
    HighFilter[0][0]= (-1.0-sqrt(3.0))/(4*sqrt(2.0));

    LowInvFilter[0][0] = 0.0;
    LowInvFilter[1][0] = 0.0;
    LowInvFilter[1][0] = (1.0+sqrt(3.0))/(4*sqrt(2.0));
    LowInvFilter[2][0] = (3.0+sqrt(3.0))/(4*sqrt(2.0));
    LowInvFilter[3][0] = (3.0-sqrt(3.0))/(4*sqrt(2.0));
    LowInvFilter[4][0]= (1.0-sqrt(3.0))/(4*sqrt(2.0));

    HighInvFilter[0][0] = 0.0;
    HighInvFilter[1][0] = 0.0;
    HighInvFilter[1][0] = (1.0-sqrt(3.0))/(4*sqrt(2.0));
    HighInvFilter[2][0] = (-3.0+sqrt(3.0))/(4*sqrt(2.0));
    HighInvFilter[3][0] = (3.0+sqrt(3.0))/(4*sqrt(2.0));
    HighInvFilter[4][0]= (-1.0-sqrt(3.0))/(4*sqrt(2.0));

    double tmp_1 = 0;
    for(int i = 0;i<HighFilter.GetNbRows();i++)
        tmp_1 += HighFilter[i][0]*HighFilter[i][0];

    LowFilterL1 = LowFilter/LowFilter.LpNorm(1);
    HighFilterL1 = HighFilter/HighFilter.LpNorm(1);
}
// end of wavelet filters


void RGB2YUVDCT(const M3Matrix& R,const M3Matrix& G,const M3Matrix& B,M3Matrix& Y,M3Matrix& U,M3Matrix& V)
{
    double r,g,b,y,u,v;
    int NbRows = R.GetNbRows();
    int NbCols = R.GetNbCols();
    Y.Reshape(NbRows,NbCols);
    U.Reshape(NbRows,NbCols);
    V.Reshape(NbRows,NbCols);

    Y.SetNbBits(R.GetNbBits()+1);
    Y.SetRangeMin(0.);
    Y.SetRangeMax(R.GetRangeMax());
    U.SetNbBits(R.GetNbBits()+1);
    U.SetRangeMin(0.);
    U.SetRangeMax(R.GetRangeMax());
    V.SetNbBits(R.GetNbBits()+1);
    V.SetRangeMin(0.);
    V.SetRangeMax(R.GetRangeMax());

    double T[3][3]   =
    {
        {   0.5773502588272094726562500000000000000000,
            0.5773502588272094726562500000000000000000,
            0.5773502588272094726562500000000000000000,     },

        {   0.7071067690849304199218750000000000000000,
            0.0000000000000000000000000000000000000000,
            -0.7071067690849304199218750000000000000000, },
        {
            0.4082483053207397460937500000000000000000,
            -0.8164966106414794921875000000000000000000,
            0.4082483053207397460937500000000000000000      }
    };


    for(int i=0;i<NbRows;i++)
        for(int j=0;j<NbCols;j++)
        {
            r = R[i][j];
            g = G[i][j];
            b = B[i][j];

            y = T[0][0] * r + T[0][1] *g + T[0][2] * b;
            u = T[1][0] * r + T[1][1] *g + T[1][2] * b;
            v = T[2][0] * r + T[2][1] *g + T[2][2] * b;

            Y[i][j] = y;
            U[i][j] = u;
            V[i][j] = v;
        }

    U += (1<<(U.GetNbBits() -1));
    V += (1<<(V.GetNbBits()-1));
}



void YUVDCT2RGB(const M3Matrix& Y, const M3Matrix& U, const M3Matrix& V, M3Matrix& R, M3Matrix& G, M3Matrix& B)
{
    double r,g,b,y,u,v;
    int NbRows = V.GetNbRows();
    int NbCols = V.GetNbCols();
    R.Reshape(NbRows,NbCols);
    G.Reshape(NbRows,NbCols);
    B.Reshape(NbRows,NbCols);

    double T[3][3]   =
    {
        {   0.5773502588272094726562500000000000000000,
            0.5773502588272094726562500000000000000000,
            0.5773502588272094726562500000000000000000,     },

        {   0.7071067690849304199218750000000000000000,
            0.0000000000000000000000000000000000000000,
            -0.7071067690849304199218750000000000000000, },
        {
            0.4082483053207397460937500000000000000000,
            -0.8164966106414794921875000000000000000000,
            0.4082483053207397460937500000000000000000      }
    };

    double bitAdjustU = (1<<(U.GetNbBits()-1));
    double bitAdjustV = (1<<(V.GetNbBits()-1));

    for(int i=0;i<NbRows;i++)
        for(int j=0;j<NbCols;j++)
        {
            y = Y[i][j];
            u = U[i][j] - bitAdjustU;
            v = V[i][j] - bitAdjustV;

            r = T[0][0] * y + T[1][0] * u + T[2][0] * v;
            g = T[0][1] * y + T[1][1] *u + T[2][1] * v;
            b = T[0][2] * y + T[1][2] *u + T[2][2] * v;

            R[i][j] = r;
            G[i][j] = g;
            B[i][j] = b;
        }

    R.SetNbBits(Y.GetNbBits()-1);
    R.SetRangeMin(0.);
    R.SetRangeMax(Y.GetRangeMax());
    G.SetNbBits(Y.GetNbBits()-1);
    G.SetRangeMin(0.);
    G.SetRangeMax(Y.GetRangeMax());
    B.SetNbBits(Y.GetNbBits()-1);
    B.SetRangeMin(0.);
    B.SetRangeMax(Y.GetRangeMax());
}


//***********************************************
// * This is the double transformation RGBtoYUV
//***********************************************
void RGBtoYUV( double r, double g, double b, double &y, double &u, double &v)
{
    y = 0.299*r + 0.587 * g + 0.144*b;
    u = (b-y)*0.493;
    v = (r-y)*0.877;
}


void YUVtoRGB( double y, double u, double v, double &r, double &g, double &b)
{
    b = y + u/0.493;
    r = y + v/0.877;
    g= (y -0.299*r -0.144*b)/0.587;
}


void RGBtoYCoCg( double r, double g, double b, double &y, double &co, double &cg)
{
    y= 0.25* r + 0.5 *g + 0.25*b;
    cg= -0.25* r + 0.5 *g - 0.25*b;
    co= 0.5* r - 0.5*b;
}


void YUVtoRGB(M3Matrix& Y, M3Matrix& U, M3Matrix& V, M3Matrix& R, M3Matrix& G, M3Matrix& B)
{
    double r,g,b,y,u,v;
    int NbRows = V.GetNbRows();
    int NbCols = V.GetNbCols();
    R.Reshape(NbRows,NbCols);
    G.Reshape(NbRows,NbCols);
    B.Reshape(NbRows,NbCols);

    for(int i=0;i<NbRows;i++)
        for(int j=0;j<NbCols;j++)
        {
            y = Y[i][j];
            u = U[i][j];
            v = V[i][j];
            YUVtoRGB(y,u,v,r,g,b);
            R[i][j] = r;
            G[i][j] = g;
            B[i][j] = b;
        }
}


void RGBtoYUV(M3Matrix& R,M3Matrix& G,M3Matrix& B,M3Matrix& Y,M3Matrix& U,M3Matrix& V)
{
    double r,g,b,y,u,v;
    int NbRows = R.GetNbRows();
    int NbCols = R.GetNbCols();
    Y.Reshape(NbRows,NbCols);
    U.Reshape(NbRows,NbCols);
    V.Reshape(NbRows,NbCols);

    for(int i=0;i<NbRows;i++)
        for(int j=0;j<NbCols;j++)
        {
            r = R[i][j];
            g = G[i][j];
            b = B[i][j];
            RGBtoYUV(r,g,b,y,u,v);
            Y[i][j] = y;
            U[i][j] = u;
            V[i][j] = v;
        }
}


vector<double> operator+( const vector<double>& v1, const M3Matrix& V2 )
{
    vector<double> v3;
    if(V2.GetNbRows() == int(v1.size()))
    {
        v3.resize(v1.size());
        for(int i = 0; i< int(v3.size()); i++)
            v3[i] = v1[i] + V2[i][0];
    }
    else
    {
        cout << "operator addition of vector<double> and M3Matrix used for incompatible sizes" << endl;
        exit(-1);
    }

    return v3;
}


vector<double> operator-( const vector<double>& v1, const M3Matrix& V2 )
{
    vector<double> v3;
    if(V2.GetNbRows() == int(v1.size()))
    {
        v3.resize(v1.size());
        for(int i = 0; i< int(v3.size()); i++)
            v3[i] = v1[i] - V2[i][0];
    }
    else
    {
        cout << "operator addition of vector<double> and M3Matrix used for incompatible sizes" << endl;
        exit(-1);
    }

    return v3;
}



//*******************************
//* Generates a synthetic image
//*******************************
void generate_synthetic_image(int option)
{
    M3Matrix Image;
    //Parameters
    int NbRows = 1024;
    int NbCols = 1024;
    string filename,filename_no_extension;
    double id,jd;

    if(option==0) //large circle with some capillary structures
    {
        double r = 0.26;
        double cx=0.5,cy=0.5; //circle center, normalized to 1
        double parabolic_factor = 512;
        double val, val2;
        Image.Reshape(NbRows,NbCols);

        //Generation of a circle with capillar structures
        for(int i = 0; i<NbRows; i++)
            for(int j = 0; j<NbCols; j++)
            {
                id = (double)i/(double)NbRows;
                jd = (double)j/(double)NbCols;

                double sqdist2center = (id-cx)*(id-cx)+(jd-cy)*(jd-cy);
                if(sqdist2center<r*r) //inside the circle
                {
                    Image[i][j] = 160 + ((id-0.45)*(id-0.45)+(jd-0.55)*(jd-0.55))*parabolic_factor;

                    //Here add a capillary structure
                    val = sin(5.2*id)*0.24+0.5;
                    val2 = 1-val;
                    if(fabs(jd-val)<0.016)
                    {
                        Image[i][j] = 160;
                    }

                    if(fabs(jd-val2)<0.012)
                    {
                        Image[i][j] = 160;
                    }

                    val2 = 0.6-0.4*id;
                    if(fabs(jd-val2)<0.01)
                    {
                        Image[i][j] = 150;
                    }

                }
                else
                {
                    Image[i][j] = 100;
                }
            }

        filename_no_extension = "SyntheticTexturedImage";
    }
    //Option = 1 corresponds to the TikZ image: source code in exemple_tikz5.tex
    else if(option==1) //small circles on a gray background, with varying line widths and contrast to the background
    {
        Image.Reshape(NbRows,NbCols);
        double local_gray,int_r,ext_r;
        double Diff = 40;

        //Generation of a circle with capillar structures
        for(int i = 0; i<NbRows; i++)
            for(int j = 0; j<NbCols; j++)
            {
                id = (double)i/(double)NbRows;
                jd = (double)j/(double)NbCols;
                id *= 2.;
                jd *= 2.;

                //Shaded Image
                Image[j][i] = 128-Diff+2*Diff*jd;
                //Unshaded Image
                Image[j][i] = 128; //-Diff+2*Diff*jd;

                //Loop on the circles
                double x,y,z;
                int Width = 5;
                for(int i_x=1;i_x<2*Width+1;i_x++)
                {
                    x = (double) i_x/((double)Width+0.5);
                    z = 0.001 + 0.008*x;// \pgfmathsetmacro{\z}{0.025*\x}

                    for(int i_y=1;i_y < 2*Width+1;i_y++)
                    {
                        y = (double) i_y/((double)Width+0.5);

                        local_gray = 128+255*(0.5*y-0.5);
                        int_r = 0.04-z/2.; //inter
                        ext_r = 0.04+z/2.; //external radius
                        if((id-x)*(id-x)+(jd-y)*(jd-y) < ext_r*ext_r && (id-x)*(id-x)+(jd-y)*(jd-y) > int_r*int_r)
                        {
                            Image[j][i] = local_gray;
                        }
                    }
                }
            }
        filename_no_extension ="Circles";
    }
    else if(option==2) //small circles on a gray background, with varying line widths and contrast to the background
    {
        Image.Reshape(NbRows,NbCols);
        double local_gray,int_r,ext_r;
        double radius;

        //Generation of a circle with capillar structures
        for(int i = 0; i<NbRows; i++)
            for(int j = 0; j<NbCols; j++)
            {
                id = (double)i/(double)NbRows;
                jd = (double)j/(double)NbCols;

                id *= 2.;
                jd *= 2.;

                Image[j][i] = 128.;

                //Loop on the circles
                double x,y,z;
                for(int i_x=0;i_x<7;i_x++)
                {
                    x = 0.125 + 0.275*(double) i_x;
                    z = 0.035*x;// \pgfmathsetmacro{\z}{0.025*\x}
                    radius = 0.04 + 0.013*x;

                    for(int i_y=0;i_y < 7;i_y++)
                    {
                        y = 0.125 + 0.275*(double) i_y;

                        local_gray = 70+59*y;
                        int_r = radius-z/2.; //inter
                        ext_r = radius+z/2.; //external radius
                        if((id-x)*(id-x)+(jd-y)*(jd-y) < ext_r*ext_r && (id-x)*(id-x)+(jd-y)*(jd-y) > int_r*int_r)
                        {
                            Image[j][i] = local_gray;
                        }
                    }
                }

                if(jd >1.09 && jd < 1.1)
                {
                    Image[j][i] = 124;
                }
                if(id > 1.07 && id < 1.08)
                {
                    Image[j][i] = 156;
                }
            }

        filename_no_extension ="Circles2";
    }
    else if(option==3) //small circles on a gray background, with varying line widths and contrast to the background
    {
        Image.Reshape(NbRows,NbCols);
        double local_gray,int_r,ext_r;
        double radius;
        double Diff = 40;

        //Generation of a circle with capillar structures
        for(int i = 0; i<NbRows; i++)
            for(int j = 0; j<NbCols; j++)
            {
                id = (double)i/(double)NbRows;
                jd = (double)j/(double)NbCols;

                id *= 2.;
                jd *= 2.;

                //Background value: here shaded
                Image[j][i] = 128-Diff+2*Diff*jd;

                //Loop on the circles
                double x,y,z;
                for(int i_x=0;i_x<7;i_x++)
                {
                    x = 0.125 + 0.275*(double) i_x;
                    z = 0.035*x;// \pgfmathsetmacro{\z}{0.025*\x}
                    radius = 0.04 + 0.013*x;

                    for(int i_y=0;i_y < 7;i_y++)
                    {
                        y = 0.125 + 0.275*(double) i_y;


                        local_gray = 70+59*y;
                        int_r = radius-z/2.; //inter
                        ext_r = radius+z/2.; //external radius
                        if((id-x)*(id-x)+(jd-y)*(jd-y) < ext_r*ext_r && (id-x)*(id-x)+(jd-y)*(jd-y) > int_r*int_r)
                        {
                            Image[j][i] = local_gray;
                        }
                    }
                }

                if(jd >1.09 && jd < 1.1)
                {
                    Image[j][i] = 124;
                }
                if(id > 1.07 && id < 1.08)
                {
                    Image[j][i] = 156;
                }
            }

        filename_no_extension = "Circles2_Shaded";
    }
    else if(option == 4)
    {
        NbRows = 256;
        NbCols = 256;

        Image.Reshape(NbRows,NbCols);

        for(int i = 0;i<NbRows;i++)
            for(int j = 0;j<NbCols;j++)
            {
                if(j<NbCols/2)
                {
                    Image[i][j] = 64;
                }
                else
                {
                    Image[i][j] = 196;
                }
            }
        filename="Jump";
    }
    else if(option == 5)
    {
        NbRows = 512;
        NbCols = 512;
        Image.Reshape(NbRows,NbCols);

        double cx = 0.5;
        double cy = 0.5;
        double r = 0.2;


        for(int i = 0;i<NbRows;i++)
            for(int j = 0;j<NbCols;j++)
            {
                double id = (double)i/(double)NbRows;
                double jd = (double)j/(double)NbCols;

                if( (id-cx)*(id-cx) +(jd-cy)*(jd-cy) < r*r )
                {
                    Image[i][j] = 64;
                }
                else
                {
                    Image[i][j] = 196;
                }
            }

        filename_no_extension = "Circle";
    }
    else if(option == 6)
    {
        NbRows = 512;
        NbCols = 512;
        Image.Reshape(NbRows,NbCols);
        Image.SetValues(64.);

        double Value = 196;
        double xStart = -0.1;
        double Width = 0.03;
        double Slope = 0.5;

        Image.DrawDiagonalBand(xStart,  Width,  Value, Slope);
        filename_no_extension = "Geometric";
    }
    else if(option == 7)
    {
        NbRows = 4001;
        NbCols = 4001;

        Image.Reshape(NbRows,NbCols);
        Image.SetValues(0);

        for(int i = 0;i<NbRows;i++)
            for(int j = 0;j<NbCols;j++)
            {
                double id = (double)i/(double)NbRows;
                double jd = (double)j/(double)NbCols;

                if(id <0.5)
                {
                    if(jd <=0.5)
                    {
                        Image[i][j] = jd*2*255;
                    }
                    else
                    {
                        Image[i][j] = (2-2*jd)*255;
                    }
                }
                else
                {
                    if(jd <=0.5)
                    {
                        Image[i][j] = (1-2*jd)*255;
                    }
                    else
                    {
                        Image[i][j] = (2*jd-1)*255;
                    }
                }
            }

        filename_no_extension = "ContrastedRamp";
    }
    else if(option == 8)
    {
        NbRows = 512;
        NbCols = 512;
        Image.Reshape(NbRows,NbCols);
        Image.SetValues(64.);

        double Value = 196;
        double radius  = 0.1;
        double theta = 0.5;
        double cx = 0.5;
        double cy = 0.5;
        double a = 1;
        double b = 2;

        Image.DrawEllipse(radius,  cx, cy, a, b, theta, Value);
        filename_no_extension = "Ellipse";
    }

    if(option == 9)
    {
        NbRows = 512;
        NbCols = 512;
        Image.Reshape(NbRows,NbCols);
        Image.SetValues(128); // gray background

        Image.DrawRectangle(0.1, 0.1, 0.4, 0.4,  192);
        Image.DrawRectangle(0.6, 0.1, 0.9, 0.4,  64);
        Image.DrawRectangle(0.1, 0.6, 0.4, 0.9,  0);
        Image.DrawRectangle(0.6, 0.6, 0.9, 0.9,  255);

        Image.DrawRectangle(0.2, 0.7, 0.3, 0.8,  255);
        Image.DrawRectangle(0.7, 0.7, 0.8, 0.8,  0);
        filename_no_extension = "Squares2";
    }
    else if(option == 10)
    {
        NbRows = 512;
        NbCols = 512;
        Image.Reshape(NbRows,NbCols);
        Image.SetValues(130); // gray background

        Image.DrawRectangle(0.1, 0.1, 0.4, 0.4,  185);
        Image.DrawRectangle(0.6, 0.1, 0.9, 0.4,  75);
        Image.DrawRectangle(0.1, 0.6, 0.4, 0.9,  20);
        Image.DrawRectangle(0.6, 0.6, 0.9, 0.9,  240);

        Image.DrawRectangle(0.2, 0.2, 0.3, 0.3,  75);
        Image.DrawRectangle(0.2, 0.7, 0.3, 0.8,  240);
        Image.DrawRectangle(0.7, 0.7, 0.8, 0.8,  20);
        filename_no_extension = "Squares3";
    }
    else if(option == 11)
    {
        NbRows = 512;
        NbCols = 512;
        Image.Reshape(NbRows,NbCols);
        Image.SetValues(130); // gray background

        //Top left
        Image.DrawRectangle(0.1, 0.05, 0.4, 0.35,  185);
        Image.DrawRectangle(0.2, 0.15, 0.3, 0.25,  75);

        //Bottom right
        Image.DrawRectangle(0.6, 0.6, 0.9, 0.9,  20);
        Image.DrawRectangle(0.7, 0.7, 0.8, 0.8,  240);

        //Cross
        Image.DrawLine(0.1,0.6,180,0.7854,0.45,0.0023);
        Image.DrawLine(0.4,0.6,180,3.14159*0.75,0.45,0.0023);

        //Vertical lines
        Image.DrawLine(0.4,0.55,180,3.14159,0.3,0.002);
        Image.DrawLine(0.4,0.5,180,3.14159,0.3,0.004);
        Image.DrawLine(0.4,0.45,180,3.14159,0.3,0.006);

        //Bottom left: some elliptic crowns
        Image.DrawEllipticCrown(0.1,0.11,0.75,0.25,1,1.5,0.5,80);
        Image.DrawEllipticCrown(0.1,0.11,0.75,0.25,1,1.4,1.8,110);

        //Image.DrawRectangle(0.7, 0.7, 0.8, 0.8,  20);
        filename_no_extension = "Geometries";
    }

    //Original Image
    filename =filename_no_extension+".pgm";
    Image.Digitalize();
    Image.SavePGM(filename);

    M3Matrix TmpImage,BlurredImage,Filter;

    Filter.MakeGaussianFilter(5);
    Image.LinearFilter(TmpImage,Filter);
    TmpImage.LinearFilter(BlurredImage,Filter);
    BlurredImage.Digitalize();
    filename = filename_no_extension +"_Blurred.pgm";
    BlurredImage.SavePGM(filename);
}


void generate_color_image()
{
    M3Matrix R1,G1,B1;
    int NbRows = 200, NbCols = 300;
    R1.Reshape(NbRows,NbCols);
    G1.Reshape(NbRows,NbCols);
    B1.Reshape(NbRows,NbCols);

    for(int i =0;i<NbRows;i++)
    {
        for(int j = 0;j<NbCols/3;j++)
        {
            R1[i][j] = 0;
            G1[i][j] = 255;
            B1[i][j] = 0;
        }
        for(int j = 10;j<2*(NbCols/3);j++)
        {
            R1[i][j] = 255;
            G1[i][j] = 255;
            B1[i][j] = 255;
        }
        for(int j = 20;j<NbCols;j++)
        {
            R1[i][j] = 255;
            G1[i][j] = 0;
            B1[i][j] = 0;
        }
    }

    M3Matrix SynthColor;
    SynthColor.CombineRGBChannels(R1,G1,B1);

    string filenameC = "Italien.ppm";
    SynthColor.SavePPM(filenameC);
}


void system_engineering_test()
{
    //Check that own files work
    M3Matrix A(2,2), B(2,2), C;
    A.SetValues(3.0),  B.SetValues(4.0);
    C = A + 3.0 * B;
    cout << "C.Dump()" << endl;
    C.Dump();
}


void test_RGB2XYZ_matrices()
{
    M3Matrix A(3,3),B(3,3),C;
    A = {0.638088, 0.214704, 0.097744,
         0.291954, 0.823841, -0.115795,
         0.002798, -0.067034, 1.153294};

    B= {1.789066, -0.482534, -0.200076,
        -0.639849, 1.396400, 0.194432,
        -0.041532, 0.082335, 0.878868};

    C = A*B;
    cout << "C : " << endl;
    C.Dump(cout,5,0.001);
}


void test_mat_vec_mult()
{
    vector<double> x, y, z;
    x = {1,2,3.3};
    y = {4,2,3.1};
    z=  {0,1,2};
    z = x+2*y-z*3;
    dump_vec(z,"z : ");
    M3Matrix D(3,1);
    D.SetValues(3.14);
    z = D + x;
    dump_vec(z,"z : ");
    z = x - D;
    dump_vec(z,"z : ");

    M3Matrix Z;
    Z = {1,2,3.2,4,5,7};
    Z.Dump(); //Z=  {}
    Z.Reshape(2,3);
    Z.Dump(); //Z=  {}
    Z.Reshape(1,6);
    Z.Dump(); //Z=  {}
    Z.Reshape(3,2);
    Z.Dump(); //Z=  {}

    M3Matrix A(2,3);
    A.SetValues({1,2,3,4,5,6});
    vector<double> b({3,5,7});

    vector<double> c;
    c = A*b;
    dump_vec(c,"c");
}

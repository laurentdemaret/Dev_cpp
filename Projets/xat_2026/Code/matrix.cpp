#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>

#include "bitstr.h"
#include "matrix.h"

#ifdef __cplusplus
extern "C" {
#endif

// Comparison function
int fncompare(const void * el1,const void * el2)
{
  double **d1 = (double **)el1;
  double **d2 = (double **)el2;

  //cout << **d1 << endl;
  //cout << **d2 << endl;

  if(**d1<**d2) return -1;
  if(**d1>**d2) return 1;
  return 0;
}


int fncompare2(const void * el1,const void * el2)
{
  double **d1 = (double **)el1;
  double **d2 = (double **)el2;

  if((*d1)[0]<(*d2)[0]) return -1;
  else
  {
    if((*d1)[0]>(*d2)[0])
      return 1;
    else
    {
      if((*d1)[1]<(*d2)[1]) return -1;
      if((*d1)[1]>(*d2)[1]) return 1;
      return 0;
    }
  }
}
	

// Comparison functions with a bloc criterion
int fncomparebloc(const void * el1, const void * el2)
{
  double **d1 = (double **)el1;
  double **d2 = (double **)el2;

  int x1 = COORD_ROUND((*d1)[0]);
  int x2 = COORD_ROUND((*d2)[0]);
  int y1 = COORD_ROUND((*d1)[1]);
  int y2 = COORD_ROUND((*d2)[1]);

  int BS = 8; // Bloc Size

  if(x1/BS < x2/BS) return -1;
  else
  {
    if(x1/BS>x2/BS)
      return 1;
    else
    {
      if(y1/BS<y2/BS) return -1;
      else
      {
        if(y1/BS>y2/BS) return 1;
        else
        {
          if(x1%BS < x2%BS ) return -1;
          else
          {
            if(x1%BS > x2%BS ) return 1;
            else
            {
              if(y1%BS < y2%BS) return -1;
              if(y1%BS > y2%BS) return -1;
              return 0;
            }

          }
        }
      }
    }
  }
}

	
int comp_int(const void *a, const void *b)
// comparaison pour qsort
{
  return *(int*) a - *(int*) b;
}


#ifdef __cplusplus
}
#endif


using namespace std;

// Default constructor
M3Matrix::M3Matrix()
{
  t_NbRows = 0;
  t_NbCols = 0;
  Data = NULL;
}


// Constructor with dimensions
M3Matrix::M3Matrix(int NbRows, int NbCols)
{
  t_NbRows = NbRows;
  t_NbCols = NbCols;
  Data = NULL;
  Alloc();
}

// Constructor with dimensions and initial value
M3Matrix::M3Matrix(int NbRows, int NbCols, const double& r_Value)
{
  t_NbRows = NbRows;
  t_NbCols = NbCols;

  Data = NULL;
  Alloc();
  SetValues(r_Value);
}

// Copy constructor : r_Source copied into this
M3Matrix::M3Matrix(const M3Matrix& r_Source)
{
  t_NbRows = r_Source.t_NbRows;
  t_NbCols = r_Source.t_NbCols;

  Data = NULL;
  Alloc();
  if(Data)
  {
    for(int i=0;i<t_NbRows;i++)
      for(int j=0;j<t_NbRows;j++)
        Data[i][j] =  r_Source.Data[i][j];
  }
}


// Destructor
M3Matrix::~M3Matrix()
{
    Free();
    Data = NULL; //No clean memory
}

// Accessors
int M3Matrix::GetNbRows()
{
  return t_NbRows;
}

int M3Matrix::GetNbCols()
{
  return t_NbCols;
}

double* M3Matrix::operator[](int Row)
{
  return Data[Row];
}


// Dumps the matrix coefficients on the output Out
void M3Matrix::Dump(ostream& Out)
{
  for (int i=0;i<GetNbRows();i++)
  {
    for (int j=0;j<GetNbCols();j++)
    {
      Out << (*this)[i][j] << " ";
    }
    Out << endl;
  }
}


// Saves the M3Matrix values into a PGM file
// Warning the values are
void M3Matrix::SavePGM(ostream& Out)
{
  Out << "P2" << endl;
  Out << "#CREATOR : MyProgram" << endl;

  Out << t_NbRows << " " << t_NbCols << endl;
  Out << 255 << endl;

  Dump(Out);
}


// Operators
M3Matrix& M3Matrix::operator=(const M3Matrix& r_Src)
{
  // Check that the two objects are initially different
  if(this != &r_Src)
  {
    Reshape(r_Src.t_NbRows,r_Src.t_NbCols);
    if ((t_NbRows)||( t_NbCols))
    {
      for(int i=0;i<t_NbRows;i++)
        for(int j=0;j<t_NbCols;j++)
          Data[i][j] = r_Src.Data[i][j];
    }
  }
  return *this;
}


// Addition of two matrices
M3Matrix* M3Matrix::operator+(M3Matrix& r_MRight)
{
  M3Matrix* somme = new M3Matrix;
  (*somme).Reshape(t_NbRows,t_NbCols);
  if(t_NbRows == r_MRight.GetNbRows() && t_NbCols == r_MRight.GetNbCols())
  {
    for(int i=0;i<t_NbRows;i++)
      for(int j=0;j<t_NbCols;j++)
        (*somme).Data[i][j] = Data[i][j] + r_MRight.Data[i][j];
  }

  return somme;
}

// Multiplication of two matrices (for non sparse matrices)
M3Matrix* M3Matrix::operator*(M3Matrix& r_MRight)
{
  M3Matrix* prod = new M3Matrix;
  double tmp = 0.;
  (*prod).Reshape(t_NbRows,r_MRight.GetNbCols());
  if(t_NbCols == r_MRight.GetNbRows())
  {
    for(int i=0;i<t_NbRows;i++)
      for(int j=0;j<t_NbCols;j++)
      {
        tmp = 0.;
        for(int k=0;k<t_NbCols;k++)
          tmp += Data[i][k] * r_MRight.Data[k][j];

        (*prod).Data[i][j]  = tmp;
      }
  }

  return prod;
}


// Transformers

// Initializes all the M3Matrix coefficients with r_Value
void M3Matrix::SetValues(const double &r_Value)
{
  for(int i=0;i<t_NbRows;i++)
    for(int j=0;j<t_NbCols;j++)
      Data[i][j] = r_Value;
}

// Deallocation
void M3Matrix::Free()
{
  if (Data != 0)
  {
    delete [] Data[0];
    delete [] Data;
    Data = 0;
  }
}

// Destroys all the elements of the array
// and reshape the dimemsions of the array
// Values are initialized to 0.
void M3Matrix::Reshape(int NbRows, int NbCols)
{
  // Memory deallocation
  if((t_NbRows * t_NbCols) !=0)
    Free();

  // Memory allocation
  t_NbRows = NbRows;
  t_NbCols = NbCols;
  Alloc();

  SetValues(0.);
}


// Memory allocation
void M3Matrix::Alloc()
{
  if(!((t_NbRows==0) || (t_NbCols==0)))
  {
    double **pt_RowPtr;
    double *p_Row;
    int RowIndex = t_NbRows;

    // Memory allocation
    Data = new double* [t_NbRows];
    //ASSERT(Data != 0);
    Data[0] = new double[t_NbRows * t_NbCols];
    //ASSERT(Data[0] != 0);

    pt_RowPtr = Data;
    p_Row = Data[0];

    while (RowIndex > 0)
    {
      *pt_RowPtr = p_Row;
      pt_RowPtr++;
      p_Row += t_NbCols;
      RowIndex--;
    }
  }
}


// Computes the minimal coefficient of the array
double M3Matrix::GetMin()
{
  double min = Data[0][0];

  for(int i=0;i<t_NbRows;i++)
    for(int j=0;j<t_NbCols;j++)
    {
      if(Data[i][j]< min)
        min = Data[i][j];
    }
  return min;
}

// Computes the maximal coefficient of the array
double M3Matrix::GetMax()
{
  double max = Data[0][0];

  for(int i=0;i<t_NbRows;i++)
    for(int j=0;j<t_NbCols;j++)
    {
      if(Data[i][j]> max)
        max = Data[i][j];
    }
  return max;
}


// Quantization of the coefficients of the M3Matrix according to the QStep
// The Symbol M3Matrix has a size of the original M3Matrix and contains the symbols
// used to be coded
void M3Matrix::Quantize(double& QStep, M3Matrix& Symbol)
{
  Symbol.Reshape(t_NbRows,t_NbCols);

  double tmp,qvalue;
  double HalfStep = QStep/2;
  double symbol;

  for(int i=0;i<t_NbRows;i++)
    for(int j=0;j<t_NbCols;j++)
    {
      tmp = Data[i][j] + HalfStep;
      symbol = DBL_ROUND(tmp/QStep);
      if (symbol < 0.) symbol -= 1;

      qvalue  = symbol * QStep;

      Data[i][j] = qvalue;
      Symbol[i][j] = symbol;
    }
}

// Sets the line i in the M3Matrix Line
void M3Matrix::SetLine(int i, M3Matrix& Line)
{
  Line.Reshape(1,t_NbCols);

  for(int j=0;j<t_NbCols;j++)
    Line[0][j] = Data[i][j];
}

// Sets the col j in the M3Matrix Col
void M3Matrix::SetCol(int j, M3Matrix& Col)
{
  Col.Reshape(t_NbRows,1);

  for(int i=0;i<t_NbRows;i++)
    Col[i][0] = Data[i][j];
}

// Computes the mean square error between two matrices
// (used for measuring the error between images)
double M3Matrix::MSE(M3Matrix& M3Matrix2)
{
  double tmp = 0.;

  for(int i=0;i<t_NbRows;i++)
    for(int j=0;j<t_NbCols;j++)
      tmp += (Data[i][j] - M3Matrix2[i][j])*(Data[i][j] - M3Matrix2[i][j]);

  tmp /= (double)(t_NbRows*t_NbCols);

  return tmp;
}


// Computes the PSNR between two matrices of chars
// (i.e. each coefficient is written on 8 bits)
// (very useful for image comparisons)
double M3Matrix::PSNR(M3Matrix& M3Matrix2)
{
  double PSNR;
  double MSE = (*this).MSE(M3Matrix2);

  PSNR = 10. * log10(256.*256./MSE);

  return PSNR;
}


// Computes the histogram of a value
// WARNING (VERY IMPORTANT):
// this function is only useful for symbol matrices
// returns the entropy of the signal
double M3Matrix::Histogram(M3Matrix& Histo)
{
  int Min = DBL2INT(GetMin());
  int Max = DBL2INT(GetMax());

  //cout << "Min : " << Min << endl;
  //cout << "Max : " << Max << endl;

  double entropy = 0.;

  Histo.Reshape(Max-Min+1,1);

  for(int i=0;i<t_NbRows;i++)
    for(int j=0;j<t_NbCols;j++)
      Histo[DBL2INT(Data[i][j])-Min][0] += 1;

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

#ifndef NO_DEBUG_MESSAGES
  double cost = entropy * NbCoeff;
  cout << "Entropy of the quantized signal :" << entropy << endl;
  cout << "Theoretical cost of the signal (with a perfect arithmetic encoding and a correct probability estimation) : " << cost << endl;
#endif

  return entropy;
}

// Loads the M3Matrix from InputFile
// WARNING the format must be on the first line : NbRows,NbCols
// and then the values
void M3Matrix::LoadFromFile(const char* InputFile)
{
  ifstream infile(InputFile);
  int NbRows;
  int NbCols;

  infile >> NbRows;
  infile >> NbCols;

  (*this).Reshape(NbRows,NbCols);


  for(int i=0;i<t_NbRows;i++)
    for(int j=0;j<t_NbCols;j++)
      infile >> Data[i][j];
}


// Loads the M3Matrix from InputFile
// WARNING the format must be coherent
// Same number of elements on each line
void M3Matrix::LoadFromFileWithoutHeader(const char* InputFile)
{
  ifstream infile(InputFile);
  int NbRows = 0;
  ifstream* pFile;
  pFile = new ifstream(InputFile);
  pFile->seekg(0);
  char Buffer;

  // TO BE ABSOLUTELY CHANGED
  int NbCols = 3;
  // ... UNTIL HERE

  while(pFile->get(Buffer))
  {
    if(Buffer == '\n')
      NbRows++;
  }


  (*this).Reshape(NbRows,NbCols);
#ifndef NO_DEBUG_MESSAGES
  cout << "NbRows : " << NbRows << endl;
  cout << "NbCols : " << NbCols << endl;
#endif

  for(int i=0;i<t_NbRows;i++)
    for(int j=0;j<t_NbCols;j++)
    {
      infile >> Data[i][j];
    }
#ifndef NO_DEBUG_MESSAGES
  cout << "fin du chargement" << endl;
#endif
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
// WARNING :
void M3Matrix::SortByBlocs(int c1, int c2)
{
  // The column i is set as first col
  if(c1 != 0)
    ExchangeCols(c1,0);
  if(c2 != 1)
    ExchangeCols(c2,1);

  qsort(Data,t_NbRows,sizeof(double *),fncomparebloc);

  // Column newly exchanged to provide a good
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

  if(Data[i1][0] > Data[i2][0])
    res = -1;
  else
  {
    if(Data[i1][0] < Data[i2][0])
      res = 1;
    else
      res = 0;
  }
  return res;
}


// Exchange columns i and j of the M3Matrix
void M3Matrix::ExchangeCols(int j1, int j2)
{
  double tmp;

  for(int i=0;i<t_NbRows;i++)
  {
    tmp = Data[i][j1];
    Data[i][j1] = Data[i][j2];
    Data[i][j2] = tmp;
  }
}


// Exchange rows i and j of the M3Matrix
void M3Matrix::ExchangeRows(int i1, int i2)
{
  double tmp;

  for(int j=0;j<t_NbCols;j++)
  {
    tmp = Data[i1][j];
    Data[i1][j] = Data[i2][j];
    Data[i2][j] = tmp;
  }
}


// Fills up the matrix with the luminance values
// read in a PGM file
void M3Matrix::LoadYFromPGM(const char* infile)
{
  ifstream* pFile;
  pFile = new ifstream(infile);
	if (!*pFile) 
	{
		cout << "Can't open output file (does the file exist ?). Abort ... " << endl;
		exit(1);
	}  
	
  pFile->seekg(0);
  char Buffer;
	
  //int row = 0;
	
	//Added : 07 april 2010
	pFile->get(Buffer); //First character
	while(Buffer == '#')
	{
		while(Buffer != '\n')
		{
			pFile->get(Buffer);
		}
		pFile->get(Buffer);
	}	
	
	//Now Buffer should be the first non commented character
	
	if(Buffer != 'P') //not correct pgm file
	{
		cout << "Error: this is not a  pgm file, abort ..." << endl;
		exit(1);
	}
	
	pFile->get(Buffer);
	if(Buffer !='2')
	{
	  cout 	<< "Error: this is not an ASCII pgm file (probably a binary pgm ?), abort ... " << endl;
		exit(1);
	}	
	
	while(Buffer != '\n')
	{
		pFile->get(Buffer);
	}
	
	Buffer = pFile->peek();
  
	while(Buffer == '#')
	{
		while(Buffer != '\n')
		{
			pFile->get(Buffer);
		}
		Buffer = pFile->peek();
	}	
	
	//while(Buffer )
	//pFile->get(Buffer);
	//end of added (07 april 2010)
	
  //while( row < 2)
  //{
  //  pFile->get(Buffer);
  //  if(Buffer == '\n')
  //    row++;
	// }
	
  int NbRows,NbCols,dum;
	
  *pFile >> NbRows;
  *pFile >> NbCols;
  *pFile >> dum;
	
	if(dum != 255)
	{
	  cout << "warning ! dynamic is not 256 !!!" << endl;
	}
	
  (*this).Reshape(NbCols,NbRows);
	for(int j=0;j<t_NbRows;j++)
    for(int i=0;i<t_NbCols;i++)
      *pFile >> Data[j][i];

  delete pFile;
}


void M3Matrix::EncodeOctTree(WBitStream& OutStream,int NbRows,int NbCols,int quantization)
{
  // WARNING : in this simple first implementation, we use
  // sizes of image : 2^N*2^N where N = 7,8,9 and quantization is
  // 2,4,8 or 16
  unsigned int quantCode; // quantization code - {1,2,3,4}

  // convert quantization
  if(quantization == 2)
    quantCode = 1;
  if(quantization == 4)
    quantCode = 2;
  if(quantization == 8)
    quantCode = 3;
  if(quantization == 16)
    quantCode = 4;

  // saves the number of the rows in 10 bits
  OutStream.SaveValue(10,NbRows);
  // saves the number of the columns in 10 bits
  OutStream.SaveValue(10,NbCols);
  //saves quantization step in 3 bits
  OutStream.SaveValue(3,quantCode);
  //saves the number of vertices (maximal value : NbRows*NbCols)
  OutStream.SaveValue(M3Matrix::getNumberOfBitsNeeded(NbRows*NbCols),t_NbRows);

  int numberOfCells = M3Matrix::getNumberOfBitsNeeded(NbRows)-1;

  //Foo *a = new Foo[n];
  //delete [] a;
  //a = NULL; // clear to prevent using invalid memory reference (dangling pointer)

  //M3Matrix cells[numberOfCells];
  M3Matrix *cells = new M3Matrix[numberOfCells];
  int local_quant = (int)(256/quantization);
  int part_cells = 1;
  for (int i=0;i<numberOfCells;i++)
  {
    if (local_quant >= 2)
      part_cells = part_cells * 8;
    else
      part_cells = part_cells * 4;

    local_quant = (int)(local_quant/2);
    cells[i].Reshape(part_cells,1);
#ifndef NO_DEBUG_MESSAGES
    //cout << part_cells << endl;
#endif
  }

  //Quantization of the values
  for(int i=0;i<(int)t_NbRows;i++)
  {
    (*this)[i][2] /= quantization;
    (*this)[i][2] = DBL_ROUND((*this)[i][2]);
  }

  int x,y,z,x1=0,y1=0,z1=0;
  unsigned int xyz;

  int localquantbits = M3Matrix::getNumberOfBitsNeeded((int)(256/quantization))-1;
//  int luminanceBits = M3Matrix::getNumberOfBitsNeeded(256)-1;
#ifndef NO_DEBUG_MESSAGES
//  cout << "localquantbits" << localquantbits << endl;
#endif
  int localcounter = numberOfCells - localquantbits;
  for(int i=0;i<(int)t_NbRows;i++)
  {
    x = COORD_ROUND((*this)[i][0]);
    y = COORD_ROUND((*this)[i][1]);
    z = DBL2INT((*this)[i][2]);

    x1 = x >> (numberOfCells-1);
    y1 = y >> (numberOfCells-1);
//    z1 = z >> (luminanceBits-1);
    z1 = z >> (numberOfCells-1-localcounter);
    xyz = (unsigned)((x1 << 2) + (y1 << 1) + z1);
    cells[0][xyz][0] += 1;
    for(int j=1;j<numberOfCells;j++)
    {
      x1 = x >> (numberOfCells-j-1);
      y1 = y >> (numberOfCells-j-1);
      if (j >= localquantbits)
      {
        xyz = (unsigned)((xyz << 2) + ((x1&1) << 1) + (y1&1));
      }
      else
      {
        z1 = z >> (numberOfCells-j-1-localcounter);
//        z1 = z >> (luminanceBits-j-1);
        xyz = (unsigned)((xyz << 3) + ((x1&1) << 2) + ((y1&1) << 1) + (z1&1));
      }
      cells[j][xyz][0] += 1;
    }
  }

  //  Transmission cost evaluation

  int HeaderCost = OutStream.GetCounter(), LastCost = HeaderCost, CurrCost;
  int CurrentMax = 0;
  M3Matrix Vector(8,1);

  for(int k=0;k<8;k++)
  {
    Vector[k][0] = cells[0][k][0];
  }
  TransmitCube(Vector,t_NbRows,OutStream);

  CurrCost = OutStream.GetCounter();
#ifndef NO_DEBUG_MESSAGES
  cout << "1 layer cost: " << CurrCost-LastCost << endl;
#endif
  LastCost = CurrCost;

  M3Matrix Vector4(4,1);
  for(int j=1;j<numberOfCells-1;j++)
  {
    for(int i=0;i<cells[j-1].GetNbRows();i++)
    {
      if(cells[j-1][i][0] > 0.)
      {
        CurrentMax = DBL2INT(cells[j-1][i][0]);
        if (j >= localquantbits)
        {
          for(int k=0;k<4;k++)
          {
              Vector4[k][0] = cells[j][4*i+k][0];
          }
          TransmitCube4(Vector4,CurrentMax,OutStream);
        }
        else
        {
          for(int k=0;k<8;k++)
          {
              Vector[k][0] = cells[j][8*i+k][0];
          }
          TransmitCube(Vector,CurrentMax,OutStream);
        }
      }
    }

    CurrCost = OutStream.GetCounter();
#ifndef NO_DEBUG_MESSAGES
    cout << j+1 << " layer cost: " << CurrCost-LastCost << endl;
#endif
    LastCost = CurrCost;
  }
  Vector4.Free();

  for(int i=0;i<cells[numberOfCells-2].GetNbRows();i++)
  {
    if(cells[numberOfCells-2][i][0] > 0.)
    {
      if(cells[numberOfCells-2][i][0] == 1.)
      {
        // Codes only the position of the vertex among the 4 positions
        for(int k=0;k<4;k++)
        {
          if(cells[numberOfCells-1][4*i+k][0] == 1.)
          {
            OutStream.SaveValue(2,k);
          }
        }
      }
      if(cells[numberOfCells-2][i][0] == 2.)
      {
        M3Matrix Vector;
        Vector.Reshape(4,1);
        for(int k=0;k<4;k++)
          Vector[k][0] = cells[numberOfCells-1][4*i+k][0];

        int Code = 0;

        if((Vector[0][0]  == 1.) && (Vector[1][0] == 1.))
          Code = 0;
        if((Vector[2][0]  == 1.) && (Vector[3][0] == 1.))
          Code = 1;
        if((Vector[0][0]  == 1.) && (Vector[2][0] == 1.))
          Code = 2;
        if((Vector[1][0]  == 1.) && (Vector[3][0] == 1.))
          Code = 3;
        if((Vector[0][0]  == 1.) && (Vector[3][0] == 1.))
          Code = 4;
        if((Vector[1][0]  == 1.) && (Vector[2][0] == 1.))
          Code = 5;

        if(Code == 4 || Code == 5)
        {
          OutStream.SaveValue(1,0);
          if(Code == 4)
             OutStream.SaveValue(1,0);
          if(Code == 5)
             OutStream.SaveValue(1,1);
        }
        else
        {
          OutStream.SaveValue(1,1);
          OutStream.SaveValue(2,Code);
        }

      }

      if(cells[numberOfCells-2][i][0] == 3.)
      {
        for(int k=0;k<4;k++)
        {
          if(cells[numberOfCells-1][4*i+k][0] != 1.)
          {
            OutStream.SaveValue(2,k);
          }
        }
      }
    }
  }


  CurrCost = OutStream.GetCounter();
#ifndef NO_DEBUG_MESSAGES
  cout << numberOfCells << " layer cost: " << CurrCost-LastCost << endl;
  cout << "Total cost: " << HeaderCost << "(header) + "<< CurrCost-HeaderCost << "(layers)" <<  endl;
#endif
  for (int i=0;i<numberOfCells;i++)
    cells[i].Free();
    
   delete [] cells;
   cells = NULL; 
}


// Transmission of the number of elements in each cube
void M3Matrix::TransmitCube(M3Matrix& Vector, unsigned int Max, WBitStream& OutStream)
{
  unsigned int CurrentMax = Max;
  int k;

  M3Matrix sum2(4,1);
  M3Matrix sum4(2,1);


  // Transmission of the values
  for(k=0;k<4;k++)
    sum2[k][0] = Vector[2*k][0] + Vector[2*k+1][0];
  for(k=0;k<2;k++)
    sum4[k][0] = sum2[2*k][0] + sum2[2*k+1][0];

  if(CurrentMax>0)
  {
    OutStream.SaveValue(M3Matrix::getNumberOfBitsNeeded(CurrentMax),DBL2INT(sum4[0][0]));
  }

  for(k=0;k<2;k++)
  {
    CurrentMax = DBL2INT(sum4[k][0]);
    if(CurrentMax>0)
    {
      OutStream.SaveValue(M3Matrix::getNumberOfBitsNeeded(CurrentMax),DBL2INT(sum2[2*k][0]));
    }
  }

  for(k=0;k<4;k++)
  {
    CurrentMax = DBL2INT(sum2[k][0]);
    if(CurrentMax>0)
    {
      OutStream.SaveValue(M3Matrix::getNumberOfBitsNeeded(CurrentMax),DBL2INT(Vector[2*k][0]));
    }
  }
}


// Transmission of the number of elements in each cube
void M3Matrix::TransmitCube4(M3Matrix& Vector, unsigned int Max, WBitStream& OutStream)
{
  unsigned int CurrentMax = Max;

  M3Matrix sum4(2,1);

  // Transmission of the values
  for(int k=0;k<2;k++)
    sum4[k][0] = Vector[2*k][0] + Vector[2*k+1][0];
  OutStream.SaveValue(M3Matrix::getNumberOfBitsNeeded(CurrentMax),DBL2INT(sum4[0][0]));

  // for each of these two parts of the cube, examine the subcubes for number of points
  // and save them in binary
  for(int k=0;k<2;k++)
  {
    CurrentMax = DBL2INT(sum4[k][0]);
    if(CurrentMax>0)
      OutStream.SaveValue(M3Matrix::getNumberOfBitsNeeded(CurrentMax),DBL2INT(Vector[2*k][0]));
  }
}


// Decoding of the stream
void M3Matrix::DecodeOctTree(int* rows,int* cols, RBitStream& InStream, int* quant)
{
  // number of points left after the adaptive thinning
  unsigned int NbRows;
  unsigned int quantCode;
  unsigned int numberRows;
  unsigned int numberCols;

  // reads the number of the rows and columns
  InStream.LoadValue(10,numberRows);
  InStream.LoadValue(10,numberCols);
  // reads the quantization step
  InStream.LoadValue(3,quantCode);
  //reads the number of vertices kept for the representation
  InStream.LoadValue(M3Matrix::getNumberOfBitsNeeded(numberRows*numberCols),NbRows);

  *rows = numberRows;
  *cols = numberCols;

  unsigned int CurrentMax = NbRows;
  int quantization = DBL2INT(pow(2.,(double)quantCode));
#ifndef NO_DEBUG_MESSAGES
  //  cout << "quantization: " << quantization << endl;
#endif
  *quant = quantization;

  unsigned int numberOfCells = M3Matrix::getNumberOfBitsNeeded(numberRows)-1;
  M3Matrix *cells = new M3Matrix[numberOfCells];

  int local_quant = (int)(256/quantization);
  int part_cells = 1;
  for (int i=0;i<(int)numberOfCells;i++)
  {
    if (local_quant >= 2)
      part_cells = part_cells * 8;
    else
      part_cells = part_cells * 4;

    local_quant = (int)(local_quant/2);
    cells[i].Reshape(part_cells,1);
#ifndef NO_DEBUG_MESSAGES
    //cout << part_cells << endl;
#endif
  }

  int localquantbits = M3Matrix::getNumberOfBitsNeeded((int)(256/quantization))-1;

  M3Matrix Vector(8,1);

  LoadCube(Vector,CurrentMax,InStream);
  for(int k=0;k<8;k++)
  {
    cells[0][k][0] = Vector[k][0];
  }

  M3Matrix Vector1(4,1);
  for(int j=1;j<(int)numberOfCells-1;j++)
  {
    for(int i=0;i<cells[j-1].GetNbRows();i++)
    {
        if (j>=localquantbits)
        {
          CurrentMax = DBL2INT(cells[j-1][i][0]);
          if(CurrentMax>0)
          {
            LoadCube4(Vector1,CurrentMax,InStream);
          }
          else
            Vector1.Reshape(4,1);

          for(int k = 0;k<4;k++)
          {
             cells[j][4*i+k][0] = Vector1[k][0];
          }
        }
        else
        {
          CurrentMax = DBL2INT(cells[j-1][i][0]);
          if(CurrentMax>0)
            LoadCube(Vector,CurrentMax,InStream);
          else
            Vector.Reshape(8,1);

          for(int k = 0;k<8;k++)
          {
             cells[j][8*i+k][0] = Vector[k][0];
          }
        }
    }
  }
  Vector1.Free();

  unsigned int tmp;
  for(int i=0;i<cells[numberOfCells-2].GetNbRows();i++)
  {
    if(cells[numberOfCells-2][i][0]>0.)
    {
      if(cells[numberOfCells-2][i][0] == 1.)
      {
        InStream.LoadValue(2,tmp);
        cells[numberOfCells-1][4*i + tmp][0] = 1.;
      }
      else if(cells[numberOfCells-2][i][0] == 2.)
      {
        InStream.LoadValue(1,tmp);
        if(tmp == 0)
        {
          InStream.LoadValue(1,tmp);
          if(tmp == 0)
          {
            cells[numberOfCells-1][4*i][0] = 1.;
            cells[numberOfCells-1][4*i+3][0] = 1.;
          }
          else
          {
            cells[numberOfCells-1][4*i+1][0] = 1.;
            cells[numberOfCells-1][4*i+2][0] = 1.;
          }
        }
        else
        {
          InStream.LoadValue(2,tmp);
          if(tmp == 0)
          {
            cells[numberOfCells-1][4*i][0] = 1.;
            cells[numberOfCells-1][4*i+1][0] = 1.;
          }

          if(tmp == 1)
          {
            cells[numberOfCells-1][4*i+2][0] = 1.;
            cells[numberOfCells-1][4*i+3][0] = 1.;
          }

          if(tmp == 2)
          {
            cells[numberOfCells-1][4*i+0][0] = 1.;
            cells[numberOfCells-1][4*i+2][0] = 1.;
          }

          if(tmp == 3)
          {
            cells[numberOfCells-1][4*i+1][0] = 1.;
            cells[numberOfCells-1][4*i+3][0] = 1.;
          }
        }
      }
      else if(cells[numberOfCells-2][i][0] == 3.)
      {
        cells[numberOfCells-1][4*i+0][0] = 1.;
        cells[numberOfCells-1][4*i+1][0] = 1.;
        cells[numberOfCells-1][4*i+2][0] = 1.;
        cells[numberOfCells-1][4*i+3][0] = 1.;

        InStream.LoadValue(2,tmp);
        cells[numberOfCells-1][4*i+tmp][0] = 0.;
      }
      else if(cells[numberOfCells-2][i][0] == 4.)
      {
        cells[numberOfCells-1][4*i+0][0] = 1.;
        cells[numberOfCells-1][4*i+1][0] = 1.;
        cells[numberOfCells-1][4*i+2][0] = 1.;
        cells[numberOfCells-1][4*i+3][0] = 1.;
      }
    }
  }


  Reshape(NbRows,3);
  SetValues(0.);
  int index = 0;
  int x, y, f;
  int addx, addx0, addy, addy0, addf, addf0;
  int localcounter = numberOfCells - localquantbits;
  for(int k=0;k<cells[numberOfCells-1].GetNbRows();k++)
  {
    if(cells[numberOfCells-1][k][0] != 0)
    {
      x = 0;
      y = 0;
      f = 0;
      for(int l = 0; l<(int)numberOfCells;l++)
      {
        if (l < localcounter)
        {
          addx = k & (1 << (1 + 2*l));
          addx0 = (addx >> (1 + l));
          x += addx0;

          addy = k & (1 << (2*l));
          addy0 = (addy >> (l));
          y += addy0;
        }
        else
        {
          addx = k & (1 << (2 + 3*(l-localcounter) + 2*localcounter));
          addx0 = (addx >> (2 + 2*l - localcounter));
          x += addx0;

          addy = k & (1 << (1 + 3*(l-localcounter) + 2*localcounter));
          addy0 = (addy >> (1 + 2*l - localcounter));
          y += addy0;

          addf = k & (1 << (3*(l-localcounter) + 2*localcounter));
          addf0 = addf >> (2*l);
          f += addf0;
        }
      }
      (*this)[index][0] = x;
      (*this)[index][1] = y;
      (*this)[index][2] = f*quantization + quantization/2;

      index++;

    }
  }

  for (unsigned int i=0;i<numberOfCells;i++)
    cells[i].Free();
    
  delete [] cells;
  cells = NULL;
    
  Vector.Free();
}


void M3Matrix::LoadCube(M3Matrix& Vector, unsigned int Max, RBitStream& InStream)
{
  unsigned int CurrentMax = Max;
  M3Matrix sum2(4,1);
  M3Matrix sum4(2,1);

  unsigned int tmp;
  InStream.LoadValue(M3Matrix::getNumberOfBitsNeeded(CurrentMax),tmp);
  sum4[0][0] = tmp;
  sum4[1][0] = (int)Max-(int)tmp;

  int k;

  for(k=0;k<2;k++)
  {
    CurrentMax = DBL2INT(sum4[k][0]);
    if(CurrentMax>0)
    {
      InStream.LoadValue(M3Matrix::getNumberOfBitsNeeded(CurrentMax),tmp);
      sum2[2*k][0] = tmp;
      sum2[2*k+1][0] = (int)CurrentMax-(int)tmp;
    }
    else
    {
      sum2[2*k][0] = 0;
      sum2[2*k+1][0] = 0;
    }
  }

  for(k=0;k<4;k++)
  {
    CurrentMax = DBL2INT(sum2[k][0]);
    if(CurrentMax>0)
    {
      InStream.LoadValue(M3Matrix::getNumberOfBitsNeeded(CurrentMax),tmp);
      Vector[2*k][0] = tmp;
      Vector[2*k+1][0] = (int)CurrentMax-(int)tmp;
    }
    else
    {
      Vector[2*k][0] = 0;
      Vector[2*k+1][0] = (int)CurrentMax-(int)0;
    }
  }
}

void M3Matrix::LoadCube4(M3Matrix& Vector,unsigned int Max,RBitStream& InStream)
{
  unsigned int CurrentMax = Max;
  M3Matrix sum4(2,1);

  unsigned int tmp;
  InStream.LoadValue(M3Matrix::getNumberOfBitsNeeded(CurrentMax),tmp);
  sum4[0][0] = tmp;
  sum4[1][0] = (int)Max-(int)tmp;

  for(int k=0;k<2;k++)
  {
    CurrentMax = DBL2INT(sum4[k][0]);
    if(CurrentMax>0)
    {
      InStream.LoadValue(M3Matrix::getNumberOfBitsNeeded(CurrentMax),tmp);
      Vector[2*k][0] = tmp;
      Vector[2*k+1][0] = (int)CurrentMax-(int)tmp;
    }
    else
    {
      Vector[2*k][0] = 0;
      Vector[2*k+1][0] = 0;
    }
  }
}

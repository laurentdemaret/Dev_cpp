#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <iostream>
#include <fstream>

#include "bitstr.h"

using std::ostream;

#define DBL2INT(x)       ((int)(x+0.000000000001))
#define DBL_ROUND(x)     ((double)(int)(x+0.000000000001))
#define COORD_ROUND(x)   ((int)(x+0.0011))

#ifdef __STDC__
typedef	signed char int8;	/* NB: non-ANSI compilers may not grok */
#else
typedef	char int8;
#endif

typedef	unsigned char uint8;
typedef	short int16;
typedef	unsigned short uint16;	/* sizeof (uint16) must == 2 */

#if defined(__alpha) || (defined(_MIPS_SZLONG) && _MIPS_SZLONG == 64)
typedef	int int32;
typedef	unsigned int uint32;	/* sizeof (uint32) must == 4 */
#else
typedef	long int32;
typedef	unsigned long uint32;	/* sizeof (uint32) must == 4 */
#endif

using std::cout;

class M3Matrix
{
 public:

  // Constructors
  M3Matrix();
  M3Matrix(int NbR, int NbC);
  M3Matrix(int NbRows, int NbCols, const double& r_Value);

  // Copy
  M3Matrix(const M3Matrix& r_Source);

  // Destructor
  virtual ~M3Matrix();

  // Accessors
  int GetNbRows();
  int GetNbCols();

  double* operator[](int Row);

  // Extraction of rows and columns
  void SetLine(int j, M3Matrix& Col);
  void SetCol(int j, M3Matrix& Col);

  // Operators
  M3Matrix& operator=(const M3Matrix& ro_ToAssign);
  M3Matrix* operator+(M3Matrix& r_MRight);
  M3Matrix* operator*(M3Matrix& r_MRight);

  // Transformers
  void SetValues(const double &r_Value);
  void Reshape(int NbRows, int NbCols);

  void ExchangeRows(int i1, int i2);
  void ExchangeCols(int i1, int i2);

  // Sorting functions
  void SortByCol(int c);
  void SortByCols(int c1, int c2);
  void SortByBlocs(int c1, int c2);

  // Standard computation tools
  double MSE(M3Matrix& M3Matrix2);
  double PSNR(M3Matrix& M3Matrix2);
  double Histogram(M3Matrix& Histo);

  // Loading the M3Matrix from a file
  void LoadFromFile(const char* InputFile);
  void LoadFromFileWithoutHeader(const char* InputFile);
  void LoadYFromPGM(const char* infile);

  // Gets the minimal value of the array
  double GetMin();
  double GetMax();

  // Quantization of the coefficients
  void Quantize(double& QStep, M3Matrix& Symbol);

  // Display
  void Dump(ostream& Out = cout);
  void SavePGM(ostream& Out);

  void EncodeOctTree(WBitStream& OutStream,int NbRows,int NbCols,int quantization);
  void TransmitCube(M3Matrix& Vector, unsigned int Max, WBitStream& OutStream);
  void TransmitCube4(M3Matrix& Vector, unsigned int Max, WBitStream& OutStream);

  void DecodeOctTree(int* rows,int* cols, RBitStream& InStream, int* quant);
  void LoadCube(M3Matrix& Vector, unsigned int Max, RBitStream& InStream);
  void LoadCube4(M3Matrix& Vector,unsigned int Max,RBitStream& InStream);


  // returns number of bits needed to represent a given number
  static unsigned int getNumberOfBitsNeeded(unsigned int number) {
    unsigned int numberOfBits=1;
    while (number>1) { number>>=1; numberOfBits++; }
    return numberOfBits;
  }

 protected:
  int t_NbRows;
  int t_NbCols;
  double **Data;

 private:
  // Memory deallocation
  void Alloc();
  void Free();

  // Sorting tools
  int vCompareFirstCol(const void *a1, const void *a2);
};

#endif /* of _MATRIX_H_ */

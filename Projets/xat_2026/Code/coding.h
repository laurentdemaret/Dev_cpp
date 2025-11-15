#ifndef _CODING_H_
#define _CODING_H_

#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>

/*#include "../triangulation/Point2D.h"
#include "../tools/matrix/matrix.h"
#include "../tools/matrix/bitstr.h"*/

#include "Point2D.h"
#include "matrix.h"
#include "bitstr.h"


using std::cout;
using std::endl;
using std::vector;

class Coding {

 private:

  Coding();

 public:

  ~Coding();

  static vector<Point2D*> Decode(int* NbRows, 
  					 int* NbCols, 
  					 const char* inputFilename) {
    RBitStream* pInStream = new RBitStream(inputFilename);
	vector<Point2D*> nodes;
    M3Matrix Decode;
    int quantization;
    Decode.DecodeOctTree(NbRows, NbCols, *pInStream, &quantization);
#ifndef NO_DEBUG_MESSAGES
    cout << "quantization: " << quantization << endl;
#endif

	for (int i=0;i<Decode.GetNbRows();i++) {
		Point2D* p = Point2D::makePoint2D(Decode[i][0],Decode[i][1],Decode[i][2]);
		nodes.push_back(p);
    }

#ifndef NO_DEBUG_MESSAGES
    cout << "size of p2da in decode is " << nodes.size() << endl;
#endif

    delete pInStream;
    
    return nodes;
  }


  static void Encode(int NbRows, 
  					 int NbCols, 
  					 vector<Point2D*> nodes, 
  					 const char* outputFilename, 
  					 int quantization) {
    M3Matrix Vertices;

    Vertices.Reshape(nodes.size(),3);

    for (unsigned int i=0;i<nodes.size();i++) {
      Vertices[i][0] = nodes[i]->x;
      Vertices[i][1] = nodes[i]->y;
      Vertices[i][2] = nodes[i]->f;
      if (Vertices[i][2]<0.)   Vertices[i][2] = 0.;
      if (Vertices[i][2]>255.) Vertices[i][2] = 255.;
    }

    WBitStream* pOutStream = new WBitStream(outputFilename);

    Vertices.EncodeOctTree(*pOutStream,NbRows,NbCols,quantization);
    for (int co=0;co<8;co++)
      pOutStream->WriteBit((unsigned char)0);

    pOutStream->Flush();

    delete pOutStream;
  }


  // Coding of the Value where we know that MaxValue is the number of possible values
  static void CodeHuffmann(int Value, int MaxValue, WBitStream& OutStream) {
    int NbBits = M3Matrix::getNumberOfBitsNeeded(MaxValue);
    int Rest =  MaxValue % ((int)pow(2.,(double)NbBits));

    int s = MaxValue - 2*Rest;
    if (Value<s)
      OutStream.SaveValue(NbBits,Value);
    else {
	  int Symbol = (int)pow(2.,(double)NbBits+1) - (MaxValue - Value);
      OutStream.SaveValue(NbBits, Symbol/2);
      OutStream.SaveValue(1, Symbol%2);
    }
  }


  static int DecodeHuffmann(int MaxValue, RBitStream& InStream) {
    int NbBits = M3Matrix::getNumberOfBitsNeeded(MaxValue);
    int Rest =  MaxValue % ((int)pow(2.,(double)NbBits));

    int s = MaxValue -2*Rest;
    unsigned int DecodedValue = 0;
    InStream.LoadValue(NbBits,DecodedValue);

    if ((int)DecodedValue >= s) {
      unsigned int Value = 0;
      InStream.LoadValue(1,Value);
      int DecodedSymbol = (unsigned int)(DecodedValue*2 + Value);
      DecodedValue = MaxValue - ((int)pow(2.,(double)NbBits + 1) - DecodedSymbol);
    }

    int Decoded = (int)DecodedValue;

    return Decoded;
  }

};

#endif /* of _CODING_H_ */

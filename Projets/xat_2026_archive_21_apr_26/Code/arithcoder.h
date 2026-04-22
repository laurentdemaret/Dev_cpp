#ifndef _ARITHCODER_H_
#define _ARITHCODER_H_

#include "qsmodel.h"
#include "rangecod.h"
#include "bitstr.h"

#define BASE_2_LOG_FREQ_COUNT   12 // base2 log of total frequency count
#define RESCALING_INTERVAL      2000


class ArithmCoder
{
 public:
  // Constructor
  ArithmCoder(unsigned int NbOfSymbols);

  // Destructor
  virtual ~ArithmCoder();

  // Accessors
  unsigned int GetNbEncodedSymb();

  // Encoding functions
  void StartEncoding();
  void EncodeValue(unsigned int value, WBitStream& OutBitStream);
  void DoneEncoding(WBitStream& OutBitStream);

  // Decoding functions
  void StartDecoding(RBitStream& ro_InBitStream);
  bool DecodeValue(unsigned int& rui_Value, RBitStream& InBitStream);
  void DoneDecoding(RBitStream& ro_InBitStream);


 protected:
  unsigned int NbSymbols;  // Does not include EOF
  unsigned int NbEncodedSymbols;
  // Data
  rangecoder rc;
  qsmodel qsm;
  int syfreq, ltfreq;

};

#endif

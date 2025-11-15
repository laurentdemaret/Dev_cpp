//  ARITHMETIC ENCODING ALGORITHM

#include "qsmodel.h"
#include "rangecod.h"
#include "bitstr.h"
#include "arithcoder.h"

using namespace std;


// Constructor with number of symbols
ArithmCoder::ArithmCoder(unsigned int NbSymb)
{
  NbSymbols = NbSymb;
  NbEncodedSymbols = 0;
}


// Destructor
ArithmCoder::~ArithmCoder()
{
}


unsigned int ArithmCoder::GetNbEncodedSymb()
{
  return NbEncodedSymbols;
}


void ArithmCoder::StartEncoding()
{
  // We need NbSymbols+1 because of EOF
  initqsmodel(&qsm,NbSymbols+1,BASE_2_LOG_FREQ_COUNT,RESCALING_INTERVAL,NULL,1);
  start_encoding(&rc,0);
}


void ArithmCoder::EncodeValue(unsigned int value, WBitStream& OutBitStream)
{
  qsgetfreq(&qsm,value,&syfreq,&ltfreq);
  encode_shift(&rc,syfreq,ltfreq,BASE_2_LOG_FREQ_COUNT,OutBitStream);
  qsupdate(&qsm,value);
  NbEncodedSymbols++;
}


void ArithmCoder::DoneEncoding(WBitStream& OutBitStream)
{
  qsgetfreq(&qsm,NbSymbols,&syfreq,&ltfreq);
  encode_shift(&rc,syfreq,ltfreq,BASE_2_LOG_FREQ_COUNT,OutBitStream);
  done_encoding(&rc,OutBitStream);
  deleteqsmodel(&qsm);
  NbEncodedSymbols++;
}


void ArithmCoder::StartDecoding(RBitStream& InBitStream)
{
  // NbSymbols+1 because of EOF
  initqsmodel(&qsm,NbSymbols+1,BASE_2_LOG_FREQ_COUNT,RESCALING_INTERVAL,NULL,0);
  start_decoding(&rc,InBitStream);
}


bool ArithmCoder::DecodeValue(unsigned int& Value, RBitStream& InBitStream)
{
  ltfreq = decode_culshift(&rc,BASE_2_LOG_FREQ_COUNT,InBitStream);
  Value = qsgetsym(&qsm,ltfreq);
  NbEncodedSymbols++;
  if (Value == NbSymbols)  // check for end-of-file 
    return false;

  qsgetfreq(&qsm,Value,&syfreq,&ltfreq);
  decode_update(&rc,syfreq,ltfreq,1<<BASE_2_LOG_FREQ_COUNT);
  qsupdate(&qsm,Value);

  return true;
}


void ArithmCoder::DoneDecoding(RBitStream& InBitStream)
{
  qsgetfreq(&qsm,NbSymbols,&syfreq,&ltfreq);
  decode_update(&rc,syfreq,ltfreq,1<<BASE_2_LOG_FREQ_COUNT);
  done_decoding(&rc,InBitStream);
  deleteqsmodel(&qsm);
}

#ifndef _IO_H_
#define _IO_H_

#include <string>
#include <vector>


#define EXT_UNKNOWN (-1)
#define EXT_NOEXTENSION (0)
#define EXT_PGM (1)
#define EXT_AT (2)
#define EXT_SEQ (3)
#define EXT_INSIDE (4)
//#define PARSE_ARGUMENT(X)

std::vector<int>* parse_multiple_iterations(int argc,char*const* argv);
int recognizeExtension(const char* filename);
std::string stripExtension(const char* filename);
char** parse_command_line_at_and_compress(int, char*const*, int*, int*, int*, int*,int*, int*, char**, int*);
void print_usage();

#endif

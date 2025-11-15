#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include "IO.h"

using namespace std;

int recognizeExtension(const char* filename) 
{
  if (filename!=NULL) 
  {
    const char* extension=strrchr(filename,'.');
    if (extension==NULL) return EXT_NOEXTENSION;
    if (strcasecmp(extension,".pgm")==0) return EXT_PGM;
    if (strcasecmp(extension,".at")==0) return EXT_AT;
    if (strcasecmp(extension,".seq")==0) return EXT_SEQ;
    if (strcasecmp(extension,".outside")==0) return EXT_SEQ;
		if (strcasecmp(extension,".inside")==0) return EXT_INSIDE;
	}
  return EXT_UNKNOWN;
}


string stripExtension(const char* filename) 
{
  if (filename!=NULL) 
  {
    const char* extension=strrchr(filename,'.');
    if (extension==NULL) return string(filename);
    else {
      // strndup is not standard :(
      char* tmp = (char*)malloc(extension-filename+1);
      memcpy(tmp,filename,extension-filename);
      tmp[extension-filename] = '\0';
      string ret(tmp);
      free(tmp);
      return ret;
    }
  }
  return "";
}


extern void print_usage();


#define PARSE_ARGUMENT(X) \
if ((*opt)=='=' || (*opt)==':') opt++; \
else if ((*opt)=='\0') opt=argv[++optind]; \
if (optind<argc) { tmp = strtol(opt,&stop_char,0); } \
if ((optind>=argc) || !(*opt!='\0' && *stop_char=='\0')) { \
  cerr << endl << "Error: option " << X << " requires an integer argument, aborting..." << endl << endl; \
  free(filenames); \
  exit(EXIT_FAILURE); \
} else


/** The following function is used for test purposes, i.g. when you have a .seq
 * file and you want to compress a picture in multiple bpp.
 * see README for AT */
vector<int>* parse_multiple_iterations(int argc,char*const* argv)
{
  char** filenames = (char**)malloc((argc-1)*sizeof(char*));
  vector<int>* result =new vector<int>;
  int optind = 1;
  char* opt,*stop_char;
  int tmp;
  while (optind<argc) 
	{
    opt=argv[optind];
    if (*opt == '-' && *(++opt) == 'i'){
      opt++;
    	PARSE_ARGUMENT("-i") result->push_back(tmp);
    }
    optind++;
  }
  return result;
}


char** parse_command_line_at_and_compress(int argc, 
											char*const* argv, 
											int* at_alg, 
											int* iterations, 
											int* exchange_iterations, 
											int* exchange_radius, 
											int* quantization, 
											int* compressing_options, 
											char** action_key,
											int* filenames_num) 
										{

  *filenames_num=0;
  *at_alg=-1;
  *iterations=-1;
  *quantization=-1;
  *compressing_options=0;
  *exchange_iterations = -1;
  *exchange_radius = 1;
  *action_key = (char*)"";

  char** filenames = (char**)malloc((argc-1)*sizeof(char*));

  int optind=1, tmp;
  char *opt, *stop_char;

  while (optind<argc) 
	{
    tmp=-1;
    stop_char=0;
    opt=argv[optind];
    if (*opt == '-') /* an option */ 
		{
      opt++; /* skip trailing '-' */
      switch (*opt) 
			{
      case 'a':
      case 'A':
		opt++; if ((*opt)=='t' || (*opt)=='T') opt++; /* recognize also -at */
		PARSE_ARGUMENT("-a[t]")	*at_alg=tmp;
      	break;
      case 'i':
      case 'I':
		opt++;
		PARSE_ARGUMENT("-i") *iterations=tmp;
		break;
      case 'S':
      case 's':
		opt++;
		PARSE_ARGUMENT("-s") *exchange_iterations=tmp;
		break;
      case 'M':
      case 'm':
		opt++;
		*action_key=(char*)"reconstruct";
		break;
      case 'r':
      case 'R':
		opt++;
		PARSE_ARGUMENT("-r") *exchange_radius=tmp;
		break;
      case 'q':
      case 'Q':
		opt++;
		PARSE_ARGUMENT("-q") *quantization=tmp;
		break;
      case '\0': cerr << "Warning: no option specified after '-', ignoring..." << endl; break;
      default:
	while (*opt!='\0') {
	  switch (*opt) {
	  case '?':
	  case 'h': print_usage(); free(filenames); exit(EXIT_SUCCESS); break;
/*	  
	  case 'y': *compressing_options|=PRINT_ALL_POINTS_L;           break;
	  case 'z': *compressing_options|=PRINT_ALL_POINTS_S;           break;
	  case 'l': *compressing_options|=PRINT_NODES_L;                break;
	  case 's': *compressing_options|=PRINT_NODES_S;                break;
	  case 'o': *compressing_options|=PRINT_THINNED;                break;
	  case 'e': *compressing_options|=PRINT_EDGES;                  break;
	  case 't': *compressing_options|=PRINT_TRIANGLES;              break;
	  case 'p': *compressing_options|=PRINT_POLY;                   break;
	  case 'C': *compressing_options|=USE_CONVEX;                   break;
	  case 'x': *compressing_options|=PRINT_PERFORMANCE;            break;
*/	  
	  default:
	    cerr << endl << "Error: unrecognized option '-" << *opt << "', aborting..." << endl << endl;
	    print_usage(); free(filenames); exit(EXIT_FAILURE);
	    break;
	  }
	  opt++;
	}
	break;
      }
    }
    else /* a filename */ 
		{
      filenames[*filenames_num]=opt;
      (*filenames_num)++;
    }
    optind++;
  }
  return filenames;
}


void print_usage() 
{
  cerr << endl
       << "Usage: xat [-a[t]INT] [-iINT] [-qINT] [-sINT] [-hyzloetpCx] <input.pgm [output.at]>..." << endl
       << "   -a[t]INT   Adaptive thinning algorithm to be used (default=AT5)" << endl
       << "   -iINT      Number of iterations to be done (default=90% of number of points)" << endl
       << "   -qINT      Quantization to be used (default=4)" << endl
      	<< "  -sINT      Number of exchange iterations" << endl
       << "   -h         Print this help and exit" << endl
       << "   -y         Print all points and their luminance values (default=off)" << endl
       << "   -z         Print all points and their significances (default=off)" << endl
       << "   -l         Print only the nodes and their luminance values (default=off)" << endl
       << "   -o         Print only the points thinned out (default=off)" << endl
       << "   -e         Print all edges of the triangulation (default=off)" << endl
       << "   -t         Print all triangles of the triangulation (default=off)" << endl
       << "   -p         Create poly file (default=off)" << endl
       << "   -C         Use convex triangulation (FIXME) (default=off)" << endl
       << "   -x         Print performance data (default=off)" << endl
	     << "Example of use: " << endl
	     << "./xat -q8 -at6  -i14000 -s5000 -r3  Data/lena128.pgm" << endl	
	<< endl;
}


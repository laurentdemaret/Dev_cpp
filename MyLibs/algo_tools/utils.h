#ifndef UTILS_H_
#define UTILS_H_

#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <map>
#include <cctype>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <fcntl.h>
#include <cstring>
#include <complex>
#include <valarray>
#include <time.h>
#include <chrono>
#include <deque>

#define TWO_PI 6.2831853071795864769252866
#define PI     3.14159265358979323846264338

#define COORD_ROUND(x)   ((int)(x+0.0011))

using namespace std;
using namespace std::chrono;

typedef high_resolution_clock Clock;
typedef Clock::time_point ClockTime;

//**************
//*** SYSTEM ***
//**************
string GetHomeDir();
string GetCurrentWorkingDir();
int system(string cmd_line);
void add_lib_path(string new_path_lib_path);


//************
//*** TIME ***
//************
void get_start_time(struct timespec& Start);
double get_elapsed_time(struct timespec Start);


//****************************
//* DOUBLE POINTER ROUTINES
//****************************
vector<double> operator+( const vector<double>& v1, const vector<double>& v2 );
vector<double> operator-( const vector<double>& v1, const vector<double>& v2 );
vector<double> operator*(  const vector<double>& v1, const double& lambda);
vector<double> operator*(  const double& lambda, const vector<double>& v1);

//computes the mean of a double ptr
double vec_mean(double* x, int Length);

//dumps a double ptr
void dump_dblptr(double* x, int Length, std::string first_line_str ="");


//********************************
//* Complex vectors and matrices
//********************************
void conjugate(vector<vector<complex<double>>>& x);
vector<complex<double>> conj(vector<complex<double>>& x);
vector<double> abs(vector<complex<double>>& x);
vector<vector<complex<double>>> conj(vector<vector<complex<double>>>& x);
vector<vector<double>> abs(vector<vector<complex<double>>>& x);
void abs(vector<vector<complex<double>>>& x, vector<vector<double>>& y);

vector<complex<double>> operator /(const vector<complex<double>>& v1, const double& lambda);

vector<vector<complex<double>>> hadamardprod(vector<vector<complex<double>>>& x, vector<vector<complex<double>>>& y);
// ... and its "efficient" variant. modifies the buffer x
void hadamardprod_na(vector<vector<complex<double>>>& x, vector<vector<complex<double>>>& y);

//Warning : this is the pointwise complex normalisation
vector<vector<complex<double>>> normalise_pw(vector<vector<complex<double>>>& x);
// ... and its "efficient" variant. modifies the buffer
void normalise_pw_na(vector<vector<complex<double>>>& x);


//***********************
//* Vector templates
//***********************
template<class T>
vector<T> subvector(vector<T>& v1, vector<int>& indices);

template<class T>
vector<T> subvector(vector<T>& v1, slice& s);

//template <class T>
void resize(vector<vector<double>>& mat, int m, int n);
void resize(vector<vector<complex<double>>>& mat, int m, int n);


//****************************
//* STRING TOOLBOX
//****************************
//adds numbers to a string
std::string append(std::string str, int i);
std::string int2str(int i);
string append(string str, double x, int precision = 2);
double string2double(char* ch);
int str2int(string s);
double str2double(string s);

string double2str(double x, int precision = 2);

string append(string str, double x, int precision);
string string_replace(string src, string const& target, string const& repl);

bool compare_begin(string str1, string str2);
bool space(char c);
bool not_space(char c);


vector<std::string> split(const string& str);
vector<string> split(string str, string delimiter);

vector<string> split(const string& str);
void remove_stop_words(vector<string>& vec, vector<string> stop_words);
void dump(vector<string> vec);
void dump2(vector<string> vec, vector<int> freq);
void set_stop_words(vector<string>& stop_words);
void compute_frequencies(vector<string> input_vec, vector<string> most_freq_str, vector<int>& frequencies);
void remove_punctuation_and_split(string str,vector<string>& str_split);
void load_text(ifstream& Text, vector<string>& text_vec);
void load_text(string filename, vector<string>& text_vec);
void clean_and_load_text(ifstream& Text, vector<string>& text_vec);
void write_text(ofstream& TextFile, vector<string> text_vec);
void write_text(string filename,vector<string> text_vec);

//Filename handling
void fileparts(string filename, string& path, string& name, string& ext);
string extract_path(string filename);
string extract_name(string filename);
string extract_ext(string filename);
void dump_extension_path_name(string filename);

string remove_extension(string filename);
string replace_extension(string filename, string new_ext);


// *******************************
// Parsing Tools
// *******************************
//Standard surface types used by generate_synthetic_OFF
enum OPTION_TYPE{stand_alone,int_params,double_params,string_params};

class Option
{
public:
  Option() {};
  virtual ~Option() {};

  //Members
  string option_name;
  vector<string> delimitators;
  int NbParams;
  bool opt_detected = false;
  OPTION_TYPE type;

  //not very beautiful but this works !
  vector<int> IntValues;
  vector<double> DoubleValues;
  vector<string> StringValues;
};


class Parser
{
public:
  Parser() {};
  virtual ~Parser() {};

  void add_option_stand_alone(string parsename, vector<string> delim_vec);
  void add_option_stand_alone(string option_name, string delim);

  void add_option_string(string parsename, vector<string> delim, int NbParams=1);
  void add_option_int(string parsename, vector<string> delim, int NbParams=1);
  void add_option_double(string parsename, vector<string> delim, int NbParams=1);

  void add_option_string(string parsename, string delim, int NbParams=1);
  void add_option_int(string parsename, string delim, int NbParams=1);
  void add_option_double(string parsename, string delim, int NbParams=1);

  bool get_int_param(string option_name, int& option_value);
  bool get_double_param(string option_name, double& option_value);
  bool get_string_param(string option_name,string& option_value);

  bool get_int_params(string option_name, vector<int>& option_vec);
  bool get_double_params(string option_name, vector<double>& option_vec);
  bool get_string_params(string option_name, vector<string>& option_vec);

  bool check_option(string option_name);

  void dump_options_and_params();

  void parse_command_line(int argc, char*const* argv);
  vector<Option> options;

  int nb_stand_alone_options = 0;
  int nb_string_options = 0;
  int nb_double_options = 0;
  int nb_int_options = 0;
};


//****************************
//* FLOATING BASIC TOOLS
//****************************
double truncate(double x, int nb_dec);
double myclip(double x, double min, double max);
double sign(double x);

double min3(double a, double b, double c);
double min4(double a, double b, double c, double d);

double max3(double a, double b, double c);
double max4(double a, double b, double c, double d);


//*****************************
//* COMPLEX MANIPULATIONS
//*****************************
void dump_complex(complex<double> , string delim="\n", int precision = 5);

//****************************
//* VECTOR MANIPULATIONS
//****************************
void dump_vec(vector<double> vec, string message = "", int precision = 5);
void dump_vec(vector<vector<double>> vec, string message = "", int precision = 5);
void dump_vec(vector<complex<double>> vec, string message ="", string delim="\n", int precision = 5);
void dump_vec(vector<vector<complex<double>>> array, string message="", string delim="\n", int precision=5);

//Minimum and maximum of vectors
double findMin(vector<double> vec,int& MinIndx);
double findMax(vector<double> vec,int& MaxIndx);
double Mean(vector<double> vec);
double Var(vector<double> vec);
double Norm(vector<double> vec, double p = 2.0);


//****************************
//* RANDOM GENERATOR
//****************************
//Generate a random variable with uniform distribution on [0,1]
double rand_unif(int seed = -1);

//Generate a random integer between 1 and Max (default : simulated ouptut of a dice)
int rand_int(int Max = 6, int seed = -1);

//Generate gaussian random noise
double randn(int seed = -1);
vector<double> randn(int size, int seed);

int my_cmp(const void *x, const void *y);
int my_cmp_descent(const void *x, const void *y);

void WriteTXT( vector<double>& Values, string& filenameOut);
void skip_comment_lines(ifstream& In, char comment_symbol =  '#');

void writeIntTo4Bytes(ofstream& OutStr, int x);
int readIntFrom4Bytes(istream& InStr);

//Comparison tools
int fncompare(const void * el1,const void * el2);
int fncompare2(const void * el1,const void * el2);
int fncomparebloc(const void * el1, const void * el2);

double spherical_angle(double theta0, double phi0,double theta1, double phi1);

//Computes manifold distance between two points on a cylinder, given in cylinder coordinates
double dist_cyl(double theta0, double z0, double theta1, double z1);

//Two useful functions to get radius and angle (polar coordinates) from cartesian coordinates
double Radius(double x, double y);
double Theta(double x, double y);

bool is_power_of_two(int n);

void printExecutionTime(ClockTime start_time, ClockTime end_time);


string date_stamp();
string date_hour_stamp();

//Test routines
void test_string();
void test_formats();
void test_fileparts();
void test_precision();
void test_map();
void test_write_char_ptr();
void test_slice();
void test_std_complex();
void test_time();
void test_deque();

#endif

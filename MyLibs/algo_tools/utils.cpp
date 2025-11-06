#include "utils.h"
#include <cmath>

//Warning: this will work only in LINUX/UNIX systems !
#include <unistd.h>
#include <ctime>

//*************************************
// Author: Laurent Demaret, 2014-2025
//*************************************
#define num_tol 0.000001
#define display_tol 0.00001

//*******************************************************
// System utilities
// Useful for command line calls inside of C++ programs
//*******************************************************
#define GetCurrentDir getcwd


int system(string cmd_line)
{
    int ret = system(cmd_line.c_str());
    return ret;
}


string GetHomeDir()
{
  char *homedir = getenv("HOME");
  string homedir_str(homedir);
  return homedir_str;
}


string GetCurrentWorkingDir()
{
    char buff[FILENAME_MAX];
    string current_working_dir("");
    if( GetCurrentDir( buff, FILENAME_MAX ) )
      current_working_dir.assign(buff);
    return current_working_dir;
}


//Adds the given path to the LD_LIBRARY_PATH environment variable
void add_lib_path(string new_path_lib_path)
{
    string ld_library_path = std::getenv("LD_LIBRARY_PATH");
    //cout << "ld_library_path : " << ld_library_path << endl;
    new_path_lib_path += ":"+ld_library_path;
    setenv("LD_LIBRARY_PATH",new_path_lib_path.c_str(),1);
    //cout << "LD_LIBRARY_PATH is: " << getenv("LD_LIBRARY_PATH") << endl;
}


//*****************************
// Vector operators
//*****************************
vector<double> operator+( const vector<double>& v1, const vector<double>& v2 )
{
    vector<double> v3;
    v3.resize(min(v1.size(),v2.size()));

    for(int i = 0; i< int(v3.size()); i++)
        v3[i] = v1[i] + v2[i];

    return v3;
}


vector<double> operator-( const vector<double>& v1, const vector<double>& v2 )
{
    vector<double> v3;
    v3.resize(min(v1.size(),v2.size()));

    for(int i = 0; i< int(v3.size()); i++)
        v3[i] = v1[i] - v2[i];

    return v3;
}


vector<double> operator*( const vector<double>& v1, const double& lambda)
{
    vector<double> v2(v1.size());
    for(int i =0; i < int(v1.size());i++)
        v2[i] = lambda * v1[i];

    return v2;
}


vector<double> operator*(  const double& lambda, const vector<double>& v1)
{
    vector<double> v2(v1.size());
    for(int i =0; i < int(v1.size());i++)
        v2[i] = lambda * v1[i];

    return v2;
}


vector<complex<double>> operator /(const vector<complex<double>>& v1, const double& lambda)
{
    vector<complex<double>> v2(v1.size());
    if(lambda != 0.)
    {
        for(int i =0; i < int(v1.size());i++)
            v2[i] =  v1[i]  / lambda;
    }

    return v2;
}


//Pointwise multiplication of two complex "matrices"
//Warning: user of this function must ensure that x and y have same sizes
vector<vector<complex<double>>> hadamardprod(vector<vector<complex<double>>>& x, vector<vector<complex<double>>>& y)
{
    vector<vector<complex<double>>> z;
    z.resize(x.size());

    for(int i = 0;i<(int)x.size();i++)
    {
        for(int j = 0;j<(int)x[i].size();j++)
        {
           z[i].push_back(x[i][j] * y[i][j]);
        }
    }

    return z;
}


//Stores the result of the pointwise multiplication of both buffers (matrices, x and y) in the first buffer matrix (x)
void hadamardprod_na(vector<vector<complex<double>>>& x, vector<vector<complex<double>>>& y)
{
   //cout << "hadamardprod_na x.size() " << x.size()  << endl;
   //cout << "hadamardprod_na y.size() " << y.size()  << endl;

    for(int i = 0;i<(int)x.size();i++)
    {
        //cout << "hadamardprod_na x["<<i<<"].size() " << x[i].size()  << endl;
        //cout << "hadamardprod_na y["<<i<<"].size() " << y[i].size()  << endl;
        for(int j = 0;j<(int)x[i].size();j++)
        {
           x[i][j] = x[i][j] * y[i][j];
        }
    }
}



// This is the same as normalise_pw but directly stores the
// adapted for efficiency and potential parallelisation
void normalise_pw_na(vector<vector<complex<double>>>& x)
{
    double abs_val;
    for(int i = 0;i<(int)x.size();i++)
    {
        for(int j = 0;j<(int)x[i].size();j++)
        {
           abs_val =  std::abs(x[i][j]);
           if(abs_val > num_tol )
             x[i][j] = x[i][j]  / abs_val;
           else
             x[i][j]= 0.;
        }
    }

}


vector<vector<complex<double>>> normalise_pw(vector<vector<complex<double>>>& x)
{
    vector<vector<complex<double>>> xn;
    xn.resize(x.size());

    for(int i = 0;i<(int)x.size();i++)
    {
        for(int j = 0;j<(int)x[i].size();j++)
        {
           if(std::abs(x[i][j]) > num_tol )
             xn[i].push_back( x[i][j]  / std::abs(x[i][j]));
           else
               xn[i].push_back(0.0);
        }
    }

    return xn;
}


vector<vector<double>> abs(vector<vector<complex<double>>>& x)
{
    vector<vector<double>> y;
    //cout << "abs : x.size() "  << x.size() << endl;
    y.resize(x.size());
    for(int i = 0;i<(int)x.size();i++)
    {
        //cout << "abs : x[i].size() "  << x[i].size() << endl;
        for(int j = 0;j<(int)x[i].size();j++)
          y[i].push_back(std::abs(x[i][j]));
    }

    return y;
}


void abs(vector<vector<complex<double>>>& x, vector<vector<double>>& y)
{
    for(int i = 0;i<(int)x.size();i++)
    {
        //cout << "abs : x[i].size() "  << x[i].size() << endl;
        for(int j = 0;j<(int)x[i].size();j++)
          y[i][j] = std::abs(x[i][j]);
    }
}


void conjugate(vector<vector<complex<double>>>& x)
{
    for(int i = 0;i<(int)x.size();i++)
    {
        for(int j = 0;j<(int)x[i].size();j++)
          x[i][j] = std::conj(x[i][j]);
    }

}


vector<vector<complex<double>>> conj(vector<vector<complex<double>>>& x)
{
    vector<vector<complex<double>>> y;
    y.resize(x.size());
    for(int i = 0;i<(int)x.size();i++)
    {
        //more simply : y.push_back(conj(x[i]));
        for(int j = 0;j<(int)x[i].size();j++)
          y[i].push_back(std::conj(x[i][j]));
    }

    return y;
}


vector<complex<double>> conj(vector<complex<double>>& x)
{
    vector<complex<double>> y;
    for(int i = 0;i<(int)x.size();i++)
    {
        y.push_back(std::conj(x[i]));
    }
    return y;
}


vector<double> abs(vector<complex<double>>& x)
{
    vector<double> y;
    for(int i = 0;i<(int)x.size();i++)
    {
        y.push_back(std::abs(x[i]));
    }
    return y;
}


void resize(vector<vector<double>>& mat, int m, int n)
{
  if((int)mat.size() != m)
    mat.resize(m);

  for(int i = 0;i<m;i++)
  {
      if((int)mat[i].size() != n)
          mat[i].resize(n);
  }
}


void resize(vector<vector<complex<double>>>& mat, int m, int n)
{
  if((int)mat.size() != m)
    mat.resize(m);

  for(int i = 0;i<m;i++)
  {
      if((int)mat[i].size() != n)
          mat[i].resize(n);
  }
}


template <class T>
vector<T> subvector(vector<T>& v1, vector<int>& indices)
{
    vector<T> v2;
    for(int i =0; i < int(indices.size());i++)
    {
        if(indices[i] >= 0 && indices[i] < (int)v1.size())
            v2.push_back(v1[indices[i]]);
    }

    return v2;
}


template <class T>
vector<T> subvector(vector<T>& v1, slice& s)
{
    vector<T> v2;

    for(int i = 0;i<(int)s.size();i++)
    {
       int indx = s.start()+i*s.stride();
       if(indx >= 0 && indx < (int)v1.size())
           v2.push_back(v1[indx]);
    }

    return v2;
}


void dump_dblptr(double* x, int Length, string first_line_str)
{
    cout << first_line_str << endl;
    cout << std::setprecision(7);
    for(int i = 0;i<Length;i++)
    {
        cout << x[i] << endl;
    }
}


void dump_vec(vector<vector<double>> vec, string message, int precision)
{
    cout << std::fixed;
    cout << std::setprecision(precision);
    if(message != "")
        cout << message << " ";

    cout << "dump_vec :: vec.size() " << vec.size() <<endl;
    for(int i =0; i < int(vec.size());i++)
    {
        cout << "dump_vec :: vec[i].size() " << vec[i].size() <<endl;
      for(int j =0; j < int(vec[i].size());j++)
      {
        cout <<  std::setprecision(precision) << vec[i][j] << "  ";
      }
      cout << endl;
    }
}



void dump_vec(vector<double> vec, string message, int precision)
{
    cout << std::fixed;
    cout << std::setprecision(precision);
    if(message != "")
        cout << message << " ";

    cout << "(";
    for(int i =0; i < int(vec.size())-1;i++)
    {
        cout <<  std::setprecision(precision) << vec[i] << ", ";
    }
    cout <<  std::setprecision(precision) << vec[vec.size()-1] << ")" << endl;
}


void dump_complex(complex<double> z, string delim, int precision)
{
    //cout << std::fixed;
    cout << std::setprecision(precision);
    if(abs(z) < display_tol)
    {
        cout << "0.0 ";
    }
    else
    {
        cout << z.real();
        if(fabs(z.imag() ) > display_tol)
        {
            if(z.imag() > 0. )
                cout << " + ";
            else
                cout << " - ";

            cout << fabs(z.imag()) << " * i ";
        }
    }

    cout << delim;
}


void dump_vec(vector<complex<double>> vec, string message, string delim, int precision)
{
    //cout << std::fixed;
    cout << std::setprecision(precision);
    if(message != "")
        cout << message;

    //cout << "(";
    for(int i =0; i < int(vec.size());i++)
    {
        dump_complex(vec[i], delim, precision);
    }
    //cout <<  std::setprecision(precision) << vec[vec.size()-1] << ")" << endl;
}


void dump_vec(vector<vector<complex<double>>> array, string message, string delim, int precision)
{
    for(int i =0; i < int(array.size());i++)
    {
        cout << "Line : " << i << endl;
        dump_vec(array[i],message, delim, precision);
        cout << endl;
    }
}




//**************************
// Routines for strings
//**************************
string append(string str, int i)
{
    stringstream s;
    s  << str << i;
    return s.str();
}


string int2str(int i)
{
    stringstream s;
    s  << i;
    return s.str();
}


string double2str(double x, int precision)
{
    stringstream s;
    s << std::setprecision(precision);
    s << x;
    return s.str();
}


double string2double(char* ch)
{
    char* dummy;
    return strtod(ch,&dummy);
}


double str2double(string s)
{
    char* dummy;
    return strtod(s.c_str(),&dummy);
}


int str2int(string s)
{
    char* dummy;
    return strtol(s.c_str(),&dummy,10);
}



// Floating point replaced by an underscore
// Precision set to 2
string append(string str, double x, int precision)
{
    stringstream s;
    bool neg = false;
    if(x<0) {x=-x; neg = true;}

    int y = (int)(floor(x));
    int z;

    precision = 2;
    //if(precision == 2)
    z = int(floor((x-y)*pow(10,precision+1) + 0.5));


    s << str;
    if(neg == true)
    {
        s << "-";
    }
    s    << y << "_";
    if(z<10)
        s <<"00";
    else if(z<100)
        s <<"0";

    s << z;

    return s.str();
}


string string_replace( string src, string const& target, string const& repl)
{
    // handle error situations/trivial cases
    if (target.length() == 0)
    {
        // searching for a match to the empty string will result in an infinite loop
        //  it might make sense to throw an exception for this case
        return src;
    }

    if (src.length() == 0)
    {
        return src;  // nothing to match against
    }

    size_t idx = 0;

    for (;;) {
        idx = src.find( target, idx);
        if (idx == string::npos)  break;

        src.replace( idx, target.length(), repl);
        idx += repl.length();
    }

    return src;
}


bool compare_begin(string str1,string str2)
{
    size_t minlength = str1.size();
    if(str2.size()<minlength)
        minlength = str2.size();

    string str1_short = string(str1,0,minlength-1);
    string str2_short = string(str2,0,minlength-1);

    bool val; //short is equal to the first of the longest
    if(str1_short==str2_short)
        val= true;
    else
        val = false;

    return val;
}


// true if the argument is whitespace, false otherwise
bool space(char c)
{
    return isspace(c) != 0;
}


// false if the argument is whitespace, true otherwise
bool not_space(char c)
{
    return !isspace(c);
}


vector<string> split(string str, string delimiter)
{
    vector<string> internal;
    stringstream ss(str); // Turn the string into a stream.
    string tok;
    char delim = *(delimiter.c_str());
    while(getline(ss, tok, delim))
    {
        internal.push_back(tok);
    }

    return internal;
}


vector<string> split(const string& str)
{
    typedef string::const_iterator iter;
    vector<string> ret;
    iter i = str.begin();

    while (i != str.end())
    {
        // ignore leading blanks
        i = find_if(i, str.end(), not_space);
        // find end of next word
        iter j = find_if(i, str.end(), space);
        // copy the characters in [i, j)
        if (i != str.end())
            ret.push_back(string(i, j));
        i = j;
    }
    return ret;
}


void remove_stop_words(vector<string>& vec, vector<string> stop_words)
{
    int i = 0;
    while(i < int(vec.size()) )
    {
        bool erase = false;
        for(int j=0; j < int(stop_words.size());j++)
        {
            if( vec[i] == stop_words[j])
            {
                vec.erase(vec.begin() +i);
                erase = true;
            }
        }
        if(erase == false)
            i++;
    }
}


void dump(vector<string> vec)
{
    for(int i = 0;i < int(vec.size());i++)
    {
        cout << "[" << i << "] : " << vec[i] << endl;
    }
}


void dump2(vector<string> vec, vector<int> freq)
{
    for(int i = 0;i < int(vec.size());i++)
    {
        cout << "[" << i << "] : " << vec[i] << " : " << freq[i] << endl;
    }
}


void set_stop_words(vector<string>& stop_words)
{
    stop_words.push_back("du");
    stop_words.push_back("des");
    stop_words.push_back("les");
    stop_words.push_back("la");
    stop_words.push_back("le");
    stop_words.push_back("Le");
    stop_words.push_back("un");
    stop_words.push_back("une");
}


void compute_frequencies(vector<string> input_vec, vector<string> most_freq_str, vector<int>& frequencies)
{
    frequencies.resize(most_freq_str.size());
    for(int i = 0; i < int(frequencies.size());i++)
    {
        frequencies[i] = 0;
    }


    for(int i = 0;i < int(input_vec.size());i++)
    {
        for(int j=0; j < int(most_freq_str.size());j++)
        {
            if( input_vec[i] == most_freq_str[j])
            {
                frequencies[j] += 1;
            }
        }
    }
}


void remove_punctuation_and_split(string str,vector<string>& str_split)
{
    str = string_replace(str,";",".");
    str = string_replace(str,"?",".");
    str = string_replace(str,"!",".");
    str = string_replace(str,"'",".");
    str = string_replace(str,",",".");

    str_split = split(str, ".");
}


void load_text(ifstream& Text, vector<string>& text_vec)
{
    text_vec.resize(0);
    string dummy;
    while(!Text.eof())
    {
        Text >> dummy;
        text_vec.push_back(dummy);
    }
}


void load_text(string filename, vector<string>& text_vec)
{
    ifstream TextStream(filename);
    load_text(TextStream,text_vec);
}


void clean_and_load_text(ifstream& Text, vector<string>& text_vec)
{
    text_vec.resize(0);

    //Load the text
    string dummy;
    Text >> dummy;
    std::transform(dummy.begin(), dummy.end(), dummy.begin(), ::tolower);

    vector<string> dummy_vec;
    remove_punctuation_and_split(dummy,dummy_vec);
    for(int i = 0; i < int(dummy_vec.size());i++)
    {
        text_vec.push_back(dummy_vec[i]);
    }

    while(!Text.eof())
    {
        Text >> dummy;
        std::transform(dummy.begin(), dummy.end(), dummy.begin(), ::tolower);

        remove_punctuation_and_split(dummy,dummy_vec);

        for(int i = 0; i < int(dummy_vec.size());i++)
        {
            text_vec.push_back(dummy_vec[i]);
        }
    }
}


void write_text(ofstream& TextFile, vector<string> text_vec)
{
    for(int i=0; i < int(text_vec.size());i++)
    {
        TextFile << text_vec[i] << endl;
    }
}


void write_text(string filename, vector<string> text_vec)
{
    ofstream TextStream(filename);
    write_text(TextStream,text_vec);
}


void fileparts(string filename, string& path, string& name, string& ext)
{
    int idx0 = filename.rfind("/");
    int idx1 = filename.rfind(".");

    path = filename.substr(0,idx0+1);
    name = filename.substr(idx0+1,idx1-idx0-1);
    ext  = filename.substr(idx1+1);
}


string replace_extension(string filename, string new_ext)
{
    string path,name,ext;
    fileparts(filename,path,name,ext);
    if(new_ext.length() != 0 && new_ext.at(0) == '.')
    {
        return path+name+new_ext;
    }
    else
    {
        return path+name+"."+new_ext;
    }
}


string remove_extension(string filename)
{
    string path,name,ext;
    fileparts(filename,path,name,ext);
    return path+name;
}


string extract_path(string filename)
{
    string path,name,ext;
    fileparts(filename,path,name,ext);
    return path;
}


string extract_name(string filename)
{
    string path,name,ext;
    fileparts(filename,path,name,ext);
    return name;
}


string extract_ext(string filename)
{
    string path,name,ext;
    fileparts(filename,path,name,ext);
    return ext;
}


void dump_extension_path_name(string filename)
{
    cout << filename << endl;
    cout << "Extension : " <<  extract_ext(filename) << endl;
    cout << "Path : " <<       extract_path(filename) << endl;
    cout << "Name : " <<       extract_name(filename) << endl;
}


// *******************************
// Parsing Tools
// *******************************
void Parser::add_option_stand_alone(string option_name, vector<string> delim_vec)
{
    Option Opt;
    Opt.option_name = option_name;
    Opt.delimitators = delim_vec;
    Opt.NbParams = 0;

    Opt.type = stand_alone; //This is an option without parameters
    options.push_back(Opt);
}


void Parser::add_option_stand_alone(string option_name, string delim)
{
    Option Opt;
    Opt.option_name = option_name;
    Opt.delimitators.resize(1);
    Opt.delimitators[0] = delim;
    Opt.NbParams = 0;

    Opt.type = stand_alone; //This is an option without parameters
    options.push_back(Opt);
    nb_stand_alone_options++;

}


void Parser::add_option_string(string option_name, vector<string> delim_vec, int NbParams)
{
    Option Opt;
    Opt.option_name = option_name;
    Opt.delimitators = delim_vec;
    Opt.type = string_params;
    Opt.NbParams = NbParams;
    Opt.StringValues.resize(NbParams);
    options.push_back(Opt);
    nb_string_options++;
}


void Parser::add_option_double(string option_name, vector<string> delim_vec, int NbParams)
{
    Option Opt;
    Opt.option_name = option_name;
    Opt.delimitators = delim_vec;
    Opt.type = double_params;
    Opt.NbParams = NbParams;
    Opt.DoubleValues.resize(NbParams);
    options.push_back(Opt);
    nb_int_options++;
}


void Parser::add_option_int(string option_name, vector<string> delim_vec, int NbParams)
{
    Option Opt;
    Opt.option_name = option_name;
    Opt.delimitators = delim_vec;
    Opt.type = int_params;
    Opt.NbParams = NbParams;
    Opt.IntValues.resize(NbParams);
    options.push_back(Opt);
}



void Parser::add_option_string(string option_name, string delim, int NbParams)
{
    Option Opt;
    Opt.option_name = option_name;
    Opt.delimitators.resize(1);
    Opt.delimitators[0] = delim;
    Opt.type = string_params;
    Opt.NbParams = NbParams;
    Opt.StringValues.resize(NbParams);
    options.push_back(Opt);
    nb_string_options++;
}


void Parser::add_option_int(string option_name, string delim, int NbParams)
{
    Option Opt;
    Opt.option_name = option_name;
    Opt.delimitators.resize(1);
    Opt.delimitators[0] = delim;
    Opt.type = int_params;
    Opt.NbParams = NbParams;
    Opt.IntValues.resize(NbParams);
    options.push_back(Opt);
}


void Parser::add_option_double(string option_name, string delim, int NbParams)
{
    Option Opt;
    Opt.option_name = option_name;
    Opt.delimitators.resize(1);
    Opt.delimitators[0] = delim;
    Opt.type = double_params;
    Opt.NbParams = NbParams;
    Opt.DoubleValues.resize(NbParams);

    options.push_back(Opt);
    nb_double_options++;
}


void Parser::parse_command_line(int argc, char*const* argv)
{
    for(int optind = 1;optind<argc;optind++)
    {
        //cout << "arg" << optind << " : " << argv[optind] << endl;

        //Loop on the different options
        for(int i = 0; i < int(options.size());i++)
        {
            for(int j=0; j < int(options[i].delimitators.size());j++)
            {
                if(options[i].delimitators[j] == argv[optind])
                {
                    //cout <<  "options[i].delimitators[j] recognized : " << options[i].delimitators[j] << endl;
                    options[i].opt_detected = true;
                    int NbParams = options[i].NbParams;
                    if(options[i].type == int_params)
                    {
                        options[i].IntValues.resize(NbParams);
                        for(int k = 0;k<NbParams;k++)
                        {
                            if(optind+k+1 < argc)
                            {
                                options[i].IntValues[k] = str2int(argv[optind+k+1]);
                                //cout << "options[i].IntValues[k] : " << options[i].IntValues[k] << endl;
                            }
                        }
                    }
                    else if(options[i].type == double_params)
                    {
                        //cout << "options[i].type : double_params. NbParams : " << NbParams << endl;
                        options[i].DoubleValues.resize(NbParams);
                        for(int k = 0;k<NbParams;k++)
                        {
                            //cout << "k  : " << k << endl;
                            if(optind+k+1 < argc)
                            {
                                //cout << "optind+k+1  : " << optind+k+1 << endl;
                                options[i].DoubleValues[k] = string2double(argv[optind+k+1]);
                                //cout << "options[i].DoubleValues[k] : " << options[i].DoubleValues[k] << endl;
                            }
                        }
                    }
                    else if(options[i].type == string_params)
                    {
                        options[i].StringValues.resize(NbParams);
                        for(int k = 0;k<NbParams;k++)
                        {
                            if(optind+k+1 < argc)
                                options[i].StringValues[k].assign(argv[optind+k+1]);
                        }
                    }
                }
            }
        }
    }
}



bool Parser::check_option(string option_name)
{
    for(int i = 0; i < int(options.size());i++)
    {
        if(options[i].option_name == option_name)
        {
            return options[i].opt_detected;
        }
    }
    
    return false;
}


bool Parser::get_string_params(string option_name,vector<string>& option_vec)
{
    for(int i = 0;i< int((*this).options.size());i++)
    {
        if(options[i].option_name == option_name)
        {
            if(options[i].type == string_params &&  options[i].opt_detected == true)
            {
                option_vec = options[i].StringValues;
            }
            return options[i].opt_detected;
        }
    }
    
    return false;
}


bool Parser::get_int_params(string option_name, vector<int>& option_vec)
{
    for(int i = 0; i < int((*this).options.size());i++)
    {
        if(options[i].option_name == option_name)
        {
            if(options[i].type == int_params &&  options[i].opt_detected == true)
            {
                option_vec = options[i].IntValues;
            }
            return options[i].opt_detected;
        }
    }
    return false;
}


bool Parser::get_double_params(string option_name,vector<double>& option_vec)
{
    for(int i = 0; i < int((*this).options.size());i++)
    {
        if(options[i].option_name == option_name)
        {
            if(options[i].type == double_params &&  options[i].opt_detected == true)
            {
                option_vec = options[i].DoubleValues;
            }
            return options[i].opt_detected;
        }
    }
    return false;
}


bool Parser::get_string_param(string option_name,string& option_value)
{
    for(int i = 0; i < int((*this).options.size());i++)
    {
        if(options[i].option_name == option_name)
        {
            if(options[i].type == string_params &&  options[i].opt_detected == true)
            {
                option_value = options[i].StringValues[0];
            }
            return options[i].opt_detected;
        }
    }
    return false;
}


bool Parser::get_int_param(string option_name, int& option_value)
{
    for(int i = 0; i < int((*this).options.size());i++)
    {
        if(options[i].option_name == option_name)
        {
            if(options[i].type == int_params &&  options[i].opt_detected == true)
            {
                option_value = options[i].IntValues[0];
            }
            return options[i].opt_detected;
        }
    }
    return false;
}


bool Parser::get_double_param(string option_name, double& option_value)
{
    for(int i = 0; i < int((*this).options.size());i++)
    {
        if(options[i].option_name == option_name)
        {
            if(options[i].type == double_params && options[i].opt_detected == true)
            {
                option_value = options[i].DoubleValues[0];
            }
            return options[i].opt_detected;
        }
    }
    return false;
}


void Parser::dump_options_and_params()
{
    cout <<  endl;
    cout << "******************************" << endl;
    cout << "*** Options and Parameters ***"<< endl;
    cout << "******************************" << endl;
    cout <<  endl;

    if(nb_stand_alone_options > 0)
        cout << "Stand Alone Options" << endl;

    for(int i = 0; i < int(options.size());i++)
    {
        if(options[i].type == stand_alone)
        {
            cout << "    " << options[i].option_name << " [";

            for(int j = 0; j < int(options[i].delimitators.size());j++)
            {
                cout << options[i].delimitators[j] << " ";
            }

            cout << "] : ";


            if(options[i].opt_detected == true )
            {
                cout << " recognized" << endl;
            }
            else
            {
                cout << " not recognized" << endl;
            }
        }
    }
    cout << endl;

    if(nb_int_options > 0)
        cout << "Integer Parameters" << endl;

    for(int i = 0; i < int(options.size());i++)
    {
        if(options[i].type == int_params)
        {
            cout << options[i].option_name << " [";

            for(int j = 0; j < int(options[i].delimitators.size());j++)
            {
                cout << options[i].delimitators[j] << " ";
            }

            cout << "] : ";

            for(int p = 0;p < int(options[i].IntValues.size());p++)
                if(options[i].opt_detected == true )
                {
                    cout << options[i].IntValues[p] << " ";
                }
                else
                {
                    cout << "not recognized";
                }

            cout <<  endl;
        }
    }
    cout << endl;

    if(nb_double_options > 0)
        cout << "Double Parameters" << endl;

    for(int i = 0; i < int(options.size());i++)
    {
        if(options[i].type == double_params)
        {
            cout << options[i].option_name << " [";

            for(int j = 0; j < int(options[i].delimitators.size());j++)
            {
                cout << options[i].delimitators[j] << " ";
            }

            cout << "] : ";
            if(options[i].opt_detected == true )
            {
                for(int p = 0;p < int(options[i].DoubleValues.size());p++)
                    cout << options[i].DoubleValues[p] << " ";
            }
            else
            {
                cout << " not assigned";
            }
            cout <<  endl;
        }
    }

    cout << endl;

    if(nb_string_options > 0)
        cout << "String Parameters" << endl;

    for(int i = 0; i < int(options.size());i++)
    {
        if(options[i].type == string_params)
        {
            cout << options[i].option_name << " [";
            for(int j = 0; j < int(options[i].delimitators.size());j++)
            {
                cout << options[i].delimitators[j] << " ";
            }

            cout << "] : ";
            if(options[i].opt_detected == true )
            {
                for(int p = 0;p < int(options[i].StringValues.size());p++)
                    cout << options[i].StringValues[p];
            }
            else
            {
                cout << "not assigned";
            }

            cout <<  endl;
        }
    }
    cout << endl;
}


//**********************************************
// * Numerical routines
//**********************************************
double truncate(double x, int nb_dec)
{
    double z = x;
    z *= pow(10,nb_dec);
    z = trunc(z);
    z /= pow(10,nb_dec);

    return z;
}

double myclip(double x, double min, double max)
{
    double y = x;
    if(min<max)
    {
        if(y<min)
            y =min;

        if(y>max)
            y = max;
    }
    return y;
}


double sign(double x)
{
    double y;
    if(x>0)
        y=1;
    else if(x<0)
        y=-1;
    else
        y = 0;
    return y;
}

double min3(double a, double b, double c)
{  return std::min(a, std::min(b, c)); }
double min4(double a, double b, double c, double d)
{  return std::min(a, min3(b, c, d)); }

double max3(double a, double b, double c)
{  return std::max(a, std::max(b, c)); }
double max4(double a, double b, double c, double d)
{  return std::max(a, max3(b, c, d)); }

double Mean(vector<double> vec)
{
    double mean_val = 0;

    if(vec.size() > 0)
    {
        for(int i = 0; i < int(vec.size());i++)
            mean_val += vec[i];

        mean_val /= vec.size();
    }

    return mean_val;
}


double Var(vector<double> vec)
{
    double var_val = 0;

    if(vec.size() > 0)
    {
        double mean_val = Mean(vec);

        for(int i = 0;i< int(vec.size());i++)
            var_val += (vec[i] - mean_val)*(vec[i] - mean_val);

        var_val /= vec.size();
    }

    return var_val;
}


double Norm(vector<double> vec, double p)
{
    double norm_val = 0.;
    if(p==2)//Default: this is the L2-norm
    {
        for(int i = 0; i < int(vec.size());i++)
            norm_val  += vec[i]*vec[i];

        norm_val = sqrt(norm_val);
    }
    else if(p==1)
    {
        for(int i = 0; i < int(vec.size());i++)
            norm_val += fabs(vec[i]);
    }
    else if(p > 0.0)
    {
        for(int i = 0; i < int(vec.size());i++)
            norm_val += pow(fabs(vec[i]),p);

        norm_val = pow(norm_val,1.0/p);
    }
    return norm_val;
}


double findMin(vector<double> vec,int& MinIndx)
{
    MinIndx = -1;
    double minval = std::numeric_limits<double>::max();
    for(int i = 0; i < int(vec.size());i++)
    {
        if(vec[i] < minval)
        {
            minval = vec[i];
            MinIndx = i;
        }
    }
    return minval;
}


double findMax(vector<double> vec,int& MaxIndx)
{
    MaxIndx = -1;
    double maxval = std::numeric_limits<double>::min();
    for(int i = 0;i < int(vec.size());i++)
    {
        if(vec[i] > maxval)
        {
            maxval = vec[i];
            MaxIndx = i;
        }
    }
    return maxval;
}


void WriteTXT( vector<double>& Values, string& filenameOut)
{
    ofstream out_str(filenameOut);
    out_str << Values.size()  << endl;
    out_str << 1 << endl;

    for(unsigned int i =0;i<Values.size();i++)
        out_str << Values[i] << endl;
}


void skip_comment_lines(ifstream& In, char comment_symbol)
{
    char Buffer;
    string line_str;
    if(In.tellg() > 0)
    {
        In.seekg(-1, std::ios::cur);
        getline(In,line_str);
    }

    Buffer = In.peek();
    while(Buffer == comment_symbol && !In.eof() ) //Skip comment lines
    {
        getline(In,line_str);
        Buffer = In.peek();
    }
}


//**************************************************
// Uniform random variable generation
//**************************************************
double rand_unif(int seed)
{
    if(seed != -1)
        srand(seed);

    double rand1 = rand() / ((double) RAND_MAX);
    return rand1;
}


//**************************************************************
// Integer from 1 to Max (default : simulated ouptut of a dice)
//**************************************************************
int rand_int(int Max, int seed)
{
    double rand1 = rand_unif(seed);
    return ceil(rand1*Max);
}


//**************************************************
// Gaussian random variable generation:
// Simple implementation of the Box Mueller method
//**************************************************
double randn(int seed)
{
    if(seed != -1)
        srand(seed);

    double rand1, rand2;
    rand1 = rand() / ((double) RAND_MAX);
    if(rand1 < 1e-100) rand1 = 1e-100;

    rand1 = -2 * log(rand1);
    rand2 = (rand() / ((double) RAND_MAX)) * TWO_PI;

    return sqrt(rand1) * cos(rand2);
}


vector<double> randn(int size, int seed)
{
    vector<double> rand_vec;
    rand_vec.push_back(randn(seed));

    for(int i=1;i<size;i++)
    {
        rand_vec.push_back(randn());
    }

    return rand_vec;
}


int my_cmp(const void *x, const void *y)
{
    double xx = *(double*)x, yy = *(double*)y;
    if (xx < yy) return -1;
    if (xx > yy) return  1;
    return 0;
}

int my_cmp_descent(const void *x, const void *y)
{
    double xx = *(double*)x, yy = *(double*)y;
    if (xx < yy) return 1;
    if (xx > yy) return  -1;
    return 0;
}


double vec_mean(double* x, int Length)
{
    //Warning: no size control !
    double mean = 0.;
    for(int i =0;i<Length;i++)
        mean += x[i];
    if(Length > 0)
        mean /= Length;

    return mean;
}


//Warning: the OutStr should be initilazed before using this function
//This follows the bytewise bigendian convention to cope with fread,fwrite
void writeIntTo4Bytes(ofstream& OutStr, int x)
{
    int x1 = ((x >> 24) & 255);
    int x2 = ((x >> 16) & 255);
    int x3 = ((x >> 8) & 255);
    int x4 = ((x)& 255);

    char c1 = (char)x1;
    char c2 = (char)x2;
    char c3 = (char)x3;
    char c4 = (char)x4;

    OutStr << c4;
    OutStr << c3;
    OutStr << c2;
    OutStr << c1;

    OutStr.flush();
}


//This follows the bytewise bigendian convention to cope with fread,fwrite
int readIntFrom4Bytes(istream& InStr)
{
    int x = 0;
    char c1, c2, c3, c4;
    InStr >> c4;
    InStr >> c3;
    InStr >> c2;
    InStr >> c1;

    x = ((c1 & 255) << 24) + ((c2 & 255) << 16) + ((c3 & 255) << 8) + (c4 & 255);
    return x;
}



int fncompare(const void * el1,const void * el2)
{
    double **d1 = (double **)el1;
    double **d2 = (double **)el2;

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


bool is_power_of_two(int n)
{
    return n > 0 && ((n & (n-1)) == 0);
}


void printExecutionTime(ClockTime start_time, ClockTime end_time)
{
    auto execution_time_ns = duration_cast<nanoseconds>(end_time - start_time).count();
    auto execution_time_ms = duration_cast<microseconds>(end_time - start_time).count();
    auto execution_time_sec = duration_cast<seconds>(end_time - start_time).count();
    auto execution_time_min = duration_cast<minutes>(end_time - start_time).count();
    auto execution_time_hour = duration_cast<hours>(end_time - start_time).count();

    cout << "\nExecution Time: ";
    if(execution_time_hour > 0)
    cout << "" << execution_time_hour << " Hours, ";
    if(execution_time_min > 0)
    cout << "" << execution_time_min % 60 << " Minutes, ";
    if(execution_time_sec > 0)
    cout << "" << execution_time_sec % 60 << " Seconds, ";
    if(execution_time_ms > 0)
    cout << "" << execution_time_ms % long(1E+3) << " MicroSeconds, ";
    if(execution_time_ns > 0)
    cout << "" << execution_time_ns % long(1E+6) << " NanoSeconds, ";
}


//Date and hour stamp
string date_stamp()
{
  // Using time point and system_clock
  std::time_t save_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

  struct tm* timeinfo;
  timeinfo = localtime(&save_time);

  string date_stamp_str;
  date_stamp_str += int2str(timeinfo->tm_year + 1900);
  date_stamp_str += "_";
  date_stamp_str += int2str(timeinfo->tm_mon +1);
  date_stamp_str += "_";
  date_stamp_str += int2str(timeinfo->tm_mday);

  return date_stamp_str;
}


//Date and hour stamp
string date_hour_stamp()
{
  // Using time point and system_clock
  std::time_t save_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  struct tm* timeinfo;
  timeinfo = localtime(&save_time);

  string date_stamp_str;
  date_stamp_str += int2str(timeinfo->tm_year + 1900);
  date_stamp_str += "_";
  date_stamp_str += int2str(timeinfo->tm_mon +1);
  date_stamp_str += "_";
  date_stamp_str += int2str(timeinfo->tm_mday);
  date_stamp_str += "_";

  if(timeinfo->tm_hour < 10)
      date_stamp_str += "0";

  date_stamp_str += int2str(timeinfo->tm_hour);
  date_stamp_str += ":";
  if(timeinfo->tm_min < 10)
      date_stamp_str += "0";

  date_stamp_str += int2str(timeinfo->tm_min);
  date_stamp_str += ":";

  if(timeinfo->tm_sec < 10)
      date_stamp_str += "0";

  date_stamp_str += int2str(timeinfo->tm_sec);

  return date_stamp_str;
}



//*******************************
//*******************************
// * Tests
//*******************************
//*******************************
void test_fileparts()
{
    cout << "test fileparts" << endl;
    string filename("Desktop/3D/svn_TLB/Photo/tavira.png.pgm");
    string path, name, ext;
    fileparts(filename,path, name, ext);

    cout << "decomposition : " << path << " , " << name << " , " << ext << endl;
}

void test_string()
{
    string test_string = "abc";
    int x = 5;
    double y = 12.376589;
    cout << "preappended string: " << test_string << endl;

    test_string = append(test_string,x)+"efg";
    cout << "appended string: " << test_string << endl;
    test_string = append(test_string,(int)y)+"hij";
    cout << "double appended: " << test_string << endl;
    //test_string = apppend(test_string,);

    double f = 3.14159;
    cout << std::setprecision(5) << f << '\n';
    cout << std::setprecision(9) << f << '\n';
    cout << std::fixed;
    cout << std::setprecision(5) << f << '\n';
    cout << std::setprecision(9) << f << '\n';
}



void test_formats()
{
    cout << "sizeof(char) : " << sizeof(char) << endl;
    cout << "sizeof(unsigned char) : " << sizeof(unsigned char) << endl;
    cout << "sizeof(short) : " << sizeof(short) << endl;
    cout << "sizeof(unsigned short) : " << sizeof(unsigned short) << endl;
    cout << "sizeof(int) : " << sizeof(int) << endl;
    cout << "sizeof(unsigned int) : " << sizeof(unsigned int) << endl;
    cout << "sizeof(float) : " << sizeof(float) << endl;
    cout << "sizeof(double) : " << sizeof(double) << endl;
    cout << endl;

    double a = 1.33;
    char b = (char)20;
    char c = (char)(a*b);
    double d = a*b;
    cout << std::setprecision(7);
    cout << "a : " << a << endl;
    cout << "(double)b : " << (double)b << endl;
    cout << "(double)c : " << (double)c << endl;
    cout <<  "d : " << std::setprecision(7)  << d << endl;

    cout << "1<<3 : " << (1<<3) << endl;
    cout << "1>>3 : " << (1>>3) << endl;

    /* int z= 0x3F;
  cout << "z : " << z << endl;
  int xi=5,yi=5;
  int z1 = xi >> 1;
  cout << "z1 : " << z1;
  int z2 = yi << 1;
  cout << "z2 : " << z2; */
}


void test_precision()
{
    double z = 0.73456712345678912345354682479;
    cout << "z : " << z << endl;
    cout << setprecision(16) << endl;
    cout << "z (prec. 16) : " << z << endl;
    int i = 11;
    z = truncate(z,i);
    /*z *= pow(10,i);
  z = trunc(z);
  z /= 10e14; */
    cout << "truncated z : " << z << endl;
    cout << setprecision(16) << "truncated z : " << z << endl;
}


void test_map()
{
    //************ test map
    cout << "test_map" << endl;
    map<int,double> x;
    x[3]=3.45;
    x[7]=5.45;
    cout << "x[7] : " << x[7] << endl;

    map<char,int> first;

    first['a'] = 12;
    first['b'] = 30;
    first['c'] = 50;
    first['d'] = 70;

    cout << "first['e'] : " << first['e'] << endl;
    auto search = first.find('c');

    cout <<  "search->first :  " <<  search->first << endl;;
    cout <<  "search->second :  " <<  search->second << endl;

    //************ test tuple
    tuple<int,double,double> ixy;
    ixy = std::make_tuple(2,3.1,4.0);
    cout << std::get<0>(ixy) << " " << std::get<1>(ixy) << endl;

    vector<tuple<int,double,double>> gt_indices_features;
    gt_indices_features.push_back(make_tuple(1,2.3,3.4));
    gt_indices_features.push_back(make_tuple(3,14.3,20.4));
    gt_indices_features.push_back(make_tuple(35,37.8,43.1));

    for(int i = 0;i<(int)gt_indices_features.size();i++)
    {
        int indx = get<0>(gt_indices_features[i]);
        double px = get<1>(gt_indices_features[i]);
        double py = get<2>(gt_indices_features[i]);
        cout << "indx : " << indx << " , px : " << px << " , py : " << py << endl;
    }
}


//Minimal function for the functionalities of unique_ptr
void test_write_char_ptr()
{
    /* char *p = "hello";
    char q[] = "hello"; // no need to count this

    printf("%zu\n", sizeof(p)); // => size of pointer to char -- 4 on x86, 8 on x86-64
    printf("%zu\n", sizeof(q)); // => size of char array in memory -- 6 on both

    // size_t strlen(const char *s) and we don't get any warnings here:
    printf("%zu\n", strlen(p)); // => 5
    printf("%zu\n", strlen(q)); // => 5

    string filename("/home/ash/Desktop/test.txt");
    int fd = open(filename.c_str(), O_CREAT | O_WRONLY | O_SYNC, 0777);

    //First with a string
    string buffer_str("konnichiwa\n");
    int sz = write(fd, buffer_str.c_str(), strlen(buffer_str.c_str()));

    //Second with a pointer to char
    char *r = "arigato\n";
    sz = write(fd, r, strlen(r));

    //Third with an array of chars
    char s[] = "gozaimasu\n";
    sz = write(fd, s, strlen(s));
    cout << "sz : " << sz << endl;

    //Fourth with a unique_ptr of char[]
    int buffer_size = 7;
    std::unique_ptr<char[]> char_Ptr(new char[buffer_size]);
    char_Ptr[0]='a';
    char_Ptr[1]='b';
    char_Ptr[2]='c';
    char_Ptr[3]='d';
    char_Ptr[4]='e';
    char_Ptr[5]='f';
    char_Ptr[6]='\n';

    sz = write(fd,char_Ptr.get(),strlen(char_Ptr.get()));
    cout << "strlen(char_Ptr.get()) : "   << strlen(char_Ptr.get()) << endl;

    //Fourth with a unique_ptr of char
    buffer_size = 9;
    std::unique_ptr<char> char_Ptr2(new char[buffer_size]);
    char_Ptr2.get()[0]='z';
    char_Ptr2.get()[1]='y';
    char_Ptr2.get()[2]='x';
    char_Ptr2.get()[3]='w';
    char_Ptr2.get()[4]='v';
    char_Ptr2.get()[5]='u';
    char_Ptr2.get()[6]='t';
    char_Ptr2.get()[7]='s';
    char_Ptr2.get()[8]='\n';

    sz = write(fd,char_Ptr2.get(),(size_t)buffer_size);

    close(fd);
    string cmd_line("/usr/bin/gedit "+filename);
    system(cmd_line.c_str()); */
}



void test_slice()
{
    valarray<double> x;
    x = {1.2,2.4,3.6,4.8, 6.2,7.55,8.0,9.0};
    int N = x.size();
    for(int i = 0;i < N;i++)
    {
        cout << "x["<<i<<"] : "<< x[i] << endl;
    }
    valarray<double> y=x[slice(0,N/2,2)];
    for(int i = 0;i < (int)y.size();i++)
    {
        cout << "y["<<i<<"] : "<< y[i] << endl;
    }


    slice s = slice(0,N/2,2);
    cout << "slice(0,N/2,2)" << endl;
    for(int i = 0;i<(int)s.size();i++)
    {
       cout <<  s.start()+i*s.stride() << " , ";
    }
    cout << endl;

    vector<double> x_vec;
    x_vec = {1.2,2.4,3.6,4.8, 6.2,7.55,8.0,9.0};

    vector<double> y_vec, y_vec1;
    y_vec = subvector(x_vec,s);
    vector<int> indices;
    indices = {0,2,4,6};
    y_vec1 = subvector(x_vec,indices);

    dump_vec(y_vec,"y_vec");
    dump_vec(y_vec1,"y_vec1");

    vector<complex<double>> z;
    z = {1,2,3,4,5,6,7,8};
    vector<complex<double>> zz;
    zz = subvector(z,s);

    cout << "z" << endl;
    for(int i =0;i<(int)z.size();i++)
    {
        cout << "z["<<i<<"] : " << z[i]<< endl;
    }

    cout << "zz" << endl;
    for(int i =0;i<(int)zz.size();i++)
    {
        cout << "zz["<<i<<"] : " << zz[i]<< endl;
    }

}


//****************
/* void test_std_complex()
{
    using namespace std::complex_literals;
    complex<double> x = 0.8 + 0.8i;
    complex<double> y = 1i;

    //complex<double> z = x*y;
    cout << "x : " << x << endl;
    cout << "y : " << y << endl;
    //cout << "z : " << z << endl;
    //cout << "z:" << z.real() << " + " << z.imag() << " * i " << endl;
    //cout << "conj(z):" << conj(z).real() << " + " << conj(z).imag() << " * i " << endl;
    //dump_complex(z);
    vector<complex<double>> vec;
    vec.push_back(0.0007);
    cout << "lab1 : " << vec.size() << endl;
    vec.push_back(complex<double>(0.3,0.2));
    cout << "lab2 : " << vec.size() << endl;
    vec.push_back(x);
    cout << "lab3 : " << vec.size() << endl;

    cout << "vec[1] : " << vec[1] << endl;
    cout << "vec[2] : " << vec[2] << endl;
    dump_vec(vec);
}*/


void get_start_time(struct timespec& Start)
{
    clock_gettime(CLOCK_MONOTONIC, &Start);
}


// Warning : This works in LINUX
// inspired from arrilib/HighPerfTimer.h
double get_elapsed_time(struct timespec Start)
{
    struct timespec Stop;
    clock_gettime(CLOCK_MONOTONIC, &Stop);
    double duration_mu = ((Stop.tv_sec * 1000000) + (Stop.tv_nsec / 1000)) - ((Start.tv_sec * 1000000) + (Start.tv_nsec / 1000));
    double duration_sec = duration_mu / 1000000.;
    //cout << "Elapse time from last Start " << duration_mu / 1000000. << " sec"<< endl;
    return duration_sec;
}


void test_deque()
{
  std::deque<int> myints;
  std::cout << "0. size: " << myints.size() << '\n';

  for (int i=0; i<5; i++) myints.push_back(i);
  std::cout << "1. size: " << myints.size() << '\n';

  myints.insert (myints.begin(),5,100);
  std::cout << "2. size: " << myints.size() << '\n';

  myints.pop_back();
  std::cout << "3. size: " << myints.size() << '\n';
}

void test_time()
{
    // Get starting timepoint
    auto start = high_resolution_clock::now();

    //write these two lines at start
    struct timespec Start; //, Stop;
    get_start_time(Start);

    // Call the function, here sort()
    //sort(values.begin(), values.end());
    double x = 0.;
    int NbIt = 1000*1000*1000;
    for(int i = 0;i<NbIt;i++)
    {
        x += 0.01;
    }

    get_elapsed_time(Start);

    // Get ending timepoint
    auto stop = high_resolution_clock::now();

    // Get duration. Substart timepoints to get durarion.
    // To cast it to proper unit use duration cast method
    auto duration = duration_cast<microseconds>(stop - start);

    cout << "Time taken by function: " << duration.count() << " microseconds" << endl;
}



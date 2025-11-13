#ifndef TEST_FUNCTION_H
#define TEST_FUNCTION_H
#include <iostream>

int my_test_add(int a, int b);
double exp2_elliptic(double x,double y,double alpha,double beta);

class Identite
{
      public:                      // begin public section
        Identite(std::string nom,std::string prenom);       // constructor
        ~Identite();                    // destructor

        void Dump();

     private:                      // begin private section
        std::string m_Nom;                // member variable
        std::string m_Prenom;                // member variable
};


// constructor of Cat,
Identite::Identite(std::string nom,std::string prenom)
{
  m_Nom = nom;
  m_Prenom = prenom;
}

// destructor, just an example
Identite::~Identite()
{
}

void Identite::Dump()
{
    std::cout << "Prenom, Nom : " << m_Prenom << " , " << m_Nom << std::endl;
}
#endif // TEST_FUNCTION_H

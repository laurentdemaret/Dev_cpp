#include <iostream>
#include <vector>
#include <map>

using namespace std;


#include <test_function.h>
#include <Point2D.h>

//A partir du
#include "matrix.h"


void test_algo_tools()
{
    std::cout << "fonctionne sous LINUX  et sous MAC 13 nov." << std::endl;
    std::cout << "13:28" << std::endl;

    M3Matrix A;
    A.Reshape(2,2,0.);

    A[0][1]=1.0;
    A[1][1]=2.0;
    A[1][0]=3.0;
    A[1][1]=4.0;
    std::cout << "test_algo_tools " << std::endl;
    std::cout << "Matrix A" << std::endl;
    A.Dump();
}

void test_map_2()
{
    std::map<int, Point2D> pointMap;

    // Create a std::map with key as int and value as string
    std::map<int, std::string> studentMap;

    // Inserting elements into the map
    studentMap.insert(std::make_pair(1, "John"));
    studentMap.insert(std::make_pair(2, "Sarah"));
    studentMap.insert(std::make_pair(4, "Michael"));
    studentMap.insert(std::make_pair(4, "Julie"));
    // Accessing elements using keys
    std::cout << "Student with key 4: " << studentMap[4] << std::endl;

    // Updating the value associated with a key
    studentMap[1] = "Alice";

    // Checking if a key exists
    if (studentMap.find(3) != studentMap.end()) {
        std::cout << "Key 3 exists!" << std::endl;
    }
    else
    {
        std::cout << "Key 3 does not exist!" << std::endl;
    }


    if (studentMap.find(4) != studentMap.end()) {
        std::cout << "Key 4 exists!" << std::endl;
    }
    else
    {
        std::cout << "Key 4 does not exist!" << std::endl;
    }

    // Iterating over the map
    std::cout << "All students: ";
    /*for (const auto& student : studentMap) {
        std::cout << student.second << " ";
    }
    std::cout << std::endl;
*/

    if (studentMap.find(2) != studentMap.end()) {
        std::cout << "Avant  effacageKey 2 existe!" << std::endl;
    }
    else
    {
        std::cout << "Avant effacage Key 2 n'existe pas!" << std::endl;
    }

    // Removing an element
    studentMap.erase(2);

    if (studentMap.find(2) != studentMap.end()) {
        std::cout << "Après effacage Key 2 existe!" << std::endl;
    }
    else
    {
        std::cout << "Après effacage Key 2 n'existe pas" << std::endl;
    }

    // Checking the size of the map
    std::cout << "Number of students: " << studentMap.size() << std::endl;
}

void bonjour()
{
    cout << "Bonjour -- 13 novembre 2025 !" << endl;
    int x=4;
    int y=5;
    std::cout << my_test_add(x,y) << std::endl;

    vector<double> a;
    a.push_back(1.0);
    a.push_back(2.0);
    a.push_back(3.0);

    cout << "a.size(): " << a.size() <<  endl;

    string prenom("Jean");
    string nom("Dupont");

    Identite personne1(nom,prenom);

    personne1.Dump();
}


void test_point2D()
{
    Point2D A(3,.4);
    Point2D B(4,.1);

    std::cout << " A.x: " << A.xm  << std::endl;
    std::cout << " A.y: " << A.ym  << std::endl;
    std::cout << std::endl;

    Point2D C(B);
    C.Move(10,10);
    C.Display();
    C = A;
    std::cout << "après le réassignement "<< std::endl;
    std::cout << "bla"<< std::endl;

    C.Display();
}


int main()
{
    test_algo_tools();
   //bonjour();
   //test_point2D();

   // ************************
   //test_map();

   return 0;
}

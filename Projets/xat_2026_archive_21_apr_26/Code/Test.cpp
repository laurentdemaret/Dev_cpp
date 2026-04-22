#include <vector>
#include <iostream>
#include "Test.h"


using namespace std;

Test::Test()
{
}

Test::~Test()
{
}

void Test::list(vector<Point2D*> l)
{
	for (unsigned int i = 0; i < l.size(); i++) l[i]->debug();
}

void Test::list(std::vector<Edge*> l)
{
	for (unsigned int i = 0; i < l.size(); i++) l[i]->debug();
}

void Test::list(std::vector<Triangle*>)
{
	//for (unsigned int i = 0; i < l.size(); i++) l[i]->debug();
}



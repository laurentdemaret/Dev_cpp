#include <cmath>
#include <vector>
#include <list>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "Point2D.h"
//#include "Edge.h"
//#include "Triangle.h"
#include "VoronoiCell.h"


//class Point2D;

using namespace std;

void VoronoiCell::CellFromCentroid()
{
//TODO: copier de Point2D.h
    vector<Triangle*> DelaunayNeighborTriangles = Centroid->getTriangles();

}

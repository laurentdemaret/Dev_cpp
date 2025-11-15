#ifndef POINTGRID_H_
#define POINTGRID_H_

#include <vector>
#include "Point2D.h"

class PointGrid
{
public:

	PointGrid(std::vector<Point2D*>, unsigned int, unsigned int);
	virtual ~PointGrid();
	
	/* attributes */
	std::vector<Point2D*> nodes;
	unsigned int nbRows;
	unsigned int nbCols;
	unsigned int xMax;
	unsigned int yMax;
	unsigned int xMin;
	unsigned int yMin;
		
};

#endif /*POINTGRID_H_*/

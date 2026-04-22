#include "PointGrid.h"

PointGrid::PointGrid(std::vector<Point2D*> n, unsigned int r, unsigned int c)
{
	this->nodes = n;
	this->nbRows = r;
	this->nbCols = c;
	this->xMin = 0;
	this->yMin = 0;
    //TODO: check this - here x is a column index. Logical ? Invert it ?
    //THese lines:
    /*this->xMax = c-1;
    this->yMax = r-1;*/
    //replaced (22/04/2026) by
    this->xMax = r-1;
    this->yMax = c-1;

	for (unsigned int i = 0; i < this->nodes.size(); i++) {
		this->nodes[i]->gridBound = ( 
										( this->nodes[i]->x == this->xMax || this->nodes[i]->x == this->xMin )
									&&  ( this->nodes[i]->y == this->yMax || this->nodes[i]->y == this->yMin )
									);
		this->nodes[i]->convexHull = ( this->nodes[i]->x == this->xMax
									|| this->nodes[i]->y == this->yMax
									|| this->nodes[i]->x == this->xMin
									|| this->nodes[i]->y == this->yMin );
										
	}
}

PointGrid::~PointGrid()
{	
}

#include <vector>
#include <iostream>
#include <math.h>
#include "Triangle.h"

using namespace std;

vector<Triangle*> Triangle::pool;
int Triangle::new_triangle_cnt = 0;

void Triangle::initTriangle(Point2D* p1, Point2D* p2, Point2D* p3)
{
	// points in counter clockwise
    this->v1 = p1;
    if(v1->ccw(p2,p3))
	{
        this->v2 = p2;
        this->v3 = p3;
	}
	else
	{
        this->v2 = p3;
        this->v3 = p2;
	}
}

Triangle::Triangle(Point2D* p1, Point2D* p2, Point2D* p3)
{
	this->initTriangle(p1,p2,p3);
}

Triangle::Triangle(Edge* edge, Point2D* point3) 
{
	this->initTriangle(edge->org,edge->dest,point3);
}

Triangle::~Triangle()
{
}

void Triangle::deleteTriangles(vector<Triangle*> triangles)
{
	for (unsigned int i = 0; i < triangles.size(); i++) 
	{
		triangles[i]->deleteTriangle();
	}
}

void Triangle::deleteTriangle()
{
	// refactor triangle 
    this->v1 = NULL;
    this->v2 = NULL;
    this->v3 = NULL;

	if(Triangle::pool.size() < TRIANGLE_POOL_SIZE)
	{
		Triangle::pool.push_back(this);
	}
	else
	{
		// free triangle 
		delete this;	
	}	
}

Triangle* Triangle::makeTriangle(Point2D* p1,Point2D* p2,Point2D* p3)
{
	if(Triangle::pool.size() > 0){
		Triangle* t = *(Triangle::pool.begin());
		Triangle::pool.erase(Triangle::pool.begin());
		t->initTriangle(p1,p2,p3);
		return t;
	}
	else{
		new_triangle_cnt++;
		return new Triangle(p1, p2, p3);	
	}	
}

Triangle* Triangle::makeTriangle(Triangle* triangle)
{
    return Triangle::makeTriangle(triangle->v1,triangle->v2,triangle->v3);
}


Triangle* Triangle::makeTriangle(Edge* edge,Point2D* point3)
{
	return Triangle::makeTriangle(edge->org,edge->dest,point3);
}

bool Triangle::isVertex(Point2D* p)
{
    return (this->v1 == p) || (this->v2 == p) || (this->v3 == p);
}


bool Triangle::insideTriangle(Point2D* point)
{
	return this->insideTriangle(point->x,point->y);
}

// int or double (LD: 18/11/25)
// attribute ??
int Triangle::xMin()
{
    int xMin = v1->x;
    if(v2->x < xMin) xMin = v2->x;
    if(v3->x < xMin) xMin = v3->x;

    return xMin;
}

int Triangle::yMin()
{
    int yMin = v1->y;
    if(v2->y < yMin) yMin = v2->y;
    if(v3->y < yMin) yMin = v3->y;

    return yMin;
}

int Triangle::xMax()
{
    int xMax= v1->x;
    if(v2->x > xMax) xMax = v2->x;
    if(v3->x > xMax) xMax = v3->x;

    return xMax;
}

int Triangle::yMax()
{
    int yMax= v1->y;
    if(v2->y > yMax) yMax = v2->y;
    if(v3->y > yMax) yMax = v3->y;

    return yMax;
}

double Triangle::Det()
{
    double det = (v2->x - v1->x)*(v3->y-v1->y) - (v2->y-v1->y)*(v3->x-v1->x);

    return det;
}

// check if thispoint inside or on one edge
// of the triangle builded by the 3 points
// parameters: point1, point2, point3
bool Triangle::insideTriangle(double px, double py)
{
    int p1x = this->v1->x;
    int p1y = this->v1->y;
	
	int p2x;
	int p2y;
	int p3x;
	int p3y;
    if(this->v1->ccw(this->v2,this->v3))
	{
        p2x = this->v2->x;
        p2y = this->v2->y;
        p3x = this->v3->x;
        p3y = this->v3->y;
	}
	else{
        p2x = this->v3->x;
        p2y = this->v3->y;
        p3x = this->v2->x;
        p3y = this->v2->y;
	}

	double ab = (py-p1y)*(p2x-p1x) - 
				(px-p1x)*(p2y-p1y);
	double bc = (py-p2y)*(p3x-p2x) -
				(px-p2x)*(p3y-p2y);
	double ca = (py-p3y)*(p1x-p3x) -
				(px-p3x)*(p1y-p3y);
	
	return (ab>=0. && bc >= 0. && ca >= 0.);
	
//TODO pruefen ob das wirklich on edge ist ??? 
//TODO wenn doch dann bei inline formell benutzen	
//		if (ab*bc>0 && bc*ca>0) {
//			return true;
//		}
//		else {
//			return this.inline(point1,point2) 
//					|| this.inline(point2,point3) 
//					|| this.inline(point1,point3);
//		}
}

void Triangle::getEdges(Edge** e1,Edge** e2,Edge** e3)
{
	int count = 0;
    Edge* entry = this->v1->entry;
    if(entry->isVertex(this->v2))
	{
		*e1 = entry;
		count++;
	}
    if(entry->isVertex(this->v3))
	{
		*e3 = entry;
		count++;
	}		
    Edge* orbit = entry->getNext(this->v1);
	while ( count < 2 && orbit != entry ) 
	{
        if(orbit->isVertex(this->v2))
		{
			*e1 = orbit;
			count++;
		}
        if( orbit->isVertex(this->v3))
		{
			*e3 = orbit;
			count++;
		}
        orbit = orbit->getNext(this->v1);
	}
    entry = this->v2->entry;
    if(entry->isVertex(this->v3))
	{
		*e2 = entry;
	}		
    orbit = entry->getNext(this->v2);
	while ( orbit != entry ) {
        if( orbit->isVertex(this->v3) )
		{
			*e2 = orbit;
			break;
		}
        orbit = orbit->getNext(this->v2);
	}
}


/*std::pair<Point2D*, Point2D*> Triangle::getOtherPoints(Point2D* a)
{
    // Check that this is a triangle
    if (a == point1) return {point2, point3};
    if (a == point2) return {point1, point3};
    if (a == point3) return {point1, point2};

    throw std::invalid_argument("Point not in triangle");
}*/


Point2D* Triangle::getOtherPoint(Point2D* a, Point2D* b)
{
    if(this->v1 != a
    && this->v1 != b )
	{
        return this->v1;
	}
    if( this->v2 != a
    &&  this->v2 != b )
	{
        return this->v2;
	}
	else
	{
        return this->v3;
	}
}


Point2D* Triangle::getOtherPoint(Edge* edge)
{
	return this->getOtherPoint(edge->org, edge->dest);
}


Point2D* Triangle::getVoronoiCentroid(int px, int py)
{
    double d1, d2, d3, d;

    int indx_nn = 1;
    d1 = ((*this).v1->x -px)*((*this).v1->x -px) + ((*this).v1->y -py)*((*this).v1->y -py);
    d2 = ((*this).v2->x -px)*((*this).v2->x -px) + ((*this).v2->y -py)*((*this).v2->y -py);
    d3 = ((*this).v3->x -px)*((*this).v3->x -px) + ((*this).v3->y -py)*((*this).v3->y -py);
    d=d1;
    if(d2<d)
    {
        d=d2;
        indx_nn = 2;
    }

    if(d3<d)
    {
        d=d3;
        indx_nn = 3;
    }

    if(indx_nn == 1)
    {
        return this->v1;
    }
    if(indx_nn == 2)
    {
        return this->v2;
    }
    if(indx_nn == 3)
    {
        return this->v3;
    }

}



Point2D* Triangle::getVoronoiCentroid(Point2D* p)
{
    return getVoronoiCentroid(p->x,p->y);
    /*double d1, d2, d3, d;

    int indx_nn = 1;
    d1 = ((*this).v1->x -p->x)*((*this).v1->x -p->x) + ((*this).v1->y -p->y)*((*this).v1->y -p->y);
    d2 = ((*this).v2->x -p->x)*((*this).v2->x -p->x) + ((*this).v2->y -p->y)*((*this).v2->y -p->y);
    d3 = ((*this).v3->x -p->x)*((*this).v3->x -p->x) + ((*this).v3->y -p->y)*((*this).v3->y -p->y);
    d=d1;
    if(d2<d)
    {
        d=d2;
        indx_nn = 2;
    }

    if(d3<d)
    {
        d=d3;
        indx_nn = 3;
    }

    if(indx_nn == 1)
    {
        return this->v1;
    }
    if(indx_nn == 2)
    {
        return this->v2;
    }
    if(indx_nn == 3)
    {
        return this->v3;
    }*/
}

Triangle* Triangle::getNeighbourByEdge(Edge *edge)
{
	Point2D* other = this->getOtherPoint(edge);
	Point2D* nextOther = edge->onext->getOtherPoint(edge->org);
	if( nextOther == other ){
		Point2D* prevOther = edge->oprev->getOtherPoint(edge->org);		
		if(edge->checkTriangle(prevOther) && prevOther != other )
		{
			return Triangle::makeTriangle(edge, prevOther);
		}
		else
		{
			return NULL;
		}
	}
	else
	{
		if(edge->checkTriangle(nextOther))
		{
			return Triangle::makeTriangle(edge, nextOther);
		}
		else
		{
			return NULL;
		}
	}
}

Point2D* Triangle::getNeighbourPointByEdge(Edge* edge)
{
	Point2D* other = this->getOtherPoint(edge);
	Point2D* nextOther = edge->onext->getOtherPoint(edge->org);
	if( nextOther == other ){
		Point2D* prevOther = edge->oprev->getOtherPoint(edge->org);		
		if( edge->checkTriangle(prevOther) && prevOther != other ){
			return prevOther;
		}
		else
		{
			return NULL;
		}
	}
	else
	{
		if(edge->checkTriangle(nextOther))
		{
			return nextOther;
		}
		else{
			return NULL;
		}
	}
}


double Triangle::GetArea()
{	
    double ax = (double)(this->v1->x);
    double ay = (double)(this->v1->y);
    double bx = (double)(this->v2->x);
    double by = (double)(this->v2->y);
    double cx = (double)(this->v3->x);
    double cy = (double)(this->v3->y);

	double area =  0.5*((ax-cx)*(by-cy)-(ay-cy)*(bx-cx)); 
	
	return area;
}


double Triangle::GetAspectRatio()
{
	double area = this->GetArea();
	double aspectratio=0.;

	
    double ax = (double)(this->v1->x);
    double ay = (double)(this->v1->y);
    double bx = (double)(this->v2->x);
    double by = (double)(this->v2->y);
    double cx = (double)(this->v3->x);
    double cy = (double)(this->v3->y);
	
	//compute the edge lengths of the triangle
	double l1 = sqrt((bx-ax)*(bx-ax)+(by-ay)*(by-ay));
	double l2 = sqrt((cx-ax)*(cx-ax)+(cy-ay)*(cy-ay));
	double l3 = sqrt((bx-cx)*(bx-cx)+(by-cy)*(by-cy));
  
	//maximal edge length of the triangle
	double lmax = l1;
	if(l2>lmax)
		lmax = l2;
	if(l3>lmax)
		lmax = l3;
						
	//computes the aspect ratio of the triangle
	aspectratio = lmax*(l1+l2+l3)/(4.0*sqrt(3)*area);
	
	return aspectratio;
}


void Triangle::debug()
{
    cout<<this->v1->x<<","<<this->v1->y<<"-"<<this->v2->x<<","<<this->v2->y<<"-"<<this->v3->x<<","<<this->v3->y<<endl;
}



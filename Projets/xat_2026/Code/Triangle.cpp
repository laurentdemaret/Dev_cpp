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
	this->point1 = p1;
	if(point1->ccw(p2,p3))
	{
		this->point2 = p2;
		this->point3 = p3;
	}
	else
	{
		this->point2 = p3;
		this->point3 = p2;
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
	this->point1 = NULL;	
	this->point2 = NULL;
	this->point3 = NULL;

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
	return Triangle::makeTriangle(triangle->point1,triangle->point2,triangle->point3);
}


Triangle* Triangle::makeTriangle(Edge* edge,Point2D* point3)
{
	return Triangle::makeTriangle(edge->org,edge->dest,point3);
}

bool Triangle::isVertex(Point2D* p)
{
	return (this->point1 == p) || (this->point2 == p) || (this->point3 == p);
}


bool Triangle::insideTriangle(Point2D* point)
{
	return this->insideTriangle(point->x,point->y);
}


// check if this point inside or on one edge 
// of the triangle builded by the 3 points
// parameters: point1, point2, point3
bool Triangle::insideTriangle(double px, double py)
{
	int p1x = this->point1->x;
	int p1y = this->point1->y;
	
	int p2x;
	int p2y;
	int p3x;
	int p3y;
	if(this->point1->ccw(this->point2,this->point3))
	{
		p2x = this->point2->x;
		p2y = this->point2->y;
		p3x = this->point3->x;
		p3y = this->point3->y;
	}
	else{
		p2x = this->point3->x;
		p2y = this->point3->y;
		p3x = this->point2->x;
		p3y = this->point2->y;
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
	Edge* entry = this->point1->entry;
	if(entry->isVertex(this->point2))
	{
		*e1 = entry;
		count++;
	}
	if(entry->isVertex(this->point3))
	{
		*e3 = entry;
		count++;
	}		
	Edge* orbit = entry->getNext(this->point1);
	while ( count < 2 && orbit != entry ) 
	{
		if(orbit->isVertex(this->point2))
		{
			*e1 = orbit;
			count++;
		}
		if( orbit->isVertex(this->point3))
		{
			*e3 = orbit;
			count++;
		}
		orbit = orbit->getNext(this->point1);
	}
	entry = this->point2->entry;
	if(entry->isVertex(this->point3))
	{
		*e2 = entry;
	}		
	orbit = entry->getNext(this->point2);
	while ( orbit != entry ) {
		if( orbit->isVertex(this->point3) )
		{
			*e2 = orbit;
			break;
		}
		orbit = orbit->getNext(this->point2);
	}
}


Point2D* Triangle::getOtherPoint(Point2D* a, Point2D* b)
{
	if(this->point1 != a 
	&& this->point1 != b )
	{
		return this->point1;
	}
	if( this->point2 != a 
	&&  this->point2 != b )
	{
		return this->point2;
	}
	else
	{
		return this->point3;
	}
}


Point2D* Triangle::getOtherPoint(Edge* edge)
{
	return this->getOtherPoint(edge->org, edge->dest);
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
	double ax = (double)(this->point1->x);
	double ay = (double)(this->point1->y);
	double bx = (double)(this->point2->x);
	double by = (double)(this->point2->y);
	double cx = (double)(this->point3->x);
	double cy = (double)(this->point3->y);

	double area =  0.5*((ax-cx)*(by-cy)-(ay-cy)*(bx-cx)); 
	
	return area;
}


double Triangle::GetAspectRatio()
{
	double area = this->GetArea();
	double aspectratio=0.;

	
	double ax = (double)(this->point1->x);
	double ay = (double)(this->point1->y);
	double bx = (double)(this->point2->x);
	double by = (double)(this->point2->y);
	double cx = (double)(this->point3->x);
	double cy = (double)(this->point3->y);
	
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
	cout<<this->point1->x<<","<<this->point1->y<<"-"<<this->point2->x<<","<<this->point2->y<<"-"<<this->point3->x<<","<<this->point3->y<<endl;
}

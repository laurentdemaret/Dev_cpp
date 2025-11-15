#include <cmath>
#include <vector>
#include <list>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "Point2D.h"
#include "Edge.h"
#include "Triangle.h"

using namespace std;


vector<Point2D*> Point2D::pool;
int Point2D::new_point2d_cnt = 0;

Point2D::Point2D(int cx, int cy)
{
	this->x = cx;
	this->y = cy;
}

Point2D::Point2D(int cx, int cy, int cv) 
{
		this->x = cx;
		this->y = cy;
		this->f = cv;
		this->thinned = false;
		this->entry = NULL;
		this->epsilon = __DBL_MAX__;
		this->clone = NULL;
		this->fh_el = NULL;
		this->exchange_fh_el = NULL;
		this->dirty = false;
		this->exchanged = false;
}

Point2D::~Point2D()
{
}


inline bool Point2D::lessThan(Point2D* b)
{
	if ( this->x < b->x ) return true;
	if ( this->x > b->x ) return  false;

	if ( this->y < b->y ) return true;
	return  false;
}

void Point2D::deletePoint2D()
{
	// refactor Point2D 
	this->entry = NULL;	
	this->bestNeighbor = NULL;
	this->thinnedTriangle = NULL;
	this->bestNeighborDiff = 0.;
	this->epsilon = 0.;
	if(this->clone != NULL){
		this->clone->clone = NULL;
		this->clone = NULL;
	}
	this->fh_el = NULL;
	
	if(Point2D::pool.size() < POINT2D_POOL_SIZE){
		Point2D::pool.push_back(this);
	}
	else{
		/* free edge */
		delete this;	
	}
}

Point2D* Point2D::makePoint2D(int cx, int cy, int cv)
{
	if(Point2D::pool.size() > 0)
	{
		Point2D* p = *(Point2D::pool.begin());
		Point2D::pool.erase(Point2D::pool.begin());
		p->x = cx;
		p->y = cy;
		p->f = cv;
		p->clone = NULL;
		p->fh_el = NULL;
		p->exchange_fh_el = NULL;
		p->gridBound = false;
		p->dirty = false;
		p->exchanged = false;
		return p;
	}
	else
	{
		Point2D::new_point2d_cnt++;
		return new Point2D(cx,cy,cv);	
	}
}
	
Point2D* Point2D::clonePoint2D()
{
	Point2D* c = Point2D::makePoint2D(this->x, this->y, this->f);
	this->clone = c;
	c->clone = this;
	return c;
}


vector<Point2D*> Point2D::clonePoint2D(vector<Point2D*> arr)
{
	vector<Point2D*> ret;
	for (unsigned int i = 0; i < arr.size(); i++) {
		Point2D* c = arr[i]->clonePoint2D();
		ret.push_back(c);
	}
	return ret;
}

/* ************************ geometric operations *********************** */
int Point2D::crossProduct(Point2D* p2, Point2D* p3)
{
	return ((p2->x - this->x) * (p3->y - this->y) 
				- (p2->y - this->y) * (p3->x - this->x));
}

int Point2D::crossProduct(Point2D* p2, Point2D* p3, Point2D* p4)
{
	int v1x = p2->x - this->x;
	int v1y = p2->y - this->y;
	int v2x = p4->x - p3->x;
	int v2y = p4->y - p3->y;
	
	return v1x*v2y - v1y*v2x;
	
}

int Point2D::dotProduct(Point2D* p2, Point2D* p3, Point2D* p4)
{
	int u1 = p2->x - this->x;
	int v1 = p2->y - this->y;
	int u2 = p4->x - p3->x;
	int v2 = p4->y - p3->y;
	
	return u1 * u2 + v1 * v2;
}
	
double Point2D::sin(Point2D* p2, Point2D* p3, Point2D* p4)
{
	int u1 = p2->x - this->x;
	int v1 = p2->y - this->y;
	int u2 = p4->x - p3->x;
	int v2 = p4->y - p3->y;
	
	return ((double)(u1 * u2 + v1 * v2))/(sqrt(u1*u1+v1*v1)*sqrt(u2*u2+v2*v2));
}

double Point2D::computeAngleCos(Point2D* a, Point2D* c)
{
	int ax = a->x;
	int ay = a->y;
	int cx = c->x;
	int cy = c->y;
	return this->computeAngleCos(ax,ay,cx,cy);
}


double Point2D::computeAngleCos(int ax, int ay, int cx, int cy)
{
	Point2D* b = this;
	int bx = b->x;
	int by = b->y;
	int bax = ax-bx;
	int bay = ay-by;
	int bcx = cx-bx;
	int bcy = cy-by;
	/* calculate angle */
//	double t1 = bax*bcx + bay*bcy;
//	double t2 = (bax*bax+bay*bay) * (bcx*bcx+bcy*bcy);
//	double t3 = sqrt( t2 );
//	double angleCos =  t1/t3;
	double angleCos =  ( bax*bcx + bay*bcy ) 
				/ (sqrt( (bax*bax+bay*bay) * (bcx*bcx+bcy*bcy)) );
	return angleCos;
}


int Point2D::inCircle(Point2D* point1, Point2D* point2, Point2D* point3)
{
	int px = this->x;
	int py = this->y;
	int p1x = point1->x;
	int p1y = point1->y;
	int p2x = -1;
	int p2y = -1;
	int p3x = -1;
	int p3y = -1;
	if(point1->ccw(point2,point3)){
		p2x = point2->x;
		p2y = point2->y;
		p3x = point3->x;
		p3y = point3->y;
	}
	else{
		p2x = point3->x;
		p2y = point3->y;
		p3x = point2->x;
		p3y = point2->y;
	}
	
	int adx = p1x - px;
	int bdx = p2x - px;
	int cdx = p3x - px;
	int ady = p1y - py;
	int bdy = p2y - py;
	int cdy = p3y - py;
	int bdxcdy = bdx * cdy;
	int cdxbdy = cdx * bdy;
	int alift  = adx * adx + ady * ady;
	
	int cdxady = cdx * ady;
	int adxcdy = adx * cdy;
	int blift  = bdx * bdx + bdy * bdy;
	
	int adxbdy = adx * bdy;
	int bdxady = bdx * ady;
	int clift  = cdx * cdx + cdy * cdy;
	
	int det = alift*(bdxcdy - cdxbdy) + 
				blift*(cdxady - adxcdy) + 
				clift*(adxbdy - bdxady);

	if(det > 0.){
		return 1;
	}
	if(det < 0.){
		return -1;
	}
	else{
		return 0;	
	}
	int ret = (det<0)? -1 : ((det>0)? 1 : 0 );
	return ret;

}

bool Point2D::ccw(Point2D* point1, Point2D* point2){
	
	int ax = this->x;
	int ay = this->y;
	int bx = point1->x;
	int by = point1->y;
	int cx = point2->x;
	int cy = point2->y;
	return ( (bx-ax)*(cy-ay) + (ax-cx)*(by-ay) ) > 0;
}

/* check if ab has score > cd */
bool Point2D::checkScore(Point2D* b, Point2D* c, Point2D* d){

	if( this->lessThan(b) ){
		return ( this->lessThan(c) && this->lessThan(d) );
	}
	else{
		return ( b->lessThan(c) && b->lessThan(d) );
	}

}

vector<Point2D*> Point2D::getNeighbors()
{
	vector<Point2D*> neighbors;
	neighbors.push_back(this->entry->getOtherPoint(this));
	Edge* orbit = this->entry->getNext(this);
	while ( this->entry != orbit ){
		neighbors.push_back(orbit->getOtherPoint(this));
		orbit = orbit->getNext(this);
	}
	return neighbors;
}

vector<Point2D*> Point2D::getNeighborsCopy()
{
	return Point2D::clonePoint2D(this->getNeighbors());
}

vector<Triangle*> Point2D::getTriangles()
{
	vector<Triangle*> triangles;
	Point2D* last_neighbor = entry->getOtherPoint(this);
	Edge* orbit = entry->getNext(this);
	while ( entry != orbit ){
		if(orbit->checkTriangle(last_neighbor)){
			triangles.push_back(Triangle::makeTriangle(orbit,last_neighbor));
		}
		last_neighbor = orbit->getOtherPoint(this);
		orbit = orbit->getNext(this);
	}
	/* last time */
	if(orbit->checkTriangle(last_neighbor)){
		triangles.push_back(Triangle::makeTriangle(orbit,last_neighbor));
	}
	return triangles;
}


vector<Edge*> Point2D::getCellHull(int xMax, int xMin, int yMax, int yMin)
{
	vector<Edge*> cell;
	Edge* entry = this->entry;
	Edge* edge = entry->getNext(entry->getOtherPoint(this)); 
	if(entry->isConvexHull(xMax,xMin,yMax,yMin)){
		cell.push_back(entry);
		if(edge->checkTriangle(this)){
			cell.push_back(edge);
		}
	}
	else{
		cell.push_back(edge);
	}
	Edge* orbit = entry->getNext(this);
	while ( entry != orbit ){
		edge = orbit->getNext(orbit->getOtherPoint(this));
		if(orbit->isConvexHull(xMax,xMin,yMax,yMin)){
			cell.push_back(orbit);
			if(edge->checkTriangle(this)){
				cell.push_back(edge);
			}
		}
		else{
			cell.push_back(edge);
		}
		
		orbit = orbit->getNext(this);
	}
	return cell;
}

bool Point2D::isConvexHull(Point2D* b, int xMax, int xMin, int yMax, int yMin)
{
	if (this->convexHull && b->convexHull){
		/* check same side */
		if( this == b ){
			return this->convexHull;
		}
		if( this->y == b->y ){
			return this->y == yMin || this->y == yMax;
		}
		if( this->x == b->x ){
			return this->x == xMin || this->x == xMax;
		}
		return false;
	}
	return false;
}


vector<Edge*> Point2D::getEdges()
{
	vector<Edge*> edges;
	edges.push_back(this->entry);
	Edge* orbit = this->entry->getNext(this);
	while ( this->entry != orbit ){
		edges.push_back(orbit);
		orbit = orbit->getNext(this);
	}
	return edges;
}

Edge* Point2D::getEdge(Point2D* b)
{
	Edge* entry = this->entry;
	if(entry->isVertex(b)){
		return entry;
	}
	Edge* orbit = entry->getNext(this);
	while ( entry != orbit ){
		if( orbit->isVertex(b) ){
			return orbit;
		}
		orbit = orbit->getNext(this);
	}
	return NULL;
}

/* ************************************* HELPER ****************************** */
vector<Point2D*> Point2D::set2vector(set<Point2D*> pSet)
{
	vector<Point2D*> ret;
	set<Point2D*>::iterator iter = pSet.begin();
	for( iter = pSet.begin(); iter != pSet.end(); iter++ ){
		ret.push_back(*iter);
	}
	return ret;
}
set<Point2D*> Point2D::vector2set(vector<Point2D*> vec)
{
	set<Point2D*> ret;
	for( unsigned int i = 0; i< vec.size(); i++){
		ret.insert(vec[i]);
	}
	return ret;
}
set<Point2D*> Point2D::makeUniqueList(set<Point2D*> pSet, set<Point2D*> pSet2)
{
	set<Point2D*> ret;
	set<Point2D*>::iterator iter = pSet.begin();
	for( iter = pSet.begin(); iter != pSet.end(); iter++ ){
		ret.insert(*iter);
	}
	iter = pSet2.begin();
	for( iter = pSet2.begin(); iter != pSet2.end(); iter++ ){
		ret.insert(*iter);
	}
	return ret;	
}

vector<Point2D*> Point2D::makeUniqueList(vector<Point2D*> vec, vector<Point2D*> vec2)
{
	set<Point2D*> s1 = Point2D::vector2set(vec);
	set<Point2D*> s2 = Point2D::vector2set(vec2);
	set<Point2D*> s3 = Point2D::makeUniqueList(s1, s2);
	return Point2D::set2vector(s3);
}

void Point2D::addAll(set<Point2D*> *s, vector<Point2D*> v)
{
	for(unsigned int i = 0; i < v.size(); i++){
		s->insert(v[i]);
	}
}

void Point2D::addAll(set<Point2D*> *s1, set<Point2D*> s2)
{
	set<Point2D*>::iterator i = s2.begin();
	for(i = s2.begin(); i != s2.end(); i++){
		s1->insert(*i);
	}
}


void Point2D::debug()
{
	cout<<this->x<<","<<this->y<<endl;
}

void Point2D::debugNeighbors(){
	vector<Point2D*> neighbors = this->getNeighbors();
	for (unsigned int i = 0; i < neighbors.size(); i++) {
		neighbors[i]->debug();
	}	
}

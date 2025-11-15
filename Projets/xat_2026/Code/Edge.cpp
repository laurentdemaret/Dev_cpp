#include <vector>
#include <set>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "Edge.h"
#include "Point2D.h"
#include "Triangle.h"


using namespace std;


vector<Edge*> Edge::pool;
int Edge::new_edge_cnt = 0;

Edge::Edge()
{
  this->initEdge();
}

Edge::Edge(Point2D* o, Point2D* d)
{
  this->org = o;
  this->dest = d;
  this->initEdge();
  if (this->org->entry == NULL)
  {
    this->org->entry = this;
  }
  if (this->dest->entry == NULL)
  {
    this->dest->entry = this;
  }
}

Edge::~Edge(){}

void Edge::initEdge()
{
	this->onext = this;
	this->oprev = this;
	this->dnext = this;
	this->dprev = this;		
    // linked list
	this->pred = this;
	this->succ = this;
}


Edge* Edge::getNext(Point2D* p)
{
	return (this->org == p)? this->onext : this->dnext;
}


Edge* Edge::getPrev(Point2D* p)
{
	return (this->org == p)? this->oprev : this->dprev;
}

void Edge::setNext(Point2D* p,Edge* e)
{
	if(this->org == p){
		this->onext = e;	
	}
	else{
		this->dnext = e;	
	}
}


void Edge::setPrev(Point2D* p,Edge* e)
{
	if(this->org == p){
		this->oprev = e;	
	}
	else{
		this->dprev = e;	
	}
}

Point2D* Edge::getOtherPoint(Point2D* p)
{
  return (p == this->org)? this->dest : this->org;
}


bool Edge::assumeUniquness(int level)
{
//	if(level > 3) return false;
	Point2D* np = this->onext->getOtherPoint(this->org);
	if( np != this->dprev->getOtherPoint(this->dest) ){
		return false;
	}
	Point2D* pp = this->dnext->getOtherPoint(this->dest);
	if( pp == np ) return false;
	if( pp != this->oprev->getOtherPoint(this->org) ){
		return false;
	}
	
	if( pp->checkScore(np, this->org, this->dest)
    &&  pp->inCircle(np, this->org, this->dest) == 0. )
    {
		this->swap();
		int l1 = level+1;
		this->onext->assumeUniquness(l1);
		this->dnext->assumeUniquness(l1);
		this->oprev->assumeUniquness(l1);
		this->dprev->assumeUniquness(l1);
		return true;
	}
	return false;
}
	
	
void Edge::swap(){
	
	Edge* a = this->oprev;
	Edge* b = this->dprev;
	Edge* c = this->onext;
	Edge* d = this->dnext;
    // store org+dest
	Point2D* nOrg = a->getOtherPoint(this->org);
	Point2D* nDest = b->getOtherPoint(this->dest);
	
    // adjust oprev
    if( a->org == this->org )
    {
		a->onext = this->onext;
		a->dprev = this;
	}
	else{
		a->dnext = this->onext;
		a->oprev = this;
	}
    // adjust dprev
    if( b->org == this->dest )
    {
		b->dprev = this;
		b->onext = this->dnext;
	}
    else
    {
      b->oprev = this;
      b->dnext = this->dnext;
	}
    // adjust onext
    if( c->org  == this->org )
    {
      c->oprev = a;
      c->dnext = this;
	}
    else
    {
      c->dprev = a;
      c->onext = this;
	}
    // adjust dnext
    if( d->org == this->dest )
    {
      d->oprev = b;
      d->dnext = this;
	}
	else{
		d->dprev = b;
		d->onext = this;
	}
	/* adjust new links */
	this->onext = a;
	this->dnext = b;
	this->oprev = d;
	this->dprev = c; 
	/* adjust point entry */
	if( this == this->org->entry ){
		this->org->entry = a;
	}
	if( this == this->dest->entry ){
		this->dest->entry = b;
	}
	this->org = nOrg;
	this->dest = nDest;
	this->reverseEdge();
}
	
	
	
void Edge::deleteEdgeInsideCell(set<Edge*> cell)
{
	Edge* on = this->onext;
	Edge* op = this->oprev;
	Edge* dn = this->dnext;
	Edge* dp = this->dprev;
	
	if(cell.find(on) != cell.end() ){
		on->setPrev(this->org, op);
		op->setNext(this->org, on);
	}
	if(cell.find(op) != cell.end() ){
		op->setNext(this->org, on);
		on->setPrev(this->org, op);
	}
	if(cell.find(dn) != cell.end() ){
		dn->setPrev(this->dest, dp);
		dp->setNext(this->dest, dn);
	}
	if(cell.find(dp) != cell.end() ){
		dp->setNext(this->dest, dn);
		dn->setPrev(this->dest, dp);
	}
}
	
		
//  @param e new edge
// @param p crossing point
void Edge::connect(Edge* e, Point2D* p)
{
	Point2D* o = e->getOtherPoint(p);
	Edge* next = this->getNext(p);
	Edge* prev = this->getPrev(p);
	if(o->ccw(p,this->getOtherPoint(p))){
		e->setNext(p,next);
		e->setPrev(p,this);
		this->setNext(p,e);
		next->setPrev(p,e);
	}
	else{
		e->setNext(p,this);
		e->setPrev(p,prev);
		this->setPrev(p,e);
		prev->setNext(p,e);
	}
}

void Edge::adjustQuadStructure(Point2D* point1, 
						 Edge* next1, 
						 Edge* prev1, 
						 Point2D* point2, 
						 Edge* next2, 
						 Edge* prev2){
	this->setNext(point1, next1);
	this->setNext(point2, next2);
	this->setPrev(point1, prev1);
	this->setPrev(point2, prev2);
	next1->setPrev(point1, this);
	next2->setPrev(point2, this);
	prev1->setNext(point1, this);
	prev2->setNext(point2, this);
}


// Add an edge to a ring of edges.
// @param b
// @param v
void Edge::splice(Edge* b, Point2D* v)
{
	Edge* next;

	  
	/* b must be the unnattached edge and a must 
	 * be the previous ccw edge to b. */
	if ( this->org == v ) {
		next = this->onext;
		this->onext = b;
	}
	else {
		next = this->dnext;
		this->dnext = b;
	}

	if ( next->org == v ){
		next->oprev = b;
	}
	else{
		next->dprev = b;
	}
	
	if ( b->org == v ) {
		b->onext = next;
		b->oprev = this;
	}
	else {
		b->dnext = next;
		b->dprev = this;
	}
}


/**
 * Creates a new edge and adds it to two rings of edges.
 * @param u
 * @param b
 * @param v
 * @param s
 * @return
 */
Edge* Edge::join(Point2D* u, Edge* b, Point2D* v, int s)
{
	/*
	 * u and v are the two vertices which are being joined. a and b are the
	 * two edges associated with u and v res.
	 */

	Edge* e = Edge::makeEdge(u, v);

	if (s == LEFT) {
		if ( this->org == u ) {
			this->oprev->splice(e, u);
		} else {
			this->dprev->splice(e, u);
		}
		b->splice(e, v);
	} else {
		this->splice(e, u);
		if ( b->org == v ){
			b->oprev->splice(e, v);
		}
		else{
			b->dprev->splice(e, v);
		}
	}
	return e;
}


vector<Point2D*> Edge::getNeighbors()
{
	vector<Point2D*> neighbors;
	/* orbit over the org */
	Edge* orbit = this->getNext(this->org);
	while ( this != orbit ){
		neighbors.push_back(orbit->getOtherPoint(this->org));
		orbit = orbit->getNext(this->org);
	}
	/* orbit over the dest: do not register next 
	 * because already registred by the org orbit
	 * */
	orbit = this->getNext(this->dest)->getNext(this->dest);
	Edge* last = this->getPrev(this->dest);
	while ( last != orbit ){
		neighbors.push_back(orbit->getOtherPoint(this->dest));
		orbit = orbit->getNext(this->dest);
	}
	return neighbors;
}


vector<Edge*> Edge::getEdges()
{
	vector<Edge*> edges;
	/* orbit over the org */
	Edge* orbit = this->onext;
	while ( this != orbit ){
		edges.push_back(orbit);
		orbit = orbit->getNext(this->org);
	}
	/* orbit over the dest: do not register next 
	 * because already registred by the org orbit
	 * */
	orbit = this->dnext->getNext(this->dest);
	Edge* last = this->dprev;
	while ( last != orbit ){
		edges.push_back(orbit);
		orbit = orbit->getNext(this->dest);
	}
	return edges;
}

vector<Edge*> Edge::getCellHull(int xMax, int xMin, int yMax, int yMin)
{
	vector<Edge*> cell;
	
	vector<Edge*> co = this->org->getCellHull(xMax,xMin,yMax,yMin);
	for (unsigned int i = 0; i < co.size(); i++) {
		if( !co[i]->isVertex(this->dest) || co[i]->isConvexHull(xMax,xMin,yMax,yMin) ){
			cell.push_back(co[i]);
		}
	}
	
	vector<Edge*> cd = this->dest->getCellHull(xMax,xMin,yMax,yMin);
	for (unsigned int i = 0; i < cd.size(); i++) {
		if( ( !cd[i]->isVertex(this->org)	|| cd[i]->isConvexHull(xMax,xMin,yMax,yMin) )
        && find(cell.begin(), cell.end(), cd[i]) == cell.end() )
        {
			cell.push_back(cd[i]);
		}
	}	

	return cell;
}



bool Edge::checkEdge(Triangle* triangle, int level)
{
	
	Triangle *to, *td;
//	if(level > 3) return false;
	Point2D* adjPoint = triangle->getNeighbourPointByEdge(this);
	if(adjPoint != NULL){
		Point2D* po = this->org;
		Point2D* pd = this->dest;
		Edge* eo = adjPoint->getEdge(po);
		Edge* ed = adjPoint->getEdge(pd);
		Point2D* point = triangle->getOtherPoint(this);
		
		switch (adjPoint->inCircle(triangle->point1,
									triangle->point2,
									triangle->point3 )) {
		case 1:
			this->swap();
			to = Triangle::makeTriangle(po,point,adjPoint);
			td = Triangle::makeTriangle(pd,point,adjPoint);
			eo->checkEdge(to, level+1);
			ed->checkEdge(td, level+1);
			to->deleteTriangle();
			td->deleteTriangle();
			return true;
		case 0:
			if(point->checkScore(adjPoint,po,pd)){
				this->swap();
				//TODO ??? recursive ??? 
				eo->assumeUniquness(level+1);
				ed->assumeUniquness(level+1);
			}
			break;
		default:
			return false;
		}
	}
	else{
		return false;
	}
	return false;
}


//  Remove an edge. 
void Edge::deleteEdge()
{
	// Cache origin and destination 
	Point2D* u = this->org;
	Point2D* v = this->dest;

	// Adjust entry points 
	if ( u->entry  == this ){
		u->entry = this->onext;
	}
	if ( v-> entry == this ){
		v->entry = this->dnext;
	}

    // Four edge links to change 
	if( this->onext->org == u ){
		this->onext->oprev = this->oprev;
	}
	else{
		this->onext->dprev = this->oprev;
	}

	if( this->oprev->org == u ){
		this->oprev->onext =this->onext;
	}
	else{
		this->oprev->dnext = this->onext;
	}

	if( this->dnext->org == v ){
		this->dnext->oprev = this->dprev;
	}
	else{
		this->dnext->dprev = this->dprev;
	}

	if( this->dprev->org == v ){
		this->dprev->onext = this->dnext;
	}
	else{
		this->dprev->dnext = this->dnext;
	}
	
	this->refactorEdge();

}


void Edge::refactorEdge()
{	
	// refactor edge 
	this->org = NULL;	
	this->dest = NULL;
	this->fh_el = NULL;
	this->onext = this;
	this->oprev = this;
	this->dnext = this;
	this->dprev = this;
	this->epsilon = 0.;
	this->succ = this;
	this->pred = this;

	if( Edge::pool.size() < EDGE_POOL_SIZE ){
		Edge::pool.push_back(this);
	}
	else{
		/* free edge */
		delete this;
	}
}

void Edge::fillEdgePool()
{
	unsigned int pool_size = Edge::pool.size();
	for (unsigned int i = pool_size; i < EDGE_POOL_SIZE ; i++) {
		Edge::pool.push_back(new Edge());
	}
}

Edge* Edge::makeDummyEdge()	
{
	if( Edge::pool.size() > 0 ){
		Edge* e = *(Edge::pool.begin());
		Edge::pool.erase(Edge::pool.begin());
		return e;
	}
	else{
		Edge::new_edge_cnt++;
		return new Edge();	
	}
}
	
Edge* Edge::makeEdge(Point2D* o, Point2D* d)	
{

	if(Edge::pool.size() > 0){
		Edge* e = *(Edge::pool.begin());
		Edge::pool.erase(Edge::pool.begin());
		e->org = o;
		e->dest = d;
		return e;
	}
	else{
		Edge::new_edge_cnt++;
		return new Edge(o,d);	
	}
}
	
	
	
	
/**
 * check if this triangle is valid true=valid, false=invalid
 * @param point
 * @return triangle is ok
 */
bool Edge::checkTriangle(Point2D* point)
{
	return /*!this->isColinear(point) && */ 
			( ( this->onext->getOtherPoint(this->org) == point
				&& this->dprev->getOtherPoint(this->dest) == point ) )
		||
			  ( this->oprev->getOtherPoint(this->org) == point 
			 	&& this->dnext->getOtherPoint(this->dest) == point );
}

bool Edge::isVertex(Point2D* p)
{
	return (this->org == p) || (this->dest == p);
}

bool Edge::isVirtualVertex(Point2D* p)
{
	return ( this->org->x == p->x && this->org->y == p->y ) 
	    || ( this->dest->x == p->x && this->dest->y == p->y );
}

bool Edge::isVirtualEqual(Edge* e)
{
	return ( this->isVirtualVertex(e->org) && this->isVirtualVertex(e->dest) );
}
	
// @return Returns the convexHull.
bool Edge::isConvexHull(unsigned int xMax, 
						unsigned int xMin, 
						unsigned int yMax, 
						unsigned int yMin)
{
	return this->org->isConvexHull(this->dest,xMax,xMin,yMax,yMin);
}	
	
void Edge::debug()
{
	cout<<this->org->x<<","<<this->org->y<<"-"<<this->dest->x<<","<<this->dest->y<<endl;
}	
void Edge::debugStruct()
{
	cout<<endl<<"this: ";
	this->debug();
	cout<<endl<<"onext: ";
	this->onext->debug();
	cout<<endl<<"dnext: ";
	this->dnext->debug();
	cout<<endl<<"oprev: ";
	this->oprev->debug();
	cout<<endl<<"dprev: ";
	this->dprev->debug();
	cout<<endl;
}

	
/* **************************************** PRIVATE ************************************ */	
bool Edge::reverseEdge(){
	Edge* on = this->onext;
	Edge* dn = this->dnext;
	Edge* op = this->oprev;
	Edge* dp = this->dprev;
	Point2D* org = this->org;
	Point2D* dest = this->dest;
	if( this != on && this != dp 
	&& ( ( on->dest == this->org && dp->org == this->dest )  
		|| ( op->dest == this->org && dn->org == this->dest ) ) ){
		this->onext = dn;
		this->oprev = dp;
		this->dnext = on;
		this->dprev = op;
		this->org = dest;
		this->dest = org;
		return true;
	}
	return false;
}	



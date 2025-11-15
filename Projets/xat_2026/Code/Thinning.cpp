#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <set>
#include <iostream>
#include <time.h>
#include <string>
#include <fstream>
#include <math.h>

#include "Thinning.h"

#include "Point2D.h"
#include "Edge.h"
#include "Triangle.h"
#include "Test.h"

using namespace std;

int Thinning::exchangeRadius;

Thinning::~Thinning()
{
}

/* **************************************** COMPARATORS **************************************** */
inline int point_significance_cmp( void *x, void *y )
{
  Point2D* a = (Point2D*)x;
  Point2D* b = (Point2D*)y;

  if ( a->epsilon < b->epsilon ) return -1;
  if ( a->epsilon > b->epsilon ) return 1;

  return 0;
}


inline int edge_significance_cmp( void *x, void *y )
{
  Edge* a = (Edge*)x;
  Edge* b = (Edge*)y;

  if ( a->epsilon < b->epsilon ) return -1;
  if ( a->epsilon > b->epsilon ) return 1;

  return 0;
}

inline int neighbor_significance_cmp( void *x, void *y )
{
  Point2D* a = (Point2D*)x;
  Point2D* b = (Point2D*)y;

  if ( a->bestNeighborDiff < b->bestNeighborDiff ) return 1;
  if ( a->bestNeighborDiff > b->bestNeighborDiff ) return -1;

  return 0;
}


Thinning::Thinning(Triangulation* tri)
{
	this->triangulation = tri;
	// init point heap 
	this->nodeHeap = fh_makeheap();
	fh_setcmp(this->nodeHeap, point_significance_cmp);
	Point2D* empty_point = Point2D::makePoint2D(0,0,0);
	empty_point->epsilon = -__DBL_MAX__;
	fh_setneginf(this->nodeHeap, (void*)empty_point);

	// init edge heap 
	this->edgeHeap = fh_makeheap();
	fh_setcmp(this->edgeHeap, edge_significance_cmp);
	Edge* empty_edge = new Edge();
	empty_edge->epsilon = -__DBL_MAX__;
	fh_setneginf(this->edgeHeap, (void*)empty_edge);
	
	// init thinnedNeighborHeap heap 
	this->thinnedNeighborHeap = fh_makeheap();
	fh_setcmp(this->thinnedNeighborHeap, neighbor_significance_cmp);
	Point2D* empty_point_n = Point2D::makePoint2D(0,0,0);
	empty_point_n->bestNeighborDiff = __DBL_MAX__;
	fh_setneginf(this->thinnedNeighborHeap, (void*)empty_point_n);	
}

double Thinning::calculateNewError(vector<Triangle*> nTriangles,	
		 							 Point2D* point)
{
	double max = 0.;
	for (unsigned int i = 0; i < nTriangles.size(); i++) 
	{
		double err = calculateError(nTriangles[i], point);
		if (max <= err) 
		{
			max = err;
			return max;
		}
	}
	return max;
}


// This function calculate the significance of an edge, for attachedPoints, old triangles (oTriangles)
//and new triangles (nTriangles)
double Thinning::calculateSignificance(vector<Triangle*> oTriangles, vector<Triangle*> nTriangles,
										 vector<Point2D*> attachedPoints, Edge* edge)
{
	double newError = this->calculateNewError(nTriangles, edge->org);
	newError += this->calculateNewError(nTriangles, edge->dest);
	
	double max = 0.;
	double formerError = 0.;

	for (unsigned int i = 0; i < attachedPoints.size(); i++) 
	{
		max = 0.;
		for (unsigned int j = 0; j < oTriangles.size(); j++) 
		{
			double err = calculateError(oTriangles[j], attachedPoints[i]);
			if (max < err) 
			{
				max = err;
				break;
			}
		}
		// increase former error
		formerError += max;
		
		max = 0.;
		for (unsigned int j = 0; j < nTriangles.size(); j++) 
		{
			double err = calculateError(nTriangles[j], attachedPoints[i]);
			if (max < err) 
			{
				max = err;
				break;
			}
		}
		newError += max;
	}
	
	return (newError>formerError)? newError-formerError : 0;
}


// @param oTriangles
// @param nTriangles
// @param attachedPoints
// @param edge
// @return
double Thinning::calculateSignificance(vector<Triangle*> oTriangles,
										 vector<Triangle*> nTriangles,
										 vector<Point2D*> attachedPoints,
										 Point2D* point)
{
	double max;
	double newError = calculateNewError(nTriangles, point);
//if(true) return newError;
	double formerError = 0.;
	for (unsigned int i = 0; i < attachedPoints.size(); i++) 
	{
		//TODO max = -1: komisch, dass die sig besser ist, wenn man alle
		//TODO points beruecksichtigt, auch die die nicht zur alten cell gehören...
		max = -1.;
		for (unsigned int j = 0; j < oTriangles.size(); j++) {
			double err = calculateError(oTriangles[j], attachedPoints[i]);
			if (max < err) {
				max = err;
				break;
			}
		}
		// attached point inside the old cell 
		if( max >= 0 ){
			// increase former error 
			formerError += max;
			
			max = 0.;
			for (unsigned int j = 0; j < nTriangles.size(); j++) 
			{
				double err = calculateError(nTriangles[j], attachedPoints[i]);
				if (max < err) 
				{
					max = err;
					break;
				}
			}
			newError += max;
		}
	}
	return (newError>formerError)? newError-formerError : 0;
}


double Thinning::calculateError(Triangle* triangle, Point2D* point)
{
	int x = point->x;
	int y = point->y;
	int v = point->f;
	
	Point2D* a = triangle->point1;
	Point2D* b = triangle->point2;
	Point2D* c = triangle->point3;
	
	int ax = a->x;
	int ay = a->y;
	
	int bx = b->x;
	int by = b->y;
	
	int cx = c->x;
	int cy = c->y;
	
	int ba = (x*(by-cy) + bx*(cy-y) + cx*(y-by));
	if(ba<0) return -1;
	
	int bb = (ax*(y-cy) + x*(cy-ay) + cx*(ay-y));
	if(bb<0) return -1;
	
	int bc = (ax*(by-y) + bx*(y-ay) + x*(ay-by));
	if(bc<0) return -1;
	
	double det = ax*(by-cy) + bx*(cy-ay) + cx*(ay-by);		
	
	int av = a->f;
	int bv = b->f;
	int cv = c->f;
	
	//This term is the quadratic error term 
	double tmp = (((ba*av)+(bb*bv)+(bc*cv))/det)-v;
	
	
	//Now computes an L1 penalisation of the solution
	double pen = 0., det2 = (bx-ax)*(cy-ay) -(cx-ax)*(by-ay);
	if(fabs(det2)>0) 
	{
		//compute the parameters dx and dy of: f(x,y) = dx*x+dy*y+k
		double dx = (( cy-ay)*(bv-av)-(by-ay)*(cv-av))/det2; //
		double dy = (( ax-cx)*(bv-av)+(bx-ax)*(cv-av))/det2; //
		pen = pow(dx*dx+dy*dy,0.5);
		//TODO: include tau as parameter and use as follows
		//tmp += pow((dx*dx+dy*dy),tau);
	}
	
	//double lambda = this->lambda;
	return tmp*tmp + (this->lambda)*pen;	
}


/*double Thinning::calculateError(Triangle* triangle, Point2D* point,double lambda)
{
	int x = point->x;
	int y = point->y;
	int v = point->f;

	Point2D* a = triangle->point1;
	Point2D* b = triangle->point2;
	Point2D* c = triangle->point3;
	
	int ax = a->x;
	int ay = a->y;
	
	int bx = b->x;
	int by = b->y;
	
	int cx = c->x;
	int cy = c->y;
	
	int ba = (x*(by-cy) + bx*(cy-y) + cx*(y-by));
	if(ba<0) return -1;
	
	int bb = (ax*(y-cy) + x*(cy-ay) + cx*(ay-y));
	if(bb<0) return -1;
	
	int bc = (ax*(by-y) + bx*(y-ay) + x*(ay-by));
	if(bc<0) return -1;
	
	double det = ax*(by-cy) + bx*(cy-ay) + cx*(ay-by);		
	
	int av = a->f;
	int bv = b->f;
	int cv = c->f;

	//This term is the quadratic error term 
	double tmp = (((ba*av)+(bb*bv)+(bc*cv))/det)-v;
  
	//Now computes an L1 penalisation of the solution
	double pen; 
	double det2 = (bx-ax)*(cy-ay) -(cx-ax)*(by-ay);
	if(fabs(det2)>0) 
	{
		//compute the parameters dx and dy of: f(x,y) = dx*x+dy*y+k
		double dx = (( cy-ay)*(bv-av)-(by-ay)*(cv-av))/det2; //
		double dy = (( ax-cx)*(bv-av)+(bx-ax)*(cv-av))/det2; //
		pen = pow(dx*dx+dy*dy,0.5);
		//TODO: include tau as parameter and use as follows
		//tmp += pow((dx*dx+dy*dy),tau);
	}
		
	return tmp*tmp + lambda*pen;
}*/


double Thinning::calculatePointSignificance(Point2D* point)
{
	// get neighbors 
	vector<Point2D*> neighborsCopy = point->getNeighborsCopy();
	vector<Triangle*> oTriangles = point->getTriangles();
	Triangulation* tri = Triangulation::makeTriangulation(neighborsCopy);
	vector<Triangle*> nTriangles = tri->getTriangles();		
	vector<Point2D*> attachedPoints = this->triangulation->getAttachedPoints(neighborsCopy);
	double sig = calculateSignificance(oTriangles,nTriangles,attachedPoints,point);

	// refactor objects 
	Triangle::deleteTriangles(oTriangles);

	tri->deleteTriangulation(nTriangles);
	return sig;
}

double Thinning::calculatePointSignificance_ForAT8(Point2D* point)
{
	// get neighbors 
	vector<Point2D*> neighborsCopy = point->getNeighborsCopy();
	vector<Triangle*> oTriangles = point->getTriangles();
	//Triangulation* tri = Triangulation::makeTriangulation(neighborsCopy);
	//vector<Triangle*> nTriangles = tri->getTriangles();		
	//vector<Point2D*> attachedPoints = this->triangulation->getAttachedPoints(neighborsCopy);
	//double sig = calculateSignificance(oTriangles,nTriangles,attachedPoints,point);
    double sig; 		
	if((int)(point->x)%3 ==0 && (int)(point->y)%3 ==0)
      sig = 12000.*rand()/RAND_MAX;
    else
	  sig = 1.*rand()/RAND_MAX; //random value

	/* refactor objects */
	Triangle::deleteTriangles(oTriangles);

	//tri->deleteTriangulation(nTriangles);
	return sig;
}


double Thinning::calculateEdgeSignificance(Edge* edge)
{
	/* get neighbors */
	vector<Point2D*> neighbors;
	vector<Point2D*> on = edge->org->getNeighbors();
	for (unsigned int i = 0; i < on.size(); i++) {
		if(on[i] != edge->dest){
			neighbors.push_back(on[i]);
		}
	}
	vector<Point2D*> dn = edge->dest->getNeighbors();
	for (unsigned int i = 0; i < dn.size(); i++ ) {
		if( find(neighbors.begin(), neighbors.end(), dn[i]) == neighbors.end() 
		 && dn[i] != edge->org){
			neighbors.push_back(dn[i]);
		}
	}
	
	vector<Triangle*> oTriangles = this->triangulation->getTriangles(edge);
	vector<Point2D*> neighborsCopy = Point2D::clonePoint2D(neighbors);
	Triangulation* tri = Triangulation::makeTriangulation(neighborsCopy);
	vector<Triangle*> nTriangles = tri->getTriangles();		
	vector<Point2D*> attachedPoints = this->triangulation->getAttachedPoints(neighbors);
	double sig = calculateSignificance(oTriangles,nTriangles,attachedPoints,edge);
	/* refactor objects */
	Triangle::deleteTriangles(oTriangles);

	tri->deleteTriangulation(nTriangles);

	return sig;
}

void Thinning::fastThinning(unsigned int tcount)
{
	/* build point heap */
	for (unsigned int i = 0; i < this->triangulation->nodes.size(); i++) 
	{
		this->calculateSignificanceForHeap(this->triangulation->nodes[i]);
	}
	for(unsigned int n = 0; n < tcount; n++) {
		/* get less significant node from point heap */
		Point2D* ph_min = (Point2D*)fh_extractmin(this->nodeHeap);
		
		/* check if current node is dirty */
		if(ph_min->dirty){
			n--;
			this->calculateSignificanceForHeap(ph_min);
			continue;
		}
		
		vector<Point2D*> neighbors = ph_min->getNeighbors();
		/* delete point */
		this->triangulation->deletePoint(ph_min);
		/* mark neighbors as dirty */
		for(unsigned int i = 0; i < neighbors.size(); i++){
			neighbors[i]->dirty = true;
		}
		
		
/* ********************** DEBUG ***************************** */		
cout<<"[deleting point] "<<n<<"/"<<tcount<<"  : "<<
ph_min->x<<","<<ph_min->y<<" -> "<<ph_min->epsilon<<": "<<endl;
/* ********************** DEBUG ***************************** */
		
	}
	
	// free mem 
	fh_deleteheap(this->nodeHeap);	
}

void Thinning::thinningAT5(unsigned int tcount)
{
	/* build point heap */
	for (unsigned int i = 0; i < this->triangulation->nodes.size(); i++) 
	{
		this->calculateSignificanceForHeap(this->triangulation->nodes[i]);
	}
	for(unsigned int n = 0; n < tcount; n++) 
	{
		/* get less significant node from point heap */
		Point2D* ph_min = (Point2D*)fh_extractmin(this->nodeHeap);
		ph_min->fh_el = NULL;
		vector<Point2D*> neighbors = ph_min->getNeighbors();
		/* delete point */
		this->triangulation->deletePoint(ph_min);
		/* recalculate significance of neighbors */
		for(unsigned int i = 0; i < neighbors.size(); i++)
		{
		  if(neighbors[i]->gridBound) continue;
		  fh_delete(this->nodeHeap, (fibheap_el*)neighbors[i]->fh_el);
		  neighbors[i]->fh_el = NULL;
		  this->calculateSignificanceForHeap(neighbors[i]);
		}
/* ********************** DEBUG ***************************** */		
//cout<<"[deleting point] "<<n<<"/"<<tcount<<"  : "<<
//ph_min->x<<","<<ph_min->y<<" -> "<<ph_min->epsilon<<endl;
/* ********************** DEBUG ***************************** */
	}
	
	/* free mem */
	fh_deleteheap(this->nodeHeap);	
	fh_deleteheap(this->edgeHeap);	
}


void Thinning::thinningAT6(unsigned int tcount)
{
    std::cout << "AT 6" << std::endl;
	/* build point heap */
	for (unsigned int i = 0; i < this->triangulation->nodes.size(); i++) 
	{
	  this->calculateSignificanceForHeap(this->triangulation->nodes[i]);
	}
	
	/* build edge heap */
	vector<Edge*> edges = this->triangulation->getEdges();
	for (unsigned int i = 0; i < edges.size(); i++) 
	{
		this->calculateSignificanceForEdgeHeap(edges[i]);
	}
	for (unsigned int n = 0; n < tcount; n++) 
	{
		/* get less significant node from point heap */
		Point2D* ph_min = (Point2D*)fh_extractmin(this->nodeHeap);
		Point2D* ph_min2 = (Point2D*)fh_min(this->nodeHeap);
		
		/* get less significant edge */
		Edge* eh_min = (Edge*)fh_min(this->edgeHeap);
		/* edge has less sig: edge based */
		if( eh_min->epsilon < ph_min->epsilon+ph_min2->epsilon ){
cout<<"[deleting edge] "<<n<<"/"<<tcount<<"  : "<<
eh_min->org->x<<","<<eh_min->org->y<<"-"<<
eh_min->dest->x<<","<<eh_min->dest->y<<
" -> "<<eh_min->epsilon<<": "<<endl;
			/* inc counter */
			n++;
			/* reinsert ph_min */
			ph_min->fh_el = (void*)fh_insert( this->nodeHeap, (void*)(ph_min) );
			this->deleteEdge(eh_min);
		}
		else{
cout<<"[deleting point] "<<n<<"/"<<tcount<<"  : "<<
ph_min->x<<","<<ph_min->y<<" -> "<<ph_min->epsilon<<": "<<endl;
			this->deletePoint(ph_min);
		}
	}

	/* free mem */
	fh_deleteheap(this->nodeHeap);	
	fh_deleteheap(this->edgeHeap);	

}


void Thinning::thinningAT8(unsigned int tcount)
{
	/* build point heap */
	for (unsigned int i = 0; i < this->triangulation->nodes.size(); i++) 
	{
		this->calculateSignificanceForHeap(this->triangulation->nodes[i]);
	}
	for(unsigned int n = 0; n < tcount; n++) 
	{
		/* get less significant node from point heap */
		Point2D* ph_min = (Point2D*)fh_extractmin(this->nodeHeap);
		ph_min->fh_el = NULL;
		vector<Point2D*> neighbors = ph_min->getNeighbors();
		/* delete point */
		this->triangulation->deletePoint(ph_min);
		/* recalculate significance of neighbors */
		for(unsigned int i = 0; i < neighbors.size(); i++){
			if(neighbors[i]->gridBound) continue;
			fh_delete(this->nodeHeap, (fibheap_el*)neighbors[i]->fh_el);
			neighbors[i]->fh_el = NULL;
			//this->calculateSignificanceForHeap(neighbors[i]);
			this->calculateSignificanceForHeap_ForAT8(neighbors[i]);
		}
/* ********************** DEBUG ***************************** */		
//cout<<"[deleting point] "<<n<<"/"<<tcount<<"  : "<<
//ph_min->x<<","<<ph_min->y<<" -> "<<ph_min->epsilon<<endl;
/* ********************** DEBUG ***************************** */
		
	}
	
	/* free mem */
	fh_deleteheap(this->nodeHeap);	
	fh_deleteheap(this->edgeHeap);	
}




void Thinning::deletePoint(Point2D* point)
{
	/* remove old edges from edde heap */
	this->updateHeaps(point->getNeighbors(),point, NULL);
}

void Thinning::deleteEdge(Edge* edge)
{
	Point2D* org = edge->org;
	Point2D* dest = edge->dest;
	/* remove points from node heap */
	fh_delete(this->nodeHeap, (fibheap_el*)(org->fh_el)); 
	org->fh_el=NULL;
	fh_delete(this->nodeHeap, (fibheap_el*)(dest->fh_el)); 
	dest->fh_el=NULL;
	/* remove edge from edge heap */
	fh_delete(this->edgeHeap, (fibheap_el*)(edge->fh_el)); 
	edge->fh_el=NULL;
	
	/* collect unique list of edge neighbor nodes */
	vector<Point2D*> neighbors;
	vector<Point2D*> orgNeighbors = org->getNeighbors();
	for (unsigned int i = 0; i < orgNeighbors.size(); i++) {
		if(orgNeighbors[i] != dest){
			neighbors.push_back(orgNeighbors[i]);
		}
	}
	vector<Point2D*> destNeighbors = dest->getNeighbors();
	for (unsigned int i = 0; i < destNeighbors.size(); i++) {
		if( find(neighbors.begin(), neighbors.end(), destNeighbors[i]) == neighbors.end() 
		&& destNeighbors[i] != org){
			neighbors.push_back(destNeighbors[i]);
		}
	}
	
	this->updateHeaps(neighbors,edge->org, edge);
}

void Thinning::updateHeaps(vector<Point2D*> neighbors, 
							Point2D* point, 
							Edge* edge)
{
	
	/* remove old edges from edge heap */
	for (unsigned int j = 0; j < neighbors.size(); j++) {

		if(neighbors[j]->gridBound) continue;
		fh_delete(this->nodeHeap, (fibheap_el*)(neighbors[j]->fh_el)); 
		neighbors[j]->fh_el=NULL;
		vector<Edge*> nEdges = neighbors[j]->getEdges();
		for (unsigned int i= 0; i < nEdges.size(); i++) {
			if( nEdges[i]->fh_el != NULL
			&& !nEdges[i]->org->gridBound 
			&& !nEdges[i]->dest->gridBound ){
				fh_delete(this->edgeHeap, (fibheap_el*)(nEdges[i]->fh_el)); 
				nEdges[i]->fh_el = NULL;
			}
		}
	}
		
	if(edge != NULL){
//		this->triangulation->deleteEdge(edge);
/* */
		Point2D* o = edge->org;
		Point2D* d = edge->dest;
		this->triangulation->deletePoint(o);
		this->triangulation->deletePoint(d);

	}
	else{
		this->triangulation->deletePoint(point);		
	}
	
	vector<Edge*> dirtryEdges;
	for (unsigned int j = 0; j < neighbors.size(); j++) {
		if(neighbors[j]->gridBound) continue;
		this->calculateSignificanceForHeap(neighbors[j]);
	}
	
	
	/* this is necessary to avoid conflicts on dirty flags of edges */
	for (unsigned int j = 0; j < neighbors.size(); j++) {
		vector<Edge*> nEdges = neighbors[j]->getEdges();
		for (unsigned int i= 0; i < nEdges.size(); i++) {
			if(find(dirtryEdges.begin(), dirtryEdges.end(), nEdges[i]) == dirtryEdges.end() ){
				dirtryEdges.push_back(nEdges[i]);
			}
		}
	}	
	
	/* calculate edge significance */
	for (unsigned int i= 0; i < dirtryEdges.size(); i++) 
	{
		this->calculateSignificanceForEdgeHeap(dirtryEdges[i]);
	}
}


void Thinning::calculateSignificanceForHeap(Point2D* point)
{
	if(!point->gridBound){
		/* calculate significance the this node */
		point->epsilon = this->calculatePointSignificance(point);
	    /* insert node into the fib heap and store pointer to the heap element */
	    point->fh_el = (void*)fh_insert( this->nodeHeap, (void*)(point) );
	    point->dirty = false;
	}
}

void Thinning::calculateSignificanceForHeap_ForAT8(Point2D* point)
{
	if(!point->gridBound){
		/* calculate significance the this node */
		//point->epsilon = this->calculatePointSignificance(point);
		/*if((int)(point->x)%3 ==0 && (int)(point->y)%3 ==0)
          point-> epsilon = 12000.*rand()/RAND_MAX;
        else
		  point->epsilon = 1.*rand()/RAND_MAX;*/ //random value
	    point->epsilon = this->calculatePointSignificance_ForAT8(point);
	    /* insert node into the fib heap and store pointer to the heap element */
	    point->fh_el = (void*)fh_insert( this->nodeHeap, (void*)(point) );
	    point->dirty = false;
	}
}


void Thinning::calculateSignificanceForEdgeHeap(Edge* edge)
{
	if(!edge->org->gridBound && !edge->dest->gridBound){
		/* calculate significance the this node */
		edge->epsilon = this->calculateEdgeSignificance(edge);
	    /* insert node into the fib heap and store pointer to the heap element */
	    edge->fh_el = (void*)fh_insert( this->edgeHeap, (void*)(edge) );
	}
}

int Thinning::exchange(unsigned int exchange_iterations)
{
	/* build thinned neighbor heap */
	for (unsigned int i = 0; i < this->triangulation->nodes.size(); i++) 
	{
		if(!this->triangulation->nodes[i]->thinned
		&& !this->triangulation->nodes[i]->gridBound)
		{
			calculatePairSignificance(this->triangulation->nodes[i]);		
		}
	}
	/* do exchange */
	int xcount = 0;
	while ( true && (/*exchange_iterations<0 || */ exchange_iterations > xcount)) {
		/* heap is empty */
		if(fh_min(this->thinnedNeighborHeap) == NULL) break;
		Point2D* point = (Point2D*)fh_extractmin(this->thinnedNeighborHeap);
		point->exchange_fh_el = NULL;
		if(point->dirty){
			calculatePairSignificance(point);
			continue;
		}
		
		Point2D* bestNeighbor = point->bestNeighbor; 
		double diff = point->bestNeighborDiff;
		if(point->thinned
		|| bestNeighbor == NULL
		|| !bestNeighbor->thinned ) continue;
		/* exchange */
		xcount++;
		
		set<Point2D*> neighbors;
		vector<Point2D*> nbs = point->getNeighbors();
		Point2D::addAll(&neighbors, nbs);
		vector<Triangle*> triangles = point->getTriangles();
		for(unsigned int i=0; i < triangles.size(); i++){
			if( triangles[i]->insideTriangle(bestNeighbor) ){
				bestNeighbor->thinnedTriangle = triangles[i];
				break;
			}
		}
		this->triangulation->insertPoint(bestNeighbor);
		Triangle::deleteTriangles(triangles);
		triangles.clear();
		nbs = point->getNeighbors();
		vector<Edge*> edges = this->triangulation->deletePoint(point);
		
		
		if(nbs.size() == 3){
			Triangle* t = Triangle::makeTriangle(nbs[0],nbs[1],nbs[2]);
			triangles.push_back(t);
		}
		else{
			triangles = this->triangulation->getTriangles(edges);		
		}
		
		point->exchanged = true;
		
		neighbors.insert(bestNeighbor);
		Point2D::addAll(&neighbors, bestNeighbor->getNeighbors());	
		for (set<Point2D*>::iterator i = neighbors.begin(); 
									 i != neighbors.end(); i++) {
			(*i)->dirty = true;
			calculatePairSignificance(*i);
		}
			
double nDiff = bestNeighbor->epsilon-point->epsilon;
cout<<"[Exchange "<<xcount<<"] "<<point->x<<","<<point->y<<
","<<point->epsilon<<" - "<<bestNeighbor->x<<","<<bestNeighbor->y<<
","<<bestNeighbor->epsilon<<" : "<<nDiff<<
" - "<<diff<<endl;
if(abs((int)(nDiff-diff)) > 1)
{
	//cout<<"FEHLER "<<nDiff<<" : "<<diff<<endl;
	if(nDiff<0 && true){
		for(unsigned int i=0; i < triangles.size(); i++){
			if( triangles[i]->insideTriangle(point) ){
				point->thinnedTriangle = triangles[i];
				break;
			}
		}
		this->triangulation->insertPoint(point);		
		this->triangulation->deletePoint(bestNeighbor);	
		if(bestNeighbor->exchange_fh_el != NULL){
			fh_delete(this->thinnedNeighborHeap, (fibheap_el*)bestNeighbor->exchange_fh_el);	
			bestNeighbor->exchange_fh_el = NULL;
		}
		point->exchanged = false;
		neighbors.insert(point);
		neighbors.erase(bestNeighbor);
		for (set<Point2D*>::iterator i = neighbors.begin(); 
									 i != neighbors.end(); i++) {
			calculatePairSignificance(*i);
			(*i)->dirty = true;
		}
		xcount--;
		//cout<<"BACK EXCHANGE "<<endl;
	}
}

		/* refactor triangle objects */
		Triangle::deleteTriangles(triangles);
	}		
		
	fh_deleteheap(this->thinnedNeighborHeap);	
	return xcount;

}

void Thinning::calculatePairSignificance(Point2D* point)
{
	/* reset point neighbor attributes */
	point->bestNeighbor = NULL;
	point->bestNeighborDiff = 0;

	if(point->exchange_fh_el != NULL){
		fh_delete(this->thinnedNeighborHeap, (fibheap_el*)point->exchange_fh_el);
		point->exchange_fh_el = NULL;
	}
	if(point->gridBound) return;

	point->dirty = false;
	/* get the next neighbors */
	vector<Point2D*> neighbors;
	vector<Triangle*> triangles = this->deletePointWithSigRecalc(point);
	
	int startX = max(0,point->x-Thinning::exchangeRadius);
	int endX = min(this->triangulation->xMax,point->x+Thinning::exchangeRadius);
	int startY = max(0,point->y-Thinning::exchangeRadius);
	int endY = min(this->triangulation->yMax,point->y+Thinning::exchangeRadius);
	for (int i = startX; i <= endX; i++) {
		int xOffset = i*this->triangulation->nbRows;
		for (int j = startY; j <= endY; j++) {
			Point2D* np = this->triangulation->nodes[xOffset+j];
			if(np->thinned && np != point){
				bool inside = false;
				// collect only thinned nodes inside fo cell hull 
				for (unsigned int k = 0; k < triangles.size(); ++k) {
					if( triangles[k]->insideTriangle(np->x, np->y) ){
						np->thinnedTriangle = triangles[k];
						inside = true;
						break;
					}
				}
				if(inside){
					neighbors.push_back(np);
				}
			}
		}
	}

	
	/* remove current point */
	Point2D* bestNeighbor = NULL;
	double bestDiff = 0.;
	double eps1 = point->epsilon;

	for( unsigned int i = 0; i < neighbors.size(); i++ ) {
		if( neighbors[i]->gridBound || neighbors[i]->exchanged ) continue;
		this->triangulation->insertPoint(neighbors[i]);
		vector<Triangle*> tt = this->deletePointWithSigRecalc(neighbors[i]);
		/* free triangle objects */
		Triangle::deleteTriangles(tt);
		
		double eps2 = neighbors[i]->epsilon;
		/* calculate significance */
		double diff = eps2 - eps1;
		if(bestDiff < diff){
			bestNeighbor = neighbors[i];
			bestDiff = diff;
			bestNeighbor->epsilon = eps2;
		}
		
	}

	if(bestNeighbor != NULL){
		point->bestNeighbor = bestNeighbor;
		point->bestNeighborDiff = bestDiff;
		point->exchange_fh_el = (void*)fh_insert( this->thinnedNeighborHeap, (void*)(point) );
	}
	
	/* set point thinned triangle */
	this->triangulation->insertPoint(point);
	// delete triangles
	Triangle::deleteTriangles(triangles);

}


//deletes point and sets the epsilon 
// DO NOT forget to free result triangles 
vector<Triangle*> Thinning::deletePointWithSigRecalc(Point2D* point)
{
	vector<Point2D*> neighbors = point->getNeighbors();
	vector<Point2D*> attachedPoints = this->triangulation->getAttachedPoints(neighbors);
	vector<Triangle*> oTriangles = point->getTriangles();
	vector<Edge*> edges = this->triangulation->deletePoint(point);
	vector<Triangle*> nTriangles;		
	if(neighbors.size() == 3){
		Triangle* t = Triangle::makeTriangle(neighbors[0],neighbors[1],neighbors[2]);
		nTriangles.push_back(t);
	}
	else{
		nTriangles = this->triangulation->getTriangles(edges);		
	}
	for (int i = nTriangles.size()-1; i >= 0; i--) {
		if( nTriangles[i]->insideTriangle(point) ){
			point->thinnedTriangle = nTriangles[i];
			break;
		}
	}
	
	double sig = calculateSignificance(oTriangles,nTriangles,attachedPoints,point);

	// refactor objects 
	Triangle::deleteTriangles(oTriangles);

	point->epsilon = sig;
	return nTriangles;
}


// IO-Functions 
void printSingleNode(ostream& out, Point2D* point, int option)
{
    out.precision(6);
    out.setf(ios::fixed, ios::floatfield);
    out.width(12);
    out << point->x << " ";
    out.width(12);
    out << point->y << " ";
    out.width(16);
    out.precision(8);
    if (option==0){
		if( point->epsilon == __DBL_MAX__ )
		{
			out << "INFINITE";
		}
      	else{
			out << point->epsilon;
      	}
    }
    else{
		out << point->f;
    }
}


void printSingleNodeEPS(ostream& out, Point2D* point/*, int option*/)
{
    //option = 1;
    
    out.precision(6);
    out.setf(ios::fixed, ios::floatfield);
    out << 0;
	out.width(12);
    out << point->x << " ";
    out.width(12);
    out << point->y << " ";
    out.width(16);
    out.precision(8);

	//Greyscale value writtenThree times because of RGB 
	out << point->f/255.0 << " ";     
	out << point->f/255.0 << " ";    
	out << point->f/255.0 << endl;    
}


void Thinning::printNodes(const string& filename, int option)
{
  vector<Point2D*> nodes = this->triangulation->getNodes();
  ofstream out(filename.c_str());
  //assert(out.is_open());
  out << nodes.size() << endl;
  out << this->triangulation->nbRows << endl;
  out << this->triangulation->nbCols << endl;
  for (unsigned int i=0; i<nodes.size(); i++) 
  {
  	printSingleNode(out,nodes[i],option);
  	out<<endl;
  }
  out.close();
}

void Thinning::printEdges(const string& filename, int option)
{
  vector<Edge*> edges = this->triangulation->getEdges();
  ofstream out(filename.c_str());
  //assert(out.is_open());
  out << edges.size() << endl;
  out << this->triangulation->nbRows << endl;
  out << this->triangulation->nbCols << endl;
  for (unsigned int i=0; i<edges.size(); i++) 
	{
  	printSingleNode(out,edges[i]->org,option);
  	printSingleNode(out,edges[i]->dest,option);
  	out<<endl;
  }
  out.close();
}


void Thinning::printTriangles(const string& filename, int option)
{
  vector<Triangle*> triangles = this->triangulation->getTriangles();
  ofstream out(filename.c_str());
  //assert(out.is_open());
  out << triangles.size() << endl;
  out << this->triangulation->nbRows << endl;
  out << this->triangulation->nbCols << endl;
  for (unsigned int i=0; i<triangles.size(); i++) 
	{
  	printSingleNode(out,triangles[i]->point1,option);
  	printSingleNode(out,triangles[i]->point2,option);
	  printSingleNode(out,triangles[i]->point3,option);
	  out<<endl;
  }
  out.close();
}


void Thinning::printTriangulationOff(const string& filename)
{
  cout << "[printTriangulationOff]" << endl;
  ofstream out(filename.c_str());
	
  vector<Triangle*> triangles = this->triangulation->getTriangles();
  vector<Point2D*> nodes = this->triangulation->getNodes();
	
  out << "OFF" << endl; 	
  out << nodes.size() << " " << triangles.size() << " " << "0"<< endl;
	
  M3Matrix Indices(this->triangulation->nbRows,this->triangulation->nbCols);
	
  //Write the nodes : coordinates and value
	for (unsigned int i=0; i<nodes.size(); i++) 
	{
		out << nodes[i]->x << " ";
		out << nodes[i]->y << " ";
		out << nodes[i]->f << " ";
		out << endl;
		
		Indices[int(nodes[i]->x)][int(nodes[i]->y)] = i;
	}
	
	//Write the nodes : coordinates and value
	for (unsigned int i=0; i<triangles.size(); i++) 
	{
    out << "3 ";
		int x1 = (int)triangles[i]->point1->x;
		int y1 = (int)triangles[i]->point1->y;
		int x2 = (int)triangles[i]->point2->x;
		int y2 = (int)triangles[i]->point2->y;
		int x3 = (int)triangles[i]->point3->x;
		int y3 = (int)triangles[i]->point3->y;
		out << Indices[x1][y1] << " ";
		out << Indices[x2][y2] << " ";
		out << Indices[x3][y3] << " ";
		out << endl;
	}
	
	out.close();
}


void Thinning::printTriangulationEps(const string& filename)
{
  cout << "[printTriangulationEps]" << endl;
	 
  vector<Edge*> edges = this->triangulation->getEdges();
  ofstream out(filename.c_str());
  //assert(out.is_open());
  
  out << "%!PS-Adobe-3.0 EPSF-3.0" << endl;
	out << "%%BoundingBox: " << 0 << " " << 0 << " " << this->triangulation->nbRows*2  << " "<< this->triangulation->nbCols*2  << endl;

  out << "0.01 setlinewidth" << endl;
  out << "0 0 0 setrgbcolor" << endl;
  //int nbrows = this->triangulation->nbRows;
  int nbcols = this->triangulation->nbCols;
	
  for (unsigned int i=0; i<edges.size(); i++) 
  {	  
  	double edge_AD = fabs((double)(edges[i]->org->f - edges[i]->dest->f));
    if(edge_AD>40.)
    {
      out << "1 0 0 setrgbcolor" << endl;
    }
    else
    {
      if(edge_AD<41.)
      {
        out << "0 1 0 setrgbcolor" << endl;
      }
      else
      {
        out << "0 0 0 setrgbcolor" << endl;
      }
    }
    
    out << (edges[i]->org->x)*2 +1<< " ";
    out << (nbcols-edges[i]->org->y-1)*2 +1 << " ";
	  
	  out << " moveto ";
	  out << (edges[i]->dest->x)*2 +1<< " ";
	  out << (nbcols-edges[i]->dest->y-1)*2 +1 << " ";
    out << " lineto stroke";
	  out << endl;
  }

  out.close();
}


void Thinning::printEps(const string& filename)
{
  vector<Triangle*> triangles = this->triangulation->getTriangles();
  ofstream out(filename.c_str());
  
  out << "%!PS-Adobe-3.0 EPSF-3.0" << endl;
  out << "%%BoundingBox: " << 0 << " " << 0 << " " << this->triangulation->nbCols*2  << " "<< this->triangulation->nbRows*2 << endl;

  out << "<<" << endl;
  out << "/ShadingType 4" << endl;
  out << "/ColorSpace [/DeviceRGB]" << endl;
  out << "/DataSource" << endl;
  out << "[" << endl;

  for (unsigned int i=0; i<triangles.size(); i++) 
  {
  	printSingleNodeEPS(out,triangles[i]->point1);
  	printSingleNodeEPS(out,triangles[i]->point2);
	printSingleNodeEPS(out,triangles[i]->point3);
	out<<endl;
  }

  out << "]" << endl;
  out << ">>" << endl;
  out << "shfill" << endl;

  out.close();
}


void Thinning::printThinnedNodes(const string& filename, int option)
{
  vector<Point2D*> nodes = this->triangulation->getThinnedNodes();
  ofstream out(filename.c_str());
  //assert(out.is_open());
  out << nodes.size() << endl;
  out << this->triangulation->nbRows << endl;
  out << this->triangulation->nbCols << endl;
  for (unsigned int i=0; i<nodes.size(); i++) 
	{
  	printSingleNode(out,nodes[i],option);
  	out<<endl;
  }
  out.close();
}




#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <set>
#include <map>
#include <string>

#include "Triangulation.h"
#include "Point2D.h"
#include "Edge.h"
#include "Triangle.h"
#include "Test.h"

using namespace std;

set<Triangulation*> Triangulation::pool;
Edge* Triangulation::dummyEdge = new Edge();

Triangulation::Triangulation(){}

Triangulation::Triangulation(PointGrid* g)
{
  this->nbRows = g->nbRows;
  this->nbCols = g->nbCols;
  this->xMax = g->xMax;
  this->yMax = g->yMax;
  this->xMin = g->xMin;
  this->yMin = g->yMin;
  this->nodes = g->nodes;
  this->triangulate();
}


Triangulation::~Triangulation()
{
}

Triangulation* Triangulation::makeTriangulation(vector<Point2D*> n)
{
	if(Triangulation::pool.size() > 0)
	{
		Triangulation* t = *(Triangulation::pool.begin());
		Triangulation::pool.erase(t);
		t->nodes = n;
		t->entry = NULL;
		t->triangulate();
		return t;
	}
	else
	{
		Triangulation* t = new Triangulation();
		t->nodes = n;
		t->triangulate();
		return t;	
	}
	return NULL;
}


void Triangulation::deleteTriangulation()
{
  // loop over edges
  Edge* e = this->entry;
  if(e != NULL)
  {
    vector<Edge*> edges = this->getEdges();
	for (unsigned int i = 0; i < edges.size(); i++ ) 
	{
	  edges[i]->refactorEdge();			
	}	
  }

  for (int i = this->nodes.size()-1; i >= 0; i--) 
  {
	this->nodes[i]->deletePoint2D();
  }
  this->nodes.clear();
	
  // refactor triangulation 
  if(Triangulation::pool.size() < TRIANGULATION_POOL_SIZE)
  {
    Triangulation::pool.insert(this);
  }
  else
  {
    // free edge 
	delete this;	
  }	
}
	
void Triangulation::deleteTriangulation(vector<Triangle*> triangles)
{
	// free triangle objects
	Triangle::deleteTriangles(triangles);
	
	this->deleteTriangulation();
}


void Triangulation::triangulate()
{
	this->entry = NULL;
	// sort nodes 
	sort( this->nodes.begin(), this->nodes.end(), positionComparator);
	this->divide(0,this->nodes.size()-1, &dummyEdge, &dummyEdge);
}

inline bool Triangulation::positionComparator(const void* a, const void* b)
{
	Point2D* pa = (Point2D*)a;
	Point2D* pb = (Point2D*)b;
	unsigned int ax = pa->x;
	unsigned int ay = pa->y;
	unsigned int bx = pb->x;
	unsigned int by = pb->y;
	
	if ( ax < bx ) return true;
	if ( ax > bx ) return  false;

	if ( ay < by ) return true;
	return  false;
}



/* ************************************* PRIVATE ***************************** */


void Triangulation::divide(int l, int r, Edge** le, Edge** re)
{
	int n,split;
	double c_p;

	Point2D* point_l = this->nodes[l];
	Point2D* point_r = this->nodes[r];
	Point2D* point_ln = this->nodes[l+1];
	
	n = r - l + 1;
	
	if (n == 2) 
	{
		// Bottom of the recursion. Make an edge 
		Edge* e = Edge::makeEdge(point_l, point_r);
		this->addEdge(e);
		
		*le = e;
		*re = e;
		return;
	} 
	else if (n == 3) 
	{
		// Bottom of the recursion. Make a triangle or two edges 
		Edge* a = Edge::makeEdge(point_l, point_ln);
		this->addEdge(a);
	    Edge* b = Edge::makeEdge(point_ln, point_r);
	    this->addEdge(b);
	    a->splice(b, point_ln);
	    c_p = point_l->crossProduct(point_ln, point_r);
	    
	    if (c_p > 0.0)
		{
	      // Make a triangle 
	      Edge* c = a->join(point_l, b, point_r, RIGHT);
	      this->addEdge(c);
		  *le = a;
		  *re = b;
		  return;
	    } 
        else if (c_p < 0.0) 
        {
	      // Make a triangle 
	      Edge* c = a->join(point_l, b, point_r, LEFT);
	      this->addEdge(c);
		  *le = c;
		  *re = c;
		  return;
	    } 
		else 
		{
	      // Points are collinear,  no triangle 
		  *le = a;
		  *re = b;
		  return;
	    }
	} 
	//  if (n  > 3)  
	else{
	    // Continue to divide 

	    // Calculate the split point 
	    split = (l + r) / 2;
	  
	    // Divide 
	    Edge *l_ccw_l, *r_cw_l, *l_ccw_r, *r_cw_r, *l_tangent;
	    this->divide(l, split, &l_ccw_l, &r_cw_l);
	    this->divide(split+1, r, &l_ccw_r, &r_cw_r);

	    // Merge 
	    l_tangent = this->merge(r_cw_l, this->nodes[split], l_ccw_r, this->nodes[split+1]);
	    
	    // The lower tangent added by merge may have invalidated 
	    //   l_ccw_l or r_cw_r. Update them if necessary. 
	    if ( l_tangent->org == point_l ){
	    	l_ccw_l = l_tangent;
	    }
	    if ( l_tangent->dest == point_r ){
	      r_cw_r = l_tangent;
	    }
	
	    // Update edge refs to be passed back  
	    *le = l_ccw_l;
	    *re = r_cw_r;

	    return;
	} 
}


Edge** Triangulation::divide(int l, int r)
{
	int n;
	int split;
	// Edge* l_ccw_l, r_cw_l, l_ccw_r, r_cw_r, l_tangent;
	// Edge* a, b, c;
	double c_p;

	Edge** ret = (Edge**)malloc(2*sizeof(Edge*));
	
	Point2D* point_l = this->nodes[l];
	Point2D* point_r = this->nodes[r];
	Point2D* point_ln = this->nodes[l+1];
	
	n = r - l + 1;
	
	if (n == 2) 
	{
	  // Bottom of the recursion. Make an edge 
	  Edge* e = Edge::makeEdge(point_l, point_r);
	  this->addEdge(e);
	  ret[0] = e;
	  ret[1] = e;
	  return ret;
	} 
	else if (n == 3) {
		// Bottom of the recursion. Make a triangle or two edges 
		Edge* a = Edge::makeEdge(point_l, point_ln);
		this->addEdge(a);
	    Edge* b = Edge::makeEdge(point_ln, point_r);
	    this->addEdge(b);
	    a->splice(b, point_ln);
	    c_p = point_l->crossProduct(point_ln, point_r);
	    
	    if (c_p > 0.0){
	    	// Make a triangle 
	    	Edge* c = a->join(point_l, b, point_r, RIGHT);
	    	this->addEdge(c);
			ret[0] = a;
			ret[1] = b;
			return ret;
	    } else if (c_p < 0.0) {
	    	// Make a triangle 
	    	Edge* c = a->join(point_l, b, point_r, LEFT);
	    	this->addEdge(c);
			ret[0] = c;
			ret[1] = c;
			return ret;
	    } else {
	    	// Points are collinear,  no triangle  
			ret[0] = a;
			ret[1] = b;
			return ret;
	    }
	} else if (n  > 3) {
	    // Continue to divide 

	    // Calculate the split point 
	    split = (l + r) / 2;
	  
	    // Divide 
	    Edge** e1 = this->divide(l, split);
	    Edge* l_ccw_l = e1[0];
	    Edge* r_cw_l = e1[1];
	    free(e1);
	    e1 = NULL;
	    
	    Edge** e2 = this->divide(split+1, r);
	    Edge* l_ccw_r = e2[0];
	    Edge* r_cw_r = e2[1];
	    free(e2);
	    e2 = NULL;

	
	    // Merge 
	    Edge* l_tangent = this->merge(r_cw_l, this->nodes[split], l_ccw_r, this->nodes[split+1]);
	    
	    // The lower tangent added by merge may have invalidated 
	    //   l_ccw_l or r_cw_r. Update them if necessary. 
	    if ( l_tangent->org == point_l ){
	    	l_ccw_l = l_tangent;
	    }
	    if ( l_tangent->dest == point_r ){
	      r_cw_r = l_tangent;
	    }
	
	    // Update edge refs to be passed back  
	    ret[0] = l_ccw_l;
	    ret[1] = r_cw_r;
	    return ret;
	} 
	return NULL;
}



Edge* Triangulation::merge(Edge* r_cw_l, Point2D* s, Edge* l_ccw_r, Point2D* u)
{
	double cot_l_cand = 0.;
	double cot_r_cand = 0.;
	
	// Create first cross edge by joining lower common tangent 
	Edge* l_lower, *r_lower;
	Point2D *org_r_lower, *org_l_lower;

	getLowerTangent(r_cw_l, s, l_ccw_r, u, 
						&org_l_lower, 
						&org_r_lower, 
						&l_lower, 
						&r_lower);

	Edge* base = l_lower->join(org_l_lower, r_lower, org_r_lower, RIGHT);
	this->addEdge(base);
	Point2D* org_base = org_l_lower;
	Point2D* dest_base = org_r_lower;
	  
	// Need to return lower tangent. 
	Edge* l_tangent = base;

	// Main merge loop 
	do{
	    // Initialise l_cand and r_cand 
	    Edge* l_cand = base->getNext(org_base);
	    Edge* r_cand = base->getPrev(dest_base);
	    Point2D* dest_l_cand = l_cand->getOtherPoint(org_base);
	    Point2D* dest_r_cand = r_cand->getOtherPoint(dest_base);

	    // Above tests
	    double c_p_l_cand = dest_l_cand->crossProduct(org_base, dest_l_cand, dest_base); 
	    double c_p_r_cand = dest_r_cand->crossProduct(org_base, dest_r_cand, dest_base);
	    bool above_l_cand = c_p_l_cand > 0.0;
	    bool above_r_cand = c_p_r_cand > 0.0;
	    if (!above_l_cand && !above_r_cand)
	      break;        /* Finished. */

	    /* Advance l_cand ccw,  deleting the old l_cand edge,  until the 
	       "in_circle" test fails. */
	    if (above_l_cand)
	    {
	      double d_p_l_cand = dest_l_cand->dotProduct(org_base, dest_l_cand, dest_base);
	      cot_l_cand = d_p_l_cand / c_p_l_cand;

	      do 
	      {
	        Edge* next = l_cand->getNext(org_base);
	        Point2D* dest_next = next->getOtherPoint(org_base);
	        double c_p_next = dest_next->crossProduct(org_base, dest_next, dest_base);
	        bool above_next = c_p_next > 0.0;

	        if (!above_next) 
	          break;    /* Finished. */

	        double d_p_next = dest_next->dotProduct(org_base, dest_next, dest_base);
	        double cot_next = d_p_next / c_p_next;

	        //TODO
	        if (cot_next > cot_l_cand)
	          break;

	        //TODO this->edges.premove(l_cand);
	        this->removeEdge(l_cand);
	        
	        l_cand = next;
	        cot_l_cand = cot_next;
	  
	      } while (true);
	    }

	    /* Now do the symmetrical for r_cand */
	    if (above_r_cand)
	    {
	      double d_p_r_cand = dest_r_cand->dotProduct(org_base, dest_r_cand, dest_base);
	      cot_r_cand = d_p_r_cand / c_p_r_cand;

	      do
	      {
	        Edge* prev = r_cand->getPrev(dest_base);
	        Point2D* dest_prev = prev->getOtherPoint(dest_base);
	        double c_p_prev = dest_prev->crossProduct(org_base, dest_prev, dest_base);
	        bool above_prev = c_p_prev > 0.0;

	        if (!above_prev) 
	          break;    /* Finished. */

	        double d_p_prev = dest_prev->dotProduct(org_base, dest_prev, dest_base);
	        double cot_prev = d_p_prev / c_p_prev;

	        //TODO
	        /* das war stric > nachvollziehen!!! */
	        if (cot_prev > cot_r_cand)
	          break;    /* Finished. */

	        //TODO this->edges.remove(r_cand);
	        this->removeEdge(r_cand);
	        r_cand = prev;
	        cot_r_cand = cot_prev;

	      } while (true);
	    }

	    // Now add a cross edge from base to either l_cand or r_cand. 
	    // If both are valid choose on the basis of the in_circle test . 
	    // Advance base and  whichever candidate was chosen.
	    dest_l_cand = l_cand->getOtherPoint(org_base);
	    dest_r_cand = r_cand->getOtherPoint(dest_base);

	    // Now add a cross edge from base to either l_cand or r_cand. 
	    // If both are valid choose on the basis of the in_circle test . 
	    // Advance base and  whichever candidate was chosen.	     
	    Edge* oldBase = base;
	    if (!above_l_cand 
	    || (above_l_cand && above_r_cand && cot_r_cand < cot_l_cand)){
	    	/* Connect to the right */
	    	base = base->join(org_base, r_cand, dest_r_cand, RIGHT);
	    	this->addEdge(base);
	    	dest_base = dest_r_cand;
	    	r_cand->assumeUniquness(0);
	    } 
	    else {
	    	if(!above_l_cand 
	    	|| !above_r_cand 
	    	|| cot_r_cand != cot_l_cand){
		    	/* Connect to the left */
		    	base = l_cand->join(dest_l_cand, base, dest_base, RIGHT);
		    	this->addEdge(base);
		      	org_base = dest_l_cand;
		      	l_cand->assumeUniquness(0);
	    	}
	    	else{
		    	// Connect to the right 
		    	base = base->join(org_base, r_cand, dest_r_cand, RIGHT);
		    	this->addEdge(base);
		    	dest_base = dest_r_cand;
		    	r_cand->assumeUniquness(0);
	    	}
	    }
    	
	    // assume uniqueness
	    // check oldbase
	    // check score and on-circle
	    oldBase->assumeUniquness(0);
	    	    	    
	} while (true);
	return l_tangent;
}


//  Determines the lower tangent of two triangulations. 
void Triangulation::getLowerTangent(Edge* r_cw_l, 
						Point2D* s, 
						Edge* l_ccw_r, 
						Point2D* u,
						Point2D **o_l,
						Point2D **o_r,
						Edge **l,
						Edge **r)
{
	*l = r_cw_l;
	*r = l_ccw_r;
	*o_l = s;
	Point2D* d_l = (*l)->getOtherPoint(s);
	*o_r = u;
	Point2D* d_r = (*r)->getOtherPoint(u);
	bool finished = false;
  
	while (!finished)
	{
	    if ((*o_l)->crossProduct(d_l, *o_r) > 0.0) {
			(*l) = (*l)->getPrev(d_l);
			*o_l = d_l;
			d_l = (*l)->getOtherPoint(*o_l);
	    } else if ((*o_r)->crossProduct(d_r, *o_l) < 0.0) {
			*r = (*r)->getNext(d_r);
			*o_r = d_r;
			d_r = (*r)->getOtherPoint(*o_r);
	    } else{
	      finished = true;
	    }
	}
}

	
vector<Point2D*> Triangulation::getAttachedPoints(unsigned int maxX, 
												  unsigned int minX, 
												  unsigned int maxY,
												  unsigned int minY)
{
	vector<Point2D*> ret;
	int rowOffsetInit = this->nbRows;
	for (unsigned int i = minX; i <= maxX; i++) 
	{
	  int rowOffset = rowOffsetInit*i;
	  for (unsigned int j = minY; j <= maxY; j++) 
	  {
	    Point2D* point = this->nodes[rowOffset+j];
		if(point->thinned)
		{
		  ret.push_back(point);
		}
	  }
	}
	return ret;
}	
	
vector<Point2D*> Triangulation::getAttachedPoints(vector<Point2D*> neighbors)
{
	unsigned int maxX = 0, minX = this->xMax, maxY = 0, minY = this->yMax;
	for (unsigned int i = 0; i < neighbors.size(); i++) {
		unsigned int currX = neighbors[i]->x;
		unsigned int currY = neighbors[i]->y;
		if(currX < minX){
			minX = currX;
		}
		if(currY < minY){
			minY = currY;
		}
		if(currX > maxX){
			maxX = currX;
		}
		if(currY > maxY){
			maxY = currY;
		}
	}
	
	return this->getAttachedPoints(maxX,minX,maxY,minY);
}

vector<Point2D*> Triangulation::getAttachedPoints(Point2D* point)
{
  return this->getAttachedPoints(point->getNeighbors());
}
	
	
vector<Edge*> Triangulation::getEdges()
{
  vector<Edge*> edges;
  Edge* e = this->entry;
	
  // loop over edges 
  while( true )
  {
    edges.push_back(e);
	if( e == e->succ )
	{
	  break;	
	}
	else
	{
	  e = e->succ;
	}
  }
  return edges;
}


vector<Triangle*> Triangulation::getTriangles(Edge* edge)
{
  return this->getTriangles(edge->getEdges());
}

vector<Triangle*> Triangulation::getTriangles()
{
  return this->getTriangles(this->getEdges());
}

vector<Triangle*> Triangulation::getTriangles(vector<Edge*> edges)
{
	vector<Triangle*> triangles;
	// to mark edge as processed 
	set<Edge*> processedEdges;
	// loop over edges 
	for( unsigned int i = 0; i < edges.size(); i++ )
    {
		Point2D* org = edges[i]->org;
		
		bool collected = false;
		Edge* ne = edges[i]->dnext;
		// check if this edges can be precessed 
		if( processedEdges.find(ne) == processedEdges.end()
  	    &&  processedEdges.find(edges[i]->oprev) == processedEdges.end() ){
			Triangle* t1 = processTrianglesByEdge(ne, org);
			if(t1 != NULL){
				triangles.push_back(t1);	
				collected = true;
			}
  	    }

		Edge* pe = edges[i]->dprev;
		// check if this edges can be precessed 
		if(  (!collected || pe != ne)
		&&  processedEdges.find(pe) == processedEdges.end()
  	    &&  processedEdges.find(edges[i]->onext) == processedEdges.end() ){
	  	    Triangle* t1 = processTrianglesByEdge(pe,org);
			if(t1 != NULL){
				triangles.push_back(t1);	
			}
  	    }
		// mark edge as processed 
		processedEdges.insert(edges[i]);
	}
	return triangles;
}

	
Triangle* Triangulation::processTrianglesByEdge(Edge* edge, Point2D* point)
{
  if(edge->checkTriangle(point))
  {
	return Triangle::makeTriangle(edge, point);
  }
  else
  {
	return NULL;	
  }
}	


void Triangulation::addEdge(Edge* e)
{
  if( this->entry == NULL )
  {
		this->entry = e;
		e->pred = e ;
		e->succ = e ;
  }
  else
  {
	this->entry->pred = e;
	e->pred = e;
	e->succ = this->entry;
	this->entry = e;	
  }
}


void Triangulation::removeEdgeFromLinkedList(Edge* e)
{
  // first edge 
  if( this->entry == e )
  {
	this->entry = e->succ;
  }
  // last edge 
  else if(e->succ == e)
  {
	Edge* pred = e->pred;
	pred->succ = pred;
  }
  // otherwise 
  else
  {
	Edge* succ = e->succ;
	Edge* pred = e->pred;
	pred->succ = succ;
	succ->pred = pred;
  }
}


void Triangulation::removeEdge(Edge* e)
{
  this->removeEdgeFromLinkedList(e);
  e->deleteEdge();
}


vector<Edge*> Triangulation::integrateSubTriangulation(vector<Edge*> cell, 
										vector<Triangle*> cellTriangles,
										vector<Point2D*> neighbors)
{
    vector<Edge*> newEdges;
	vector<Point2D*> neighborsCopy = Point2D::clonePoint2D(neighbors);

	// make a virtual triangulation 
	Triangulation* tri = Triangulation::makeTriangulation(neighborsCopy);
	vector<Edge*> edges = tri->getEdges();
	// map old cell edges to the cell edges in the copy triangulation 
	map<Edge*, Edge*> edgeMap;

	// loop over edges 
	for (unsigned int i = 0; i < cell.size(); i++) 
	{
	  cell[i]->org->entry = cell[i];
	  cell[i]->dest->entry = cell[i];

		bool do_continue = false;
		for (unsigned int j = 0; j < edges.size(); j++) {
			Point2D* o = edges[j]->org->clone;
			Point2D* d = edges[j]->dest->clone;
			if(cell[i]->isVertex(o)
			&& cell[i]->isVertex(d)){
				edgeMap.insert(make_pair(edges[j],cell[i]));
				do_continue = true;
				break;
			}
		}
		if(do_continue) continue;
	
		// edge not found: we must swap an edge 
		for (unsigned int j = 0; j < edges.size(); j++) {	
			Point2D* p1 = edges[j]->dnext->getOtherPoint(edges[j]->dest)->clone;
			Point2D* p2 = edges[j]->onext->getOtherPoint(edges[j]->org)->clone;
			if(cell[i]->isVertex(p1)
			&& cell[i]->isVertex(p2)){
				edges[j]->swap();
				edgeMap.insert(make_pair(edges[j],cell[i]));
				break;
			}
		}
		
	}

	map<Edge*,bool> edgeCollected;
			
	// insert new edges 
	for (unsigned int j = 0; j < edges.size(); j++) 
	{
		// check if this is an old cell edge 
		map<Edge*, Edge*>::iterator iter = edgeMap.find(edges[j]);
		// edge found 
		if(iter != edgeMap.end()) continue;
		
		// edge not found
		Point2D* o = edges[j]->org->clone;
		Point2D* d = edges[j]->dest->clone;
		
		// check if this is an edge which is outside of the old cell: 
		// edge is on the convex hull of the new cell but the old cell is not convex
		
		bool insideHull = false;
		
		// compute edge center point 
		double cx = (edges[j]->org->x+edges[j]->dest->x) / 2.0;
		double cy = (edges[j]->org->y+edges[j]->dest->y) / 2.0;
		for (unsigned int i = 0; i < cellTriangles.size(); i++) 
		{
			if(cellTriangles[i]->insideTriangle(cx,cy))
			{
				insideHull = true;
				break;
			}
		}
		if(!insideHull) continue;
		
		edges[j]->org = o;
		edges[j]->dest = d;
		this->addEdge(edges[j]);
		// mark this edge as collected 
		edgeCollected.insert(make_pair(edges[j],true));
		
		o->entry = edges[j];
		d->entry = edges[j];
		newEdges.push_back(edges[j]);
		/* update next/prev */
		iter = edgeMap.find(edges[j]->onext);
		if(iter != edgeMap.end()){
			Edge* onext = iter->second;
			edges[j]->onext = onext;
			onext->setPrev(o, edges[j]);
		}
		iter = edgeMap.find(edges[j]->dnext);
		if(iter != edgeMap.end()){
			Edge* dnext = iter->second;
			edges[j]->dnext = dnext;
			dnext->setPrev(d, edges[j]);
		}
		iter = edgeMap.find(edges[j]->oprev);
		if(iter != edgeMap.end()){
			Edge* oprev = iter->second;
			edges[j]->oprev = oprev;
			oprev->setNext(o,edges[j]);
		}
		iter = edgeMap.find(edges[j]->dprev);
		if(iter != edgeMap.end()){
			Edge* dprev = iter->second;
			edges[j]->dprev = dprev;
			dprev->setNext(d, edges[j]);
		}
	}
	
	// free mem 
	for (unsigned int i = 0; i < edges.size(); i++) 
    {
	  // if edge not used by main triangulation: free/refactor edge 
	  if( edgeCollected.find(edges[i]) == edgeCollected.end() )
      {
		edges[i]->refactorEdge();				
	  }
	}

	edgeCollected.clear();
	// to avoid trying free edges
	tri->entry = NULL;
	tri->deleteTriangulation();
			
	return newEdges;
}


vector<Edge*> Triangulation::insertPoint(Point2D* point)
{
	Triangle* thinned = point->thinnedTriangle;
	Edge* thinnedEdge1;
	Edge* thinnedEdge2;
	Edge* thinnedEdge3;
	thinned->getEdges(&thinnedEdge1, &thinnedEdge2, &thinnedEdge3);
	
	set<Point2D*> neighbors;
	set<Edge*> cell;
	cell.insert(thinnedEdge1);
	cell.insert(thinnedEdge2);
	cell.insert(thinnedEdge3);
	vector<Triangle*> cellTriangles;
	cellTriangles.push_back(Triangle::makeTriangle(thinned));
	neighbors.insert(thinned->point1);
	neighbors.insert(thinned->point2);
	neighbors.insert(thinned->point3);

	// check if point is on an edge 
	Triangle* adj = NULL;
	Edge* adjEdge = NULL;

	// edge 1 
	if( point->crossProduct(thinned->point1,thinned->point2) == 0 )
	{
		adj = thinned->getNeighbourByEdge(thinnedEdge1); 
		if(adj != NULL){
			adjEdge = thinnedEdge1;
		}
		else{
			if(thinnedEdge1->isConvexHull(this->xMax,
											this->xMin,
											this->yMax,
											this->yMin)){
				Edge *ne1, *ne2;
				this->splitEdge(point,thinnedEdge1,&ne1, &ne2);
				this->addEdge(ne1);
				this->addEdge(ne2);
				cell.insert(ne1);
				cell.insert(ne2);
				cell.erase(thinnedEdge1);
			}
		}
	}			
	// edge 2 
	else if( point->crossProduct(thinned->point2,thinned->point3) == 0 ){
		adj = thinned->getNeighbourByEdge(thinnedEdge2); 
		if(adj != NULL){
			adjEdge = thinnedEdge2;
		}
		else{
			if(thinnedEdge2->isConvexHull(this->xMax,
										  this->xMin,
										  this->yMax,
										  this->yMin)){
				Edge *ne1, *ne2;
				this->splitEdge(point, thinnedEdge2, &ne1, &ne2);
				this->addEdge(ne1);
				this->addEdge(ne2);
				cell.insert(ne1);
				cell.insert(ne2);
				cell.erase(thinnedEdge2);
			}
		}
	}			
	// edge 3 
	else if( point->crossProduct(thinned->point1,thinned->point3) == 0 )
	{
		adj = thinned->getNeighbourByEdge(thinnedEdge3); 
		if(adj != NULL)
		{
			adjEdge = thinnedEdge3;
		}
		else
		{
			if(thinnedEdge3->isConvexHull(this->xMax,this->xMin,
										  this->yMax,this->yMin))
			{
				Edge *ne1, *ne2;
				this->splitEdge(point, thinnedEdge3, &ne1, &ne2);
				this->addEdge(ne1);
				this->addEdge(ne2);
				cell.insert(ne1);
				cell.insert(ne2);
				cell.erase(thinnedEdge3);
			}
		}
	}			
	
	if(adj != NULL){
		cellTriangles.push_back(adj);
		Edge *edge1, *edge2, *edge3;
		adj->getEdges(&edge1, &edge2, &edge3);
		cell.insert(edge1);
		cell.insert(edge2);
		cell.insert(edge3);
		neighbors.insert(adj->point1);
		neighbors.insert(adj->point2);
		neighbors.insert(adj->point3);
		cell.erase(adjEdge);
		adjEdge->deleteEdgeInsideCell(cell);
		this->removeEdgeFromLinkedList(adjEdge);
	}
	
	neighbors.insert(point);
	
	/* copy cell 2 vector */
	vector<Edge*> v_cell;
	set<Edge*>::iterator cIter = cell.begin();
	for(cIter = cell.begin(); cIter != cell.end(); cIter++){
		v_cell.push_back(*cIter);
	}
	/* copy neighbors 2 vector */
	vector<Point2D*> v_neighbors;
	set<Point2D*>::iterator nIter = neighbors.begin();
	for(nIter = neighbors.begin(); nIter != neighbors.end(); nIter++)
    {
		v_neighbors.push_back(*nIter);
	}

	
	vector<Edge*> ret = integrateSubTriangulation(v_cell, cellTriangles, v_neighbors);
	
	// check edge for Delaunay 
	for (unsigned int i = 0; i < v_cell.size(); i++)
	{
		if(v_cell[i]->checkTriangle(point)){
			Triangle *t = Triangle::makeTriangle(v_cell[i],point);
			v_cell[i]->checkEdge(t,0);
			t->deleteTriangle();
		}
	}
	
	// free triangle objects 
	Triangle::deleteTriangles(cellTriangles);
	
	point->thinned = false;

	return ret;
}


void Triangulation::splitEdge(Point2D* point, Edge* edge, Edge **e1, Edge **e2)
{
	this->removeEdgeFromLinkedList(edge);
	Edge *neo = Edge::makeEdge(edge->org,point);
	neo->adjustQuadStructure(edge->org,edge->onext,edge->oprev,point,neo,neo);
	Edge *ned = Edge::makeEdge(edge->dest,point);
	point->entry = ned;
	ned->adjustQuadStructure(edge->dest,edge->dnext,edge->dprev,point,neo,neo);
	if( edge->org->entry == edge ){
		edge->org->entry = neo;
	}
	if( edge->dest->entry == edge ){
		edge->dest->entry = ned;
	}
	*e1 = neo;
	*e2 = ned;
}


void Triangulation::removePoint(Point2D* point, vector<Edge*> *cell)
{
	vector<Edge*> edges2delete;
	vector<Edge*> nEdges = point->getEdges();

	set<Edge*> edgeCollected;

	// first remove point 
	for ( unsigned int i = 0;	i < nEdges.size(); i++ ) 
	{
		if( edgeCollected.find(nEdges[i]) != edgeCollected.end() ) continue;
		if( nEdges[i]->isConvexHull(this->xMax,
									this->xMin,
									this->yMax,
									this->yMin) ){
			/* get the other hull edge */
			Edge* e2 = NULL;
			for ( unsigned int i2 = 0;	i2 < nEdges.size(); i2++ ) {
				if(i2 != i
				&& edgeCollected.find(nEdges[i2]) == edgeCollected.end()
				&& nEdges[i2]->isConvexHull(this->xMax,
											this->xMin,
											this->yMax,
											this->yMin)
				&& ( nEdges[i]->isVertex(nEdges[i2]->org) || nEdges[i]->isVertex(nEdges[i2]->dest) ) ){
					e2 = nEdges[i2];
					break;
				}
			}
			Point2D* other = e2->getOtherPoint(point);
			Edge* next = e2->getNext(other);
			Edge* prev = e2->getPrev(other);
			if( nEdges[i]->org == point ){
				nEdges[i]->org = other;
				nEdges[i]->onext = next;
				nEdges[i]->oprev = prev;
			}
			else{
				nEdges[i]->dest = other;
				nEdges[i]->dnext = next;
				nEdges[i]->dprev = prev;
			}
			next->setPrev(other, nEdges[i]);
			prev->setNext(other, nEdges[i]);
			if( other->entry == e2 ){
				other->entry = nEdges[i];
			}
			vector<Edge*>::iterator pos = remove(cell->begin(),cell->end(),e2);
			cell->erase(pos,cell->end());
			// mark e2 as processed 
			edgeCollected.insert(e2);
			this->removeEdgeFromLinkedList(e2);
			edges2delete.push_back(e2);
		}
		else{
			Point2D* other = nEdges[i]->getOtherPoint(point);
			/* splice at other */
			Edge* next = nEdges[i]->getNext(other);
			Edge* prev = nEdges[i]->getPrev(other);
			if( other->entry == nEdges[i] ){
				other->entry = next;
			}
			next->setPrev(other, prev);
			prev->setNext(other, next);
			this->removeEdgeFromLinkedList(nEdges[i]);
			edges2delete.push_back(nEdges[i]);
		}
		// mark e2 as processed 	
		edgeCollected.insert(nEdges[i]);
	}
	
	edgeCollected.clear();
	
	// refactor edges 
	for (unsigned int i = 0; i < edges2delete.size(); i++) 
	{
		// remove from edge list
		edges2delete[i]->refactorEdge();
	}

	edges2delete.clear();
	nEdges.clear();
}


vector<Edge*> Triangulation::deleteEdge(Edge* edge)
{
	Point2D* org = edge->org;
	Point2D* dest = edge->dest;
	vector<Point2D*> neighbors = edge->getNeighbors();
	// triangles of the cell 
	vector<Triangle*> cellTriangles = this->getTriangles(edge);
	// get cell hull of this point
	vector<Edge*> cell = edge->getCellHull(this->xMax,
											this->xMin,
											this->yMax,
											this->yMin);
	// remove point from cell 
	removePoint(org, &cell);
	removePoint(dest, &cell);
	
	vector<Edge*> newEdges = this->integrateSubTriangulation(cell, cellTriangles, neighbors);
	org->thinned = true;
	dest->thinned = true;
	org->entry = NULL;
	dest->entry = NULL;

	return newEdges;
}


vector<Edge*> Triangulation::deletePoint(Point2D* point)
{
	
	vector<Point2D*> neighbors = point->getNeighbors();
	/* triangles of the cell */
	vector<Triangle*> cellTriangles = point->getTriangles();
	/* get cell hull of this point*/
	vector<Edge*> cell = point->getCellHull(this->xMax,
											this->xMin,
											this->yMax,
											this->yMin);
	/* remove point from cell */
	vector<Edge*> newEdges;
	this->removePoint(point, &cell);
	
	vector<Edge*> ne = this->integrateSubTriangulation(cell,cellTriangles,neighbors);
	for( unsigned int i = 0; i < ne.size(); i++ ){
		newEdges.push_back(ne[i]);
	}
	
	point->thinned = true;
	point->entry = NULL;

	// free triangle objects 
	Triangle::deleteTriangles(cellTriangles);

	return newEdges;
} 


void Triangulation::debug()
{
	Edge* nEdge = this->entry;
	while( true )
	{
	  nEdge->debug();
	  // iteration end 
	  if( nEdge == nEdge->succ )
	  {
	    break;	
	  }
	  else
	  {
		nEdge = nEdge->succ;
	  }
	}
}


vector<Point2D*> Triangulation::getNodes()
{
  vector<Point2D*> ret;
  for (unsigned int i = 0; i < this->nodes.size(); ++i) 
  {
    if(!this->nodes[i]->thinned) ret.push_back(this->nodes[i]);
  }
  return ret;
}


void Triangulation::optimize()
{
	unsigned int i,j,k;
	int x,y,count,nl;

	M3Matrix pointVisited(this->nbCols, this->nbRows);
	
	// double eps = -0.000000000001;
	// To get the positions in the sparse matrix
	M3Matrix PointPosition(this->nbCols,this->nbRows);
  
	M3Matrix imgM(this->nbCols, this->nbRows);
	for(i=0; i<this->nodes.size(); i++) {
		imgM[this->nodes[i]->x][this->nodes[i]->y] = this->nodes[i]->f;
	}
  
	vector<Point2D*> triNodes = this->getNodes();
	unsigned int plen = triNodes.size();
	M3Matrix IPosition(plen+1,1);

	int maxneighbors = 0;

	count = 0;
	for (i=0;i<plen;i++) {
		nl = triNodes[i]->getNeighbors().size();
		if (nl > maxneighbors){
			maxneighbors = nl;
		}
		PointPosition[triNodes[i]->x][triNodes[i]->y] = (double)count;
		count++;
	}
		
	M3Matrix JPosition(plen,maxneighbors);

	count = 0;
	for (i=0;i<plen;i++) {
		vector<Point2D*> neighbors = triNodes[i]->getNeighbors();
		nl = neighbors.size();
		IPosition[i][0] = count;
		count += nl+1;
		for(j=0;j<neighbors.size();j++) {
			JPosition[i][j] = PointPosition[neighbors[j]->x][neighbors[j]->y];
		}
	}
  
	IPosition[plen][0] = (double)count;

	M3Matrix SparseMatrix(count,1);

	M3Matrix fS(plen,1);

	// now fills the Phi matrix
	unsigned int NbN;

	double ab,ac,bc;
	double aa,bb,cc;
	double fa,fb,fc;
	int xa, xb, xc, ya, yb, yc;
  
//			  SortedSet sortedTriangles = new TreeSet();
//			  sortedTriangles.addAll(this.getTriangles());
	vector<Triangle*> triangles = this->getTriangles();
	for(unsigned int t = 0;t<triangles.size();t++) {

		xa = triangles[t]->point1->x;
		ya = triangles[t]->point1->y;
		xb = triangles[t]->point2->x;
		yb = triangles[t]->point2->y;
		xc = triangles[t]->point3->x;
		yc = triangles[t]->point3->y;

		ab =0.;
		ac =0.;
		bc =0.;
		
		fa = 0.;
		fb = 0.;
		fc = 0.;
				
		int ia = (int)PointPosition[xa][ya];
		int ib = (int)PointPosition[xb][yb];
		int ic = (int)PointPosition[xc][yc];
		
		aa = 0.;
		bb = 0.;
		cc = 0.;
		
		int xmin = xa;
		int ymin = ya;
		if (xb<xmin) xmin = xb;
		if (xc<xmin) xmin = xc;
		if (yb<ymin) ymin = yb;
		if (yc<ymin) ymin = yc;
		
		int xmax = xa;
		int ymax = ya;
		if (xb>xmax) xmax = xb;
		if (xc>xmax) xmax = xc;
		if (yb>ymax) ymax = yb;
		if (yc>ymax) ymax = yc;
		
		double dxa = (double)xa;
		double dxb = (double)xb;
		double dxc = (double)xc;
		double dya = (double)ya;
		double dyb = (double)yb;
		double dyc = (double)yc;
		
		// det = determinant; ba, bb, bc barycentric coordinates of p
		double det = (dxa*(dyb-dyc)) + (dxb*(dyc-dya)) + (dxc*(dya-dyb));

		// loop on all the points of the triangle
		for (x=xmin;x<=xmax;x++) {
			for (y=ymin;y<=ymax;y++) {
				if (pointVisited[x][y] != 1.) {
					double dx = (double)x;
					double dy = (double)y;

					double da = (dx*(dyb-dyc) + dxb*(dyc-dy) + dxc*(dy-dyb))/det;
					double db = (dxa*(dy-dyc) + dx*(dyc-dya) + dxc*(dya-dy))/det;
					double dc = (dxa*(dyb-dy) + dxb*(dy-dya) + dx*(dya-dyb))/det;
					if (da>=0.0 && db>=0.0 && dc>=0.0) 
					{ 
						// Filling of the Gram matrix : <Phi_i,Phi_j>

						// Non diagonal elements
						ab += da*db;
						ac += da*dc;
						bc += db*dc;
						
						// Filling of the second member vector <Phi_i,I>
						int v = (int)imgM[x][y];
						fa += da*v;
						fb += db*v;
						fc += dc*v;
						
						// Diagonal elements
						aa += da*da;
						bb += db*db;
						cc += dc*dc;
						
						pointVisited[x][y] = 1.;
					}
				}
			}
		}


		fS[ia][0] += fa;
		fS[ib][0] += fb;
		fS[ic][0] += fc;
		
		int jab = 0;
		int jac = 0;
		int jba = 0;
		int jbc = 0;
		int jca = 0;
		int jcb = 0;
		
		// Point a
		NbN = (int)(IPosition[ia+1][0]) - (int)(IPosition[ia][0]) - 1;
		for (k=0;k<NbN;k++) {
			if (JPosition[ia][k] == ib) jab = k;
			if (JPosition[ia][k] == ic) jac = k;
		}

		// Point b
		NbN = (int)(IPosition[ib+1][0]) - (int)(IPosition[ib][0]) - 1;
		for (k=0;k<NbN;k++) {
			if ((int)JPosition[ib][k] == ia) jba = k;
			if ((int)JPosition[ib][k] == ic) jbc = k;
		}
		
		// Point c
		NbN = (int)(IPosition[ic+1][0]) - (int)(IPosition[ic][0]) - 1;
		for (k=0;k<NbN;k++) {
			if ((int)JPosition[ic][k] == ia) jca = k;
			if ((int)JPosition[ic][k] == ib) jcb = k;
		}

		// Increments SparseMatrix (which is in fact a vector)
		// diagonal elements
		SparseMatrix[(int)(IPosition[ia][0])][0] += aa;
		SparseMatrix[(int)(IPosition[ib][0])][0] += bb;
		SparseMatrix[(int)(IPosition[ic][0])][0] += cc;

		// Other elements
		SparseMatrix[(int)(IPosition[ia][0]) + jab + 1][0] += ab;
		SparseMatrix[(int)(IPosition[ia][0]) + jac + 1][0] += ac;
		SparseMatrix[(int)(IPosition[ib][0]) + jba + 1][0] += ab;
		SparseMatrix[(int)(IPosition[ib][0]) + jbc + 1][0] += bc;
		SparseMatrix[(int)(IPosition[ic][0]) + jca + 1][0] += ac;
		SparseMatrix[(int)(IPosition[ic][0]) + jcb + 1][0] += bc;

	}

	// Resolution of the system : PCG method
	M3Matrix X(plen,1);
	M3Matrix r(plen,1);
	M3Matrix p(plen,1);
	M3Matrix q(plen,1);
	M3Matrix s(plen,1);

	double gamma;
	double tmp;

	// Initialization of X,r,p,s,gamma with the image values
	for (i=0;i<plen;i++) 
	{
    	X[i][0] = triNodes[i]->f;
	}

	for (i=0;i<plen;i++) {
		tmp = 0.;
		NbN = (int)(IPosition[i+1][0]) - (int)(IPosition[i][0]);
		tmp = SparseMatrix[(int)(IPosition[i][0])][0] * X[i][0];
		for (j=1;j<NbN;j++){
			tmp += SparseMatrix[(int)(IPosition[i][0])+j][0] * X[(int)(JPosition[i][j-1])][0];
		}
		r[i][0] = fS[i][0] - tmp;
	}


	for (i=0;i<plen;i++) {
		tmp = 0.;
		NbN = (int)(IPosition[i+1][0]) - (int)(IPosition[i][0]);
		tmp = SparseMatrix[(int)(IPosition[i][0])][0] * r[i][0];
		for (j=1;j<NbN;j++){
			tmp += SparseMatrix[(int)(IPosition[i][0])+j][0] * r[(int)(JPosition[i][j-1])][0];
		}
		p[i][0] = tmp;
		s[i][0] = tmp;
	}

	gamma = 0.;
	for (i=0;i<plen;i++) gamma += s[i][0]*s[i][0];
	
	// Iterations of the method
	unsigned int NbIt = 100;
	
	double alpha,beta;
	double tmp1=0.;
	double Tol = 0.0000000001;
	
	if (gamma>Tol) 
	{
		for (k=0;k<NbIt;k++) 
		{
		     // q_k = Phi * p_k
			for (i=0;i<plen;i++) {
				tmp = 0.;
				NbN = (int)(IPosition[i+1][0]) - (int)(IPosition[i][0]);
				tmp = SparseMatrix[(int)(IPosition[i][0])][0] * p[i][0];
				for (j=1;j<NbN;j++){
					tmp += SparseMatrix[(int)(IPosition[i][0])+j][0] * p[(int)(JPosition[i][j-1])][0];
				}
				q[i][0] = tmp;
			}

			// alpha_k = gamma_k/N_2(q_k)
			tmp1 = 0.;
			for (i=0;i<plen;i++){
				tmp1 += q[i][0]*q[i][0];
			}

			alpha = gamma/tmp1;
			
			// X_{k+1} =X_k + alpha_k * p_k
			// r_{k+1} =r_k - alpha_k * q_k
			for (i=0;i<plen;i++) {
				X[i][0] = X[i][0] + alpha * p[i][0];
				r[i][0] = r[i][0] - alpha * q[i][0];
			}

			// s_{k+1}=A*r_{k+1}
			for (i=0;i<plen;i++) {
				tmp = 0.;
				NbN = (int)(IPosition[i+1][0]) - (int)(IPosition[i][0]);
				tmp = SparseMatrix[(int)(IPosition[i][0])][0] * r[i][0];
				for (j=1;j<NbN;j++){
					tmp += SparseMatrix[(int)(IPosition[i][0])+j][0] * r[(int)(JPosition[i][j-1])][0];
				}
				s[i][0] = tmp;
			}

			// gamma_{k+1} = N_2(s_{k+1})
			tmp1 = 0.;
			for (i=0;i<plen;i++){
				tmp1 += s[i][0]*s[i][0];
			}

			//
			beta = tmp1/gamma;
			gamma = tmp1;
			
			// p_{k+1} = s_{k+1} + beta * p_k
			for (i=0;i<plen;i++){
				p[i][0] = s[i][0] + beta * p[i][0];
			}
			if (gamma < Tol){
				k = NbIt;
			}
		}
	}


	for(unsigned int t = 0; t<triangles.size();t++) 
	{
		int ia = (int)(PointPosition[triangles[t]->point1->x][triangles[t]->point1->y]);
		int ib = (int)(PointPosition[triangles[t]->point2->x][triangles[t]->point2->y]);
		int ic = (int)(PointPosition[triangles[t]->point3->x][triangles[t]->point3->y]);

		int av = (int)(X[ia][0]);
		int bv = (int)(X[ib][0]);
		int cv = (int)(X[ic][0]);

		triangles[t]->point1->f = (av<0)? 0 : (av>255)? 255 : av;
		triangles[t]->point2->f = (bv<0)? 0 : (bv>255)? 255 : bv;
		triangles[t]->point3->f = (cv<0)? 0 : (cv>255)? 255 : cv;
	}
}	


void Triangulation::writeNodes(string fn)
{
	ofstream out(fn.c_str());
	vector<Point2D*> nds = this->getNodes();
	for (int i = nds.size()-1; i>=0; --i) 
	{
	  out<<nds[i]->x<<","<<nds[i]->y<<endl;
	}
	out.close();
}


void Triangulation::writeEdges(string fn)
{
  ofstream out(fn.c_str());
  vector<Edge*> el = this->getEdges();
  for ( unsigned int i = 0; i < el.size(); i++ ) 
  {
    out<<el[i]->org->x<<","<<el[i]->org->y<<","<<
		el[i]->dest->x<<","<<el[i]->dest->y<<","<<endl;
  }
  out.close();
}


void Triangulation::checkTri()
{
	Edge* nEdge = this->entry;
	while( true )
	{
		if(nEdge->pred == nEdge && nEdge != this->entry)
		{
			this->entry->debug();	
		}
		// iteration end 
		if( nEdge == nEdge->succ )
		{
			break;	
		}
		else
		{
			nEdge = nEdge->succ;
		}
	}
	
}

vector<Point2D*> Triangulation::getThinnedNodes()
{
	vector<Point2D*> ret;
	for (unsigned int i = 0; i < this->nodes.size(); ++i) 
	{
		if(this->nodes[i]->thinned) ret.push_back(this->nodes[i]);
	}
	return ret;
}

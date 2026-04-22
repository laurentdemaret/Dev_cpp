#ifndef TRIANGULATION_H_
#define TRIANGULATION_H_

#include <vector>
#include <set>
#include <string>

#include "Point2D.h"
#include "Edge.h"
#include "PointGrid.h"
#include "matrix.h"

#define	TRIANGULATION_POOL_SIZE (10)

class Triangulation
{

public:
	// attributes 
	std::vector<Point2D*> nodes;
	int nbRows, nbCols, xMax, yMax, xMin, yMin;
	// edge linked list
	Edge* entry;
	
	// constructors/destructor 
	Triangulation();
	Triangulation(PointGrid*);
	virtual ~Triangulation();
	
	static bool positionComparator(const void*, const void*);
	void triangulate();
	std::vector<Point2D*> getAttachedPoints(Point2D*);
	std::vector<Point2D*> getAttachedPoints(std::vector<Point2D*>);
	std::vector<Point2D*> getAttachedPoints(unsigned int, 
														  unsigned int, 
														  unsigned int,
														  unsigned int);
	std::vector<Triangle*> getTriangles(std::vector<Edge*>);														  
	std::vector<Triangle*> getTriangles(Edge*);
	std::vector<Triangle*> getTriangles();
	Triangle* processTrianglesByEdge(Edge*, 
													Point2D*);
	std::vector<Edge*> deleteEdge(Edge*);
	std::vector<Edge*> insertPoint(Point2D*);
	
	// constructor from pool 
	static Triangulation* makeTriangulation(std::vector<Point2D*>);
	void deleteTriangulation();
	void deleteTriangulation(std::vector<Triangle*>);
	std::vector<Edge*> deletePoint(Point2D* point);
	std::vector<Edge*> getEdges();
	std::vector<Point2D*> getNodes();
	std::vector<Point2D*> getThinnedNodes();
	void optimize();
	void debug();
	
	void writeNodes(std::string);
	void writeEdges(std::string);

private:
	void divide(int,int,Edge**,Edge**);
	
	Edge** divide(int,int);
	
	Edge* merge(Edge*,Point2D*,Edge*,Point2D*);
	void getLowerTangent(Edge*, Point2D*, Edge*, Point2D*,
										Point2D**, Point2D**, Edge**, Edge**);
	std::vector<Edge*> integrateSubTriangulation(std::vector<Edge*>, std::vector<Triangle*>, std::vector<Point2D*>);
	void removeEdgeFromLinkedList(Edge*);
	void removeEdge(Edge*);
	void addEdge(Edge*);
	void splitEdge(Point2D*,Edge*,Edge**,Edge**);
	void removePoint(Point2D*,std::vector<Edge*>*);
	
	// attributes 
	static std::set<Triangulation*> pool;
	
	// dummy edge 
	static Edge* dummyEdge;

	// tests 
	void checkTri();
};

#endif /*TRIANGULATION_H_*/

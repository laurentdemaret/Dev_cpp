#ifndef EDGE_H_
#define EDGE_H_

#include <vector>
#include <set>
#include "Point2D.h"

#define LEFT (0)
#define RIGHT (1)
#define	EDGE_POOL_SIZE (100)



class Edge
{
public:

	// constructors/destructors 
	Edge();
	Edge(Point2D*, Point2D*);
	virtual ~Edge();

    // attributes
	Point2D* org;
	Point2D* dest;
	Edge* onext;
	Edge* oprev;
	Edge* dnext;
	Edge* dprev;
	double epsilon;	
	/* linked list of edges: used to get edges of */
	Edge* pred;
	Edge* succ;
	void* fh_el;

	// methods 
	Point2D* getOtherPoint(Point2D*);
	bool assumeUniquness(int);
	void swap();
	void connect(Edge*, Point2D*);
	void splice(Edge*, Point2D*);
	Edge* getNext(Point2D*);
	Edge* getPrev(Point2D*);
	void setNext(Point2D*,Edge*);
	void setPrev(Point2D*,Edge*);
	Edge* join(Point2D*, Edge*, Point2D*, int);
	void refactorEdge();
	void deleteEdge();
	bool checkTriangle(Point2D*);
	/* constructor from pool */
	static Edge* makeEdge(Point2D*, Point2D*);
	static Edge* makeDummyEdge();
	bool isVertex(Point2D*);
	
	bool isVirtualVertex(Point2D*);
	bool isVirtualEqual(Edge*);
	
	bool isConvexHull(unsigned int, 
							unsigned int, 
							unsigned int, 
							unsigned int);
	void adjustQuadStructure(Point2D*, 
									 Edge*, 
									 Edge*, 
									 Point2D*, 
									 Edge*, 
									 Edge*);
	void deleteEdgeInsideCell(std::set<Edge*>);
	bool checkEdge(Triangle*, int);
	
	std::vector<Point2D*> getNeighbors();
	std::vector<Edge*> getEdges();
	std::vector<Edge*> getCellHull(int, int, int, int);
	
	static void fillEdgePool();
	
	// debug 
	void debug();
	void debugStruct();

	static int new_edge_cnt;
	static std::vector<Edge*> pool;
	
private:
	// methods 
	bool reverseEdge();
	void initEdge();
	
	// attributes	
};

#endif /*EDGE_H_*/

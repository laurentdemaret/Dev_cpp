#ifndef TRIANGLE_H_
#define TRIANGLE_H_

#include <vector>
#include <set>

#include "Point2D.h"
#include "Edge.h"

#define	TRIANGLE_POOL_SIZE (1000)

class Triangle
{
public:

	// constructors/destructors 
	Triangle(Point2D*,Point2D*,Point2D*);
	Triangle(Edge*,Point2D*);
	virtual ~Triangle();
	
	// attributes
	Point2D* point1;
	Point2D* point2;
	Point2D* point3;
	
	// methods 
	bool isVertex(Point2D*);
	void deleteTriangle();
	
	// constructor from pool 
	static Triangle* makeTriangle(Point2D*,Point2D*,Point2D*);
	static Triangle* makeTriangle(Edge*,Point2D*);
	static Triangle* makeTriangle(Triangle*);
	static void deleteTriangles(std::vector<Triangle*>);
	bool insideTriangle(Point2D*);
	bool insideTriangle(double, double);
	void getEdges(Edge**,Edge**,Edge**);
	Point2D* getOtherPoint(Point2D*, Point2D*);
	Point2D* getOtherPoint(Edge*);
	Triangle* getNeighbourByEdge(Edge*);
	Point2D* getNeighbourPointByEdge(Edge* edge);
	
	double GetArea();
    double GetAspectRatio();
	
	static int new_triangle_cnt;
	
	// debug 
	void debug();

	static std::vector<Triangle*> pool;
	
private:	
	// methods 
	void initTriangle(Point2D*, Point2D*, Point2D*);	
};

#endif /*TRIANGLE_H_*/

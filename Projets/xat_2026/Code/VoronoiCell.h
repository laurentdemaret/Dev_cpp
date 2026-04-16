#ifndef VORONOICELL_H_
#define VORONOICELL_H_

#include <vector>
#include <set>

// Forward declaration of Point2D, Edge,Triangle for cyclic dependency
class Point2D;
class Edge;
class Triangle;

#define	POINT2D_POOL_SIZE (100)

class VoronoiCell
{
public:
	// constructors/destructors
    //Point2D();
    //Point2D(int, int);
    //Point2D(int, int, int);
    //virtual ~Point2D();
	
	// attributes
    //int x, y, f; //x,y: are the coordinates of the point / f is an associated

    Point2D* Centroid;

    std::vector<Point2D*> pixels;
	
	std::vector<Triangle*> getTriangles();
	std::vector<Edge*> getEdges();
	Edge* getEdge(Point2D* b);
	
		
    // constructor from pool
    //static Point2D* makePoint2D(int, int, int);


    //void VoronoiCell();

//private:
    // attributes
//	static std::vector<Point2D*> pool;
};

#endif //VoronoiCell_H_

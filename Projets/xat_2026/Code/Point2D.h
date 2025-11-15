#ifndef Point2D_H_
#define Point2D_H_


#include <vector>
#include <set>

// Forward declaration of Edge,Triangle for cyclic dependency
class Edge;
class Triangle;

#define	POINT2D_POOL_SIZE (100)

class Point2D
{
public:
	// constructors/destructors
	Point2D(int, int);
	Point2D(int, int, int);
	virtual ~Point2D();
	
	// attributes
	int x, y, f;
	double epsilon;
	bool gridBound;
	bool convexHull;
	Edge* entry;
	bool thinned;
	bool dirty;
	bool exchanged;
	/* a pointer to the corresponding fibheap element */
	void* fh_el;
	Triangle* thinnedTriangle;
	Point2D* bestNeighbor;
	double bestNeighborDiff;
	void* exchange_fh_el;
	// copy: needed producing temp triangulations
	Point2D* clone;
	

	// methods
	int crossProduct(Point2D*, Point2D*);
	int crossProduct(Point2D*, Point2D*, Point2D*);
	int dotProduct(Point2D*, Point2D*, Point2D*);
	double sin(Point2D*, Point2D*, Point2D*);
	double computeAngleCos(int, int, int, int);
	double computeAngleCos(Point2D*, Point2D*);
	int inCircle(Point2D*, Point2D*, Point2D*);
	bool ccw(Point2D*, Point2D*);
	bool checkScore(Point2D*, Point2D*, Point2D*);
	static std::vector<Point2D*> clonePoint2D(std::vector<Point2D*>);
	Point2D* clonePoint2D();
	
	bool isConvexHull(Point2D*, int, int, int, int);
	std::vector<Edge*> getCellHull(int, int, int, int);
	
	/* */
	std::vector<Point2D*> getNeighbors();
	std::vector<Point2D*> getNeighborsCopy();
	std::vector<Triangle*> getTriangles();
	std::vector<Edge*> getEdges();
	Edge* getEdge(Point2D* b);
	
	void deletePoint2D();
	
	
	/* some helper functions */
	static std::vector<Point2D*> set2vector(std::set<Point2D*>);
	static std::set<Point2D*> vector2set(std::vector<Point2D*>);
	static std::set<Point2D*> makeUniqueList(std::set<Point2D*>, std::set<Point2D*>);
	static std::vector<Point2D*> makeUniqueList(std::vector<Point2D*>, std::vector<Point2D*>);
	static void addAll(std::set<Point2D*>*,std::vector<Point2D*>);
	static void addAll(std::set<Point2D*>*,std::set<Point2D*>);
	
	/* constructor from pool */
	static Point2D* makePoint2D(int, int, int);

	bool lessThan(Point2D*);

	/* debug */
	void debug();
	void debugNeighbors();
	static int new_point2d_cnt;

private:
	/* attributes*/
	static std::vector<Point2D*> pool;
};

#endif /*Point2D_H_*/

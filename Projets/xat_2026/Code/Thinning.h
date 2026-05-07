#ifndef THINNING_H_
#define THINNING_H_

#include <string>

#include "Triangulation.h"
#include "fib.h"

class Thinning
{
public:
	Thinning(Triangulation*);
	virtual ~Thinning();
	
	static bool positionComparator(const void*, const void*);
	static bool rankComparator( const void* x, const void* y);
	// attributes 
	Triangulation* triangulation;
	static int exchangeRadius;
	double lambda; //this parameter is used for the L_tau penalisation
	double tau;
	
	//Matrix with the original pixel greyvalues 
	M3Matrix OriginalImage;

	//Matrix with the current approximated greyvalues 
	M3Matrix ApproximatedImage;

    // methods
    void thinningAT9(unsigned int tcount);
	void thinningAT5(unsigned int);
	void thinningAT6(unsigned int);
	void thinningAT8(unsigned int);
	void fastThinning(unsigned int);
	int exchange(unsigned int);
	// IO
	void printNodes(const std::string& filename, int option);
	void printThinnedNodes(const std::string& filename, int option);
	void printEdges(const std::string& filename, int option);
	void printTriangles(const std::string& filename, int option);
	void printEps(const std::string& filename);
	void printTriangulationEps(const std::string& filename);
	void printTriangulationOff(const std::string& filename);

private:
	fibheap* nodeHeap;	
	fibheap* edgeHeap;
	fibheap* thinnedNeighborHeap;
	
    // methods
    double calculateErrorVoronoi(Triangle* triangle, Point2D* point);

    double calculatePointSignificanceVoronoi(Point2D* point);
	double calculatePointSignificance(Point2D*);
	double calculatePointSignificance_ForAT8(Point2D*);
	double calculateEdgeSignificance(Edge* edge);
    double calculateSignificanceVoronoi(std::vector<Triangle*> oTriangles,
                                        std::vector<Triangle*> nTriangles,
                                        std::vector<Point2D*> attachedPoints,
                                        Point2D* point);
	double calculateSignificance(std::vector<Triangle*>,
											 std::vector<Triangle*>,
											 std::vector<Point2D*>,
											 Point2D*);
	double calculateSignificance(std::vector<Triangle*>,
											 std::vector<Triangle*>,
											 std::vector<Point2D*>,
											 Edge*);


	double calculateNewError(std::vector<Triangle*>, Point2D*);

	double calculateError(Triangle*, Point2D*);
	double calculateError(Triangle*, Point2D*,double);

    //Fonctions rajoutées pour "corriger l'erreur ? " (avril 2026)
    /*double calculateSignificance26(std::vector<Triangle*> oTriangles,
                                   std::vector<Triangle*> nTriangles,
                                   std::vector<Point2D*> attachedPoints,
                                   Point2D* point);*/
    /*double calculateNewError26(std::vector<Triangle*> nTriangles,
                               Point2D* point);*/


    //Gebraucht um den Voronoi Fehler zu berechnen !
    double calculateNewError_Voronoi(std::vector<Triangle*> nTriangles,
                                     Point2D* point);
    double calculateError_Voronoi(Triangle* triangle, Point2D* point);

    void calculateSignificanceVoronoiForHeap(Point2D* point);
    void calculateSignificanceForHeap(Point2D*);
	void calculateSignificanceForHeap_ForAT8(Point2D*);
	void calculateSignificanceForEdgeHeap(Edge*);
	void deletePoint(Point2D*);
	std::vector<Triangle*> deletePointWithSigRecalc(Point2D*);
	void deleteEdge(Edge*);
	void updateHeaps(std::vector<Point2D*>,Point2D*,Edge*);
	void calculatePairSignificance(Point2D*);
	void setThinnedTriangle(Point2D*, std::set<Point2D*>);
	std::vector<Point2D*> getExchangeCandidats(Point2D* point);	
	
	//point2D_array getThinnedPoints();
};

#endif /*THINNING_H_*/

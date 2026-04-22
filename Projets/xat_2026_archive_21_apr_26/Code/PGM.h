#ifndef PGM_H_
#define PGM_H_

#include <vector>
#include "PointGrid.h"
#include "matrix.h"
#include "Triangulation.h"

class PGM
{
public:
	PGM();
	virtual ~PGM();
	
	// methods 
	static PointGrid* readFile(const char*);
    static void renderTriangulation(M3Matrix&,
                                         int NbRows,
                                         int NbCols,
                                         const char*  outputFilename,
                                         std::vector<Triangle*>);

    static void renderVoronoi(M3Matrix& image,
                            int NbRows,
                            int NbCols,
                            const char*  outputFilename,
                              std::vector<Triangle*> dtta);

	static void optimize(int, 
								int, 
								M3Matrix&, 
								std::vector<Triangle*>,
								std::vector<Point2D*> nodes);
};

#endif /*PGM_H_*/

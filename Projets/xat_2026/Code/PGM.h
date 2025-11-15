#ifndef PGM_H_
#define PGM_H_

#include <vector>
/*#include "../triangulation/PointGrid.h"
#include "../tools/matrix/matrix.h"
#include "../triangulation/Triangulation.h"*/
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
										 int, 
										 int, 
										 const char* , 
										 std::vector<Triangle*>);	
	static void optimize(int, 
								int, 
								M3Matrix&, 
								std::vector<Triangle*>,
								std::vector<Point2D*> nodes);
};

#endif /*PGM_H_*/

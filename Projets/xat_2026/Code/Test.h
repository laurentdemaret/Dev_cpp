#ifndef TEST_H_
#define TEST_H_

#include <vector>
/*#include "../triangulation/Point2D.h"
#include "../triangulation/Edge.h"
#include "../triangulation/Triangle.h"*/
#include "Point2D.h"
#include "Edge.h"
#include "Triangle.h"


class Test
{
public:
	Test();
	virtual ~Test();
	
	static void list(std::vector<Point2D*>);
	static void list(std::vector<Edge*>);
	static void list(std::vector<Triangle*>);
	
};

#endif /*TEST_H_*/

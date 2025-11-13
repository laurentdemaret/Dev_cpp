#ifndef EDGE2D_H
#define EDGE2D_H

#include <iostream>
#include <map>
#include <set>

class Edge2D
{
  public:
  Edge2D();
  //Point2D(double x, double y);       // constructor
  //Point2D(const Point2D& p_src); //copy constructor


 /* Point2D& operator=(const Point2D& p_src)
  {
      if (this != &p_src) {  // protection contre l'auto-affectation
          xm = p_src.xm;
          ym = p_src.ym;
          // copier d'autres attributs ici si nécessaire
      }
      return *this;
  };*/

  //void Move(double dx,double dy);

  // destructor
  ~Edge2D(){};

  //void Display();

  //void Move

  //Coordinates of the point

  //Identificator of the point
  int id;

  int p1;
  int p2;

  //
  std::set<int> triangle_ids;
};


#endif // EDGE2D_H



#ifndef POINT2D_H
#define POINT2D_H

#include <iostream>
#include <map>
#include <set>

class Point2D
{
  public:
  Point2D();
  Point2D(double x, double y);       // constructor
  Point2D(const Point2D& p_src); //copy constructor


  Point2D& operator=(const Point2D& p_src)
  {
      if (this != &p_src) {  // protection contre l'auto-affectation
          xm = p_src.xm;
          ym = p_src.ym;
          // copier d'autres attributs ici si nécessaire
      }
      return *this;
  };

  void Move(double dx,double dy);

  // destructor
  ~Point2D(){};

  void Display();

  //void Move

  //Coordinates of the point
  double xm;
  double ym;

  //Identificator of the point
  int id;


  //
  std::set<int> edge_ids;
  std::set<int> triangle_ids;
};


#endif // POINT2D_H



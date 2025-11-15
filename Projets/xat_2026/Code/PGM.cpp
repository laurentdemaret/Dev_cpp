#include <cstdio>
#include <fstream>
#include <vector>
#include <set>
#include <iostream>
#include <cmath>

#include "PGM.h"
/*#include "../tools/matrix/matrix.h"
#include "../triangulation/Triangle.h"
#include "../triangulation/Triangulation.h"*/

#include "matrix.h"
#include "Triangle.h"
#include "Triangulation.h"

using namespace std;

PGM::PGM()
{
}

PGM::~PGM()
{
}


PointGrid* PGM::readFile(const char* filename)
{
	vector<Point2D*> nodes;
    int b,h;
	ifstream input(filename);
	//assert(input.is_open());
	char buf[1000];
 	input.getline(buf, 1000);
 	char buf1[1];
	input.read(buf1,1);
	while (*buf1 == '#') {
		input.getline(buf, 1000);
		input.read(buf1,1);
	}

	input.unget();
	input >> b;
	input >> h;

	int f;
	input >> f;
    for (int i=0;i<b;i++)
    {
      for (int j=0;j<h;j++)
      {
        input >> f;
        Point2D* p = Point2D::makePoint2D(j,i,f);
        nodes.push_back(p);
      }
	}
  
	PointGrid* grid = new PointGrid(nodes, h, b);  
	return grid;
}


/* static void anisotropicDiffusion(M3Matrix& image,bool** fixed, int NBRows,int NBCols,int diffusionIterations)
{
	// gradients
	double grad1;
	double grad2;
	// structure tensor
	double g1;
	double g2;
	double g3;
	// tensor field
	double t1;
	double t2;
	double t3;
	
	double** cube = new double*[NBCols];
	for(int ix=0;ix<NBCols;++ix)
	{
		cube[ix] = new double[NBRows];
		for(int iy=0;iy<NBRows;++iy){
			cube[ix][iy] = image[ix][iy];
		}
	}
	
	// consts
	double dt        	= .1;  // Adapting time step
	double alpha  		= .001;
	double sqrt2		= sqrt(2.);
	
	
	for(int i=0;i<diffusionIterations;++i){
		// get initial stats for later normalizing
		double initial_max=cube[0][0], initial_min=cube[0][0];
		for(int x=0;x<NBCols;++x){
			for(int y=0;y<NBRows;++y){
				double pix = cube[x][y];
				if (pix>initial_max) initial_max=pix;
				if (pix<initial_min) initial_min=pix;
			}
		}
		
		// compute gradients
		double Ipp,Icp,Inp=0,Ipc,Icc,Inc=0,Ipn,Icn,Inn=0;
		for(int x=0;x<NBCols;++x){
			int px=x-1; if (px<0) px=0;
			int nx=x+1; if (nx==NBCols) nx--;
			for(int y=0;y<NBRows;++y){
				// escape fixed points 
				//				if( fixed[x][y] == 1 ) continue;
				int py=y-1; if (py<0) py=0;
				int ny=y+1; if (ny==NBRows) ny--;
				Ipp=cube[px][py];
				Ipc=cube[px][y] ;
				Ipn=cube[px][ny];
				Icp=cube[x] [py];
				Icn=cube[x] [ny];
				Inp=cube[nx][py];
				Inc=cube[nx][y] ;
				Inn=cube[nx][ny];
				grad1 = (float)(sqrt2*Inn+Inc+sqrt2*Inp-sqrt2*Ipp-Ipc-sqrt2*Ipn)/(2+4.*sqrt2);
				grad2 = (float)(sqrt2*Inn+Icn+sqrt2*Ipn-sqrt2*Ipp-Icp-sqrt2*Inp)/(2+4.*sqrt2);
				
				
				// compute structure tensor field G
				double fx = grad1;
				double fy = grad2;
				g1 = fx*fx;
				g2 = fx*fy;
				g3 = fy*fy;
				
				// compute the tensor field T, used to drive the diffusion
				
				// eigenvalues:
				double a = g1, b = g2, c = g2, d = g3, e = a+d;
				double f = sqrt(e*e-4*(a*d-b*c));
				double l1 = 0.5*(e-f), l2 = 0.5*(e+f);
				// more precise computing of quadratic equation
				if (e>0) { if (l1!=0) l2 = (a*d - b*c)/l1; }
				else     { if (l2!=0) l1 = (a*d - b*c)/l2; }
				
				
				//						double val1 = (l2 / drange2);
				//						double val2 = (l1 / drange2);
				double f1,f2;
				// slight cheat speedup for default a1 value
				//						f1 = (a1==.5) ? (1/sqrt(1.0f+val1+val2)) : (pow(1.0f+val1+val2,-a1));
				//						f2 = pow(1.0+val1+val2,-a2);
				// coherence 
				f1 = alpha + (1 - alpha)*exp(-3.31488/pow(l1 - l2,8));
				f2 = alpha;
				// edge enhanced 
				//						f1 = 1-exp(-3.31488/pow(l1/3.6,8));
				//						f2 = 1;
				
				// eigenvectors:
				double u, v, n;
				if (fabs(b) > fabs(a-l1)) { u = 1; v = (l1-a)/b; }
				else { if (a-l1!=0) { u = -b/(a-l1); v = 1; }
				else { u = 1; v = 0; }
				}
				n = sqrt(u*u+v*v);
				if( n > 0 ){
					u/=n; v/=n;
				}
				//cout<<"--------------------------------------------"<<endl;
				//cout<<"\t\t"<<u<<"\t\t"<<v<<endl;
				//cout<<"\t\t"<<l1<<"\t\t"<<l2<<endl;
				
				double vec1 = u;
				double vec2 = v;
				float vec11 = vec1*vec1, vec12 = vec1*vec2, vec22 = vec2*vec2;
				t1 = f1*vec11 + f2*vec22;
				t2 = (f1-f2)*vec12;
				t3 = f1*vec22 + f2*vec11;
				//cout<<"--------------------------------------------"<<endl;
				//cout<<"\t\t"<<g1<<"\t\t"<<g2<<endl;
				//cout<<"\t\t"<<g2<<"\t\t"<<g3<<endl;
				
				// compute the PDE velocity and update the iterated image
				Ipp=cube[px][py];
				Ipc=cube[px][y] ;
				Ipn=cube[px][ny];
				Icp=cube[x] [py];
				Icc=cube[x] [y] ;
				Icn=cube[x] [ny];
				Inp=cube[nx][py];
				Inc=cube[nx][y] ;
				Inn=cube[nx][ny];
				//					float ixx = Inc+Ipc-2*Icc,
				//					iyy = Icn+Icp-2*Icc,
				//					ixy = 0.5f*(Ipp+Inn-Ipn-Inp);
				double ixx = (Inc+Ipc-2*Icc),
				iyy = (Icn+Icp-2*Icc),
				ixy = (Ipp+Inn-Ipn-Inp)/sqrt2;
				// update image
				cube[x][y] += ( t1*ixx + t2*ixy + t3*iyy )*dt;
				// normalize image to the original range
				if (cube[x][y] < initial_min) cube[x][y] = initial_min;
				if (cube[x][y] > initial_max) cube[x][y] = initial_max;
			}
		}
		
	}
	
	for(int ix=0;ix<NBCols;++ix){
		for(int iy=0;iy<NBRows;++iy){
			double value = (int)(cube[ix][iy]+.5);
			image[ix][iy]  = (value > 255)? 255 : (value < 0)? 0 : value;
		}
	}
	
	for(int x=0;x<NBCols;++x){
		delete [] cube[x];
	}
	delete [] cube;
}*/


void PGM::renderTriangulation(M3Matrix& image, 
							  int NbRows, 
							  int NbCols, 
							  const char*  outputFilename, 
							  vector<Triangle*> dtta)
{
  image.Reshape(NbRows,NbCols);
  image.SetValues(-1);

  double val;
  double det;
  double ba, bb , bc;
  int xmin, ymin, xmax, ymax;

//Following lines added for anisotropic diffusion 19/10/2010
	/* init fixed */
  int** fixed = new int*[NbCols];
  for(int x=0;x<NbCols;++x){
	  fixed[x] = new int[NbRows];
	  for(int y=0;y<NbRows;++y)
		{
		  fixed[x][y] = false;
	  }
  }
	//end of addition 19/10/2010
	
	
  for ( unsigned int i = 0; i < dtta.size(); i++) 
	{
		//Following lines added for anisotropic diffusion 19/10/2010
		fixed[dtta[i]->point1->x][dtta[i]->point1->y] = true;
		fixed[dtta[i]->point1->x][dtta[i]->point2->y] = true;
		fixed[dtta[i]->point1->x][dtta[i]->point3->y] = true;
		//end of addition 19/10/2010
		
    // Warning : This is true only with integer coordinates
    xmin = dtta[i]->point1->x;
    ymin = dtta[i]->point1->y;
    if ( dtta[i]->point2->x < xmin ) xmin = dtta[i]->point2->x;
    if ( dtta[i]->point3->x < xmin ) xmin = dtta[i]->point3->x;
    if ( dtta[i]->point2->y < ymin ) ymin = dtta[i]->point2->y;
    if ( dtta[i]->point3->y < ymin ) ymin = dtta[i]->point3->y;

    xmax = dtta[i]->point1->x;
    ymax = dtta[i]->point1->y;
    if ( dtta[i]->point2->x > xmax ) xmax = dtta[i]->point2->x;
    if ( dtta[i]->point3->x > xmax ) xmax = dtta[i]->point3->x;
    if ( dtta[i]->point2->y > ymax ) ymax = dtta[i]->point2->y;
    if ( dtta[i]->point3->y > ymax ) ymax = dtta[i]->point3->y;

    // det = determinant; ba, bb, bc barycentric coordinates of p
    det = (dtta[i]->point1->x*(dtta[i]->point2->y-dtta[i]->point3->y)) + 
    		(dtta[i]->point2->x*(dtta[i]->point3->y-dtta[i]->point1->y)) + 
    			(dtta[i]->point3->x*(dtta[i]->point1->y-dtta[i]->point2->y));

	// loop on all the points of the triangle
	for (int dx=xmin;dx<=xmax;dx++) {
		for (int dy=ymin;dy<=ymax;dy++) {
			
	        ba = (dx*(dtta[i]->point2->y-dtta[i]->point3->y) + dtta[i]->point2->x*(dtta[i]->point3->y-dy) + dtta[i]->point3->x*(dy-dtta[i]->point2->y))/det;
	        if(ba<0) continue;
	        bb = (dtta[i]->point1->x*(dy-dtta[i]->point3->y) + dx*(dtta[i]->point3->y-dtta[i]->point1->y) + dtta[i]->point3->x*(dtta[i]->point1->y-dy))/det;
	        if(bb<0) continue;
	        bc = (dtta[i]->point1->x*(dtta[i]->point2->y-dy) + dtta[i]->point2->x*(dy-dtta[i]->point1->y) + dx*(dtta[i]->point1->y-dtta[i]->point2->y))/det;
			if(bc<0) continue;
			// get L(f,Ty) and epsilon
			val = (ba*dtta[i]->point1->f)+(bb*dtta[i]->point2->f)+(bc*dtta[i]->point3->f);
			if (val<0.)   val=0.;
			if (val>255.) val=255.;
			// WARNING : to be tested on non square images
			// to check that NbRows and NbCols are not inverted
			image[dy][dx] = DBL_ROUND(val+0.5);
		}
    }
  }

	//anisotropic diffusion: added the 19/10/2010
	//anisotropicDiffusion(image,fixed,NbRows,NbCols,1);
	//end of addition 19/10/2010

	
  dtta.clear();
  ofstream out(outputFilename);
  image.SavePGM(out);
}





void PGM::optimize(int NbRows, 
				   int NbCols, 
				   M3Matrix& image, 
				   vector<Triangle*> dtta, 
				   vector<Point2D*> nodes)
{
  M3Matrix PointVisited;
  PointVisited.Reshape(NbRows,NbCols);
  PointVisited.SetValues(0.);

  //double eps = -0.000000000001;
  // To get the positions in the sparse matrix
  M3Matrix PointPosition;
  PointPosition.Reshape(NbRows,NbCols);
  unsigned int x,y,count,i,j,k;

  unsigned int plen = nodes.size();
  M3Matrix IPosition;
  IPosition.Reshape(plen+1,1);
  IPosition[0][0] = 0.;

  unsigned int nl;
  unsigned int maxneighbors = 0;

  count = 0;
  for (i = 0; i < nodes.size(); i++) 
	{
  	vector<Point2D*> neighbors = nodes[i]->getNeighbors();
    nl = neighbors.size();
    if (nl > maxneighbors)
      maxneighbors = nl;
    x = nodes[i]->x;
    y = nodes[i]->y;

    PointPosition[x][y] = (double)count;
    count++;
  }

  M3Matrix JPosition;
  JPosition.Reshape(plen,maxneighbors);

  count = 0;
  for (i=0;i<plen;i++) 
	{
	  vector<Point2D*> neighbors = nodes[i]->getNeighbors();
	  nl = neighbors.size();
    IPosition[i][0] = count;
    count += nl+1;

    for (j=0;j<nl;j++) 
		{
      x = neighbors[j]->x;
      y = neighbors[j]->y;
      JPosition[i][j] = PointPosition[x][y];
    }
  }
  IPosition[plen][0] = (double)count;

  M3Matrix SparseMatrix;
  SparseMatrix.Reshape(count,1);
  SparseMatrix.SetValues(0.);

  M3Matrix fS;
  fS.Reshape(plen,1);
  fS.SetValues(0.);

  // now fills the Phi matrix
  int NbN;
  int trlen = dtta.size();
  
  double ab,ac,bc;
  double aa,bb,cc;
  double fa,fb,fc;
  int xa, xb, xc, ya, yb, yc;
  for (int i=0;i<trlen;i++) 
	{
    xa = dtta[i]->point1->x;
    ya = dtta[i]->point1->y;
    xb = dtta[i]->point2->x;
    yb = dtta[i]->point2->y;
    xc = dtta[i]->point3->x;
    yc = dtta[i]->point3->y;

    ab =0.;
    ac =0.;
    bc =0.;

    fa = 0.;
    fb = 0.;
    fc = 0.;

    int ia = DBL2INT(PointPosition[xa][ya]);
    int ib = DBL2INT(PointPosition[xb][yb]);
    int ic = DBL2INT(PointPosition[xc][yc]);

    aa = 0.;
    bb = 0.;
    cc = 0.;

    int xmin = xa;
    int ymin = ya;
    if (xb<xmin) xmin = xb;
    if (xc<xmin) xmin = xc;
    if (yb<ymin) ymin = yb;
    if (yc<ymin) ymin = yc;

    int xmax = xa;
    int ymax = ya;
    if (xb>xmax) xmax = xb;
    if (xc>xmax) xmax = xc;
    if (yb>ymax) ymax = yb;
    if (yc>ymax) ymax = yc;

    double dxa = (double)xa;
    double dxb = (double)xb;
    double dxc = (double)xc;
    double dya = (double)ya;
    double dyb = (double)yb;
    double dyc = (double)yc;

    // det = determinant; ba, bb, bc barycentric coordinates of p
    double det = (dxa*(dyb-dyc)) + (dxb*(dyc-dya)) + (dxc*(dya-dyb));

    // loop on all the points of the triangle
    for (int x=xmin;x<=xmax;x++) {
      for (int y=ymin;y<=ymax;y++) {
        if (PointVisited[x][y] == 0.) {
          double dx = (double)x;
          double dy = (double)y;

          double da = (dx*(dyb-dyc) + dxb*(dyc-dy) + dxc*(dy-dyb))/det;
          double db = (dxa*(dy-dyc) + dx*(dyc-dya) + dxc*(dya-dy))/det;
          double dc = (dxa*(dyb-dy) + dxb*(dy-dya) + dx*(dya-dyb))/det;
          if (da>=0.0 && db>=0.0 && dc>=0.0) { // FIXME: >=epsilon?
            // Filling of the Gram matrix : <Phi_i,Phi_j>

            // Non diagonal elements
            ab += da*db;
            ac += da*dc;
            bc += db*dc;

            // Filling of the second member vector <Phi_i,I>
            fa += da*image[x][y];
            fb += db*image[x][y];
            fc += dc*image[x][y];

            // Diagonal elements
            aa += da*da;
            bb += db*db;
    	    cc += dc*dc;

            PointVisited[x][y] = 1.0;
	  }
       	}
      }
    }


    fS[ia][0] += fa;
    fS[ib][0] += fb;
    fS[ic][0] += fc;

    int jab = 0;
    int jac = 0;
    int jba = 0;
    int jbc = 0;
    int jca = 0;
    int jcb = 0;

    // Point a
    NbN = DBL2INT(IPosition[ia+1][0]) - DBL2INT(IPosition[ia][0]) - 1;
    for (k=0;k<NbN;k++) 
		{
      if (JPosition[ia][k] == ib)
        jab = k;
      if (JPosition[ia][k] == ic)
        jac = k;
    }

    // Point b
    NbN = DBL2INT(IPosition[ib+1][0]) - DBL2INT(IPosition[ib][0]) - 1;
    for (k=0;k<NbN;k++) {
      if (JPosition[ib][k] == ia)
        jba = k;
      if (JPosition[ib][k] == ic)
        jbc = k;
    }

    // Point c
    NbN = DBL2INT(IPosition[ic+1][0]) - DBL2INT(IPosition[ic][0]) - 1;
    for (k=0;k<NbN;k++) {
      if (JPosition[ic][k] == ia)
        jca = k;
      if (JPosition[ic][k] == ib)
        jcb = k;
    }

    // Increments SparseMatrix (which is in fact a vector)
    // diagonal elements
    SparseMatrix[DBL2INT(IPosition[ia][0])][0] += aa;
    SparseMatrix[DBL2INT(IPosition[ib][0])][0] += bb;
    SparseMatrix[DBL2INT(IPosition[ic][0])][0] += cc;

    // Other elements
    SparseMatrix[DBL2INT(IPosition[ia][0]) + jab + 1][0] += ab;
    SparseMatrix[DBL2INT(IPosition[ia][0]) + jac + 1][0] += ac;
    SparseMatrix[DBL2INT(IPosition[ib][0]) + jba + 1][0] += ab;
    SparseMatrix[DBL2INT(IPosition[ib][0]) + jbc + 1][0] += bc;
    SparseMatrix[DBL2INT(IPosition[ic][0]) + jca + 1][0] += ac;
    SparseMatrix[DBL2INT(IPosition[ic][0]) + jcb + 1][0] += bc;
  }

  /* was zum T* ist das? Raus?
  int jp;
  for (i = 0;i<IPosition.GetNbRows()-1;i++) {
    NbN = DBL2INT(IPosition[i+1][0]) - DBL2INT(IPosition[i][0]) - 1;
    for (j=0;j<NbN;j++)
      jp = (int)JPosition[i][j]; //??
  }
  */


  // Resolution of the system : PCG method
  M3Matrix X;
  X.Reshape(fS.GetNbRows(),1);
  M3Matrix r;
  r.Reshape(fS.GetNbRows(),1);
  M3Matrix p;
  p.Reshape(fS.GetNbRows(),1);
  M3Matrix q;
  q.Reshape(fS.GetNbRows(),1);
  M3Matrix s;
  s.Reshape(fS.GetNbRows(),1);

  double gamma;
  double tmp;

  // Initialization of X,r,p,s,gamma with the image values
  for (i=0;i<X.GetNbRows();i++) {
    X[i][0] = nodes[i]->f;
  }

  for (i=0;i<X.GetNbRows();i++) {
    tmp = 0.;
    NbN = DBL2INT(IPosition[i+1][0]) - DBL2INT(IPosition[i][0]);
    tmp = SparseMatrix[DBL2INT(IPosition[i][0])][0] * X[i][0];
    for (j=1;j<NbN;j++)
      tmp += SparseMatrix[DBL2INT(IPosition[i][0])+j][0] * X[DBL2INT(JPosition[i][j-1])][0];
    r[i][0] = fS[i][0] - tmp;
  }


  for (i=0;i<p.GetNbRows();i++) 
	{
    tmp = 0.;
    NbN = DBL2INT(IPosition[i+1][0]) - DBL2INT(IPosition[i][0]);
    tmp = SparseMatrix[DBL2INT(IPosition[i][0])][0] * r[i][0];
    for (j=1;j<NbN;j++)
      tmp += SparseMatrix[DBL2INT(IPosition[i][0])+j][0] * r[DBL2INT(JPosition[i][j-1])][0];
    p[i][0] = tmp;
    s[i][0] = tmp;
  }

  gamma = 0.;
  for (i=0;i<p.GetNbRows();i++)
    gamma += s[i][0]*s[i][0];

  // Iterations of the method
  unsigned int NbIt = 100;

  double alpha,beta;
  double tmp1=0.;
  double Tol = 0.0000000001;

  if (gamma>Tol) {
   for (k=0;k<NbIt;k++) {
     // q_k = Phi * p_k
     for (i=0;i<p.GetNbRows();i++) {
       tmp = 0.;
       NbN = DBL2INT(IPosition[i+1][0]) - DBL2INT(IPosition[i][0]);
       tmp = SparseMatrix[DBL2INT(IPosition[i][0])][0] * p[i][0];
       for (j=1;j<NbN;j++){
		tmp += SparseMatrix[DBL2INT(IPosition[i][0])+j][0] * p[DBL2INT(JPosition[i][j-1])][0];
       }
       q[i][0] = tmp;
     }



     // alpha_k = gamma_k/N_2(q_k)
     tmp1 = 0.;
     for (i=0;i<p.GetNbRows();i++) tmp1 += q[i][0]*q[i][0];

     alpha = gamma/tmp1;

     // X_{k+1} =X_k + alpha_k * p_k
     // r_{k+1} =r_k - alpha_k * q_k
     for (i=0;i<p.GetNbRows();i++) {
       X[i][0] = X[i][0] + alpha * p[i][0];
       r[i][0] = r[i][0] - alpha * q[i][0];
     }

     // s_{k+1}=A*r_{k+1}
     for (i=0;i<p.GetNbRows();i++) {
       tmp = 0.;
       NbN = DBL2INT(IPosition[i+1][0]) - DBL2INT(IPosition[i][0]);
       tmp = SparseMatrix[DBL2INT(IPosition[i][0])][0] * r[i][0];
       for (j=1;j<NbN;j++)
         tmp += SparseMatrix[DBL2INT(IPosition[i][0])+j][0] * r[DBL2INT(JPosition[i][j-1])][0];
       s[i][0] = tmp;
     }

     // gamma_{k+1} = N_2(s_{k+1})
     tmp1 = 0.;
     for (i=0;i<p.GetNbRows();i++) tmp1 += s[i][0]*s[i][0];

     //
     beta = tmp1/gamma;
     gamma = tmp1;

     // p_{k+1} = s_{k+1} + beta * p_k
     for (i=0;i<p.GetNbRows();i++) p[i][0] = s[i][0] + beta * p[i][0];
     if (gamma < Tol) k = NbIt;
   }
  }


  for (int i=0;i<trlen;i++) 
	{
    int ia = DBL2INT(PointPosition[dtta[i]->point1->x][dtta[i]->point1->y]);
    int ib = DBL2INT(PointPosition[dtta[i]->point2->x][dtta[i]->point2->y]);
    int ic = DBL2INT(PointPosition[dtta[i]->point3->x][dtta[i]->point3->y]);

    dtta[i]->point1->f = X[ia][0];
    dtta[i]->point2->f = X[ib][0];
    dtta[i]->point3->f = X[ic][0];
  }
  
	cout << "[End of optimize. XDump()]" << endl;
	X.Dump();
}




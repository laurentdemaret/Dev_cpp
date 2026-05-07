#include <ctime>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <fstream>
#include <sstream>
#include <set>
#include <vector>
#include <unistd.h>

#include "IO.h"
#include "matrix.h"
#include "PGM.h"
#include "PointGrid.h"
#include "Point2D.h"
//#include "Edge.h"
#include "Triangle.h"
#include "Triangulation.h"
#include "Thinning.h"
#include "coding.h"

using namespace std;

// **********************************************************************************************
// * Liste von Todos (Stand: 20.04.26)
// 1 - Check rows and cols - erledigt: 27.40.26
// 2 - Wo double Thinning::calculateError_Voronoi(Triangle* triangle, Point2D* point) einsetzen
// 3 - Berechnung des Aktualisierungsfehlers ()
// 4 - Liste der Voronoi Zelle Pixeln -> direkt in welche Klasse rein (aktuell noch nicht)
// 5 - Vereinfachung der I/O Funktionalitäten - präziser angeben
// 6 - Print (in eps Format) die Voronoi-Zellen
// **********************************************************************************************


class TestTriangle;


class TestPoint
{
public:
    TestPoint(){};
    TestPoint(double xd, double yd){
        x = xd;
        y = yd;
    };
    //virtual ~TestPoint() {}

    void Dump()
    {
        std::cout << "Dump (TestPoint) :" << x << ", " << y  << std::endl;
    }

    double x=0.;
    double y=0.;

    vector<TestTriangle*> neighbour_triangles;
};


class TestTriangle
{
public:
    TestTriangle(){};
    TestTriangle(TestPoint* a,TestPoint* b,TestPoint* c)
    {
        p1 = a;
        p2 = b;
        p3 = c;
    };
    virtual ~TestTriangle() {}

    TestPoint* p1 = nullptr;
    TestPoint* p2 = nullptr;
    TestPoint* p3 = nullptr;

    /*void SetPoints(TestPoint* a,TestPoint* b,TestPoint* c)
    {
        p1 = a;
        p2 = b;
        p3 = c;
    };*/
    double GetArea()
    {
        double p12x = p2->x-p1->x;
        double p12y = p2->y-p1->y;
        double p13x = p3->x-p1->x;
        double p13y = p3->y-p1->y;

        double area = fabs(0.5*(p12x*p13y-p12y*p13x));
        return area;
    }

};


//Classe pour tester pas vraiment utile
class ImageMesh
{
public:
    //std::vector<std::unique_ptr<TestPoint>> pixels;
    //std::vector<std::unique_ptr<TestPoint>> nodes;
    //std::vector<std::unique_ptr<TestTriangle>> triangles;
    std::vector<TestPoint*> pixels;
    std::vector<TestPoint*> nodes;
    std::vector<TestTriangle*> triangles;

    TestPoint* addPixel(double x, double y)
    {
        TestPoint* p = new TestPoint(0,0);
        pixels.push_back(p);
        //auto* p = pixels.back().get();
        p->x = x;
        p->y = y;
        return p;
    }
};


void test_Pointers()
{
    std::cout << "test Pointers " << std::endl;
    double xd = 0.4;
    double yd = 0.5;

    TestPoint* p = new TestPoint(xd,yd);
    {
        vector<TestPoint*> points;

      points.push_back(p);
      points[0]->Dump();
      std::cout << "1. p        " <<  p << std::endl;
      std::cout << "2. points[0] " <<  points[0] << std::endl;
    }
    std::cout << "3. p        " <<  p << std::endl;
    //std::cout << "4. points[0] " <<  points[0] << std::endl;

    //points[0]->Dump();

    ImageMesh m;
    m.addPixel(0.,0.);
    m.addPixel(0.,1.);
    m.addPixel(1.,1.);

    /*TestTriangle* t1 = new TestTriangle();
    TestTriangle* t2 = new TestTriangle();

    TestPoint* a = new TestPoint(0.,0.0);
    TestPoint* b = new TestPoint(0.,1.0);
    TestPoint* c = new TestPoint(1.,0.0);
    TestPoint* d = new TestPoint(1.,1.0); */

      /*t1->SetPoints(a,b,c);
      a->neighbour_triangles.push_back(t1);
      b->neighbour_triangles.push_back(t1);
      c->neighbour_triangles.push_back(t1);

      t2->SetPoints(b,c,d);
      b->neighbour_triangles.push_back(t2);
      c->neighbour_triangles.push_back(t2);
      d->neighbour_triangles.push_back(t2);*/


      //std::cout << "Area of the triangle before change: " << t1->GetArea() << std::endl;

}


// *************************************
// Test Triangulations: ajouté nov. 25
// *************************************
void test_Triangulations()
{
    std::cout << "test Triangulations" << std::endl;

    //Load an image
    M3Matrix OriginalImage;
    std::string OriginalImageName("Original/church_gray_ascii.pgm");
    //std::string OriginalImageName("Data/lena128.pgm");
    OriginalImage.LoadFromPGM(OriginalImageName);

    std::cout << "OriginalImage.GetNbRows():" <<  OriginalImage.GetNbRows() << " OriginalImage.GetNbCols():" << OriginalImage.GetNbCols() <<  std::endl;
    std::string SaveName("image_test_a.pgm");
    ofstream outStream_Save(SaveName);

    OriginalImage.Reshape(3,4);
    OriginalImage.SetValues(0.0);

    OriginalImage[0][3]=255;

    OriginalImage.SavePGM(outStream_Save);

    exit(-1);

    M3Matrix TestImage;
    std::string TestImageName("Out/test.pgm");


    TestImage = OriginalImage;
    for(int i=0;i<TestImage.GetNbRows();i++)
    {
        for(int j=0;j<TestImage.GetNbCols();j++)
        {
            TestImage[i][j] = 255-TestImage[i][j];
        }
    }

    ofstream outStream_Test(TestImageName);
    TestImage.SavePGM(outStream_Test);


    //Define some points
    vector<Point2D*> nodes;

    //Point2D* p1 = Point2D::makePoint2D(0,0,0);

    // Four corners
    nodes.push_back(Point2D::makePoint2D(0,0,OriginalImage[0][0]));
    nodes.push_back(Point2D::makePoint2D(127,0,OriginalImage[127][0]));
    nodes.push_back(Point2D::makePoint2D(0,127,OriginalImage[0][127]));
    nodes.push_back(Point2D::makePoint2D(127,127,OriginalImage[127][127]));

    for(int i = 1;i<22;i++)
        for(int j = 1;j<22;j++)
        {
            nodes.push_back(Point2D::makePoint2D(5*i,5*j,OriginalImage[5*i][5*j]));
        }
/*    nodes.push_back(Point2D::makePoint2D(10,20,OriginalImage[10][20]));
    nodes.push_back(Point2D::makePoint2D(30,70,OriginalImage[30][70]));
    nodes.push_back(Point2D::makePoint2D(10,10,OriginalImage[30][70]));
    nodes.push_back(Point2D::makePoint2D(100,60,OriginalImage[100][60]));
    nodes.push_back(Point2D::makePoint2D(40,40,OriginalImage[40][40]));
    nodes.push_back(Point2D::makePoint2D(60,60,OriginalImage[60][60])); */

    Triangulation* tri = Triangulation::makeTriangulation(nodes);
    vector<Edge*> edges = tri->getEdges();
    vector<Triangle*> triangles = tri->getTriangles(edges);


    vector<Point2D*> test_attached_points = tri->getAttachedPoints(nodes[0]);

    cout << endl;
    cout << endl;
    cout << "********************************"<< endl;
    cout << "********************************"<< endl;
    cout << "nodes[0] : " << nodes[0]->x << ", " << nodes[0]->y << ":" << nodes[0]->f  << endl;
    cout << "voici les coordonnées des points attaches à nodes[0]:"<< endl;

    for(int i = 0;i <(int) test_attached_points.size();i++)
    {
        cout << "test_attached_points["<<i<<"] : " << test_attached_points[i]->x << ", " << test_attached_points[i]->y<< ":" << test_attached_points[i]->f  << endl;

    }


    for(int i = 0;i<triangles.size();i++)
    {
        Triangle* t= triangles.at(i);
        std::cout << "Triangle :" << i<< std::endl;
        std::cout << "(" << t->v1->x << " , " << t->v1->y << ")" << std::endl;
        std::cout << "(" << t->v2->x << " , " << t->v2->y << ")" << std::endl;
        std::cout << "(" << t->v3->x << " , " << t->v3->y << ")" << std::endl;
        std::cout << std::endl;
    }


    //std::stream
    //std::string();

//    tri->writeNodes();


    // Premier test Voronoi sur un triangle et un pixel sélectionnés
    //Point2D* pixel = new Point2D(9,7,0);
    //Point2D* pixel = Point2D::makePoint2D(9,7,0);

    //std::size_t NbRows = 10;
    //std::size_t NbCols = 10;


    int NbRows = 127, NbCols = 127;
    std::vector<std::vector<Point2D*> > pixels(
        NbRows,
        std::vector<Point2D*>(NbCols, nullptr)
        );


    //vector<Point2D*> pixels;
    for(int i = 0;i < NbRows;i++)
        for(int j = 0;j < NbCols;j++)
        {
            //pixels.push_back(Point2D::makePoint2D(i,j,0));
            pixels[i][j] = new Point2D(i, j, 0);

            for(int t = 0; t < (int)triangles.size();t++)
            {
                Triangle* triangle = triangles.at(t);
                if(triangle->insideTriangle(pixels[i][j]))
                {
                    std::cout << "Le point "<< i << "," << j  <<" est dans le triangle " << t << std::endl;

                    Point2D* centroid = triangle->getVoronoiCentroid(pixels[i][j]);
                    pixels[i][j]->setVoronoiCentroid(centroid);
                    pixels[i][j]->f = centroid->f;
                    std::cout << "pixel : " << pixels[i][j] << std::endl;

                    //std::cout << " pixel->clone: " << pixels[i][j]->clone << std::endl;
                    //std::cout << " pixel->VoronoiCentroid: " << pixels[i][j]->VoronoiCentroid << std::endl;
                }
            }
        }

    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << "Test après la boucle " << std::endl;


    int i = 70;
    int j = 40;

    std::cout << "i,j" << i << "," << j << std::endl;

    int vc_x = pixels[i][j]->VoronoiCentroid->x;
    int vc_y = pixels[i][j]->VoronoiCentroid->y;
    int vc_f = pixels[i][j]->VoronoiCentroid->f;

    std::cout << "vc_x:" << vc_x << std::endl;
    std::cout << "vc_y:" << vc_y << std::endl;
    std::cout << "vc_f:" << vc_f << std::endl;
    std::cout << "f:" << pixels[i][j]->f << std::endl;


    M3Matrix Image(NbRows,NbCols);
    for(int i =0;i<Image.GetNbRows();i++)
        for(int j =0;j<Image.GetNbCols();j++)
        {
            Image[i][j] = pixels[i][j]->f;
        }

    ofstream outStream("image.pgm");
    Image.SavePGM(outStream);

    //Implémentation naive revoir l'indexation des pixels ? classe matrix ?


    /*pixels.push_back(Point2D::makePoint2D(9,7,0));
    Point2D* pixel = pixels.at(0);
    for(int i = 0;i < (int)triangles.size();i++)
    {
        Triangle* triangle = triangles.at(i);
        if(triangle->insideTriangle(pixel))
        {
           std::cout << "Le point est dans le triangle " << i << std::endl;
           Point2D* centroid = triangle->getVoronoiCentroid(pixel);
           pixel->setVoronoiCentroid(centroid);
           std::cout << "pixel : " << pixel << std::endl;
           std::cout << " pixel->clone: " << pixel->clone << std::endl;
           std::cout << " pixel->VoronoiCentroid: " << pixel->VoronoiCentroid << std::endl;

           //std::cout << "centroid : " << centroid->x  << " "  << centroid->y << std::endl;
        }
        else
        {
          std::cout << "Le point est en dehors du triangle" << i <<  std::endl;
        }
        std::cout << "SET centroid ptr dans la boucle: " << pixel->VoronoiCentroid << std::endl;
    }*/

    /*std::cout << "SET centroid ptr après la boucle: " << pixel->VoronoiCentroid << std::endl;
    std::cout << "pixel->centroid " << pixel->VoronoiCentroid->x << ", " << pixel->VoronoiCentroid->y << std::endl;*/

    //Point2D* p1 = triangle->point1;
    //std::pair<Point2D*, Point2D*>  p_op = triangle->getOtherPoints(p1);

    /*vector<Point2D*> nodes_from_tri = tri->getNodes();
    for(int i = 0;i<nodes_from_tri;i++)
    {}*/

    //Triangulation()

    std::cout << "push_back successful ?" << std::endl;

   //exit(1);
}

double my_function_2D(double x)
{
    double y = x*x;
    return y;
}


//
void test()
{
    std::cout << "This is a simple test for short testing of functionalities"  << std::endl;
    double x = 3.5;

    double y = my_function_2D(x);
    cout << "x : " << x << endl;
    cout << "y : x^2 " << y << endl;
}


//this is the main for ./xat
int main(int argc, char* argv[]) 
{
    //test();
    //test_Pointers();
    //test_Triangulations();

	int at_alg, iterations, exchange_iterations,exchange_radius,quantization, compressing_options;
	int filenames_num;
    M3Matrix image, image_voronoi;
	char* action_key;

    std::cout << "Version portable de at (20/04/2026)" << std::endl;

	char** filenames = parse_command_line_at_and_compress(argc,argv,
														&at_alg,
														&iterations,
														&exchange_iterations,
														&exchange_radius,
														&quantization,
														&compressing_options, 
														&action_key,
														&filenames_num);


    if( strcasecmp(action_key,"test") == 0 )
    {
        cout<<"[test]"<<endl;
        test_Triangulations();
        std::cout << "test_Triangulations done" << std::endl;
    }

	if ( filenames_num == 0 ) 
	{
		cerr << "Error: no file to process, aborting..." << endl;
		print_usage();
		free(filenames);
		return EXIT_FAILURE;
	}
	
	// output 
	string filename = stripExtension(filenames[0]);
	
    // action = reconstruct
    // TODO: add a typical command line
	if( strcasecmp(action_key,"reconstruct") == 0 )
	{
        cout<<"[reconstruct]"<<endl;
        // read nodes
        //TODO: which formate is expected ?
		ifstream input(filenames[0]);
		//assert(input.is_open());
		int r,c,s;
		char buf[10];
	 	input.getline(buf, 10);
	 	s = atoi(buf);
		input.unget();
		input >> r;
		input >> c;
		
		int x,y,f;
		vector<Point2D*> tnds;
		for (int j=0;j<s;j++) 
		{
			input >> x;
			input >> y;
			input >> f;
			Point2D* p = Point2D::makePoint2D(x,y,f);
			tnds.push_back(p);		
		}
		
        cout<<"[Node Count] "<<tnds.size()<<endl;
		Triangulation* ttri = Triangulation::makeTriangulation(tnds);
		ttri->nbRows = r;
		ttri->nbCols = c;
		string fn = filename+".rec.pgm";
		PGM::renderTriangulation(image, r, c, 
														 fn.c_str(), ttri->getTriangles());
        string fn_voronoi = filename+"_Voronoi.pgm";
        PGM::renderTriangulation(image_voronoi, r, c,
                                 fn_voronoi.c_str(), ttri->getTriangles());
		Thinning* tthinning = new Thinning(ttri);
		tthinning->printEdges(filename+".edge",1);
		tthinning->printTriangles(filename+".triangle",1);
		return (0);
	}
		
	clock_t start, start0;
	start  = clock();
	start0 = clock();
	PointGrid* grid = PGM::readFile(filenames[0]);

	cout<<"[read] "<<(clock()-start)<<endl;
	start  = clock();
	Triangulation* tri = new Triangulation(grid);
	cout<<"[triangulation] "<<(clock()-start)<<endl;

    std::cout <<   " tri->nbRows :" << tri->nbRows   << "  tri->nbCols : "  << tri->nbCols << std::endl;
    // *************************
    //* repousser cette borne
    // *************************

    //This is the thinning algorithm
	start  = clock();
	Thinning* thinning = new Thinning(tri);
  //thinning->lambda = 200.0; //TODO: this should be a parameter !!!!
	thinning->lambda = 0.0; //TODO: this should be a parameter !!!!

	switch (at_alg) 
	{
        //TODO: implement an algorithm with Voronoi cells
        case 9:
            thinning->thinningAT9(iterations);
            break;
        case 5:
			thinning->thinningAT5(iterations);		
			break;
		case 6:
			thinning->thinningAT6(iterations);		
			break;
		case 8:
			thinning->thinningAT8(iterations);		
			break;
		default:
			thinning->fastThinning(iterations);		
			break;
	}

    cout<<"[Thinning] "<<(clock()-start)<<endl;

	// original image 
	M3Matrix image0, optImage;
    //image0.Reshape(thinning->triangulation->nbCols, thinning->triangulation->nbRows);
	optImage.Reshape(thinning->triangulation->nbRows, thinning->triangulation->nbCols );

	vector<Triangle*> dtta = thinning->triangulation->getTriangles();
	vector<Point2D*> nodes = thinning->triangulation->getNodes();

	int row, col;
	for (unsigned int n = 0; n < thinning->triangulation->nodes.size(); n++) 
	{
      row = thinning->triangulation->nodes[n]->x;
      col = thinning->triangulation->nodes[n]->y;
      optImage[row][col] = thinning->triangulation->nodes[n]->f;
      //image0[row][col] = thinning->triangulation->nodes[n]->f;
	}

	ostringstream oss_temp;
	oss_temp << "_AT" << at_alg << "_I";
	if (iterations==-1) oss_temp << "default";
	else oss_temp << iterations;
	oss_temp <<"_Q" << quantization ;
	
	string filenameComplete = filename;
	filenameComplete.append(oss_temp.str());
		
	string filenameTriangulationEps = filename;
	filenameTriangulationEps.append(oss_temp.str());
	string filenameTriangulationOff = filename;
	filenameTriangulationOff.append(oss_temp.str());
	
	string filenameRightAfterThinning = filenameComplete+"_RightAfterThinning.pgm";
	
    filenameTriangulationEps+="_TriangulationRAT.eps";
	filenameTriangulationOff+="_TriangulationRAT.off";
	
    thinning->printTriangulationEps(filenameTriangulationEps);
	thinning->printTriangulationOff(filenameTriangulationOff);

    std::cout << "[print Triangulation functionalities done]" << std::endl;

    PGM::renderTriangulation(image, thinning->triangulation->nbRows,
                                    thinning->triangulation->nbCols,
                                    filenameRightAfterThinning.c_str(),dtta);

    std::cout << "after render Triangulation" << std::endl;
    string fn_voronoi = filename+"_Voronoi.pgm";
    PGM::renderVoronoi( image_voronoi, row, col,
                             fn_voronoi.c_str(), thinning->triangulation->getTriangles());

    std::cout << "after render Voronoi" << std::endl;
    exit(-1);

	char psnr[10];
	sprintf(psnr,"%f",image.PSNR(image0));
	string psnr_str = psnr;
	cout<<"[PSNR] RightAfterThinning: "<<image.PSNR(image0)<<endl;
	// exchange 
	start  = clock();
	// build thinned neighbor heap 
	Thinning::exchangeRadius = exchange_radius;
	thinning->exchange(exchange_iterations);
	cout<<"[Exchange] "<<(clock()-start)<<endl;
	int time = (clock()-start0)/1000;
	
	Triangle::deleteTriangles(dtta);
	dtta = thinning->triangulation->getTriangles();
	nodes = thinning->triangulation->getNodes();
	
	int min = time / 60;
	int sec = time % 60;
	cout<<"[Total Time] "<<min<<"min:"<<sec<<"s "<<endl;
	
	string filenameExchanged = filenameComplete+"_Exchanged.pgm";
	PGM::renderTriangulation(image,
													 thinning->triangulation->nbRows, 
													 thinning->triangulation->nbCols, 
													 filenameExchanged.c_str(), 
													 dtta);
	sprintf(psnr,"%f",image.PSNR(image0));
	psnr_str = psnr;
	cout<<"[PSNR] Exchanged: "<<image.PSNR(image0)<<endl;
	
	// L2 optimization 	
	PGM::optimize(thinning->triangulation->nbRows, thinning->triangulation->nbCols, optImage, dtta, nodes );
	
	string filenameOutside = filename;
	filenameOutside.append(oss_temp.str());
	filenameOutside+=".outside";
	
	string filenameInside = filename;
	filenameInside.append(oss_temp.str());
	filenameInside+=".inside";
	
	string filenameEdge = filename;
	filenameEdge.append(oss_temp.str());
	filenameEdge+=".edge";
	
	string filenameTriangle = filename;
	filenameTriangle.append(oss_temp.str());
	filenameTriangle+=".triangle";
	
	string filenameEps = filename;
	filenameEps.append(oss_temp.str());
	filenameEps+=".eps";
	
	string filenameTriangulationEpsExchanged = filename;
	filenameTriangulationEpsExchanged.append(oss_temp.str());
	filenameTriangulationEpsExchanged+="_TriangulationExchanged.eps";
	
	string filenameTriangulationOffExchanged = filename;
	filenameTriangulationOffExchanged.append(oss_temp.str());
	filenameTriangulationOffExchanged+="_TriangulationExchanged.off";
		
	//print removed nodes
	thinning->printThinnedNodes(filename+".outside",1);
	// print thinning nodes  
	thinning->printNodes(filename+".sig.inside",0);
	thinning->printNodes(filenameInside,1);
	// print edges  
	thinning->printEdges(filenameEdge,1);
	// print triangles  
	thinning->printTriangles(filenameTriangle,1);
	
	//print eps file
	thinning->printEps(filenameEps);
	thinning->printTriangulationEps(filenameTriangulationEpsExchanged);
	thinning->printTriangulationOff(filenameTriangulationOffExchanged);
	
	string filenameOptimized = filenameComplete+"_Optimized.pgm";
	
	PGM::renderTriangulation(image,thinning->triangulation->nbRows,thinning->triangulation->nbCols, 
							 filenameOptimized.c_str(),dtta);
	
	sprintf(psnr,"%f",image.PSNR(image0));
	psnr_str = psnr;
	cout<<"[PSNR] Optimized: "<<image.PSNR(image0)<<endl;
	
	// encode 
	string filenameEncoded = filename;
	filenameEncoded.append(oss_temp.str());
	filenameEncoded+=".at";
	
	Coding::Encode(thinning->triangulation->nbRows, thinning->triangulation->nbCols, thinning->triangulation->getNodes(), filenameEncoded.c_str(),8);
	string filenameDecoded = filename+"_decoded_Q8.pgm";
	vector<Point2D*> nds = Coding::Decode(&(thinning->triangulation->nbRows), &(thinning->triangulation->nbCols), filenameEncoded.c_str());
	Triangulation* tri0 = Triangulation::makeTriangulation(nds);
	PGM::renderTriangulation(image,
													 thinning->triangulation->nbRows, 
													 thinning->triangulation->nbCols, 
													 filenameDecoded.c_str(), 
													 tri0->getTriangles());
	
	sprintf(psnr,"%f",image.PSNR(image0));
	psnr_str = psnr;
	cout<<"[PSNR] Decoded: "<<image.PSNR(image0)<<endl;
	tri->deleteTriangulation(dtta);
	tri0->deleteTriangulation();
	FILE * pFile;
	pFile = fopen (filenameEncoded.c_str(),"r+");
	fseek (pFile,0,SEEK_END);
	double fs = ftell(pFile);
	double fsize = fs*8;
	fclose (pFile);
	double bpp = fsize/(tri->nbCols*tri->nbRows);
	cout<<"[Coded Size] "<<fs<<endl;
	cout<<"[BPP] "<<bpp<<endl;
	
	return (0);	
}











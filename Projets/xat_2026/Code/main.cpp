#include <time.h>
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
#include "Edge.h"
#include "Triangle.h"
#include "Triangulation.h"
#include "Thinning.h"
#include "coding.h"

using namespace std;

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

    exit(-1);
}

//this is the main for ./xat
int main(int argc, char* argv[]) 
{
    test();
	std::cout << "This is the new version 2025 of at" << std::endl;
	
	int at_alg, iterations, exchange_iterations,exchange_radius,quantization, compressing_options;
	int filenames_num;
	M3Matrix image;
	char* action_key;
	
	char** filenames = parse_command_line_at_and_compress(argc,argv,
														&at_alg,
														&iterations,
														&exchange_iterations,
														&exchange_radius,
														&quantization,
														&compressing_options, 
														&action_key,
														&filenames_num);
	
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
	if( strcasecmp(action_key,"reconstruct") == 0 )
	{
        cout<<"[reconstruct]"<<endl;
        // read nodes
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
	
  //This is the thinning algorithm
	start  = clock();
	Thinning* thinning = new Thinning(tri);
  //thinning->lambda = 200.0; //TODO: this should be a parameter !!!!
	thinning->lambda = 0.0; //TODO: this should be a parameter !!!!
	
	switch (at_alg) 
	{
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
	image0.Reshape(thinning->triangulation->nbCols, thinning->triangulation->nbRows);
	optImage.Reshape(thinning->triangulation->nbRows, thinning->triangulation->nbCols );
	
	vector<Triangle*> dtta = thinning->triangulation->getTriangles();
	vector<Point2D*> nodes = thinning->triangulation->getNodes();
	
	int row, col;
	for (unsigned int n = 0; n < thinning->triangulation->nodes.size(); n++) 
	{
	  col = thinning->triangulation->nodes[n]->x;
	  row = thinning->triangulation->nodes[n]->y;
	  optImage[col][row] = thinning->triangulation->nodes[n]->f;
	  image0[row][col] = thinning->triangulation->nodes[n]->f;
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
	
	PGM::renderTriangulation(image,
													 thinning->triangulation->nbRows, 
													 thinning->triangulation->nbCols, 
													 filenameRightAfterThinning.c_str(),dtta);
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











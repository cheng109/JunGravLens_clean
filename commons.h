/*
 * commons.h
 *
 *  Created on: Sep 30, 2015
 *      Author: cheng109
 */

#ifndef COMMONS_H_
#define COMMONS_H_


#include <string>
#include <vector>
#include "Image.h"
#include "nanoflann.hpp"
//#include "Model.h"
#include <map>
#include <Eigen/Sparse>

using namespace std;
using namespace nanoflann; 

typedef Eigen::SparseMatrix<double> sp_mat;
typedef Eigen::SparseMatrix<double,Eigen::RowMajor> sp_mat_row; 
typedef Eigen::VectorXd vec;

struct normVec{ 
	double n0;  
	double n1; 
	double n2; 
	normVec() {};
	normVec(double n0, double n1, double n2):n0(n0), n1(n1), n2(n2) {} ;
}; 

struct Point{
	double x; 
	double y; 
	double z;
	Point(double a, double b, double c):x(a), y(b), z(c) {};
}; 


template <typename T>
struct PointCloud
{
	struct Point
	{
		T  x,y,z;
	};

	vector<Point>  pts;

	// Must return the number of data points
	inline size_t kdtree_get_point_count() const { return pts.size(); }

	// Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
	inline T kdtree_distance(const T *p1, const size_t idx_p2,size_t /*size*/) const
	{
		const T d0=p1[0]-pts[idx_p2].x;
		const T d1=p1[1]-pts[idx_p2].y;
		//const T d2=p1[2]-pts[idx_p2].z;
		return d0*d0+d1*d1 ; 
	}

	inline T kdtree_get_pt(const size_t idx, int dim) const
	{
		if (dim==0) return pts[idx].x;
		else if (dim==1) return pts[idx].y;
		else return pts[idx].z;
	}
	template <class BBOX>
	bool kdtree_get_bbox(BBOX& /*bb*/) const { return false; }
};






class Conf{
public:
	// Cosmic Constants:
	double omega;
	double lambda;
	double weos;

	double hubble;
	double srcZ;
	double lenZ;

	int causticLevel; 

	double length;
	size_t srcSize[2];
	size_t imgSize[2];
	size_t potSize[2];
	double srcRes;
	double imgRes;
	double potRes;
	double srcXCenter;
	double srcYCenter;
	double imgXCenter;
	double imgYCenter;
	double potXCenter;
	double potYCenter;

	// Flags: 
	int verbose; 
	int usingRegion; 

	int outputImages; 
	int	srcBackground ; 

	double back_mean; 
	double back_std; 
	int numSources; 

	double srcRegLevel; 
	string srcRegType; 

	long potN;
	int bitpix;
	double beta; 

	string imageFileName; 
	string criticalName; 
	string causticName; 
	string contourCritName; 
	string contourCausName; 

	Conf(Image* image, map<string, string> confMap);
	void printConfList();
};


inline double dist(Point A, Point B) ;
inline double area(Point A, Point B, Point C);
vector<double> getTriWeight(Point A, Point B, Point C, Point P);

void printerror( int status);

bool pnpoly(size_t nvert, vector<double> *vertx, vector<double> *verty, double testx, double testy);
void updateConf(string confFileName);
vector<string> parseRegionFile(string regionFileName, vector<vector<double> > *xpos, vector<vector<double> > *ypos); 

map<string, string> parseConfigure(string confFileName) ;
//double getPenalty(sp_mat* M, vec* r, vec* d, sp_mat* C);
double lm_arctanh(double x);
normVec getNormVector(Point A, Point B, Point C);
void getLinearInterpolate(Point A, Point B,  Point C,  Point *P,  char direction);
vector<double> getPentWeigth(Point A, Point B, Point C, Point D, Point E);
double getEisteinRadius(Conf* conf, double Mtot);
void getAngularSizeDistance(Conf* conf, double z, double* comoveD, double* angularD) ;
int sign(double x);
normVec meanNormVector(vector<normVec>  normList);
sp_mat generatePSFoperator(string psfFileName, int axis1, int axis2);
sp_mat generatePSFoperator(string psfFileName, Image* image);
vector<double> eigenV_to_cV(vec * v) ;
vec cV_to_eigenV(vector<double>* s);
vector<string> splitString(string s);
double	lm_arccosh(double x); 
double	lm_nfw_mass(double x); 

double getMassLuminosity(Image* lensImage, Image* dataImage,  string regionFileName, double background); 


Image * magDiffMap(string img1FileName, string img2FileName, double back1, double back2,
					double std1, double std2,  string regionFile , int pixelCombine); 
//vector<double> getCritCaustic(Conf* conf, Model* model); 



#endif

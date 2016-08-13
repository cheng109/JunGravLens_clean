/*
 * Model.h
 *
 *  Created on: Oct 9, 2015
 *      Author: juncheng
 */

#ifndef MODEL_H_
#define MODEL_H_
#include <string>
#include <vector>
#include "commons.h"
#include "Image.h"
#include <map>
#include <Eigen/Sparse>
#include <fstream>
#include <sstream>
//#include <boost/algorithm/string.hpp>
#define NUM_PTMASS_PARAM 	3
#define NUM_SIE_PARAM 		5
#define NUM_NFW_PARAM 		6
#define NUM_SPEMD_PARAM 	6
#define NUM_SERSIC_PARAM 	7

typedef Eigen::SparseMatrix<double> sp_mat;
typedef Eigen::SparseMatrix<double,Eigen::RowMajor> sp_mat_row; 
typedef Eigen::VectorXd vec;

using namespace std;

struct mixModels {
	string name; 
	vector<double> paraList; 
	mixModels(string name): name(name), paraList(8, 0.0) {

	}
}; 



struct SingleModelParam {
	string name;
	double mass;
	double nParam;
	// All shared parameters: 
	double centerX, centerXFrom, centerXTo, centerXInc;
	double centerY, centerYFrom, centerYTo, centerYInc;
	
	double e, eFrom, eTo, eInc;
	double q, qFrom, qTo, qInc;
	double PA, PAFrom, PATo, PAInc;

	// For PTMASS model and SIE model: 
	double critRad, critRadFrom, critRadTo, critRadInc;
	double power, powerFrom, powerTo, powerInc; 
	double core, coreFrom, coreTo, coreInc;
	// For NFW model
	double massScale, massScaleFrom, massScaleTo, massScaleInc;   
	double radScale, radScaleFrom, radScaleTo, radScaleInc; 	
	// For sersic model: 
	double kap, kapFrom, kapTo, kapInc; 
	double sersicScale, sersicScaleFrom, sersicScaleTo, sersicScaleInc; 
	double m, mFrom, mTo, mInc; 

	


	SingleModelParam() { 
		name = "none"; 
		mass = 0;  nParam = 0; 
		centerX = centerXFrom = centerXTo = centerXInc = 0 ;
		centerY = centerYFrom = centerYTo = centerYInc = 0 ;

		e =  eFrom =  eTo =  eInc = 0 ; 
	 	q =  qFrom =  qTo =  qInc = 0 ; 
	 	PA =  PAFrom =  PATo =  PAInc = 0 ; 

	
	 	critRad =  critRadFrom =  critRadTo =  critRadInc = 0 ; 
	 	power =  powerFrom =  powerTo =  powerInc = 0 ;  
	 	core =  coreFrom =  coreTo =  coreInc = 0 ; 
	
	 	massScale =  massScaleFrom =  massScaleTo =  massScaleInc = 0 ;    
	 	radScale =  radScaleFrom =  radScaleTo =  radScaleInc = 0 ;  	
	
	 	kap =  kapFrom =  kapTo =  kapInc = 0 ;  
	 	sersicScale =  sersicScaleFrom =  sersicScaleTo =  sersicScaleInc = 0 ;  
	 	m =  mFrom =  mTo =  mInc = 0 ;  
	}
	bool reachEnd() { // only for SIE; 
		return critRad >= critRadTo and centerX>=centerXTo and centerY >= centerYTo
				and e>=eTo and PA >= PATo and coreTo>=coreTo; 
	}
};

class MultModelParam{


public:
	vector<SingleModelParam> parameter;  // not huge one; 
	int nLens ;
	long long int nComb ; 
	vector<int> nParam;
	vector<vector<mixModels> > mixAllModels; 

	MultModelParam(map<string,string> confMap) ; 
	void printModels() ; 
	void mix(int opt);
	vector<string> printCurrentModels(int curr); 

};

class Model {
	


public:
	int length;
	int k;   // number of sources; 
	
	// For regularization;
	double chi2;
	double srcR;
	double phiR;
	double penalty;
	// Number of lens models;
	double nLens;
	double totalParam;

	double occupation; 

	MultModelParam param;

	vector<double> srcPosXList;	  // Source position after deflection in X direction, in arcsec;
	vector<double> srcPosYList;	  // Source position after deflection in Y direction, in arcsec;

	vector<double> srcPosXListPixel;	  // Source position after deflection in X direction, in pixel;
	vector<double> srcPosYListPixel;	  // Source position after deflection in Y direction, in pixel;

	vector<int> group; 


	vector<double> srcPosXListPixel_increase; 
	vector<double> srcPosYListPixel_follow; 
	
	vector<double> srcPosXListPixel_follow; 
	vector<double> srcPosYListPixel_increase; 





	vector<double> pDeltaX;  	// Deflection angle in X direction;
	vector<double> pDeltaY; 	// Deflection angle in Y direction;
	vector<double> invMag;
	vector<double> dSy1; 
	vector<double> dSy2; 

	vector<double> res_img;    // Residual brightness
	vector<double> res_full_img;  //  the full residual map; 
	vector<double> simple_res_img;


	vector<double> mod_img;
	vector<double> critical;


	//vector<double> s;

	vector<double> something;
	sp_mat_row L;
	sp_mat M;
	vec r;
	vec new_r;
	vec s;

	vec phi;
	//vec square_s;
	//
	sp_mat Ds;
	sp_mat Dphi;


	sp_mat_row Hs1;
	sp_mat_row Hs2;

	sp_mat_row Hessian_grad; 

	sp_mat H_zero; 
	sp_mat H_grad; 
	sp_mat H_curv; 




	sp_mat H0H; 

	sp_mat Hphi;
	sp_mat HtH;
	sp_mat HphiH;
	sp_mat T;

	sp_mat RtR;

	sp_mat H0;  // zeroth-order regularization;
	sp_mat H1;  // gradient-order regularization;
	sp_mat H2;  // curvature-order regularization;
	double lambdaS;
	double lambdaC; 
	double lambdaPhi;


	sp_mat REG; 

	vector<vector<normVec> > normV;
	vector<normVec> meanNormV;
	map<pair<int, int>,int> posMap;



typedef	struct	_lmDeflCache {
	short			iType;
	double			param[3];		/* axratio, scale, r_power */
	double			pix_size;
	double			dimension[2];
	double			*pValsX;
	double			*pValsY;
	double			*pDeflX;
	double			*pDeflY;
	struct  _lmDeflCache *pNext;
}	lmDeflCache;


public:
	Model();
	Model(Conf* conList, MultModelParam multModelParam, double lambdaS);
	void update_H_zero(Conf* conf); 
	void update_H_grad(Conf* conf); 
	void update_H_curv(Conf* conf); 

	static vector<double> getDeflectionAngle(Conf* conList, double pfX, double pfY, double *pDeltaX, double *pDeltaY,  MultModelParam * param);
	void updatePosMapping(Image* image,  Conf* conList);
	void updateLensAndRegularMatrix(Image* dataImage,  Conf* constList);
	void updateGradient(Image* dataImage);
	void updateSource(Conf* conf); 
	void Logging(Image* dataImage, Conf* conList, string outFileName);
	void updateRegularMatrix();
	void solveSource( sp_mat* invC, vec* d, string R_type);

	void writeSrcImage(string outFileName, Conf* conList);
	//void updateCritCaustic(Image* dataImage,  Conf* constList);
	virtual ~Model();
	Image* getFullResidual(Image* dataImage);
	double getRegularizationSrcValue (vec d);


	double getScatterReg(); 
	double getKmeansScatter(Image* dataImage) ; 
	double getZerothOrderReg  	(Conf* conf, vector<double> briList);
	double getGradientOrderReg	(Conf* conf, vector<double> briList); 
	double getCurvatureOrderReg	(Conf* conf, vector<double> briList); 

	vector<vector<double> > getCritCausticFine(vector<double> xPosListArc, vector<double> yPosListArc, Conf* conf, MultModelParam * param, int level); 

	//vector<string> &split(string &s, char delim, vector<string> &elems) ; 

	void resetVectors(Conf* conf);
    void copyParam(Conf* conf, int i);
    void copyParam(int i1, int i2);
    void updateVegettiRegularization(); 

};

vector<Image* > getCritCaustic(Conf* conf, MultModelParam * param); 
void createDs9Contour(vector<double>* xList, vector<double>* yList, double level, string contourCritName) ; 
Image* createLensImage(Conf* conf, MultModelParam * param) ; 
void writeSrcModResImage(Model* model, Image* dataImage, Conf* conf, string fileName, string dir) ; 

#endif /* MODEL_H_ */

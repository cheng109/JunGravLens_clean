#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <map>
#include <ctime>
#include "Image.h"
#include "commons.h"
#include "fitsio.h"
#include "Model.h"
#include "parafit.h"
//#include "gsl/gsl_multimin.h"
//#include <armadillo>
#include <tuple>
#include <map>
//#include "QuadTree.h"
#include <Eigen/SparseCholesky>
#include <Eigen/Dense>

using namespace std;
//using namespace arma;


int main(int argc, char* argv[]) {


	
	if (argc<4)  { 
		cout << "***************Usage***************" << endl; 
		cout << "3 arguements are required. For example: " << endl; 
		cout << "./junGL horseshoe_test/ conf.txt output.txt" << endl; 
		return 1; 
	}


	// Assign arguments: 

	string dir 		= argv[1]; 
	string strConf 	= argv[2]; 
	string output 	= argv[3]; 
    int opt(0);
    if (argc > 4) opt = atoi(argv[4]);          //0: grid, 1: mcmc

	/***** prepare*****/
	string confFileName = dir+ strConf; 
	cout << confFileName << endl; 
	map<string, string> mapConf = parseConfigure(confFileName);

	
	cout << mapConf["imageFileName"] << endl; 

	Image* dataImage = new Image(mapConf["imageFileName"]);
	dataImage->updateFilterImage(mapConf["regionFileName"], stoi(mapConf["usingRegion"]));


	
	Conf *conf = new Conf(dataImage, mapConf);
	dataImage->updateVarList(90, conf->back_mean, conf->back_std); // (threshold, var);
	dataImage->updateBackSubtract(conf->back_mean, conf->back_std); 
	dataImage->updateGridPointType();
	
	dataImage->invC = dataImage->getVarMatrix();
	
 
	conf->printConfList();

	//cout << "start " << endl; 
	MultModelParam param = MultModelParam(mapConf);
	
	param.printModels();
	//param.mix(opt); 	 // update 'AllMixModels'; 

	// cout << "=====================" << endl; 
	// for(int i=0; i<11; ++i) {
	// 	string filterName = "f110" ; 
	// 	string endfix = "_set06"; 
	// 	double factor =0 ; 

	// 	map<string, double> backMap, backStdMap; 
	// 	backStdMap["f475"] = 6.4  ;
	// 	backStdMap["f606"] = 5.93 ; 
	// 	backStdMap["f814"] = 6.53 ; 
	// 	backStdMap["f110"] = 14.5 ; 
	// 	backStdMap["f160"] = 14.8 ; 

	// 	backMap["f475"] = 69.6  + factor * backStdMap["f475"] ;  // 6.4
	// 	backMap["f606"] = 68.2  + factor * backStdMap["f475"];  // 5.93
	// 	backMap["f814"] = 84.2  + factor * backStdMap["f814"];  // 6.53
	// 	backMap["f110"] = 665.5 + factor * backStdMap["f110"];  // 14.5
	// 	backMap["f160"] = 644.5 + factor * backStdMap["f160"];  // 14.8

		

	
	// 	string prefix = "regFiles/reg_"; 
		
	// 	map<string, string> regFileName; 
	// 	regFileName["f475"] = prefix + to_string(i) + endfix + ".reg"; 
	// 	regFileName["f606"] = prefix + to_string(i) + endfix + ".reg"; 
	// 	regFileName["f814"] = prefix + to_string(i) + endfix + ".reg"; 
	// 	regFileName["f110"] = prefix + to_string(i) + endfix + "_image_IR.reg"; 
	// 	regFileName["f160"] = prefix + to_string(i) + endfix + "_image_IR.reg"; 


		
	// 	Image* lensImage1 = new Image("horseshoe_test/color_test/f110_clean.fits"); 
	// 	string imgName = "galfit_work/" + filterName + "_ADU.fits"; 
	// 	Image* dataImage1 = new Image(imgName);
	// 	double luminosity =  getMassLuminosity(lensImage1, dataImage1, regFileName[filterName], backMap[filterName]) ; 
	// 	cout << luminosity << "\n"; 
	// } 
	// cout << endl; 
		if(0) {
		map<string, double> backMap, backStdMap; 
		backStdMap["f475"] = 6.4  ;
		backStdMap["f606"] = 5.93 ; 
		backStdMap["f814"] = 6.53 ; 
		backStdMap["f110"] = 14.5 ; 
		backStdMap["f160"] = 14.8 ; 

		double factor =0 ; 
		backMap["f475"] = 69.6  + factor * backStdMap["f475"] ;  // 6.4
		backMap["f606"] = 68.2  + factor * backStdMap["f475"];  // 5.93
		backMap["f814"] = 84.2  + factor * backStdMap["f814"];  // 6.53
		backMap["f110"] = 665.5 + factor * backStdMap["f110"];  // 14.5
		backMap["f160"] = 644.5 + factor * backStdMap["f160"];  // 14.8

		backStdMap["src"] = 0; 
		backMap["src"] = 0; 

		string img1Name = "galfit_work/clean_img/f475_clean.fits";
		string img2Name = "galfit_work/clean_img/f606_clean.fits";
		string img3Name = "galfit_work/clean_img/f814_clean.fits";
		string img4Name = "galfit_work/clean_img/f110_clean.fits";
		string img5Name = "galfit_work/clean_img/f160_clean.fits";

		vector<Image*> mapList1; 
		vector<Image*> mapList2; 
		int numRegion = 13; 
		// for(int i=0; i<numRegion; ++i) {

		// 	string regionFile = "new_region/part" + to_string(i) +".reg"; 
		// 	mapList1.push_back(magDiffMap(img1Name, img2Name, backMap["f475"], backMap["f606"],backStdMap["f475"], backStdMap["f606"], regionFile)); 
		// 	mapList2.push_back(magDiffMap(img2Name, img3Name, backMap["f606"], backMap["f814"],backStdMap["f606"], backStdMap["f814"], regionFile)); 
		// }
		
		// for(int i=0; i<mapList1.size(); ++i) {
		// 	ofstream myfile; 
		// 	myfile.open("map_" + to_string(i) + ".txt"); 
		// 	for(int j = 0; j<mapList1[i]->dataList.size(); ++j){
		// 		myfile << mapList1[i]->dataList[j] << "\t" << mapList2[i]->dataList[j] << endl;
		// 	}
		// 	myfile.close(); 
		// }

		

		int pixelCombine = 4;    // combine 3x3 pixels together; 

		Image* img0 =  magDiffMap("f475_source.fits", "f606_source.fits", backMap["src"], backMap["src"],backStdMap["src"], backStdMap["src"], "whatever", pixelCombine); 

		Image* img1 =  magDiffMap(img1Name, img2Name, backMap["f475"], backMap["f606"],backStdMap["f475"], backStdMap["f606"], "whatever", pixelCombine); 
		Image* img2 =  magDiffMap(img2Name, img3Name, backMap["f606"], backMap["f814"],backStdMap["f606"], backStdMap["f814"], "whatever", pixelCombine); 
		Image* img3 =  magDiffMap(img4Name, img5Name, backMap["f110"], backMap["f160"],backStdMap["f110"], backStdMap["f160"], "whatever", pixelCombine); 
		Image* img4 =  magDiffMap(img1Name, img3Name, backMap["f475"], backMap["f814"],backStdMap["f475"], backStdMap["f814"], "whatever", pixelCombine); 

		
		img0->writeFilterImage("img0_color.fits"); 

		img1->writeFilterImage("img1_color.fits"); 
		img2->writeFilterImage("img2_color.fits"); 
		img3->writeFilterImage("img3_color.fits"); 
		img4->writeFilterImage("img4_color.fits"); 

		// create  Mask image1;  

		Image* imgLeft = new Image(*img1); 
		Image* imgRight = new Image(*img1); 

		for(int i=0; i< img1->dataList.size(); ++i) {
			double val = img1->dataList[i]; 
			if(val  < -0.5 and val > -1.0 ) {
			 	imgLeft->dataList[i] = 0.0; 
			 	imgRight->dataList[i] = 1.0; 
			}

			else if (val  < 0.0 and val > -0.5 ) {
				imgLeft->dataList[i] = 1.0; 
				imgRight->dataList[i] = 0.0; 
			}
			else {
				imgLeft->dataList[i] = 0.0; 
				imgRight->dataList[i] = 0.0; 
			}


		}

		imgLeft->writeFilterImage("left.fits"); 

		// now we have two 'mask' images ; 

		Image* newDataImage = new Image(mapConf["imageFileName"]); 
		newDataImage->updateFilterImage(mapConf["regionFileName"], 0); // no filter at all; 

		Image* maskedLeft = new Image(*newDataImage); 
		Image* maskedRight = new Image(*newDataImage);
		maskedLeft->multiple(imgLeft); 
		maskedLeft->writeFilterImage("maskedLeft.fits"); 
		maskedRight->multiple(imgRight); 
		maskedRight->writeFilterImage("maskedRight.fits"); 



		delete maskedLeft, maskedRight, imgLeft, imgRight; 
		delete img0, img1, img2, img3, img4; 
	}

	// Inititate dataImage1 ; 
	Image* dataImage1 = new Image(mapConf["imageFileName"]);
	dataImage1->updateFilterImage("horseshoe_test/north1_points_sparse_extended.reg", 1);
	dataImage1->updateBackSubtract(conf->back_mean, conf->back_std); 

	Image* dataImage2 = new Image(mapConf["imageFileName"]);
	dataImage2->updateFilterImage("horseshoe_test/north2_points_sparse.reg", 1);
	dataImage2->updateBackSubtract(conf->back_mean, conf->back_std); 

	Image* dataImageSmall = new Image(mapConf["imageFileName"]);
	dataImageSmall->updateFilterImage("horseshoe_test/trial_arc.reg", 1);
	dataImageSmall->updateBackSubtract(conf->back_mean, conf->back_std); 


	vector<Image* > dataImageList; 
	//dataImageList.push_back(dataImage1); 
	//dataImageList.push_back(dataImage2); 

	dataImageList.push_back(dataImage);
	//dataImageList.push_back(dataImageSmall); 
	gridSearchVegetti(conf, param,  dataImageList, dir, output);	


	delete conf;
	delete dataImage, dataImage1, dataImage2;

	return 0;
}






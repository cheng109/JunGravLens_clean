/*
 * parafit.cpp
 *
 *  Created on: Dec 24, 2015
 *      Author: cheng109
 */


//#include "gsl/gsl_multimin.h"
#include "Model.h"
#include "Image.h"
#include <iostream>
#include "commons.h"
#include <fstream>
#include <limits>
#include <ctime> 
#include "parafit.h"

using namespace std;

void gridSearchVegetti(Conf* conf, MultModelParam param_old, vector<Image*> dataImageList, string dir, string outputFileName) {
	double lambdaS = conf->srcRegLevel;  


	conf->length = dataImageList[0]->length; 
	Model *model = new Model(conf, param_old, lambdaS);
	//cout << "nLen: " << param_old.parameter.size() << endl; 
		
	int minIndex = 0; 
	double minPenalty = std::numeric_limits<double>::max(); 

	ofstream output; 
	output.open(outputFileName); 



	vector<vector<double>> minParam(3, vector<double> (6, 0.0)); 
	long long count = 0; 
	vector<SingleModelParam> parameter = model->param.parameter; // make a copy; 
	//for(int i=0 ; i< model->param.nComb ; ++i) {
	//for(int i=0; i<3; ++i) {   //  3 models combine; 

		//model->length = dataImageList[0]->data.size(); 
		
		//model->updateReserve(conf); 

		
		clock_t begin = clock(); 
		model->param.parameter[0].name = "SIE"; 
		model->param.parameter[1].name = "SIE"; 
		model->param.parameter[2].name = "SIE"; 
		for (double critRad0 = parameter[0].critRadFrom;	critRad0 <= parameter[0].critRadTo; critRad0 += parameter[0].critRadInc) {
	        for (double centerX0 = parameter[0].centerXFrom;	centerX0 <= parameter[0].centerXTo; centerX0 += parameter[0].centerXInc) {
	        	for (double centerY0 = parameter[0].centerYFrom;	centerY0 <= parameter[0].centerYTo; centerY0 += parameter[0].centerYInc) {
	        		for (double e0 = parameter[0].eFrom;	e0 <= parameter[0].eTo; e0 += parameter[0].eInc) {
	        			for (double PA0 = parameter[0].PAFrom; PA0 <= parameter[0].PATo; PA0 += parameter[0].PAInc) {
	        				for(double core0 = parameter[0].coreFrom; core0 <= parameter[0].coreTo; core0 += parameter[0].coreInc) {


		for (double critRad1 = parameter[1].critRadFrom;	critRad1 <= parameter[1].critRadTo; critRad1 += parameter[1].critRadInc) {
	        for (double centerX1 = parameter[1].centerXFrom;	centerX1 <= parameter[1].centerXTo; centerX1 += parameter[1].centerXInc) {
	        	for (double centerY1 = parameter[1].centerYFrom;	centerY1 <= parameter[1].centerYTo; centerY1 += parameter[1].centerYInc) {
	        		for (double e1 = parameter[1].eFrom;	e1 <= parameter[1].eTo; e1 += parameter[1].eInc) {
	        			for (double PA1 = parameter[1].PAFrom; PA1 <= parameter[1].PATo; PA1 += parameter[1].PAInc) {
	        				for(double core1 = parameter[1].coreFrom; core1 <= parameter[1].coreTo; core1 += parameter[1].coreInc) {

	    

	    for (double critRad2 = parameter[2].critRadFrom;	critRad2 <= parameter[2].critRadTo; critRad2 += parameter[2].critRadInc) {
	        for (double centerX2 = parameter[2].centerXFrom;	centerX2 <= parameter[2].centerXTo; centerX2 += parameter[2].centerXInc) {
	        	for (double centerY2 = parameter[2].centerYFrom;  centerY2 <= parameter[2].centerYTo; centerY2 += parameter[2].centerYInc) {
	        		for (double e2 = parameter[2].eFrom;	e2 <= parameter[2].eTo; e2 += parameter[2].eInc) {
	        			for (double PA2 = parameter[2].PAFrom; PA2 <= parameter[2].PATo; PA2 += parameter[2].PAInc) {
	        				for(double core2 = parameter[2].coreFrom; core2 <= parameter[2].coreTo; core2 += parameter[2].coreInc) {
	        					
	        					

	        					model->resetVectors (conf); 
	        					//SingleModelParam s0, s1, s2; 

		 						
	        					model->param.parameter[0].critRad = critRad0;
	        					model->param.parameter[0].centerX = centerX0;
	        					model->param.parameter[0].centerY = centerY0;
	        					model->param.parameter[0].e 	  = e0;
	        					model->param.parameter[0].PA 	  = PA0;
	        					model->param.parameter[0].core    = core0;  

	        					model->param.parameter[1].critRad = critRad1;
	        					model->param.parameter[1].centerX = centerX1;
	        					model->param.parameter[1].centerY = centerY1;
	        					model->param.parameter[1].e 	  = e1;
	        					model->param.parameter[1].PA 	  = PA1;
	        					model->param.parameter[1].core    = core1;  

	        					model->param.parameter[2].critRad = critRad2;
	        					model->param.parameter[2].centerX = centerX2;
	        					model->param.parameter[2].centerY = centerY2;
	        					model->param.parameter[2].e 	  = e2;
	        					model->param.parameter[2].PA 	  = PA2;
	        					model->param.parameter[2].core    = core2;  

	        					
	        					double sum = 0 ; 
								for(auto & dataImage: dataImageList) {
									conf->length = dataImage->length; 
									model->length = dataImage->data.size(); 
									sum += getPenalty(model,  dataImage, conf)[2] ; 
								}
								if(sum < minPenalty) {
									minPenalty = sum; 
									minIndex = count; 
									minParam[0][0]= critRad0;
	        						minParam[0][1]= centerX0;
	        						minParam[0][2]= centerY0;
	        						minParam[0][3]= e0;
	        						minParam[0][4]= PA0;
	        						minParam[0][5]= core0;  

	        						minParam[1][0]= critRad1;
	        						minParam[1][1]= centerX1;
	        						minParam[1][2]= centerY1;
	        						minParam[1][3]= e1;
	        						minParam[1][4]= PA1;
	        						minParam[1][5]= core1;  

	        						minParam[2][0]= critRad2;
	        						minParam[2][1]= centerX2;
	        						minParam[2][2]= centerY2;
	        						minParam[2][3]= e2;
	        						minParam[2][4]= PA2;
	        						minParam[2][5]= core2;  
	        	
								}

								count ++; 
								cout << "[" + to_string(count) + "/" + to_string(model->param.nComb) + "] " << sum << endl; 
								if(conf->outputImages)
									writeSrcModResImage(model,dataImageList[0], conf, to_string(0), dir) ; 

	        				}
	        			}
	        		}
	        	}
	        }
	    }
	    					}
	        			}
	        		}
	        	}
	        }
	    }
	    					}
	        			}
	        		}
	        	}
	        }
	    }

	    
	    delete model; 
	


	// print the best solutions; 
	cout << "Best Models: [" << minIndex+1 << "/" + to_string(model->param.nComb) + "] " << endl; 
	for(int i=0 ;i<3; ++i) {
		cout << "Model[" << i << "]: ";  
		for(int j=0; j<6; ++j) 
			cout << minParam[i][j] << "\t" ; 
		cout << endl; 


	}


}

vector<double> getPenalty(Model* model, Image* dataImage, Conf* conf) {
	int begin = clock(); 
	vector<double> penalty(3); 
	vec d = cV_to_eigenV (&dataImage->dataList); 
	model->updatePosMapping(dataImage, conf);  // time used: 0.03s; 
	model->updateCompactMatrix(dataImage);  // 

	if(1) {
		model->update_H_zero(conf); 
		model->updateLensAndRegularMatrix(dataImage, conf);  // get matrix 'L' and 'RTR'; most time consuming part; 
		model->solveSource(&dataImage->invC, &dataImage->d, conf->srcRegType); 

		vec &s = model->s; 
		vec res = ( model->L * s - dataImage->d) ; 
		vec chi2 = res.transpose() *  dataImage->invC * res * model->lambdaC* model->lambdaC  ; 
		vec srcR = s  .transpose() *  model->REG      * s   * model->lambdaS* model->lambdaS  ; 
		penalty[0] = chi2[0]; 
		penalty[1] = srcR[0]; 
		penalty[2] = chi2[0] + srcR[0]; 
	}
	
	//penalty[2] = model->getScatterReg() ; 	
	//penalty[2] = model->getKmeansScatter(dataImage) ;
	//cout << "Scatter used : " << double(clock()-begin)/CLOCKS_PER_SEC << " seconds !" << endl; 
	return penalty; 
}


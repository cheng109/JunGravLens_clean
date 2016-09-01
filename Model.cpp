/*
 * Model.cpp
 *
 *  Created on: Oct 9, 2015
 *      Author: juncheng
 */

#include "Model.h"
#include "commons.h"
//#include <armadillo>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <ctime>
#include <cmath>
#include <fstream>
#include "fastell.h"
#include "nanoflann.hpp"
#include <Eigen/SparseCholesky>
#include <Eigen/Dense>
#include <assert.h>
#include <stdlib.h>
#include <random>
#include "Kmeans.h"
#include <memory>
using namespace std;
//using namespace arma;
using namespace Eigen;
using namespace nanoflann; 
#define DIFF 0.000000001



/* Model::Model() {
	// TODO Auto-generated constructor stub
} */

Model::Model(Conf* conf, MultModelParam param, double lambdaS):
		length(conf->length),
		k(conf->numSources),
		Kmeans_group(conf->length, 0),
		Kmeans_centers(conf->numSources), 
		param(param),
		beta(conf->beta),
		mod_img(length),
		L(length,length),
		r(2*length),
		new_r(2*length),
		s(length),
		phi(length),
		Ds(length, 2*length),
		Dphi(2*length,length),
		Hs1(length, length),
		Hs2(length, length),
		Hessian_grad(conf->srcSize[0]*conf->srcSize[1], conf->srcSize[0]*conf->srcSize[1]),
		H_zero(conf->srcSize[0]*conf->srcSize[1], length),
		H_grad(conf->srcSize[0]*conf->srcSize[1], length),
		H_curv(conf->srcSize[0]*conf->srcSize[1], length),
		H0H(length, length),
		Hphi(length, length),
		HtH(length, length),
		HcH(length, length),
		HphiH(length, length),
		T(conf->srcSize[0]*conf->srcSize[1], length),
		RtR(2*length, 2*length),
		H0(conf->srcSize[0]*conf->srcSize[1], conf->srcSize[0]*conf->srcSize[1]),
		lambdaS(lambdaS)
		{
	// initial s;
	nLens = param.parameter.size();


	occupation  = 0 ; 

	lambdaC = 1.0; 

	for(int i=0; i<nLens; ++i) {

		totalParam += param.nParam[i];
	}
	
	normVec n(0, 0, 0);
	vector<normVec> temp;
	temp.push_back(n);

	for(int i=0; i<conf->length; ++i) {
		//s.push_back(0);
		normV.push_back(temp);
		dSy1.push_back(0);
		dSy2.push_back(0);
	}

	for(int i=0; i<conf->srcSize[0]*conf->srcSize[1]; ++i) {
			H0.insert(i, i)= 1;
	}

	Hessian_grad.reserve(Eigen::VectorXi::Constant(conf->srcSize[0]*conf->srcSize[1], 6));  
	update_H_grad(conf); 
	 

}



vector<double> Model::getDeflectionAngle(Conf* conf, double pfX, double pfY, double *pDeltaX, double *pDeltaY, MultModelParam * param) {
	double fDenom = 0 ; 
	double srcX   = 0;  
	double srcY   = 0 ; 
	double fX  = 0; 
	double fY  = 0; 

	vector<double> srcPos;
	int nLens = param->parameter.size(); 
	


	//cout << "nLens: " <<  param->parameter.size() << endl; 
	for(int i=0; i<nLens; ++i) {
		fX =  pfX - (param->parameter[i].centerX * conf->imgRes);   // lens center frame, 
		fY =  pfY - (param->parameter[i].centerY * conf->imgRes);
		// Unit:  aresecond.

		if(param->parameter[i].name.compare("PTMASS")==0) {
			fDenom = fX*fX+fY*fY;
			double fMult = param->parameter[i].critRad*param->parameter[i].critRad/fDenom;
			*pDeltaX +=  fX*fMult;
			*pDeltaY +=  fY*fMult;
		}
		if(param->parameter[i].name.compare("SIE")==0) {
			double phi,root1mq,fq,fac,fCosTheta,fSinTheta,x1,y1,deltax1,deltay1;
			double fCore= param->parameter[i].core; 
			if (fX == 0 && fY == 0)  {
				*pDeltaX += param->parameter[i].critRad; 
				*pDeltaY += param->parameter[i].critRad;

			}
			//pre-calculate constants
			fCosTheta = cos(param->parameter[i].PA*M_PI/180 + 0.5* M_PI);
			fSinTheta = sin(param->parameter[i].PA*M_PI/180 + 0.5* M_PI);
			fq = 1-param->parameter[i].e;
			if (fq>1.0) cout << "Axis ratio should be smaller than 1. " << endl;
			if (fq==1.0) fq = 0.999;

			//rotate reference frame to x-axis
			x1 = fX*fCosTheta + fY*fSinTheta;
			y1 = -fX*fSinTheta + fY*fCosTheta;

			root1mq = sqrt(1.0-fq*fq);
			phi = sqrt(fq*fq*(fCore*fCore + x1*x1) + y1*y1);
			fac = param->parameter[i].critRad*sqrt(fq)/root1mq;
            if (phi==0. && fCore==0.) {
                deltax1 = 0.;
                deltay1=0.;
            } else {
                deltax1 = fac*atan(root1mq*x1/(phi + fCore));
                deltay1 = fac*lm_arctanh(root1mq*y1/(phi+ fCore*fq*fq));
            }

			//cout << root1mq << "\t" << y1 << "\t " << phi << "\t" << fq << "\t" << deltay1 << endl; 
			*pDeltaX += (deltax1*fCosTheta - deltay1*fSinTheta);
			*pDeltaY += (deltay1*fCosTheta + deltax1*fSinTheta);
		}
		
		
		if(param->parameter[i].name.compare("NFW")==0) {
		
			double fEllip,fCosTheta,fSinTheta,x1,y1,fPhi,fAngRadius,fTempResult,fCosPhi,fSinPhi,fScale;
  			fCosTheta = cos(param->parameter[i].PA*M_PI/180 + 0.5*M_PI);
			fSinTheta = sin(param->parameter[i].PA*M_PI/180 + 0.5*M_PI);

			fEllip = param->parameter[i].e; 
			fScale = param->parameter[i].radScale;

            if (fEllip >= 1.0 || fEllip < 0 || fScale < 0) 
                cout << "Bad parameters of 'e' or 'scale'. " << endl; 
                    
      		// create elliptical co-ords still in angle units from rotated frame sky coords 
			x1 = sqrt(1.0 - fEllip)*(fX*fCosTheta + fY*fSinTheta);
			y1 = sqrt(1.0 + fEllip)*(-fX*fSinTheta + fCosTheta*fY);
			fPhi = atan2(y1,x1);

			// angular radius is in dimensionless units 
			fAngRadius = sqrt(x1*x1 + y1*y1)/fScale;

			if (fAngRadius > 0.0) {
				double	deflx,defly;

				fCosPhi = cos(fPhi);
				fSinPhi = sin(fPhi);
					
				fTempResult = param->parameter[i].massScale * fScale * lm_nfw_mass(fAngRadius)/(fAngRadius);
				//cout << param.parameter[i].radScale << "\t" << fScale << "\t" << fAngRadius << "\t" << fTempResult  << endl; 
				deflx = sqrt(1.-fEllip)*fTempResult*fCosPhi;
				defly = sqrt(1.+fEllip)*fTempResult*fSinPhi;
				*pDeltaX += (deflx*fCosTheta - defly*fSinTheta);
				*pDeltaY += (deflx*fSinTheta + defly*fCosTheta);
			}
			else {
				*pDeltaX += 0.0;
				*pDeltaY += 0.0;
			} 
		}
		

		if(param->parameter[i].name.compare("SPEMD")==0)  {  

			double fTempKappa = 0.0; 
			double fTempCoreSqu=0.0; 
			double fTempAxratio =1.0; 
			double fTempDefl[2]; 
			double fTempGamma = 0.0; 
			double fTempCenterX=0; 
			double fTempCenterY=0;
			double	x1, y1;
			double 	fCosTheta, fSinTheta;

            
			fTempAxratio = 1.0 - param->parameter[i].e;
			fTempCoreSqu = param->parameter[i].core*param->parameter[i].core;
			fTempGamma = param->parameter[i].power;
            if (fTempAxratio > 1.0 || fTempAxratio<0 ) {
                cout << "Ellipticity is out of range; " << endl; 
            }
			if (fX == 0 && fY == 0 && fTempGamma >= 0.5) {
				*pDeltaX = param->parameter[i].critRad;
				*pDeltaY = param->parameter[i].critRad;
				//iStatus =LM_IGNORE_PROJ;
				//break;
			}
			fCosTheta = cos(param->parameter[i].PA*M_PI/180 + 0.5*M_PI);
			fSinTheta = sin(param->parameter[i].PA*M_PI/180 + 0.5*M_PI);

                
            x1 = (fX*fCosTheta + fY*fSinTheta);
            y1 = (-fX*fSinTheta + fY*fCosTheta);

            //cout << x1 << "\t" << y1 << endl; 
			fTempKappa = 0.5 * param->parameter[i].critRad * pow((2.0-2.0*fTempGamma)/fTempAxratio,fTempGamma);
			
			//cout << "ftempkapp: " << fTempKappa << "\t" << param->parameter[i].critRad  << endl; 
			fastelldefl_(&x1,&y1,&fTempKappa,&fTempGamma,&fTempAxratio,&fTempCoreSqu,fTempDefl);

			*pDeltaX = fTempDefl[0]*fCosTheta - fTempDefl[1]*fSinTheta;
			*pDeltaY = fTempDefl[1]*fCosTheta + fTempDefl[0]*fSinTheta;  
		}

		/* parameters: kappa,axis ratio (0-1),angle,scale len, M parameter */
		if(param->parameter[i].name.compare("SERSIC")==0)  {
		
				// double fCosTheta,fSinTheta,x1,y1,deflx,defly;
				// double fEllip = param->parameter[i].e;
    //             if (fEllip > 1.0 || fEllip < 0) {
    //                 cout << "Ellipticity is out of range; " << endl; 
    //                 break;
    //             }
				// /* rotate so that major axis of mass in x direction */
				// fCosTheta = cos(param->parameter[i].PA*M_PI/180 + 0.5*M_PI);
				// fSinTheta = sin(param->parameter[i].PA*M_PI/180 + 0.5*M_PI);
				

				// x1 = fX*fCosTheta + fY*fSinTheta;
				// y1 = -fX*fSinTheta + fY*fCosTheta;

				// iStatus = lm_deflCacheLookup(LM_SERSIC,&deVaucCache,x1,y1,pLensComp->fParameter[1],pLensComp->fParameter[3],pLensComp->fParameter[4], &deflx, &defly);
				// if (iStatus == LM_CACHE_MISS) {
				// 	iStatus = lm_CalcSersicDefl(x1,y1,pLensComp->fParameter[1],pLensComp->fParameter[3],pLensComp->fParameter[4], &deflx, &defly);
				// 	if (iStatus != 0) {
				// 		goto EXIT;
				// 	}
				// }

				// /* rotate back again */
				// *pDeltaX = (deflx*fCosTheta - defly*fSinTheta)*pLensComp->fParameter[0];
				// *pDeltaY = (defly*fCosTheta + deflx*fSinTheta)*pLensComp->fParameter[0];

		}


		/*

		


		if(param.parameter[i].name.compare("SERSIC")==0)  {  
			double fCosTheta,fSinTheta,x1,y1,deflx,defly;
			double fEllip = param.parameter[i].e; 
            if (fEllip > 1.0 || fEllip <= 0) 
            	  cout << "Bad parameters of 'e' . " << endl; 
			
			// rotate so that major axis of mass in x direction 
			fCosTheta = cos(param.parameter[i].PA*M_PI/180 + 0.5*M_PI);
			fSinTheta = sin(param.parameter[i].PA*M_PI/180 + 0.5*M_PI);

			x1 = fX*fCosTheta + fY*fSinTheta;
			y1 = -fX*fSinTheta + fY*fCosTheta;
			
			//iStatus = lm_deflCacheLookup(LM_SERSIC,&deVaucCache,x1,y1,pLensComp->fParameter[1],pLensComp->fParameter[3],pLensComp->fParameter[4], &deflx, &defly);
			
			if (iStatus == LM_CACHE_MISS) {
					iStatus = lm_CalcSersicDefl(x1,y1,pLensComp->fParameter[1],pLensComp->fParameter[3],pLensComp->fParameter[4], &deflx, &defly);
					if (iStatus != 0) {
						goto EXIT;
					}
				}
			
				
			*pDeltaX = (deflx*fCosTheta - defly*fSinTheta)*param.parameter[i].kap;
			*pDeltaY = (defly*fCosTheta + deflx*fSinTheta)*param.parameter[i].kap;

			

		}



		
		*/
		}

		
	// SrcX and SrcY are in unit aresec;   
	// (SrcX, SrcY) = (0, 0) is the center point of source plane; 
	srcX = pfX - (*pDeltaX);     
	srcY = pfY - (*pDeltaY);
	// In arcsecond
	srcPos.push_back(srcX);
	srcPos.push_back(srcY);
	return srcPos;

}




void Model::updatePosMapping(Image* image, Conf* conf) {

	length = conf->length;
	vector<double> srcPos;

	
	for(int i=0; i<length; ++i) {
		int imgX = image->xList[i];
		int imgY = image->yList[i];
		//cout << param.critRad << endl;
		double defX = 0 ;
		double defY = 0 ;

		double pfX = (imgX - conf->imgXCenter ) * conf->imgRes; 			// image center frame;
		double pfY = (imgY - conf->imgYCenter ) * conf->imgRes;

		srcPos = getDeflectionAngle(conf,pfX, pfY, &defX, &defY, &param);

		pDeltaX.push_back(defX);
		pDeltaY.push_back(defY);
		srcPosXList.push_back(srcPos[0]);
		srcPosYList.push_back(srcPos[1]);
		srcPosXListPixel.push_back(srcPos[0]/conf->srcRes+conf->srcXCenter);
		srcPosYListPixel.push_back(srcPos[1]/conf->srcRes+conf->srcYCenter);
		posMap[make_pair(imgX, imgY)] = i;
	}

}


/***************************
Function:   	update_H_zero
Description:  	(1) Member function of class Model;
				(2) H_zero is matrix convert an adaptive grid to a regular grid;
				(3) H_zero*S_adpative = S_regular ;
				(4) It will use the position of adaptive grid nodes:  srcPosXListPixel, srcPosYListPixel.
Arguments:		Conf*
Returns:		None
****************************/
void Model::update_H_zero(Conf* conf) {

	// T size:  [srcSize[0]*srcSize[1],  2*length]
	int x, y, iList;

	for(int i=0; i<length; ++i) {
		x = nearbyint(srcPosXListPixel[i]);
		y = nearbyint(srcPosYListPixel[i]);
		if(x>0 && x< conf->srcSize[0] && y>0 && y<conf->srcSize[1]) {
			iList = conf->srcSize[0]*y+x;
			H_zero.insert(iList, i) = 1.0;
		}
	}

}


void Model::update_H_grad(Conf* conf) {

	int src_length = conf->srcSize[0]*conf->srcSize[1]; 
	//sp_mat_row Hessian_grad(src_length, src_length); 

	for(int i=0; i<src_length; ++i) {
		int y = i/conf->srcSize[0]; 
		int x = i%conf->srcSize[0]; 
		vector<int> neighbors; 
		neighbors.push_back(i-1);  			// Left
		neighbors.push_back(i+1);  			// Right
		neighbors.push_back(i+conf->srcSize[0]);  // Up
		neighbors.push_back(i-conf->srcSize[0]);  // Down

		for(int n: neighbors) {
			if(n > 0 and n <src_length)  
				Hessian_grad.insert(i, n) = -1; 

		}
		Hessian_grad.insert(i, i) = 4; 
	}	
}




void Model::updateLensAndRegularMatrix(Image* dataImage,  Conf* conf) {

	map<pair<int, int>,int>::iterator left, right, up, down;
	vector<double> w;


	// Make pair; 
	
	// vector<pair<double, double>> coordPair(conf->length); 
	// for(int i=0; i<conf->length; ++i) 
	// 	coordPair[i] = make_pair(srcPosXListPixel[i], srcPosYListPixel[i]); 

	// sort(coordPair.begin(), coordPair.end(), [=](pair<double, double>& a,  pair<double, double> & b) {
	// 	return a.first < b.first; 
	// }); 

	for (int i=0; i<conf->length; ++i) {


		if(dataImage->type[i]==1) {
			L.insert(i,i)=1;
		}

        else if(dataImage->type[i]==0) {
            left  = posMap.find(make_pair(dataImage->xList[i]-1, dataImage->yList[i]));
            right = posMap.find(make_pair(dataImage->xList[i]+1, dataImage->yList[i]));
            up    = posMap.find(make_pair(dataImage->xList[i], dataImage->yList[i]+1));
            down  = posMap.find(make_pair(dataImage->xList[i], dataImage->yList[i]-1));

			if(left!=posMap.end() && up!=posMap.end() && down!=posMap.end()) {
				int iLeft = left->second;
				int iUp   = up  ->second;
				int iDown = down->second;

				Point A(srcPosXList[iLeft], srcPosYList[iLeft], s(iLeft));
				Point B(srcPosXList[iUp  ], srcPosYList[iUp	 ], s(iUp  ));
				Point C(srcPosXList[iDown], srcPosYList[iDown], s(iDown));
				Point P(srcPosXList[i	 ], srcPosYList[i	 ], s(i    ));

				w = getTriWeight( A, B, C, P);
				L.insert(i, iLeft) 	= w[0];
				L.insert(i, iUp	 )  = w[1];
				L.insert(i, iDown) 	= w[2];

				normVec n = getNormVector(A, B, C);
				normV[iLeft].push_back(n);
				normV[iUp  ].push_back(n);
				normV[iDown].push_back(n);
			}
			else L.insert(i, i) = 1;
		}
	}
}

void Model::updateGradient(Image* dataImage) {
	vec vSy1(length);  vSy1.fill(0);
	vec vSy2(length);  vSy2.fill(0);

	for (int i=0; i<length; ++i) {
		if(dataImage->type[i]==1) {
			normVec mean =  meanNormVector(normV[i]);
			if (mean.n2==0) {
				vSy1[i] = 0;
				vSy2[i] = 0;
			}
			else {
				vSy1[i] = -mean.n0/mean.n2;
				vSy2[i] = -mean.n1/mean.n2;
			}
		}
	}

	vSy1 = L*vSy1;
	vSy2 = L*vSy2;
	//update Ds and Dphi

	for(int i=0; i<length; ++i){
		Ds.coeffRef(i, 2*i+0) = vSy1(i);
		Ds.coeffRef(i, 2*i+1) = vSy2(i);
	}

}


void Model::updateVegettiRegularization() {
	
	vector<double> w5;
	// vector<TriData> TriDataListX ; 

	// for(int i=0; i<length; ++i) {
	// 	TriDataListX.push_back(TriData(srcPosXListPixel[i], srcPosYListPixel[i], i)); 
	// }
	// vector<TriData> TriDataListY(TriDataListX); 
	// sort(TriDataListX.begin(), TriDataListX.end(), [](const TriData & left, const TriData & right){
	// 		return left.srcX < right.srcX; }); 
	// sort(TriDataListY.begin(), TriDataListY.end(), [](const TriData & left, const TriData & right){
	// 		return left.srcY < right.srcY; }); 

	// Build cloud; 
	
	PointCloud<double> cloud;
	cloud.pts.resize(length);
    KDTreeAdaptor index(2, cloud, KDTreeSingleIndexAdaptorParams(20 ));
	for(size_t j=0; j<length; ++j) {
		cloud.pts[j].x = srcPosXListPixel[j];
		cloud.pts[j].y = srcPosYListPixel[j];
	}
	index.buildIndex();


	int indexRegion1 = -1;
	int indexRegion2 = -1;
	int indexRegion3 = -1;
	int indexRegion4 = -1;


	if(1){
	for (int i=0; i<length; ++i) {
		
        double query_pt[2] = {srcPosXListPixel[i], srcPosYListPixel[i]};
		const size_t num_results =10;
		std::vector<size_t>   ret_index(num_results);
		std::vector<double> out_dist_sqr(num_results);
		index.knnSearch(&query_pt[0], num_results, &ret_index[0], &out_dist_sqr[0]);
        double delta(1e-6);

		indexRegion1 = -1;
		indexRegion2 = -1;
		indexRegion3 = -1;
		indexRegion4 = -1;

		for(int k=1; k<num_results; ++k) {

			size_t j = ret_index[k] ;
			//cout << i<< "\t" << j << "\t" << srcPosXListPixel[j] << "\t" << srcPosXListPixel[i] << "\t" << srcPosYListPixel[j] << "\t" << srcPosYListPixel[i] << endl; 
			if( indexRegion1 == -1 and srcPosXListPixel[j] > srcPosXListPixel[i] + delta and srcPosYListPixel[j] > srcPosYListPixel[i]    ) {   // in region1; 
					indexRegion1 = j ;
			}
			else if(indexRegion2 == -1 and srcPosXListPixel[j] < srcPosXListPixel[i] - delta and srcPosYListPixel[j] > srcPosYListPixel[i]  ) {
					indexRegion2 = j ;
			}
			else if(indexRegion3 == -1 and srcPosXListPixel[j] < srcPosXListPixel[i] - delta and srcPosYListPixel[j] < srcPosYListPixel[i] ) {   // in region3; 
					indexRegion3 = j ;
			}
			else if(indexRegion4 == -1 and srcPosXListPixel[j] > srcPosXListPixel[i] + delta and srcPosYListPixel[j] < srcPosYListPixel[i]  ) {   // in region4; 
					indexRegion4 = j ;
			}
			else if(indexRegion1 != -1 and indexRegion2 != -1 and indexRegion3 != -1 and indexRegion4 != -1)
				break;
		}

		if (indexRegion1 > 0 and indexRegion2 > 0  and indexRegion3 > 0 and indexRegion4 > 0) {
			Point A(srcPosXList[indexRegion3 ], srcPosYList[indexRegion3 ], s(indexRegion3 ));
			Point B(srcPosXList[indexRegion2 ], srcPosYList[indexRegion2 ], s(indexRegion2 ));
			Point C(srcPosXList[i            ], srcPosYList[i            ], s(i            ));
			Point D(srcPosXList[indexRegion4 ], srcPosYList[indexRegion4 ], s(indexRegion4 ));
			Point E(srcPosXList[indexRegion1 ], srcPosYList[indexRegion1 ], s(indexRegion1 ));

			w5 = getPentWeigth(A, B, C, D, E);
			Hs1.insert(i, indexRegion3	) 	= w5[0];
			Hs1.insert(i, indexRegion2	) 	= w5[1];
			Hs1.insert(i, i				) 	= w5[2];
			Hs1.insert(i, indexRegion4	) 	= w5[3];
			Hs1.insert(i, indexRegion1	) 	= w5[4];

			Hs2.insert(i, indexRegion3	) 	= w5[5];
			Hs2.insert(i, indexRegion2	) 	= w5[6];
			Hs2.insert(i, i				) 	= w5[7];
			Hs2.insert(i, indexRegion4	) 	= w5[8];
			Hs2.insert(i, indexRegion1	) 	= w5[9];
		}
	}
	}
	

	if(0) {  // naive vegetti way to solve the problem; 

	for (int i=0; i<length; ++i) {
		
		double largeNumber = 10000; 
		double distRegion1 = largeNumber; 
		double distRegion2 = largeNumber; 
		double distRegion3 = largeNumber; 
		double distRegion4 = largeNumber; 

		indexRegion1 = -1; 
		indexRegion2 = -1; 
		indexRegion3 = -1; 
		indexRegion4 = -1; 
		for(int j=0; j<length; ++j) {
			if (j == i ) 	continue; 
			double dx = srcPosXListPixel[j] - srcPosXListPixel[i] ; 
			double dy = srcPosYListPixel[j] - srcPosYListPixel[i] ; 
			double dist = dx * dx + dy * dy; 

			if(srcPosXListPixel[j] > srcPosXListPixel[i] and srcPosYListPixel[j] > srcPosYListPixel[i] and dist<distRegion1 ) {   // in region1; 
					distRegion1 = dist ; 
					indexRegion1 = j ;  
			}
			else if(srcPosXListPixel[j] < srcPosXListPixel[i] and srcPosYListPixel[j] > srcPosYListPixel[i]  and dist<distRegion2 ) {   // in region2; 
					distRegion2 = dist ; 
					indexRegion2 = j ;  
			}
			else if(srcPosXListPixel[j] < srcPosXListPixel[i] and srcPosYListPixel[j] < srcPosYListPixel[i] and dist<distRegion3 ) {   // in region3; 
					distRegion3 = dist ; 
					indexRegion3 = j ;  
			}
			else if(srcPosXListPixel[j] > srcPosXListPixel[i] and srcPosYListPixel[j] < srcPosYListPixel[i]  and dist<distRegion4 ) {   // in region4; 
					distRegion4 = dist ; 
					indexRegion4 = j ;  
			}
		}


			if (indexRegion1 > 0 and indexRegion2 > 0  and indexRegion3 > 0 and indexRegion4 > 0
			and indexRegion1 != i and indexRegion2 != i  and indexRegion3 != i and indexRegion4 != i  )  {
			Point A(srcPosXList[indexRegion3 ], srcPosYList[indexRegion3 ], s(indexRegion3 ));
			Point B(srcPosXList[indexRegion2 ], srcPosYList[indexRegion2 ], s(indexRegion2 ));
			Point C(srcPosXList[i            ], srcPosYList[i            ], s(i            ));
			Point D(srcPosXList[indexRegion4 ], srcPosYList[indexRegion4 ], s(indexRegion4 ));
			Point E(srcPosXList[indexRegion1 ], srcPosYList[indexRegion1 ], s(indexRegion1 ));

			w5 = getPentWeigth(A, B, C, D, E);
			Hs1.insert(i, indexRegion3	) 	= w5[0];
			Hs1.insert(i, indexRegion2	) 	= w5[1];
			Hs1.insert(i, i				) 	= w5[2];
			Hs1.insert(i, indexRegion4	) 	= w5[3];
			Hs1.insert(i, indexRegion1	) 	= w5[4];

			Hs2.insert(i, indexRegion3	) 	= w5[5];
			Hs2.insert(i, indexRegion2	) 	= w5[6];
			Hs2.insert(i, i				) 	= w5[7];
			Hs2.insert(i, indexRegion4	) 	= w5[8];
			Hs2.insert(i, indexRegion1	) 	= w5[9]; 
		}
	}
}

	HtH = Hs1.transpose()*Hs1 + Hs2.transpose()*Hs2;

}



void Model::updateCompactMatrix(Image* dataImage) {
	// update CompactMatrix HcH; 
	
	//Suppose we know the 'srcPosXListPixel' and 'srcPosYListPixel': 


	//vector<double> w(srcPosXListPixel.size(), 1); 
	vector<double> w = dataImage->dataList; 
	Kmeans kmeans(k, srcPosXListPixel, srcPosYListPixel, Kmeans_group, w); 
	kmeans.initiate(); 
	
	// using the weight to partition the source images; 

	double diff = INT_MAX;
	vector<pair<double, double>> pre_centers ; 
	int iterNum = 0; 
	while( diff > 1.0e-6  && iterNum < 20 ) {  // unit pixel; 
		iterNum ++; 
		kmeans.updateGroup(); 
		vector<int> count = kmeans.countGroup(); 
		pre_centers = kmeans.centers; 
		kmeans.updateCenters(); 
		diff = 0; 
		for(int i=0; i<k; ++i) {
			double xdiff = pre_centers[i].first - kmeans.centers[i].first; 
			double ydiff = pre_centers[i].second -kmeans.centers[i].second; 
			diff += (xdiff * xdiff + ydiff * ydiff );
		}
	}

	// Now the centers are found; 
	vector<int> countGroup  = kmeans.countGroup();   // countGroup.size() == k;

	// save group and cetners information; 
	Kmeans_group = kmeans.group; 
	Kmeans_centers = kmeans.centers;

	//getKmeansScatter(dataImage);   // Kmeans_group, Kmeans_center; 
	for(int i=0; i<length; ++i) {
		double xdiff = srcPosXListPixel[i] - Kmeans_centers[Kmeans_group[i]].first; 
		double ydiff = srcPosYListPixel[i] - Kmeans_centers[Kmeans_group[i]].second; 
		HcH.insert(i, i) = xdiff * xdiff + ydiff * ydiff ; 
	}
}

void Model::Logging(Image* dataImage, Conf* conList, string outFileName) {
	ofstream f(outFileName);
	string tab = "\t";
	string entry;
	f << "#1 index\n" << "#2 imgX\n" << "#3 imgY\n" << "#4 imgBri\n";
	f << "#5 srcX\n"  << "#6 srcY\n" << "#7 mapIndex" << "#8 ";
	for(int i=0; i<conList->length; ++i) {
		entry =  to_string(dataImage->iList[i]) + "\t"
				+to_string(dataImage->xList[i]) + "\t"
				+to_string(dataImage->yList[i]) + "\t"
				+to_string(dataImage->dataList[i]) + "\t"
				+to_string(srcPosXList[i]) + "\t"
				+to_string(srcPosYList[i]) + "\t"
				+to_string(posMap[make_pair(dataImage->xList[i], dataImage->yList[i])]) + "\t" ;  //  Matched position
				//+to_string()
		f << entry << endl ;
	}
	f.close();
}



void Model::solveSource(sp_mat* invC, vec* d , string R_type) {

	//cout << "Beta: " <<  beta << endl; 
	if (R_type == "zero") 
		REG = H_zero.transpose()*H_zero; 
	else if(R_type == "grad")
		REG = H_zero.transpose() * Hessian_grad * H_zero; 

	else if(R_type == "vege")  {
		updateVegettiRegularization(); 
		REG = HtH + beta * HcH ;
	}
	else {
		cout << "Regularization type is not supported yet!" << endl; 
		exit(1); 
	} 

	
	sp_mat 	A = lambdaC * lambdaC * L.transpose()*(*invC)*L + lambdaS * lambdaS  * REG; // Hessian_grad.transpose();
	vec 	b = lambdaC * lambdaC * L.transpose()*(*invC)*(*d);
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > chol(A);
	s = chol.solve(b);
}


void Model::updateSource(Conf* conf) {

	std::default_random_engine generator; 
	std::normal_distribution<double> distribution(conf->back_mean, conf->back_std); 
	for(int i=0; i<conf->length; ++i) {
		double randNum = distribution(generator); 
		s[i]  += randNum; 
	}

}




Image* Model::getFullResidual(Image* dataImage) {
	// Assume "mod_image" is known;  's' is the source brightness after solve the linear equaion; 
	mod_img = eigenV_to_cV(&s); 

	Image* fullResidualImage = new Image(* dataImage); 
	//map<pair<int, int>, double> modMap; 
	for(int i=0; i< mod_img.size(); ++i)  {
		int pos = dataImage->yList[i] * dataImage->naxis1 + dataImage->xList[i]; 
		fullResidualImage->data[pos] -= mod_img[i]; 
	}
	return fullResidualImage; 	
}

void Model::writeSrcImage(string outFileName, Conf* conList) {
	vector<double> sBright = eigenV_to_cV(&s);
	Image* srcImg = new Image(srcPosXListPixel, srcPosYListPixel, &sBright, conList->srcSize[0], conList->srcSize[1], conList->bitpix);
	srcImg->writeToFile(outFileName);
	delete srcImg;

}


double Model::getRegularizationSrcValue (vec d) {
	// HtH and d is known
	return  d.transpose()*HtH*d;
}


double Model::getScatterReg() {

	// Suppose we know the 'srcPosXListPixel' and 'srcPosYListPixel': 

	double scatter = 0;   // STD of position 'x' and 'y'
	double sumX = 0; 
	double sumY = 0;  
	for (int i=0; i<srcPosXList.size(); ++i) {
		//cout << srcPosXList[i] << "\t" ; 
		sumX += srcPosXList[i]; 
		sumY += srcPosYList[i]; 
	}

	

	double xPosMean = sumX / srcPosXList.size(); 
	double yPosMean = sumY / srcPosYList.size();

	for (int i=0; i<srcPosXList.size(); ++i) {
		scatter += (srcPosXList[i]-xPosMean) * (srcPosXList[i]-xPosMean) ;
		scatter += (srcPosYList[i]-yPosMean) * (srcPosYList[i]-yPosMean) ;
	}
	scatter = scatter / srcPosXList.size();  
	return scatter; 
	
}


double Model::getKmeansScatter(Image* dataImage) {
	//Suppose we know the 'srcPosXListPixel' and 'srcPosYListPixel': 


	//vector<double> w(srcPosXListPixel.size(), 1); 
	vector<double> w = dataImage->dataList; 
	Kmeans kmeans(k, srcPosXListPixel, srcPosYListPixel, Kmeans_group, w); 
	kmeans.initiate(); 
	
	// using the weight to partition the source images; 

	double diff = INT_MAX;
	vector<pair<double, double>> pre_centers ; 
	int iterNum = 0; 
	while( diff > 1.0e-6  && iterNum < 20 ) {
		iterNum ++; 
		kmeans.updateGroup(); 
		vector<int> count = kmeans.countGroup(); 
		pre_centers = kmeans.centers; 
		kmeans.updateCenters(); 
		diff = 0; 
		for(int i=0; i<k; ++i) {
			double xdiff = pre_centers[i].first - kmeans.centers[i].first; 
			double ydiff = pre_centers[i].second -kmeans.centers[i].second; 
			diff += (xdiff * xdiff + ydiff * ydiff );  
		}
	}

	// Now the centers are found; 
	vector<int> countGroup  = kmeans.countGroup();   // countGroup.size() == k; 

	// save group and cetners information; 
	Kmeans_group = kmeans.group; 
	Kmeans_centers = kmeans.centers; 



	double scatter  = 0; 
	for(int i=0; i<srcPosXListPixel.size(); ++i) {
		double xdiff = srcPosXListPixel[i] - kmeans.centers[kmeans.group[i]].first; 
		double ydiff = srcPosYListPixel[i] - kmeans.centers[kmeans.group[i]].second; 
		scatter += kmeans.w[i] * (xdiff * xdiff + ydiff * ydiff) ; 
	}
	return scatter/srcPosXListPixel.size(); 

}


double Model::getZerothOrderReg (Conf* conf, vector<double> briList) {
	//s is known;
	double sum = 0;
	long naxis1 = conf->srcSize[0];
	long naxis2 = conf->srcSize[1];
	Image* srcImg =new Image(srcPosXListPixel, srcPosYListPixel, &briList, naxis1, naxis2, conf->bitpix);

	for (int i=0; i< naxis1 ; ++i) {
		for (int j=0; j< naxis2 ; ++j) {
			int index = i+j*naxis1;

			double val = srcImg->data[index] * srcImg->data[index]; 
			if(val!=0) {
				sum += val;
				occupation += 1; 
			}

		}
	}
	delete srcImg;
	return sum ;
}


double Model::getGradientOrderReg(Conf* conf, vector<double> briList) {
	double sum = 0;
	double diff = 0;  
	long naxis1 = conf->srcSize[0];
	long naxis2 = conf->srcSize[1];
	int index = 0 ; 
	int index_edge = 0; 
	int index_next = 0; 
	Image* srcImg =new Image(srcPosXListPixel, srcPosYListPixel, &briList, naxis1, naxis2, conf->bitpix);

	for (int i=0; i< naxis1 ; ++i) {
		index_edge = (naxis2 - 1)* naxis1 + i; 
		sum +=  srcImg->data[index_edge] * srcImg->data[index_edge]; 		
		for (int j=0; j< naxis2-1 ; ++j) {
			index  		= i + j 	* naxis1;
			index_next 	= i + (j+1) * naxis1; 
			diff = srcImg->data[index] - srcImg->data[index_next]; 
			sum += diff*diff; 
		}
	}

	for (int j=0; j< naxis2 ; ++j) {
		index_edge = j * naxis2 + (naxis1-1) ; 
		sum +=  srcImg->data[index_edge] *  srcImg->data[index_edge]; 
		for (int i=0; i< naxis1-1 ; ++i) {
			index  		= i 	+ j * naxis1;
			index_next 	= (i+1) + j * naxis1; 
			diff = srcImg->data[index] - srcImg->data[index_next]; 
			sum += diff*diff; 
		}
	}
	delete srcImg; 
	return sum ; 
}


double Model::getCurvatureOrderReg(Conf* conf, vector<double> briList) {

	double sum = 0;
	double diff = 0;  
	double diff_edge = 0; 
	long naxis1 = conf->srcSize[0];
	long naxis2 = conf->srcSize[1];

	int index1 = 0 ;
	int index2 = 0 ;
	int index3 = 0 ;

	int index_edge1 = 0; 
	int index_edge2 = 0;	 

	Image* srcImg =new Image(srcPosXListPixel, srcPosYListPixel, &briList, naxis1, naxis2, conf->bitpix);

	for (int i=0; i< naxis1 ; ++i) {
		index_edge1 = (naxis2 - 2)* naxis1 + i;
		index_edge2 = (naxis2 - 1)* naxis1 + i; 
		diff_edge = srcImg->data[index_edge1] -  srcImg->data[index_edge2]; 
		sum +=  diff_edge * diff_edge +  srcImg->data[index_edge2]*  srcImg->data[index_edge2]; 

		for ( int j=0 ; j< naxis2-2; ++j) {
			index1 = i + j  * naxis1; 
			index2 = index1 + naxis1; 
			index3 = index2 + naxis1; 
			diff = srcImg->data[index1] - 2* srcImg->data[index2] + srcImg->data[index3]; 
			sum += diff * diff; 

		}

	}	
	
	for ( int j=0 ; j< naxis2; ++j) {
		index_edge1 = j * naxis1 + naxis1-2;
		index_edge2 = j * naxis1 + naxis1-1; 
		diff_edge = srcImg->data[index_edge1] -  srcImg->data[index_edge2]; 
		sum +=  diff_edge * diff_edge +  srcImg->data[index_edge2]*  srcImg->data[index_edge2];

		for (int i=0; i< naxis1-2 ; ++i) {
			index1 = i + j  * naxis1; 
			index2 = index1 + 1; 
			index3 = index2 + 1; 
			diff = srcImg->data[index1] - 2* srcImg->data[index2] + srcImg->data[index3]; 
			sum += diff * diff; 
		}
	}	
	delete srcImg; 
	return sum ; 
}



MultModelParam::MultModelParam(map<string,string> confMap) {
		

		nLens = 0;
		map<string, string>::iterator itPTMASS	= confMap.find("PTMASS") ;
		map<string, string>::iterator itSIE    	= confMap.find("SIE");
		map<string, string>::iterator itNFW   	= confMap.find("NFW");
		map<string, string>::iterator itSPEMD 	= confMap.find("SPEMD");
		map<string, string>::iterator itSERSIC 	= confMap.find("SERSIC");

		if(itSERSIC != confMap.end()) {
			vector<string> strs;
			std::string s = itSIE->second; 
			string delimiter = "_&&_";

			size_t pos = s.find(delimiter) ; 
			while( pos!=std::string::npos) {
				strs.push_back(s.substr(0, pos)); 
				s = s.substr(pos+4); 
				pos = s.find(delimiter);
			}; 
			strs.push_back(s); 
			for(int i=0; i<strs.size(); ++i) {

				

				vector<string> items = splitString(itSERSIC->second);
				SingleModelParam tempParam;
				tempParam.name = "SERSIC";
				tempParam.centerXFrom 		= stof(items[0]);
				tempParam.centerXTo 		= stof(items[1]);
				tempParam.centerXInc 		= stof(items[2]);
				tempParam.centerYFrom 		= stof(items[3]);
				tempParam.centerYTo 		= stof(items[4]);
				tempParam.centerYInc 		= stof(items[5]);
				tempParam.kapFrom 			= stof(items[6]);
				tempParam.kapTo   			= stof(items[7]);
				tempParam.kapInc  			= stof(items[8]);
				tempParam.eFrom 			= stof(items[9]);
				tempParam.eTo   			= stof(items[10]);
				tempParam.eInc  			= stof(items[11]);
				tempParam.PAFrom			= stof(items[12]);
				tempParam.PATo 				= stof(items[13]);
				tempParam.PAInc 			= stof(items[14]);
				tempParam.sersicScaleFrom 	= stof(items[15]);
				tempParam.sersicScaleTo		= stof(items[16]);
				tempParam.sersicScaleInc	= stof(items[17]);
				tempParam.mFrom 			= stof(items[18]);
				tempParam.mTo   			= stof(items[19]);
				tempParam.mInc  			= stof(items[20]);
			
				assert (tempParam.centerXInc 	 >DIFF and 
						tempParam.centerYInc 	 >DIFF and 
						tempParam.kapInc 		 >DIFF and 
						tempParam.eInc       	 >DIFF and  
						tempParam.PAInc		 	 >DIFF and 
						tempParam.sersicScaleInc >DIFF and 
						tempParam.mInc 		 	 >DIFF); 

				parameter.push_back(tempParam);
				nParam.push_back(NUM_SERSIC_PARAM);

				nLens +=1;
			}
			
		}

		if(itPTMASS != confMap.end()) {
			vector<string> strs;
			std::string s = itPTMASS->second; 
			string delimiter = "_&&_";
			size_t pos = s.find(delimiter) ; 
			while( pos!=std::string::npos) {
				strs.push_back(s.substr(0, pos)); 
				s = s.substr(pos+4); 
				pos = s.find(delimiter);
				break; 
			}; 
			strs.push_back(s); 
			for(int i=0; i<strs.size(); ++i) {
				vector<string> items = splitString(strs[i]);
				SingleModelParam tempParam;
				tempParam.name = "PTMASS";
				tempParam.centerXFrom   = stof(items[0]);
				tempParam.centerXTo 	= stof(items[1]);
				tempParam.centerXInc 	= stof(items[2]);
				tempParam.centerYFrom 	= stof(items[3]);
				tempParam.centerYTo 	= stof(items[4]);
				tempParam.centerYInc 	= stof(items[5]);
				tempParam.critRadFrom	= stof(items[6]);
				tempParam.critRadTo   	= stof(items[7]);
				tempParam.critRadInc  	= stof(items[8]);

				parameter.push_back(tempParam);
				nParam.push_back(NUM_PTMASS_PARAM);
				nLens +=1;
			}

		}
		if(itSIE != confMap.end()) {
			
			vector<string> strs;
			std::string s = itSIE->second; 
			string delimiter = "_&&_";

			size_t pos = s.find(delimiter) ; 
			while( pos!=std::string::npos) {
				strs.push_back(s.substr(0, pos)); 
				s = s.substr(pos+4); 
				pos = s.find(delimiter);
			}; 
			strs.push_back(s); 
			for(int i=0; i<strs.size(); ++i) {

				vector<string> items = splitString(strs[i]);
				SingleModelParam tempParam;

				tempParam.name = "SIE";

				tempParam.centerXFrom 	= stof(items[0]);
				tempParam.centerXTo 	= stof(items[1]);
				tempParam.centerXInc 	= stof(items[2]);
				tempParam.centerYFrom 	= stof(items[3]);
				tempParam.centerYTo 	= stof(items[4]);
				tempParam.centerYInc 	= stof(items[5]);
				tempParam.critRadFrom 	= stof(items[6]);
				tempParam.critRadTo   	= stof(items[7]);
				tempParam.critRadInc  	= stof(items[8]);
				tempParam.eFrom			= stof(items[9]);
				tempParam.eTo 			= stof(items[10]);
				tempParam.eInc 			= stof(items[11]);
				tempParam.PAFrom 		= stof(items[12]);
				tempParam.PATo			= stof(items[13]);
				tempParam.PAInc			= stof(items[14]);
				tempParam.coreFrom		= stof(items[15]); 
				tempParam.coreTo		= stof(items[16]); 
				tempParam.coreInc		= stof(items[17]); 
				assert (tempParam.centerXInc >DIFF and 
						tempParam.centerYInc >DIFF and 
						tempParam.critRadInc >DIFF and 
						tempParam.eInc       >DIFF and  
						tempParam.PAInc		 >DIFF and 
						tempParam.coreInc 	 >DIFF); 


				parameter.push_back(tempParam);
				nParam.push_back(NUM_SIE_PARAM);

				nLens +=1;
			}
		}

		if(itNFW != confMap.end()) {

			vector<string> strs;
			std::string s = itSIE->second; 
			string delimiter = "_&&_";

			size_t pos = s.find(delimiter) ; 
			while( pos!=std::string::npos) {
				strs.push_back(s.substr(0, pos)); 
				s = s.substr(pos+4); 
				pos = s.find(delimiter);
			}; 
			strs.push_back(s); 
			for(int i=0; i<strs.size(); ++i) {

				vector<string> items = splitString(itNFW->second);
				SingleModelParam tempParam;

				tempParam.name = "NFW";

				tempParam.centerXFrom 	= stof(items[0]);
				tempParam.centerXTo 	= stof(items[1]);
				tempParam.centerXInc 	= stof(items[2]);
				tempParam.centerYFrom 	= stof(items[3]);
				tempParam.centerYTo 	= stof(items[4]);
				tempParam.centerYInc 	= stof(items[5]);
				tempParam.massScaleFrom = stof(items[6]);
				tempParam.massScaleTo   = stof(items[7]);
				tempParam.massScaleInc  = stof(items[8]);
				tempParam.radScaleFrom 	= stof(items[9]);
				tempParam.radScaleTo   	= stof(items[10]);
				tempParam.radScaleInc  	= stof(items[11]);
				tempParam.eFrom			= stof(items[12]);
				tempParam.eTo 			= stof(items[13]);
				tempParam.eInc 			= stof(items[14]);
				tempParam.PAFrom 		= stof(items[15]);
				tempParam.PATo			= stof(items[16]);
				tempParam.PAInc			= stof(items[17]);

				assert (tempParam.centerXInc   >DIFF and 
						tempParam.centerYInc   >DIFF and 
						tempParam.massScaleInc >DIFF and 
						tempParam.radScaleInc  >DIFF and 
						tempParam.eInc         >DIFF and  
						tempParam.PAInc		   >DIFF); 

				parameter.push_back(tempParam);
				nParam.push_back(NUM_NFW_PARAM);

				nLens +=1;
			}
		}

		if(itSPEMD != confMap.end()) {

			vector<string> strs;
			std::string s = itSIE->second; 
			string delimiter = "_&&_";

			size_t pos = s.find(delimiter) ; 
			while( pos!=std::string::npos) {
				strs.push_back(s.substr(0, pos)); 
				s = s.substr(pos+4); 
				pos = s.find(delimiter);
			}; 
			strs.push_back(s); 
			for(int i=0; i<strs.size(); ++i) {
				vector<string> items = splitString(itSPEMD->second);
				SingleModelParam tempParam;
				tempParam.name = "SPEMD";
				tempParam.centerXFrom 	= stof(items[0]);
				tempParam.centerXTo 	= stof(items[1]);
				tempParam.centerXInc 	= stof(items[2]);
				tempParam.centerYFrom 	= stof(items[3]);
				tempParam.centerYTo 	= stof(items[4]);
				tempParam.centerYInc 	= stof(items[5]);
				tempParam.critRadFrom 	= stof(items[6]);
				tempParam.critRadTo   	= stof(items[7]);
				tempParam.critRadInc  	= stof(items[8]);
				tempParam.eFrom			= stof(items[9]);
				tempParam.eTo 			= stof(items[10]);
				tempParam.eInc 			= stof(items[11]);
				tempParam.PAFrom 		= stof(items[12]);
				tempParam.PATo			= stof(items[13]);
				tempParam.PAInc			= stof(items[14]);
				tempParam.powerFrom		= stof(items[15]);
				tempParam.powerTo 		= stof(items[16]);
				tempParam.powerInc 		= stof(items[17]);
				tempParam.coreFrom 		= stof(items[18]);
				tempParam.coreTo		= stof(items[19]);
				tempParam.coreInc		= stof(items[20]);

				parameter.push_back(tempParam);
				nParam.push_back(NUM_SPEMD_PARAM);
				nLens +=1;
			}
		}



		// compute nComb; 
		nComb = 0;
        /*
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
	        					nComb ++; 


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
        } */
        nComb = 1;
		for(size_t i=0; i<parameter.size(); ++i) {
			nComb *= (int((parameter[i].centerXTo-parameter[i].centerXFrom )/ parameter[i].centerXInc) +1);
			nComb *= (int((parameter[i].centerYTo-parameter[i].centerYFrom )/ parameter[i].centerYInc) +1);
			nComb *= (int((parameter[i].critRadTo-parameter[i].critRadFrom )/ parameter[i].critRadInc) +1);
			nComb *= (int((parameter[i].eTo-parameter[i].eFrom)/ parameter[i].eInc) +1);
			nComb *= (int((parameter[i].PATo-parameter[i].PAFrom )/ parameter[i].PAInc) +1);
			nComb *= (int((parameter[i].coreTo-parameter[i].coreFrom )/ parameter[i].coreInc) +1);

		}
        cout << "nComb "<< nComb << endl;


	}


void MultModelParam::printModels() {
	cout << "\n*********** Models *********" << endl;
	cout << "nLens: 	 " << nLens << endl; 
	for(int i=0; i<nLens; ++i) {
		if (parameter[i].name=="PTMASS") {

			cout << "[" << parameter[i].name  << "]:" 
				<< "\n_From\t"
				<< parameter[i].centerXFrom << "\t" 
				<< parameter[i].centerYFrom << "\t"
				<< parameter[i].critRadFrom << "\t"
				<< "\n_To\t"
				<< parameter[i].centerXTo<< "\t" 
				<< parameter[i].centerYTo << "\t"
				<< parameter[i].critRadTo << "\t"
				<< "\n_Inc\t"
				<< parameter[i].centerXInc << "\t" 
				<< parameter[i].centerYInc << "\t"
				<< parameter[i].critRadInc << "\t"
				<< endl; 
				
		}

		if (parameter[i].name=="SIE") {

			cout << "[" << parameter[i].name  << "]:" 
				<< "\n_From\t"
				<< parameter[i].centerXFrom << "\t" 
				<< parameter[i].centerYFrom << "\t"
				<< parameter[i].critRadFrom << "\t"
				<< parameter[i].eFrom << "\t" 
				<< parameter[i].PAFrom << "\t"
				<< parameter[i].coreFrom << "\t"
				<< "\n_To\t"
				<< parameter[i].centerXTo<< "\t" 
				<< parameter[i].centerYTo << "\t"
				<< parameter[i].critRadTo << "\t"
				<< parameter[i].eTo << "\t" 
				<< parameter[i].PATo << "\t" 
				<< parameter[i].coreTo << "\t"
				<< "\n_Inc\t"
				<< parameter[i].centerXInc << "\t" 
				<< parameter[i].centerYInc << "\t"
				<< parameter[i].critRadInc << "\t"
				<< parameter[i].eInc << "\t" 
				<< parameter[i].PAInc << "\t"
				<< parameter[i].coreInc << "\t"
				<< endl; 
				
		}

	

		if (parameter[i].name=="SERSIC") {

			cout << "[" << parameter[i].name  << "]:" 
				<< "\n_From\t"
				<< parameter[i].centerXFrom << "\t" 
				<< parameter[i].centerYFrom << "\t"
				<< parameter[i].kapFrom << "\t"
				<< parameter[i].eFrom << "\t" 
				<< parameter[i].PAFrom << "\t"
				<< parameter[i].sersicScaleFrom << "\t"
				<< parameter[i].mFrom << "\t"
				<< "\n_To\t"
				<< parameter[i].centerXTo << "\t" 
				<< parameter[i].centerYTo << "\t"
				<< parameter[i].kapTo << "\t"
				<< parameter[i].eTo << "\t" 
				<< parameter[i].PATo << "\t"
				<< parameter[i].sersicScaleTo << "\t"
				<< parameter[i].mTo << "\t"
				<< "\n_Inc\t"
				<< parameter[i].centerXInc << "\t" 
				<< parameter[i].centerYInc << "\t"
				<< parameter[i].kapInc << "\t"
				<< parameter[i].eInc << "\t" 
				<< parameter[i].PAInc << "\t"
				<< parameter[i].sersicScaleInc << "\t"
				<< parameter[i].mInc << "\t"
				<< endl; 
				
		}

		if (parameter[i].name=="NFW") {

			cout << "[" << parameter[i].name  << "]:" 
				<< "\n_From\t"
				<< parameter[i].centerXFrom << "\t" 
				<< parameter[i].centerYFrom << "\t"
				<< parameter[i].massScaleFrom << "\t"
				<< parameter[i].radScaleFrom << "\t"
				<< parameter[i].eFrom << "\t" 
				<< parameter[i].PAFrom << "\t" 
				<< "\n_To\t"
				<< parameter[i].centerXTo<< "\t" 
				<< parameter[i].centerYTo << "\t"
				<< parameter[i].massScaleTo << "\t"
				<< parameter[i].radScaleTo << "\t"
				<< parameter[i].eTo << "\t" 
				<< parameter[i].PATo << "\t" 
				<< "\n_Inc\t"
				<< parameter[i].centerXInc << "\t" 
				<< parameter[i].centerYInc << "\t"
				<< parameter[i].massScaleInc << "\t"
				<< parameter[i].radScaleInc << "\t"
				<< parameter[i].eInc << "\t" 
				<< parameter[i].PAInc << "\t" 
				<< endl; 
		}

		if (parameter[i].name=="SPEMD") {
	// 		double power, powerFrom, powerTo, powerInc; 
	// double core, coreFrom, coreTo, coreInc;
			cout << "[" << parameter[i].name  << "]:" 
				<< "\n_From\t"
				<< parameter[i].centerXFrom << "\t" 
				<< parameter[i].centerYFrom << "\t"
				<< parameter[i].critRadFrom << "\t"
				<< parameter[i].powerFrom << "\t"
				<< parameter[i].eFrom << "\t" 
				<< parameter[i].PAFrom << "\t" 
				<< parameter[i].coreFrom << "\t" 
				<< "\n_To\t"
				<< parameter[i].centerXTo << "\t" 
				<< parameter[i].centerYTo << "\t"
				<< parameter[i].critRadTo << "\t"
				<< parameter[i].powerTo << "\t"
				<< parameter[i].eTo << "\t" 
				<< parameter[i].PATo << "\t" 
				<< parameter[i].coreTo << "\t" 
				<< "\n_Inc\t"
				<< parameter[i].centerXInc << "\t" 
				<< parameter[i].centerYInc << "\t"
				<< parameter[i].critRadInc << "\t"
				<< parameter[i].powerInc << "\t"
				<< parameter[i].eInc << "\t" 
				<< parameter[i].PAInc << "\t" 
				<< parameter[i].coreInc << "\t" 
				<< endl; 
		}
	}
	cout << "*******************************" << endl;
}


void MultModelParam::mix(int opt) { 
	
	
	/*	mixModel structure:    {name, paraList[0, 0, 0, 0, 0, 0, 0, 0]}
	*
	*/

	vector<vector<mixModels> > mix;
	for(int i=0; i<nLens; ++i) {

		if (parameter[i].name=="PTMASS") {
			vector<mixModels> v1;
            if (opt == 0) {
	        	for (double critRad = parameter[i].critRadFrom;	critRad <= parameter[i].critRadTo; critRad += parameter[i].critRadInc) {
	        		for (double centerX = parameter[i].centerXFrom;	centerX <= parameter[i].centerXTo; centerX += parameter[i].centerXInc) {
	        			for (double centerY = parameter[i].centerYFrom;	centerY <= parameter[i].centerYTo; centerY += parameter[i].centerYInc) {
	        						mixModels sModel("PTMASS");
	        						sModel.paraList[0] = critRad;
	        						sModel.paraList[1] = centerX;
	        						sModel.paraList[2] = centerY;
	        						sModel.paraList[3] = 0;
	        						sModel.paraList[4] = 0;
	        						sModel.paraList[5] = 0;  //parameter[i].core;
	        						v1.push_back(sModel);
	        			}
	        		}
	        	}
            } else if (opt >= 1) {
                mixModels sModel("PTMASS");
                sModel.paraList[0] = 0.;
                sModel.paraList[1] = 0.;
                sModel.paraList[2] = 0.;
                sModel.paraList[3] = 0.;
                sModel.paraList[4] = 0.;
                sModel.paraList[5] = 0.;
                v1.push_back(sModel);
                v1.push_back(sModel);
                v1.push_back(sModel);
                sModel.paraList[0] = parameter[i].critRadFrom;
                sModel.paraList[1] = parameter[i].centerXFrom;
                sModel.paraList[2] = parameter[i].centerYFrom;
                v1.push_back(sModel);
                v1.push_back(sModel);
                v1.push_back(sModel);
                v1.push_back(sModel);
                sModel.paraList[0] = parameter[i].critRadTo;
                sModel.paraList[1] = parameter[i].centerXTo;
                sModel.paraList[2] = parameter[i].centerYTo;
                v1.push_back(sModel);
                sModel.paraList[0] = 0.1 * (parameter[i].critRadTo - parameter[i].critRadFrom);
                sModel.paraList[1] = 0.1 * (parameter[i].centerXTo - parameter[i].centerXFrom);
                sModel.paraList[2] = 0.1 * (parameter[i].centerYTo - parameter[i].centerYFrom);
                v1.push_back(sModel);
            }
			mix.push_back(v1);

		} else if (parameter[i].name=="SIE") {
			vector<mixModels> v1;
            if (opt == 0) {
	        	for (double critRad = parameter[i].critRadFrom;	critRad <= parameter[i].critRadTo; critRad += parameter[i].critRadInc) {
	        		for (double centerX = parameter[i].centerXFrom;	centerX <= parameter[i].centerXTo; centerX += parameter[i].centerXInc) {
	        			for (double centerY = parameter[i].centerYFrom;	centerY <= parameter[i].centerYTo; centerY += parameter[i].centerYInc) {
	        				for (double e = parameter[i].eFrom;	e <= parameter[i].eTo; e += parameter[i].eInc) {
	        					for (double PA = parameter[i].PAFrom; PA <= parameter[i].PATo; PA += parameter[i].PAInc) {
	        						for(double core = parameter[i].coreFrom; core <= parameter[i].coreTo; core += parameter[i].coreInc) {
	        							mixModels sModel("SIE");
	        							sModel.paraList[0] = critRad;
	        							sModel.paraList[1] = centerX;
	        							sModel.paraList[2] = centerY;
	        							sModel.paraList[3] = e;
	        							sModel.paraList[4] = PA;
	        							sModel.paraList[5] = core;  //parameter[i].core;
	        							v1.push_back(sModel);
	        						}
	        					}
	        				}
	        			}
	        		}
	        	}
            } else if (opt >= 1) {
                mixModels sModel("SIE");
                sModel.paraList[0] = 0.;
                sModel.paraList[1] = 0.;
                sModel.paraList[2] = 0.;
                sModel.paraList[3] = 0.;
                sModel.paraList[4] = 0.;
                sModel.paraList[5] = 0.;
                v1.push_back(sModel);
                v1.push_back(sModel);
                v1.push_back(sModel);
                sModel.paraList[0] = 0.5 * (parameter[i].critRadTo + parameter[i].critRadFrom);
                sModel.paraList[1] = 0.5 * (parameter[i].centerXTo + parameter[i].centerXFrom);
                sModel.paraList[2] = 0.5 * (parameter[i].centerYTo + parameter[i].centerYFrom);
                sModel.paraList[3] = 0.5 * (parameter[i].eTo       + parameter[i].eFrom      );
                sModel.paraList[4] = 0.5 * (parameter[i].PATo      + parameter[i].PAFrom     );
                sModel.paraList[5] = 0.5 * (parameter[i].coreTo    + parameter[i].coreFrom   );
                v1.push_back(sModel);
                v1.push_back(sModel);
                v1.push_back(sModel);
                sModel.paraList[0] = parameter[i].critRadFrom;
                sModel.paraList[1] = parameter[i].centerXFrom;
                sModel.paraList[2] = parameter[i].centerYFrom;
                sModel.paraList[3] = parameter[i].eFrom;
                sModel.paraList[4] = parameter[i].PAFrom;
                sModel.paraList[5] = parameter[i].coreFrom;
                v1.push_back(sModel);
                sModel.paraList[0] = parameter[i].critRadTo;
                sModel.paraList[1] = parameter[i].centerXTo;
                sModel.paraList[2] = parameter[i].centerYTo;
                sModel.paraList[3] = parameter[i].eTo;
                sModel.paraList[4] = parameter[i].PATo;
                sModel.paraList[5] = parameter[i].coreTo;
                v1.push_back(sModel);
                sModel.paraList[0] = parameter[i].critRadInc;
                sModel.paraList[1] = parameter[i].centerXInc;
                sModel.paraList[2] = parameter[i].centerYInc;
                sModel.paraList[3] = parameter[i].eInc;
                sModel.paraList[4] = parameter[i].PAInc;
                sModel.paraList[5] = parameter[i].coreInc;
                v1.push_back(sModel);
            }

			mix.push_back(v1);
		}

	// double massScale, massScaleFrom, massScaleTo, massScaleInc;   
	// double radScale, radScaleFrom, radScaleTo, radScaleInc; 
		else if (parameter[i].name=="NFW") {
			vector<mixModels> v1;
            if (opt == 0) {
	        	for (double massScale = parameter[i].massScaleFrom;	massScale <= parameter[i].massScaleTo; massScale += parameter[i].massScaleInc) {
	        		for (double centerX = parameter[i].centerXFrom;	centerX <= parameter[i].centerXTo; centerX += parameter[i].centerXInc) {
	        			for (double centerY = parameter[i].centerYFrom;	centerY <= parameter[i].centerYTo; centerY += parameter[i].centerYInc) {
	        				for (double e = parameter[i].eFrom;	e <= parameter[i].eTo; e += parameter[i].eInc) {
	        					for (double PA = parameter[i].PAFrom; PA <= parameter[i].PATo; PA += parameter[i].PAInc) {
	        						for(double radScale = parameter[i].radScaleFrom; radScale <= parameter[i].radScaleTo; radScale += parameter[i].radScaleInc) {
	        							mixModels sModel("NFW");
	        							sModel.paraList[0] = massScale;
	        							sModel.paraList[1] = centerX;
	        							sModel.paraList[2] = centerY;
	        							sModel.paraList[3] = e;
	        							sModel.paraList[4] = PA;
	        							sModel.paraList[5] = radScale;  //parameter[i].core;
	        							v1.push_back(sModel);
	        						}
	        					}
	        				}
	        			}
	        		}
	        	}
            } 

			mix.push_back(v1);
		}

		else if (parameter[i].name=="SPEMD") {
			vector<mixModels> v1;
            if (opt == 0) {
            	for (double critRad = parameter[i].critRadFrom;	critRad <= parameter[i].critRadTo; critRad += parameter[i].critRadInc) {
	        		for (double centerX = parameter[i].centerXFrom;	centerX <= parameter[i].centerXTo; centerX += parameter[i].centerXInc) {
	        			for (double centerY = parameter[i].centerYFrom;	centerY <= parameter[i].centerYTo; centerY += parameter[i].centerYInc) {
	        				for (double e = parameter[i].eFrom;	e <= parameter[i].eTo; e += parameter[i].eInc) {
	        					for (double PA = parameter[i].PAFrom; PA <= parameter[i].PATo; PA += parameter[i].PAInc) {
	        						for(double core = parameter[i].coreFrom; core <= parameter[i].coreTo; core += parameter[i].coreInc) {
	        							for(double power = parameter[i].powerFrom; power <= parameter[i].powerTo; power += parameter[i].powerInc) {

	        								mixModels sModel("SPEMD");
	        								sModel.paraList[0] = critRad;
	        								sModel.paraList[1] = centerX;
	        								sModel.paraList[2] = centerY;
	        								sModel.paraList[3] = e;
	        								sModel.paraList[4] = PA;
	        								sModel.paraList[5] = core;  //parameter[i].core;
	        								sModel.paraList[6] = power; 
	        								v1.push_back(sModel);
	        							}
	        						}
	        					}
	        				}
	        			}
	        		}
	        	}
	        	
            } 

			mix.push_back(v1);
		}


		//     	tempParam.centerXFrom 		= stof(items[0]);
				// tempParam.centerXTo 		= stof(items[1]);
				// tempParam.centerXInc 		= stof(items[2]);
				// tempParam.centerYFrom 		= stof(items[3]);
				// tempParam.centerYTo 		= stof(items[4]);
				// tempParam.centerYInc 		= stof(items[5]);
				// tempParam.kapFrom 			= stof(items[6]);
				// tempParam.kapTo   			= stof(items[7]);
				// tempParam.kapInc  			= stof(items[8]);
				// tempParam.eFrom 			= stof(items[9]);
				// tempParam.eTo   			= stof(items[10]);
				// tempParam.eInc  			= stof(items[11]);
				// tempParam.PAFrom			= stof(items[12]);
				// tempParam.PATo 				= stof(items[13]);
				// tempParam.PAInc 			= stof(items[14]);
				// tempParam.sersicScaleFrom 	= stof(items[15]);
				// tempParam.sersicScaleTo		= stof(items[16]);
				// tempParam.sersicScaleInc	= stof(items[17]);
				// tempParam.mFrom 			= stof(items[18]);
				// tempParam.mTo   			= stof(items[19]);
				// tempParam.mInc  			= stof(items[20]);

		else if (parameter[i].name=="SERSIC") {
			vector<mixModels> v1;
            if (opt == 0) {
            	for (double kap = parameter[i].kapFrom;	kap <= parameter[i].kapTo; kap += parameter[i].kapInc) {
	        		for (double centerX = parameter[i].centerXFrom;	centerX <= parameter[i].centerXTo; centerX += parameter[i].centerXInc) {
	        			for (double centerY = parameter[i].centerYFrom;	centerY <= parameter[i].centerYTo; centerY += parameter[i].centerYInc) {
	        				for (double e = parameter[i].eFrom;	e <= parameter[i].eTo; e += parameter[i].eInc) {
	        					for (double PA = parameter[i].PAFrom; PA <= parameter[i].PATo; PA += parameter[i].PAInc) {
	        						for(double sersicScale = parameter[i].sersicScaleFrom; sersicScale <= parameter[i].sersicScaleTo; sersicScale += parameter[i].sersicScaleInc) {
	        							for(double  m= parameter[i].mFrom; m <= parameter[i].mTo; m += parameter[i].mInc) {

	        								mixModels sModel("SERSIC");
	        								sModel.paraList[0] = kap;
	        								sModel.paraList[1] = centerX;
	        								sModel.paraList[2] = centerY;
	        								sModel.paraList[3] = e;
	        								sModel.paraList[4] = PA;
	        								sModel.paraList[5] = sersicScale;  
	        								sModel.paraList[6] = m; 
	        								v1.push_back(sModel);
	        							}
	        						}
	        					}
	        				}
	        			}
	        		}
	        	}
	        	
            } 

			mix.push_back(v1);
		}


	}
	mix.resize(3);
    size_t ms1 = (mix[1].size() > 0) ? mix[1].size():1;
    size_t ms2 = (mix[2].size() > 0) ? mix[2].size():1;

    /* for maximum 3 models:  j, k, m */
    if (opt == 0) {
        for(size_t j=0; j<mix[0].size(); ++j) {
            for(size_t k=0; k<ms1; ++k ) {
                for(size_t m=0; m<ms2; ++m) {
                    vector<mixModels> v2;
                    v2.push_back(mix[0][j]);
                    if (nLens>1) v2.push_back(mix[1][k]);
                    if (nLens>2) v2.push_back(mix[2][m]);
                    mixAllModels.push_back(v2);
                }
            }
        }
    } else if (opt >= 1) {
        for(size_t j=0; j<mix[0].size(); ++j) {
            vector<mixModels> v2;
            for(int i=0; i<nLens; ++i) v2.push_back(mix[i][j]);
            mixAllModels.push_back(v2);
        }
    }
	nComb = mixAllModels.size();
    cout<< "nComb "<< nComb <<" "<<ms1<<" "<<ms2<<" "<< mix[0].size()<<endl;
}


vector<string> MultModelParam::printCurrentModels(int curr) {
		
	//cout <<  mixAllModels[curr].size() << endl; 
	vector<string> ret; 
	string modelsInRow; 
	string modelsInCol; 

	// return a string for "output.txt"; 
	modelsInCol += ( "[" + to_string(curr+1)+"/"+to_string(nComb)  + "]\n" ); 
	for(int i=0; i<mixAllModels[curr].size(); ++i) {
		modelsInCol += (mixAllModels[curr][i].name  + ":\t") ; 
		for (int j=0; j< 8; ++j) {
			modelsInCol += (to_string(mixAllModels[curr][i].paraList[j]) + "\t") ; 
			modelsInRow += (to_string(mixAllModels[curr][i].paraList[j]) + "\t")  ; 

		}
		modelsInCol += "\n";  
	} 
	modelsInCol += "\n";
	ret.push_back(modelsInRow); 
	ret.push_back(modelsInCol); 

	return ret; 
}


void Model::resetVectors(Conf* conf) {

	L.		setZero();
	M.		setZero(); 
	r.		setZero(); 
	new_r.	setZero(); 
	s.		setZero(); 
	phi.	setZero(); 
	//square_s.setZero(); 
	Ds.		setZero(); 
	Dphi.	setZero(); 
	Hs1.	setZero(); 
	Hs2.	setZero(); 
	HcH. 	setZero(); 
	Hphi.	setZero(); 
	HphiH.	setZero(); 
	T.		setZero(); 
	RtR.	setZero(); 
	H0.		setZero(); 
	H1.		setZero(); 
	H2.		setZero(); 
	H0H.	setZero(); 
	H_zero.	setZero(); 



	vector<vector<normVec> > normV;
	vector<normVec> meanNormV;



	//param.parameter.	clear(); 
	srcPosXListPixel.	clear(); 
	srcPosYListPixel.	clear(); 
	srcPosXList.		clear(); 
	srcPosYList.		clear(); 
	pDeltaX.			clear(); 
	pDeltaY.			clear(); 
	critical.			clear();

	res_img.			clear(); 
	res_full_img.		clear(); 
	simple_res_img.		clear(); 
	mod_img.			clear(); 

	L.	 reserve(Eigen::VectorXi::Constant(length,3));
	Ds.  reserve(2*length);
	Dphi.reserve(2*length);
	
	Hs1. reserve(Eigen::VectorXi::Constant(length,10));
	Hs2. reserve(Eigen::VectorXi::Constant(length,10));
	HcH. reserve(Eigen::VectorXi::Constant(length,1 )); 
	RtR. reserve(200*length);
	T.   reserve(length);

	
	H_zero.reserve(Eigen::VectorXi::Constant(length, 2)); 
	H_grad.reserve(Eigen::VectorXi::Constant(length, 5)); 
	H_curv.reserve(Eigen::VectorXi::Constant(length, 10)); 
	
	Hessian_grad.reserve(Eigen::VectorXi::Constant(conf->srcSize[0]*conf->srcSize[1], 6)); 

}

void Model::copyParam(Conf* conf, int i) {
    resetVectors(conf);
    for (int j=0; j<param.nLens; ++j) {
        if (param.mixAllModels[i][j].name == "PTMASS") {
            param.parameter[j].critRad = param.mixAllModels[i][j].paraList[0];
            param.parameter[j].centerX = param.mixAllModels[i][j].paraList[1];
            param.parameter[j].centerY = param.mixAllModels[i][j].paraList[2];
        } else if (param.mixAllModels[i][j].name == "SIE") {
            param.parameter[j].critRad = param.mixAllModels[i][j].paraList[0];
            param.parameter[j].centerX = param.mixAllModels[i][j].paraList[1];
            param.parameter[j].centerY = param.mixAllModels[i][j].paraList[2];
            param.parameter[j].e       = param.mixAllModels[i][j].paraList[3];
            param.parameter[j].PA      = param.mixAllModels[i][j].paraList[4];
            param.parameter[j].core    = param.mixAllModels[i][j].paraList[5];
        } else if (param.mixAllModels[i][j].name == "NFW") {
            param.parameter[j].massScale = param.mixAllModels[i][j].paraList[0];
            param.parameter[j].centerX   = param.mixAllModels[i][j].paraList[1];
            param.parameter[j].centerY   = param.mixAllModels[i][j].paraList[2];
            param.parameter[j].e         = param.mixAllModels[i][j].paraList[3];
            param.parameter[j].PA        = param.mixAllModels[i][j].paraList[4];
            param.parameter[j].radScale  = param.mixAllModels[i][j].paraList[5];
        } else if (param.mixAllModels[i][j].name == "SERSIC") {
            param.parameter[j].kap 	  		= param.mixAllModels[i][j].paraList[0];
            param.parameter[j].centerX 		= param.mixAllModels[i][j].paraList[1];
            param.parameter[j].centerY 		= param.mixAllModels[i][j].paraList[2];
            param.parameter[j].e       		= param.mixAllModels[i][j].paraList[3];
            param.parameter[j].PA      		= param.mixAllModels[i][j].paraList[4];
            param.parameter[j].sersicScale  = param.mixAllModels[i][j].paraList[5];
            param.parameter[j].m   			= param.mixAllModels[i][j].paraList[6];
        } else if (param.mixAllModels[i][j].name == "SPEMD") {
            param.parameter[j].critRad = param.mixAllModels[i][j].paraList[0];
            param.parameter[j].centerX = param.mixAllModels[i][j].paraList[1];
            param.parameter[j].centerY = param.mixAllModels[i][j].paraList[2];
            param.parameter[j].e       = param.mixAllModels[i][j].paraList[3];
            param.parameter[j].PA      = param.mixAllModels[i][j].paraList[4];
            param.parameter[j].core    = param.mixAllModels[i][j].paraList[5];
            param.parameter[j].power   = param.mixAllModels[i][j].paraList[6];
        }
    }
}

void Model::copyParam(int i1, int i2) {
    for(int j=0; j<param.nLens; ++j) {
        for (size_t k=0; k<param.mixAllModels[0][j].paraList.size(); ++k) {
            param.mixAllModels[i2][j].paraList[k] = param.mixAllModels[i1][j].paraList[k];
        }
    }
}

vector<vector<double> > getCritCausticFine(vector<double> xPosListArc, vector<double> yPosListArc, Conf* conf, MultModelParam * param, int level) {
	// level = 5; 
	vector<double> invMag;
	map<pair<int, int>,int> posMap;
	vector<double> srcPos, pDeltaX, pDeltaY; 
	vector<double> srcPosXListArc; 
	vector<double> srcPosYListArc; 

	double curr_res = conf->imgRes/level; 

	vector<vector<double> > ret; 

	vector<double> imgXList; 
	vector<double> imgYList; 

	double side = 1.0/level;   // unit is pixel; 

	int index = 0 ; 
	for (int i=0; i<xPosListArc.size(); ++i) {
		int old_imgX = int(round(xPosListArc[i] / conf->imgRes + conf->imgXCenter )); 
		int old_imgY = int(round(yPosListArc[i] / conf->imgRes + conf->imgYCenter )); 

		for(int j=0; j<level; ++j) {		
			for(int k=0; k<level; ++k) {
			double new_posX = (j * side) * curr_res + xPosListArc[i];   // in arcsecond
			double new_posY = (k * side) * curr_res + yPosListArc[i];   // in arcsecond

			double defX = 0 ;
			double defY = 0 ;  
			srcPos = Model::getDeflectionAngle(conf, new_posX, new_posY, &defX, &defY, param);
			pDeltaX.push_back(defX);    // in arcsec; 
			pDeltaY.push_back(defY);		
			srcPosXListArc.push_back(srcPos[0]); 
			srcPosYListArc.push_back(srcPos[1]);

			int imgX = old_imgX * level + k; 
			int imgY = old_imgY * level + j; 
			imgXList.push_back(imgX); 
			imgYList.push_back(imgY); 
			posMap[make_pair(imgX, imgY)] = index;
			++index; 

			}
		}
	}		

	map<pair<int, int>,int>::iterator left, right, up, down;
	vector<double> w, w5;

	vector<double> a11, a12, a21, a22;
	vector<double> new_xPosListArc; 
	vector<double> new_yPosListArc; 
	vector<double> new_xSrcPosListArc; 
	vector<double> new_ySrcPosListArc; 


	double h = curr_res;
	for (int i=0; i< index; ++i) {
		left  = posMap.find(make_pair(imgXList[i]-1, imgYList[i]));
		right = posMap.find(make_pair(imgXList[i]+1, imgYList[i]));
		up    = posMap.find(make_pair(imgXList[i], imgYList[i]+1));
		down  = posMap.find(make_pair(imgXList[i], imgYList[i]-1));

		if(left!=posMap.end() && up!=posMap.end() && down!=posMap.end() && right!=posMap.end()) {

			int iLeft = left->second;
			int iUp   = up  ->second;
			int iDown = down->second;
			int iRight= right->second;

			a11.push_back((pDeltaX[iRight]-pDeltaX[iLeft])/(2*h));
			a12.push_back((pDeltaX[iUp]-pDeltaX[iDown])/(2*h));
			a21.push_back((pDeltaY[iRight]-pDeltaY[iLeft])/(2*h));
			a22.push_back((pDeltaY[iUp]-pDeltaY[iDown])/(2*h));
		}
		else {
			a11.push_back(0);
			a12.push_back(0);
			a21.push_back(0);
			a22.push_back(0);
		}
		invMag.push_back(1.0-(a11[i]+a22[i])+a11[i]*a22[i]-a12[i]*a21[i]);
	}
	// update inverse magnification
	int sign_t = 0;
	

	for (int i=0; i<index; ++i) {
			
			left  = posMap.find(make_pair(imgXList[i]-1, imgYList[i]));
			right = posMap.find(make_pair(imgXList[i]+1, imgYList[i]));
			up    = posMap.find(make_pair(imgXList[i], imgYList[i]+1));
			down  = posMap.find(make_pair(imgXList[i], imgYList[i]-1));

			if(left!=posMap.end() && up!=posMap.end() && down!=posMap.end() && right!=posMap.end()) {
				int iLeft = left->second;
				int iUp   = up  ->second;
				int iDown = down->second;
				int iRight= right->second;
				sign_t = sign(invMag[i])*(sign(invMag[iLeft]) + sign(invMag[iRight]) + sign(invMag[iUp]) + sign(invMag[iDown]));
	
			}
			else
				sign_t = 5;  // Assign a value bigger than 4;

			if(sign_t<4 ) { //} && distSqure > 50*50) {
				new_xPosListArc.push_back((imgXList[i]-conf->imgXCenter*level)* curr_res); 
				new_yPosListArc.push_back((imgYList[i]-conf->imgYCenter*level)* curr_res);
				new_xSrcPosListArc.push_back(srcPosXListArc[i]); 
				new_ySrcPosListArc.push_back(srcPosYListArc[i]); 


			}
	}

	// Center:  (0, 0)  in arcsecond
	ret.push_back(new_xPosListArc); 
	ret.push_back(new_yPosListArc); 
	ret.push_back(new_xSrcPosListArc); 
	ret.push_back(new_ySrcPosListArc); 
	return ret; 
}


vector<Image* > getCritCaustic(Conf* conf, MultModelParam * param) {
	// automatic using  full Image ( no region files ); 
	
	vector<Image*> ret; 
	vector<double> critical, xList, yList;
	vector<double> caustic, srcXList, srcYList; 
	vector<double> srcPosXListPixel; 
	vector<double> srcPosYListPixel; 
	int level =  conf->causticLevel; 

	Image* dataImage = new Image(conf->imageFileName); 
	dataImage->updateFilterImage("whatever..", 0);   // without any region file --- using all data points; 
	int length = dataImage->data.size();	

	vector<double> xPosListArc; 
	vector<double> yPosListArc; 

	for(int i=0; i<length; ++i) {
		xPosListArc.push_back((dataImage->xList[i]-conf->imgXCenter)*conf->imgRes); 
		yPosListArc.push_back((dataImage->yList[i]-conf->imgYCenter)*conf->imgRes); 
	}

	// Basic critical searching; 
	vector<vector<double> >  critXY  = getCritCausticFine(xPosListArc, yPosListArc, conf, param, 1);
	// Finer critical searching with higher level; 
	if (level > 1 ) 
		critXY = getCritCausticFine(critXY[0], critXY[1], conf, param, level); 

	double newResolution = conf->imgRes/level; 
	double newSrcResolution = conf->srcRes/level; 

	for(int i=0; i<critXY[0].size(); ++i) {


		// Image coordinate (in pixel)
		xList.push_back(critXY[0][i]/newResolution + conf->imgXCenter*level +1); 
		yList.push_back(critXY[1][i]/newResolution + conf->imgYCenter*level +1); 

		srcXList.push_back(critXY[2][i]/newSrcResolution + conf->srcXCenter*level +1); 
		srcYList.push_back(critXY[3][i]/newSrcResolution + conf->srcYCenter*level +1); 
		critical.push_back(1); 

	}
	createDs9Contour(&xList,    &yList,    level, conf->contourCritName) ; 
	createDs9Contour(&srcXList, &srcYList, level, conf->contourCausName) ; 

	/// modify ending: 
	Image* critImg = new Image(xList, yList, &critical, conf->imgSize[0]*level, conf->imgSize[1]*level, conf->bitpix);
	Image* causImg = new Image(srcXList, srcYList, &critical, conf->srcSize[0]*level, conf->srcSize[1]*level, conf->bitpix);

	ret.push_back(critImg); 
	ret.push_back(causImg); 
 	return ret; 
}


void createDs9Contour(vector<double>* xList, vector<double>* yList, double level, string contourFileName) {
	//   contour in image coordinate; 
	ofstream contourFile(contourFileName);
	map<pair<int, int>, int> posMap;
	for (int i=0; i<xList->size(); ++i) {
		posMap[make_pair(int((*xList)[i]), int((*yList)[i]))] = i; 
	}
	vector <map<pair<int, int>,int>::iterator> n(9); //  n1, n2, n3, n4, n5, n6, n7, n8, n9; 
	/* find all the neighbors of each points; 
	 *  n0, || n1,  n2
	 * 	n3, || n4,  n5
	 *	n6, || n7,  n8
	 */
	for (int i=0; i<xList->size(); ++i) {
			int x = int((*xList)[i]); 
			int y = int((*yList)[i]); 

			n[0] = posMap.find(make_pair(x-1, y+1));
			n[1] = posMap.find(make_pair(x  , y+1));
			n[2] = posMap.find(make_pair(x+1, y+1));
			n[3] = posMap.find(make_pair(x-1, y  ));
			n[4] = posMap.find(make_pair(x  , y  ));    // self; 
			n[5] = posMap.find(make_pair(x+1, y  ));
			n[6] = posMap.find(make_pair(x-1, y-1));
			n[7] = posMap.find(make_pair(x  , y-1));
			n[8] = posMap.find(make_pair(x+1, y-1));


			for (int j=0; j<9; ++j) {
				if(n[j]!=posMap.end()) {
					int ind = n[j]->second; 
					contourFile << to_string((*xList)[i]/level)   << "\t" << to_string((*yList)[i]/level)   << endl; 
 					contourFile << to_string((*xList)[ind]/level) << "\t" << to_string((*yList)[ind]/level) << endl; 
 					contourFile << endl; 
				}

			}

		}
 	contourFile.close(); 
}



Image* createLensImage(Conf* conf, MultModelParam * param) {

	double level = 1; 
	vector<int> xList;  // in pixel, the same size as 'dataImage'; 
	vector<int> yList;
	vector<double> val;
	for(int i=0; i<conf->imgSize[1]*level; ++i) {
		for(int j=0; j<conf->imgSize[0]*level; ++j) {
			xList.push_back(j); 
			yList.push_back(i); 
			val.push_back(0); 
		}
	}






	for(int i=0; i< param->nLens; ++i) {

		if(param->parameter[i].name =="PTMASS") {
		}
		if (param->parameter[i].name =="SIE") {

			double centerX = param->parameter[i].centerX;
			double centerY = param->parameter[i].centerY; 
			double critRad = param->parameter[i].critRad; 
			double q = 1-param->parameter[i].e; 
			double PA = param->parameter[i].PA/180 * M_PI + 0.5 * M_PI; 
			double core = param->parameter[i].core; 
			if (q==1)  q=0.999; 
			for(int j=0; j<xList.size(); ++j) {
				double x = (xList[j] - (centerX + conf->imgXCenter))*conf->imgRes  ; 
				double y = (yList[j] - (centerY + conf->imgYCenter))*conf->imgRes  ; 
				double new_x = x*cos(PA) + y*sin(PA); 
				double new_y = x*sin(PA) - y*cos(PA);  
				double coeff =1 ; // sqrt(q/(1-q*q)); 
				val[j] += critRad*coeff /sqrt(new_x*new_x*q + new_y*new_y/q + core*core); 
			}
		}
		if(param->parameter[i].name =="NFW") {
		}

		if(param->parameter[i].name =="SERSIC") {
		}
		




	}


	if(0) {
				// model1 ; 
				double PA = 115.0/180.0*M_PI + 0.5*M_PI; 
				double q = 1 - 0.1; 
				double kappa = 1.0; 
				double Re = 108*0.04; 
				double m = 2.74; 

				double bn = 1.9992*m - 0.3271; 

				for(int j=0; j<xList.size(); ++j) {
					double x = (xList[j] - (0 + conf->imgXCenter))*conf->imgRes  ; 
					double y = (yList[j] - (0 + conf->imgYCenter))*conf->imgRes  ; 
					double new_x = x*cos(PA) + y*sin(PA); 
					double new_y = x*sin(PA) - y*cos(PA); 
					double new_r = sqrt(new_x*new_x*q + new_y*new_y/q) ; 
					val[j] += 7.8* exp(-bn*(pow(new_r/Re, 1.0/m)-1)) / 163116.0; 
				}

				double sum = 0; 
				for(int j=0; j<xList.size(); ++j) {
					sum += val[j]; 

				}	
				cout << "sum of model 1: " << sum  << endl; 


				PA = 175.0/180.0*M_PI + 0.5*M_PI; 
				q = 1 - 0.4; 
				kappa = 1.0; 
				Re = 6.74*0.04; 
				m = 0.85; 

				bn = 1.9992*m - 0.3271; 

				for(int j=0; j<xList.size(); ++j) {
					double x = (xList[j] - (0 + conf->imgXCenter))*conf->imgRes  ; 
					double y = (yList[j] - (0 + conf->imgYCenter))*conf->imgRes  ; 
					double new_x = x*cos(PA) + y*sin(PA); 
					double new_y = x*sin(PA) - y*cos(PA); 
					double new_r = sqrt(new_x*new_x*q + new_y*new_y/q) ; 
					val[j] += 1.0* exp(-bn*(pow(new_r/Re, 1.0/m)-1)) / 507.0; 
				}

				sum = 0; 
				for(int j=0; j<xList.size(); ++j) {
					sum += val[j]; 

				}	
				cout << "sum of model 1 & 2: " << sum << endl ; 

			
		}	
	Image* lensImg = new Image(xList, yList, &val, conf->imgSize[0]*level, conf->imgSize[1]*level, conf->bitpix);

	return lensImg; 


}



void writeSrcModResImage(Model* model, Image* dataImage, Conf* conf, string fileName, string dir) {

	// Part I:  output 'src', 'mod', and 'res' images: 
	// Part II: output 'crit', 'caus' and 'lens' images; 
	vector<double> s = eigenV_to_cV(&model->s); 

	////////    output image by separating sourse; 
	// src1 and src2: 


	//  For each individual sources; 

	int k = model->k ;
	vector<vector<double>> sList(k, vector<double>(s.size(), 0)); 



	vector<double> s0(s.size(), 0); 
	vector<double> s1(s.size(), 0); 
	for(int i=0; i<s.size(); ++i) {
		for(int j=0; j<k; ++j) {
			if(model->Kmeans_group[i]==j) 
				sList[j][i] = s[i]; 
			
		}
	}

	// write src1: 

	for(int i=0; i<k; ++i) {
		unique_ptr<Image> individualSrc (new Image(model->srcPosXListPixel, model->srcPosYListPixel, &sList[i], conf->srcSize[0], conf->srcSize[1], conf->bitpix));		
		unique_ptr<Image> individualMod (new Image(dataImage->xList, dataImage->yList, &sList[i], conf->imgSize[0], conf->imgSize[1], conf->bitpix));
		individualSrc -> writeToFile (dir + "img_src_" + fileName + "_" + to_string(i) + ".fits"); //, conf->back_mean, conf->back_std);
		individualMod -> writeToFile (dir + "img_mod_" + fileName + "_" + to_string(i) + ".fits");
	}

	// unique_ptr<Image> srcImg0 (new Image(model->srcPosXListPixel, model->srcPosYListPixel, &s0, conf->srcSize[0], conf->srcSize[1], conf->bitpix));		
	// unique_ptr<Image> modImg0 (new Image(dataImage->xList, dataImage->yList, &s0, conf->imgSize[0], conf->imgSize[1], conf->bitpix));


	// unique_ptr<Image> srcImg1 (new Image(model->srcPosXListPixel, model->srcPosYListPixel, &s1, conf->srcSize[0], conf->srcSize[1], conf->bitpix));		
	// unique_ptr<Image> modImg1 (new Image(dataImage->xList, dataImage->yList, &s1, conf->imgSize[0], conf->imgSize[1], conf->bitpix));

	// srcImg0 -> writeToFile (dir + "img_src_" + fileName + "_0.fits", conf->back_mean, conf->back_std);
	// modImg0 -> writeToFile (dir + "img_mod_" + fileName + "_0.fits");

	// srcImg1 -> writeToFile (dir + "img_src_" + fileName + "_1.fits", conf->back_mean, conf->back_std);
	// modImg1 -> writeToFile (dir + "img_mod_" + fileName + "_1.fits");







	// Output fits for overall images ; 
	unique_ptr<Image> srcImg (new Image(model->srcPosXListPixel, model->srcPosYListPixel, &s, conf->srcSize[0], conf->srcSize[1], conf->bitpix));		
	unique_ptr<Image> modImg  (new Image(dataImage->xList, dataImage->yList, &s, conf->imgSize[0], conf->imgSize[1], conf->bitpix));
	Image* resImg = model->getFullResidual(dataImage);
	//srcImg -> writeToFile (dir + "img_src_" + fileName + ".fits", conf->back_mean, conf->back_std ) ;
	srcImg -> writeToFile (dir + "img_src_" + fileName + ".fits"); //, conf->back_mean, conf->back_std);
	resImg -> writeToFile (dir + "img_res_" + fileName + ".fits");
	modImg -> writeToFile (dir + "img_mod_" + fileName + ".fits");
	//delete srcImg, modImg, resImg; 

	vector<Image* > curve =  getCritCaustic(conf, &model->param); 
	Image* critImg = curve[0]; 
	Image* causImg = curve[1]; 
	Image* lensImg = createLensImage(conf, &model->param); 
	lensImg -> writeToFile(dir + "img_lens_" + fileName + ".fits");
	critImg -> writeToFile(dir + "img_crit_" + fileName + ".fits");
	causImg -> writeToFile(dir + "img_caus_" + fileName + ".fits");
	//delete critImg, causImg, lensImg ; 
}







Model::~Model() {
	// TODO Auto-generated destructor stub
}


#if 0
#endif 







/*
 * parafit2.cpp
 *
 *  Created on: Apr 22, 2016
 *      Author: En-Hsin Peng
 */


#include "Model.h"
#include "Image.h"
#include <iostream>
#include "commons.h"
#include <fstream>
#include <limits>
#include <iomanip>
#include "parafit.h"
#include "mc.h"
#include <omp.h>

using namespace std;



Model *model;
#pragma omp threadprivate(model)
void mcFitGW(Conf* conf, MultModelParam param_old, vector<Image*> dataImageList, string dir, string outputFileName) {
    size_t nLoops(conf->nLoops), lag(10), thin(15), writeChkpt(10), iter(0);
    double lambdaS = conf->srcRegLevel;
    size_t nWalkers(conf->nWalkers);
    auto objective = std::bind(getLogProb, std::placeholders::_1, dataImageList[0], conf);

    #pragma omp parallel
    model = new Model(conf, param_old, lambdaS);

    MC mc(model, conf, objective, nWalkers, outputFileName, iter);


    nLoops += iter;
    for (size_t loop=iter; loop<nLoops; ++loop) {
        #pragma omp parallel for
        for (size_t m=0; m<nWalkers; ++m) {
            mc.stretchMove(model,m);
        }
        if (loop % lag == 0) {
            mc.writeOutput(loop, thin);
            if (loop % writeChkpt == 0) mc.checkPoint(loop);
        }
    }
    mc.writeOutput(nLoops);
    mc.checkPoint(nLoops);
    mc.copyParam(model->param);
    model->copyParam(conf, 5);
    vector<double> bestChi = getPenalty(model, dataImageList[0], conf);
    cout << "best chi: " << bestChi[0] << " dof: " << dataImageList[0]->dataList.size() << " reg: "<< bestChi[1] << " total:" << bestChi[2] << endl;

    writeSrcModResImage(model, dataImageList[0], conf, "mcgw_" + to_string(conf->seed), dir) ;

    // Print out the best model :
    cout << "************************\nThe best models : "<< endl;
    cout << model->param.printCurrentModels(5).at(1);
    cout << "************************\n" << endl;
    cout << model->param.printCurrentModels(6).at(1);
    cout << model->param.printCurrentModels(7).at(1);
    cout << model->param.printCurrentModels(8).at(1);
    #pragma omp parallel
    delete model;
}


double getLogProb(Model* model, Image* dataImage, Conf* conf) {

    model->copyParam(conf, 3);
    model->updatePosMapping(dataImage, conf);
    model->updateCompactMatrix(dataImage);
    model->update_H_zero(conf);
    model->updateLensAndRegularMatrix(dataImage, conf);
    model->solveSource(&dataImage->invC, &dataImage->d, conf->srcRegType);

    vec &s = model->s;
    vec res = ( model->L * s - dataImage->d);
    vec chi2 = res.transpose() *  dataImage->invC * res * model->lambdaC* model->lambdaC  ;
    vec srcR = s  .transpose() *  model->REG      * s   * model->lambdaS* model->lambdaS  ;

    //double chi2 = (res.cwiseProduct(dataImage->invSigma)).squaredNorm()*model->lambdaC* model->lambdaC;
    double lp = srcR[0] + chi2[0];
    if (conf->verbose || std::isnan(lp)) {
        cout << "Penalty " << chi2[0] << " " << srcR[0] << " " << s.norm() << " " << model->REG.norm() << " " << res.norm() << endl;
        cout << model->L.norm() << " " << dataImage->d.norm() << " " << dataImage->invC.norm() << endl;
        if (std::isnan(lp)) return -1.0;
    }
    return lp;

}


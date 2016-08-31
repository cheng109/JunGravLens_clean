/*
 * mc.h
 *
 *  Created on: Apr, 2016
 *      Author: En-Hsin Peng
 */

#ifndef MC_H_
#define MC_H_

#include <vector>
#include <random>
#include <functional>
#include <cmath>
#include "Model.h"

class MC {
public:
    MC(MultModelParam &param, unsigned seed);
    MC(Model* model, Conf* conf, const std::function<double(Model*)> &objective,
            size_t n, string outputFileName, size_t &iter);
    MC(Model* model, Conf* conf, const std::function<double(Model*)> &objective,
            size_t n, string outputFileName);
    double gauss(double sigma, double mu, double x);
    void makeCgauss();
    double cgauss();
    double random();
    double stepPar(MultModelParam &param, int &j, int &k);
    void stepPar(MultModelParam &param);
    double stepPar(vector<vec> &src, double cfac, size_t &iter);
    void checkPoint(string outputFileName, MultModelParam &param, bool moveAll, double L, double LMax);
    void checkPoint(size_t loop);
    void load(string fileName, MultModelParam &param, bool moveAll, double &L, double &LMax);
    void load(size_t &loop);
    void printIterNum(bool moveAll);
    void stretchMove(Model *model, size_t kk);
    void writeOutput(size_t loop);
    void writeOutput(size_t loop, int thin);
    void copyParam(MultModelParam &param);
    void startGA();
    void startNM();
    void evaluate(Model *model, size_t m);
    void elitism();
    double getRMin(){return RMin;};
    void nelmin(Model *model, size_t m);

private:
    std::vector<double> cgArr;
    std::mt19937 rng_engine;
    //std::default_random_engine rng_engine;
    std::function<double()> rng;
    std::function<double()> normal;
    double cfac;
    vector<vector<size_t>> freePar, iter;
    size_t iters, nAccept;

    void setRng(unsigned seed);
    void setParam(MultModelParam &param);
    std::function<double(Model*)> objective;

    //GW
    size_t nFreePar;
    size_t nWalker;
    double stepA;
    double RMin;
    vector<vector<double>> par, bound;
    vector<double> bestPar, R0;
    void setupGW(MultModelParam &param, size_t n);
    ofstream output;
    string chkptFileName;
    bool runGA;
    //GA
    size_t nGens;
    vector<size_t> index;
    double pCrossOver, pMutation;
    vector<double> bestParPrev;
    vector<vector<double>> parSum;
    double RMinPrev;
    void setupGA();
    //Simplex
    double *step;
    double evaluate(Model *model, double par[]);
    void nelmin (Model *model, double start[], double xmin[],
            double *ynewlo, double reqmin, int konvge, int kcount,
            int *icount, int *numres, int *ifault);
};

#endif /* MC_H_ */

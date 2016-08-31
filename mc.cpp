/*
 * mc.cpp
 *
 *  Created on: Apr 23, 2016
 *      Author: En-Hsin Peng
 */

#include "mc.h"
#include <iostream>
#include <iomanip>

void MC::setRng(unsigned seed) {
    rng_engine.seed(seed);
    rng = std::bind(std::uniform_real_distribution<double>(0.,1.), std::ref(rng_engine));
    normal = std::bind(std::normal_distribution<double>(0.,1.), std::ref(rng_engine));
    //makeCgauss();
}

void MC::setParam(MultModelParam &param) {
    freePar.resize(param.nLens);
    iter.resize(param.nLens);
    nFreePar=0;
    for(int j=0; j<param.nLens; ++j) {
        iter[j].resize(param.mixAllModels[0][j].paraList.size(),0);
        for (size_t k=0; k<param.mixAllModels[0][j].paraList.size(); ++k) {
            if (param.mixAllModels[6][j].paraList[k] < param.mixAllModels[7][j].paraList[k]) {
                param.mixAllModels[4][j].paraList[k] = param.mixAllModels[6][j].paraList[k] +
                                random()*(param.mixAllModels[7][j].paraList[k]-param.mixAllModels[6][j].paraList[k]);
                freePar[j].push_back(k);
                nFreePar++;
            }
        }
    }
    cout << "total par " << nFreePar<< endl;
    iters=0;
    nAccept=0;
}
MC::MC(MultModelParam &param, unsigned seed) {
    setRng(seed);
    cfac=1.0;
    setParam(param);
}

//Simplex
MC::MC(Model* model, Conf* conf, const std::function<double(Model*)> &objective, size_t n, string outputFileName) {
    setRng(conf->seed);
    setParam(model->param);
    setupGW(model->param, n);
    this->objective = objective;
    step = new double[nFreePar];
    for (size_t c=0; c<nFreePar; ++c) step[c] = 0.1*(bound[1][c]-bound[0][c]);
    setupGA();
    pMutation = 0.4;

    if (conf->resume) {
        ifstream input(outputFileName.c_str());
        string line, token;
        if (input) {
            istringstream ss;
            while (getline(input, line)) {
                ss.str(line);
            }
            for (size_t l=0; l<nFreePar; ++l) {
                getline(ss, token, ' ');
                bestPar[l] = stod(token);
                cout << bestPar[l] << endl;
            }
            getline(ss, token, ' ');
            RMin = stod(token);
        }
        output.open(outputFileName, std::ofstream::out | std::ofstream::app);
    } else {
        output.open(outputFileName);
    }
}

//GW & GA
MC::MC(Model* model, Conf* conf,
     const std::function<double(Model*)> &objective, size_t n, string outputFileName, size_t &iter) {
    setRng(conf->seed);
    setParam(model->param);
    setupGW(model->param, n);
    if (conf->GA) {
        runGA = true;
        setupGA();
        chkptFileName = "mcga_chkpt_"+to_string(conf->seed)+".txt";
    } else {
        runGA = false;
        chkptFileName = "mcgw_chkpt_"+to_string(conf->seed)+".txt";
    }
    if (conf->resume) {
        load(iter);
        output.open(outputFileName, std::ofstream::out | std::ofstream::app);
    } else {
        output.open(outputFileName);
    }
    this->objective = objective;
}

void MC::setupGW(MultModelParam &param, size_t n) {
    nWalker = n;
    stepA = 2.;
    par.resize(n);
    par[0].resize(nFreePar);
    bound.resize(2);
    bound[0].resize(nFreePar);
    bound[1].resize(nFreePar);
    bestPar.resize(nFreePar,0);
    R0.resize(n,std::numeric_limits<double>::max());
    RMin = std::numeric_limits<double>::max();
    size_t c(0);
    double ss(5e-3);
    for (size_t j=0; j<freePar.size(); ++j) {
        for (auto k: freePar[j]) {
            bound[0][c] = param.mixAllModels[6][j].paraList[k];
            bound[1][c] = param.mixAllModels[7][j].paraList[k];
            par[0][c] = 0.5 * (bound[0][c] + bound[1][c]);
            par[0][c] += (random()-0.5) * 0.5 * (bound[1][c]-bound[0][c]);
            c++;
        }
    }
    for (size_t m=1; m<nWalker; ++m) {
        par[m].resize(nFreePar);
        for (size_t c=0; c<nFreePar; ++c) par[m][c] = par[m-1][c] + ss*(random()-0.5)*(bound[1][c]-bound[0][c]);
    }
}

void MC::setupGA() {
    pCrossOver = 0.6;
    pMutation = 0.4;
    index.resize(nWalker);
    parSum.resize(4);
    for (size_t m=0; m<4; ++m) parSum[m].resize(nFreePar,0);
    for (size_t m=0; m<nWalker; ++m) {
        index[m] = m;
        for (size_t c=0; c<nFreePar; ++c) par[m][c] = bound[0][c] + random()*(bound[1][c]-bound[0][c]);
    }
}

double MC::random() {
    return rng();
}

double MC::gauss(double sigma, double mu, double x) {
    double dx = x- mu;
    return exp(-dx*dx/(2*sigma*sigma))/(sigma*sqrt(2*M_PI));
}

void MC::makeCgauss() {
    cgArr.resize(100,0.);
    cgArr[0] = gauss(1., 0.,  -5.)*0.1;
    for (size_t i=1; i<100; ++i) {
        cgArr[i] = cgArr[i-1] + gauss(1., 0.,  -5. + i*0.1)*0.1;
    }
}

double MC::cgauss() {
    //double xv = 0.;
    //double r = rng();
    //for (size_t i=0; i<100; ++i) {
    //    if (cgArr[i]>r) {
    //        xv = (-5. + i*0.1) - 0.1*rng();
    //        break;
    //    }
    //}
    //return xv;
    return normal();
}


double MC::stepPar(MultModelParam &param, int &j, int &k) {
    double minSig = 1e-6;
    double eps=0.01;
    cfac *= (1+eps);
    if (cfac > 1e20) cfac = 1.0;

    // mixAllModels[0]: weight * par^2
    // mixAllModels[1]: weight
    // mixAllModels[2]: weight * par
    // mixAllModels[3]: new parameter
    // mixAllModels[4]: original parameter
    // mixAllModels[5]: best parameter
    // mixAllModels[6]: lower bound
    // mixAllModels[7]: upper bound
    // mixAllModels[8]: step
    j = random() * param.nLens;
    k = freePar[j][(int) (random() * freePar[j].size())];
    //cout << "parameter " << j << " " << k << endl;

    double stepSig = param.mixAllModels[8][j].paraList[k];
    double par0 = param.mixAllModels[4][j].paraList[k];

    param.mixAllModels[0][j].paraList[k] += cfac*par0*par0;
    param.mixAllModels[1][j].paraList[k] += cfac;
    param.mixAllModels[2][j].paraList[k] += cfac*par0;
    if (iter[j][k] > 3) {
        double sig = sqrt((1+1e-7)*param.mixAllModels[0][j].paraList[k]/param.mixAllModels[1][j].paraList[k]
                - pow(param.mixAllModels[2][j].paraList[k]/param.mixAllModels[1][j].paraList[k],2));
        //stepSig = 3.*2.38*sig/sqrt(freePar[j].size());
        stepSig = sig;
        if (std::isnan(stepSig)) {
            iter[j][k] = 0;
            param.mixAllModels[0][j].paraList[k] = 0.;
            param.mixAllModels[1][j].paraList[k] = 0.;
            param.mixAllModels[2][j].paraList[k] = 0.;
            stepSig = param.mixAllModels[8][j].paraList[k];
        }
    }
    if (stepSig < minSig) stepSig = minSig;
    double r = cgauss();
    param.mixAllModels[3][j].paraList[k] = par0 + r*stepSig;
    //std::cout << r << " " << stepSig << " " << std::endl;
    if (param.mixAllModels[3][j].paraList[k] < param.mixAllModels[6][j].paraList[k]
            || param.mixAllModels[3][j].paraList[k] > param.mixAllModels[7][j].paraList[k])
        param.mixAllModels[3][j].paraList[k] = param.mixAllModels[6][j].paraList[k] +
            random()*(param.mixAllModels[7][j].paraList[k]-param.mixAllModels[6][j].paraList[k]);

    iter[j][k]++;
    return stepSig;
}

void MC::stepPar(MultModelParam &param) {
    double minSig = 1e-6;
    double eps=0.01;
    cfac *= (1+eps);
    if (cfac > 1e20) cfac = 1.0;

    for(int j=0; j<param.nLens; ++j) {
        for (auto k: freePar[j]) {
            double stepSig = param.mixAllModels[8][j].paraList[k];
            double par0 = param.mixAllModels[4][j].paraList[k];
            param.mixAllModels[0][j].paraList[k] += cfac*par0*par0;
            param.mixAllModels[1][j].paraList[k] += cfac;
            param.mixAllModels[2][j].paraList[k] += cfac*par0;
            if (iters > 3) {
                double sig = sqrt((1+1e-7)*param.mixAllModels[0][j].paraList[k]/param.mixAllModels[1][j].paraList[k]
                        - pow(param.mixAllModels[2][j].paraList[k]/param.mixAllModels[1][j].paraList[k],2));
                stepSig = sig/sqrt(freePar[j].size());
                if (std::isnan(stepSig)) {
                    iters = 0;
                    param.mixAllModels[0][j].paraList[k] = 0.;
                    param.mixAllModels[1][j].paraList[k] = 0.;
                    param.mixAllModels[2][j].paraList[k] = 0.;
                    stepSig = param.mixAllModels[8][j].paraList[k];
                }
            }
            if (stepSig < minSig) stepSig = minSig;
            double r = cgauss();
            param.mixAllModels[3][j].paraList[k] = par0 + r*stepSig;
            //std::cout << r << " " << stepSig << " " << std::endl;
            if (param.mixAllModels[3][j].paraList[k] < param.mixAllModels[6][j].paraList[k]
                    || param.mixAllModels[3][j].paraList[k] > param.mixAllModels[7][j].paraList[k])
                param.mixAllModels[3][j].paraList[k] = param.mixAllModels[6][j].paraList[k] +
                    random()*(param.mixAllModels[7][j].paraList[k]-param.mixAllModels[6][j].paraList[k]);
        }
    }
    iters++;
}

double MC::stepPar(vector<vec> &src, double cfac, size_t &iter) {
    double minSig = 1e-6;
    double eps=0.0;
    cfac *= (1+eps);
    if (cfac > 1e20) {
        cfac = 1.0;
        for (int k=0; k<src[0].size(); ++k) {
            src[0](k) = 0.;
            src[1](k) = 0.;
            src[2](k) = 0.;
        }
    }

    for (int k=0; k<src[0].size(); ++k) {
        double stepSig = src[8](k);
        double par0 = src[4](k);

        src[0](k) += cfac*par0*par0;
        src[1](k) += cfac;
        src[2](k) += cfac*par0;
        if (iter > 3) {
            double sig = sqrt(src[0](k)/src[1](k)
                    - pow(src[2](k)/src[1](k),2));
            stepSig = 3.*2.38*sig/sqrt(src[0].size());
            if (std::isnan(stepSig)) {
                iter=0;
                src[0](k) = 0.;
                src[1](k) = 0.;
                src[2](k) = 0.;
                stepSig = src[8](k);
            }
        }
        if (stepSig < minSig) stepSig = minSig;
        double r = cgauss();
        src[3][k] = par0 + r*stepSig;
        if (src[3](k) < src[6](k))
            src[3](k) = src[6](k);
        if (src[3](k) > src[7](k))
            src[3](k) = src[7](k);
    }
    return cfac;
}

void MC::checkPoint(string outputFileName, MultModelParam &param, bool moveAll, double L, double LMax) {
    ofstream chkpt;
    chkpt.open(outputFileName);
    for (int j=0; j<param.nLens; ++j) {
        for (auto k: freePar[j]) {
            for (size_t l=0; l<6; ++l) chkpt << std::scientific << std::setprecision(7) << param.mixAllModels[l][j].paraList[k] << " ";
            double par0 = param.mixAllModels[4][j].paraList[k];
            double s0 = param.mixAllModels[0][j].paraList[k] + cfac*par0*par0;
            double s1 = param.mixAllModels[1][j].paraList[k] + cfac;
            double s2 = param.mixAllModels[2][j].paraList[k] + cfac*par0;
            double sig = sqrt((1+1e-7)*s0/s1
                        - pow(s2/s1,2));
            chkpt << sig << " ";
            chkpt << endl;
        }
    }
    chkpt << cfac << endl;
    if (moveAll) {
        chkpt << iters << endl;
    } else {
        for (int j=0; j<param.nLens; ++j) {
            for (auto k: freePar[j]) {
                chkpt << iter[j][k] << " ";
            }
            chkpt << endl;
        }
    }
    chkpt << L << endl;
    chkpt << LMax << endl;
    chkpt.close();
}

void MC::load(string fileName, MultModelParam &param, bool moveAll, double &L, double &LMax) {
    ifstream input(fileName.c_str());
    string line, token;
    if (input) {
        for (int j=0; j<param.nLens; ++j) {
            for (auto k: freePar[j]) {
                getline(input, line);
                istringstream ss(line);
                for (size_t l=0; l<6; ++l) {
                    getline(ss, token, ' ');
                    param.mixAllModels[l][j].paraList[k] = stod(token);
                }
            }
        }
        getline(input, line);
        cfac = stod(line);
        if (moveAll) {
            getline(input, line);
            iters = stod(line);
        } else {
            for (int j=0; j<param.nLens; ++j) {
                if (freePar[j].size() > 0) {
                    getline(input, line);
                    istringstream ss(line);
                    for (auto k: freePar[j]) {
                        getline(ss, token, ' ');
                        iter[j][k] = stod(token);
                    }
                }
            }
        }
        getline(input, line);
        L = stod(line);
        LMax = stod(line);
        checkPoint("copy_"+fileName, param, moveAll, L, LMax);
    }
}

void MC::printIterNum(bool moveAll) {
    cout << "Iteration #: " << endl;
    if (moveAll) {
        cout << iters << endl;
    } else {
       for(size_t j=0; j<freePar.size(); ++j) {
           cout << "Lens " << j << endl;
           for (auto k: freePar[j]) {
               cout << iter[j][k] << " ";
           }
           cout << endl;
       }
    }
}

void MC::stretchMove(Model *model, size_t kk) {
    size_t jj;
    do {
        jj = random()*nWalker;
    } while (jj == kk);

    size_t c(0);
    double zz = pow((stepA - 1)*random() + 1,2)/stepA;
    for(int j=0; j<model->param.nLens; ++j) {
        for (auto k: freePar[j]) {
            model->param.mixAllModels[3][j].paraList[k] = par[jj][c] + zz*(par[kk][c]-par[jj][c]);
            if (model->param.mixAllModels[3][j].paraList[k]<bound[0][c] || model->param.mixAllModels[3][j].paraList[k]>bound[1][c])
                return;
            c++;
        }
    }
    iters++;
    double R = objective(model);
    if (R < 0) return;
    double weight(0.5);
    if (R < R0[kk] || log(random()) <= (nFreePar-1)*log(zz) - (R-R0[kk])*weight) {
        R0[kk] = R;
        size_t c(0);
        for(int j=0; j<model->param.nLens; ++j) {
            for (auto k: freePar[j]) {
                par[kk][c] = model->param.mixAllModels[3][j].paraList[k];
                if (R < RMin) bestPar[c] = par[kk][c];
                c++;
            }
        }
        nAccept++;
    }
    if (R < RMin) RMin = R;
}

void MC::load(size_t &loop) {
    ifstream input(chkptFileName.c_str());
    string line, token;
    if (input) {
        for (size_t j=0; j<nWalker; ++j) {
            getline(input, line);
            istringstream ss(line);
            for (size_t m=0; m<nFreePar; ++m) {
                getline(ss, token, ' ');
                par[j][m] = stod(token);
            }
            getline(ss, token, ' ');
            R0[j] = stod(token);
        }
        getline(input, line);
        istringstream ss(line);
        for (size_t m=0; m<nFreePar; ++m) {
            getline(ss, token, ' ');
            bestPar[m] = stod(token);
        }
        getline(ss, token, ' ');
        RMin = stod(token);
        getline(input, line);
        loop = stoi(line);
    }
}



void MC::checkPoint(size_t loop) {
    ofstream chkpt;
    chkpt.open(chkptFileName);
    chkpt << std::scientific << std::setprecision(4);
    for (size_t j=0; j<nWalker; ++j) {
        for (size_t m=0; m<nFreePar; ++m) {
            chkpt << par[j][m] << " ";
        }
        chkpt << R0[j] << endl;
    }
    for (size_t m=0; m<nFreePar; ++m) {
        chkpt << bestPar[m] << " ";
    }
    chkpt << RMin << endl;
    if (runGA) {
        for (size_t k=0; k<4; ++k) {
            for (size_t m=0; m<nFreePar; ++m) {
                chkpt << parSum[k][m] << " ";
                parSum[k][m] = 0.;
            }
            chkpt << endl;
        }
    }
    chkpt << loop << endl;
    chkpt.close();
}

void MC::writeOutput(size_t loop) {
    output << std::scientific << std::setprecision(4);
    for (size_t m=0; m<nFreePar; ++m) {
        output << bestPar[m] << " ";
    }
    output << RMin << " ";
    if (iters>0) output << std::fixed << (double) nAccept/iters << " ";
    output << std::setw(10) << loop << endl;
    //output.close();  //let destructor close it
}

void MC::writeOutput(size_t loop, int thin) {
    double rate = (double) nAccept/iters;
    for (size_t k=0; k<nWalker/thin; ++k) {
        size_t j = random()*nWalker;
        output << std::scientific << std::setprecision(4);
        for (size_t m=0; m<nFreePar; ++m) {
            output << par[j][m] << " ";
        }
        output << R0[j] << " ";
        output << std::fixed << rate << " ";
        output << std::setw(10) << loop << endl;
    }
    nAccept = 0;
    iters = 0;
}

void MC::copyParam(MultModelParam &param) {
    size_t c(0);
    for(int j=0; j<param.nLens; ++j) {
        for (auto k: freePar[j]) {
            param.mixAllModels[5][j].paraList[k] = bestPar[c];
            c++;
        }
    }
}


void MC::evaluate(Model *model, size_t m) {
    size_t c(0);
    for(int j=0; j<model->param.nLens; ++j) {
        for (auto k: freePar[j]) {
            model->param.mixAllModels[3][j].paraList[k] = par[m][c];
            parSum[0][c] += par[m][c]*par[m][c];
            parSum[1][c] += 1.;
            parSum[2][c] += par[m][c];
            c++;
        }
    }
    R0[m] = objective(model);
}

double MC::evaluate(Model *model, double par[]) {
    size_t c(0);
    for(int j=0; j<model->param.nLens; ++j) {
        for (auto k: freePar[j]) {
            model->param.mixAllModels[3][j].paraList[k] = par[c];
            c++;
        }
    }
    return objective(model);
}

void MC::startGA() {
    //selection
    vector<vector<double>> parPrev(par);
    for (size_t m=0; m<nWalker; ++m) {
        size_t jj = pow(random(),1.3) * nWalker;
        for (size_t k=0; k<nFreePar; ++k) par[m][k] = parPrev[index[jj]][k];
    }

    double minSig = 5e-2;
    for (size_t k=0; k<nFreePar; ++k) {
        parSum[3][k] = sqrt(parSum[0][k]/parSum[1][k]-pow(parSum[2][k]/parSum[1][k],2));
        if (std::isnan(parSum[3][k]) || parSum[3][k] < minSig) parSum[3][k] = minSig;
    }

    //crossover
    int first(0);
    size_t i;
    for (size_t m=0; m<nWalker; ++m) {
        double p = random();
        if (p < pCrossOver) {
            ++first;
            if (first & 1) {  //odd
                i = m;
            } else {
                for (size_t k=0; k<random()*nFreePar; ++k) {
                    size_t kk = random()*nFreePar;
                    double t = par[i][kk];
                    par[i][kk] = par[m][kk];
                    par[m][kk] = t;
                }
            }
        } else if (p < pCrossOver + pMutation) { //mutation
            for (size_t k=0; k<nFreePar; ++k) {
                double tpar = par[m][k] + 0.5*cgauss()*parSum[3][k];
                if (tpar >= bound[0][k] && tpar <= bound[1][k]) par[m][k] = tpar;
                //else par[m][k] = bound[0][k] + random()*(bound[1][k] - bound[0][k]);
            }
        }
    }

    /*
    for (size_t t=0; t<pMutation*nFreePar*nWalker; ++t) {
        int p = random()*nFreePar*nWalker;
        size_t j = p / nFreePar;
        size_t k = p % nFreePar;
        double tpar = par[j][k] + cgauss()*parSum[3][k];
        if (tpar >= bound[0][k] && tpar <= bound[1][k]) par[j][k] = tpar;
        else par[j][k] = bound[0][k] + random()*(bound[1][k] - bound[0][k]);
    }*/
}

void MC::elitism() {
    index = sort_indexes(R0);
    if (R0[index[0]] < RMin) {
        RMin = R0[index[0]];
        for (size_t k=0; k<nFreePar; ++k) bestPar[k] = par[index[0]][k];
    } else {
        for (size_t k=0; k<nFreePar; ++k) par[index[0]][k] = bestPar[k];
        R0[index[0]] = RMin;
    }
}

void MC::startNM() {
    //selection
    vector<vector<double>> parPrev(par);
    for (size_t m=0; m<nWalker; ++m) {
        size_t jj = pow(random(),1.4) * nWalker;
        for (size_t k=0; k<nFreePar; ++k) par[m][k] = parPrev[index[jj]][k];
    }
    vector<double> psum1(nFreePar,0.);
    vector<double> psum2(nFreePar,0.);
    vector<double> sig(nFreePar);
    for (size_t m=0; m<nWalker; ++m) {
        for (size_t k=0; k<nFreePar; ++k) {
            psum2[k] += par[m][k]*par[m][k];
            psum1[k] += par[m][k];
        }
    }
    double minSig = 1e-2;
    for (size_t k=0; k<nFreePar; ++k) {
        sig[k] = sqrt(psum2[k]/nWalker - pow(psum1[k]/nWalker,2));
        if (std::isnan(sig[k]) || sig[k] < minSig) sig[k] = minSig;
    }

    //mutation
    for (size_t m=0; m<nWalker; ++m) {
        double p = random();
        if (p < pMutation) {
            for (size_t k=0; k<nFreePar; ++k) {
                par[m][k] = bound[0][k] + random()*(bound[1][k] - bound[0][k]);
            }
        } else {
            for (size_t k=0; k<nFreePar; ++k) {
                double tpar = par[m][k] + cgauss()*sig[k];
                if (tpar >= bound[0][k] && tpar <= bound[1][k]) par[m][k] = tpar;
            }
        }
    }
}

void MC::nelmin(Model *model, size_t m) {
    double reqmin(1e-3);
    double ynewlo(std::numeric_limits<double>::max());
    int konvge(10), kcount(50), icount, numres, ifault;
    double *xmin = new double[nFreePar];

    nelmin(model, &par[m][0], xmin, &ynewlo, reqmin, konvge, kcount, &icount, &numres, &ifault);
    if (ynewlo >= 0 && ynewlo < R0[m]) {
        R0[m] = ynewlo;
        for (size_t k=0; k<nFreePar; ++k) par[m][k] = xmin[k];
    }
}

//****************************************************************************80

void MC::nelmin (Model *model, double start[], double xmin[],
  double *ynewlo, double reqmin, int konvge, int kcount,
  int *icount, int *numres, int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    NELMIN minimizes a function using the Nelder-Mead algorithm.
//
//  Discussion:
//
//    This routine seeks the minimum value of a user-specified function.
//
//    Simplex function minimisation procedure due to Nelder+Mead(1965),
//    as implemented by O'Neill(1971, Appl.Statist. 20, 338-45), with
//    subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
//    25, 97) and Hill(1978, 27, 380-2)
//
//    The function to be minimized must be defined by a function of
//    the form
//
//      function fn ( x, f )
//      double fn
//      double x(*)
//
//    and the name of this subroutine must be declared EXTERNAL in the
//    calling routine and passed as the argument FN.
//
//    This routine does not include a termination test using the
//    fitting of a quadratic surface.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by R ONeill.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    John Nelder, Roger Mead,
//    A simplex method for function minimization,
//    Computer Journal,
//    Volume 7, 1965, pages 308-313.
//
//    R ONeill,
//    Algorithm AS 47:
//    Function Minimization Using a Simplex Procedure,
//    Applied Statistics,
//    Volume 20, Number 3, 1971, pages 338-345.
//
//  Parameters:
//
//    Input, double FN ( double x[] ), the name of the routine which evaluates
//    the function to be minimized.
//
//    Input, int N, the number of variables.
//
//    Input/output, double START[N].  On input, a starting point
//    for the iteration.  On output, this data may have been overwritten.
//
//    Output, double XMIN[N], the coordinates of the point which
//    is estimated to minimize the function.
//
//    Output, double YNEWLO, the minimum value of the function.
//
//    Input, double REQMIN, the terminating limit for the variance
//    of function values.
//
//    Input, double STEP[N], determines the size and shape of the
//    initial simplex.  The relative magnitudes of its elements should reflect
//    the units of the variables.
//
//    Input, int KONVGE, the convergence check is carried out 
//    every KONVGE iterations.
//
//    Input, int KCOUNT, the maximum number of function 
//    evaluations.
//
//    Output, int *ICOUNT, the number of function evaluations 
//    used.
//
//    Output, int *NUMRES, the number of restarts.
//
//    Output, int *IFAULT, error indicator.
//    0, no errors detected.
//    1, REQMIN, N, or KONVGE has an illegal value.
//    2, iteration terminated because KCOUNT was exceeded without convergence.
//
{
  int n = nFreePar;
  double ccoeff = 0.5;
  double del;
  double dn;
  double dnn;
  double ecoeff = 2.0;
  double eps = 0.001;
  int i;
  int ihi;
  int ilo;
  int j;
  int jcount;
  int l;
  int nn;
  double *p;
  double *p2star;
  double *pbar;
  double *pstar;
  double rcoeff = 1.0;
  double rq;
  double x;
  double *y;
  double y2star;
  double ylo;
  double ystar;
  double z;
//
//  Check the input parameters.
//
  if ( reqmin <= 0.0 )
  {
    *ifault = 1;
    return;
  }

  if ( n < 1 )
  {
    *ifault = 1;
    return;
  }

  if ( konvge < 1 )
  {
    *ifault = 1;
    return;
  }

  p = new double[n*(n+1)];
  pstar = new double[n];
  p2star = new double[n];
  pbar = new double[n];
  y = new double[n+1];

  *icount = 0;
  *numres = 0;

  jcount = konvge; 
  dn = ( double ) ( n );
  nn = n + 1;
  dnn = ( double ) ( nn );
  del = 1.0;
  rq = reqmin * dn;
//
//  Initial or restarted loop.
//
  for ( ; ; )
  {
    for ( i = 0; i < n; i++ )
    { 
      p[i+n*n] = start[i];
    }
    //y[n] = fn ( start );
    y[n] = evaluate(model, start);
    *icount = *icount + 1;

    for ( j = 0; j < n; j++ )
    {
      x = start[j];
      start[j] = start[j] + step[j] * del;
      for ( i = 0; i < n; i++ )
      {
        p[i+j*n] = start[i];
      }
      //y[j] = fn ( start );
      y[j] = evaluate(model, start);
      *icount = *icount + 1;
      start[j] = x;
    }
//                    
//  The simplex construction is complete.
//                    
//  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
//  the vertex of the simplex to be replaced.
//                
    ylo = y[0];
    ilo = 0;

    for ( i = 1; i < nn; i++ )
    {
      if ( y[i] < ylo )
      {
        ylo = y[i];
        ilo = i;
      }
    }
//
//  Inner loop.
//
    for ( ; ; )
    {
      if ( kcount <= *icount )
      {
        break;
      }
      *ynewlo = y[0];
      ihi = 0;

      for ( i = 1; i < nn; i++ )
      {
        if ( *ynewlo < y[i] )
        {
          *ynewlo = y[i];
          ihi = i;
        }
      }
//
//  Calculate PBAR, the centroid of the simplex vertices
//  excepting the vertex with Y value YNEWLO.
//
      for ( i = 0; i < n; i++ )
      {
        z = 0.0;
        for ( j = 0; j < nn; j++ )
        { 
          z = z + p[i+j*n];
        }
        z = z - p[i+ihi*n];  
        pbar[i] = z / dn;
      }
//
//  Reflection through the centroid.
//
      for ( i = 0; i < n; i++ )
      {
        pstar[i] = pbar[i] + rcoeff * ( pbar[i] - p[i+ihi*n] );
      }
      //ystar = fn ( pstar );
      ystar = evaluate(model, pstar);
      *icount = *icount + 1;
//
//  Successful reflection, so extension.
//
      if ( ystar < ylo )
      {
        for ( i = 0; i < n; i++ )
        {
          p2star[i] = pbar[i] + ecoeff * ( pstar[i] - pbar[i] );
        }
        //y2star = fn ( p2star );
        y2star = evaluate(model, p2star);
        *icount = *icount + 1;
//
//  Check extension.
//
        if ( ystar < y2star )
        {
          for ( i = 0; i < n; i++ )
          {
            p[i+ihi*n] = pstar[i];
          }
          y[ihi] = ystar;
        }
//
//  Retain extension or contraction.
//
        else
        {
          for ( i = 0; i < n; i++ )
          {
            p[i+ihi*n] = p2star[i];
          }
          y[ihi] = y2star;
        }
      }
//
//  No extension.
//
      else
      {
        l = 0;
        for ( i = 0; i < nn; i++ )
        {
          if ( ystar < y[i] )
          {
            l = l + 1;
          }
        }

        if ( 1 < l )
        {
          for ( i = 0; i < n; i++ )
          {
            p[i+ihi*n] = pstar[i];
          }
          y[ihi] = ystar;
        }
//
//  Contraction on the Y(IHI) side of the centroid.
//
        else if ( l == 0 )
        {
          for ( i = 0; i < n; i++ )
          {
            p2star[i] = pbar[i] + ccoeff * ( p[i+ihi*n] - pbar[i] );
          }
          //y2star = fn ( p2star );
          y2star = evaluate(model, p2star);
          *icount = *icount + 1;
//
//  Contract the whole simplex.
//
          if ( y[ihi] < y2star )
          {
            for ( j = 0; j < nn; j++ )
            {
              for ( i = 0; i < n; i++ )
              {
                p[i+j*n] = ( p[i+j*n] + p[i+ilo*n] ) * 0.5;
                xmin[i] = p[i+j*n];
              }
              //y[j] = fn ( xmin );
              y[j] = evaluate(model, xmin);
              *icount = *icount + 1;
            }
            ylo = y[0];
            ilo = 0;

            for ( i = 1; i < nn; i++ )
            {
              if ( y[i] < ylo )
              {
                ylo = y[i];
                ilo = i;
              }
            }
            continue;
          }
//
//  Retain contraction.
//
          else
          {
            for ( i = 0; i < n; i++ )
            {
              p[i+ihi*n] = p2star[i];
            }
            y[ihi] = y2star;
          }
        }
//
//  Contraction on the reflection side of the centroid.
//
        else if ( l == 1 )
        {
          for ( i = 0; i < n; i++ )
          {
            p2star[i] = pbar[i] + ccoeff * ( pstar[i] - pbar[i] );
          }
          //y2star = fn ( p2star );
          y2star = evaluate(model, p2star);
          *icount = *icount + 1;
//
//  Retain reflection?
//
          if ( y2star <= ystar )
          {
            for ( i = 0; i < n; i++ )
            {
              p[i+ihi*n] = p2star[i];
            }
            y[ihi] = y2star;
          }
          else
          {
            for ( i = 0; i < n; i++ )
            {
              p[i+ihi*n] = pstar[i];
            }
            y[ihi] = ystar;
          }
        }
      }
//
//  Check if YLO improved.
//
      if ( y[ihi] < ylo )
      {
        ylo = y[ihi];
        ilo = ihi;
      }
      jcount = jcount - 1;

      if ( 0 < jcount )
      {
        continue;
      }
//
//  Check to see if minimum reached.
//
      if ( *icount <= kcount )
      {
        jcount = konvge;

        z = 0.0;
        for ( i = 0; i < nn; i++ )
        {
          z = z + y[i];
        }
        x = z / dnn;

        z = 0.0;
        for ( i = 0; i < nn; i++ )
        {
          z = z + pow ( y[i] - x, 2 );
        }

        if ( z <= rq )
        {
          break;
        }
      }
    }
//
//  Factorial tests to check that YNEWLO is a local minimum.
//
    for ( i = 0; i < n; i++ )
    {
      xmin[i] = p[i+ilo*n];
    }
    *ynewlo = y[ilo];

    if ( kcount < *icount )
    {
      *ifault = 2;
      break;
    }

    *ifault = 0;

    for ( i = 0; i < n; i++ )
    {
      del = step[i] * eps;
      xmin[i] = xmin[i] + del;
      //z = fn ( xmin );
      z = evaluate(model, xmin);
      *icount = *icount + 1;
      if ( z < *ynewlo )
      {
        *ifault = 2;
        break;
      }
      xmin[i] = xmin[i] - del - del;
      //z = fn ( xmin );
      z = evaluate(model, xmin);
      *icount = *icount + 1;
      if ( z < *ynewlo )
      {
        *ifault = 2;
        break;
      }
      xmin[i] = xmin[i] + del;
    }

    if ( *ifault == 0 )
    {
      break;
    }
//
//  Restart the procedure.
//
    for ( i = 0; i < n; i++ )
    {
      start[i] = xmin[i];
    }
    del = eps;
    *numres = *numres + 1;
  }
  delete [] p;
  delete [] pstar;
  delete [] p2star;
  delete [] pbar;
  delete [] y;

  return;
}


//
//  biexponential.h
//  
//
//  Created by для интернета on 12.04.18.
//
//

#ifndef biexponential_h
#define biexponential_h

#include <octave-4.2.2/octave/oct.h>
#include <octave-4.2.2/octave/parse.h>

class Sweep {
    
private:
    
    Matrix sweepA, sweepB;
    double* p;
    double* q;
    int n;

public:
    
    Sweep(Matrix sweepA, Matrix sweepB, const int n);
    double* findSolution();
};

class Exponential {
    
private:
    
    Matrix f;
    double startFrom, endAt;
    int n;
    
    double h, plotHx;
    
    double* b;
    double* e;
    double* d;
    double* s;
    double* c;
    double* tau;
    
    void initializeSpline();
    void initB();
    void initS();
    void initC();
    void initE();
    void initD();
    void initTau();
    
public:
    
    Exponential(Matrix uNth, const double startFrom, const int endAt, const int n, const double plotHx);
    double valueAt(double x);
    double getTau(int i);
    double getF(int i);
};

int findIndex(double xy, double startFrom, double h);
int findNumberOfBlocks(int indexX, int indexY);

#endif /* biexponential_h */

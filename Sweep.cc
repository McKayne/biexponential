//
//  Sweep.cpp
//  
//
//  Created by для интернета on 12.04.18.
//
//

#include <stdio.h>

#include "biexponential.h"

Sweep::Sweep(Matrix sweepA, Matrix sweepB, const int n) {
    this->sweepA = sweepA;
    this->sweepB = sweepB;
    this->n = n;
    
    p = new double[n];
    q = new double[n];
    
    p[0] = -sweepA(0, 1) / sweepA(0, 0);
    q[0] = sweepB(0) / sweepA(0, 0);
    
    for (int i = 1; i < n - 1; i++) {
        p[i] = -sweepA(i, i + 1) / (sweepA(i, i) + sweepA(i, i - 1) * p[i - 1]);
        q[i] = (sweepB(i) - sweepA(i, i - 1) * q[i - 1]) / (sweepA(i, i) + sweepA(i, i - 1) * p[i - 1]);
        //cout << q[i] << endl;
    }
    
    p[n - 1] = 0.0;
    q[n - 1] = (sweepB(n - 1) - sweepA(n - 1, n - 2) * q[n - 2]) / (sweepA(n - 1, n - 1) + sweepA(n - 1, n - 2) * p[n - 2]);
    
    //cout << q[n - 1] << sweepA(n - 1, n - 1) << endl;
}

double* Sweep::findSolution() {
    double* x = new double[n];
    
    x[n - 1] = q[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        x[i] = p[i] * x[i + 1] + q[i];
    }
    
    return x;
}

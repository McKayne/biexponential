//Compile with mkoctfile Exponential.cc

#include <octave-4.2.2/octave/oct.h>
#include <octave-4.2.2/octave/parse.h>

#include <iostream>

#include "biexponential.h"

using namespace std;

Exponential::Exponential(Matrix uNth, const double startFrom, const int endAt, const int n, const double plotHx) {
    f = uNth;
    this->startFrom = startFrom;
    this->endAt = endAt;
    this->n = n;
    
    this->plotHx = plotHx;
    
    h = (endAt - startFrom) / n;
    cout << h << endl;
    
    initializeSpline();
}

void Exponential::initializeSpline() {
    initB();
    initS();
    initC();
    initE();
    initD();
    initTau();
}

void Exponential::initB() {
    b = new double[n + 1];
    
    b[0] = (f(1) - f(0)) / h - 0.0;
    
    for (int i = 1; i < n; i++) {
        b[i] = (f(i + 1) - f(i)) / h - (f(i) - f(i - 1)) / h;
        //cout << b[i] << endl;
        //cout << f(i) << endl;
    }
    
    b[n] = 0.0 - (f(n) - f(n - 1)) / h;
}

void Exponential::initS() {
    s = new double[n];
    
    for (int i = 0; i < n; i++) {
        s[i] = sinh(1.0 * h);
    }
}

void Exponential::initC() {
    c = new double[n];
    
    for (int i = 0; i < n; i++) {
        c[i] = cosh(1.0 * h);
    }
}

void Exponential::initE() {
    e = new double[n];
    
    for (int i = 0; i < n; i++) {
        e[i] = (1.0 / h - 1.0 / s[i]) / 1.0;
    }
}

void Exponential::initD() {
    d = new double[n];
    
    for (int i = 0; i < n; i++) {
        d[i] = (1.0 * c[i] / s[i] - 1 / h) / 1.0;
    }
}

void Exponential::initTau() {
    Matrix sweepA(n + 1, n + 1);
    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= n; j++) {
            sweepA(i, j) = 0.0;
        }
    }
    Matrix sweepB(n + 1, 1);
    
    sweepA(0, 0) = d[0];
    sweepA(0, 1) = e[0];
    sweepB(0) = b[0];
    
    for (int i = 1; i < n; i++) {
        sweepA(i, i - 1) = e[i - 1];
        sweepA(i, i) = d[i - 1] + d[i];
        sweepA(i, i + 1) = e[i];
        sweepB(i) = b[i];
    }
    
    sweepA(n, n - 1) = e[n - 1];
    sweepA(n, n) = d[n - 1];
    cout << d[n - 1] << endl;
    sweepB(n) = b[n];
    
    Sweep sw(sweepA, sweepB, n + 1);
    tau = sw.findSolution();
    
    cout << tau[0] << endl;
    cout << tau[1] << endl;
    cout << tau[2] << endl;
    cout << tau[3] << endl;
    cout << tau[4] << endl;
}

double Exponential::valueAt(double x) {
    int i = findIndex(x, startFrom, h);
    cout << x << "\t" << i << endl;
    
    double prevX = startFrom + i * h;
    double nextX = startFrom + (i + 1) * h;
    cout << ">> " << prevX << endl;
    
    return (tau[i] * sinh(nextX - x) + tau[i + 1] * sinh(x - prevX)) / s[i] + (f(i) - tau[i]) * (nextX - x) / h + (f(i + 1) - tau[i + 1]) * (x - prevX) / h;
}

double Exponential::getTau(int i) {
    return tau[i];
}

double Exponential::getF(int i) {
    return f(i);
}

Matrix getOctaveMatrix(Exponential e, const double startFrom, const int plotX, const double plotHx) {
	Matrix v(1, plotX + 1);
	for (int j = 0; j <= plotX; j++) {
		double x = startFrom + j * plotHx;
		v(0, j) = e.valueAt(x);
	}
	return v;
}

DEFUN_DLD(Exponential, args, nargout, "") {
	Matrix uNth = args(0).matrix_value();
	const double startFrom = args(1).double_value();
	const int endAt = args(2).int_value();
	const int n = args(3).int_value();
	const int plotX = args(4).int_value();

	const double range = endAt - startFrom;
	const double plotHx = range / plotX;

	Exponential e(uNth, startFrom, endAt, n, plotHx);
	Matrix v = getOctaveMatrix(e, startFrom, plotX, plotHx);
	return octave_value(v);
}

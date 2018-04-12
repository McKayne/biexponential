//Compile with mkoctfile Biexponential.cc

#include <octave-4.2.2/octave/oct.h>
#include <octave-4.2.2/octave/parse.h>

#include "biexponential.h"

using namespace std;

/*Matrix getOctaveMatrix(Exponential e, const double xStartFrom, const int n, const int k, const int plotX, const int plotY, const double plotHx) {
	Matrix v(plotY, plotX);
	for (int i = 0; i < plotY; i++) {
		for (int j = 0; j < plotX; j++) {
            double x = xStartFrom + j * plotHx;
			v(i, j) = e.valueAt(x);
		}
	}
	return v;
}*/

Matrix getOctaveMatrix(Exponential** ex, Exponential** ey, const double xStartFrom, const double yStartFrom, const int n, const int k, const int plotX, const int plotY, const double plotHx, const double plotHy) {
    Matrix v(plotY, plotX, 0.0);
    
    for (int i = 0; i <= k; i++) {
        cout << "+ " << plotY << endl;
        for (int j = 0; j <= plotX; j++) {
            double x = xStartFrom + j * plotHx;
            v(i * (plotY / k), j) = ex[i]->valueAt(x);
        }
    }
    
    for (int i = 0; i <= plotY; i++) {
        cout << "+ " << plotY << endl;
        for (int j = 0; j <= n; j++) {
            double y = yStartFrom + i * plotHy;
            v(i, j * (plotX / n)) = ey[j]->valueAt(y);
        }
    }
    
    return v;
}

DEFUN_DLD(Biexponential, args, nargout, "") {
	Matrix u = args(0).matrix_value();
    
    const double xStartFrom = args(1).double_value();
    const double xEndAt = args(2).double_value();
    const double yStartFrom = args(3).double_value();
    const double yEndAt = args(4).double_value();
    
	const int n = args(5).int_value();
	const int k = args(6).int_value();
    
	const int plotX = args(7).int_value();
	const int plotY = args(8).int_value();
    
    const double xRange = xEndAt - xStartFrom;
    const double yRange = yEndAt - yStartFrom;
    const double plotHx = xRange / plotX;
    const double plotHy = yRange / plotY;

    Exponential** ex = new Exponential*[k + 1];
    Exponential** ey = new Exponential*[n + 1];
    for (int i = 0; i <= k; i++) {
        Matrix uNth = (Matrix) u.row(i);
        ex[i] = new Exponential(uNth, xStartFrom, xEndAt, n, plotX);
    }
    for (int i = 0; i <= n; i++) {
        Matrix uNth = (Matrix) u.column(i);
        ey[i] = new Exponential(uNth, yStartFrom, yEndAt, k, plotY);
    }
    
    //Exponential e(uNth, xStartFrom, xEndAt, n, plotX);
	Matrix v = getOctaveMatrix(ex, ey, xStartFrom, yStartFrom, n, k, plotX, plotY, plotHx, plotHy);
	return octave_value(v);
}

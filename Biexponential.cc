//Compile with mkoctfile Biexponential.cc

#include <octave-4.2.2/octave/oct.h>
#include <octave-4.2.2/octave/parse.h>

#include "biexponential.h"

using namespace std;

int findIndex(double xy, double startFrom, double h) {
    for (int i = 1;; i++) {
        if (xy - startFrom <= i * h) {
            return i - 1;
        }
    }
    return -1;
}

double findAngle(double x, double y) {
    if (y == 0) {
        return 0;
    } else if (x == 0) {
        return 90;
    } else if (x == y) {
        return 45;
    } else {
        return atan(y / x) * (180 / M_PI);
    }
}

double findStep(double phi, double hx, double hy) {
    //cout << "c : " << cos(30 * M_PI / 180) << endl;
    //TODO improve
    double eta = hx / cos(phi * M_PI / 180);
    
    if (phi > 45) {
        return hy / sin(phi * M_PI / 180);
    }
    
    return eta;
}

double xi(Exponential* ex, Exponential* ey, const int indexX, const int indexY, const double xStartFrom, const double yStartFrom, const double hx, const double hy, double phi) {
    if (phi <= 45.0) {
        double eta = phi * hy / 45.0;
        
        return ey->valueAt(eta + indexY * hy + yStartFrom);
    } else {
        
        double eta = hx - (phi - 45) * hx / 45.0;
        
        return ex->valueAt(eta + indexX * hx + xStartFrom);
    }
}

double startTau(Exponential** ex, Exponential** ey, const int n, const int k, const double hx, const double hy, int indexX, int indexY, double phi, double xStartFrom, double yStartFrom) {
    
    //return 0.123;
    Exponential* ex1 = ex[indexY];
    Exponential* ex2 = ex[indexY + 1];
    Exponential* ey1 = ey[indexX];
    Exponential* ey2 = ey[indexX + 1];
    
    
    if (phi == 0) {
        return ex1->getTau(indexX);
    } else if (phi == 90) {
        return ey1->getTau(indexY);
    } else if (phi == 45) {
        int blocks = findNumberOfBlocks(n, k, indexX, indexY);
        //double hPhi = sqrt(8.0);
        double hPhi = sqrt(pow(hx, 2) + pow(hy, 2));
        
        int currentIndexX = indexX;
        int currentIndexY = indexY;
        int diagNth = 0;
        while (true) {
            if (currentIndexX > 0 && currentIndexY > 0) {
                currentIndexX--;
                currentIndexY--;
                diagNth++;
                //n++;
            } else {
                break;
            }
        }
        Matrix diagSpline(1, blocks + 1);
        diagSpline(0) = ex[currentIndexY]->getF(currentIndexX);
        for (int i = 1; i <= blocks; i++) {
            diagSpline(i) = ex[currentIndexY + 1]->getF(currentIndexX + 1);
            
            currentIndexX++;
            currentIndexY++;
        }
        Exponential e(diagSpline, 0, blocks * hPhi, blocks, -1);
        
        return e.getTau(diagNth);
    } else {
        int nth;
        cout << endl << "Now processing " << indexX << "; " << indexY << "; " << phi << endl << endl;
        Exponential2 e = buildDiagSpline(indexX, indexY, phi, ex, ey, xStartFrom, yStartFrom, hx, hy, &nth);
        cout << endl << "Nth = " << nth << endl << endl;
    
        return e.getTau(nth);
    }
    
    /*Matrix tau(1, 3);
    tau(0) = ex1->getTau(indexX);
    tau(1) = e.getTau(diagNth);
    tau(2) = ey1->getTau(indexY);
    
    Exponential eTau(tau, 0, 90, 2, -1);

    if (phi == 0) {
        return ex1->getTau(indexX);
    } else if (phi == 90) {
        return ey1->getTau(indexY);
    } else if (phi == 45) {
        return e.getTau(diagNth);
    }
    
    return eTau.valueAt(phi);*/
    
    
    
    
}

double endTau(Exponential** ex, Exponential**ey, const int n, const int k, const double hx, const double hy, int indexX, int indexY, double phi) {
    
    Exponential* ex1 = ex[indexY];
    Exponential* ex2 = ex[indexY + 1];
    Exponential* ey1 = ey[indexX];
    Exponential* ey2 = ey[indexX + 1];
    
    if (phi == 0) {
        return ex1->getTau(indexX + 1);
    } else if (phi == 45) {
        int blocks = findNumberOfBlocks(n, k, indexX, indexY);
        //double hPhi = sqrt(8.0);
        double hPhi = sqrt(pow(hx, 2) + pow(hy, 2));
        
        int currentIndexX = indexX;
        int currentIndexY = indexY;
        int diagNth = 0;
        while (true) {
            if (currentIndexX > 0 && currentIndexY > 0) {
                currentIndexX--;
                currentIndexY--;
                diagNth++;
                //n++;
            } else {
                break;
            }
        }
        Matrix diagSpline(1, blocks + 1);
        diagSpline(0) = ex[currentIndexY]->getF(currentIndexX);
        for (int i = 1; i <= blocks; i++) {
            diagSpline(i) = ex[currentIndexY + 1]->getF(currentIndexX + 1);
            
            currentIndexX++;
            currentIndexY++;
        }
        Exponential e(diagSpline, 0, blocks * hPhi, blocks, -1);
        
        return e.getTau(diagNth + 1);
    } else if (phi == 90) {
        return ey1->getTau(indexY + 1);
    } else {
        int nth;
        cout << endl << "Now processing " << indexX << "; " << indexY << "; " << phi << endl << endl;
        Exponential2 e = buildDiagSpline(indexX, indexY, phi, ex, ey, -5, -5, hx, hy, &nth);
        cout << endl << "Nth = " << nth << endl << endl;
        
        return e.getTau(nth + 1);
    }
}

int findNumberOfBlocks(const int n, const int k, int indexX, int indexY) {
    int num = 1;
    
    int currentIndexX = indexX;
    int currentIndexY = indexY;
    while (true) {
        if (currentIndexX < n - 1 && currentIndexY < k - 1) {
            currentIndexX++;
            currentIndexY++;
            num++;
        } else {
            break;
        }
    }
    
    currentIndexX = indexX;
    currentIndexY = indexY;
    while (true) {
        if (currentIndexX > 0 && currentIndexY > 0) {
            currentIndexX--;
            currentIndexY--;
            num++;
        } else {
            break;
        }
    }
    
    return num;
}

double valueAt(Exponential** ex, Exponential** ey, double x, double y, const double xStartFrom, const double yStartFrom, const int n, const int k, const double hx, const double hy) {
    int indexX = findIndex(x, xStartFrom, hx), indexY = findIndex(y, yStartFrom, hy);
    
    Exponential* ex1 = ex[indexY];
    Exponential* ex2 = ex[indexY + 1];
    Exponential* ey1 = ey[indexX];
    Exponential* ey2 = ey[indexX + 1];
    
    double currentX = x - (xStartFrom + indexX * hx), currentY = y - (yStartFrom + indexY * hy);
    //cout << "CX = " << currentX << "; " << currentY << endl;

    double r = sqrt(pow(currentX, 2) + pow(currentY, 2));
    //cout << "r = " << r << endl;
    
    double phi = findAngle(currentX, currentY);
    double hPhi = findStep(phi, hx, hy);
    

    double a = (startTau(ex, ey, n, k, hx, hy, indexX, indexY, phi, xStartFrom, yStartFrom) * sinh(hPhi - r) + endTau(ex, ey, n, k, hx, hy, indexX, indexY, phi) * sinh(r)) / sinh(hPhi);
    
    double b = (ex1->getF(indexX) - startTau(ex, ey, n, k, hx, hy, indexX, indexY, phi, xStartFrom, yStartFrom)) * (hPhi - r) / hPhi;
    
    double c = (xi(ex2, ey2, indexX, indexY, xStartFrom, yStartFrom, hx, hy, phi) - endTau(ex, ey, n, k, hx, hy, indexX, indexY, phi)) * r / hPhi;
    
    return a + b + c;
    return 0;
}

Matrix getOctaveMatrix(Exponential** ex, Exponential** ey, const double xStartFrom, const double yStartFrom, const int n, const int k, const int plotX, const int plotY, const double hx, const double hy, const double plotHx, const double plotHy) {
    Matrix v(plotY + 1, plotX + 1, 0.0);
    
    for (int i = 0; i <= plotY; i++) {
        for (int j = 0; j <= plotX; j++) {
            double x = xStartFrom + j * plotHx;
            double y = yStartFrom + i * plotHy;
            
            v(i, j) = valueAt(ex, ey, x, y, xStartFrom, yStartFrom, n, k, hx, hy);
            if (isnan(v(i, j))) {
                cout << "NaN at x = " << (double) x << ", y = " << y << endl;
            }
            //cout << "++++ " << v(i, j) << endl;
        }
    }
    for (int i = 0; i < plotY; i++) {
        for (int j = 20; j < plotX; j++) {
            double x = xStartFrom + j * plotHx;
            double y = yStartFrom + i * plotHy;
            
            cout << "(" << x << "; " << y << ")" << endl;
            cout << "Plot " << plotHx << endl;
        }
    }
    
    /*for (int i = 0; i <= k; i++) {
        //cout << "+ " << plotY << endl;
        for (int j = 0; j <= plotX; j++) {
            double x = xStartFrom + j * plotHx;
            double y = yStartFrom + i * plotHy;
            //cout << -5 + i * (plotY / k) << endl;
            v(i * (plotY / k), j) = valueAt(ex, ey, x, yStartFrom + i * hy, xStartFrom, yStartFrom, hx, hy);
        }
    }
    
    for (int i = 0; i <= plotY; i++) {
        cout << "+ " << plotY << endl;
        for (int j = 0; j <= n; j++) {
            double x = xStartFrom + j * plotHx;
            double y = yStartFrom + i * plotHy;
            v(i, j * (plotX / n)) = valueAt(ex, ey, xStartFrom + j * hx, y, xStartFrom, yStartFrom, hx, hy);
        }
    }*/

    return v;
}

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

/*Matrix getOctaveMatrix(Exponential** ex, Exponential** ey, const double xStartFrom, const double yStartFrom, const int n, const int k, const int plotX, const int plotY, const double plotHx, const double plotHy) {
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
}*/

void toPrevBlock(double xStartFrom, double yStartFrom, int indexX, int indexY, double phi, double hx, double hy, int* nth, double* pointX, double* pointY) {
    double xStart = 0, yStart = 0;
    //double pointX = 0, pointY = 0;

    *nth = 0;
    
    xStart = 0;
    yStart = 0;
    
    *pointX = xStartFrom;
    *pointY = yStartFrom;
    while (indexX > 0 && indexY > 0) {
    //for (int i = 0; i < 7; i++) {
        double eta = phi * (hy - xStart) / 45.0 + yStart;
        
        //cout << "Phi = " << phi << "; eta = " << eta << endl;
        
        if (eta < hy) {
            //cout << "Goto x - 1" << endl;
            indexX--;
            
            xStart = 0;
            yStart = eta;//abs(phi - 45) * hx / 45.0;
        } else if (eta > hy) {
            //cout << "Goto y - 1" << endl;
            indexY--;
            
            xStart = hx - abs(phi - 45) * hx / 45.0;
            yStart = 0;
        } else {
            //cout << "Goto x - 1 and y - 1" << endl;
            indexX--;
            indexY--;
            
            xStart = 0;
            yStart = 0;
        }
        
        //pointX -= xStart;
        //pointY -= yStart;
        
        *pointX = 5 - ((5 - indexX) * hx + xStart) + 0;
        *pointY = 5 - ((5 - indexY) * hy + yStart) + 0;
        //cout << "New point = " << pointX << "; " << pointY << endl;
        //cout << "New center = " << xStart << "; " << yStart << endl;
        //cout << "New index = " << indexX << "; " << indexY << endl;
        
        *nth = *nth + 1;
    }
    
    //cout << "Start point = " << pointX << "; " << pointY << endl;
}

double* toNextBlock(Exponential** ex, Exponential** ey, double xStartFrom, double yStartFrom, int indexX, int indexY, double phi, double hx, double hy, double xStart, double yStart, int* totalValues, double* h, bool debug) {
    double pointX = xStartFrom, pointY = yStartFrom;
    
    double* values = NULL;
    if (*totalValues > 0) {
        values = new double[*totalValues];
    }
    
    *totalValues = 1;
    
    pointX = indexX * hx + xStart + xStartFrom, pointY = indexY * hy + yStart + yStartFrom;
    if (values != NULL) {
        values[0] = ey[indexX]->valueAt(yStart);
    }
    if (debug) {
        cout << "New center = " << xStart << "; " << yStart << endl;
        cout << "New index = " << indexX << "; " << indexY << endl;
        cout << "New point = " << pointX << "; " << pointY << ", diff " << 0 << endl;
        cout << endl;
    }
    while (indexX < 5 && indexY < 5) {
        double eta = phi * (hy - xStart) / 45.0 + yStart;
    
        //cout << "Phi = " << phi << "; eta = " << eta << endl;
    
        if (eta < hy) {
            //cout << "Goto x + 1" << endl;
            indexX++;
        
            xStart = 0;
            yStart = eta;
            
            if (values != NULL) {
                values[*totalValues] = ey[indexX]->valueAt(yStart);
            }
            
            if (debug) {
                cout << "New center = " << xStart << "; " << yStart << endl;
                cout << "New index = " << indexX << "; " << indexY << endl;
            }
        } else if (eta > hy) {
            //cout << "Goto y + 1" << endl;
            indexY++;
        
            xStart = hx - abs(phi - 45) * hx / 45.0;
            yStart = 0;
            
            if (values != NULL) {
                values[*totalValues] = ex[indexY]->valueAt(xStart);
            }
            
            if (debug) {
                cout << "New center = " << xStart << "; " << yStart << endl;
                cout << "New index = " << indexX << "; " << indexY << endl;
            }
        } else {
            cout << "Goto x + 1 and y + 1" << endl;
            indexX++;
            indexY++;
        
            xStart = 0;
            yStart = 0;
            
            if (values != NULL) {
                values[*totalValues] = ey[indexX]->valueAt(0);
            }
                
            if (debug) {
                cout << "New center = " << xStart << "; " << yStart << endl;
                cout << "New index = " << indexX << "; " << indexY << endl;
            }
        }
        
        double oldX = pointX, oldY = pointY;
        pointX = indexX * hx + xStart + xStartFrom, pointY = indexY * hy + yStart + yStartFrom;
        double diff = sqrt(pow(pointX - oldX, 2) + pow(pointY - oldY, 2));
        
        if (debug) {
            cout << "New point = " << pointX << "; " << pointY << ", diff " << diff << endl;
            cout << endl;
        }
        
        if (values != NULL) {
            h[*totalValues - 1] = diff;
        }
        
        *totalValues = *totalValues + 1;
        //
        //
    }
    
    cout << "End point = " << pointX << "; " << pointY << ", " << *totalValues << " total" << endl;
    
    return values;
}

Exponential2 buildDiagSpline( int indexX,  int indexY, const double phi, Exponential** ex, Exponential** ey, const double xStartFrom, const double yStartFrom, const double hx, const double hy, int* nth) {
    double pointX, pointY;
    
    int totalValues = 0;
    
    toPrevBlock(xStartFrom, yStartFrom, indexX, indexY, phi, hx, hy, nth, &pointX, &pointY);
    //cout << "Nth = " << nth << endl;
    cout << "Start point = " << pointX << "; " << pointY << endl;
    
    int newX = findIndex(pointX, xStartFrom, hx), newY = findIndex(pointY, yStartFrom, hy);
    double startX = pointX - (newX * hx + xStartFrom), startY = pointY - (newY * hy + yStartFrom);
    
    //TODO dirty hack
    /*if (startX == hx) {
        newX++;
        startX = 0;
    }
    if (startY == hy) {
        newY++;
        startY = 0;
    }*/
    
    cout << "Index = " << newX << "; " << newY << endl;
    cout << "Start = " << startX << "; " << startY << endl;
    
    toNextBlock(ex, ey, xStartFrom, yStartFrom, newX, newY, phi, hx, hy, startX, startY, &totalValues, NULL, true);
    cout << "Total " << totalValues << endl;
    
    double* h = new double[totalValues];
    double* values = toNextBlock(ex, ey, xStartFrom, yStartFrom, newX, newY, phi, hx, hy, startX, startY, &totalValues, h, false);
    
    for (int i = 0; i < totalValues; i++) {
        cout << ">> " << values[i] << ", " << h[i] << endl;
    }
    //cout << values[0] << ", " << values[3] << ", " << endl;
    cout << endl;
    
    return Exponential2(values, h, 0, totalValues - 1);
}

/*int newIndexAt(double val) {
    for (int i = 1; i < 5; i++) {
        if (i * 2 > val + 5) {
            return i - 1;
        }
    }
}*/

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
    
    const double hx = xRange / n;
    const double hy = yRange / k;
    
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
    Matrix v = getOctaveMatrix(ex, ey, xStartFrom, yStartFrom, n, k, plotX, plotY, hx, hy, plotHx, plotHy);

    /*cout << "A: " << findAngle(1, 0) << "; " << findStep(findAngle(1, 0), hx, hy) << endl;
    cout << "A: " << findAngle(1, 0.2) << "; " << findStep(findAngle(1, 0.2), hx, hy) << endl;
    cout << "A: " << findAngle(1, 0.4) << "; " << findStep(findAngle(1, 0.4), hx, hy) << endl;
    cout << "A: " << findAngle(1, 0.6) << "; " << findStep(findAngle(1, 0.6), hx, hy) << endl;
    cout << "A: " << findAngle(1, 0.8) << "; " << findStep(findAngle(1, 0.8), hx, hy) << endl;
    cout << "A: " << findAngle(1, 1) << "; " << findStep(findAngle(1, 1), hx, hy) << endl;
    cout << "A: " << findAngle(1, 1.2) << "; " << findStep(findAngle(1, 1.2), hx, hy) << endl;
    cout << "A: " << findAngle(1, 1.5) << "; " << findStep(findAngle(1, 1.5), hx, hy) << endl;
    cout << "A: " << findAngle(0.2, 1.9) << "; " << findStep(findAngle(0.2, 1.9), hx, hy) << endl;*/
    
    /*valueAt(ex, ey, -5, 1, xStartFrom, yStartFrom, hx, hy);
    valueAt(ex, ey, -4, 1, xStartFrom, yStartFrom, hx, hy);
    valueAt(ex, ey, -3, 1, xStartFrom, yStartFrom, hx, hy);
    valueAt(ex, ey, -2, 1, xStartFrom, yStartFrom, hx, hy);
    valueAt(ex, ey, -1, 1, xStartFrom, yStartFrom, hx, hy);
    valueAt(ex, ey, 0, 1, xStartFrom, yStartFrom, hx, hy);
    valueAt(ex, ey, 1, 1, xStartFrom, yStartFrom, hx, hy);
    valueAt(ex, ey, 2, 1, xStartFrom, yStartFrom, hx, hy);
    valueAt(ex, ey, 3, 1, xStartFrom, yStartFrom, hx, hy);*/
    
    /*cout << "Xi: " << 0.0 * 2 / 45.0 << endl;
    cout << "Xi: " << 10.0 * 2 / 45.0 << endl;
    cout << "Xi: " << 30.0 * 2 / 45.0 << endl;
    cout << "Xi: " << 40.0 * 2 / 45.0 << endl;
    cout << "Xi: " << 45.0 * 2 / 45.0 << endl;*/
    
    /*for (int i = 0; i <= 4; i++) {
        cout << "index: " << i << "; " << 4 << ", blocks: " << findNumberOfBlocks(i, 4) << endl;
    }
    cout << endl;
    
    for (int i = 0; i <= 4; i++) {
        cout << "index: " << i << "; " << 3 << ", blocks: " << findNumberOfBlocks(i, 3) << endl;
    }
    cout << endl;
    
    for (int i = 0; i <= 4; i++) {
        cout << "index: " << i << "; " << 2 << ", blocks: " << findNumberOfBlocks(i, 2) << endl;
    }
    cout << endl;
    
    for (int i = 0; i <= 4; i++) {
        cout << "index: " << i << "; " << 1 << ", blocks: " << findNumberOfBlocks(i, 1) << endl;
    }
    cout << endl;
    
    for (int i = 0; i <= 4; i++) {
        cout << "index: " << i << "; " << 0 << ", blocks: " << findNumberOfBlocks(i, 0) << endl;
    }
    cout << endl;*/
    
    
    
    //Now processing 0; 0; 36.8699
    cout << "Start debug" << endl << endl;
    
    int nth;
    Exponential2 e = buildDiagSpline(0, 0, 36.8699, ex, ey, xStartFrom, yStartFrom, hx, hy, &nth);
    cout << endl << "Nth = " << nth << endl << endl;
    
    cout << "Test " << valueAt(ex, ey, 5, -1, xStartFrom, yStartFrom, n, k, hx, hy) << endl << endl;
    
    //toPrevBlock(xStartFrom, yStartFrom, 5, 5, 60, hx, hy, 0, 0);
    //toNextBlock(xStartFrom, yStartFrom, 0, 0, 60, hx, hy, 0, 0);
    //cout << endl;
    //toPrevBlock(xStartFrom, yStartFrom, 4, 4, 60, hx, hy, 0, 0);
    //toNextBlock(xStartFrom, yStartFrom, 0, 0, 45, hx, hy, 0, 0);
    
    return octave_value(v);
}

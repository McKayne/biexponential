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

double l(double value1, double value2, double value3, double phi) {
    if (phi == 0) {
        return value1;
    } else if (phi == 45) {
        return value2;
    } else if (phi == 90) {
        return value3;
    } else if (phi < 45) {
        return value1 + (value2 - value1) / 45 * phi;
    } else {
        return value2 + (value3 - value2) / 45 * (phi - 45);
    }
}

double l2(double value1, double value2, double phi) {
    return value1 + (value2 - value1) / 90 * phi;
}

double xi(Exponential* ex, Exponential* ey, const int indexX, const int indexY, const double xStartFrom, const double yStartFrom, double hx, double hy, double phi) {
    if (phi <= 45.0) {
        double eta = phi * hy / 45.0;
        
        return ey->valueAt(eta + indexY * 2 + yStartFrom);
    } else {
        
        double eta = 2 - (phi - 45) * hx / 45.0;
        
        return ex->valueAt(eta + indexX * 2 + xStartFrom);
    }
}

double startTau(Exponential** ex, Exponential** ey, int indexX, int indexY, double phi) {
    Exponential* ex1 = ex[indexY];
    Exponential* ex2 = ex[indexY + 1];
    Exponential* ey1 = ey[indexX];
    Exponential* ey2 = ey[indexX + 1];
    
    
    int blocks = findNumberOfBlocks(indexX, indexY);
    double hPhi = sqrt(8.0);
    
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
    
    Matrix tau(1, 3);
    tau(0) = ex1->getTau(indexX);
    tau(1) = e.getTau(diagNth);
    tau(2) = ey1->getTau(indexY);
    
    Exponential eTau(tau, 0, 90, 2, -1);
    
    return eTau.valueAt(phi);
    
    
    /*if (phi == 0) {
        return ex1->getTau(indexX);
    } else if (phi == 45) {
        Matrix tau(1, 2);
        tau(0) = ex1->getTau(indexX);
        tau(1) = ey1->getTau(indexY);
        
        Exponential eTau(tau, 0, 90, 1, -1);
        
        return eTau.valueAt(phi);
    } else if (phi == 90) {
        return ey1->getTau(indexY);
    } else if (phi < 45) {
        Matrix tau(1, 2);
        tau(0) = ex1->getTau(indexX);
        tau(1) = ey1->getTau(indexY);
        
        Exponential eTau(tau, 0, 90, 1, -1);
        
        return eTau.valueAt(phi);
    } else {
        Matrix tau(1, 2);
        tau(0) = ex1->getTau(indexX);
        tau(1) = ey1->getTau(indexY);
        
        Exponential eTau(tau, 0, 90, 1, -1);
        
        return eTau.valueAt(phi);
    }*/
}

double endTau(Exponential** ex, Exponential**ey, int indexX, int indexY, double phi) {
    Exponential* ex1 = ex[indexY];
    Exponential* ex2 = ex[indexY + 1];
    Exponential* ey1 = ey[indexX];
    Exponential* ey2 = ey[indexX + 1];
    
    
    int blocks = findNumberOfBlocks(indexX, indexY);
    double hPhi = sqrt(8.0);
    
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
    
    Matrix tau(1, 3);
    tau(0) = ex1->getTau(indexX + 1);
    tau(1) = e.getTau(diagNth + 1);
    tau(2) = ey1->getTau(indexY + 1);
    
    Exponential eTau(tau, 0, 90, 2, -1);
    
    return eTau.valueAt(phi);
    
    
    if (phi == 0) {
        return ex1->getTau(indexX + 1);
    } else if (phi == 45) {
        Matrix tau(1, 2);
        tau(0) = ex1->getTau(indexX + 1);
        tau(1) = ey1->getTau(indexY + 1);
        
        Exponential eTau(tau, 0, 90, 1, -1);
        
        return eTau.valueAt(phi);
    } else if (phi == 90) {
        return ey1->getTau(indexY + 1);
    } else if (phi < 45) {
        Matrix tau(1, 2);
        tau(0) = ex1->getTau(indexX + 1);
        tau(1) = ey1->getTau(indexY + 1);
        
        Exponential eTau(tau, 0, 90, 1, -1);
        
        return eTau.valueAt(phi);
    } else {
        Matrix tau(1, 2);
        tau(0) = ex1->getTau(indexX + 1);
        tau(1) = ey1->getTau(indexY + 1);
        
        Exponential eTau(tau, 0, 90, 1, -1);
        
        return eTau.valueAt(phi);
    }
}

int findNumberOfBlocks(int indexX, int indexY) {
    int n = 1;
    
    int currentIndexX = indexX;
    int currentIndexY = indexY;
    while (true) {
        if (currentIndexX < 5 - 1 && currentIndexY < 5 - 1) {
            currentIndexX++;
            currentIndexY++;
            n++;
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
            n++;
        } else {
            break;
        }
    }
    
    return n;
}

double valueAt(Exponential** ex, Exponential** ey, double x, double y, const double xStartFrom, const double yStartFrom, const double hx, const double hy) {
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
    
    
    
    
    
    double a = (l2(ex1->getTau(indexX), ey1->getTau(indexY), phi) * sinh(hPhi - r) + l2(ex1->getTau(indexX + 1), ey1->getTau(indexY + 1), phi) * sinh(r)) / sinh(hPhi);
    
    double b = (ex1->getF(indexX) - l2(ex1->getTau(indexX), ey1->getTau(indexY), phi)) * (hPhi - r) / hPhi;
    
    double c = (xi(ex2, ey2, indexX, indexY, xStartFrom, yStartFrom, hx, hy, phi) - l2(ex1->getTau(indexX + 1), ey1->getTau(indexY + 1), phi)) * r / hPhi;
    
    a = (startTau(ex, ey, indexX, indexY, phi) * sinh(hPhi - r) + endTau(ex, ey, indexX, indexY, phi) * sinh(r)) / sinh(hPhi);
    b = (ex1->getF(indexX) - startTau(ex, ey, indexX, indexY, phi)) * (hPhi - r) / hPhi;
    c = (xi(ex2, ey2, indexX, indexY, xStartFrom, yStartFrom, hx, hy, phi) - endTau(ex, ey, indexX, indexY, phi)) * r / hPhi;
    
    /*double prevX = xStartFrom + indexX * hx;
    double nextX = xStartFrom + (indexX + 1) * hx;
    
    //r = currentX;
    nextX = 2;
    prevX = 0;
    
    a = (ex1->getTau(indexX) * sinh(nextX - r) + ex1->getTau(indexX + 1) * sinh(r - prevX)) / sinh(2);
    b = (ex1->getF(indexX) - ex1->getTau(indexX)) * (nextX - r) / hx;
    c = (ex1->getF(indexX + 1) - ex1->getTau(indexX + 1)) * (r - prevX) / hx;*/
    
    return a + b + c;
    return 0;
}

Matrix getOctaveMatrix(Exponential** ex, Exponential** ey, const double xStartFrom, const double yStartFrom, const int n, const int k, const int plotX, const int plotY, const double hx, const double hy, const double plotHx, const double plotHy) {
    Matrix v(plotY, plotX, 0.0);
    
    for (int i = 0; i < plotY; i++) {
        for (int j = 0; j < plotX; j++) {
            double x = xStartFrom + j * plotHx;
            double y = yStartFrom + i * plotHy;
            
            v(i, j) = valueAt(ex, ey, x, y, xStartFrom, yStartFrom, hx, hy);
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

    cout << "A: " << findAngle(1, 0) << "; " << findStep(findAngle(1, 0), hx, hy) << endl;
    cout << "A: " << findAngle(1, 0.2) << "; " << findStep(findAngle(1, 0.2), hx, hy) << endl;
    cout << "A: " << findAngle(1, 0.4) << "; " << findStep(findAngle(1, 0.4), hx, hy) << endl;
    cout << "A: " << findAngle(1, 0.6) << "; " << findStep(findAngle(1, 0.6), hx, hy) << endl;
    cout << "A: " << findAngle(1, 0.8) << "; " << findStep(findAngle(1, 0.8), hx, hy) << endl;
    cout << "A: " << findAngle(1, 1) << "; " << findStep(findAngle(1, 1), hx, hy) << endl;
    cout << "A: " << findAngle(1, 1.2) << "; " << findStep(findAngle(1, 1.2), hx, hy) << endl;
    cout << "A: " << findAngle(1, 1.5) << "; " << findStep(findAngle(1, 1.5), hx, hy) << endl;
    cout << "A: " << findAngle(0.2, 1.9) << "; " << findStep(findAngle(0.2, 1.9), hx, hy) << endl;
    
    valueAt(ex, ey, -5, 1, xStartFrom, yStartFrom, hx, hy);
    valueAt(ex, ey, -4, 1, xStartFrom, yStartFrom, hx, hy);
    valueAt(ex, ey, -3, 1, xStartFrom, yStartFrom, hx, hy);
    valueAt(ex, ey, -2, 1, xStartFrom, yStartFrom, hx, hy);
    valueAt(ex, ey, -1, 1, xStartFrom, yStartFrom, hx, hy);
    valueAt(ex, ey, 0, 1, xStartFrom, yStartFrom, hx, hy);
    valueAt(ex, ey, 1, 1, xStartFrom, yStartFrom, hx, hy);
    valueAt(ex, ey, 2, 1, xStartFrom, yStartFrom, hx, hy);
    valueAt(ex, ey, 3, 1, xStartFrom, yStartFrom, hx, hy);
    
    cout << "R: " << l(1, 2, 3, 0) << endl;
    cout << "R: " << l(1, 2, 3, 30) << endl;
    cout << "R: " << l(1, 2, 3, 45) << endl;
    cout << "R: " << l(1, 2, 3, 75) << endl;
    cout << "R: " << l(1, 2, 3, 90) << endl;
    
    cout << "Xi: " << 0.0 * 2 / 45.0 << endl;
    cout << "Xi: " << 10.0 * 2 / 45.0 << endl;
    cout << "Xi: " << 30.0 * 2 / 45.0 << endl;
    cout << "Xi: " << 40.0 * 2 / 45.0 << endl;
    cout << "Xi: " << 45.0 * 2 / 45.0 << endl;
    
    for (int i = 0; i <= 4; i++) {
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
    cout << endl;
    
    return octave_value(v);
}

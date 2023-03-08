#ifndef VARIABLES_H
#define VARIABLES_H
#include "Variables.h"


//CreateGrid
int** Grid(int fx, int fy, int lx, int ly, int m)
{
    int delx = (lx - fx) / m;
    int dely = (ly - fy) / m;

    int n = m + 1;
    int** ptr = 0;
    ptr = new int* [n ^ 2];
    int x = fx, y = fy;
    for (int i = 0; i < n; i++)
    {
        y = fy;
        ptr[i] = new int[2];
        for (int j = 0; j < n; j++)
        {
            ptr[i][0] = x;
            ptr[i][1] = y;
            y += dely;
        }
        x += delx;
    }
    return ptr;
}

//rotate/flip a quadrant appropriately
void rot(int n, int* x, int* y, int rx, int ry) {
    if (ry == 0) {
        if (rx == 1) {
            *x = n - 1 - *x;
            *y = n - 1 - *y;
        }

        //Swap x and y
        int t = *x;
        *x = *y;
        *y = t;
    }
}

//convert (x,y) to d
int xy2d(int n, int x, int y) {
    int rx, ry, s, d = 0;
    for (s = n / 2; s > 0; s /= 2) {
        rx = (x & s) > 0;
        ry = (y & s) > 0;
        d += s * s * ((3 * rx) ^ ry);
        rot(n, &x, &y, rx, ry);
    }
    return d;
}

//convert d to (x,y)
void d2xy(int n, int d, int* x, int* y) {
    int rx, ry, s, t = d;
    *x = *y = 0;
    for (s = 1; s < n; s *= 2) {
        rx = 1 & (t / 2);
        ry = 1 & (t ^ rx);
        rot(s, x, y, rx, ry);
        *x += s * rx;
        *y += s * ry;
        t /= 4;
    }
}

//Position Vector Reference configuration
Eigen::VectorXd X(Eigen::VectorXd x1, Eigen::VectorXd x2, double exi)
{
    //Shape function matrix
    Eigen::MatrixXd N = Eigen::MatrixXd::Zero(2, 4);
    N(0, 0) = (1 - exi) / 2;
    N(0, 2) = (1 + exi) / 2;

    N(1, 1) = (1 - exi) / 2;
    N(1, 3) = (1 + exi) / 2;

    Eigen::VectorXd xi = Eigen::VectorXd::Zero(4);
    xi(0) = x1(0);
    xi(1) = x1(1);
    xi(2) = x2(0);
    xi(3) = x2(1);

    Eigen::VectorXd y(2);
    y = N * xi;

    return y;
}

//Shape functions for displacement vectors
Eigen::VectorXd u(Eigen::VectorXd x1, Eigen::VectorXd x2, double exi, double h)
{
    //Shape function matrix
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(2, 6);

    H(0, 0) = (1 - exi) / 2;
    H(0, 3) = (1 + exi) / 2;;

    H(1, 1) = (1.0 / 4) * (2 - 3 * exi + 2 * pow(exi, 3));
    H(1, 4) = (1.0 / 4) * (2 + 3 * exi - 2 * pow(exi, 3));

    H(1, 2) = (h / 8) * (pow(exi, 2) - 1) * (exi - 1);
    H(1, 5) = (h / 8) * (pow(exi, 2) - 1) * (exi + 1);

    Eigen::VectorXd xi = Eigen::VectorXd::Zero(6);
    xi(0) = x1(0);
    xi(1) = x1(1);
    xi(2) = x1(2);
    xi(3) = x2(0);
    xi(4) = x2(1);
    xi(5) = x2(2);

    Eigen::VectorXd y;

    y = H * xi;

    return y;
}

//Position vector current configuration
Eigen::VectorXd x(Eigen::VectorXd x1, Eigen::VectorXd x2, double exi)
{
    return (X(x1, x2, exi) + u(x1, x2, exi, x2(0) - x1(0)));
}

//First derivative of position vector Reference configuration
Eigen::VectorXd dXdexi(Eigen::VectorXd x1, Eigen::VectorXd x2, double exi)
{
    Eigen::MatrixXd N = Eigen::MatrixXd::Zero(2, 4);
    N(0, 0) = -(1.0 / 2);
    N(0, 2) = (1.0 / 2);

    N(1, 1) = -(1.0 / 2);
    N(1, 3) = (1.0 / 2);

    Eigen::VectorXd xi = Eigen::VectorXd::Zero(4);
    xi(0) = x1(0);
    xi(1) = x1(1);
    xi(2) = x2(0);
    xi(3) = x2(1);

    Eigen::VectorXd y(2);
    y = N * xi;

    return y;
}

//First derivative of displacement function
Eigen::VectorXd dudexi(Eigen::VectorXd x1, Eigen::VectorXd x2, double exi, double h)
{
    //Shape function matrix
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(2, 6);

    H(0, 0) = - 1 / 2;
    H(0, 3) = 1 / 2;

    H(1, 1) = (1.0 / 4 * h) * 6 * (pow(exi, 2) - 1);
    H(1, 4) = (1.0 / 4 * h) * 6 * (1 - pow(exi, 2));

    H(1, 2) = (1 / 4) * (3 * pow(exi, 2) - 2 * exi - 1);
    H(1, 5) = (1 / 4) * (3 * pow(exi, 2) + 2 * exi - 1);

    Eigen::VectorXd xi = Eigen::VectorXd::Zero(6);
    xi(0) = x1(0);
    xi(1) = x1(1);
    xi(2) = x1(2);
    xi(3) = x2(0);
    xi(4) = x2(1);
    xi(5) = x2(2);

    Eigen::VectorXd y;

    y = H * xi;

    return y;
}

//Second derivative of displacement function
Eigen::VectorXd dudexi2(Eigen::VectorXd x1, Eigen::VectorXd x2, double exi, double h)
{
    //Shape function matrix
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(2, 6);

    H(0, 0) = 0;
    H(0, 3) = 0;

    H(1, 1) = -(1.0 / pow(h, 2)) * 6 * exi;
    H(1, 4) = -(1.0 / h) * (3 * exi - 1);

    H(1, 2) = (1.0 / pow(h, 2)) * 6 * exi;
    H(1, 5) = (1.0 / h) * (3 * exi + 1);

    Eigen::VectorXd xi = Eigen::VectorXd::Zero(6);
    xi(0) = x1(0);
    xi(1) = x1(1);
    xi(2) = x1(2);
    xi(3) = x2(0);
    xi(4) = x2(1);
    xi(5) = x2(2);

    Eigen::VectorXd y;

    y = H * xi;

    return y;
}

//First derivative of position vector current configuration
Eigen::VectorXd dxdexi(Eigen::VectorXd x1, Eigen::VectorXd x2, double exi)
{
    return (dXdexi(x1, x2, exi) + dudexi(x1, x2, exi, x2(0) - x1(0)));
}

//Double derivative of position vector current configuration
Eigen::VectorXd dxdexi2(Eigen::VectorXd x1, Eigen::VectorXd x2, double exi)
{
    return dudexi2(x1, x2, exi, x2(0) - x1(0));
}

//Local contact search
Eigen::VectorXd ContactPoints(Eigen::VectorXd x1, Eigen::VectorXd x2, Eigen::VectorXd y1, Eigen::VectorXd y2, int nbodies)
{
    Eigen::VectorXd bb = x2 + x1;
    Eigen::VectorXd tb = x2 - x1;
    Eigen::VectorXd b = y2 + y1;
    Eigen::VectorXd t = y2 - y1;
    std::cout << "b Bar" << std::endl;
    std::cout << bb << std::endl;
    std::cout << "t Bar" << std::endl;
    std::cout << tb << std::endl;
    std::cout << "b" << std::endl;
    std::cout << b << std::endl;
    std::cout << "t" << std::endl;
    std::cout << t << std::endl;
    Eigen::VectorXd exi = Eigen::VectorXd::Zero(nbodies);
    std::cout << (tb.dot(tb)) * (t.dot(t)) - pow(tb.dot(t), 2) << std::endl;
    exi(0) = -(bb - b).dot((tb * (tb.dot(t)) - t * (tb.dot(tb))) / (tb.dot(tb) * (t.dot(t)) - pow(tb.dot(t), 2)));
    exi(1) = (bb - b).dot((t * (tb.dot(t)) - tb * (t.dot(t))) / (tb.dot(tb) * (t.dot(t)) - pow(tb.dot(t), 2)));
    std::cout << "exi and exibar" << std::endl;
    std::cout << exi << std::endl;

    return exi;
}

double MinimumDistance(Eigen::VectorXd x1, Eigen::VectorXd x2, Eigen::VectorXd y1, Eigen::VectorXd y2, Eigen::VectorXd exi)
{
    Eigen::VectorXd xm(2);
    xm = x(x1, x2, exi(0));
    Eigen::VectorXd xs(2);
    xs = x(y1, y2, exi(1));
    double d = (xm - xs).norm();
    std::cout << d << std::endl;
    
    return d;
}

Eigen::MatrixXd EvaluateHTildeexiMatrix(Eigen::VectorXd x1, Eigen::VectorXd x2, Eigen::VectorXd y1, Eigen::VectorXd y2, Eigen::VectorXd exi)
{
    //Shape function matrix
    Eigen::MatrixXd Hexi = Eigen::MatrixXd::Zero(2, 6);
    double h = x2(0) - x1(0);

    Hexi(0, 0) = -1 / 2;
    Hexi(0, 3) = 1 / 2;;

    Hexi(1, 1) = (1.0 / 4 * h) * 6 * (pow(exi(0), 2) - 1);
    Hexi(1, 4) = (1.0 / 4 * h) * 6 * (1 - pow(exi(0), 2));

    Hexi(1, 2) = (1.0 / 4) * (3 * pow(exi(0), 2) - 2 * exi(0) - 1);
    Hexi(1, 5) = (1.0 / 4) * (3 * pow(exi(0), 2) + 2 * exi(0) - 1);


    Eigen::MatrixXd Hexib = Eigen::MatrixXd::Zero(2, 6);
    h = y2(0) - y1(0);

    Hexib(0, 0) = -1.0 / 2;
    Hexib(0, 3) = 1.0 / 2;

    Hexib(1, 1) = (1.0 / 4 * h) * 6 * (pow(exi(1), 2) - 1);
    Hexib(1, 4) = (1.0 / 4 * h) * 6 * (1 - pow(exi(1), 2));

    Hexib(1, 2) = (1.0 / 4) * (3 * pow(exi(1), 2) - 2 * exi(1) - 1);
    Hexib(1, 5) = (1.0 / 4) * (3 * pow(exi(1), 2) + 2 * exi(1) - 1);

    Eigen::MatrixXd HTilde(2, 12);
    for (int i = 0; i < 6; i++)
    {
        HTilde(0, i) = Hexi(0, i);
        HTilde(1, i) = Hexi(1, i);

        HTilde(0, i + 6) = Hexib(0, i);
        HTilde(1, i + 6) = Hexib(1, i);
    }

    return HTilde;
}

Eigen::MatrixXd EvaluateAMatrix(Eigen::VectorXd x1, Eigen::VectorXd x2, Eigen::VectorXd y1, Eigen::VectorXd y2, Eigen::VectorXd exi)
{
    Eigen::MatrixXd A(2, 2);
    Eigen::VectorXd dxdxib(2);
    dxdxib = dxdexi(x1, x2, exi(0));
    Eigen::VectorXd dxdxi(2);
    dxdxi = dxdexi(y1, y2, exi(1));

    A(0, 0) = -dxdxi.dot(dxdxi);
    A(0, 1) = dxdxib.dot(dxdxi);
    A(1, 0) = -A(0, 1);
    A(1, 1) = dxdxib.dot(dxdxib);

    std::cout << A << std::endl;

    return A;
}

Eigen::MatrixXd EvaluateBMatrix(Eigen::VectorXd x1, Eigen::VectorXd x2, Eigen::VectorXd y1, Eigen::VectorXd y2, Eigen::VectorXd exi)
{
    Eigen::MatrixXd B(2, 4);
    Eigen::VectorXd dxdxib(2);
    dxdxib = dxdexi(x1, x2, exi(0));
    Eigen::VectorXd dxdxi(2);
    dxdxi = dxdexi(y1, y2, exi(1));

    B(0, 0) = dxdxi(0);
    B(0, 1) = dxdxi(1);
    B(0, 2) = -dxdxi(0);
    B(0, 3) = -dxdxi(1);

    B(1, 0) = dxdxib(0);
    B(1, 1) = dxdxib(1);
    B(1, 2) = -dxdxib(0);
    B(1, 3) = -dxdxib(1);

    return B;
}

Eigen::MatrixXd EvaluateCMatrix(Eigen::VectorXd x1, Eigen::VectorXd x2, Eigen::VectorXd y1, Eigen::VectorXd y2, Eigen::VectorXd exi)
{
    Eigen::VectorXd xm(2);
    xm = x(x1, x2, exi(0));
    Eigen::VectorXd xs(2);
    xs = x(y1, y2, exi(1));

    Eigen::MatrixXd C(2, 4);
    C(0, 0) = -(xm(0) - xs(0));
    C(0, 1) = -(xm(0) - xs(0));
    C(0, 2) = 0;
    C(0, 3) = 0;

    C(1, 0) = 0;
    C(1, 1) = 0;
    C(1, 2) = -(xm(0) - xs(0));
    C(1, 3) = -(xm(0) - xs(0));

    return C;

}

Eigen::MatrixXd EvaluateDMatrix(Eigen::MatrixXd A, Eigen::MatrixXd B, Eigen::MatrixXd C, Eigen::VectorXd x1, Eigen::VectorXd x2, Eigen::VectorXd y1, Eigen::VectorXd y2, Eigen::VectorXd exi)
{
    //Shape function matrix
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(2, 6);
    double h = x2(0) - x1(0);

    H(0, 0) = -1 / 2;
    H(0, 3) = 1 / 2;;

    H(1, 1) = (1.0 / 4 * h) * 6 * (pow(exi(0), 2) - 1);
    H(1, 4) = (1.0 / 4 * h) * 6 * (1 - pow(exi(0), 2));

    H(1, 2) = (1 / 4) * (3 * pow(exi(0), 2) - 2 * exi(0) - 1);
    H(1, 5) = (1 / 4) * (3 * pow(exi(0), 2) + 2 * exi(0) - 1);

    Eigen::MatrixXd Hb = Eigen::MatrixXd::Zero(2, 6);
    h = y2(0) - y1(0);

    Hb(0, 0) = -1 / 2;
    Hb(0, 3) = 1 / 2;;

    Hb(1, 1) = (1.0 / 4 * h) * 6 * (pow(exi(1), 2) - 1);
    Hb(1, 4) = (1.0 / 4 * h) * 6 * (1 - pow(exi(1), 2));

    Hb(1, 2) = (1 / 4) * (3 * pow(exi(0), 2) - 2 * exi(1) - 1);
    Hb(1, 5) = (1 / 4) * (3 * pow(exi(0), 2) + 2 * exi(1) - 1);

    Eigen::MatrixXd HH = Eigen::MatrixXd::Zero(4, 12);
    for (int i = 0; i < 6; i++)
    {
        HH(0, i) = H(0, i);
        HH(1, i) = H(1, i);
        HH(2, i + 6) = Hb(0, i);
        HH(3, i + 6) = Hb(0, i);
    }

    Eigen::MatrixXd Hexi = Eigen::MatrixXd::Zero(2, 6);

    Hexi(0, 0) = -1 / 2;
    Hexi(0, 3) = 1 / 2;;

    Hexi(1, 1) = (1.0 / 4 * h) * 6 * (pow(exi(0), 2) - 1);
    Hexi(1, 4) = (1.0 / 4 * h) * 6 * (1 - pow(exi(0), 2));

    Hexi(1, 2) = (1.0 / 4) * (3 * pow(exi(0), 2) - 2 * exi(0) - 1);
    Hexi(1, 5) = (1.0 / 4) * (3 * pow(exi(0), 2) + 2 * exi(0) - 1);

    Eigen::MatrixXd Hexib = Eigen::MatrixXd::Zero(2, 6);

    Hexib(0, 0) = -1 / 2;
    Hexib(0, 3) = 1 / 2;;

    Hexib(1, 1) = (1.0 / 4 * h) * 6 * (pow(exi(1), 2) - 1);
    Hexib(1, 4) = (1.0 / 4 * h) * 6 * (1 - pow(exi(1), 2));

    Hexib(1, 2) = (1 / 4) * (3 * pow(exi(1), 2) - 2 * exi(1) - 1);
    Hexib(1, 5) = (1 / 4) * (3 * pow(exi(1), 2) + 2 * exi(1) - 1);

    Eigen::MatrixXd HHexi = Eigen::MatrixXd::Zero(4, 12);
    for (int i = 0; i < 6; i++)
    {
        HHexi(0, i) = Hexi(0, i);
        HHexi(1, i) = Hexi(1, i);
        HHexi(2, i + 6) = Hexib(0, i);
        HHexi(3, i + 6) = Hexib(0, i);
    }

    Eigen::MatrixXd D(2, 12);
    D = (A.inverse()) * (B * HH + C * HHexi);

    return D;
}

Eigen::VectorXd EvaluateNormalVector(Eigen::VectorXd x1, Eigen::VectorXd x2, Eigen::VectorXd y1, Eigen::VectorXd y2, Eigen::VectorXd exi)
{
    Eigen::VectorXd xm(2);
    xm = x(x1, x2, exi(0));
    Eigen::VectorXd xs(2);
    xs = x(y1, y2, exi(1));
    Eigen::VectorXd n(2);
    n = (xm - xs) / (xm - xs).norm();
    std::cout << n << std::endl;

    return n;
}

Eigen::MatrixXd EvaluateEMatrix(Eigen::VectorXd x1, Eigen::VectorXd x2, Eigen::VectorXd y1, Eigen::VectorXd y2, Eigen::VectorXd exi, Eigen::VectorXd n, Eigen::VectorXd d, Eigen::VectorXd db)
{
    Eigen::MatrixXd Hexi = Eigen::MatrixXd::Zero(2, 6);
    double h = x2(0) - x1(0);

    Hexi(0, 0) = -1.0 / 2;
    Hexi(0, 3) = 1.0 / 2;;

    Hexi(1, 1) = (1.0 / 4 * h) * 6 * (pow(exi(0), 2) - 1);
    Hexi(1, 4) = (1.0 / 4 * h) * 6 * (1 - pow(exi(0), 2));

    Hexi(1, 2) = (1.0 / 4) * (3 * pow(exi(0), 2) - 2 * exi(0) - 1);
    Hexi(1, 5) = (1.0 / 4) * (3 * pow(exi(0), 2) + 2 * exi(0) - 1);

    Eigen::MatrixXd Hexib = Eigen::MatrixXd::Zero(2, 6);
    h = y2(0) - y2(0);

    Hexib(0, 0) = -1.0 / 2;
    Hexib(0, 3) = 1.0 / 2;

    Hexib(1, 1) = (1.0 / 4 * h) * 6 * (pow(exi(1), 2) - 1);
    Hexib(1, 4) = (1.0 / 4 * h) * 6 * (1 - pow(exi(1), 2));

    Hexib(1, 2) = (1.0 / 4) * (3 * pow(exi(1), 2) - 2 * exi(1) - 1);
    Hexib(1, 5) = (1.0 / 4) * (3 * pow(exi(1), 2) + 2 * exi(1) - 1);

    Eigen::MatrixXd E(12, 12), Esub1(6, 12), Esub2(6, 12);
    Esub1 = -Hexi.transpose() * n * d;
    Esub2 = Hexib.transpose() * n * db;
    E << Esub1,
        Esub2;

    return E;

}

Eigen::MatrixXd EvaluateFMatrix(Eigen::VectorXd x1, Eigen::VectorXd x2, Eigen::VectorXd y1, Eigen::VectorXd y2, Eigen::VectorXd exi, Eigen::VectorXd n, Eigen::VectorXd d, Eigen::VectorXd db)
{
    Eigen::MatrixXd F(12, 12);
    F = db.transpose() * dudexi2(x1, x2, exi(0), x2(0) - x1(0)) * n * db - d.transpose() * n.transpose() * dudexi2(y1, y2, exi(1), y2(0) - y1(0)) * d;

    return F;
}

Eigen::MatrixXd EvaluateGMatrix(Eigen::VectorXd x1, Eigen::VectorXd x2, Eigen::VectorXd y1, Eigen::VectorXd y2, Eigen::VectorXd exi, Eigen::VectorXd n, Eigen::VectorXd d, Eigen::VectorXd db, Eigen::MatrixXd HTilde)
{
    Eigen::MatrixXd HTildeexi(2, 12);
    HTildeexi = EvaluateHTildeexiMatrix(x1, x2, y1, y2, exi);

    Eigen::MatrixXd Gsub1(12, 2), Gsub2(2, 12), G(12, 12);
    G = (HTildeexi.transpose() + db.transpose() * dxdexi(x1, x2, exi(0)).transpose() - d.transpose() * dxdexi(y1, y2, exi(1)).transpose()) * (Eigen::MatrixXd::Identity(2, 2) - n * n.transpose()) * (HTildeexi + dxdexi(x1, x2, exi(0)) * db -
        dxdexi(y1, y2, exi(1)) * d);

    return G;
}

Eigen::MatrixXd EvaluateHTildeMatrix(Eigen::VectorXd x1, Eigen::VectorXd x2, Eigen::VectorXd y1, Eigen::VectorXd y2, Eigen::VectorXd exi)
{
    //Shape function matrix
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(2, 6);
    double h = x2(0) - x1(0);

    H(0, 0) = (1 - exi(0)) / 2;
    H(0, 3) = (1 + exi(0)) / 2;;

    H(1, 1) = (1.0 / 4) * (2 - 3 * exi(0) + 2 * pow(exi(0), 3));
    H(1, 4) = (1.0 / 4) * (2 + 3 * exi(0) - 2 * pow(exi(0), 3));

    H(1, 2) = (h / 8) * (pow(exi(0), 2) - 1) * (exi(0) - 1);
    H(1, 5) = (h / 8) * (pow(exi(0), 2) - 1) * (exi(0) + 1);

    Eigen::MatrixXd Hb = Eigen::MatrixXd::Zero(2, 6);
    h = y2(0) - y1(0);

    H(0, 0) = (1 - exi(1)) / 2;
    H(0, 3) = (1 + exi(1)) / 2;;

    H(1, 1) = (1.0 / 4) * (2 - 3 * exi(1) + 2 * pow(exi(1), 3));
    H(1, 4) = (1.0 / 4) * (2 + 3 * exi(1) - 2 * pow(exi(1), 3));

    H(1, 2) = (h / 8) * (pow(exi(1), 2) - 1) * (exi(1) - 1);
    H(1, 5) = (h / 8) * (pow(exi(1), 2) - 1) * (exi(1) + 1);

    Eigen::MatrixXd HTilde(2, 12);
    for (int i = 0; i < 6; i++)
    {
        HTilde(0, i) = H(0, i);
        HTilde(1, i) = H(1, i);

        HTilde(0, i + 6) = H(0, i);
        HTilde(1, i + 6) = H(1, i);
    }

    return HTilde;
}



#endif

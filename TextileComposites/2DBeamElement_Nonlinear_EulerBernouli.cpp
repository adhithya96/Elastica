#ifndef VARIABLES_H
#define VARIABLES_H

#include "Variables.h"

//INP file 
NLEBBE2D::NLEBBE2D()
{
    NLEBBE2D::NNODE = 17;
    NLEBBE2D::NELEM = 16;
    NLEBBE2D::NMAT = 2;
    NLEBBE2D::NDOF = 3;
    NLEBBE2D::NLS = 12;

    NLEBBE2D::NODE = Eigen::MatrixXd::Zero(NLEBBE2D::NNODE, 3);
    NLEBBE2D::ELEM = Eigen::MatrixXd::Zero(NLEBBE2D::NELEM, 3);

    NLEBBE2D::E = 30e6;
    NLEBBE2D::nu = 0.33;

    NLEBBE2D::h = 1;
    NLEBBE2D::b = 1;

    NLEBBE2D::NODE(0, 0) = 0;
    NLEBBE2D::NODE(0, 1) = 0;
    NLEBBE2D::NODE(0, 2) = 0;
    NLEBBE2D::NODE(1, 0) = 6.25;
    NLEBBE2D::NODE(1, 1) = 0;
    NLEBBE2D::NODE(1, 2) = 0;
    NLEBBE2D::NODE(2, 0) = 12.5;
    NLEBBE2D::NODE(2, 1) = 0;
    NLEBBE2D::NODE(2, 2) = 0;
    NLEBBE2D::NODE(3, 0) = 18.75;
    NLEBBE2D::NODE(3, 1) = 0;
    NLEBBE2D::NODE(3, 2) = 0;
    NLEBBE2D::NODE(4, 0) = 25;
    NLEBBE2D::NODE(4, 1) = 0;
    NLEBBE2D::NODE(4, 2) = 0;
    NLEBBE2D::NODE(5, 0) = 31.25;
    NLEBBE2D::NODE(5, 1) = 0;
    NLEBBE2D::NODE(5, 2) = 0;
    NLEBBE2D::NODE(6, 0) = 37.5;
    NLEBBE2D::NODE(6, 1) = 0;
    NLEBBE2D::NODE(6, 2) = 0;
    NLEBBE2D::NODE(7, 0) = 43.75;
    NLEBBE2D::NODE(7, 1) = 0;
    NLEBBE2D::NODE(7, 2) = 0;
    NLEBBE2D::NODE(8, 0) = 50;
    NLEBBE2D::NODE(8, 1) = 0;
    NLEBBE2D::NODE(8, 2) = 0;
    NLEBBE2D::NODE(9, 0) = 56.25;
    NLEBBE2D::NODE(9, 1) = 0;
    NLEBBE2D::NODE(9, 2) = 0;
    NLEBBE2D::NODE(10, 0) = 62.5;
    NLEBBE2D::NODE(10, 1) = 0;
    NLEBBE2D::NODE(10, 2) = 0;
    NLEBBE2D::NODE(11, 0) = 68.75;
    NLEBBE2D::NODE(11, 1) = 0;
    NLEBBE2D::NODE(11, 2) = 0;
    NLEBBE2D::NODE(12, 0) = 75;
    NLEBBE2D::NODE(12, 1) = 0;
    NLEBBE2D::NODE(12, 2) = 0;
    NLEBBE2D::NODE(13, 0) = 81.25;
    NLEBBE2D::NODE(13, 1) = 0;
    NLEBBE2D::NODE(13, 2) = 0;
    NLEBBE2D::NODE(14, 0) = 87.5;
    NLEBBE2D::NODE(14, 1) = 0;
    NLEBBE2D::NODE(14, 2) = 0;
    NLEBBE2D::NODE(15, 0) = 93.75;
    NLEBBE2D::NODE(15, 1) = 0;
    NLEBBE2D::NODE(15, 2) = 0;
    NLEBBE2D::NODE(16, 0) = 100;
    NLEBBE2D::NODE(16, 1) = 0;
    NLEBBE2D::NODE(16, 2) = 0;

    NLEBBE2D::ELEM(0, 0) = 1;
    NLEBBE2D::ELEM(0, 1) = 1;
    NLEBBE2D::ELEM(0, 2) = 2;
    NLEBBE2D::ELEM(1, 0) = 1;
    NLEBBE2D::ELEM(1, 1) = 2;
    NLEBBE2D::ELEM(1, 2) = 3;
    NLEBBE2D::ELEM(2, 0) = 1;
    NLEBBE2D::ELEM(2, 1) = 3;
    NLEBBE2D::ELEM(2, 2) = 4;
    NLEBBE2D::ELEM(3, 0) = 1;
    NLEBBE2D::ELEM(3, 1) = 4;
    NLEBBE2D::ELEM(3, 2) = 5;
    NLEBBE2D::ELEM(4, 0) = 1;
    NLEBBE2D::ELEM(4, 1) = 5;
    NLEBBE2D::ELEM(4, 2) = 6;
    NLEBBE2D::ELEM(5, 0) = 1;
    NLEBBE2D::ELEM(5, 1) = 6;
    NLEBBE2D::ELEM(5, 2) = 7;
    NLEBBE2D::ELEM(6, 0) = 1;
    NLEBBE2D::ELEM(6, 1) = 7;
    NLEBBE2D::ELEM(6, 2) = 8;
    NLEBBE2D::ELEM(7, 0) = 1;
    NLEBBE2D::ELEM(7, 1) = 8;
    NLEBBE2D::ELEM(7, 2) = 9;
    NLEBBE2D::ELEM(8, 0) = 1;
    NLEBBE2D::ELEM(8, 1) = 9;
    NLEBBE2D::ELEM(8, 2) = 10;
    NLEBBE2D::ELEM(9, 0) = 1;
    NLEBBE2D::ELEM(9, 1) = 10;
    NLEBBE2D::ELEM(9, 2) = 11;
    NLEBBE2D::ELEM(10, 0) = 1;
    NLEBBE2D::ELEM(10, 1) = 11;
    NLEBBE2D::ELEM(10, 2) = 12;
    NLEBBE2D::ELEM(11, 0) = 1;
    NLEBBE2D::ELEM(11, 1) = 12;
    NLEBBE2D::ELEM(11, 2) = 13;
    NLEBBE2D::ELEM(12, 0) = 1;
    NLEBBE2D::ELEM(12, 1) = 13;
    NLEBBE2D::ELEM(12, 2) = 14;
    NLEBBE2D::ELEM(13, 0) = 1;
    NLEBBE2D::ELEM(13, 1) = 14;
    NLEBBE2D::ELEM(13, 2) = 15;
    NLEBBE2D::ELEM(14, 0) = 1;
    NLEBBE2D::ELEM(14, 1) = 15;
    NLEBBE2D::ELEM(14, 2) = 16;
    NLEBBE2D::ELEM(15, 0) = 1;
    NLEBBE2D::ELEM(15, 1) = 16;
    NLEBBE2D::ELEM(15, 2) = 17;

    NLEBBE2D::CNODE = Eigen::MatrixXd::Zero(4, 2);

    //pinned
    NLEBBE2D::CNODE(0, 0) = 1;
    NLEBBE2D::CNODE(0, 1) = 1;
    NLEBBE2D::CNODE(1, 0) = 1;
    NLEBBE2D::CNODE(1, 1) = 2;
    NLEBBE2D::CNODE(2, 0) = 17;
    NLEBBE2D::CNODE(2, 1) = 1;
    NLEBBE2D::CNODE(3, 0) = 17;
    NLEBBE2D::CNODE(3, 1) = 2;
    //clamped
    /*NLEBBE2D::CNODE(0, 0) = 1;
    NLEBBE2D::CNODE(0, 1) = 1;
    NLEBBE2D::CNODE(1, 0) = 1;
    NLEBBE2D::CNODE(1, 1) = 2;
    NLEBBE2D::CNODE(2, 0) = 1;
    NLEBBE2D::CNODE(2, 1) = 3;
    NLEBBE2D::CNODE(3, 0) = 17;
    NLEBBE2D::CNODE(3, 1) = 1;
    NLEBBE2D::CNODE(4, 0) = 17;
    NLEBBE2D::CNODE(4, 1) = 2;
    NLEBBE2D::CNODE(5, 0) = 17;
    NLEBBE2D::CNODE(5, 1) = 3;*/

    NLEBBE2D::vf = 0;
    NLEBBE2D::af = 0;
}






int NLEBBE2D::get_nen()
{
    return NLEBBE2D::NEN;
}

int NLEBBE2D::get_ndof()
{
    return NLEBBE2D::NDOF;
}

int NLEBBE2D::get_nnode()
{
    return NLEBBE2D::NNODE;
}

int NLEBBE2D::get_nelem()
{
    return NLEBBE2D::NELEM;
}

int NLEBBE2D::get_nls()
{
    return NLEBBE2D::NLS;
}

double NLEBBE2D::get_coordinates(int i, int j)
{
    return NLEBBE2D::NODE(i, j);
}

int NLEBBE2D::get_connectivity(int i, int j)
{
    return NLEBBE2D::ELEM(i, j);
}

int NLEBBE2D::get_cnode(int i, int j)
{
    return NLEBBE2D::CNODE(i, j);
}

double NLEBBE2D::get_modelprop(std::string str)
{
    if (str == "E")
        return NLEBBE2D::E;
    else if (str == "nu")
        return NLEBBE2D::nu;
    else if (str == "b")
        return NLEBBE2D::b;
    else if (str == "h")
        return NLEBBE2D::h;
}

double NLEBBE2D::get_loadprop(std::string str)
{
    if (str == "vf")
        return NLEBBE2D::vf;
    else if (str == "af")
        return NLEBBE2D::af;
}

void NLEBBE2D::set_loadprop(std::string str, double value)
{
    if (str == "vf")
        NLEBBE2D::vf = value;
    else if (str == "af")
        NLEBBE2D::af = value;
}


//Jacobian
double Jacobian(double xa, double xb)
{
    double J = (xb - xa) / 2;
    assert(J > 0);
    return J;
}

//Hermite Cubics
double HC(double x, int num, double xa, double xb)
{
    double h = xb - xa;
    if (num == 1)
        return 1 - 3 * pow((x - xa) / h, 2) + 2 * pow((x - xa) / h, 3);
    else if (num == 2)
        return -(x - xa) * pow((1 - (x - xa) / h), 2);
    else if (num == 3)
        return 3 * pow((x - xa) / h, 2) - 2 * pow((x - xa) / h, 3);
    else if (num == 4)
        return -(x - xa) * (pow((x - xa) / h, 2) - (x - xa) / h);
    else
        return -99;
}

//Function to calculate double derivative of Hermite Cubic Functions
double DDHC(double x, int num, double xa, double xb)
{
    double h = xb - xa;
    if (num == 1)
        return -(6 / pow(h, 2)) * (1 - 2 * (x - xa) / h);
    else if (num == 2)
        return -(2 / h) * (3 * (x - xa) / h - 2);
    else if (num == 3)
        return (6 / pow(h, 2)) * (1 - 2 * (x - xa) / h);
    else if (num == 4)
        return -(2 / h) * (3 * (x - xa) / h - 1);
    else
        return -99;
}

//Function to calculate single derivative of Hermite Cubic functions
double SDHC(double x, int num, double xa, double xb)
{
    double h = xb - xa;
    if (num == 1)
        return -(6 * (x - xa) / pow(h, 2)) * (1 - (x - xa) / h);
    else if (num == 2)
        return -(1 + 3 * pow((x - xa) / h, 2) - 4 * (x - xa) / h);
    else if (num == 3)
        return (6 * (x - xa) / pow(h, 2)) * (1 - (x - xa) / h);
    else if (num == 4)
        return -((x - xa) / h) * (3 * (x - xa) / h - 2);
    else
        return -99;
}

//Hat functions
double HF(double x, double num, double xa, double xb)
{
    double h = xb - xa;
    if (num == 1)
        return (xb - x) / h;
    else if (num == 2)
        return (x - xa) / h;
    else
        return -99;
}

//Function to evaluate single derivative of hat functions
double SDHF(double num, double xa, double xb)
{
    double h = xb - xa;
    if (num == 1)
        return -1 / h;
    else if (num == 2)
        return 1 / h;
    else
        return -99;
}

//To rearrange the elements of stiffness matrix and force vector
void NLEBBE2D::RearrangeElementStiffness_NLEBBE(Eigen::MatrixXd& k, Eigen::MatrixXd& t, Eigen::VectorXd& f)
{
    Eigen::MatrixXd temp(6, 6), temp2(6, 6);
    Eigen::VectorXd temp3(6);
    //Rearranging rows
    for (int j = 0; j < 6; j++)
    {
        //Stiffness
        temp(0, j) = k(0, j);
        temp(1, j) = k(2, j);
        temp(2, j) = k(3, j);
        temp(3, j) = k(1, j);
        temp(4, j) = k(4, j);
        temp(5, j) = k(5, j);
        //Tangent
        temp2(0, j) = t(0, j);
        temp2(1, j) = t(2, j);
        temp2(2, j) = t(3, j);
        temp2(3, j) = t(1, j);
        temp2(4, j) = t(4, j);
        temp2(5, j) = t(5, j);
        //Force
        temp3(j) = f(j);
    }
    //    Rearranging columns
    for (int i = 0; i < 6; i++)
    {
        //Stiffness
        k(i, 0) = temp(i, 0);
        k(i, 1) = temp(i, 2);
        k(i, 2) = temp(i, 3);
        k(i, 3) = temp(i, 1);
        k(i, 4) = temp(i, 4);
        k(i, 5) = temp(i, 5);
        //Tangent
        t(i, 0) = temp2(i, 0);
        t(i, 1) = temp2(i, 2);
        t(i, 2) = temp2(i, 3);
        t(i, 3) = temp2(i, 1);
        t(i, 4) = temp2(i, 4);
        t(i, 5) = temp2(i, 5);
    }
    f(0) = temp3(0);
    f(1) = temp3(2);
    f(2) = temp3(3);
    f(3) = temp3(1);
    f(4) = temp3(4);
    f(5) = temp3(5);
}



//How to calculate dw0/dx from previous iteration
double dw0dx(Eigen::VectorXd& U, double x, int a, int b, double xa, double xb)
{
    double dwdx;
    double wa = U(3 * a + 1);
    double wb = U(3 * b + 1);
    double theta_a = U(3 * a + 2);
    double theta_b = U(3 * b + 2);
    dwdx = wa * SDHC(x, 1, xa, xb) + theta_a * SDHC(x, 2, xa, xb) + wb * SDHC(x, 3, xa, xb) + theta_b * SDHC(x, 4, xa, xb);

    return dwdx;
}

double dU0dx(Eigen::VectorXd& U, int a, int b, double xa, double xb)
{
    double dudx;
    double ua = U(3 * a);
    double ub = U(3 * b);
    dudx = ua * SDHF(1, xa, xb) + ub * SDHF(2, xa, xb);

    return dudx;
}

double ddw0dx(Eigen::VectorXd& U, double x, int a, int b, double xa, double xb)
{
    double ddwdx;
    double wa = U(3 * a + 1);
    double wb = U(3 * b + 1);
    double theta_a = U(3 * a + 2);
    double theta_b = U(3 * b + 2);
    ddwdx = wa * DDHC(x, 1, xa, xb) + theta_a * DDHC(x, 2, xa, xb) + wb * DDHC(x, 3, xa, xb) + theta_b * DDHC(x, 4, xa, xb);

    return ddwdx;
}

Eigen::MatrixXd NLEBBE2D::TangentStiffnessMatrix_NLEBBE(Eigen::MatrixXd& k, double xa, double xb, double E, double nu, double A, double I, Eigen::VectorXd& U, int a, int b)
{
    //One point Gauss Quadrature rule
    double ow1 = 2;
    double ox1 = 0;
    //Calculate Axx, Dxx
    double Axx, Dxx;
    Axx = E * A;
    Dxx = E * I;

    double x0 = (xb - xa) / 2 * ox1 + (xb + xa) / 2;
    Eigen::MatrixXd t = Eigen::MatrixXd::Zero(6, 6);

    //T11
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            t(i, j) = k(i, j);
    //T12
    for (int i = 0; i < 2; i++)
        for (int j = 2; j < 6; j++)
            t(i, j) = 2 * k(i, j);
    //T21
    for (int i = 2; i < 6; i++)
        for (int j = 0; j < 2; j++)
            t(i, j) = k(i, j);
    //T22
    for (int i = 2; i < 6; i++)
        for (int j = 2; j < 6; j++)
        {
            //Integration points
            //This is used to indicate the map between x and exi
            //A subparametric map (in this case a linear map) is used to convert x to exi.
            double fxa = Axx * (dU0dx(U, a, b, xa, xb) + pow(dw0dx(U, x0, a, b, xa, xb), 2)) * SDHC(x0, i - 1, xa, xb) * SDHC(x0, j - 1, xa, xb);
            t(i, j) = k(i, j) + ow1 * fxa * Jacobian(xa, xb);
        }

    return t;
}


Eigen::MatrixXd NLEBBE2D::StiffnessMatrix_NLEBBE(double xa, double xb, double E, double nu, double A, double I, Eigen::VectorXd& U, int a, int b)
{
    //Two Point Gauss Quadrature
    double tx1 = -0.577350269;
    double tx2 = 0.577350269;
    double tw1 = 1;
    double tw2 = 1;
    //One point Gauss Quadrature
    double ox1 = 0;
    double ow1 = 2;
    //Calculate Axx, Dxx
    double Axx, Dxx;
    Axx = E * A;
    Dxx = E * I;

    Eigen::MatrixXd k = Eigen::MatrixXd::Zero(6, 6);

    double x1 = (xb - xa) / 2 * tx1 + (xb + xa) / 2;
    double x2 = (xb - xa) / 2 * tx2 + (xb + xa) / 2;
    double x0 = (xb - xa) / 2 * ox1 + (xb + xa) / 2;
    //Stiffness matrix elements
    //K11
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            //constant function
            //Exact integral can be evaluated since Axx is a constant
            k(i, j) = 2 * Axx * SDHF(i + 1, xa, xb) * SDHF(j + 1, xa, xb) * Jacobian(xa, xb);
    //K12
    for (int i = 0; i < 2; i++)
        for (int j = 2; j < 6; j++)
        {
            //Quadrature
            double fxa = 0.5 * Axx * dw0dx(U, x0, a, b, xa, xb) * SDHF(i + 1, xa, xb) * SDHC(x0, j - 1, xa, xb);
            k(i, j) = fxa * ow1 * Jacobian(xa, xb);
        }
    //K21
    for (int i = 2; i < 6; i++)
        for (int j = 0; j < 2; j++)
        {
            //Quadrature
            double fxa = Axx * dw0dx(U, x0, a, b, xa, xb) * SDHF(j + 1, xa, xb) * SDHC(x0, i - 1, xa, xb);
            k(i, j) = fxa * ow1 * Jacobian(xa, xb);
        }
    //K22
    for (int i = 2; i < 6; i++)
        for (int j = 2; j < 6; j++)
        {
            double gxa = Dxx * DDHC(x1, i - 1, xa, xb) * DDHC(x1, j - 1, xa, xb);
            double gxb = Dxx * DDHC(x2, i - 1, xa, xb) * DDHC(x2, j - 1, xa, xb);
            double fxa = 0.5 * Axx * pow(dw0dx(U, x0, a, b, xa, xb), 2) * SDHC(x0, i - 1, xa, xb) * SDHC(x0, j - 1, xa, xb);
            k(i, j) = (gxa * tw1 + gxb * tw2 + fxa * ow1) * Jacobian(xa, xb);
        }
    //TangentStiffnessEulerBernoulli(t,k,xa,xb,E,nu,base,height,U,a,b);

    return k;
}

//Local Force Vector
Eigen::VectorXd NLEBBE2D::LocalFoceVec_NLEBBE(double xa, double xb, int a, int b, double vf, double af, Eigen::VectorXd& U, double E, double nu)
{
    //Two Point Gauss Quadrature
    double tx1 = -0.577350269;
    double tx2 = 0.577350269;
    double tw1 = 1;
    double tw2 = 1;

    Eigen::VectorXd f(6);
    for (int i = 0; i < 6; i++)
        f(i) = 0;
    //Two things have to be done. Input axial load and transverse load in ReadInp File and use it in this function
    //For now q is taken as 1
    //2 point gaussian quadrature.
    //Axial load acting
    double x1 = ((xb - xa) / 2) * tx1 + (xb + xa) / 2;
    double x2 = ((xb - xa) / 2) * tx2 + (xb + xa) / 2;
    for (int i = 0; i < 2; i++)
    {
        //        std::cout<<x1<<std::endl;
        //        std::cout<<x2<<std::endl;
        //        std::cout<<x3<<std::endl;
        //        std::cout<<w1<<std::endl;
        //        std::cout<<w2<<std::endl;
        //        std::cout<<w3<<std::endl;
        double fxa = af * HF(x1, i + 1, xa, xb);
        double fxb = af * HF(x2, i + 1, xa, xb);
        f(i) = (tw1 * fxa + tw2 * fxb) * Jacobian(xa, xb);
    }

    //Transverse load
    for (int i = 2; i < 6; i++)
    {
        double fxa = vf * HC(x1, i - 1, xa, xb);
        double fxb = vf * HC(x2, i - 1, xa, xb);
        f(i) = (tw1 * fxa + tw2 * fxb) * Jacobian(xa, xb);
    }

    return f;
}

void NLEBBE2D::ApplyConstraints_NLEBBE(Eigen::SparseMatrix<double, Eigen::ColMajor>& T, Eigen::VectorXd& U, int NNODE, Eigen::VectorXd& R)
{
    for (int j = 0; j < NLEBBE2D::CNODE.rows(); j++)
    {
        int nodenum = NLEBBE2D::CNODE(j, 0);
        int dir = NLEBBE2D::CNODE(j, 1) - 1;
        for (int i = 0; i < 3 * NNODE; i++)
            T.coeffRef(3 * (nodenum - 1) + dir, i) = 0;
        R(3 * (nodenum - 1) + dir) = 0;
        T.coeffRef(3 * (nodenum - 1) + dir, 3 * (nodenum - 1) + dir) = 1;
    }
}

/*void PostProcessing(Eigen::VectorXd U, NonLinearEulerBernouliBeamElement NLEBBE, int fiter)
{
    std::fstream file1;
    std::string filename, temp; 
    std::stringstream s1;
    s1 << fiter;
    temp = s1.str();
    filename = "E:/Adhithya/MTech_Res_Thesis/Cpp/ThesisCode/TextileComposites/TextileComposites/2DBeam/Displacement" + temp + ".vtk";
    std::cout << filename << std::endl;
    file1.open(filename.c_str(), std::ios::out);
    if (file1.is_open())
    {
        file1 << "# vtk DataFile Version 2.0" << std::endl;
        file1 << "# Cantilever Beam " << std::endl;
        file1 << "ASCII" << std::endl;
        file1 << "DATASET UNSTRUCTURED_GRID" << std::endl;
        file1 << "POINTS " << NLEBBE.NNODE << " float" << std::endl;
        for (int i = 0; i < NLEBBE.NNODE; i++)
        {
            file1 << NLEBBE.NODE(i, 0) << " " << NLEBBE.NODE(i, 1) << std::endl;
        }
        file1 << std::endl;
        file1 << "CELLS " << NLEBBE.NELEM << " " << NLEBBE.NELEM * 3 << std::endl;
        for (int i = 0; i < NLEBBE.NELEM; i++)
        {
            file1 << "2 " << NLEBBE.ELEM(i, 1) - 1 << " " << NLEBBE.ELEM(i, 2) - 1 << std::endl;
        }
        file1 << std::endl;
        file1 << "CELL_TYPES " << NLEBBE.NELEM << std::endl;
        for (int i = 0; i < NLEBBE.NELEM; i++)
            file1 << "3 " << std::endl;
        file1 << std::endl;
        file1 << "POINT_DATA" << " " << NLEBBE.NNODE << std::endl;
        file1 << "VECTORS " << "U " <<"float " << std::endl;
        for (int i = 0; i < NLEBBE.NNODE; i++)
        {
            file1 << U(i) << " " << U(i + 1) << std::endl;
        }
        file1 << "SCALARS " << "Theta " << "float " << std::endl;
        for (int i = 0; i < NLEBBE.NNODE; i++)
        {
            file1 << U(i + 2) << std::endl;
        }
    }
    else
        std::cout << " File is not open " << std::endl;

}*/

#endif


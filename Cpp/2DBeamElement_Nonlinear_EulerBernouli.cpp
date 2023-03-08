#ifndef VARIABLES_H
#define VARIABLES_H

#include "Variables.h"

struct ForceVector_IntConstants
{
    double Q1;
    double Q2;
    double Q3;
    double Q4;
    double Q5;
    double Q6;
};

//Jacobian
double Jacobian(double xa, double xb)
{
    double J = (xb - xa) / 2;
    assert(J > 0);
    return J;
}

double Area(struct CrossSection* CS)
{
    if ((*CS).choice == "RECT")
    {
        double h = (*CS).Rect.height;
        double b = (*CS).Rect.width;
        return (b * h);
    }
    else if ((*CS).choice == "CIRCLE")
    {
        double r = (*CS).Cir.radius;
        return (M_PI * pow(r, 2));
    }
    else
        return -99;
}

double MomentOfInertia(struct CrossSection* CS)
{
    if ((*CS).choice == "RECT")
    {
        double h = (*CS).Rect.height;
        double b = (*CS).Rect.width;
        return (b * pow(h, 3) / 12);
    }
    else if ((*CS).choice == "CIRCLE")
    {
        double r = (*CS).Cir.radius;
        return (M_PI * pow(r, 4) / 4);
    }
    else
        return -99;
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
void RearrangeElementStiffness_NLEBBE(Eigen::MatrixXd& k, Eigen::MatrixXd& t, Eigen::VectorXd& f)
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

Eigen::MatrixXd TangentStiffnessMatrix_NLEBBE(Eigen::MatrixXd& k, double xa, double xb, double E, double nu, double A, double I, Eigen::VectorXd& U, int a, int b)
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


Eigen::MatrixXd StiffnessMatrix_NLEBBE(double xa, double xb, double E, double nu, double A, double I, Eigen::VectorXd& U, int a, int b)
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
Eigen::VectorXd LocalFoceVec_NLEBBE(double xa, double xb, int a, int b, double vf, double af, Eigen::VectorXd& U, double E, double nu)
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

void ApplyConstraints_NLEBBE(Eigen::SparseMatrix<double, Eigen::ColMajor>& T, Eigen::VectorXd& U, Eigen::MatrixXd& CNODE, int NNODE, Eigen::VectorXd& R)
{
    for (int j = 0; j < CNODE.rows(); j++)
    {
        int nodenum = CNODE(j, 0);
        int dir = CNODE(j, 1) - 1;
        for (int i = 0; i < 3 * NNODE; i++)
            T.coeffRef(3 * (nodenum - 1) + dir, i) = 0;
        R(3 * (nodenum - 1) + dir) = 0;
        T.coeffRef(3 * (nodenum - 1) + dir, 3 * (nodenum - 1) + dir) = 1;
    }
}

void PostProcessing(Eigen::VectorXd U, NonLinearEulerBernouliBeamElement NLEBBE, int fiter)
{
    std::fstream file1;
    std::string filename, temp; 
    std::stringstream s1;
    s1 << fiter;
    temp = s1.str();
    filename = "G:/MTech_Aerospace/MTech_Res_Thesis/Cpp/ThesisCode/TextileComposites/Displacement.vtu." + temp;
    std::cout << filename << std::endl;
    file1.open(filename.c_str(), std::ios::out);
    if (file1.is_open())
    {
        file1 << "# vtk DataFile Version 2.0" << std::endl;
        file1 << "Cantilever Beam " << std::endl;
        file1 << "ASCII" << std::endl;
        file1 << "DATASET STRUCTURED_POINTS" << std::endl;
        file1 << "POINTS " << NLEBBE.NNODE * 2 << " float" << std::endl;
        for (int i = 0; i < NLEBBE.NNODE; i++)
        {
            file1 << NLEBBE.NODE(i, 0) << " " << NLEBBE.NODE(i, 1) << std::endl;
        }
        file1 << "CELLS " << NLEBBE.NELEM * 3 << std::endl;
        for (int i = 0; i < NLEBBE.NELEM; i++)
        {
            file1 << "3 " << NLEBBE.ELEM(i, 1) - 1 << " " << NLEBBE.ELEM(i, 2) - 1 << std::endl;
        }
        file1 << "POINTDATA " << std::endl;
        for (int i = 0; i < NLEBBE.NNODE; i++)
        {
            file1 << U(i) << U(i + 1) << U(i + 2) << std::endl;
        }
    }
    else
        std::cout << " File is not open " << std::endl;

}

#endif


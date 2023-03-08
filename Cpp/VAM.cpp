
//Couple of things to figure out

//What units to use for k1
//FDM code for first and last node
//How to validate? Experimental data is not given. Just the graphs.
#ifndef VARIABLES_H
#define VARIABLES_H


#include "Variables.h"
#include <corecrt_math_defines.h>


//Q matrix
Eigen::MatrixXd ReducedStiffnessMatrix(double E1, double E2, double nu12, double nu21, double G12)
{
    Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(3, 3);
    Q(0, 0) = E1 / (1 - nu12 * nu21);
    Q(0, 1) = nu12 * E2 / (1 - nu12 * nu21);
    Q(1, 1) = E2 / (1 - nu12 * nu21);
    Q(1, 0) = Q(0, 1);
    Q(2, 2) = G12;

    //    std::cout<<Q<<std::endl;

    return Q;
}

Eigen::MatrixXd TransformationMatrix(double theta)
{
    double m = cos(M_PI * theta / 180);
    double n = sin(M_PI * theta / 180);

    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(3, 3);
    T(0, 0) = pow(m, 2);
    T(0, 1) = pow(n, 2);
    T(0, 2) = -2 * m * n;
    T(1, 0) = pow(n, 2);
    T(1, 1) = pow(m, 2);
    T(1, 2) = 2 * m * n;
    T(2, 0) = m * n;
    T(2, 1) = -m * n;
    T(2, 2) = pow(m, 2) - pow(n, 2);

    //    std::cout<<T<<std::endl;

    return T;
}

//A matrix
void AMatrix(Eigen::MatrixXd Q, double z0, double z1, Eigen::MatrixXd* ptr)
{
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            (*ptr)(i, j) += Q(i, j) * (z1 - z0);
}

//B matrix
void BMatrix(Eigen::MatrixXd Q, double z0, double z1, Eigen::MatrixXd* ptr)
{
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            (*ptr)(i, j) += (1.0 / 2.0) * Q(i, j) * (pow(z1, 2) - pow(z0, 2));
}

//D matrix
void DMatrix(Eigen::MatrixXd Q, double z0, double z1, Eigen::MatrixXd* ptr)
{
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            (*ptr)(i, j) += (1.0 / 3.0) * Q(i, j) * (pow(z1, 3) - pow(z0, 3));
}

//Variables needed for Calculating VAM simplified Stiffness matrix
/*Variables StiffnessVariableCalc(Eigen::MatrixXd A, Eigen::MatrixXd B, Eigen::MatrixXd D)
{
    Variables V;



    return V;
}*/

//Stiffness Matrices after simplification using VAM
/*Variables StiffnessMatrix_VAM(double E1, double E2, double nu12, double  G12, double nu21, double width, double height,
    Eigen::VectorXd Orient, int np, std::fstream& file1)
{
    struct Variables V;

    double b = width;

    //Midplane thickness and calculate hi
    double mt = np * height / 2;

    file1 << "            Midplane thickness  " << mt << std::endl;
    Eigen::VectorXd h(np + 1);
    h(0) = -mt;
    for (int i = 1; i < np + 1; i++)
        h(i) = h(i - 1) + height;
    //assert(h(np)==mt);
    file1 << "            Height of nth composite lamina  " << h(np) << std::endl;

    //Reduced stiffness matrix at different angles
    Eigen::MatrixXd Q = ReducedStiffnessMatrix(E1, E2, nu12, nu21, G12);
    file1 << "            Reduced Stiffness Matrix " << std::endl;
    for (int j = 0; j < 3; j++)
    {
        file1 << "                ";
        for (int k = 0; k < 3; k++)
        {
            file1 << Q(j, k) << "  ";
        }
        file1 << std::endl;
    }
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(3, 3);
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, 3);
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(3, 3);
    for (int i = 0; i < Orient.size(); i++)
    {
        Eigen::MatrixXd T = TransformationMatrix(Orient(i));
        /*file1 << "                Evaluating Transformation matrix at " << Orient(i) << std::endl;
        for (int j = 0; j < 3; j++)
        {
            file1 << "                    ";
            for (int k = 0; k < 3; k++)
            {
                file1 << T(j, k) << "  ";
            }
            file1 << std::endl;
        }
        Eigen::MatrixXd Qdash = T * Q * (T.transpose());
        /*file1 << "                Evaluating Q matrix at " << Orient(i) << std::endl;
        for (int j = 0; j < 3; j++)
        {
            file1 << "                    ";
            for (int k = 0; k < 3; k++)
            {
                file1 << Qdash(j, k) << "  ";
            }
            file1 << std::endl;
        }
        double z0 = h(i);
        double z1 = h(i + 1);
        AMatrix(Qdash, z0, z1, &A);
        BMatrix(Qdash, z0, z1, &B);
        DMatrix(Qdash, z0, z1, &D);
    }
    file1 << "            A matrix" << std::endl;
    for (int j = 0; j < 3; j++)
    {
        file1 << "                ";
        for (int k = 0; k < 3; k++)
        {
            file1 << A(j, k) << "  ";
        }
        file1 << std::endl;
    }
    file1 << "            B matrix" << std::endl;
    for (int j = 0; j < 3; j++)
    {
        file1 << "                ";
        for (int k = 0; k < 3; k++)
        {
            file1 << B(j, k) << "  ";
        }
        file1 << std::endl;
    }
    file1 << "            D matrix" << std::endl;
    for (int j = 0; j < 3; j++)
    {
        file1 << "                ";
        for (int k = 0; k < 3; k++)
        {
            file1 << D(j, k) << "  ";
        }
        file1 << std::endl;
    }

    V.Amat = Eigen::MatrixXd::Zero(3, 3);
    V.Bmat = Eigen::MatrixXd::Zero(3, 3);
    V.Dmat = Eigen::MatrixXd::Zero(3, 3);

    V.Amat = A;
    V.Bmat = B;
    V.Dmat = D;

    //Variables required for analytical cross sectional stiffness matrix
    double A11 = A(0, 0);
    double A12 = A(0, 1);
    double A16 = A(0, 2);
    double A22 = A(1, 1);
    double A26 = A(1, 2);
    double A66 = A(2, 2);

    double B11 = B(0, 0);
    double B12 = B(0, 1);
    double B16 = B(0, 2);
    double B22 = B(1, 1);
    double B26 = B(1, 2);
    double B66 = B(2, 2);

    double D11 = D(0, 0);
    double D12 = D(0, 1);
    double D16 = D(0, 2);
    double D22 = D(1, 1);
    double D26 = D(1, 2);
    double D66 = D(2, 2);

    V.B.A11 = A11 + (pow(A16, 2) * A22 - 2 * A12 * A16 * A26 + pow(A12, 2) * A66) / (pow(A26, 2) - A22 * A66);
    //std::cout << V.B.A11 << std::endl;
    V.B.B11 = B11 + (A12 * A66 * B12 + A16 * A22 * B16 - A26 * (A16 * B12 + A12 * B16)) / (pow(A26, 2) - A22 * A66);
    //std::cout << V.B.B11 << std::endl;
    V.B.B12 = B12 + (A12 * A66 * B22 + A16 * A22 * B26 - A26 * (A16 * B22 + A12 * B26)) / (pow(A26, 2) - A22 * A66);
    //std::cout << V.B.B12 << std::endl;
    V.B.B16 = B16 + (A12 * A66 * B26 + A16 * A22 * B66 - A26 * (A16 * B26 + A12 * B66)) / (pow(A26, 2) - A22 * A66);
    //std::cout << V.B.B16 << std::endl;

    V.B.D11 = D11 + (A66 * pow(B12, 2) - 2 * A26 * B12 * B16 + A22 * pow(B16, 2)) / (pow(A26, 2) - A22 * A66);
    V.B.D12 = D12 + (A66 * B12 * B22 + A22 * B16 * B26 - A26 * (B16 * B22 + B12 * B26)) / (pow(A26, 2) - A22 * A66);
    V.B.D22 = D22 + (A66 * pow(B22, 2) - 2 * A26 * B22 * B26 + A22 * pow(B26, 2)) / (pow(A26, 2) - A22 * A66);
    V.B.D16 = D16 + (A66 * B12 * B26 + A22 * B16 * B66 - A26 * (B16 * B26 + B12 * B66)) / (pow(A26, 2) - A22 * A66);
    V.B.D26 = D26 + (A66 * B22 * B26 + A22 * B26 * B66 - A26 * (pow(B26, 2) + B22 * B66)) / (pow(A26, 2) - A22 * A66);
    V.B.D66 = D66 + (A66 * pow(B26, 2) - 2 * A26 * B26 * B66 + A22 * pow(B66, 2)) / (pow(A26, 2) - A22 * A66);

    V.BB.A11 = V.B.A11 - pow(V.B.B12, 2) / V.B.D22;
    V.BB.B11 = V.B.B11 - V.B.B12 * V.B.D12 / V.B.D22;
    V.BB.B16 = V.B.B16 - V.B.B12 * V.B.D26 / V.B.D22;
    V.BB.D11 = V.B.D11 - pow(V.B.D12, 2) / V.B.D22;
    V.BB.D16 = V.B.D16 - V.B.D12 * V.B.D26 / V.B.D22;
    V.BB.D66 = V.B.D66 - pow(V.B.D26, 2) / V.B.D22;

    /*file1 << "            Values of variables from appendix" << std::endl;
    file1 << "                A11:  " << V.B.A11 << std::endl;
    file1 << "                B11:  " << V.B.B11 << std::endl;
    file1 << "                B12:  " << V.B.B12 << std::endl;
    file1 << "                B16:  " << V.B.B16 << std::endl;
    file1 << "                D11:  " << V.B.D11 << std::endl;
    file1 << "                D12:  " << V.B.D12 << std::endl;
    file1 << "                D22:  " << V.B.D22 << std::endl;
    file1 << "                D16:  " << V.B.D16 << std::endl;
    file1 << "                D26:  " << V.B.D26 << std::endl;
    file1 << "                D66:  " << V.B.D66 << std::endl;
    file1 << "            Values of variables with double bar from appendix" << std::endl;
    file1 << "                A11:  " << V.BB.A11 << std::endl;
    file1 << "                B11:  " << V.BB.B11 << std::endl;
    file1 << "                B16:  " << V.BB.B16 << std::endl;
    file1 << "                D11:  " << V.BB.D11 << std::endl;
    file1 << "                D16:  " << V.BB.D16 << std::endl;
    file1 << "                D66:  " << V.BB.D66 << std::endl;*/

    //Eigen::MatrixXd Seq = Eigen::MatrixXd::Zero(4, 4);
    /*S.Sl = Eigen::MatrixXd::Zero(4, 4);
    S.Sln = Eigen::MatrixXd::Zero(4, 5);
    S.Sn = Eigen::MatrixXd::Zero(5, 5);

    S.Sl(0, 0) = b * V.BB.A11;
    S.Sl(0, 1) = (-2) * b * V.BB.B16 + pow(b, 3) * V.BB.A11 * k1 / 12;
    S.Sl(0, 2) = b * V.BB.B11;
    S.Sl(1, 0) = S.Sl(0, 1);
    S.Sl(1, 1) = 4 * b * V.BB.D66 - pow(b, 3) * V.BB.B16 * k1 / 3 + pow(b, 5) * V.BB.A11 * pow(k1, 2) / 80;
    S.Sl(1, 2) = -2 * b * V.BB.D16 + pow(b, 3) * V.BB.B11 * k1 / 12;
    S.Sl(2, 0) = S.Sl(0, 2);
    S.Sl(2, 1) = S.Sl(1, 2);
    S.Sl(2, 2) = b * V.BB.D11;
    S.Sl(3, 3) = pow(b, 3) * V.BB.A11 / 12;

    S.Sln(0, 0) = pow(b, 3) * V.BB.A11 / 24;
    S.Sln(1, 0) = -pow(b, 3) * V.BB.B16 / 12 + pow(b, 5) * V.BB.A11 * k1 / 160;
    S.Sln(1, 1) = pow(b, 5) * V.BB.A11 * V.B.D12 * k1 / (360 * V.B.D22);
    S.Sln(1, 2) = pow(b, 5) * V.BB.A11 * V.B.B12 * k1 / (360 * V.B.D22);
    S.Sln(2, 0) = (pow(b, 3) * V.BB.B11 / 24 + pow(b, 5) * V.BB.A11 * V.B.D26 * k1 / (180 * V.B.D22)
        + pow(b, 7) * V.BB.A11 * V.B.B12 * pow(k1, 2) / (10080 * V.B.D22));
    S.Sln(3, 3) = -pow(b, 5) * V.BB.A11 * V.B.B12 / (720 * V.B.D22);

    S.Sn(0, 0) = pow(b, 5) * V.BB.A11 / 320;
    S.Sn(0, 2) = pow(b, 5) * V.BB.A11 * V.B.B12 / (720 * V.B.D22);
    S.Sn(0, 4) = -pow(b, 5) * V.BB.A11 * V.B.D26 / (360 * V.B.D22) + pow(b, 7) * V.BB.A11 * V.B.B12 * k1 / (10080 * V.B.D22);
    S.Sn(1, 1) = pow(b, 5) * V.BB.A11 * pow(V.B.D12, 2) / (720 * pow(V.B.D22, 2));
    S.Sn(1, 2) = pow(b, 5) * V.BB.A11 * V.B.B12 * V.B.D12 / (720 * pow(V.B.D22, 2));
    S.Sn(1, 4) = -pow(b, 5) * V.BB.A11 * V.B.D12 * V.B.D26 / (360 * pow(V.B.D22, 2))
        - pow(b, 7) * V.BB.A11 * V.B.B12 * V.B.D12 * k1 / (60480 * pow(V.B.D22, 2));
    S.Sn(2, 2) = pow(b, 5) * V.BB.A11 * pow(V.B.B12, 2) / (720 * pow(V.B.D22, 2));
    S.Sn(2, 4) = -pow(b, 5) * V.BB.A11 * V.B.B12 * V.B.D26 / (360 * pow(V.B.D22, 2))
        - pow(b, 7) * V.BB.A11 * pow(V.B.B12, 2) * k1 / (60480 * pow(V.B.D22, 2));
    S.Sn(3, 3) = (pow(b, 7) * V.BB.A11 * pow(V.B.B12, 2) / (10080 * pow(V.B.D22, 2)) - pow(b, 7) * pow(V.BB.A11, 2) / (30240 * V.B.D22));
    S.Sn(4, 4) = (pow(b, 5) * V.BB.A11 * pow(V.B.D26, 2) / (180 * pow(V.B.D22, 2)) + pow(b, 5) * V.BB.A11 * V.B.D12 / (360 * V.B.D22)
        + pow(b, 7) * V.BB.A11 * V.B.B12 * V.B.D26 * k1 / (15120 * pow(V.B.D22, 2))
        - (pow(b, 9) * pow(V.BB.A11, 2) / (90720 * V.B.D22) + pow(b, 9) * V.BB.A11 * pow(V.B.B12, 2) / (403200 * pow(V.B.D22, 2))) * pow(k1, 2));

//    std::cout<<S<<std::endl;

    return V;
}*/

Eigen::MatrixXd Equivalent_StiffnessMatrix_ZerothOrder()
{
    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(6, 6);

    S(0, 0) = 0.16330186e7;
    S(0, 1) = 0.82478352e-7;
    S(0, 2) = -0.34645950e-8;
    S(0, 3) = -0.34119601e6;
    S(0, 4) = -0.26268760e-7;
    S(0, 5) = 0.10627862e-6;

    S(1, 0) = 0.82478352e-7;
    S(1, 1) = 0.24102720e6;
    S(1, 2) = -0.58356855e-8;
    S(1, 3) = 0.63811455e-7;
    S(1, 4) = 0.12876931e6;
    S(1, 5) = -0.60489689e-6;

    S(2, 0) = -0.34645950e-8;
    S(2, 1) = -0.58184268e-8;
    S(2, 2) = 0.75401574e3;
    S(2, 3) = -0.27501371e-06;
    S(2, 4) = -0.29004439E-08;
    S(2, 5) = -0.28277805E+05;

    S(3, 0) = -0.34119601E+06;
    S(3, 1) = 0.63811455E-07;
    S(3, 2) = -0.27501371E-06;
    S(3, 3) = 0.18988820E+06;
    S(3, 4) = -0.26554143E-07;
    S(3, 5) = 0.10274357E-04;

    S(4, 0) = -0.41817119E-07;
    S(4, 1) = 0.12876931E+06;
    S(4, 2) = -0.29004439E-08;
    S(4, 3) = -0.28599112E-07;
    S(4, 4) = 0.21183083E+06;
    S(4, 5) = 0.62027751E-07;

    S(5, 0) = 0.10673869E-06;
    S(5, 1) = -0.60489689E-06;
    S(5, 2) = -0.28277805E+05;
    S(5, 3) = 0.10274757E-04;
    S(5, 4) = 0.93949918E-07;
    S(5, 5) = 0.85205269E+08;


    return S;
}

Eigen::MatrixXd Equivalent_StiffnessMatrix_FirstOrder(Eigen::VectorXd Strain, Eigen::VectorXd inittwist, double b, std::fstream& file1)
{
    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(6, 6);

    //Strains
    double gamma11, kappa1, kappa2, kappa3;
    gamma11 = Strain(0);
    kappa1 = Strain(3);
    kappa2 = Strain(4);
    kappa3 = Strain(5);

    //Initial twist and bend
    double k1, k2;
    k1 = inittwist(0);
    k2 = inittwist(1);

    /* S(0, 0) = (1.0 / 720) * b * V.BB.A11 * (720 * pow(V.B.D22, 2) + pow(kappa2, 2) * pow(b, 4) * pow(V.B.B12, 2) + pow(k2, 2) * pow(b, 4) * pow(V.B.B12, 2) +
         2 * k2 * pow(b, 4) * pow(V.B.B12, 2) * kappa2) / pow(V.B.D22, 2);

     S(0, 3) = -(1.0 / 60480) * b * (120960 * V.BB.B16 * pow(V.B.D22, 2) - 5040 * pow(b, 2) * V.BB.A11 * k1 * pow(V.B.D22, 2) - 5040 * pow(b, 2) * V.BB.A11 * kappa1 * pow(V.B.D22, 2)
         - 168 * pow(b, 4) * V.BB.A11 * V.B.B12 * k1 * kappa2 * V.B.D22 - 168 * pow(b, 4) * V.BB.A11 * V.B.B12 * kappa1 * kappa2 * V.B.D22 + 168 * pow(b, 4) * V.BB.A11 * V.B.B12 * pow(kappa2, 2) * V.B.D26
         + pow(kappa2, 2) * pow(b, 6) * V.BB.A11 * pow(V.B.B12, 2) * k1 - 168 * k2 * pow(b, 4) * V.BB.A11 * V.B.B12 * k1 * V.B.D22 + 4 * pow(k2, 2) * pow(b, 6) * pow(V.B.B12, 2) * k1 * V.B.A11 + 168 * pow(k2, 2) * pow(b, 4) * V.B.B12 * V.BB.A11 * V.B.D26
         - 168 * k2 * pow(b, 4) * V.B.B12 * V.BB.A11 * kappa1 * V.B.D22 + 8 * k2 * pow(b, 4) * pow(V.B.B12, 2) * kappa2 * k1 * V.B.A11 + 336 * k2 * pow(b, 4) * V.B.B12 * kappa2 * V.BB.A11 * V.B.D26) / pow(V.B.D22, 2);

     S(0, 4) = -(1.0 / 30240) * b * (-30240 * V.BB.B11 * pow(V.B.D22, 2) - 84 * pow(b, 4) * V.BB.A11 * V.B.B12 * k1 * kappa1 * V.B.D22 - 42 * pow(b, 4) * V.BB.A11 * V.B.B12 * pow(kappa1, 2) * V.B.D22
         - 126 * pow(b, 4) * V.BB.A11 * V.B.B12 * V.B.D12 * pow(kappa2, 2) - 84 * pow(b, 4) * V.BB.A11 * pow(V.B.B12,2) * kappa2 * gamma11 + 168 * pow(b, 4) * V.BB.A11 * V.B.B12 * kappa2 * kappa1 * V.B.D26 + pow(b, 6) * V.BB.A11 * pow(V.B.B12, 2) * kappa2 * kappa1 * k1
         - 42 * pow(k2, 2) * pow(b, 4) * V.BB.A11 * V.B.B12 * V.B.D12 - 84 * k2 * pow(b, 4) * pow(V.B.B12, 2) * V.BB.A11 * gamma11 + 4 * k2 * pow(b, 4) * pow(V.B.B12, 2) * kappa1 * k1 * V.B.A11 + 168 * k2 * pow(b, 4) * V.B.B12 * kappa1 * V.BB.A11 * V.B.D26
         - 168 * k2 * pow(b, 4) * V.B.B12 * V.BB.A11 * V.B.D12 * kappa2) / pow(V.B.D22, 2);

     S(1, 1) = 1e20;

     S(2, 2) = 1e20;

     S(3, 0) = -(1.0 / 60480) * b * (120960 * V.BB.B16 * pow(V.B.D22, 2) - 5040 * pow(b, 2) * V.BB.A11 * k1 * pow(V.B.D22, 2) - 5040 * pow(b, 2) * V.BB.A11 * kappa1 * pow(V.B.D22, 2) - 168 * pow(b, 4) * V.BB.A11 * V.B.B12 * k1 * kappa2 * V.B.D22
         - 168 * kappa1 * pow(b, 4) * V.BB.A11 * V.B.B12 * kappa2 * V.B.D22 + 168 * pow(kappa2, 2) * pow(b, 4) * V.BB.A11 * V.B.B12 * V.B.D26 + pow(kappa2, 2) * pow(b, 6) * V.BB.A11 * pow(V.B.B12, 2) * k1 - 168 * k2 * pow(b, 4) * V.BB.A11 * V.B.B12 * k1 * V.B.D22
         + 4 * pow(k2, 2) * pow(b, 6) * pow(V.B.B12, 2) * k1 * V.B.A11 + 168 * pow(k2, 2) * pow(b, 4) * V.B.B12 * V.BB.A11 * V.B.D26 - 168 * k2 * pow(b, 4) * V.B.B12 * V.BB.A11 * kappa1 * V.B.D22 + 8 * k2 * pow(b, 4) * pow(V.B.B12, 2) * kappa2 * k1 * V.B.A11
         + 336 * k2 * pow(b, 4) * V.B.B12 * kappa2 * V.BB.A11 * V.B.D26) / pow(V.B.D22, 2);

     S(3, 3) = -(1.0 / 3628800) * b * (-10080 * k2 * pow(b, 4) * V.BB.A11 * V.B.B12 * gamma11 * V.B.D22 + 40320 * kappa2 * pow(b, 4) * V.B.D22 * V.BB.A11 * k1 * V.B.D26 - 14515200 * V.BB.D66 * pow(V.B.D22, 2)
         - 960 * pow(k2, 2) * pow(b, 6) * k1 * V.B.A11 * V.B.D26 * V.B.B12 + 1209600 * pow(b, 2) * V.BB.B16 * k1 * pow(V.B.D22, 2) - 45360 * pow(b, 4) * V.BB.A11 * pow(k1, 2) * pow(V.B.D22, 2) - 302400 * gamma11 * pow(b, 2) * V.BB.A11 * pow(V.B.D22, 2)
         - 68040 * pow(b, 4) * V.BB.A11 * pow(kappa1, 2) * pow(V.B.D22, 2) + 1814400 * pow(b, 2) * kappa1 * pow(V.B.D22, 2) * V.BB.B16 - 302400 * kappa2 * pow(b, 2) * pow(V.B.D22, 2) * V.BB.B11 - 20160 * pow(kappa2, 2) * pow(b, 4) * V.BB.A11 * pow(V.B.D26, 2)
         - 20160 * pow(k2, 2) * pow(b, 4) * V.BB.A11 * pow(V.B.D26, 2) - 10080 * pow(b, 4) * V.BB.A11 * V.B.D12 * kappa2 * gamma11 * V.B.D22 - 136080 * pow(b, 4) * kappa1 * pow(V.B.D22, 2) * V.BB.A11 * k1
         - 720 * kappa2 * pow(b, 6) * V.B.D22 * V.BB.A11 * V.B.B12 * pow(k1, 2) - 2160 * pow(b, 6) * V.BB.A11 * kappa2 * kappa1 * V.B.D22 * V.B.B12 * k1 + 60480 * pow(b, 4) * V.BB.A11 * kappa2 * kappa1 * V.B.D22 * V.B.D26
         - 10080 * pow(kappa2, 2) * pow(b, 4) * V.BB.A11 * V.B.D12 * V.B.D22 - 240 * pow(kappa2, 2) * pow(b, 6) * V.BB.A11 * V.B.B12 * V.B.D26 * k1 + 40 * pow(kappa2, 2) * pow(b, 8) * pow(V.BB.A11, 2) * pow(k1, 2) * V.B.D22 + 9 * pow(kappa2, 2) * pow(b, 8) * V.BB.A11 * pow(k1, 2) * pow(V.B.B12, 2)
         - 720 * k2 * pow(b, 6) * pow(k1, 2) * V.B.D22 * V.B.A11 * V.B.B12 + 60480 * k2 * pow(b, 4) * kappa1 * V.B.D22 * V.BB.A11 * V.B.D26 - 10080 * k2 * pow(b, 4) * kappa2 * V.BB.A11 * V.B.D12 * V.B.D22
         - 1080 * k2 * pow(b, 6) * kappa1 * V.B.D22 * k1 * V.B.A11 * V.B.B12 + 60480 * k2 * pow(b, 4) * kappa1 * V.B.D22 * V.BB.A11 * V.B.D26 - 10080 * k2 * pow(b, 4) * kappa2 * V.BB.A11 * V.B.D12 * V.B.D22
         - 40320 * k2 * pow(b, 4) * kappa2 * V.BB.A11 * pow(V.B.D26, 2) - 1920 * k2 * pow(b, 6) * kappa2 * k1 * V.B.A11 * V.B.D26 * V.B.B12 + 80 * k2 * pow(b, 8) * kappa2 * pow(k1, 2) * pow(V.B.A11, 2) * V.B.D22) / pow(V.B.D22, 2);

     S(3, 4) = -(1.0 / 1814400) * b * (20160 * pow(b, 4) * kappa1 * V.B.D22 * V.BB.A11 * k1 * V.B.D26 - 151200 * pow(b, 2) * V.BB.B11 * k1 * pow(V.B.D22, 2) - 151200 * pow(b, 2) * kappa1 * pow(V.B.D22, 2) * V.BB.B11
         + 3628800 * V.BB.D16 * pow(V.B.D22, 2) - 10080 * pow(b, 4) * V.BB.A11 * V.B.D12 * k1 * kappa2 * V.B.D22 - 5040 * pow(b, 4) * V.BB.A11 * V.B.B12 * k1 * gamma11 * V.B.D22
         - 5040 * k2 * pow(b, 4) * V.BB.A11 * V.B.D12 * k1 * V.B.D22 - 360 * pow(b, 6) * kappa1 * V.B.D22 * V.BB.A11 * V.B.B12 * pow(k1, 2) - 5040 * kappa1 * pow(b, 4) * V.BB.A11 * V.B.D22 * V.B.B12 * gamma11
         + 15120 * pow(kappa1, 2) * pow(b, 4) * V.BB.A11 * V.B.D22 * V.B.D26 - 540 * pow(kappa1, 2) * pow(b, 6) * V.BB.A11 * V.B.D22 * V.B.B12 * k1 + 90 * pow(b, 6) * V.BB.A11 * V.B.D12 * pow(kappa2, 2) * V.B.B12 * k1
         + 15120 * pow(b, 4) * V.BB.A11 * V.B.D12 * pow(kappa2, 2) * V.B.D26 + 10080 * gamma11 * pow(b, 4) * V.BB.A11 * V.B.B12 * kappa2 * V.B.D26 + 60 * gamma11 * pow(b, 6) * V.BB.A11 * pow(V.B.B12, 2) * kappa2 * k1
         - 20160 * pow(b, 4) * V.BB.A11 * kappa2 * kappa1 * pow(V.B.D26, 2) - 10080 * gamma11 * pow(b, 4) * V.BB.A11 * kappa2 * kappa1 * V.B.D12 * V.B.D22 - 240 * pow(b, 6) * V.BB.A11 * kappa1 * kappa2 * V.B.B12 * V.B.D26 * k1
         + 40 * pow(b, 8) * pow(V.BB.A11, 2) * kappa2 * kappa1 * pow(k1, 2) * V.B.D22 + 9 * pow(b, 8) * V.BB.A11 * kappa2 * kappa1 * pow(k1, 2) * pow(V.B.B12, 2) + 120 * pow(k2, 2) * pow(b, 6) * V.B.D12 * k1 * V.B.A11 * V.B.B12 + 5040 * pow(k2, 2) * pow(b, 4) * V.B.D12 * V.BB.A11 * V.B.D26
         + 480 * k2 * pow(b, 6) * V.B.D12 * kappa2 * k1 * V.B.A11 * V.B.B12 + 20160 * k2 * pow(b, 4) * V.B.D12 * kappa2 * V.BB.A11 * V.B.D26 + 240 * k2 * pow(b, 4) * pow(V.B.B12, 2) * gamma11 * k1 * V.B.A11
         + 10080 * k2 * pow(b, 4) * V.B.B12 * gamma11 * V.BB.A11 * V.B.D26 - 5040 * k2 * pow(b, 4) * kappa1 * V.BB.A11 * V.B.D12 * V.B.D22 - 20160 * k2 * pow(b, 4) * kappa1 * V.BB.A11 * pow(V.B.D26, 2)
         - 960 * k2 * pow(b, 6) * kappa1 * k1 * V.B.A11 * V.B.D26 * V.B.B12 + 40 * k2 * pow(b, 8) * kappa1 * pow(k1, 2) * pow(V.B.A11, 2) * V.B.D22 - 180 * k2 * pow(b, 6) * V.BB.A11 * kappa3 * pow(V.B.B12, 2)
         + 60 * k2 * pow(b, 6) * pow(V.BB.A11, 2) * kappa3 * V.B.D22) / pow(V.B.D22, 2);

     S(3, 5) = -(1.0 / 30240) * k2 * pow(b, 7) * V.BB.A11 * (-3 * pow(V.B.B12, 2) + V.BB.A11 * V.B.D22) * kappa2 / pow(V.B.D22, 2);

     S(4, 0) = -(1.0 / 30240) * b * (-30240 * V.BB.B11 * pow(V.B.D22, 2) - 84 * pow(b, 4) * V.BB.A11 * V.B.B12 * k1 * kappa1 * V.B.D22 - 42 * pow(b, 4) * V.BB.A11 * V.B.B12 * pow(kappa1, 2) * V.B.D22 - 126 * pow(b, 4) * V.BB.A11 * V.B.B12 * V.B.D12 * pow(kappa2, 2)
         - 84 * pow(b, 4) * V.BB.A11 * pow(V.B.B12, 2) * kappa2 * gamma11 + 168 * pow(b, 4) * V.BB.A11 * V.B.B12 * kappa2 * kappa1 * V.B.D26 + pow(b, 6) * V.BB.A11 * pow(V.B.B12, 2) * kappa2 * kappa1 * k1 - 42 * pow(k2, 2) * pow(b, 4) * V.BB.A11 * V.B.B12 * V.B.D12
         - 84 * k2 * pow(b, 4) * pow(V.B.B12, 2) * V.BB.A11 * gamma11 + 4 * k2 * pow(b, 4) * pow(V.B.B12, 2) * kappa1 * k1 * V.B.A11 + 168 * k2 * pow(b, 4) * V.B.B12 * kappa1 * V.BB.A11 * V.B.D26
         - 168 * k2 * pow(b, 4) * V.B.B12 * V.BB.A11 * V.B.D12 * kappa2) / pow(V.B.D22, 2);

     S(4, 3) = -(1.0 / 1814400) * b * (20160 * pow(b, 4) * kappa1 * V.B.D22 * V.BB.A11 * k1 * V.B.D26 - 151200 * pow(b, 2) * V.BB.B11 * k1 * pow(V.B.D22, 2) - 151200 * pow(b, 2) * kappa1 * pow(V.B.D22, 2) * V.BB.B11
         + 3628800 * V.BB.D16 * pow(V.B.D22, 2) - 10080 * pow(b, 4) * V.BB.A11 * V.B.D12 * k1 * kappa2 * V.B.D22 - 5040 * pow(b, 4) * V.BB.A11 * V.B.B12 * k1 * gamma11 * V.B.D22
         - 5040 * k2 * pow(b, 4) * V.BB.A11 * V.B.D12 * k1 * V.B.D22 - 360 * pow(b, 6) * kappa1 * V.B.D22 * V.BB.A11 * V.B.B12 * pow(k1, 2) - 5040 * kappa1 * pow(b, 4) * V.BB.A11 * V.B.D22 * V.B.B12 * gamma11
         + 15120 * pow(kappa1, 2) * pow(b, 4) * V.BB.A11 * V.B.D22 * V.B.D26 - 540 * pow(kappa1, 2) * pow(b, 6) * V.BB.A11 * V.B.D22 * V.B.B12 * k1 + 90 * pow(b, 6) * V.BB.A11 * V.B.D12 * pow(kappa2, 2) * V.B.B12 * kappa1
         + 15120 * pow(b, 4) * V.BB.A11 * V.B.D12 * pow(kappa2, 2) * V.B.D26 + 10080 * gamma11 * pow(b, 4) * V.BB.A11 * V.B.B12 * kappa2 * V.B.D26 + 60 * gamma11 * pow(b, 6) * V.BB.A11 * pow(V.B.B12, 2) * kappa2 * k1
         - 20160 * pow(b, 4) * V.BB.A11 * kappa2 * kappa1 * pow(V.B.D26, 2) - 10080 * pow(b, 4) * V.BB.A11 * kappa2 * kappa1 * V.B.D12 * V.B.D22 - 240 * pow(b, 6) * V.BB.A11 * kappa2 * kappa1 * V.B.B12 * V.B.D26 * k1
         + 40 * pow(b, 8) * pow(V.BB.A11, 2) * kappa2 * kappa1 * pow(k1, 2) * V.B.D22 + 9 * pow(b, 8) * V.BB.A11 * kappa2 * kappa1 * pow(k1, 2) * pow(V.B.B12, 2) + 120 * pow(k2, 2) * pow(b, 6) * V.B.D12 * k1 * V.B.A11 * V.B.B12 + 5040 * pow(k2, 2) * pow(b, 4) * V.B.D12 * V.BB.A11 * V.B.D26
         + 480 * k2 * pow(b, 6) * V.B.D12 * kappa2 * k1 * V.B.A11 * V.B.B12 + 20160 * k2 * pow(b, 4) * V.B.D12 * kappa2 * V.BB.A11 * V.B.D26 + 240 * k2 * pow(b, 4) * pow(V.B.B12, 2) * gamma11 * k1 * V.B.A11
         + 10080 * k2 * pow(b, 4) * V.B.B12 * gamma11 * V.BB.A11 * V.B.D26 - 5040 * k2 * pow(b, 4) * kappa1 * V.BB.A11 * V.B.D12 * V.B.D22 - 20160 * k2 * pow(b, 4) * kappa1 * V.BB.A11 * pow(V.B.D26, 2)
         - 960 * k2 * pow(b, 6) * kappa1 * k1 * V.B.A11 * V.B.D26 * V.B.B12 + 40 * k2 * pow(b, 8) * kappa1 * pow(k1, 2) * pow(V.B.A11, 2) * V.B.D22 - 180 * k2 * pow(b, 6) * V.BB.A11 * kappa3 * pow(V.B.B12, 2)
         + 60 * k2 * pow(b, 6) * pow(V.BB.A11, 2) * kappa3 * V.B.D22) / pow(V.B.D22, 2);

     S(4, 4) = -(1.0 / 3628800) * b * (60480 * pow(b, 4) * V.BB.A11 * V.B.D12 * kappa2 * kappa1 * V.B.D26 - 30240 * pow(b, 4) * V.BB.A11 * pow(V.B.D12, 2) * pow(kappa2, 2) - 5040 * pow(k2, 2) * pow(b, 4) * V.BB.A11 * pow(V.B.D12, 2)
         - 5040 * pow(gamma11, 2) * pow(b, 4) * V.BB.A11 * pow(V.B.B12, 2) - 360 * pow(kappa3, 2) * pow(b, 6) * V.BB.A11 * pow(V.B.B12, 2) + 120 * pow(kappa3, 2) * pow(b, 6) * pow(V.BB.A11, 2) * V.B.D22 - 20160 * pow(kappa1, 2) * pow(b, 4) * V.BB.A11 * pow(V.B.D26, 2)
         - 3628800 * V.BB.D11 * pow(V.B.D22, 2) - 20160 * pow(b, 4) * V.BB.A11 * V.B.D12 * k1 * kappa1 * V.B.D22 - 30240 * pow(b, 4) * V.BB.A11 * V.B.B12 * V.B.D12 * kappa2 * gamma11
         + 360 * pow(b, 6) * V.BB.A11 * V.B.D12 * kappa2 * kappa1 * V.B.B12 * k1 + 20160 * gamma11 * pow(b, 4) * V.BB.A11 * V.B.B12 * kappa1 * V.B.D26 + 120 * gamma11 * pow(b, 6) * V.BB.A11 * pow(V.B.B12, 2) * kappa1 * k1
         - 10080 * pow(kappa1, 2) * pow(b, 4) * V.BB.A11 * V.B.D12 * V.B.D22 - 240 * pow(kappa1, 2) * pow(b, 6) * V.BB.A11 * V.B.B12 * V.B.D26 * k1 + 40 * pow(kappa1, 2) * pow(b, 8) * pow(V.BB.A11, 2) * pow(k1, 2) * V.B.D22 + 9 * pow(kappa1, 2) * pow(b, 8) * V.BB.A11 * pow(k1, 2) * pow(V.B.B12, 2)
         + 960 * k2 * pow(b, 6) * V.B.D12 * kappa1 * k1 * V.B.A11 * V.B.B12 + 40320 * k2 * pow(b, 4) * V.B.D12 * kappa1 * V.BB.A11 * V.B.D26 - 30240 * k2 * pow(b, 4) * pow(V.B.D12, 2) * V.BB.A11 * kappa2
         - 20160 * k2 * pow(b, 4) * V.B.D12 * V.BB.A11 * V.B.B12 * gamma11) / pow(V.B.D22, 2);

     S(4, 5) = -(1.0 / 30240) * pow(b, 5) * V.BB.A11 * (84 * V.B.B12 * kappa3 * V.B.D22 - 6 * pow(b, 2) * kappa2 * kappa3 * pow(V.B.B12, 2) + 2 * pow(b, 2) * kappa2 * kappa3 * V.BB.A11 * V.B.D22 - 3 * k2 * pow(b, 2) * kappa1 * pow(V.B.B12, 2)
         + k2 * pow(b, 2) * kappa1 * V.BB.A11 * V.B.D22) / pow(V.B.D22, 2);

     S(5, 3) = -(1.0 / 30240) * k2 * pow(b, 7) * V.BB.A11 * (-3 * pow(V.B.B12, 2) + V.BB.A11 * V.B.D22) * kappa2 / pow(V.B.D22, 2);

     S(5, 4) = -(1.0 / 30240) * pow(b, 5) * V.BB.A11 * (84 * V.B.B12 * kappa3 * V.B.D22 - 6 * pow(b, 2) * kappa2 * kappa3 * pow(V.B.B12, 2) + 2 * pow(b, 2) * kappa2 * kappa3 * V.BB.A11 * V.B.D22 - 3 * k2 * pow(b, 2) * kappa1 * pow(V.B.B12, 2)
         + k2 * pow(b, 2) * kappa1 * V.BB.A11 * V.B.D22) / pow(V.B.D22, 2);

     S(5, 5) = -(1.0 / 30240) * pow(b, 3) * V.BB.A11 * (-2520 * pow(V.B.D22, 2) + 84 * pow(b, 2) * V.B.B12 * kappa2 * V.B.D22 - 3 * pow(kappa2, 2) * pow(b, 4) * pow(V.B.B12, 2) + pow(kappa2, 2) * pow(b, 4) * V.BB.A11 * V.B.D22 + 84 * k2 * pow(b, 2) * V.B.B12 * V.B.D22
         - 3 * pow(k2, 2) * pow(b, 4) * pow(V.B.B12, 2) + pow(k2, 2) * pow(b, 4) * V.BB.A11 * V.B.D22) / pow(V.B.D22, 2);*/

    // Antisymmetric Laminate Cross-Sectional Stiffness Matrix Values
    /*S(0, 0) = 1.67158e6;
    S(0, 3) = -332.582 + 89.8695 * kappa1;

    S(1, 1) = 1e20;

    S(2, 2) = 1e20;

    S(3, 0) = S(0, 3);
    S(3, 3) = 0.19022 + 89.8695 * gamma11 - 0.0536422 * kappa1 + 0.0130455 * pow(kappa1, 2) + 0.00122216 * pow(kappa2, 2);
    S(3, 4) = 0.00244432 * kappa1 * kappa2;

    S(4, 3) = S(3, 4);
    S(4, 4) = 0.144537 + 0.00122216 * pow(kappa1, 2) + 0.00231857 * pow(kappa2, 2) - 0.178883 * pow(kappa3, 2);
    S(4, 5) = -0.357766 * kappa2 * kappa3;

    S(5, 4) = S(4, 5);
    S(5, 5) = 89.8695 - 0.178883 * pow(kappa2, 2);*/

    //Cross-Sectional Stiffness Matrix for CAS-1 from Wenbin Yu
    /*S(0, 0) = 0.137e7;
    S(0, 1) = -0.184e6;
    S(0, 2) = -0.150e3;
    S(1, 1) = 0.885e5;
    S(1, 2) = 0.803e2;
    S(2, 2) = 0.387e5;
    S(3, 3) = 0.170e5;
    S(3, 4) = 0.176e5;
    S(3, 5) = -0.349e3;
    S(4, 4) = 0.591e5;
    S(4, 5) = -0.371e3;
    S(5, 5) = 0.141e6;*/

    //Cross Sectional Stiffness Matrix for BT beam from Cesnik
    /*S(0, 0) = 0.8115 * pow(10, 6);
    S(0, 1) = -0.4655 * pow(10, 5);
    S(1, 0) = S(0, 1);
    S(1, 1) = 0.9368 * pow(10, 5);
    S(2, 2) = 0.6882 * pow(10, 4);
    S(3, 3) = 0.1251 * pow(10, 3);
    S(3, 4) = 0.3455 * pow(10, 2);
    S(4, 3) = S(3, 4);
    S(4, 4) = 0.1852 * pow(10, 3);
    S(5, 5) = 0.9178 * pow(10, 5);*/

    //4*4 matrix from cesnik
    /*S(0, 0) = 0.3607 * pow(10, 7);
    S(1, 1) = 1e7;
    S(2, 2) = 1e7;
    S(3, 3) = 0.4670 * pow(10, 0);
    S(3, 4) = 0.9864 * pow(10, -1);
    S(4, 3) = S(3, 4);
    S(4, 4) = 0.5297 * pow(10, 0);
    S(5, 5) = 0.2628 * pow(10, 4);*/

    //From Xia Xaoyang thesis
    S(0, 0) = 0.7884e6;
    S(1, 1) = 1e7;
    S(2, 2) = 1e7;
    S(3, 3) = 0.1290e3;
    S(3, 4) = 0.3653e2;
    S(4, 3) = S(3, 4);
    S(4, 4) = 0.1864e3;
    S(5, 5) = 0.9179e5;

    // Symmetric Matrix Cross-Sectional Stiffness Matrix Values
    /*S(0, 0) = 1.13132e7;
    S(0, 3) = 18478.2 * kappa1;

    S(1, 1) = 1e20;

    S(2, 2) = 1e20;

    S(3, 0) = S(0, 3);
    S(3, 3) = 0.311569 + 18478.2 * gamma11 + 81.4889 * pow(kappa1, 2) + 2.75785 * pow(kappa2, 2);
    S(3, 4) = 5.5157 * kappa1 * kappa2;

    S(4, 0) = 0;
    S(4, 3) = S(3, 4);
    S(4, 4) = 1.41032 + 2.75785 * pow(kappa1, 2) + 1.89003 * pow(kappa2, 2) - 279387 * pow(kappa3, 2);
    S(4, 5) = -558773 * kappa2 * kappa3;

    S(5, 4) = S(4, 5);
    S(5, 5) = 18478.2 - 279387 * pow(kappa2, 2);*/


    //std::cout << S.inverse() << std::endl;

    return S;
}


Eigen::MatrixXd Equivalent_StiffnessMatrix_VAM_Numerical(Variables V, Eigen::VectorXd Strain, Eigen::VectorXd inittwist, double b, std::fstream& file1)
{
    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(6, 6);

    Eigen::MatrixXd Slin = Eigen::MatrixXd::Zero(6, 6);

    //Strains
    double gamma11 = Strain(0);
    double kappa1 = Strain(3);
    double kappa2 = Strain(4);
    double kappa3 = Strain(5);

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(4, 4);
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(4, 4);
    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(4, 4);
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(4, 4);

    //linear stiffness matrix
    Slin(0, 0) = 1.633e6;
    Slin(0, 3) = -0.341e6;
    Slin(1, 1) = 0.241e6;
    Slin(1, 4) = 0.129e6;
    Slin(2, 2) = 0.754e6;
    Slin(2, 5) = -0.283e5;
    Slin(3, 3) = 0.190e6;
    Slin(3, 5) = 0.103e-4;
    Slin(4, 4) = 0.212e6;
    Slin(5, 5) = 85.205e6;

    //Values from A matrix
    A(0, 0) = 0;
    A(0, 1) = -0.002e6;
    A(1, 0) = A(0, 1);
    A(1, 1) = 0.001e6;
    A(2, 2) = 0.143e6;
    A(3, 3) = 84e6;

    //Values from B matrix
    B(0, 0) = 0.194e6;
    B(0, 1) = 21.120e6;
    B(1, 0) = B(0, 1);
    B(1, 1) = -7.539e6;
    B(2, 2) = 0.013e6;
    B(3, 3) = 9.4e6;

    //Values from C matrix
    C(0, 2) = -0.577e6;
    C(2, 0) = C(0, 2);
    C(1, 2) = 0.305e6;
    C(2, 1) = C(1, 2);

    //Values from D matrix
    D(0, 3) = 52.256e6;
    D(3, 0) = D(0, 3);
    D(1, 3) = -2.835e6;
    D(3, 1) = D(1, 3);

    //nonlinear stiffness matrix
    S(0, 0) = Slin(0, 0) + gamma11 * 6 * A(0, 0) + kappa1 * (4 * A(0, 1) + 2 * B(0, 0)) + kappa2 * (4 * A(0, 2) + 2 * C(0, 0)) + kappa3 * (4 * A(0, 3) + 2 * D(0, 0));
    S(0, 1) = Slin(0, 1);
    S(0, 2) = Slin(0, 2);
    S(0, 3) = Slin(0, 3) + gamma11 * (4 * A(0, 1) + 2 * B(0, 0)) + kappa1 * (2 * A(1, 1) + 4 * B(0, 1)) + kappa2 * (2 * A(1, 2) + 2 * B(0, 2) + 2 * C(0, 1)) + kappa3 * (2 * A(1, 3) + 2 * B(0, 3) + 2 * D(0, 1));
    S(0, 4) = Slin(0, 4) + gamma11 * (4 * A(0, 2) + 2 * C(0, 0)) + kappa1 * (2 * A(1, 2) + 2 * B(0, 2) + 2 * C(0, 1)) + kappa2 * (2 * A(2, 2) + 3 * C(0, 2)) + kappa3 * (2 * A(2, 3) + 2 * C(0, 3) + 2 * D(0, 2));
    S(0, 5) = Slin(0, 5) + gamma11 * (4 * A(0, 3) + 2 * D(0, 0)) + kappa1 * (2 * A(1, 3) + 2 * B(0, 3) + 2 * D(0, 1)) + kappa2 * (2 * A(2, 3) + 2 * C(0, 3) + 2 * D(0, 2)) + kappa3 * (2 * A(3, 3) + 4 * D(0, 3));
    S(1, 0) = Slin(0, 1);
    S(1, 1) = Slin(1, 1);
    S(1, 2) = Slin(1, 2);
    S(1, 3) = Slin(1, 3);
    S(1, 4) = Slin(1, 4);
    S(1, 5) = Slin(1, 5);

    S(2, 0) = Slin(0, 2);
    S(2, 1) = Slin(1, 2);
    S(2, 2) = Slin(2, 2);
    S(2, 3) = Slin(2, 3);
    S(2, 4) = Slin(2, 4);
    S(2, 5) = Slin(2, 5);

    S(3, 0) = Slin(0, 3);
    S(3, 1) = Slin(1, 3);
    S(3, 2) = Slin(2, 3);
    S(3, 3) = Slin(3, 3) + gamma11 * (2 * A(1, 1) + 4 * B(0, 1)) + kappa1 * 6 * B(1, 1) + kappa2 * (4 * B(1, 2) + 2 * C(1, 1)) + kappa3 * (4 * B(1, 3) + 2 * D(1, 1));
    S(3, 4) = Slin(3, 4) + gamma11 * (2 * A(1, 2) + 2 * B(0, 2) + 2 * C(0, 1)) + kappa1 * (4 * B(1, 2) + 2 * C(1, 1)) + kappa2 * (2 * B(2, 2) + 4 * C(1, 2)) + kappa3 * (2 * B(2, 3) + 2 * C(1, 3) + 2 * D(1, 2));
    S(3, 5) = Slin(3, 5) + gamma11 * (2 * A(1, 3) + 2 * B(0, 3) + 2 * D(0, 1)) + kappa1 * (4 * B(1, 3) + 2 * D(1, 1)) + kappa2 * (2 * B(2, 3) + 2 * C(1, 3) + 2 * D(1, 2)) + kappa3 * (2 * B(3, 3) + 4 * D(1, 3));

    S(4, 0) = Slin(0, 4);
    S(4, 1) = Slin(1, 4);
    S(4, 2) = Slin(2, 4);
    S(4, 3) = Slin(3, 4);
    S(4, 4) = Slin(4, 4) + gamma11 * (2 * A(2, 2) + 4 * C(0, 2)) + kappa1 * (2 * B(2, 2) + 4 * C(1, 2)) + kappa2 * 6 * C(2, 2) + kappa3 * (4 * C(2, 3) + 2 * D(2, 2));
    S(4, 5) = Slin(4, 5) + gamma11 * (2 * A(2, 3) + 2 * C(0, 3) + 2 * D(0, 2)) + kappa1 * (2 * B(2, 3) + 2 * C(1, 3) + 2 * D(1, 2)) + kappa2 * (4 * C(2, 3) + 2 * D(1, 1)) + kappa3 * (2 * C(3, 3) + 4 * D(2, 3));

    S(5, 0) = Slin(0, 5);
    S(5, 1) = Slin(1, 5);
    S(5, 2) = Slin(2, 5);
    S(5, 3) = Slin(3, 5);
    S(5, 4) = Slin(4, 5);
    S(5, 5) = Slin(5, 5) + gamma11 * (2 * A(3, 3) + 4 * D(0, 3)) + kappa1 * (2 * B(3, 3) + 4 * D(1, 3)) + kappa2 * (2 * C(3, 3) + 4 * D(2, 3)) + kappa3 * 6 * D(3, 3);

    return S;
}

//This function gives the equivalent stiffness matrix in 4*4 form after d
//differentiating twice the one dimensional strain energy with respect to strains
/*Eigen::MatrixXd Equivalent_ClassicalStiffnessModel_VAM(Stiffness S, Eigen::VectorXd Strain, std::fstream& file1)
{
    Eigen::MatrixXd Seq = Eigen::MatrixXd::Zero(4, 4);

    //Strains
    double gamma11, kappa1, kappa2, kappa3;
    gamma11 = Strain(0);
    kappa1 = Strain(3);
    kappa2 = Strain(4);
    kappa3 = Strain(5);

    /*file1 << "          Strains" << std::endl;
    file1 << "                gamma11: " << gamma11 << std::endl;
    file1 << "                kappa1: " << kappa1 << std::endl;
    file1 << "                kappa2: " << kappa2 << std::endl;
    file1 << "                kappa3: " << kappa3 << std::endl;

    //d^2U/dgamma11^2
    Seq(0, 0) = 2 * (0.5 * S.Sl(0, 0) + kappa2 * (S.Sln(0, 2) + 0.5 * S.Sn(2, 2) * kappa2));
    //d^2U/dgamma11*dkappa1
    Seq(0, 1) = 0.5 * S.Sl(0, 1) + S.Sln(0, 0) * kappa1 + kappa2 * (0.5 * S.Sln(0, 4) + 0.5 * S.Sln(1, 2) + S.Sn(0, 2) * kappa1 +
        0.5 * S.Sn(2, 4) * kappa2);
    //d^2U/dgamma11*dkappa2
    Seq(0, 2) = 0.5 * S.Sl(0, 2) + S.Sln(0, 2) * gamma11 + 0.5 * S.Sln(0, 4) * kappa1 + 0.5 * S.Sln(1, 2) * kappa1 + 0.5 * S.Sn(0, 2) * pow(kappa1, 2)
        + S.Sln(0, 1) * kappa2 + S.Sln(2, 2) * kappa2 + S.Sn(2, 2) * gamma11 * kappa2 +
        S.Sn(2, 4) * kappa1 * kappa2 + 1.5 * S.Sn(1, 2) * pow(kappa2, 2) +
        0.5 * S.Sln(0, 3) * kappa3 + 0.5 * S.Sln(2, 3) * kappa3 +
        S.Sn(2, 3) * kappa2 * kappa3;
    //d^2U/dgamma11*dkappa3
    Seq(0, 3) = 0.5 * S.Sl(0, 3) + 0.5 * kappa2 * S.Sln(0, 3) + S.Sln(2, 3) +
        S.Sn(2, 3) * kappa2;
    //d^U/dkappa1*dgamma11
    Seq(1, 0) = Seq(0, 1);
    //d^U/dkappa1^2
    Seq(1, 1) = 2 * (0.5 * S.Sl(1, 1) + S.Sln(0, 0) * gamma11 + 3 * S.Sln(0, 1) * kappa1 + 3 * S.Sn(0, 0) * pow(kappa1, 2) +
        S.Sln(0, 2) * kappa2 + S.Sln(1, 4) * kappa2 + S.Sn(0, 2) * gamma11 * kappa2 +
        3 * S.Sn(0, 4) * kappa1 * kappa2 + S.Sn(0, 1) * pow(kappa2, 2) +
        0.5 * S.Sn(4, 4) * pow(kappa2, 2) +
        S.Sln(0, 3) * kappa3 + S.Sn(0, 3) * kappa2 * kappa3);
    //d^U/dkappa1*dkappa2
    Seq(1, 2) = 0.5 * S.Sl(1, 2) + 0.5 * S.Sln(0, 4) * gamma11 + 0.5 * S.Sln(1, 2) * gamma11 +
        S.Sln(0, 2) * kappa1 + S.Sln(1, 4) * kappa1 + S.Sn(0, 2) * gamma11 * kappa1 + 1.5 * S.Sn(0, 4) * pow(kappa1, 2) +
        S.Sln(1, 1) * kappa2 + S.Sln(2, 4) * kappa2 + S.Sn(2, 4) * gamma11 * kappa2 +
        2 * S.Sn(0, 1) * kappa1 * kappa2 + S.Sn(4, 4) * kappa1 * kappa2 +
        1.5 * S.Sn(1, 4) * pow(kappa2, 2) + 0.5 * S.Sln(1, 3) * kappa3 + 0.5 * S.Sln(3, 4) * kappa3 +
        S.Sn(0, 3) * kappa1 * kappa3 + S.Sn(3, 4) * kappa2 * kappa3;
    //d^U/dkappa1*dkappa3
    Seq(1, 3) = 0.5 * S.Sl(1, 3) + S.Sln(0, 3) * kappa1 +
        kappa2 * (0.5 * S.Sln(1, 3) + 0.5 * S.Sln(3, 4) + S.Sn(0, 3) * kappa1 +
            0.5 * S.Sn(3, 4) * kappa2);
    //d^U/dkappa2*dgamma11
    Seq(2, 0) = Seq(0, 2);
    //d^U/dkappa2*dkappa1
    Seq(2, 1) = Seq(1, 2);
    //d^U/dkappa2*dkappa2
    Seq(2, 2) = 2 * (0.5 * S.Sl(2, 2) + S.Sln(0, 1) * gamma11 + S.Sln(2, 2) * gamma11 +
        0.5 * S.Sn(2, 2) * pow(gamma11, 2) + S.Sln(1, 1) * kappa1 +
        S.Sln(2, 4) * kappa1 + S.Sn(2, 4) * gamma11 * kappa1 +
        S.Sn(0, 1) * pow(kappa1, 2) + 0.5 * S.Sn(4, 4) * pow(kappa1, 2) +
        3 * S.Sln(1, 2) * kappa2 + 3 * S.Sn(1, 2) * gamma11 * kappa2 +
        3 * S.Sn(1, 4) * kappa1 * kappa2 + 3 * S.Sn(1, 1) * pow(kappa2, 2) +
        S.Sln(1, 3) * kappa3 + S.Sln(2, 3) * kappa3 + S.Sn(2, 3) * gamma11 * kappa3 +
        S.Sn(3, 4) * kappa1 * kappa3 + 3 * S.Sn(1, 3) * kappa2 * kappa3 +
        0.5 * S.Sn(3, 3) * pow(kappa3, 2));
    //d^U/dkappa2*dkappa3
    Seq(2, 3) = 0.5 * S.Sl(2, 3) + 0.5 * S.Sln(0, 3) * gamma11 + 0.5 * S.Sln(2, 3) * gamma11 +
        0.5 * S.Sln(1, 3) * kappa1 + 0.5 * S.Sln(3, 4) * kappa1 + 0.5 * S.Sn(0, 3) * pow(kappa1, 2) +
        S.Sln(1, 3) * kappa2 + S.Sln(2, 3) * kappa2 + S.Sn(2, 3) * gamma11 * kappa2 +
        S.Sn(3, 4) * kappa1 * kappa2 + 1.5 * S.Sn(1, 3) * pow(kappa2, 2) + S.Sln(3, 3) * kappa3 +
        S.Sn(3, 3) * kappa2 * kappa3;
    //d^U/dkappa3*dgamma11
    Seq(3, 0) = Seq(0, 3);
    //d^U/dkappa3*dkappa1
    Seq(3, 1) = Seq(1, 3);
    //d^U/dkappa3*dkappa2
    Seq(3, 2) = Seq(2, 3);
    //d^U/dkappa3*dkappa3
    Seq(3, 3) = 2 * (0.5 * S.Sl(3, 3) + kappa2 * S.Sln(3, 3) + 0.5 * S.Sn(3, 3) * kappa2);

    file1 << "        Printing equivalent cross-sectional stiffness matrix " << std::endl;
    for (int j = 0; j < 4; j++)
    {
        file1 << "            ";
        for (int k = 0; k < 4; k++)
        {
            file1 << Seq(j, k) << "         ";
        }
        file1 << std::endl;
    }

    /*Eigen::MatrixXd S_final = Eigen::MatrixXd::Zero(6, 6);

    for (int j = 0; j < 6; j++)
    {
        if (j != 1 && j != 2)
        {
            if (j == 0)
            {
                S_final(j, 0) = Seq(j, 0);
                S_final(j, 3) = Seq(j, 1);
                S_final(j, 4) = Seq(j, 2);
                S_final(j, 5) = Seq(j, 3);
            }
            else
            {
                S_final(j, 0) = Seq(j - 2, 0);
                S_final(j, 3) = Seq(j - 2, 1);
                S_final(j, 4) = Seq(j - 2, 2);
                S_final(j, 5) = Seq(j - 2, 3);
            }
        }
        else
            S_final(j, j) = 1e20;
    }


    return Seq;
}*/


/*Eigen::MatrixXd EquivalentStiffnessMatrix_VAM(Stiffness S, Eigen::VectorXd Strain, std::fstream& file1)
{
    Eigen::MatrixXd Seq = Eigen::MatrixXd::Zero(4, 4);

    //Strains
    double gamma11, kappa1, kappa2, kappa3;
    gamma11 = Strain(0);
    kappa1 = Strain(3);
    kappa2 = Strain(4);
    kappa3 = Strain(5);

    Seq(0, 0) = S.Sl(0, 0) + S.Sn(2, 2) * pow(kappa2, 2);
    Seq(0, 1) = S.Sl(0, 1) + kappa1 * S.Sln(0, 0) + 0.5 * pow(kappa2, 2) * S.Sn(2, 4) + kappa1 * kappa2 * S.Sn(0, 2);
    Seq(0, 2) = S.Sl(0, 2) + gamma11 * S.Sln(0, 2) + kappa2 * S.Sln(0, 1) +
                kappa2 * S.Sln(2, 2) + kappa1 * (0.5 * S.Sln(0, 4) +
                0.5 * S.Sln(1, 2) - 2 * kappa1 * S.Sn(0, 2)) +
                kappa2 * (-0.5 * gamma11 * S.Sn(2, 2) + kappa1 * S.Sn(2, 4) -
                2 * kappa2 * S.Sn(1, 2) + kappa3 * S.Sn(2, 3));
    Seq(0, 3) = S.Sl(0, 3) + 0.5 * kappa2 * (S.Sln(0, 3) + S.Sln(2, 3) +
                kappa2 * S.Sn(2, 3));

    Seq(1, 0) = Seq(0, 1);
    Seq(1, 1) = S.Sl(1, 1) + 2 * kappa1 * S.Sln(0, 1) + 2 * gamma11 * kappa2 * S.Sn(0, 2)
        + pow(kappa1, 2) * S.Sn(0, 0) + kappa2 * (6 * kappa1 * S.Sn(0, 4) +
            kappa2 * (S.Sn(4, 4) - 4 * S.Sn(0, 1)) + 2 * kappa3 * S.Sn(0, 3));
    Seq(1, 2) = S.Sl(1, 2) + 0.5 * gamma11 * S.Sln(0, 4) + 0.5 * gamma11 * S.Sln(1, 2) +
        kappa1 * (S.Sln(0, 2) + S.Sln(1, 4)) + gamma11 * kappa1 * S.Sn(0, 2) - 2 * pow(kappa1, 2) *
        S.Sn(0, 4) + kappa2 * S.Sln(1, 1) + kappa2 * S.Sln(2, 4) +
         + gamma11 * kappa2 * S.Sn(2, 4)  + 2 * kappa2 * kappa1 * S.Sn(0, 1) - 0.5 * kappa2 * kappa1 * S.Sn(4, 4)
        -2 * pow(kappa2, 2) * S.Sn(1, 4) + kappa3 * kappa1 * S.Sn(0, 3)  +
        kappa2 * kappa3 * S.Sn(3, 4);
    Seq(1, 3) = S.Sl(1, 3) + kappa1 * S.Sln(0, 3) + kappa2 * (0.5 * S.Sln(1, 3) +
        0.5 * S.Sln(3, 4) + kappa1 * S.Sn(0, 3) + 0.5 * kappa2 * S.Sn(3, 4));

    Seq(2, 0) = Seq(0, 2);
    Seq(2, 1) = Seq(1, 2);
    Seq(2, 2) = S.Sl(2, 2) + pow(gamma11, 2) * S.Sn(2, 2)
        - 3 * gamma11 * kappa1 * S.Sn(2, 4) + 2 * pow(kappa1, 2) * S.Sn(0, 1)
        + pow(kappa1, 2) * S.Sn(4, 4) + 2 * kappa2 * S.Sln(1, 2) +
        6 * kappa1 * kappa2 * S.Sn(1, 4) + pow(kappa2, 2) * S.Sn(1, 1) -
        3 * kappa1 * kappa3 * S.Sn(3, 4) + 6 * kappa2 * kappa3 * S.Sn(1, 3) +
        pow(kappa3, 2) * S.Sn(3, 3) +
        gamma11 * (6 * kappa2 * S.Sn(1, 2) - 3 * kappa3 * S.Sn(2, 3));
    Seq(2, 3) = S.Sl(2, 3) + 0.5 * gamma11 * S.Sln(0, 3) + 0.5 * kappa1 * S.Sln(3, 4)
        - 2 * pow(kappa1, 2) * S.Sn(0, 3) + S.Sln(2, 3) * (0.5 * gamma11 + kappa2) +
        (0.5 * kappa1 + kappa2) * S.Sln(1, 3) + kappa3 * S.Sln(3, 3) +
        kappa2 * (gamma11 * S.Sn(2, 3) + kappa1 * S.Sn(3, 4) - 2 * kappa2 * S.Sn(1, 3) -
            0.5 * kappa3 * S.Sn(3, 3)) ;

    Seq(3, 0) = Seq(0, 3);
    Seq(3, 1) = Seq(1, 3);
    Seq(3, 2) = Seq(2, 3);
    Seq(3, 3) = S.Sl(3, 3) + pow(kappa2, 2) * S.Sn(3, 3);

    file1 << "        Printing equivalent cross-sectional stiffness matrix for node " << std::endl;
    for (int j = 0; j < 4; j++)
    {
        file1 << "            ";
        for (int k = 0; k < 4; k++)
        {
            file1 << Seq(j, k) << "         ";
        }
        file1 << std::endl;
    }

    Eigen::MatrixXd S_final = Eigen::MatrixXd::Zero(6, 6);

    for (int j = 0; j < 6; j++)
    {
        if (j != 1 && j != 2)
        {
            if (j == 0)
            {
                S_final(j, 0) = Seq(j, 0);
                S_final(j, 3) = Seq(j, 1);
                S_final(j, 4) = Seq(j, 2);
                S_final(j, 5) = Seq(j, 3);
            }
            else
            {
                S_final(j, 0) = Seq(j - 2, 0);
                S_final(j, 3) = Seq(j - 2, 1);
                S_final(j, 4) = Seq(j - 2, 2);
                S_final(j, 5) = Seq(j - 2, 3);
            }
        }
        else
            S_final(j, j) = 1e20;
    }


    return Seq;
}*/
#endif
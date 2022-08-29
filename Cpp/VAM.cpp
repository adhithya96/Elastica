
//Couple of things to figure out

//What units to use for k1
//FDM code for first and last node
//How to validate? Experimental data is not given. Just the graphs.
#ifndef VARIABLES_H
#define VARIABLES_H


#include "Variables.h"
#include<iostream>
#include <corecrt_math_defines.h>
#define M_PI 3.14159265358979323846

struct PostProcessing
{
    Eigen::MatrixXd Jacobian;
    Eigen::VectorXd R;
};

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
void AMatrix(Eigen::MatrixXd Q, double z0, double z1, Eigen::MatrixXd &A)
{
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            A(i, j) += Q(i, j) * (z1 - z0);
}

//B matrix
void BMatrix(Eigen::MatrixXd Q, double z0, double z1, Eigen::MatrixXd &B)
{
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            B(i, j) += (1.0 / 2.0) * Q(i, j) * (pow(z1, 2) - pow(z0, 2));
}

//D matrix
void DMatrix(Eigen::MatrixXd Q, double z0, double z1, Eigen::MatrixXd &D)
{
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            D(i, j) += (1.0 / 3.0) * Q(i, j) * (pow(z1, 3) - pow(z0, 3));
}

//Variables needed for Calculating VAM simplified Stiffness matrix
Variables StiffnessVariableCalc(Eigen::MatrixXd A, Eigen::MatrixXd B, Eigen::MatrixXd D)
{
    Variables V;

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

    return V;
}

//Stiffness Matrices after simplification using VAM
Stiffness StiffnessMatrix_VAM(double E1, double E2, double nu12, double  G12, double nu21, double width, double height,
    double k1, Eigen::VectorXd Orient, int np, std::fstream &file1)
{
    double b = width;

    //Midplane thickness and calculate hi
    double mt = np * height / 2;

    file1 << "            Midplane thickness  " << mt << std::endl;
    Eigen::VectorXd h(np + 1);
    h(0) = -mt;
    for (int i = 1; i < np + 1; i++)
        h(i) = h(i - 1) + height;
    //assert(h(np)==mt);
//    std::cout<<h(np)<<std::endl;
//    std::cout<<mt<<std::endl;

    //Reduced stiffness matrix at different angles
    Eigen::MatrixXd Q = ReducedStiffnessMatrix(E1, E2, nu12, nu21, G12);
    file1 << "            Reduced Stiffness Matrix" << std::endl;
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
        }*/
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
        }*/
        double z0 = h(i);
        double z1 = h(i + 1);
        AMatrix(Qdash, z0, z1, A);
        BMatrix(Qdash, z0, z1, B);
        DMatrix(Qdash, z0, z1, D);
    }
    /*file1 << "            A matrix" << std::endl;
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
    }*/

    Variables V = StiffnessVariableCalc(A, B, D);

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
    file1 << "            Valeus of variables with double bar from appendix" << std::endl;
    file1 << "                A11:  " << V.BB.A11 << std::endl;
    file1 << "                B11:  " << V.BB.B11 << std::endl;
    file1 << "                B16:  " << V.BB.B16 << std::endl;
    file1 << "                D11:  " << V.BB.D11 << std::endl;
    file1 << "                D16:  " << V.BB.D16 << std::endl;
    file1 << "                D66:  " << V.BB.D66 << std::endl;*/

    Stiffness S;
    S.Sl = Eigen::MatrixXd::Zero(4, 4);
    S.Sln = Eigen::MatrixXd::Zero(4, 5);
    S.Sn = Eigen::MatrixXd::Zero(5, 5);

    S.Sl(0, 0) = b * V.BB.A11;
    S.Sl(0, 1) = (-2) * b * V.BB.B16 + pow(b, 3) * V.BB.A11 * k1 / 12;
    S.Sl(0, 2) = b * V.BB.B11;
    S.Sl(1, 0) = S.Sl(0, 1);
    S.Sl(1, 1) = (4 * b * V.BB.D66 - pow(b, 3) * V.BB.B16 * k1 / 3 + pow(b, 5) * V.BB.A11 * pow(k1, 2) / 80);
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
    S.Sln(3, 4) = -pow(b, 5) * V.BB.A11 * V.B.B12 / (720 * V.B.D22);

    S.Sn(0, 0) = pow(b, 5) * V.BB.A11 / 320;
    S.Sn(0, 2) = pow(b, 5) * V.BB.A11 * V.B.B12 / (720 * V.B.D22);
    S.Sn(0, 4) = -pow(b, 5) * V.BB.A11 * V.B.D26 / (360 * V.B.D22) + pow(b, 7) * V.BB.A11 * V.B.B12 * k1 / (10080 * V.B.D22);
    S.Sn(1, 1) = pow(b, 5) * V.BB.A11 * pow(V.B.D12, 2) / (720 * pow(V.B.D22, 2));
    S.Sn(1, 2) = pow(b, 5) * V.BB.A11 * V.B.B12 * V.B.D12 / (720 * pow(V.B.D22, 2));
    S.Sn(1, 4) = -pow(b, 5) * V.BB.A11 * V.B.D12 * V.B.D26 / (360 * pow(V.B.D22, 2))
        - pow(b, 7) * V.BB.A11 * V.B.B12 * V.B.D12 * k1 / (60480 * pow(V.B.D22, 2));
    S.Sn(2, 2) = pow(b, 5) * V.BB.A11 * pow(V.B.B12, 2) / (720 * pow(V.B.D22, 2));
    S.Sn(2, 3) = -pow(b, 5) * V.BB.A11 * V.B.B12 * V.B.D26 / (360 * pow(V.B.D22, 2))
        - pow(b, 7) * V.BB.A11 * pow(V.B.B12, 2) * k1 / (60480 * pow(V.B.D22, 2));
    S.Sn(3, 3) = (pow(b, 7) * V.BB.A11 * pow(V.B.B12, 2) / (10080 * pow(V.B.D22, 2)) - pow(b, 7) * pow(V.BB.A11, 2) / (30240 * V.B.D22));
    S.Sn(4, 4) = (pow(b, 5) * V.BB.A11 * pow(V.B.D26, 2) / (180 * pow(V.B.D22, 2)) + pow(b, 5) * V.BB.A11 * V.B.D12 / (360 * V.B.D22)
        + pow(b, 7) * V.BB.A11 * V.B.B12 * V.B.D26 * k1 / (15120 * pow(V.B.D22, 2))
        - (pow(b, 9) * pow(V.BB.A11, 2) / (90720 * V.B.D22) + pow(b, 9) * V.BB.A11 * pow(V.B.B12, 2) / (403200 * pow(V.B.D22, 2))) * pow(k1, 2));

    //    std::cout<<S<<std::endl;

    return S;
}

//This function gives the equivalent stiffness matrix in 4*4 form after d
//differentiating twice the one dimensional strain energy with respect to strains
Eigen::MatrixXd Equivalent_ClassicalStiffnessModel_VAM(Stiffness S, Eigen::VectorXd Strain, std::fstream &file1)
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
    file1 << "                kappa3: " << kappa3 << std::endl;*/

    //d^2U/dgamma11^2
    Seq(0, 0) = 0.5 * S.Sl(0, 0) + kappa2 * (S.Sln(0, 2) + 0.5 * S.Sn(2, 2) * kappa2);
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
    Seq(1, 1) = 0.5 * S.Sl(1, 1) + S.Sln(0, 0) * gamma11 + 3 * S.Sln(0, 1) * kappa1 + 3 * S.Sn(0, 0) * pow(kappa1, 2) +
        S.Sln(0, 2) * kappa2 + S.Sln(1, 4) * kappa2 + S.Sn(0, 2) * gamma11 * kappa2 +
        3 * S.Sn(0, 4) * kappa1 * kappa2 + S.Sn(0, 1) * pow(kappa2, 2) +
        0.5 * S.Sn(4, 4) * pow(kappa2, 2) +
        S.Sln(0, 3) * kappa3 + S.Sn(0, 3) * kappa2 * kappa3;
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
    Seq(2, 2) = 0.5 * S.Sl(2, 2) + S.Sln(0, 1) * gamma11 + S.Sln(2, 2) * gamma11 +
        0.5 * S.Sn(2, 2) * pow(gamma11, 2) + S.Sln(1, 1) * kappa1 +
        S.Sln(2, 4) * kappa1 + S.Sn(2, 4) * gamma11 * kappa1 +
        S.Sn(0, 1) * pow(kappa1, 2) + 0.5 * S.Sn(4, 4) * pow(kappa1, 2) +
        3 * S.Sln(1, 2) * kappa2 + 3 * S.Sn(1, 2) * gamma11 * kappa2 +
        3 * S.Sn(1, 4) * kappa1 * kappa2 + 3 * S.Sn(1, 1) * pow(kappa2, 2) +
        S.Sln(1, 3) * kappa3 + S.Sln(2, 3) * kappa3 + S.Sn(2, 3) * gamma11 * kappa3 +
        S.Sn(3, 4) * kappa1 * kappa3 + 3 * S.Sn(1, 3) * kappa2 * kappa3 +
        0.5 * S.Sn(3, 3) * pow(kappa3, 2);
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
    Seq(3, 3) = 0.5 * S.Sl(3, 3) + kappa2 * S.Sln(3, 3) + 0.5 * S.Sn(3, 3) * kappa2;

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


    return S_final;
}

Eigen::MatrixXd Equivalent_TangentStiffnessMatrix_VAM(Stiffness S, Eigen::VectorXd Strain, std::fstream &file1)
{
    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(6, 6);

    double gamma11, kappa1, kappa2, kappa3;
    gamma11 = Strain(0);
    kappa1 = Strain(3);
    kappa2 = Strain(4);
    kappa3 = Strain(5);

    /*T(0, 0) = 0.5 * S.Sl(0, 0) + kappa2.Subscript[Sln, 13] +
            1.5 Subscript[Sn, 33] Subscript[\[Kappa], 2]),
        0.5 Subscript[Sl, 12] + 2. Subscript[Sln, 11] Subscript[\[Kappa], 1] +
        Subscript[\[Kappa],
        2](1. Subscript[Sln, 15] + 1. Subscript[Sln, 23] +
            3. Subscript[Sn, 13] Subscript[\[Kappa], 1] +
            1.5 Subscript[Sn, 35] Subscript[\[Kappa], 2]),
        0.5 Subscript[Sl, 13] +
        2. Subscript[Sln, 13] Subscript[\[Gamma], 11] +
        1. Subscript[Sln, 15] Subscript[\[Kappa], 1] +
        1. Subscript[Sln, 23] Subscript[\[Kappa], 1] + 1.5 Subscript[Sn, 13]
        \!\(\ * SubsuperscriptBox[\(\[Kappa]\), \(1\), \(2\)]\) +
        2. Subscript[Sln, 12] Subscript[\[Kappa], 2] +
        2. Subscript[Sln, 33] Subscript[\[Kappa], 2] +
        3. Subscript[Sn, 33] Subscript[\[Gamma], 11] Subscript[\[Kappa],
        2] + 3. Subscript[Sn, 35] Subscript[\[Kappa], 1] Subscript[\[Kappa],
        2] + 4.5 Subscript[Sn, 23]
        \!\(\ * SubsuperscriptBox[\(\[Kappa]\), \(2\), \(2\)]\) +
        1. Subscript[Sln, 14] Subscript[\[Kappa], 3] +
        1. Subscript[Sln, 34] Subscript[\[Kappa], 3] +
        3. Subscript[Sn, 34] Subscript[\[Kappa], 2] Subscript[\[Kappa], 3],
        0.5 Subscript[Sl, 14] +
        Subscript[\[Kappa],
        2](1. Subscript[Sln, 14] + 1. Subscript[Sln, 34] +
            1.5 Subscript[Sn, 34] Subscript[\[Kappa], 2])*/

    return T;
}

#endif
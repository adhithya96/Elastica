
//Couple of things to figure out

//What units to use for k1
//FDM code for first and last node
//How to validate? Experimental data is not given. Just the graphs.
#ifndef VARIABLES_H
#define VARIABLES_H

#include "Variables.h"
#include<iostream>


struct PostProcessing
{
    Eigen::MatrixXd Jacobian;
    Eigen::VectorXd R;
};

//Q matrix
Eigen::MatrixXd ReducedStiffnessMatrix(double E1, double E2, double nu12, double nu21, double G12)
{
    Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(3,3);
    Q(0,0) = E1/(1-nu12*nu21);
    Q(0,1) = nu12*E2/(1-nu12*nu21);
    Q(1,1) = E2/(1-nu12*nu21);
    Q(1,0) = Q(0,1);
    Q(2,2) = G12;

//    std::cout<<Q<<std::endl;

    return Q;
}

Eigen::MatrixXd TransformationMatrix(double theta)
{
    double m = cos(M_PI*theta/180);
    double n = sin(M_PI*theta/180);

    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(3,3);
    T(0,0) = pow(m,2);
    T(0,1) = pow(n,2);
    T(0,2) = -2*m*n;
    T(1,0) = pow(n,2);
    T(1,1) = pow(m,2);
    T(1,2) = 2*m*n;
    T(2,0) = m*n;
    T(2,1) = -m*n;
    T(2,2) = pow(m,2) - pow(n,2);

//    std::cout<<T<<std::endl;

    return T;
}

//A matrix
Eigen::MatrixXd AMatrix(Eigen::MatrixXd Q, double z0, double z1)
{
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(3,3);
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            A(i,j) = Q(i,j)*(z1-z0);

//    std::cout<<A<<std::endl;

    return A;
}

//B matrix
Eigen::MatrixXd BMatrix(Eigen::MatrixXd Q, double z0, double z1)
{
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3,3);
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            B(i,j) = (1.0/2.0)*Q(i,j)*(pow(z1,2)-pow(z0,2));

//    std::cout<<B<<std::endl;

    return B;
}

//D matrix
Eigen::MatrixXd DMatrix(Eigen::MatrixXd Q, double z0, double z1)
{
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(3,3);
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            D(i,j) = (1.0/3.0)*Q(i,j)*(pow(z1,3)-pow(z0,3));

//    std::cout<<D<<std::endl;

    return D;
}

//Variables needed for Calculating VAM simplified Stiffness matrix
Variables StiffnessVariableCalc(Eigen::MatrixXd A, Eigen::MatrixXd B, Eigen::MatrixXd D)
{
    Variables V;

    double A11 = A(0,0);
    double A12 = A(0,1);
    double A16 = A(0,2);
    double A22 = A(1,1);
    double A26 = A(1,2);
    double A66 = A(2,2);

    double B11 = B(0,0);
    double B12 = B(0,1);
    double B16 = B(0,2);
    double B22 = B(1,1);
    double B26 = B(1,2);
    double B66 = B(2,2);

    double D11 = D(0,0);
    double D12 = D(0,1);
    double D16 = D(0,2);
    double D22 = D(1,1);
    double D26 = D(1,2);
    double D66 = D(2,2);

    V.B.A11 = A11 + (pow(A16,2)*A22 - 2*A12*A16*A26 + pow(A12,2)*A66)/(pow(A26,2)-A22*A66);
    V.B.B11 = B11 + (A12*A66*B12 + A16*A22*B16 - A26*(A16*B12 + A12*B16))/(pow(A26,2) - A22*A66);
    V.B.B12 = B12 + (A12*A66*B22 + A16*A22*B26 - A26*(A16*B22 + A12*B26))/(pow(A26,2) - A22*A66);
    V.B.B16 = B16 + (A12*A66*B26 + A16*A22*B66 - A26*(A16*B26 + A12*B66))/(pow(A26,2)-A22*A66);

    V.B.D11 = D11 + (A66*pow(B12,2) - 2*A26*B12*B16 + A22*pow(B16,2))/(pow(A26,2) - A22*A66);
    V.B.D12 = D12 + (A66*B12*B22 + A22*B16*B26 - A26*(B16*B22 + B12*B26))/(pow(A26,2) - A22*A66);
    V.B.D22 = D22 + (A66*pow(B22,2) - 2*A26*B22*B26 + A22*pow(B26,2))/(pow(A26,2) - A22*A66);
    V.B.D16 = D16 + (A66*B12*B26 + A22*B16*B66 - A26*(B16*B26 + B12*B66))/(pow(A26,2) - A22*A66);
    V.B.D26 = D26 + (A66*B22*B26 + A22*B26*B66 - A26*(pow(B26,2) + B22*B66))/(pow(A26,2) - A22*A66);
    V.B.D66 = D66 + (A66*pow(B26,2) - 2*A26*B26*B66 + A22*pow(B66,2))/(pow(A26,2) - A22*A66);

    V.BB.A11 = V.B.A11 - pow(V.B.B12,2)/V.B.D22;
    V.BB.B11 = V.B.B11 - V.B.B12*V.B.D12/V.B.D22;
    V.BB.B16 = V.B.B16 - V.B.B12*V.B.D26/V.B.D22;
    V.BB.D11 = V.B.D11 - pow(V.B.D12,2)/V.B.D22;
    V.BB.D16 = V.B.D16 - V.B.D12*V.B.D26/V.B.D22;
    V.BB.D66 = V.B.D66 - pow(V.B.D26,2)/V.B.D22;

    return V;
}

//Stiffness Matrices after simplification using VAM
Eigen::MatrixXd StiffnessMatrix_VAM(double E1,double E2,double nu12,double  G12, double nu21, double width, double height,
                                    double k1, Eigen::VectorXd Orient, int np)
{
    double b = width;
    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(9,9);

    //Midplane thickness and calculate hi
    double mt = np*height/2;
    Eigen::VectorXd h(np+1);
    h(0) = -mt;
    for(int i=1;i<np+1;i++)
        h(i) = h(i-1) + height;
    //assert(h(np)==mt);
//    std::cout<<h(np)<<std::endl;
//    std::cout<<mt<<std::endl;

    //Reduced stiffness matrix at different angles
    Eigen::MatrixXd Q = ReducedStiffnessMatrix(E1,E2,nu12,nu21,G12);
    Eigen::MatrixXd A(3,3), B(3,3), D(3,3);
    for(int i=0;i<Orient.size();i++)
    {
        Eigen::MatrixXd T = TransformationMatrix(Orient(i));
        Eigen::MatrixXd Qdash = T*Q*(T.transpose());
//        std::cout<<Qdash<<std::endl;
        double z0 = h(i);
        double z1 = h(i+1);
        A += AMatrix(Qdash,z0,z1);
//        std::cout<<A<<std::endl;
        B += BMatrix(Qdash,z0,z1);
        D += DMatrix(Qdash,z0,z1);
    }
//    std::cout<<A<<std::endl;
//    std::cout<<B<<std::endl;
//    std::cout<<D<<std::endl;
    //assert(A.transpose()==A&&B.transpose()==B&&D.transpose()==D);


    Variables V = StiffnessVariableCalc(A,B,D);

    S(0,0) = b*V.BB.A11;
    S(0,1) = (-2)*b*V.BB.B16 + pow(b,3)*V.BB.A11*k1/12;
    S(0,2) = b*V.BB.B11;
    S(1,0) = S(0,1);
    S(1,1) = (4*b*V.BB.D66 - pow(b,3)*V.BB.B16*k1/3 + pow(b,5)*V.BB.A11*pow(k1,2)/80);
    S(1,2) = -2*b*V.BB.D16 + pow(b,3)*V.BB.B11*k1/12;
    S(2,0) = S(0,2);
    S(2,1) = S(1,2);
    S(2,2) = b*V.BB.D11;
    S(3,3) = pow(b,3)*V.BB.A11/12;

    S(0,4) = pow(b,3)*V.BB.A11/24;
    S(1,4) = -pow(b,3)*V.BB.B16/12 + pow(b,5)*V.BB.A11*k1/160;
    S(1,5) = pow(b,5)*V.BB.A11*V.B.D12*k1/(360*V.B.D22);
    S(1,6) = pow(b,5)*V.BB.A11*V.B.B12*k1/(360*V.B.D22);
    S(2,4) = (pow(b,3)*V.BB.B11/24 + pow(b,5)*V.BB.A11*V.B.D26*k1/(180*V.B.D22)
            + pow(b,7)*V.BB.A11*V.B.B12*pow(k1,2)/(10080*V.B.D22));
    S(3,7) = -pow(b,5)*V.BB.A11*V.B.B12/(720*V.B.D22);

    S(4,0) = S(0,4);
    S(4,1) = S(1,4);
    S(5,1) = S(1,5);
    S(6,1) = S(1,6);
    S(4,2) = S(2,4);
    S(7,3) = S(3,7);

    S(4,4) = pow(b,5)*V.BB.A11/320;
    S(4,6) = pow(b,5)*V.BB.A11*V.B.B12/(720*V.B.D22);
    S(4,8) = -pow(b,5)*V.BB.A11*V.B.D26/(360*V.B.D22) + pow(b,7)*V.BB.A11*V.B.B12*k1/(10080*V.B.D22);
    S(5,5) = pow(b,5)*V.BB.A11*pow(V.B.D12,2)/(720*pow(V.B.D22,2));
    S(5,6) = pow(b,5)*V.BB.A11*V.B.B12*V.B.D12/(720*pow(V.B.D22,2));
    S(5,8) = -pow(b,5)*V.BB.A11*V.B.D12*V.B.D26/(360*pow(V.B.D22,2))
            - pow(b,7)*V.BB.A11*V.B.B12*V.B.D12*k1/(60480*pow(V.B.D22,2));
    S(6,6) = pow(b,5)*V.BB.A11*pow(V.B.B12,2)/(720*pow(V.B.D22,2));
    S(6,8) = -pow(b,5)*V.BB.A11*V.B.B12*V.B.D26/(360*pow(V.B.D22,2))
            - pow(b,7)*V.BB.A11*pow(V.B.B12,2)*k1/(60480*pow(V.B.D22,2));
    S(7,7) = (pow(b,7)*V.BB.A11*pow(V.B.B12,2)/(10080*pow(V.B.D22,2)) - pow(b,7)*pow(V.BB.A11,2)/(30240*V.B.D22));
    S(8,8) = (pow(b,5)*V.BB.A11*pow(V.B.D26,2)/(180*pow(V.B.D22,2)) + pow(b,5)*V.BB.A11*V.B.D12/(360*V.B.D22)
            + pow(b,7)*V.BB.A11*V.B.B12*V.B.D26*k1/(15120*pow(V.B.D22,2))
            - (pow(b,9)*pow(V.BB.A11,2)/(90720*V.B.D22) + pow(b,9)*V.BB.A11*pow(V.B.B12,2)/(403200*pow(V.B.D22,2)))*pow(k1,2));

    S(6,4) = S(4,6);
    S(8,4) = S(4,8);
    S(6,5) = S(5,6);
    S(8,5) = S(5,8);
    S(8,6) = S(6,8);

//    std::cout<<S<<std::endl;

    return S;

}

/*
Eigen::MatrixXd RearrangeStiffnessMatrix_VAM(Eigen::MatrixXd S, double gamma11, double kappa1, double kappa2, double kappa3)
{
    Eigen::MatrixXd Sf = Eigen::MatrixXd::Zero(4,4);

 /*   Sf(0,0) = S(0,0) + S(1,0) + S(2,0) + (S(0,6) + S(1,6) + S(2,6))*kappa2;
    Sf(0,1) = S(0,1) + S(1,1) + S(2,1) + (S(0,4) + S(1,4) + S(2,4))*kappa1 + (S(0,8) + S(1,8) + S(2,8))*kappa2;
    Sf(0,2) = S(0,2) + S(1,2) + S(2,2) + (S(0,5) + S(1,5) + S(2,5))*kappa2 + (S(0,7) + S(1,7) + S(2,7))*kappa3;
    Sf(0,3) = S(0,3) + S(1,3) + S(2,3);

    Sf(1,0) = S(3,0) + S(4,0) + (S(3,6) + S(4,6))*kappa2;
    Sf(1,1) = S(3,1) + S(4,1) +(S(3,4) + S(4,4))*kappa1 + (S(3,8) + S(4,8))*kappa2;
    Sf(1,2) = S(3,2) + S(4,2) + (S(3,5) + S(4,5))*kappa2 + (S(3,7) + S(4,7))*kappa3;
    Sf(1,3) = S(3,3) + S(4,3);

    Sf(2,0) = S(5,0) + S(6,0) + (S(5,6) + S(6,6))*kappa2;
    Sf(2,1) = S(5,1) + S(6,1) +(S(5,4) + S(6,4))*kappa1 + (S(5,8) + S(6,8))*kappa2;
    Sf(2,2) = S(5,2) + S(6,2) + (S(5,5) + S(6,5))*kappa2 + (S(5,7) + S(6,7))*kappa3;
    Sf(2,3) = S(5,3) + S(6,3);

    Sf(3,0) = S(7,0) + S(8,0) + (S(7,6) + S(8,6))*kappa2;
    Sf(3,1) = S(7,1) + S(8,1) +(S(7,4) + S(8,4))*kappa1 + (S(7,8) + S(8,8))*kappa2;
    Sf(3,2) = S(7,2) + S(8,2) + (S(7,5) + S(8,5))*kappa2 + (S(7,7) + S(8,7))*kappa3;
    Sf(3,3) = S(7,3) + S(8,3);*/
/*
    for(int i=0;i<4;i++)
    {
        Sf(i,0) = S(i,0) + S(i,6)*kappa2;
        Sf(i,1) = S(i,1) + S(i,4)*kappa1 + S(i,8)*kappa2;
        Sf(i,2) = S(i,2) + S(i,5)*kappa2 + S(i,7)*kappa3;
        Sf(i,3) = S(i,3);
    }*/

    /*
    Sf(0,0) = S(0,0);
    Sf(0,1) = S(0,1) + S(0,4)*kappa1 + S(0,8)*kappa2;
    Sf(0,2) = S(0,2) + S(0,5)*kappa2 + S(0,6)*gamma11;
    Sf(0,3) = S(0,3) + S(0,7)*kappa3;

    Sf(1,0) = ;*/
    /*
    //S_linear
    for(int i=0;i<4;i++)
        for(int j=0;j<4;j++)
            Sf(i,j) = S(i,j);

    //S_ln
    for(int i=0;i<4;i++)
    {
        Sf(i,0) += 2*S(i,2)*kappa2;
        Sf(i,1) += 2*S(i,4)*kappa2 + S(i,0)*kappa1 ;
        Sf(i,2) += 2*S(i,1)*kappa2;
        Sf(i,3) += 2*S(i,3)*kappa2;
    }

    //S_n
    //Rearranging Elements of Sn matrix
    //First by rows
    Eigen::MatrixXd Temp = Eigen::MatrixXd::Zero(5,4);
    for(int i=0;i<5;i++)
    {
        Temp(i,0) += S(i+4,2+4)*kappa2;
        Temp(i,1) += S(i+4,4+4)*kappa2 + S(i+4,0+4)*kappa1 ;
        Temp(i,2) += S(i+4,1+4)*kappa2;
        Temp(i,3) += S(i+4,3+4)*kappa2;
    }
    //Rearranging columns and adding it to the final 4*4 stiffness matrix
    for(int j=0;j<4;j++)
    {
        Sf(0,j) += Temp(2,j)*kappa2;
        Sf(1,j) += Temp(4,j)*kappa2 + Temp(0,j)*kappa1;
        Sf(2,j) += Temp(1,j)*kappa2;
        Sf(3,j) += Temp(3,j)*kappa2;
    }
    Eigen::MatrixXd Temp(9,4);
    for(int i=0;i<9;i++)
    {
        Temp(i,0) = S(i,0) + S(i,6)*kappa2;
        Temp(i,1) = S(i,1) + S(i,4)*kappa1 + S(i,8)*kappa2;
        Temp(i,2) = S(i,2) + S(i,5)*kappa2 + S(i,7)*kappa3;
        Temp(i,3) = S(i,3);
    }
    for(int i=0;i<4;i++)
    {
        Sf(0,i) = Temp(0,i) + Temp(6,i)*kappa2;
        Sf(1,i) = Temp(1,i) + Temp(4,i)*kappa1 + Temp(8,i)*kappa2;
        Sf(2,i) = Temp(2,i) + Temp(5,i)*kappa2 + Temp(7,i)*kappa3;
        Sf(3,i) = Temp(3,i);
    }

    std::cout<<Sf<<std::endl;

    return Sf;

}

Eigen::MatrixXd TangentStiffnessMatrix_VAM(Eigen::MatrixXd S, double gamma11, double kappa1, double kappa2, double kappa3)
{
    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(4,4);
/*
    T(0,0) = S(0,0) + S(1,0) + S(2,0) + (S(0,6) + S(1,6) + S(2,6))*kappa2;
    T(0,1) = S(0,1) + S(1,1) + S(2,1) + 2*(S(0,4) + S(1,4) + S(2,4))*kappa1 + (S(0,8) + S(1,8) + S(2,8))*kappa2;
    T(0,2) = S(0,2) + S(1,2) + S(2,2) + 2*(S(0,5) + S(1,5) + S(2,5))*kappa2 + (S(0,7) + S(1,7) + S(2,7))*kappa3
            + (S(0,6) + S(1,6) + S(2,6))*gamma11;
    T(0,3) = (S(0,7) + S(1,7) + S(2,7))*kappa2 + S(0,3) + S(1,3) + S(2,3);

    T(1,0) = S(3,0) + S(4,0) + (S(3,6) + S(4,6))*kappa2;
    T(1,1) = S(3,1) + S(4,1) + 2*(S(3,4) + S(4,4))*kappa1 + (S(3,8) + S(4,8))*kappa2;
    T(1,2) = S(3,2) + S(4,2) + 2*(S(3,5) + S(4,5))*kappa2 + (S(3,7) + S(4,7))*kappa3
            + (S(3,6) + S(4,6))*gamma11;
    T(1,3) = (S(3,7) + S(4,7))*kappa2 + S(3,3) + S(4,3);

    T(2,0) = S(5,0) + S(6,0) + (S(5,6) + S(6,6))*kappa2;
    T(2,1) = S(5,1) + S(6,1) + 2*(S(5,4) + S(6,4))*kappa1 + (S(5,8) + S(6,8))*kappa2;
    T(2,2) = S(5,2) + S(6,2) + 2*(S(5,5) + S(6,5))*kappa2 + (S(5,7) + S(6,7))*kappa3
            + (S(5,6) + S(6,6))*gamma11;
    T(2,3) = (S(5,7) + S(6,7))*kappa2 + S(5,3) + S(6,3);

    T(3,0) = S(7,0) + S(8,0) + (S(7,6) + S(8,6))*kappa2;
    T(3,1) = S(7,1) + S(8,1) + 2*(S(7,4) + S(8,4))*kappa1 + (S(7,8) + S(8,8))*kappa2;
    T(3,2) = S(7,2) + S(8,2) + 2*(S(7,5) + S(8,5))*kappa2 + (S(7,7) + S(8,7))*kappa3
            + (S(7,6) + S(8,6))*gamma11;
    T(3,3) = (S(7,7) + S(8,7))*kappa2 + S(7,3) + S(8,3);
*/
/*
    for(int i=0;i<4;i++)
    {
        T(i,0) = S(i,0) + S(i,6)*kappa2;
        T(i,1) = S(i,1) + 2*S(i,4)*kappa1 + S(i,8)*kappa2;
        T(i,2) = S(i,2) + 2*S(i,5)*kappa2 + S(i,6)*gamma11;
        T(i,3) = S(i,7)*kappa2 + S(i,3);
    }*/
/*
    for(int i=0;i<4;i++)
        for(int j=0;j<4;j++)
            T(i,j) = S(i,j);

    //1st row
    T(0,0) = S(0,0) + S(0,6)*kappa2 + (S(6,0) + S(6,6)*kappa2)*kappa2;
    T(0,1) = S(0,1) + 2*S(0,4)*kappa1 + S(0,8)*kappa2 + S(6,6)*kappa2 + 2*S(6,4)*kappa1*kappa2 + S(6,8)*pow(kappa2,2);
    T(0,2) = (S(0,6) + S(6,0) + 2*S(6,6))*gamma11 + (S(0,8) + S(6,1) + S(6,4)*kappa1 + 2*S(6,8)*kappa2)*kappa1 +
            S(0,2) + S(0,5)*kappa2 + S(0,7)*kappa3 + (S(6,2) + S(6,5)*kappa2 + S(6,7)*kappa3)*kappa2
            + (S(0,5) + 2*S(6,5)*kappa2 + S(6,7)*kappa3)*kappa2 + S(6,3)*kappa3;
    T(0,3) = (S(0,7) + S(6,7)*kappa2)*kappa2 + S(0,3) + S(6,3)*kappa2;

    //2nd row
    T(1,0) = S(1,0) + S(1,6)*kappa2 + (S(4,0) + S(4,6)*kappa2)*kappa1 + (S(8,0) + S(8,6)*kappa2)*kappa2;
    T(1,1) = (S(4,0) + S(4,6)*kappa2)*gamma11 + S(1,1) + 2*S(1,4)*kappa1 + S(1,8)*kappa2 + S(4,1) + 2*S(4,4)*kappa1
            + S(4,8)*kappa2 + S(8,1)*kappa2 + 2*S(8,4)*kappa2*kappa1 + S(8,8)*kappa2*kappa2
            + (S(4,2) + S(4,5)*kappa2 + S(4,7)*kappa3)*kappa2 + S(4,3)*kappa3;
    T(1,2) = (S(1,6) + S(4,6)*kappa1 + S(8,0) + 2*S(8,6)*kappa2)*gamma11
            + (S(1,8) + S(4,8)*kappa1 + 2*S(8,8)*kappa2)*kappa1
            + (S(1,2) + S(1,5)*kappa2 + S(1,7)*kappa3) + (S(4,2) + S(4,5)*kappa2 + S(4,7)*kappa3)*kappa1 + (S(8,2) + S(8,4)*kappa2 + S(8,7)*kappa3)*kappa2
            + (S(1,5) + S(4,5) + S(8,2) + 2*S(8,4)*kappa2 + S(8,7)*kappa3)*kappa2
            + S(8,3)*kappa3;
    T(1,3) = (S(1,7) + S(4,7)*kappa1 + S(8,7)*kappa2)*kappa2
            + S(1,3) + S(4,3)*kappa1 + S(8,3)*kappa2;

    //3rd row
    T(2,0) = S(2,0) + S(2,6)*kappa2 + (S(5,0) + S(5,6)*kappa2)*kappa2 + (S(7,0) + S(7,6)*kappa2)*kappa3;
    T(2,1) = S(2,1) + S(2,4)*kappa1 + S(2,8)*kappa2 + (S(5,1) + S(5,4)*kappa1 + S(5,8)*kappa2)*kappa1
            + (S(7,1) + S(7,4)*kappa1 + S(7,8)*kappa2)*kappa2 +
            (S(2,4) + (S(5,1) + S(5,4)*kappa1 + S(5,8)*kappa2) + S(5,4)*kappa1 + S(7,4)*kappa2)*kappa1;
    T(2,2) = (S(2,6) + 2*S(5,6)*kappa2 + S(7,6)*kappa3)*gamma11
            + (S(2,8) + S(5,8)*kappa1 + 2*S(7,8)*kappa2)*kappa1
            + S(2,2) + S(2,5)*kappa2 + S(2,7)*kappa3 + (S(5,2) + S(5,5)*kappa2 + S(5,7)*kappa3)*kappa2 + (S(7,2) + S(7,5)*kappa2 + S(7,7)*kappa3)*kappa3
            + (S(2,4) + S(5,2) + S(5,5)*kappa2 + S(5,7)*kappa3 + S(5,5)*kappa2 + S(7,5)*kappa3)*kappa2
            + S(2,5)*kappa3;
    T(2,3) = S(7,6)*kappa2*gamma11
            + S(2,7) + S(5,7)*kappa2 + S(7,2) + S(7,5)*kappa2 + 2*S(7,7)*kappa3
            + S(2,3) + S(2,5)*kappa2 + 2*S(2,7)*kappa3;

    //4th row
    T(3,0) = S(3,0) + S(3,6)*kappa2;
    T(3,1) = S(3,1) + 2*S(3,4)*kappa1 + S(3,8)*kappa2;
    T(3,2) = S(3,6)*gamma11 + S(3,8)*kappa1 + S(3,2) + 2*S(3,5)*kappa2 + S(3,7)*kappa3;
    T(3,3) = S(3,7)*kappa2 + S(3,3);

    std::cout<<T<<std::endl;

    return T;
}

Eigen::VectorXd LocalForceVec_VAM(Eigen::MatrixXd LOAD, int a, int b)
{
    Eigen::VectorXd f = Eigen::VectorXd::Zero(8);

    if(a==LOAD(0,0)-1)
    {
        int dir = LOAD(0,1)-1;
        f(dir) = LOAD(0,2);
    }
    else if(b==LOAD(0,0)-1)
    {
        int dir = LOAD(0,1)-1;
        f(dir+4) = LOAD(0,2);
    }

    return f;
}



//Post Processing
//ApplyConstraints in FDM code
void ApplyConstraints_VAM(Eigen::SparseMatrix<double, Eigen::ColMajor> &J, Eigen::VectorXd &R, Eigen::MatrixXd CNODE, int NNODE)
{
    for(int j=0;j<CNODE.rows();j++)
    {
        int nodenum = CNODE(j,0);
        int dir = CNODE(j,1);
        for(int i=0;i<3*NNODE;i++)
            J.coeffRef(3*(nodenum-1)+dir,i) = 0;
        R(3*(nodenum-1)+dir)=0;
        J.coeffRef(3*(nodenum-1)+dir,3*(nodenum-1)+dir) = 1;
    }

}


Eigen::MatrixXd RotationMatrix(Eigen::VectorXd theta)
{
    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(3,3);

    double theta1 = theta(0);
    double theta2 = theta(1);
    double theta3 = theta(2);

    C(0,0) = 1 + (1.0/4)*pow(theta1,2);
    C(0,1) = theta3 + (1.0/2)*theta1*theta2;
    C(0,2) = -theta2 + (1.0/2)*theta1*theta3;
    C(1,0) = -theta3 + (1.0/2)*theta1*theta2;
    C(1,1) = 1 + (1.0/4)*pow(theta2,2);
    C(1,2) = theta1 + (1.0/2)*theta2*theta3;
    C(2,0) = theta2 + (1.0/2)*theta3*theta1;
    C(2,1) = -theta1 + (1.0/2)*theta2*theta3;
    C(2,2) = 1 + (1.0/4)*pow(theta3,2);

    return C;
}

/*
void StrainToDisplacement(Eigen::VectorXd &Gamma, Eigen::VectorXd &U_1, Eigen::VectorXd &U_2, Eigen::VectorXd &U_0,
                          double xa, double xb)
{
    Eigen::VectorXd theta(3);
    //theta = CurvatureToRotation(Gamma, U_1, U_2, U_0, xa, xb);

    //Calculate rotation matrix
    Eigen::MatrixXd rotmat(3,3);
    double theta1 = theta[0];
    double theta2 = theta[1];
    double theta3 = theta[2];
    double psi = 1+(1/4)*(pow(theta1,2)+pow(theta2,2)+pow(theta3,2));
    double phi = 1-(1/4)*(pow(theta1,2)+pow(theta2,2)+pow(theta3,2));
    Eigen::MatrixXd Delta = Eigen::MatrixXd::Zero(3,3);
    for(int i=0;i<3;i++)
        Delta(i,i) = phi;
    Eigen::MatrixXd Tildetheta(3,3);
    Tildetheta(0,0) = 0;
    Tildetheta(0,1) = -theta3;
    Tildetheta(0,2) = theta2;
    Tildetheta(1,0) = theta3;
    Tildetheta(1,1) = 0;
    Tildetheta(1,2) = -theta1;
    Tildetheta(2,0) = -theta2;
    Tildetheta(2,1) = theta1;
    Tildetheta(2,2) = 0;
    rotmat = (Delta - Tildetheta + (1/2)*theta*theta.transpose())/psi;

    //Solve for u
    double gamma11 = Gamma(0);
    double gamma12 = Gamma(1)/2;
    double gamma13 = Gamma(2)/2;

    double u1 = U_1(0);
    double u2 = U_1(1);
    double u3 = U_1(2);

    double u_prev1 = U_0(0);
    double u_prev2 = U_0(1);
    double u_prev3 = U_0(2);

    double u_next1 = U_2(0);
    double u_next2 = U_2(1);
    double u_next3 = U_2(2);

    double h = xb-xa;
    double C11 = rotmat(0,0);
    double C12 = rotmat(0,1);
    double C13 = rotmat(0,3);
    double C21 = rotmat(1,0);
    double C22 = rotmat(1,1);
    double C23 = rotmat(1,2);
    double C31 = rotmat(2,0);
    double C32 = rotmat(2,1);
    double C33 = rotmat(2,2);

    double u_new1 = ((gamma11 + 1 - C11)*h - C12*(u2-u_prev2)-C13*(u3-u_prev3)+C11*u_prev1)/C11;
    double u_new2 = ((2*gamma12 - C21)*h - C22*(u2-u_prev2)-C23*(u3-u_prev3)+C21*u_prev1)/C21;
}


PostProcessing CurvatureToRotation_Jacobian(Eigen::VectorXd Tau, Eigen::VectorXd U_1, Eigen::VectorXd U_2, Eigen::VectorXd U_0,
                                            double xa, double xb, Eigen::VectorXd k, int i, int NNODE)
{
    struct PostProcessing PP;
    //Kappa values
    double kappa1 = Tau(1);
    double kappa2 = Tau(2);
    double kappa3 = Tau(3);

    //values of node i
    double theta1, theta2, theta3;
    theta1 = U_1(0);
    theta2 = U_1(1);
    theta3 = U_1(2);

    //values of node i+1
    double Theta1, Theta2, Theta3;
    Theta1 = U_2(0);
    Theta2 = U_2(1);
    Theta3 = U_2(2);

    //values of node i-1
    double phi1, phi2, phi3;
    phi1 = U_0(0);
    phi2 = U_0(1);
    phi3 = U_0(2);

    double h = xb-xa;

    //Initial twist vector
    double k1 = k(0);
    double k2 = k(1);
    double k3 = k(2);

    //Jacobian Matrix
    PP.Jacobian = Eigen::MatrixXd::Zero(3,3);
    double Co = 1/(1+(1/4)*(pow(theta1,2)+pow(theta2,2)+pow(theta3,2)));
    double exi1=0, exi2=0, exi3=0, psi1=0, psi2=0, psi3=0;
    if(i>0&&i<NNODE-1)
    {
        exi1 = (Theta1 - phi1)/(2*h) + (Theta2 - phi2)*theta3/(2*h*2) - (Theta3 - phi3)*theta2/(2*h*2) ;
        exi2 = -(Theta1 - phi1)*theta3/(2*2*h) + (Theta2 - phi2)/(h*2) + (Theta3 - phi3)*theta1/(2*h*2);
        exi3 = (Theta1 - phi1)*theta2/(2*h*2) - (Theta2 - phi2)*theta1/(2*h*2) + (Theta3 - phi3)/(h*2);

        psi1 = (1 + 0.25*pow(theta1,2))*k1 + (theta3 + 0.5*theta1*theta2)*k2 + (0.5*theta1*theta3 - theta2)*k3;
        psi2 = (0.5*theta1*theta2 - theta3)*k1 + (1 + 0.5*pow(theta2,2))*k2 + (theta1 + 0.5*theta2*theta3)*k3;
        psi3 = (theta2 + 0.5*theta3*theta1)*k1 + (0.5*theta3*theta2 - theta1)*k2 + (1 + 0.25*pow(theta3,2))*k3;

        PP.Jacobian(0,0) = (1/pow(Co,2))*(-0.5*theta1*exi1) - (1/pow(Co,2))*(-0.5*theta1*psi1)
                + (1/Co)*(0.5*theta1*k1 + 0.5*theta2*k2 + 0.5*theta3*k3);
        PP.Jacobian(0,1) = (0.5/Co)*(Theta3 - phi3)/(2*h) - (0.5*theta2/pow(Co,2))*exi1 - (0.5*theta2/pow(Co,2))*psi1
                + (1/Co)*(0.5*theta1*k2 - k3);
        PP.Jacobian(0,2) = (0.5/Co)*(Theta2 - phi2)/(2*h) - (0.5*theta3/pow(Co,2))*exi1 - (0.5*theta3/pow(Co,2))*psi1
                + (1/Co)*(0.5*theta1*k3 + k2);
        PP.Jacobian(1,0) = (0.5/Co)*(Theta3 - phi3)/(2*h) - (0.5*theta1/pow(Co,2))*exi2 - (0.5*theta1/pow(Co,2))*psi2
                + (1/Co)*(0.5*theta2*k1 + k3);
        PP.Jacobian(1,1) = (1/pow(Co,2))*(-0.5*theta2*exi2) - (1/pow(Co,2))*(-0.5*theta2*psi2)
                + (1/Co)*(0.5*theta1*k1 + 0.5*theta2*k2 + 0.5*theta3*k3);
        PP.Jacobian(1,2) = -(0.5/Co)*(Theta1 - phi1)/(2*h) - (0.5*theta3/pow(Co,2))*exi2 - (0.5*theta3/pow(Co,2))*psi2
                + (1/Co)*(0.5*theta2*k3 - k1);
        PP.Jacobian(2,0) = (0.5/Co)*(Theta2 - phi2)/(2*h) - (0.5*theta1/pow(Co,2))*exi3 - (0.5*theta1/pow(Co,2))*psi3
                + (1/Co)*(-0.5*theta3*k1 - k2);
        PP.Jacobian(2,1) = (0.5/Co)*(Theta1 - phi1)/(2*h) - (0.5*theta2/pow(Co,2))*exi3 - (0.5*theta2/pow(Co,2))*psi3
                + (1/Co)*(0.5*theta3*k2 - k1);
        PP.Jacobian(2,2) = (1/pow(Co,2))*(-0.5*theta3*exi3) - (1/pow(Co,2))*(-0.5*theta3*psi3)
                + (1/Co)*(0.5*theta1*k1 + 0.5*theta2*k2 + 0.5*theta3*k3);
    }
    else if(i==0)
    {
        exi1 = (Theta1 - theta1)/h + (Theta2 - theta2)*theta3/(h*2) - (Theta3 - theta3)*theta2/(h*2) ;
        exi2 = -(Theta1 - theta1)*theta3/(2*h) + (Theta2 - theta2)/h + (Theta3 - theta3)*theta1/(h*2);
        exi3 = (Theta1 - theta1)*theta2/(h*2) - (Theta2 - theta2)*theta1/(h*2) + (Theta3 - theta3)/h;

        psi1 = (1 + 0.25*pow(theta1,2))*k1 + (theta3 + 0.5*theta1*theta2)*k2 + (0.5*theta1*theta3 - theta2)*k3;
        psi2 = (0.5*theta1*theta2 - theta3)*k1 + (1 + 0.5*pow(theta2,2))*k2 + (theta1 + 0.5*theta2*theta3)*k3;
        psi3 = (theta2 + 0.5*theta3*theta1)*k1 + (0.5*theta3*theta2 - theta1)*k2 + (1 + 0.25*pow(theta3,2))*k3;

        PP.Jacobian(0,0) = (1/pow(Co,2))*(-0.5*theta1*exi1) - (1/pow(Co,2))*(-0.5*theta1*psi1)
                + (1/Co)*(0.5*theta1*k1 + 0.5*theta2*k2 + 0.5*theta3*k3);
        PP.Jacobian(0,1) = (0.5/Co)*(Theta3 - theta3)/h - (0.5*theta2/pow(Co,2))*exi1 - (0.5*theta2/pow(Co,2))*psi1
                + (1/Co)*(0.5*theta1*k2 - k3);
        PP.Jacobian(0,2) = (0.5/Co)*(Theta2 - theta2)/h - (0.5*theta3/pow(Co,2))*exi1 - (0.5*theta3/pow(Co,2))*psi1
                + (1/Co)*(0.5*theta1*k3 + k2);
        PP.Jacobian(1,0) = (0.5/Co)*(Theta3 - theta3)/h - (0.5*theta1/pow(Co,2))*exi2 - (0.5*theta1/pow(Co,2))*psi2
                + (1/Co)*(0.5*theta2*k1 + k3);
        PP.Jacobian(1,1) = (1/pow(Co,2))*(-0.5*theta2*exi2) - (1/pow(Co,2))*(-0.5*theta2*psi2)
                + (1/Co)*(0.5*theta1*k1 + 0.5*theta2*k2 + 0.5*theta3*k3);
        PP.Jacobian(1,2) = -(0.5/Co)*(Theta1 - theta1)/h - (0.5*theta3/pow(Co,2))*exi2 - (0.5*theta3/pow(Co,2))*psi2
                + (1/Co)*(0.5*theta2*k3 - k1);
        PP.Jacobian(2,0) = (0.5/Co)*(Theta2 - theta2)/h - (0.5*theta1/pow(Co,2))*exi3 - (0.5*theta1/pow(Co,2))*psi3
                + (1/Co)*(-0.5*theta3*k1 - k2);
        PP.Jacobian(2,1) = (0.5/Co)*(Theta1 - theta1)/h - (0.5*theta2/pow(Co,2))*exi3 - (0.5*theta2/pow(Co,2))*psi3
                + (1/Co)*(0.5*theta3*k2 - k1);
        PP.Jacobian(2,2) = (1/pow(Co,2))*(-0.5*theta3*exi3) - (1/pow(Co,2))*(-0.5*theta3*psi3)
                + (1/Co)*(0.5*theta1*k1 + 0.5*theta2*k2 + 0.5*theta3*k3);
    }
    else if(i==NNODE-1)
    {
        exi1 = (theta1 - phi1)/h + (theta2 - phi2)*theta3/(h*2) - (theta3 - phi3)*theta2/(h*2) ;
        exi2 = -(theta1 - phi1)*theta3/(2*h) + (theta2 - phi2)/h + (theta3 - phi3)*theta1/(h*2);
        exi3 = (theta1 - phi1)*theta2/(h*2) - (theta2 - phi2)*theta1/(h*2) + (theta3 - phi3)/h;

        psi1 = (1 + 0.25*pow(theta1,2))*k1 + (theta3 + 0.5*theta1*theta2)*k2 + (0.5*theta1*theta3 - theta2)*k3;
        psi2 = (0.5*theta1*theta2 - theta3)*k1 + (1 + 0.5*pow(theta2,2))*k2 + (theta1 + 0.5*theta2*theta3)*k3;
        psi3 = (theta2 + 0.5*theta3*theta1)*k1 + (0.5*theta3*theta2 - theta1)*k2 + (1 + 0.25*pow(theta3,2))*k3;

        PP.Jacobian(0,0) = (1/pow(Co,2))*(-0.5*theta1*exi1) - (1/pow(Co,2))*(-0.5*theta1*psi1)
                + (1/Co)*(0.5*theta1*k1 + 0.5*theta2*k2 + 0.5*theta3*k3);
        PP.Jacobian(0,1) = (0.5/Co)*(theta3 - phi3)/h - (0.5*theta2/pow(Co,2))*exi1 - (0.5*theta2/pow(Co,2))*psi1
                + (1/Co)*(0.5*theta1*k2 - k3);
        PP.Jacobian(0,2) = (0.5/Co)*(theta2 - phi2)/h - (0.5*theta3/pow(Co,2))*exi1 - (0.5*theta3/pow(Co,2))*psi1
                + (1/Co)*(0.5*theta1*k3 + k2);
        PP.Jacobian(1,0) = (0.5/Co)*(theta3 - phi3)/h - (0.5*theta1/pow(Co,2))*exi2 - (0.5*theta1/pow(Co,2))*psi2
                + (1/Co)*(0.5*theta2*k1 + k3);
        PP.Jacobian(1,1) = (1/pow(Co,2))*(-0.5*theta2*exi2) - (1/pow(Co,2))*(-0.5*theta2*psi2)
                + (1/Co)*(0.5*theta1*k1 + 0.5*theta2*k2 + 0.5*theta3*k3);
        PP.Jacobian(1,2) = -(0.5/Co)*(theta1 - phi1)/h - (0.5*theta3/pow(Co,2))*exi2 - (0.5*theta3/pow(Co,2))*psi2
                + (1/Co)*(0.5*theta2*k3 - k1);
        PP.Jacobian(2,0) = (0.5/Co)*(theta2 - phi2)/h - (0.5*theta1/pow(Co,2))*exi3 - (0.5*theta1/pow(Co,2))*psi3
                + (1/Co)*(-0.5*theta3*k1 - k2);
        PP.Jacobian(2,1) = (0.5/Co)*(theta1 - phi1)/h - (0.5*theta2/pow(Co,2))*exi3 - (0.5*theta2/pow(Co,2))*psi3
                + (1/Co)*(0.5*theta3*k2 - k1);
        PP.Jacobian(2,2) = (1/pow(Co,2))*(-0.5*theta3*exi3) - (1/pow(Co,2))*(-0.5*theta3*psi3)
                + (1/Co)*(0.5*theta1*k1 + 0.5*theta2*k2 + 0.5*theta3*k3);
    }
    double f1 = exi1*(1/Co) + psi1*(1/Co) - k1 - kappa1;
    double f2 = exi2*(1/Co) + psi2*(1/Co) - k2 - kappa2;
    double f3 = exi3*(1/Co) + psi3*(1/Co) - k3 - kappa3;

    PP.R = Eigen::VectorXd::Zero(3);
    PP.R(0) = f1;
    PP.R(1) = f2;
    PP.R(2) = f3;

    return PP;
}

//This refines the mesh and gives back the nodal and connectivity information for refined mesh.
//How to change the loading conditions for refined mesh?
//How to change the boundary conditions for refined mesh?
VAMBeamElement RefineMesh(int refine, VAMBeamElement M)
{
    VAMBeamElement M2;

    int NELEM = M.NELEM*refine;
    M2.NODE = Eigen::MatrixXd(NELEM+1,3);
    M2.ELEM = Eigen::MatrixXd(NELEM,3);
    for(int i=0;i<M.NELEM;i++)
    {
        int a = M.ELEM(i,1)-1;
        int b = M.ELEM(i,2)-1;
        double xa = M.NODE((int)a,0);
        double xb = M.NODE((int)b,0);

        for(int k=0;k<refine;k++)
        {
            if(k==0)
            {
                M2.ELEM(i+k,1) = a + 1;
                M2.ELEM(i+k,2) = M.NNODE + i + k;
                M2.NODE((int)M2.ELEM(i+k,1)-1,0) = xa;
                M2.NODE((int)M2.ELEM(i+k,2)-1,1) = xb/refine;
            }
            else
            {
                M2.ELEM(i+k,1) = M2.ELEM(i+k-1,2);
                M2.ELEM(i+k,2) = M.NNODE + i + k;
                M2.NODE((int)M2.ELEM(i+k,1)-1,0) = M2.NODE((int)M2.ELEM(i+k-1)-1,1);
                M2.NODE((int)M2.ELEM(i+k,2)-1,1) = xb*k/refine;
            }
        }


    }
    return M2;
}

//Newton raphson to find rotations from curvatures
Eigen::VectorXd FDM(Eigen::VectorXd &Tau, int NNODE, int NELEM, Eigen::VectorXd k, VAMBeamElement M)
{
    Eigen::VectorXd Theta(NNODE*3);
    Eigen::VectorXd Theta_new = Eigen::VectorXd::Zero(NNODE*3);
    Eigen::VectorXd dTheta = Eigen::VectorXd::Zero(NNODE*3);
    Eigen::VectorXd error = Eigen::VectorXd::Zero(NNODE*3);
    Eigen::VectorXd Theta_0(4), Theta_1(4), Theta_2(4), Strain(4);
    Eigen::SparseMatrix<double, Eigen::ColMajor> J(NNODE*3,NNODE*3);
    Eigen::VectorXd R = Eigen::VectorXd::Zero(NNODE*3);
    struct PostProcessing PP;
    double max;
    int iter = 0;
    int maxiter = 100;
    for(int i=0;i<NNODE*3;i++)
        Theta(i) = 1;
    //Newton Raphson
    do
    {
        Theta_0.setZero();
        Theta_1.setZero();
        Theta_2.setZero();
        for(int i=0;i<NNODE;i++)
        {
            //int a = M.ELEM(i,1)-1;
            //int b = M.ELEM(i,2)-1;
            double xa = M.NODE(0,0);
            double xb = M.NODE(1,0);


            for(int j=0;j<3;j++)
            {
                if(i>0&&i<NNODE-1)
                {
                    Theta_0(j) = Theta(3*(i-1)+j);
                    Theta_1(j) = Theta(3*i+j);
                    Theta_2(j) = Theta(3*(i+1)+j);
                    Strain(j) = Tau(4*i+j);
                }
                else if(i==0)
                {
                    Theta_1(j) = Theta(3*i+j);
                    Theta_2(j) = Theta(3*(i+1)+j);
                    Strain(j) = Tau(4*i+j);
                }
                else if(i==NNODE-1)
                {
                    Theta_0(j) = Theta(3*(i-1)+j);
                    Theta_1(j) = Theta(3*i+j);
                    Strain(j) = Tau(4*i+j);
                }

            }

            //Get the jacobian and residual for each node
            PP = CurvatureToRotation_Jacobian(Strain,Theta_1,Theta_2,Theta_0,xa,xb,k,i,NNODE);

            //Assembly
            for(int j=0;j<3;j++)
            {
                for(int k=0;k<3;k++)
                    J.coeffRef(3*i+j,3*i+k) += PP.Jacobian(j,k);
                R(3*i+j) += PP.R(j);
            }

        }

        //ApplyConstraints
        ApplyConstraints_VAM(J,R,M.CNODE,NNODE);

//        for(int i=0;i<3*NNODE;i++)
//            if(J.coeffRef(i,i)==0)
//                J.coeffRef(i,i) = pow(10,6);
//        std::cout<<J<<std::endl;

        Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>> solver;
        solver.analyzePattern(J);
        solver.factorize(J); //LU decomposition

        assert(solver.info()==Eigen::Success);

        dTheta = solver.solve(-R);

        for(int j=0;j<NNODE*3;j++)
            Theta_new(j) = Theta(j) + dTheta(j);

        //Error calculation
        for(int j=0;j<NNODE;j++)
            for(int k=0;k<3;k++)
                error(k+j) = abs((Theta_new(k+j)-Theta(k+j))/Theta(k+j));
        max = error(0);
        for(int j=0;j<M.NNODE*3;j++)
            if(max<error(j))
                max = error(j);

        //Assignment for next iteration
        for(int j=0;j<M.NNODE;j++)
            for(int k=0;k<3;k++)
                Theta(j+k) = Theta_new(j+k);
        iter++;
    }while(max>pow(10,-5)&&iter<maxiter);

    if(iter<maxiter)
        return Theta;
    else
    {
        std::cout<<"Solution did not converge"<<std::endl;
        std::exit(0);
    }
}
*/
#endif

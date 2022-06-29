#include<Eigen/Sparse>
#include<Eigen/Dense>
#include "variables.h"
#include<iostream>
#include<math.h>

//Jacobian
double Jacobian(double xa, double xb)
{
    double J = (xb-xa)/2;
    assert(J>0);
    return J;
}

//Hermite Cubics
double HC(double exi, int num, double xa, double xb)
{
    double h = 2;
    double N1_x = 4-6*exi+2*pow(exi,3);
    double N2_x = (h*(pow(exi,2)-1)*(exi-1))*Jacobian(xa,xb);
    double N3_x = 4+6*exi-2*pow(exi,3);
    double N4_x = h*(pow(exi,2)-1)*(exi+1)*Jacobian(xa,xb);
    if(num==1)
        return N1_x/8.0;
    else if(num==2)
        return N2_x/8.0;
    else if (num==3)
        return N3_x/8.0;
    else if (num==4)
        return  N4_x/8.0;
    else
        return -99;
}

//Function to calculate double derivative of Hermite Cubic Functions
double DDHC(double exi, int num, double xa, double xb)
{
    double h = 2;
    if(num==1)
        return 6*exi/pow(h,2);
    else if(num==2)
        return (3*exi-1)*Jacobian(xa,xb)/h;
    else if(num==3)
        return -6*exi/pow(h,2);
    else if(num==4)
        return (3*exi+1)*Jacobian(xa,xb)/h;
    else
        return -99;
}

//Function to calculate single derivative of Hermite Cubic functions
double SDHC(double exi, int num, double xa, double xb)
{
    double h = 2;
    if(num==1)
        return 6*(pow(exi,2)-1)/(4*h);
    else if(num==2)
        return (3*pow(exi,2)-2*exi-1)*Jacobian(xa,xb)/4;
    else if(num==3)
        return 6*(1-pow(exi,2))/(4*h);
    else if(num==4)
        return (3*pow(exi,2)+2*exi-1)*Jacobian(xa,xb)/4;
    else
        return -99;
}

//Hat functions
double HF(double exi, double num)
{
    if(num==1)
        return (1-exi)/2;
    else if(num==2)
        return (1+exi)/2;
    else
        return -99;
}

//Function to evaluate single derivative of hat functions
double SDHF(double num)
{
    double h = 2;
    if(num==1)
        return -1/h;
    else if(num==2)
        return 1/h;
    else
        return -99;
}

//To rearrange the elements of stiffness matrix and force vector
void RearrangeElementMatrices(Eigen::MatrixXd &k, Eigen::MatrixXd &t, Eigen::VectorXd &f)
{
    Eigen::MatrixXd temp(6,6), temp2(6,6);
    Eigen::VectorXd temp3(6);
    //Rearranging rows
    for(int j=0;j<6;j++)
    {
        //Stiffness
        temp(0,j) = k(0,j);
        temp(1,j) = k(2,j);
        temp(2,j) = k(3,j);
        temp(3,j) = k(1,j);
        temp(4,j) = k(4,j);
        temp(5,j) = k(5,j);
        //Tangent
        temp2(0,j) = t(0,j);
        temp2(1,j) = t(2,j);
        temp2(2,j) = t(3,j);
        temp2(3,j) = t(1,j);
        temp2(4,j) = t(4,j);
        temp2(5,j) = t(5,j);
        //Force
        temp3(j) = f(j);
    }
//    Rearranging columns
    for(int i=0;i<6;i++)
    {
        //Stiffness
        k(i,0) = temp(i,0) ;
        k(i,1) = temp(i,2) ;
        k(i,2) = temp(i,3) ;
        k(i,3) = temp(i,1) ;
        k(i,4) = temp(i,4) ;
        k(i,5) = temp(i,5) ;
        //Tangent
        t(i,0) = temp2(i,0) ;
        t(i,1) = temp2(i,2) ;
        t(i,2) = temp2(i,3) ;
        t(i,3) = temp2(i,1) ;
        t(i,4) = temp2(i,4) ;
        t(i,5) = temp2(i,5) ;
    }
    f(0) = temp3(0);
    f(1) = temp3(2);
    f(2) = temp3(3);
    f(3) = temp3(1);
    f(4) = temp3(4);
    f(5) = temp3(5);
}



//How to calculate dw0/dx from previous iteration
double dw0dx(Eigen::VectorXd &U, double x, int a, int b,double xa,double xb)
{
    double dwdx;
    double wa = U(3*a+1);
    double wb = U(3*b+1);
    double theta_a = U(3*a+2);
    double theta_b = U(3*b+2);
    dwdx = wa*SDHC(x,1,xa,xb)+theta_a*SDHC(x,2,xa,xb)+wb*SDHC(x,3,xa,xb)+theta_b*SDHC(x,4,xa,xb);
    return dwdx;
}

double dU0dx(Eigen::VectorXd &U, int a, int b)
{
    double dudx;
    double ua = U(3*a);
    double ub = U(3*b);
    dudx = ua*SDHF(1)+ub*SDHF(2);
    return dudx;
}

void TangentStiffnessEulerBernoulli(Eigen::MatrixXd &t, Eigen::MatrixXd &k, double xa, double xb,double E, double nu, double base, double height, Eigen::VectorXd& U, int a, int b)
{
    //GaussPoints
    double x1 = -sqrt(3.0/5);
    double x2 = 0;
    double x3 = sqrt(3.0/5);
    double w1 = 5.0/9;
    double w2 = 8.0/9;
    double w3 = 5.0/9;
    //Calculate Axx, Dxx
    double Axx, Dxx;
    Axx = E*base*height;
    Dxx = E*base*pow(height,3)/12;

    //T11
    for(int i=0;i<2;i++)
        for(int j=0;j<2;j++)
            t(i,j)=k(i,j);
    //T12
    for(int i=0;i<2;i++)
        for(int j=2;j<6;j++)
            t(i,j)=2*k(i,j);
    //T21
    for(int i=2;i<6;i++)
        for(int j=0;j<2;j++)
            t(i,j) = k(i,j);
    //T22
    for(int i=2;i<6;i++)
        for(int j=2;j<6;j++)
        {
            //Integration points
            //This is used to indicate the map between x and exi
            //A subparametric map (in this case a linear map) is used to convert x to exi.
            double fxa = Axx*(dU0dx(U,a,b)+pow(dw0dx(U,x1,a,b,xa,xb),2))*SDHC(x1,i-1,xa,xb)*SDHC(x1,j-1,xa,xb);
            double fxb = Axx*(dU0dx(U,a,b)+pow(dw0dx(U,x2,a,b,xa,xb),2))*SDHC(x2,i-1,xa,xb)*SDHC(x2,j-1,xa,xb);
            double fxc = Axx*(dU0dx(U,a,b)+pow(dw0dx(U,x3,a,b,xa,xb),2))*SDHC(x3,i-1,xa,xb)*SDHC(x3,j-1,xa,xb);
            t(i,j) = k(i,j) + w1*fxa+w2*fxb+w3*fxc;
        }
}


void StiffnessEulerBernoulli(Eigen::MatrixXd &k, double xa, double xb,double E, double nu, double base, double height, Eigen::VectorXd &U, int a, int b)
{
    //GaussPoints
    double x1 = -sqrt(3.0/5);
    double x2 = 0;
    double x3 = sqrt(3.0/5);
    double w1 = 5.0/9;
    double w2 = 8.0/9;
    double w3 = 5.0/9;
    //Calculate Axx, Dxx
    double Axx, Dxx;
    Axx = E*base*height;
    Dxx = E*base*pow(height,3)/12;

    //Stiffness matrix elements
    //K11
    for(int i=0;i<2;i++)
        for(int j=0;j<2;j++)
            //constant function
            //Exact integral can be evaluated since Axx is a constant
            k(i,j) =  2*Axx*SDHF(i+1)*SDHF(j+1)/Jacobian(xa,xb);
    //K12
    for(int i=0;i<2;i++)
        for(int j=2;j<6;j++)
        {
            //Quadrature
            double fxa = 0.5*Axx*dw0dx(U,x1,a,b,xa,xb)*SDHF(i+1)*SDHC(x1,j-1,xa,xb)/pow(Jacobian(xa,xb),2);
            double fxb = 0.5*Axx*dw0dx(U,x2,a,b,xa,xb)*SDHF(i+1)*SDHC(x2,j-1,xa,xb)/pow(Jacobian(xa,xb),2);
            double fxc = 0.5*Axx*dw0dx(U,x3,a,b,xa,xb)*SDHF(i+1)*SDHC(x3,j-1,xa,xb)/pow(Jacobian(xa,xb),2);
            k(i,j) = fxa*w1+fxb*w2+fxc*w3;
        }
    //K21
    for(int i=2;i<6;i++)
        for(int j=0;j<2;j++)
            k(i,j) = 2*k(j,i);
    //K22
    for(int i=2;i<6;i++)
        for(int j=2;j<6;j++)
        {
            double gxa = Dxx*DDHC(x1,i-1,xa,xb)*DDHC(x1,j-1,xa,xb)/pow(Jacobian(xa,xb),3);
            double gxb = Dxx*DDHC(x2,i-1,xa,xb)*DDHC(x2,j-1,xa,xb)/pow(Jacobian(xa,xb),3);
            double gxc = Dxx*DDHC(x3,i-1,xa,xb)*DDHC(x3,j-1,xa,xb)/pow(Jacobian(xa,xb),3);
            double fxa = 0.5*Axx*pow(dw0dx(U,x1,a,b,xa,xb),2)*SDHC(x1,i-1,xa,xb)*SDHC(x1,j-1,xa,xb)/pow(Jacobian(xa,xb),3);
            double fxb = 0.5*Axx*pow(dw0dx(U,x2,a,b,xa,xb),2)*SDHC(x2,i-1,xa,xb)*SDHC(x2,j-1,xa,xb)/pow(Jacobian(xa,xb),3);
            double fxc = 0.5*Axx*pow(dw0dx(U,x3,a,b,xa,xb),2)*SDHC(x3,i-1,xa,xb)*SDHC(x3,j-1,xa,xb)/pow(Jacobian(xa,xb),3);
            k(i,j) = (gxa+fxa)*w1+(gxb+fxb)*w2+(gxc+fxc)*w3;
        }
    //TangentStiffnessEulerBernoulli(t,k,xa,xb,E,nu,base,height,U,a,b);
}

//Local Force Vector
void ForceVec(Eigen::VectorXd &f, double xa, double xb, int a, int b, double vf, double af)
{
    //GaussPoints
    double x1 = -sqrt(3.0/5);
    double x2 = 0;
    double x3 = sqrt(3.0/5);
    double w1 = 5.0/9;
    double w2 = 8.0/9;
    double w3 = 5.0/9;
    //Two things have to be done. Input axial load and transverse load in ReadInp File and use it in this function
    //For now q is taken as 1
    //3 point gaussian quadrature.
    //Axial load acting
    for(int i=0;i<2;i++)
    {
//        std::cout<<x1<<std::endl;
//        std::cout<<x2<<std::endl;
//        std::cout<<x3<<std::endl;
//        std::cout<<w1<<std::endl;
//        std::cout<<w2<<std::endl;
//        std::cout<<w3<<std::endl;
        double fxa = af*HF(x1,i-1);
        double fxb = af*HF(x2,i-1);
        double fxc = af*HF(x3,i-1);
        f(i) = (w1*fxa + w2*fxb + w3*fxc)*Jacobian(xa,xb);
    }

    //Transverse load
    for(int i=2;i<6;i++)
    {
        double fxa = vf*HC(x1,i-1,xa,xb);
        double fxb = vf*HC(x2,i-1,xa,xb);
        double fxc = vf*HC(x3,i-1,xa,xb);
        f(i) = (w1*fxa + w2*fxb + w3*fxc)*Jacobian(xa,xb);
    }
}

void ApplyConstraints2DEB(Eigen::SparseMatrix<double, Eigen::ColMajor>& T, Eigen::VectorXd& U, Eigen::MatrixXd &CNODE, int NNODE, Eigen::VectorXd& R)
{
    for(int j=0;j<CNODE.rows();j++)
    {
        int nodenum = CNODE(j,0);
        int dir = CNODE(j,1);
        for(int i=0;i<3*NNODE;i++)
            T.coeffRef(3*(nodenum-1)+dir,i) = 0;
        R(3*(nodenum-1)+dir)=0;
        T.coeffRef(3*(nodenum-1)+dir,3*(nodenum-1)+dir) = 1;
    }
}

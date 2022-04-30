#include<Eigen/Sparse>
#include<Eigen/Dense>

//Hermite Cubics
double HC(double xa, double xb, double x, int num)
{
    double h = (xb-xa);
    double N1_x = -pow((x-xb),2)*(-h+2*(xa - x))/pow(h,3);
    double N2_x = (x - xa)*pow((x-xb),2)/pow(h,2);
    double N3_x = pow((x-xb),2)*(h + 2*(xb - x))/pow(h,3);
    double N4_x = pow((x-xb),2)*(x - xb)/pow(h,2);
    if(num==1)
        return N1_x;
    else if(num==2)
        return N2_x;
    else if (num==3)
        return N3_x;
    else if (num==4)
        return  N4_x;
    else
        return -99;
}

//Function to calculate double derivative of Hermite Cubic Functions
double DDHC(double xa, double xb, double x, int num)
{
    double h = (xb-xa);
    if(num==1)
        return -2*(-h+2*(xa-x))/pow(h,3);
    else if(num==2)
        return (4*(x-xb)+2*(x-xa))/pow(h,2);
    else if(num==3)
        return (2*(h+2*(xb-x))-6*(x-xa))/pow(h,3);
    else if(num==4)
        return (4*(x-xa)+2*(x-xb))/pow(h,2);
    else
        return -99;
}

//Function to calculate single derivative of Hermite Cubic functions
double SDHC(double xa, double xb, double x, int num)
{
    double h = (xb-xa);
    if(num==1)
        return -2*(x-xb)*(-h + 2*(xa-x)) + 2*pow((x-xb),2)/pow(h,3);
    else if(num==2)
        return (pow(x-xb,2) + 2*(x-xb)*(x-xa))/pow(h,2);
    else if(num==3)
        return (2*(x-xa)*(h+2*(xb-x)) - 2*pow((x-xa),2))/pow(h,3);
    else if(num==4)
        return (2*(x-xb)*(x-xa) + pow(x-xa, 2))/pow(h,2);
    else
        return -99;
}

//Hat functions
double HF(double xa, double xb, double x, double num)
{
    double h = (xb-xa);
    if(num==2)
        return (xb-x)/h;
    else if(num==1)
        return (x-xa)/h;
    else
        return -99;
}

//Function to evaluate single derivative of hat functions
double SDHF(double xa, double xb, double num)
{
    double h = (xb-xa);
    if(num==2)
        return -1/h;
    else if(num==1)
        return 1/h;
    else
        return -99;
}

//To rearrange the elements of stiffness matrix and force vector
void RearrangeElementStiffness(Eigen::MatrixXd &k)
{
    Eigen::MatrixXd temp(6,6);
    //Rearranging rows
    for(int j=0;j<6;j++)
    {
        temp(1,j) = k(2,j);
        temp(2,j) = k(3,j);
        temp(3,j) = k(1,j);
        temp(0,j) = k(0,j);
        temp(4,j) = k(4,j);
        temp(5,j) = k(5,j);
    }
//    Rearranging columns
    for(int i=0;i<6;i++)
    {
        k(i,1) = temp(i,2) ;
        k(i,2) = temp(i,3) ;
        k(i,3) = temp(i,1) ;
        k(i,0) = temp(i,0) ;
        k(i,4) = temp(i,4) ;
        k(i,5) = temp(i,5) ;
    }
}

//How to calculate dw0/dx from previous iteration
double dw0dx(Eigen::VectorXd U, double x, int a, int b)
{
    double dwdx;
    double wa = U(3*a+1);
    double wb = U(3*b+1);
    double theta_a = U(3*a+2);
    double theta_b = U(3*b+2);
    dwdx = wa*SDHC(-1,1,x,1)+theta_a*SDHC(-1,1,x,2)+wb*SDHC(-1,1,x,3)+theta_b*SDHC(-1,1,x,4);
    return dwdx;
}

double dU0dx(Eigen::VectorXd U, int a, int b)
{
    double dudx;
    double ua = U(3*a);
    double ub = U(3*b);
    dudx = ua*SDHF(-1,1,1)+ub*SDHF(-1,1,2);
    return dudx;
}

void TangentStiffnessEulerBernoulli(Eigen::MatrixXd &t, Eigen::MatrixXd &k, double xa, double xb,double E, double nu, double base, double height, Eigen::VectorXd U, int a, int b)
{
    double h = (xb-xa);
    //Weights for gaussian  quadrature
    double w1 = 0.555555556;
    double w2 = 0.888888889;
    double w3 = 0.555555556;
    //Function evaluation points for quadrature
    double x1 = - 0.774596669;
    double x2 = 0;
    double x3 = 0.774596669;
    //Calculate Axx, Dxx
    double Axx, Dxx;
    Axx = E*base*height;
    Dxx = E*base*pow(height,3)/12;

    //T11, T12
    for(int i=0;i<2;i++)
        for(int j=0;j<6;j++)
            t(i,j)=k(i,j);
    //T21
    for(int i=2;i<6;i++)
        for(int j=0;j<2;j++)
            t(i,j) = k(i,j);
    //T22
    for(int i=2;i<6;i++)
        for(int j=2;j<6;j++)
        {
            double fxa = Axx*(dU0dx(U,a,b)*dw0dx(U,x1,a,b)*dw0dx(U,x1,a,b))*SDHC(-1,1,x1,i-1)*SDHC(-1,1,x1,j-1);
            double fxb = Axx*(dU0dx(U,a,b)*dw0dx(U,x2,a,b)*dw0dx(U,x2,a,b))*SDHC(-1,1,x2,i-1)*SDHC(-1,1,x2,j-1);
            double fxc = Axx*(dU0dx(U,a,b)*dw0dx(U,x3,a,b)*dw0dx(U,x3,a,b))*SDHC(-1,1,x3,i-1)*SDHC(-1,1,x3,j-1);
            t(i,j) = k(i,j) + w1*fxa+w2*fxb+w3*fxc;
        }
}


void StiffnessEulerBernoulli(Eigen::MatrixXd &k, double xa, double xb,double E, double nu, double base, double height, Eigen::VectorXd U, int a, int b)
{
    double h = (xb-xa);
    //Weights for gaussian  quadrature
    double w1 = 0.555555556;
    double w2 = 0.888888889;
    double w3 = 0.555555556;
    //Function evaluation points for quadrature
    double x1 = - 0.774596669;
    double x2 = 0;
    double x3 = 0.774596669;
    //Calculate Axx, Dxx
    double Axx, Dxx;
    Axx = E*base*height;
    Dxx = E*base*pow(height,3)/12;

    //Stiffness matrix elements
    //K11
    for(int i=0;i<2;i++)
        for(int j=0;j<2;j++)
            //3 point Gaussian Quadrature
            k(i,j) =  2*(2/h)*Axx*SDHF(-1,1,i+1)*SDHF(-1,1,j+1);
    //K12
    for(int i=0;i<2;i++)
        for(int j=2;j<6;j++)
        {
            //Quadrature
            double fxa = 0.5*Axx*dw0dx(U,x1,a,b)*pow(2/h,2)*SDHF(-1, 1,i+1)*SDHC(-1,1,x1,j-1);
            double fxb = 0.5*Axx*dw0dx(U,x2,a,b)*pow(2/h,2)*SDHF(-1, 1,i+1)*SDHC(-1,1,x2,j-1);
            double fxc = 0.5*Axx*dw0dx(U,x3,a,b)*pow(2/h,2)*SDHF(-1, 1,i+1)*SDHC(-1,1,x3,j-1);
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
//            Quadrature
            double gxa = Dxx*pow(2/h,3)*DDHC(-1,1,x1,i-1)*DDHC(-1,1,x1,j-1);
            double gxb = Dxx*pow(2/h,3)*DDHC(-1,1,x2,i-1)*DDHC(-1,1,x2,j-1);
            double gxc = Dxx*pow(2/h,3)*DDHC(-1,1,x3,i-1)*DDHC(-1,1,x3,j-1);
            double fxa = 0.5*Axx*pow(dw0dx(U,x1,a,b),2)*SDHC(-1,1,x1,i-1)*SDHC(-1,1,x1,j-1);
            double fxb = 0.5*Axx*pow(dw0dx(U,x2,a,b),2)*SDHC(-1,1,x2,i-1)*SDHC(-1,1,x2,j-1);
            double fxc = 0.5*Axx*pow(dw0dx(U,x3,a,b),2)*SDHC(-1,1,x3,i-1)*SDHC(-1,1,x3,j-1);
            k(i,j) = (gxa+fxa)*w1+(gxb+fxb)*w2+(gxc+fxc)*w3;
        }
    //TangentStiffnessEulerBernoulli(t,k,xa,xb,E,nu,base,height,U,a,b);
    //RearrangeElementStiffness(k);
}


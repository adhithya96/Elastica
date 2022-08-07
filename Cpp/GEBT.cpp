//Functions required to solve GEBT
#ifndef VARIABLES_H
#define VARIABLES_H
#include"Variables.h"

Eigen::MatrixXd Rot_Mat(Eigen::VectorXd theta)
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

    return C/(1+(1.0/4)*(pow(theta1,2) + pow(theta2,2) + pow(theta3,2)));
}

Eigen::MatrixXd Transform_Mat(double theta)
{
    Eigen::MatrixXd C_ab = Eigen::MatrixXd::Zero(3,3);

    C_ab(0,0) = sin(theta);
    C_ab(1,0) = cos(theta);
    C_ab(0,1) = -cos(theta);
    C_ab(1,1) = sin(theta);
    C_ab(2,2) = 1;

    return C_ab;
}

Eigen::MatrixXd Dyadic(Eigen::VectorXd x)
{
    Eigen::MatrixXd x_tilde = Eigen::MatrixXd::Zero(3,3);

    double theta1 = x(0);
    double theta2 = x(1);
    double theta3 = x(2);

    x_tilde(0,1) = -theta3;
    x_tilde(0,2) = theta2;
    x_tilde(1,0) = theta3;
    x_tilde(1,2) = -theta1;
    x_tilde(2,0) = -theta2;
    x_tilde(2,1) = theta1;

    return x_tilde;
}

Eigen::MatrixXd dCdtheta(Eigen::VectorXd theta, double k)
{
    Eigen::MatrixXd C = Rot_Mat(theta);
    Eigen::VectorXd ek(3);
    for(int i=0;i<3;i++)
        ek(i) = 0;
    if(k==1)
        ek(0) = 1;
    else if(k==2)
        ek(1) = 1;
    else if(k==3)
        ek(2) = 3;
    Eigen::MatrixXd dCdt = Eigen::MatrixXd::Zero(3,3);
    dCdt = -(0.5*(ek.transpose()*theta*(Eigen::MatrixXd::Identity(3,3) + C) - Dyadic(ek) +
                  0.5*(ek*theta.transpose() + theta*ek.transpose())));
    return dCdt/(1 + (1.0/4)*theta.transpose()*theta);
}

//h - distance between two nodes in an element or mesh size
//U1 - Degrees of freedom of node i
//U2 - Degrees of freedom of node i+1
//F1 - Force and Moments on node i
//F2 - Force and Moments on node i+1
//External or Internal?
//S - Cross-Sectional Stiffness Matrix
//Theta -
//This function calculates the element residual vector needed for Newton Raphson
Eigen::VectorXd Element_Residual(Eigen::VectorXd U1, Eigen::VectorXd U2, Eigen::VectorXd F1, Eigen::VectorXd F2, double h,
                                 Eigen::MatrixXd S1, Eigen::MatrixXd S2, double Theta)
{
    Eigen::VectorXd F_1(3);
    F_1 << F1(0),F1(1),F1(2);
    Eigen::VectorXd F_2(3);
    F_2 << F2(0),F2(1),F2(2);
    Eigen::VectorXd M_1(3);
    M_1 << F1(3),F1(4),F1(5);
    Eigen::VectorXd M_2(3);
    M_2 << F2(3),F2(4),F2(5);

    Eigen::VectorXd u1(3);
    u1 << U1(0),U1(1),U1(2);
    Eigen::VectorXd u2(3);
    u2 << U2(0),U2(1),U2(2);
    Eigen::VectorXd theta1(3);
    theta1 << U1(3),U1(4),U1(5);
    Eigen::VectorXd theta2(3);
    theta2 << U2(3),U2(4),U2(5);

    Eigen::VectorXd Strain1(6), Strain2(6);
    Eigen::MatrixXd Sinv1(6,6), Sinv2(6,6);
    Sinv1 = S1.inverse();
    Sinv2 = S2.inverse();
    Strain1 = Sinv1*F1;
    Strain2 = Sinv2*F2;
    Eigen::VectorXd gamma1(3);
    gamma1 << Strain1(0),Strain1(1),Strain1(2);
    Eigen::VectorXd kappa1(3);
    kappa1 << Strain1(3),Strain1(4),Strain1(5);
    Eigen::VectorXd gamma2(3);
    gamma2 << Strain2(0),Strain2(1),Strain2(2);
    Eigen::VectorXd kappa2;
    kappa2 << Strain2(3),Strain2(4),Strain2(5);
    \
    Eigen::VectorXd e_1(3);
    e_1 << 1,0,0;

    //f_u
    Eigen::VectorXd fu_minus(3), fu_plus(3);
    //fu_i+1
    fu_minus = -(Rot_Mat(theta2).transpose())*Transform_Mat(Theta)*F_2;
    //fu_i
    fu_plus = (Rot_Mat(theta1).transpose())*Transform_Mat(Theta)*F_1;

    //f_psi
    Eigen::VectorXd fpsi_minus(3), fpsi_plus(3);
    //fpsi_i+1
    fpsi_minus = -(Rot_Mat(theta2).transpose())*Transform_Mat(Theta)*M_2 -
            (h/2)*(Rot_Mat(theta2).transpose())*Transform_Mat(Theta)*(Dyadic(e_1)+Dyadic(gamma2)*F_2);
    //fpsi_i
    fpsi_plus = (Rot_Mat(theta1).transpose())*Transform_Mat(Theta)*M_1 -
            (h/2)*(Rot_Mat(theta1).transpose())*Transform_Mat(Theta)*(Dyadic(e_1)+Dyadic(gamma1)*F_1);

    //f_F
    Eigen::VectorXd fF_minus(3), fF_plus(3);
    //fF_i+1
    fF_minus = u2 - (h/2)*(Rot_Mat(theta2)*Transform_Mat(Theta)*(Dyadic(e_1)+Dyadic(gamma2))-Transform_Mat(Theta)*e_1);
    //fF_i
    fF_plus = -u1 - (h/2)*(Rot_Mat(theta1)*Transform_Mat(Theta)*(Dyadic(e_1)+Dyadic(gamma1))-Transform_Mat(Theta)*e_1);

    //f_M
    Eigen::VectorXd fM_minus, fM_plus;
    //fM_i+1
    fM_minus = theta2 - (h/2)*(Eigen::MatrixXd::Identity(3,3) + Dyadic(theta2)/2 + (theta2*(theta2.transpose()))/4)
            *Transform_Mat(Theta)*kappa2;
    //fM_i
    fM_plus =  -theta1 - (h/2)*(Eigen::MatrixXd::Identity(3,3) + Dyadic(theta1)/2 + (theta1*(theta1.transpose()))/4)
            *Transform_Mat(Theta)*kappa1;

    Eigen::VectorXd R(12);
    for(int i=0;i<3;i++)
    {
        R(i) = fu_plus(i) + fu_minus(i);
        R(i+3) = fpsi_plus(i) + fpsi_minus(i);
        R(i+6) = fF_plus(i) + fF_minus(i);
        R(i+9) = fM_plus(i) + fM_minus(i);
    }

    return R;
}


Eigen::VectorXd Element_Residual(Eigen::VectorXd U1, Eigen::VectorXd F1, double h,
                                 Eigen::MatrixXd S1, double Theta, int nodenum, Boundary B)
{
    if(nodenum==0)
    {
        Eigen::VectorXd F_1(3);
        F_1 << F1(0),F1(1),F1(2);
        Eigen::VectorXd M_1(3);
        M_1 << F1(3),F1(4),F1(5);
        Eigen::VectorXd u1(3);
        u1 << U1(0),U1(1),U1(2);
        Eigen::VectorXd theta1(3);
        theta1 << U1(3), U1(4), U1(5);

        Eigen::VectorXd Strain1(6);
        Eigen::MatrixXd Sinv(6,6);
        Sinv = S1.inverse();
        Strain1 = Sinv*F1;

        Eigen::VectorXd gamma1(3);
        gamma1 << Strain1(0),Strain1(1),Strain1(2);
        Eigen::VectorXd kappa1;
        kappa1 << Strain1(3),Strain1(4),Strain1(5);
        \
        Eigen::VectorXd e_1(3);
        e_1 << 1,0,0;

        //f_u
        Eigen::VectorXd fu_minus(3);
        //fu_i+1
        fu_minus = -(Rot_Mat(theta1).transpose())*Transform_Mat(Theta)*F_1;

        //f_psi
        Eigen::VectorXd fpsi_minus(3);
        //fpsi_i+1
        fpsi_minus = -(Rot_Mat(theta1).transpose())*Transform_Mat(Theta)*M_1 -
                (h/2)*(Rot_Mat(theta1).transpose())*Transform_Mat(Theta)*(Dyadic(e_1)+Dyadic(gamma1)*F_1);

        //f_F
        Eigen::VectorXd fF_minus(3);
        //fF_i+1
        fF_minus = u1 - (h/2)*(Rot_Mat(theta1)*Transform_Mat(Theta)*(Dyadic(e_1)+Dyadic(gamma1))-Transform_Mat(Theta)*e_1);

        //f_M
        Eigen::VectorXd fM_minus;
        //fM_i+1
        fM_minus = theta1 - (h/2)*(Eigen::MatrixXd::Identity(3,3) + Dyadic(theta1)/2 + (theta1*(theta1.transpose()))/4)
                *Transform_Mat(Theta)*kappa1;

        Eigen::VectorXd R(12);
        for(int i=0;i<3;i++)
        {
            R(i) = fu_minus(i) - B.F1(i);
            R(i+3) = fpsi_minus(i) - B.M1(i);
            R(i+6) = fF_minus(i) - B.u1(i);
            R(i+9) = fM_minus(i) - B.theta1(i);
        }

        return R;
    }
    else
    {
        Eigen::VectorXd F_1(3);
        F_1 << F1(0),F1(1),F1(2);
        Eigen::VectorXd M_1(3);
        M_1 << F1(3),F1(4),F1(5);
        Eigen::VectorXd u1(3);
        u1 << U1(0),U1(1),U1(2);
        Eigen::VectorXd theta1(3);
        theta1 << U1(3),U1(4),U1(5);

        Eigen::VectorXd Strain1(6), Strain2(6);
        Eigen::MatrixXd Sinv(6,6);
        Sinv = S1.inverse();
        Strain1 = Sinv*F1;
        Eigen::VectorXd gamma1(3);
        gamma1 << Strain1(0),Strain1(1),Strain1(2);
        Eigen::VectorXd kappa1(3);
        kappa1 << Strain1(3),Strain1(4),Strain1(5);

        Eigen::VectorXd e_1(3);
        e_1 << 1,0,0;

        //f_u
        Eigen::VectorXd fu_plus(3);
        //fu_i
        fu_plus = (Rot_Mat(theta1).transpose())*Transform_Mat(Theta)*F_1;

        //f_psi
        Eigen::VectorXd fpsi_plus(3);
        //fpsi_i
        fpsi_plus = (Rot_Mat(theta1).transpose())*Transform_Mat(Theta)*M_1 -
                (h/2)*(Rot_Mat(theta1).transpose())*Transform_Mat(Theta)*(Dyadic(e_1)+Dyadic(gamma1)*F_1);

        //f_F
        Eigen::VectorXd fF_plus(3);
        //fF_i
        fF_plus = -u1 - (h/2)*(Rot_Mat(theta1)*Transform_Mat(Theta)*(Dyadic(e_1)+Dyadic(gamma1))-Transform_Mat(Theta)*e_1);

        //f_M
        Eigen::VectorXd fM_minus, fM_plus;
        //fM_i
        fM_plus =  -theta1 - (h/2)*(Eigen::MatrixXd::Identity(3,3) + Dyadic(theta1)/2 + (theta1*(theta1.transpose()))/4)
                *Transform_Mat(Theta)*kappa1;

        Eigen::VectorXd R(12);
        for(int i=0;i<3;i++)
        {
            R(i) = fu_plus(i) - B.FN(i);
            R(i+3) = fpsi_plus(i) - B.MN(i);
            R(i+6) = fF_plus(i) + B.uN(i);
            R(i+9) = fM_plus(i) + B.thetaN(i);
        }

        return R;
    }
}


//h - distance between two nodes in an element or mesh size
//U1 - Degrees of freedom of node i
//U2 - Degrees of freedom of node i+1
//F1 - Force and Moments on node i
//F2 - Force and Moments on node i+1
//External or Internal?
//S - Cross-Sectional Stiffness Matrix
//Theta -
//This function calculates the element Jacobian Matrix needed for Newton Raphson
Eigen::MatrixXd Element_Jacobian(Eigen::VectorXd U1, Eigen::VectorXd U2, Eigen::VectorXd F1, Eigen::VectorXd F2, double h,
                                 Eigen::MatrixXd S1, Eigen::MatrixXd S2, double Theta)
{
    //This portion of the code calculates the jacobian or tangent matrix to solve the newton raphson
    //Refer Equation 47 of "GEBT: A general-purpose nonlinear analysis tool for composite beams"
    //This is section is specfically for intermediate nodes. Boundary nodes will be dealt with seperately.
    //The equation contains a total of 12 unknowns. (u,theta,F,M) for nodes i and i+1
    //The jacobian matrix will be a 12*24 matrix
    Eigen::VectorXd F_1(3), M_1(3), u1(3), theta1(3), Strain1(6), gamma1(3), kappa1(3), e_1(3);
    Eigen::MatrixXd dgammadF(3,3), dgammadM(3,3), dkappadF(3,3), dkappadM(3,3);
    Eigen::MatrixXd dfu_dFi(3,3), dfu_dui(3,3), dfu_dthetai(3,3), dfu_dMi(3,3);
    Eigen::MatrixXd dfpsi_dFi(3,3), dfpsi_dui(3,3), dfpsi_dthetai(3,3), dfpsi_dMi(3,3);
    Eigen::VectorXd dfF_dFi(3), dfF_dui(3), dfF_dthetai(3), dfF_dMi(3);
    Eigen::VectorXd dfM_dFi(3), dfM_dui(3), dfM_dthetai(3), dfM_dMi(3);

    Eigen::VectorXd F_2(3), M_2(3), u2, theta2, Strain2(6), gamma2(3), kappa2(3);
    Eigen::MatrixXd dfu_dFj(3,3), dfu_duj(3,3), dfu_dthetaj(3,3), dfu_dMj(3,3);
    Eigen::MatrixXd dfpsi_dFj(3,3), dfpsi_duj(3,3), dfpsi_dthetaj(3,3), dfpsi_dMj(3,3);
    Eigen::VectorXd dfF_dFj(3), dfF_duj(3), dfF_dthetaj(3), dfF_dMj(3);
    Eigen::VectorXd dfM_dFj(3), dfM_duj(3), dfM_dthetaj(3), dfM_dMj(3);

    F_1 << F1(0),F1(1),F1(2);
    M_1 << F1(3),F1(4),F1(5);
    u1 << U1(0),U1(1),U1(2);
    theta1 << U1(3),U1(4),U1(5);
    Eigen::MatrixXd Sinv1(6,6);
    Sinv1 = S1.inverse();
    Strain1 = Sinv1*F1;
    gamma1 << Strain1(0),Strain1(1),Strain1(2);
    kappa1 << Strain1(3),Strain1(4),Strain1(5);
    e_1 << 1,0,0;

    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
        {
            dgammadF(i,j) = Sinv1(i,j);
            dgammadM(i,j) = Sinv1(i,j+3);
            dkappadF(i,j) = Sinv1(i+3,j);
            dkappadM(i,j) = Sinv1(i+3,j+3);
        }

    //Node i
    //f_u
    dfu_dFi = (Rot_Mat(theta1).transpose())*Transform_Mat(Theta);
    dfu_dui = Eigen::MatrixXd::Zero(3,3);
    dfu_dMi = Eigen::MatrixXd::Zero(3,3);

    //f_psi
    dfpsi_dFi = -(h/2)*(Rot_Mat(theta1).transpose()*Transform_Mat(Theta)*(Dyadic(e_1) + Dyadic(gamma1) - Dyadic(F_1)*dgammadF));
    dfpsi_dui = Eigen::MatrixXd::Zero(3,3);
    dfpsi_dMi = Rot_Mat(theta1)*Transform_Mat(Theta)*(Eigen::MatrixXd::Identity(3,3) - 0.5*h*Dyadic(F_1)*dgammadM);

    //f_F
    dfF_dFi = -(h/2)*(Rot_Mat(theta1)*Transform_Mat(Theta)*dgammadF);
    dfF_dui = Eigen::MatrixXd::Identity(3,3);
    dfF_dMi = -(h/2)*(Rot_Mat(theta1)*Transform_Mat(Theta)*dgammadM);

    //f_M
    dfM_dFi = -(h/2)*(Eigen::MatrixXd::Identity(3,3) + Dyadic(theta1)/2 + theta1*theta1.transpose()/4)*Transform_Mat(Theta)*dkappadF;
    dfM_dui = Eigen::MatrixXd::Zero(3,3);
    dfM_dMi = -(h/2)*(Eigen::MatrixXd::Identity(3,3) + Dyadic(theta1)/2 + theta1*theta1.transpose()/4)*Transform_Mat(Theta)*dkappadM;

    //dfdtheta
    Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(3,1);
    for(int i=0;i<2;i++)
    {
        //dfu_dtheta
        temp = dCdtheta(theta1,i+1).transpose()*Transform_Mat(Theta)*F_1;
        for(int j=0;j<3;j++)
            dfu_dthetai(i,j) = temp(j);
        //dfpsi_dtheta
        temp = dCdtheta(theta1,i+1).transpose()*Transform_Mat(Theta)*(M_1 - (h/2)*(Dyadic(e_1) + Dyadic(gamma1))*F_1);
        for(int j=0;j<3;j++)
            dfpsi_dthetai(i,j) = temp(j);
        //dff_dtheta
        temp = -(h/2)*(dCdtheta(theta1,i+1).transpose()*Transform_Mat(Theta)*(e_1 + gamma1));
        for(int j=0;j<3;j++)
            dfF_dthetai(i,j) = temp(j);
        //dfM_dtheta
        Eigen::VectorXd e(3);
        e << 0,0,0;
        if(i==1)
            e(0) = 1;
        else if(i==2)
            e(1) = 1;
        else if(i==3)
            e(2) = 1;
        temp = -e-(h/2)*(0.5*Dyadic(e) + 0.25*(e*theta1.transpose() + e.transpose()*theta1))*Transform_Mat(Theta)*kappa1;
        for(int j=0;j<3;j++)
            dfM_dthetai(i,j) = temp(j);
    }

    //Node i+1
    F_2 << F2(0),F2(1),F2(2);
    M_2 << F2(3),F2(4),F2(5);
    u2 << U2(0),U2(1),U2(2);
    theta2 << U2(3),U2(4),U2(5);
    Eigen::MatrixXd Sinv2(6,6);
    Sinv2 = S2.inverse();
    Strain2 = Sinv2*F2;
    gamma2 << Strain2(0),Strain2(1),Strain2(2);
    kappa2 << Strain2(3),Strain2(4),Strain2(5);

    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
        {
            dgammadF(i,j) = Sinv2(i,j);
            dgammadM(i,j) = Sinv2(i,j+3);
            dkappadF(i,j) = Sinv2(i+3,j);
            dkappadM(i,j) = Sinv2(i+3,j+3);
        }

    //f_u
    dfu_dFj = -(Rot_Mat(theta2).transpose())*Transform_Mat(Theta);
    dfu_duj = Eigen::MatrixXd::Zero(3,3);
    dfu_dMj = Eigen::MatrixXd::Zero(3,3);

    //f_psi
    dfpsi_dFj = -(h/2)*(Rot_Mat(theta2).transpose()*Transform_Mat(Theta)*(Dyadic(e_1) + Dyadic(gamma2) - Dyadic(F_2)*dgammadF));
    dfpsi_duj = Eigen::MatrixXd::Zero(3,3);
    dfpsi_dMj = -Rot_Mat(theta2)*Transform_Mat(Theta)*(Eigen::MatrixXd::Identity(3,3) - 0.5*h*Dyadic(F_2)*dgammadM);

    //f_F
    dfF_dFj = -(h/2)*(Rot_Mat(theta2)*Transform_Mat(Theta)*dgammadF);
    dfF_duj = Eigen::MatrixXd::Identity(3,3);
    dfF_dMj = -(h/2)*(Rot_Mat(theta2)*Transform_Mat(Theta)*dgammadM);

    //f_M
    dfM_dFj = -(h/2)*(Eigen::MatrixXd::Identity(3,3) + Dyadic(theta2)/2 + (theta2*theta2.transpose())/4)*Transform_Mat(Theta)*dkappadF;
    dfM_duj = Eigen::MatrixXd::Zero(3,3);
    dfM_dMj = -(h/2)*(Eigen::MatrixXd::Identity(3,3) + Dyadic(theta2)/2 + (theta2*theta2.transpose())/4)*Transform_Mat(Theta)*dkappadM;

    //dfdtheta
    temp = Eigen::MatrixXd::Zero(3,1);
    for(int i=0;i<2;i++)
    {
        //dfu_dtheta
        temp = -(dCdtheta(theta2,i+1).transpose())*Transform_Mat(Theta)*F_2;
        for(int j=0;j<3;j++)
            dfu_dthetaj(i,j) = temp(j);
        //dfpsi_dtheta
        temp = -dCdtheta(theta2,i+1).transpose()*Transform_Mat(Theta)*(M_2 - (h/2)*(Dyadic(e_1) + Dyadic(gamma2))*F_2);
        for(int j=0;j<3;j++)
            dfpsi_dthetaj(i,j) = temp(j);
        //dff_dtheta
        temp = -(h/2)*(dCdtheta(theta2,i+1).transpose()*Transform_Mat(Theta)*(e_1 + gamma2));
        for(int j=0;j<3;j++)
            dfF_dthetaj(i,j) = temp(j);
        //dfM_dtheta
        Eigen::VectorXd e(3);
        e << 0,0,0;
        if(i==1)
            e(0) = 1;
        else if(i==2)
            e(1) = 1;
        else if(i==3)
            e(2) = 1;
        temp = e-(h/2)*(0.5*Dyadic(e) + 0.25*(e*theta2.transpose() + theta2*e.transpose()))*Transform_Mat(Theta)*kappa2;
        for(int j=0;j<3;j++)
            dfM_dthetaj(i,j) = temp(j);
    }

    //Assembly
    Eigen::VectorXd J(12,24);
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
        {
            // Node i
            J(i,j) = dfu_dui(i,j);
            J(i+3,j) = dfu_dthetai(i,j);
            J(i+6,j) = dfu_dFi(i,j);
            J(i+9,j) = dfu_dMi(i,j);

            J(i,j+3) = dfpsi_dui(i,j);
            J(i+3,j+3) = dfpsi_dthetai(i,j);
            J(i+6,j+3) = dfpsi_dFi(i,j);
            J(i+9,j+3) = dfpsi_dMi(i,j);

            J(i,j+6) = dfF_dui(i,j);
            J(i+3,j+6) = dfF_dthetai(i,j);
            J(i+6,j+6) = dfF_dFi(i,j);
            J(i+9,j+6) = dfF_dMi(i,j);

            J(i,j+9) = dfM_dui(i,j);
            J(i+3,j+9) = dfM_dthetai(i,j);
            J(i+6,j+9) = dfM_dFi(i,j);
            J(i+9,j+9) = dfM_dMi(i,j);

            //Node i+1
            J(i,j+12) = dfu_duj(i,j);
            J(i+3,j+12) = dfu_dthetaj(i,j);
            J(i+6,j+12) = dfu_dFj(i,j);
            J(i+9,j+12) = dfu_dMj(i,j);

            J(i,j+15) = dfpsi_duj(i,j);
            J(i+3,j+15) = dfpsi_dthetaj(i,j);
            J(i+6,j+15) = dfpsi_dFj(i,j);
            J(i+9,j+15) = dfpsi_dMj(i,j);

            J(i,j+18) = dfF_duj(i,j);
            J(i+3,j+18) = dfF_dthetaj(i,j);
            J(i+6,j+18) = dfF_dFj(i,j);
            J(i+9,j+18) = dfF_dMj(i,j);

            J(i,j+21) = dfM_duj(i,j);
            J(i+3,j+21) = dfM_dthetaj(i,j);
            J(i+6,j+21) = dfM_dFj(i,j);
            J(i+9,j+21) = dfM_dMj(i,j);
        }
    return J;
}

Eigen::MatrixXd Element_Jacobian(Eigen::VectorXd U1, Eigen::VectorXd F1, double h, Eigen::MatrixXd S1, double Theta, int nodenum)
{
    Eigen::VectorXd F_1(3), M_1(3), u1(3), theta1(3), Strain1(6), gamma1(3), kappa1(3), e_1(3);
    Eigen::MatrixXd dgammadF(3,3), dgammadM(3,3), dkappadF(3,3), dkappadM(3,3);
    Eigen::MatrixXd dfu_dFi(3,3), dfu_dui(3,3), dfu_dthetai(3,3), dfu_dMi(3,3);
    Eigen::MatrixXd dfpsi_dFi(3,3), dfpsi_dui(3,3), dfpsi_dthetai(3,3), dfpsi_dMi(3,3);
    Eigen::VectorXd dfF_dFi(3), dfF_dui(3), dfF_dthetai(3), dfF_dMi(3);
    Eigen::VectorXd dfM_dFi(3), dfM_dui(3), dfM_dthetai(3), dfM_dMi(3);

    if(nodenum==0)
    {

        F_1 << F1(0),F1(1),F1(2);
        M_1 << F1(3),F1(4),F1(5);
        u1 << U1(0),U1(1),U1(2);
        theta1 << U1(3),U1(4),U1(5);
        Eigen::MatrixXd Sinv1(6,6);
        Sinv1 = S1.inverse();
        Strain1 = Sinv1*F1;
        gamma1 << Strain1(0),Strain1(1),Strain1(2);
        kappa1 << Strain1(3),Strain1(4),Strain1(5);
        e_1 << 1,0,0;

        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
            {
                dgammadF(i,j) = Sinv1(i,j);
                dgammadM(i,j) = Sinv1(i,j+3);
                dkappadF(i,j) = Sinv1(i+3,j);
                dkappadM(i,j) = Sinv1(i+3,j+3);
            }

        //f_u
        dfu_dFi = -(Rot_Mat(theta1).transpose())*Transform_Mat(Theta);
        dfu_dui = Eigen::MatrixXd::Zero(3,3);
        dfu_dMi = Eigen::MatrixXd::Zero(3,3);

        //f_psi
        dfpsi_dFi = -(h/2)*(Rot_Mat(theta1).transpose()*Transform_Mat(Theta)*(Dyadic(e_1) + Dyadic(gamma1) - Dyadic(F_1)*dgammadF));
        dfpsi_dui = Eigen::MatrixXd::Zero(3,3);
        dfpsi_dMi = -Rot_Mat(theta1)*Transform_Mat(Theta)*(Eigen::MatrixXd::Identity(3,3) - 0.5*h*Dyadic(F_1)*dgammadM);

        //f_F
        dfF_dFi = -(h/2)*(Rot_Mat(theta1)*Transform_Mat(Theta)*dgammadF);
        dfF_dui = Eigen::MatrixXd::Identity(3,3);
        dfF_dMi = -(h/2)*(Rot_Mat(theta1)*Transform_Mat(Theta)*dgammadM);

        //f_M
        dfM_dFi = -(h/2)*(Eigen::MatrixXd::Identity(3,3) + Dyadic(theta1)/2 + (theta1*theta1.transpose())/4)*Transform_Mat(Theta)*dkappadF;
        dfM_dui = Eigen::MatrixXd::Zero(3,3);
        dfM_dMi = -(h/2)*(Eigen::MatrixXd::Identity(3,3) + Dyadic(theta1)/2 + (theta1*theta1.transpose())/4)*Transform_Mat(Theta)*dkappadM;

        //dfdtheta
        Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(3,1);
        for(int i=0;i<2;i++)
        {
            //dfu_dtheta
            temp = -(dCdtheta(theta1,i+1).transpose())*Transform_Mat(Theta)*F_1;
            for(int j=0;j<3;j++)
                dfu_dthetai(i,j) = temp(j);
            //dfpsi_dtheta
            temp = -dCdtheta(theta1,i+1).transpose()*Transform_Mat(Theta)*(M_1 - (h/2)*(Dyadic(e_1) + Dyadic(gamma1))*F_1);
            for(int j=0;j<3;j++)
                dfpsi_dthetai(i,j) = temp(j);
            //dff_dtheta
            temp = -(h/2)*(dCdtheta(theta1,i+1).transpose()*Transform_Mat(Theta)*(e_1 + gamma1));
            for(int j=0;j<3;j++)
                dfF_dthetai(i,j) = temp(j);
            //dfM_dtheta
            Eigen::VectorXd e(3);
            e << 0,0,0;
            if(i==1)
                e(0) = 1;
            else if(i==2)
                e(1) = 1;
            else if(i==3)
                e(2) = 1;
            temp = e-(h/2)*(0.5*Dyadic(e) + 0.25*(e*theta1.transpose() + theta1*e.transpose()))*Transform_Mat(Theta)*kappa1;
            for(int j=0;j<3;j++)
                dfM_dthetai(i,j) = temp(j);
        }

        //Assembly
        Eigen::VectorXd J(12,12);
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
            {
                //Node i
                J(i,j) = dfu_dui(i,j);
                J(i+3,j) = dfu_dthetai(i,j);
                J(i+6,j) = dfu_dFi(i,j);
                J(i+9,j) = dfu_dMi(i,j);

                J(i,j+3) = dfpsi_dui(i,j);
                J(i+3,j+3) = dfpsi_dthetai(i,j);
                J(i+6,j+3) = dfpsi_dFi(i,j);
                J(i+9,j+3) = dfpsi_dMi(i,j);

                J(i,j+6) = dfF_dui(i,j);
                J(i+3,j+6) = dfF_dthetai(i,j);
                J(i+6,j+6) = dfF_dFi(i,j);
                J(i+9,j+6) = dfF_dMi(i,j);

                J(i,j+9) = dfM_dui(i,j);
                J(i+3,j+9) = dfM_dthetai(i,j);
                J(i+6,j+9) = dfM_dFi(i,j);
                J(i+9,j+9) = dfM_dMi(i,j);
            }
        return J;
    }
    else
    {

        F_1 << F1(0),F1(1),F1(2);
        M_1 << F1(3),F1(4),F1(5);
        u1 << U1(0),U1(1),U1(2);
        theta1 << U1(3),U1(4),U1(5);
        Eigen::MatrixXd Sinv1(6,6);
        Sinv1 = S1.inverse();
        Strain1 = Sinv1*F1;
        gamma1 << Strain1(0),Strain1(1),Strain1(2);
        kappa1 << Strain1(3),Strain1(4),Strain1(5);
        e_1 << 1,0,0;

        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
            {
                dgammadF(i,j) = Sinv1(i,j);
                dgammadM(i,j) = Sinv1(i,j+3);
                dkappadF(i,j) = Sinv1(i+3,j);
                dkappadM(i,j) = Sinv1(i+3,j+3);
            }

        //Node i
        //f_u
        dfu_dFi = (Rot_Mat(theta1).transpose())*Transform_Mat(Theta);
        dfu_dui = Eigen::MatrixXd::Zero(3,3);
        dfu_dMi = Eigen::MatrixXd::Zero(3,3);

        //f_psi
        dfpsi_dFi = -(h/2)*(Rot_Mat(theta1).transpose()*Transform_Mat(Theta)*(Dyadic(e_1) + Dyadic(gamma1) - Dyadic(F_1)*dgammadF));
        dfpsi_dui = Eigen::MatrixXd::Zero(3,3);
        dfpsi_dMi = Rot_Mat(theta1)*Transform_Mat(Theta)*(Eigen::MatrixXd::Identity(3,3) - 0.5*h*Dyadic(F_1)*dgammadM);

        //f_F
        dfF_dFi = -(h/2)*(Rot_Mat(theta1)*Transform_Mat(Theta)*dgammadF);
        dfF_dui = Eigen::MatrixXd::Identity(3,3);
        dfF_dMi = -(h/2)*(Rot_Mat(theta1)*Transform_Mat(Theta)*dgammadM);

        //f_M
        dfM_dFi = -(h/2)*(Eigen::MatrixXd::Identity(3,3) + Dyadic(theta1)/2 + theta1*theta1.transpose()/4)*Transform_Mat(Theta)*dkappadF;
        dfM_dui = Eigen::MatrixXd::Zero(3,3);
        dfM_dMi = -(h/2)*(Eigen::MatrixXd::Identity(3,3) + Dyadic(theta1)/2 + theta1*theta1.transpose()/4)*Transform_Mat(Theta)*dkappadM;
        //dfdtheta
        Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(3,1);
        for(int i=0;i<2;i++)
        {
            //dfu_dtheta
            temp = dCdtheta(theta1,i+1).transpose()*Transform_Mat(Theta)*F_1;
            for(int j=0;j<3;j++)
                dfu_dthetai(i,j) = temp(j);
            //dfpsi_dtheta
            temp = dCdtheta(theta1,i+1).transpose()*Transform_Mat(Theta)*(M_1 - (h/2)*(Dyadic(e_1) + Dyadic(gamma1))*F_1);
            for(int j=0;j<3;j++)
                dfpsi_dthetai(i,j) = temp(j);
            //dff_dtheta
            temp = -(h/2)*(dCdtheta(theta1,i+1).transpose()*Transform_Mat(Theta)*(e_1 + gamma1));
            for(int j=0;j<3;j++)
                dfF_dthetai(i,j) = temp(j);
            //dfM_dtheta
            Eigen::VectorXd e(3);
            e << 0,0,0;
            if(i==1)
                e(0) = 1;
            else if(i==2)
                e(1) = 1;
            else if(i==3)
                e(2) = 1;
            temp = -e-(h/2)*(0.5*Dyadic(e) + 0.25*(e*theta1.transpose() + e.transpose()*theta1))*Transform_Mat(Theta)*kappa1;
            for(int j=0;j<3;j++)
                dfM_dthetai(i,j) = temp(j);
        }
        //How to assemble?
        Eigen::VectorXd J(12,12);
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
            {
                // Node i
                J(i,j) = dfu_dui(i,j);
                J(i+3,j) = dfu_dthetai(i,j);
                J(i+6,j) = dfu_dFi(i,j);
                J(i+9,j) = dfu_dMi(i,j);

                J(i,j+3) = dfpsi_dui(i,j);
                J(i+3,j+3) = dfpsi_dthetai(i,j);
                J(i+6,j+3) = dfpsi_dFi(i,j);
                J(i+9,j+3) = dfpsi_dMi(i,j);

                J(i,j+6) = dfF_dui(i,j);
                J(i+3,j+6) = dfF_dthetai(i,j);
                J(i+6,j+6) = dfF_dFi(i,j);
                J(i+9,j+6) = dfF_dMi(i,j);

                J(i,j+9) = dfM_dui(i,j);
                J(i+3,j+9) = dfM_dthetai(i,j);
                J(i+6,j+9) = dfM_dFi(i,j);
                J(i+9,j+9) = dfM_dMi(i,j);
            }
        return J;
    }
}


//This function gives the equivalent stiffness matrix in 4*4 form after d
//differentiating twice the one dimensional strain energy with respect to strains
Eigen::MatrixXd ClassicalStiffnessModel_VAM(Eigen::MatrixXd S, Eigen::VectorXd Strain)
{
    Eigen::MatrixXd Seq = Eigen::MatrixXd::Zero(4,4);
    Eigen::MatrixXd Sl(4,4), Sln(4,5), Sn(5,5);
    for(int i=0;i<4;i++)
        for(int j=0;j<4;j++)
            Sl(i,j) = S(i,j);
    for(int i=0;i<4;i++)
        for(int j=0;j<5;j++)
        {
            Sln(i,j) = S(i+4,j);
        }
    for(int i=0;i<5;i++)
        for(int j=0;j<5;j++)
            Sn(i,j) = S(i+4,j+4);
    //Strains
    double gamma11, kappa1, kappa2, kappa3;
    gamma11 = Strain(0);
    kappa1 = Strain(1);
    kappa2 = Strain(2);
    kappa3 = Strain(3);

    //d^2U/dgamma11^2
    Seq(0,0) = 0.5*Sl(0,0) + kappa2*(Sln(0,2) + 0.5*Sn(2,2)*kappa2);
    //d^2U/dgamma11*dkappa1
    Seq(0,1) = 0.5*Sl(0,1) + Sln(0,0)*kappa1 + kappa2*(0.5*Sln(0,4) + 0.5*Sln(1,2) + Sn(0,2)*kappa1 +
              0.5*Sn(2,4)*kappa2);
    //d^2U/dgamma11*dkappa2
    Seq(0,2) = 0.5*Sl(0,2) + Sln(0,2)*gamma11 + 0.5*Sln(0,4)*kappa1 + 0.5*Sln(1,2)*kappa1 + 0.5*Sn(0,2)*pow(kappa1,2)
            + Sln(0,1)*kappa2 + Sln(2,2)*kappa2 + Sn(2,2)*gamma11*kappa2 +
            Sn(2,4)*kappa1*kappa2 + 1.5*Sn(1,2)*pow(kappa2,2) +
            0.5*Sln(0,3)*kappa3 + 0.5*Sln(2,3)*kappa3 +
            Sn(2,3)*kappa2*kappa3;
    //d^2U/dgamma11*dkappa3
    Seq(0,3) = 0.5*Sl(0,3) + 0.5*kappa2*Sln(0,3) + Sln(2,3) +
               Sn(2,3)*kappa2;
    //d^U/dkappa1*dgamma11
    Seq(1,0) = Seq(0,1);
    //d^U/dkappa1^2
    Seq(1,1) = 0.5*Sl(1,1) + Sln(0,0)*gamma11 + 3*Sln(0,1)*kappa1 + 3*Sn(0,0)*pow(kappa1,2) +
            Sln(0,2)*kappa2 + Sln(1,4)*kappa2 + Sn(0,2)*gamma11*kappa2 +
            3*Sn(0,4)*kappa1*kappa2 + Sn(0,1)*pow(kappa2,2) +
            0.5*Sn(4,4)*pow(kappa2,2) +
            Sln(0,3)*kappa3 + Sn(0,3)*kappa2*kappa3;
    //d^U/dkappa1*dkappa2
    Seq(1,2) = 0.5*Sl(1,2) + 0.5*Sln(0,4)*gamma11 + 0.5*Sln(1,2)*gamma11 +
            Sln(0,2)*kappa1 + Sln(1,4)*kappa1 + Sn(0,2)*gamma11*kappa1 + 1.5*Sn(0,4)*pow(kappa1,2) +
            Sln(1,1)*kappa2 + Sln(2,4)*kappa2 + Sn(2,4)*gamma11*kappa2 +
            2*Sn(0,1)*kappa1*kappa2 + Sn(4,4)*kappa1*kappa2 +
            1.5*Sn(1,4)*pow(kappa2,2) + 0.5*Sln(1,3)*kappa3 + 0.5*Sln(3,4)*kappa3 +
            Sn(0,3)*kappa1*kappa3 + Sn(3,4)*kappa2*kappa3;
    //d^U/dkappa1*dkappa3
    Seq(1,3) = 0.5*Sl(1,3) + Sln(0,3)*kappa1 +
               kappa2*(0.5*Sln(1,3) + 0.5*Sln(3,4) + Sn(0,3)*kappa1 +
                0.5*Sn(3,4)*kappa2);
    //d^U/dkappa2*dgamma11
    Seq(2,0) = Seq(0,2);
    //d^U/dkappa2*dkappa1
    Seq(2,1) = Seq(1,2);
    //d^U/dkappa2*dkappa2
    Seq(2,2) = 0.5*Sl(2,2) + Sln(0,1)*gamma11 + Sln(2,2)*gamma11 +
            0.5*Sn(2,2)*pow(gamma11,2) + Sln(1,1)*kappa1 +
            Sln(2,4)*kappa1 + Sn(2,4)*gamma11*kappa1 +
            Sn(0,1)*pow(kappa1,2) + 0.5*Sn(4,4)*pow(kappa1,2) +
            3*Sln(1,2)*kappa2 + 3*Sn(1,2)*gamma11*kappa2 +
            3*Sn(1,4)*kappa1*kappa2 + 3*Sn(1,1)*pow(kappa2,2) +
            Sln(1,3)*kappa3 + Sln(2,3)*kappa3 + Sn(2,3)*gamma11*kappa3 +
            Sn(3,4)*kappa1*kappa3 + 3*Sn(1,3)*kappa2*kappa3 +
            0.5*Sn(3,3)*pow(kappa3,2);
    //d^U/dkappa2*dkappa3
    Seq(2,3) = 0.5*Sl(2,3) + 0.5*Sln(0,3)*gamma11 + 0.5*Sln(2,3)*gamma11 +
            0.5*Sln(1,3)*kappa1 + 0.5*Sln(3,4)*kappa1 + 0.5*Sn(0,3)*pow(kappa1,2) +
            Sln(1,3)*kappa2 + Sln(2,3)*kappa2 + Sn(2,3)*gamma11*kappa2 +
            Sn(3,4)*kappa1*kappa2 + 1.5*Sn(1,3)*pow(kappa2,2) + Sln(3,3)*kappa3 +
            Sn(3,3)*kappa2*kappa3;
    //d^U/dkappa3*dgamma11
    Seq(3,0) = Seq(0,3);
    //d^U/dkappa3*dkappa1
    Seq(3,1) = Seq(1,3);
    //d^U/dkappa3*dkappa2
    Seq(3,2) = Seq(2,3);
    //d^U/dkappa3*dkappa3
    Seq(3,3) = 0.5*Sl(3,3) + kappa2*Sln(3,3) + 0.5*Sn(3,3)*kappa2;

    //The matrix is returned in the form of a 6*6 matrix to generalize for Timoshenko
    //This also helps with the generalization of the code as we can assume each node of the beam element has 3 dofs
    //For the present case u1, theta1, theta2, and theta3 are the dofs
    //u2 and u3 are taken as zero
    Eigen::MatrixXd Stiffness(6,6);

    for(int i=0;i<6;i++)
        for(int j=0;j<6;j++)
        {
            if(i==1||i==2)
                Stiffness(i,j) = 0;
            else
                Stiffness(i,j) = Seq(i,j);
        }

    return Stiffness;
}

//After the 1D non linear beam analysis is done, the values of unknowns which are
//Forces and displacements are obtained
//New strains are found from the constitutive relation. Since the stiffness matrix is a function
//of Strains an iterative loop is used to find the converged value of strains.
void Update_Strains(Eigen::VectorXd &Tau, int nnode, Eigen::MatrixXd S, Eigen::VectorXd U)
{
    for(int i=0;i<nnode;i++)
    {
        Eigen::VectorXd Force(6);
        for(int j=0;j<6;j++)
            Force(j) = U(12*j);
        Eigen::VectorXd Strain(6), Strain_new(6);
        for(int j=0;j<6;j++)
            Strain(j) = 0;
        Eigen::MatrixXd Seq(6,6);
        Seq = ClassicalStiffnessModel_VAM(S,Strain);
        double maxerror, iter = 0;
        do
        {
             Strain_new = Seq.colPivHouseholderQr().solve(Force);
             Seq = ClassicalStiffnessModel_VAM(S,Strain_new);
             Eigen::VectorXd error(6);
             for(int j=0;j<6;j++)
                error(j) = abs(Strain_new(j) - Strain(j))/Strain(j);
             maxerror = error(0);
             for(int j=1;j<6;j++)
                 if(error(j)>maxerror)
                     maxerror = error(j);
             for(int j=0;j<6;j++)
                 Strain(j) = Strain_new(j);
             iter++;
        }while(maxerror>pow(10,-6)&&iter<20);
        for(int j=0;j<6;j++)
            Tau(12*i+j) = Strain(j);
    }
}

#endif

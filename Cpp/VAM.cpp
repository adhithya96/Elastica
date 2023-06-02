
#ifndef VARIABLES_H
#define VARIABLES_H


#include "Variables.h"
//2D Cross-Section Analysis
//S is the cross sectional stiffness matrix obtained from the cross sectional analysis
//9*9 stiffness matrix formula is borrowed from "Non-classical effects in non-linear analysis of pretwisted
//anisotropic strips" by Harursampath et. al
//It is converted to 4*4 by taking double derivative of one dimensional strain energy with respect to the strains.
//The matrix is in general converted to 6*6 to generalize the code such that there's 3 dofs at each node accounting
//for a total of 6 dofs for each element.
//This can be extended to Timoshenko formulation as well.
//Formula for 1D strain energy can be found in Eq. (31) in the paper.
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

    //Cross-Sectional Stiffness for 0/0 ply

    S(0, 0) = 502859;
    S(0, 1) = 0;
    S(0, 2) = 0;
    S(0, 3) = 27.0354*kappa1;
    S(0, 4) = 0;
    S(0, 5) = 0;
    
    S(1, 0) = 0;
    S(1, 1) = 1e20;
    S(1, 2) = 0;
    S(1, 3) = 0;
    S(1, 4) = 0;
    S(1, 5) = 0;

    S(2, 0) = 0;
    S(2, 1) = 0;
    S(2, 2) = 1e20;
    S(2, 3) = 0;
    S(2, 4) = 0;
    S(2, 5) = 0;

    S(3, 0) = 27.0354 * kappa1;
    S(3, 1) = 0;
    S(3, 2) = 0;
    S(3, 3) = 0.000110668 + 27.0354 * kappa1 + 0.0039245 * pow(kappa1, 2) + 0.000174421 * pow(kappa2, 2);
    S(3, 4) = 0.000348843 * kappa1 * kappa2;
    S(3, 5) = 0;

    S(4, 0) = 0;
    S(4, 1) = 0;
    S(4, 2) = 0;
    S(4, 3) = 0.000348843 * kappa1 * kappa2;
    S(4, 4) = 0.000893245 + 0.000174421 * pow(kappa1, 2) + 0.0001567979 * pow(kappa2, 2) - 34.2061 * pow(kappa3, 2);
    S(4, 5) = -68.4122 * kappa2 * kappa3;


    S(5, 0) = 0;
    S(5, 1) = 0;
    S(5, 2) = 0;
    S(5, 3) = 0;
    S(5, 4) = -68.4122 * kappa2 * kappa3;
    S(5, 5) = 27.0354 - 34.2061 * pow(kappa2, 2);

    return S;
}

#endif
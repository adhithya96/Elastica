#pragma once

#include "Variables.h"

//Inp file
//Single beam case
//Isotropic validation with elastica
VAMBeamElement::VAMBeamElement()
{
    this->NNODE = 17;
    this->NELEM = 16;
    this->NDOF = 12;
    this->NLS = 10;

    this->NODE = Eigen::MatrixXd::Zero(this->NNODE, 3);
    this->ELEM = Eigen::MatrixXd::Zero(this->NELEM, 3);

    //    Nodal Information
    this->NODE(0, 0) = 0;
    this->NODE(0, 1) = 0;
    this->NODE(0, 2) = 0;
    this->NODE(1, 0) = 6.25;
    this->NODE(1, 1) = 0;
    this->NODE(1, 2) = 0;
    this->NODE(2, 0) = 12.5;
    this->NODE(2, 1) = 0;
    this->NODE(2, 2) = 0;
    this->NODE(3, 0) = 18.75;
    this->NODE(3, 1) = 0;
    this->NODE(3, 2) = 0;
    this->NODE(4, 0) = 25.0;
    this->NODE(4, 1) = 0;
    this->NODE(4, 2) = 0;
    this->NODE(5, 0) = 31.25;
    this->NODE(5, 1) = 0;
    this->NODE(5, 2) = 0;
    this->NODE(6, 0) = 37.5;
    this->NODE(6, 1) = 0;
    this->NODE(6, 2) = 0;
    this->NODE(7, 0) = 43.75;
    this->NODE(7, 1) = 0;
    this->NODE(7, 2) = 0;
    this->NODE(8, 0) = 50.0;
    this->NODE(8, 1) = 0;
    this->NODE(8, 2) = 0;
    this->NODE(9, 0) = 56.25;
    this->NODE(9, 1) = 0;
    this->NODE(9, 2) = 0;
    this->NODE(10, 0) = 62.5;
    this->NODE(10, 1) = 0;
    this->NODE(10, 2) = 0;
    this->NODE(11, 0) = 68.75;
    this->NODE(11, 1) = 0;
    this->NODE(11, 2) = 0;
    this->NODE(12, 0) = 75.0;
    this->NODE(12, 1) = 0;
    this->NODE(12, 2) = 0;
    this->NODE(13, 0) = 81.25;
    this->NODE(13, 1) = 0;
    this->NODE(13, 2) = 0;
    this->NODE(14, 0) = 87.5;
    this->NODE(14, 1) = 0;
    this->NODE(14, 2) = 0;
    this->NODE(15, 0) = 93.75;
    this->NODE(15, 1) = 0;
    this->NODE(15, 2) = 0;
    this->NODE(16, 0) = 100;
    this->NODE(16, 1) = 0;
    this->NODE(16, 2) = 0;

    //Connectivity information
    this->ELEM(0, 0) = 1;
    this->ELEM(0, 1) = 1;
    this->ELEM(0, 2) = 2;
    this->ELEM(1, 0) = 1;
    this->ELEM(1, 1) = 2;
    this->ELEM(1, 2) = 3;
    this->ELEM(2, 0) = 1;
    this->ELEM(2, 1) = 3;
    this->ELEM(2, 2) = 4;
    this->ELEM(3, 0) = 1;
    this->ELEM(3, 1) = 4;
    this->ELEM(3, 2) = 5;
    this->ELEM(4, 0) = 1;
    this->ELEM(4, 1) = 5;
    this->ELEM(4, 2) = 6;
    this->ELEM(5, 0) = 1;
    this->ELEM(5, 1) = 6;
    this->ELEM(5, 2) = 7;
    this->ELEM(6, 0) = 1;
    this->ELEM(6, 1) = 7;
    this->ELEM(6, 2) = 8;
    this->ELEM(7, 0) = 1;
    this->ELEM(7, 1) = 8;
    this->ELEM(7, 2) = 9;
    this->ELEM(8, 0) = 1;
    this->ELEM(8, 1) = 9;
    this->ELEM(8, 2) = 10;
    this->ELEM(9, 0) = 1;
    this->ELEM(9, 1) = 10;
    this->ELEM(9, 2) = 11;
    this->ELEM(10, 0) = 1;
    this->ELEM(10, 1) = 11;
    this->ELEM(10, 2) = 12;
    this->ELEM(11, 0) = 1;
    this->ELEM(11, 1) = 12;
    this->ELEM(11, 2) = 13;
    this->ELEM(12, 0) = 1;
    this->ELEM(12, 1) = 13;
    this->ELEM(12, 2) = 14;
    this->ELEM(13, 0) = 1;
    this->ELEM(13, 1) = 14;
    this->ELEM(13, 2) = 15;
    this->ELEM(14, 0) = 1;
    this->ELEM(14, 1) = 15;
    this->ELEM(14, 2) = 16;
    this->ELEM(15, 0) = 1;
    this->ELEM(15, 1) = 16;
    this->ELEM(15, 2) = 17;

    this->inittwist = Eigen::VectorXd::Zero(3);
    this->inittwist(0) = 0;
    this->inittwist(1) = 0;
    this->inittwist(2) = 0;

    this->a1 = Eigen::VectorXd::Zero(3);
    this->b1 = Eigen::VectorXd::Zero(3);
    this->aN = Eigen::VectorXd::Zero(3);
    this->bN = Eigen::VectorXd::Zero(3);

    this->aN(2) = 400;

}

//2D Cross-Section Analysis
//S is the cross sectional stiffness matrix obtained from the cross sectional analysis
//9*9 stiffness matrix formula is borrowed from "Non-classical effects in non-linear analysis of pretwisted
//anisotropic strips" by Harursampath et. al
//It is converted to 4*4 by taking double derivative of one dimensional strain energy with respect to the strains.
//The matrix is in general converted to 6*6 to generalize the code such that there's 3 dofs at each node accounting
//for a total of 6 dofs for each element.
//This can be extended to Timoshenko formulation as well.
//Formula for 1D strain energy can be found in Eq. (31) in the paper.
Eigen::MatrixXd VAMBeamElement::Equivalent_StiffnessMatrix_FirstOrder(std::fstream& file1)
{
    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(6, 6);

    //cross-sectional matrix for isotropic case
    S(0, 0) = 1.20024245e+07;
    S(0, 1) = 1.21226616e+05;
    S(0, 2) = 1.21226616e+05;
    S(0, 3) = -8.61796472e-13;
    S(0, 4) = -4.91340022e-14;
    S(0, 5) = -1.16220209e-10;

    S(1, 0) = 1.21226616e+05;
    S(1, 1) = 1.20017403e+07;
    S(1, 2) = 1.20921284e+05;
    S(1, 3) = -1.07935374e-10;
    S(1, 4) = 9.23952828e-12;
    S(1, 5) = -9.70156484e-12;

    S(2, 0) = 1.21226616e+05;
    S(2, 1) = 1.20921284e+05;
    S(2, 2) = 1.20017403e+07;
    S(2, 3) = -1.06141966e-10;
    S(2, 4) = -3.05978023e-11;
    S(2, 5) = 1.86570753e-11;
    
    S(3, 0) = 1.82405507e-12;
    S(3, 1) = -1.07935374e-10;
    S(3, 2) = -1.06141966e-10;
    S(3, 3) = 2.00025105e+06;
    S(3, 4) = -1.00265129e+04;
    S(3, 5) = -1.00265129e+04;

    S(4, 0) = -3.11269079e-14;
    S(4, 1) = 9.23952828e-12;
    S(4, 2) = -3.05978023e-11;
    S(4, 3) = -1.00265129e+04;
    S(4, 4) = 1.00007526e+06;
    S(4, 5) = 5.02591713e+01;

    S(5, 0) = -1.16408747e-10;
    S(5, 1) = -9.70156484e-12;
    S(5, 2) = 1.86570753e-11;
    S(5, 3) = -1.00265129e+04;
    S(5, 4) = 5.02591713e+01;
    S(5, 5) = 1.00007526e+06;

    //std::cout << S << std::endl;

    //Cross-Sectional Stiffness for 0/0 ply

    /*S(0, 0) = 502859;
    S(0, 1) = 0;
    S(0, 2) = 0;
    S(0, 3) = 27.0354 * kappa1;
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
    S(5, 5) = 27.0354 - 34.2061 * pow(kappa2, 2);*/

    return S;
}
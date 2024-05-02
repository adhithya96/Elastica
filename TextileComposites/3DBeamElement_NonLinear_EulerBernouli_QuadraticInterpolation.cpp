#ifndef VARIABLES_H
#define VARIABLES_H
#include "Variables.h"
#include<iomanip>



/*************************************************************
* AceGen    7.505 Windows (16 Aug 22)                        *
*           Co. J. Korelc  2020           18 May 23 15:11:18 *
**************************************************************
User     : Full professional version
Notebook : 3DBeamElement_EulerBernouli_Quadratic_User
Evaluation time                 : 24 s    Mode  : Debug
Number of formulae              : 1015    Method: Automatic
Subroutine                      : RKt size: 21520
Total size of Mathematica  code : 21520 subexpressions
Total size of C code            : 103503 bytes */
#include<iostream>


//INP file
//Validation case for Euler Bernouli Beam element linear elastic model(6 dofs per node)
/*NLEBBE3D::NLEBBE3D(std::string str)
{
    if (choice == 1)
    {
        NLEBBE3D::MAT = str;

        NLEBBE3D::NNODE = 17;
        NLEBBE3D::NELEM = 8;
        NLEBBE3D::NDOF = 6;
        NLEBBE3D::NLS = 10;
        NLEBBE3D::NEN = 3;
        NLEBBE3D::NDM = 3;

        NLEBBE3D::NODE = Eigen::MatrixXd::Zero(NLEBBE3D::NNODE, NLEBBE3D::NDM);
        NLEBBE3D::ELEM = Eigen::MatrixXd::Zero(NLEBBE3D::NELEM, NLEBBE3D::NEN + 1);

        //Nodal Information
        NLEBBE3D::NODE(0, 0) = 0;
        NLEBBE3D::NODE(0, 1) = 0;
        NLEBBE3D::NODE(0, 2) = 0;
        NLEBBE3D::NODE(1, 0) = 0.1204;
        NLEBBE3D::NODE(1, 1) = 4.9067;
        NLEBBE3D::NODE(1, 2) = 0;
        NLEBBE3D::NODE(2, 0) = 0.4815;
        NLEBBE3D::NODE(2, 1) = 9.8017;
        NLEBBE3D::NODE(2, 2) = 0;
        NLEBBE3D::NODE(3, 0) = 1.082;
        NLEBBE3D::NODE(3, 1) = 14.673;
        NLEBBE3D::NODE(3, 2) = 0;
        NLEBBE3D::NODE(4, 0) = 1.9214;
        NLEBBE3D::NODE(4, 1) = 19.509;
        NLEBBE3D::NODE(4, 2) = 0;
        NLEBBE3D::NODE(5, 0) = 2.9968;
        NLEBBE3D::NODE(5, 1) = 24.298;
        NLEBBE3D::NODE(5, 2) = 0;
        NLEBBE3D::NODE(6, 0) = 4.306;
        NLEBBE3D::NODE(6, 1) = 29.03;
        NLEBBE3D::NODE(6, 2) = 0;
        NLEBBE3D::NODE(7, 0) = 5.845;
        NLEBBE3D::NODE(7, 1) = 33.69;
        NLEBBE3D::NODE(7, 2) = 0;
        NLEBBE3D::NODE(8, 0) = 7.612;
        NLEBBE3D::NODE(8, 1) = 38.27;
        NLEBBE3D::NODE(8, 2) = 0;
        NLEBBE3D::NODE(9, 0) = 9.601;
        NLEBBE3D::NODE(9, 1) = 42.75;
        NLEBBE3D::NODE(9, 2) = 0;
        NLEBBE3D::NODE(10, 0) = 11.807;
        NLEBBE3D::NODE(10, 1) = 47.14;
        NLEBBE3D::NODE(10, 2) = 0;
        NLEBBE3D::NODE(11, 0) = 14.23;
        NLEBBE3D::NODE(11, 1) = 51.41;
        NLEBBE3D::NODE(11, 2) = 0;
        NLEBBE3D::NODE(12, 0) = 16.853;
        NLEBBE3D::NODE(12, 1) = 55.56;
        NLEBBE3D::NODE(12, 2) = 0;
        NLEBBE3D::NODE(13, 0) = 19.68;
        NLEBBE3D::NODE(13, 1) = 59.57;
        NLEBBE3D::NODE(13, 2) = 0;
        NLEBBE3D::NODE(14, 0) = 22.698;
        NLEBBE3D::NODE(14, 1) = 63.44;
        NLEBBE3D::NODE(14, 2) = 0;
        NLEBBE3D::NODE(15, 0) = 25.9;
        NLEBBE3D::NODE(15, 1) = 67.15;
        NLEBBE3D::NODE(15, 2) = 0;
        NLEBBE3D::NODE(16, 0) = 29.3;
        NLEBBE3D::NODE(16, 1) = 70.71;
        NLEBBE3D::NODE(16, 2) = 0;

        //Element connectivity
        NLEBBE3D::ELEM(0, 0) = 1;
        NLEBBE3D::ELEM(0, 1) = 1;
        NLEBBE3D::ELEM(0, 2) = 2;
        NLEBBE3D::ELEM(0, 3) = 3;

        NLEBBE3D::ELEM(1, 0) = 1;
        NLEBBE3D::ELEM(1, 1) = 3;
        NLEBBE3D::ELEM(1, 2) = 4;
        NLEBBE3D::ELEM(1, 3) = 5;


        NLEBBE3D::ELEM(2, 0) = 1;
        NLEBBE3D::ELEM(2, 1) = 5;
        NLEBBE3D::ELEM(2, 2) = 6;
        NLEBBE3D::ELEM(2, 3) = 7;

        NLEBBE3D::ELEM(3, 0) = 1;
        NLEBBE3D::ELEM(3, 1) = 7;
        NLEBBE3D::ELEM(3, 2) = 8;
        NLEBBE3D::ELEM(3, 3) = 9;

        NLEBBE3D::ELEM(4, 0) = 1;
        NLEBBE3D::ELEM(4, 1) = 9;
        NLEBBE3D::ELEM(4, 2) = 10;
        NLEBBE3D::ELEM(4, 3) = 11;

        NLEBBE3D::ELEM(5, 0) = 1;
        NLEBBE3D::ELEM(5, 1) = 11;
        NLEBBE3D::ELEM(5, 2) = 12;
        NLEBBE3D::ELEM(5, 3) = 13;

        NLEBBE3D::ELEM(6, 0) = 1;
        NLEBBE3D::ELEM(6, 1) = 13;
        NLEBBE3D::ELEM(6, 2) = 14;
        NLEBBE3D::ELEM(6, 3) = 15;

        NLEBBE3D::ELEM(7, 0) = 1;
        NLEBBE3D::ELEM(7, 1) = 15;
        NLEBBE3D::ELEM(7, 2) = 16;
        NLEBBE3D::ELEM(7, 3) = 17;

        NLEBBE3D::E = 1e7;
        NLEBBE3D::nu = 0.0001;
        NLEBBE3D::Bp = 1;
        NLEBBE3D::Hp = 1;
        NLEBBE3D::Zx = 0;
        NLEBBE3D::Zy = 0;
        NLEBBE3D::Zz = 1;

        NLEBBE3D::DIA = 2.0 / sqrt(M_PI);
    }
    else
        std::cout << "Wrong input" << std::endl;
}*/

NLEBBE3D::NLEBBE3D(std::string str)
{
    NLEBBE3D::MAT = str;

    NLEBBE3D::NNODE = 21;
    NLEBBE3D::NELEM = 10;
    NLEBBE3D::NDOF = 6;
    NLEBBE3D::NLS = 10;
    NLEBBE3D::NEN = 3;
    NLEBBE3D::NDM = 3;

    NLEBBE3D::NODE = Eigen::MatrixXd::Zero(NLEBBE3D::NNODE, NLEBBE3D::NDM);
    NLEBBE3D::ELEM = Eigen::MatrixXd::Zero(NLEBBE3D::NELEM, NLEBBE3D::NEN + 1);

    //Nodal Information
    NLEBBE3D::NODE(0, 0) = -250;
    NLEBBE3D::NODE(0, 1) = -83.33;
    NLEBBE3D::NODE(0, 2) = 0;
    NLEBBE3D::NODE(1, 0) = -225;
    NLEBBE3D::NODE(1, 1) = -83.33;
    NLEBBE3D::NODE(1, 2) = 0;
    NLEBBE3D::NODE(2, 0) = -200;
    NLEBBE3D::NODE(2, 1) = -83.33;
    NLEBBE3D::NODE(2, 2) = 0;
    NLEBBE3D::NODE(3, 0) = -175;
    NLEBBE3D::NODE(3, 1) = -83.33;
    NLEBBE3D::NODE(3, 2) = 0;
    NLEBBE3D::NODE(4, 0) = -150;
    NLEBBE3D::NODE(4, 1) = -83.33;
    NLEBBE3D::NODE(4, 2) = 0;
    NLEBBE3D::NODE(5, 0) = -125;
    NLEBBE3D::NODE(5, 1) = -83.33;
    NLEBBE3D::NODE(5, 2) = 0;
    NLEBBE3D::NODE(6, 0) = -100;
    NLEBBE3D::NODE(6, 1) = -83.33;
    NLEBBE3D::NODE(6, 2) = 0;
    NLEBBE3D::NODE(7, 0) = -75;
    NLEBBE3D::NODE(7, 1) = -83.33;
    NLEBBE3D::NODE(7, 2) = 0;
    NLEBBE3D::NODE(8, 0) = -50;
    NLEBBE3D::NODE(8, 1) = -83.33;
    NLEBBE3D::NODE(8, 2) = 0;
    NLEBBE3D::NODE(9, 0) = -25;
    NLEBBE3D::NODE(9, 1) = -83.33;
    NLEBBE3D::NODE(9, 2) = 0;
    NLEBBE3D::NODE(10, 0) = 0;
    NLEBBE3D::NODE(10, 1) = -83.33;
    NLEBBE3D::NODE(10, 2) = 0;
    NLEBBE3D::NODE(11, 0) = 25;
    NLEBBE3D::NODE(11, 1) = -83.33;
    NLEBBE3D::NODE(11, 2) = 0;
    NLEBBE3D::NODE(12, 0) = 50;
    NLEBBE3D::NODE(12, 1) = -83.33;
    NLEBBE3D::NODE(12, 2) = 0;
    NLEBBE3D::NODE(13, 0) = 75;
    NLEBBE3D::NODE(13, 1) = -83.33;
    NLEBBE3D::NODE(13, 2) = 0;
    NLEBBE3D::NODE(14, 0) = 100;
    NLEBBE3D::NODE(14, 1) = -83.33;
    NLEBBE3D::NODE(14, 2) = 0;
    NLEBBE3D::NODE(15, 0) = 125;
    NLEBBE3D::NODE(15, 1) = -83.33;
    NLEBBE3D::NODE(15, 2) = 0;
    NLEBBE3D::NODE(16, 0) = 150;
    NLEBBE3D::NODE(16, 1) = -83.33;
    NLEBBE3D::NODE(16, 2) = 0;
    NLEBBE3D::NODE(17, 0) = 175;
    NLEBBE3D::NODE(17, 1) = -83.33;
    NLEBBE3D::NODE(17, 2) = 0;
    NLEBBE3D::NODE(18, 0) = 200;
    NLEBBE3D::NODE(18, 1) = -83.33;
    NLEBBE3D::NODE(18, 2) = 0;
    NLEBBE3D::NODE(19, 0) = 225;
    NLEBBE3D::NODE(19, 1) = -83.33;
    NLEBBE3D::NODE(19, 2) = 0;
    NLEBBE3D::NODE(20, 0) = 250;
    NLEBBE3D::NODE(20, 1) = -83.33;
    NLEBBE3D::NODE(20, 2) = 0;

    //Element connectivity
    NLEBBE3D::ELEM(0, 0) = 1;
    NLEBBE3D::ELEM(0, 1) = 1;
    NLEBBE3D::ELEM(0, 2) = 2;
    NLEBBE3D::ELEM(0, 3) = 3;

    NLEBBE3D::ELEM(1, 0) = 1;
    NLEBBE3D::ELEM(1, 1) = 3;
    NLEBBE3D::ELEM(1, 2) = 4;
    NLEBBE3D::ELEM(1, 3) = 5;


    NLEBBE3D::ELEM(2, 0) = 1;
    NLEBBE3D::ELEM(2, 1) = 5;
    NLEBBE3D::ELEM(2, 2) = 6;
    NLEBBE3D::ELEM(2, 3) = 7;

    NLEBBE3D::ELEM(3, 0) = 1;
    NLEBBE3D::ELEM(3, 1) = 7;
    NLEBBE3D::ELEM(3, 2) = 8;
    NLEBBE3D::ELEM(3, 3) = 9;

    NLEBBE3D::ELEM(4, 0) = 1;
    NLEBBE3D::ELEM(4, 1) = 9;
    NLEBBE3D::ELEM(4, 2) = 10;
    NLEBBE3D::ELEM(4, 3) = 11;

    NLEBBE3D::ELEM(5, 0) = 1;
    NLEBBE3D::ELEM(5, 1) = 11;
    NLEBBE3D::ELEM(5, 2) = 12;
    NLEBBE3D::ELEM(5, 3) = 13;

    NLEBBE3D::ELEM(6, 0) = 1;
    NLEBBE3D::ELEM(6, 1) = 13;
    NLEBBE3D::ELEM(6, 2) = 14;
    NLEBBE3D::ELEM(6, 3) = 15;

    NLEBBE3D::ELEM(7, 0) = 1;
    NLEBBE3D::ELEM(7, 1) = 15;
    NLEBBE3D::ELEM(7, 2) = 16;
    NLEBBE3D::ELEM(7, 3) = 17;

    NLEBBE3D::ELEM(8, 0) = 1;
    NLEBBE3D::ELEM(8, 1) = 17;
    NLEBBE3D::ELEM(8, 2) = 18;
    NLEBBE3D::ELEM(8, 3) = 19;

    NLEBBE3D::ELEM(9, 0) = 1;
    NLEBBE3D::ELEM(9, 1) = 19;
    NLEBBE3D::ELEM(9, 2) = 20;
    NLEBBE3D::ELEM(9, 3) = 21;

    NLEBBE3D::E = 2.4e3;
    NLEBBE3D::nu = 0.39;
    NLEBBE3D::Bp = 1;
    NLEBBE3D::Hp = 1;
    NLEBBE3D::Zx = 0;
    NLEBBE3D::Zy = 0;
    NLEBBE3D::Zz = 1;

    NLEBBE3D::DIA = 30;
}


/*NLEBBE3D::NLEBBE3D(std::string str)
{
    NLEBBE3D::MAT = str;

    NLEBBE3D::NNODE = 21;
    NLEBBE3D::NELEM = 10;
    NLEBBE3D::NDOF = 6;
    NLEBBE3D::NLS = 10;
    NLEBBE3D::NEN = 3;
    NLEBBE3D::NDM = 3;

    NLEBBE3D::NODE = Eigen::MatrixXd::Zero(NLEBBE3D::NNODE, NLEBBE3D::NDM);
    NLEBBE3D::ELEM = Eigen::MatrixXd::Zero(NLEBBE3D::NELEM, NLEBBE3D::NEN + 1);

    //Nodal Information
    NLEBBE3D::NODE(0, 0) = -250;
    NLEBBE3D::NODE(0, 1) = -83.33;
    NLEBBE3D::NODE(0, 2) = 0;
    NLEBBE3D::NODE(1, 0) = -225;
    NLEBBE3D::NODE(1, 1) = -83.33;
    NLEBBE3D::NODE(1, 2) = 0;
    NLEBBE3D::NODE(2, 0) = -200;
    NLEBBE3D::NODE(2, 1) = -83.33;
    NLEBBE3D::NODE(2, 2) = 0;
    NLEBBE3D::NODE(3, 0) = -175;
    NLEBBE3D::NODE(3, 1) = -83.33;
    NLEBBE3D::NODE(3, 2) = 0;
    NLEBBE3D::NODE(4, 0) = -150;
    NLEBBE3D::NODE(4, 1) = -83.33;
    NLEBBE3D::NODE(4, 2) = 0;
    NLEBBE3D::NODE(5, 0) = -125;
    NLEBBE3D::NODE(5, 1) = -83.33;
    NLEBBE3D::NODE(5, 2) = 0;
    NLEBBE3D::NODE(6, 0) = -100;
    NLEBBE3D::NODE(6, 1) = -83.33;
    NLEBBE3D::NODE(6, 2) = 0;
    NLEBBE3D::NODE(7, 0) = -75;
    NLEBBE3D::NODE(7, 1) = -83.33;
    NLEBBE3D::NODE(7, 2) = 0;
    NLEBBE3D::NODE(8, 0) = -50;
    NLEBBE3D::NODE(8, 1) = -83.33;
    NLEBBE3D::NODE(8, 2) = 0;
    NLEBBE3D::NODE(9, 0) = -25;
    NLEBBE3D::NODE(9, 1) = -83.33;
    NLEBBE3D::NODE(9, 2) = 0;
    NLEBBE3D::NODE(10, 0) = 0;
    NLEBBE3D::NODE(10, 1) = -83.33;
    NLEBBE3D::NODE(10, 2) = 0;
    NLEBBE3D::NODE(11, 0) = 25;
    NLEBBE3D::NODE(11, 1) = -83.33;
    NLEBBE3D::NODE(11, 2) = 0;
    NLEBBE3D::NODE(12, 0) = 50;
    NLEBBE3D::NODE(12, 1) = -83.33;
    NLEBBE3D::NODE(12, 2) = 0;
    NLEBBE3D::NODE(13, 0) = 75;
    NLEBBE3D::NODE(13, 1) = -83.33;
    NLEBBE3D::NODE(13, 2) = 0;
    NLEBBE3D::NODE(14, 0) = 100;
    NLEBBE3D::NODE(14, 1) = -83.33;
    NLEBBE3D::NODE(14, 2) = 0;
    NLEBBE3D::NODE(15, 0) = 125;
    NLEBBE3D::NODE(15, 1) = -83.33;
    NLEBBE3D::NODE(15, 2) = 0;
    NLEBBE3D::NODE(16, 0) = 150;
    NLEBBE3D::NODE(16, 1) = -83.33;
    NLEBBE3D::NODE(16, 2) = 0;
    NLEBBE3D::NODE(17, 0) = 175;
    NLEBBE3D::NODE(17, 1) = -83.33;
    NLEBBE3D::NODE(17, 2) = 0;
    NLEBBE3D::NODE(18, 0) = 200;
    NLEBBE3D::NODE(18, 1) = -83.33;
    NLEBBE3D::NODE(18, 2) = 0;
    NLEBBE3D::NODE(19, 0) = 225;
    NLEBBE3D::NODE(19, 1) = -83.33;
    NLEBBE3D::NODE(19, 2) = 0;
    NLEBBE3D::NODE(20, 0) = 250;
    NLEBBE3D::NODE(20, 1) = -83.33;
    NLEBBE3D::NODE(20, 2) = 0;

    //Element connectivity
    NLEBBE3D::ELEM(0, 0) = 1;
    NLEBBE3D::ELEM(0, 1) = 1;
    NLEBBE3D::ELEM(0, 2) = 2;
    NLEBBE3D::ELEM(0, 3) = 3;

    NLEBBE3D::ELEM(1, 0) = 1;
    NLEBBE3D::ELEM(1, 1) = 3;
    NLEBBE3D::ELEM(1, 2) = 4;
    NLEBBE3D::ELEM(1, 3) = 5;


    NLEBBE3D::ELEM(2, 0) = 1;
    NLEBBE3D::ELEM(2, 1) = 5;
    NLEBBE3D::ELEM(2, 2) = 6;
    NLEBBE3D::ELEM(2, 3) = 7;

    NLEBBE3D::ELEM(3, 0) = 1;
    NLEBBE3D::ELEM(3, 1) = 7;
    NLEBBE3D::ELEM(3, 2) = 8;
    NLEBBE3D::ELEM(3, 3) = 9;

    NLEBBE3D::ELEM(4, 0) = 1;
    NLEBBE3D::ELEM(4, 1) = 9;
    NLEBBE3D::ELEM(4, 2) = 10;
    NLEBBE3D::ELEM(4, 3) = 11;

    NLEBBE3D::ELEM(5, 0) = 1;
    NLEBBE3D::ELEM(5, 1) = 11;
    NLEBBE3D::ELEM(5, 2) = 12;
    NLEBBE3D::ELEM(5, 3) = 13;

    NLEBBE3D::ELEM(6, 0) = 1;
    NLEBBE3D::ELEM(6, 1) = 13;
    NLEBBE3D::ELEM(6, 2) = 14;
    NLEBBE3D::ELEM(6, 3) = 15;

    NLEBBE3D::ELEM(7, 0) = 1;
    NLEBBE3D::ELEM(7, 1) = 15;
    NLEBBE3D::ELEM(7, 2) = 16;
    NLEBBE3D::ELEM(7, 3) = 17;

    NLEBBE3D::ELEM(8, 0) = 1;
    NLEBBE3D::ELEM(8, 1) = 17;
    NLEBBE3D::ELEM(8, 2) = 18;
    NLEBBE3D::ELEM(8, 3) = 19;

    NLEBBE3D::ELEM(9, 0) = 1;
    NLEBBE3D::ELEM(9, 1) = 19;
    NLEBBE3D::ELEM(9, 2) = 20;
    NLEBBE3D::ELEM(9, 3) = 21;

    NLEBBE3D::E = 2.4e3;
    NLEBBE3D::nu = 0.39;
    NLEBBE3D::Bp = 1;
    NLEBBE3D::Hp = 1;
    NLEBBE3D::Zx = 0;
    NLEBBE3D::Zy = 0;
    NLEBBE3D::Zz = 1;

    NLEBBE3D::DIA = 30;
}*/


//debugging inp file
/*NLEBBE3D::NLEBBE3D(int choice, std::string str)
{
    if (choice == 1)
    {
        NLEBBE3D::MAT = str;

        NLEBBE3D::NNODE = 5;
        NLEBBE3D::NELEM = 2;
        NLEBBE3D::NDOF = 6;
        NLEBBE3D::NLS = 10;
        NLEBBE3D::NEN = 3;
        NLEBBE3D::NDM = 3;

        NLEBBE3D::NODE = Eigen::MatrixXd::Zero(NLEBBE3D::NNODE, NLEBBE3D::NDM);
        NLEBBE3D::ELEM = Eigen::MatrixXd::Zero(NLEBBE3D::NELEM, NLEBBE3D::NEN + 1);

        //Nodal Information
        //Nodal Information
        NLEBBE3D::NODE(0, 0) = -2;
        NLEBBE3D::NODE(0, 1) = -1.5;
        NLEBBE3D::NODE(0, 2) = 0;
        NLEBBE3D::NODE(1, 0) = -1;
        NLEBBE3D::NODE(1, 1) = -1.5;
        NLEBBE3D::NODE(1, 2) = 0;
        NLEBBE3D::NODE(2, 0) = 0;
        NLEBBE3D::NODE(2, 1) = -1.5;
        NLEBBE3D::NODE(2, 2) = 0;
        NLEBBE3D::NODE(3, 0) = 1;
        NLEBBE3D::NODE(3, 1) = -1.5;
        NLEBBE3D::NODE(3, 2) = 0;
        NLEBBE3D::NODE(4, 0) = 2;
        NLEBBE3D::NODE(4, 1) = -1.5;
        NLEBBE3D::NODE(4, 2) = 0;

        //Element connectivity
        NLEBBE3D::ELEM(0, 0) = 1;
        NLEBBE3D::ELEM(0, 1) = 1;
        NLEBBE3D::ELEM(0, 2) = 2;
        NLEBBE3D::ELEM(0, 3) = 3;

        NLEBBE3D::ELEM(1, 0) = 1;
        NLEBBE3D::ELEM(1, 1) = 3;
        NLEBBE3D::ELEM(1, 2) = 4;
        NLEBBE3D::ELEM(1, 3) = 5;

        NLEBBE3D::E = 2.4e9;
        NLEBBE3D::nu = 0.39;
        NLEBBE3D::Bp = 1;
        NLEBBE3D::Hp = 1;
        NLEBBE3D::Zx = 0;
        NLEBBE3D::Zy = 0;
        NLEBBE3D::Zz = 1;

        NLEBBE3D::DIA = 0.1;
    }
    else if (choice == 3)
    {
        NLEBBE3D::MAT = str;

        NLEBBE3D::NNODE = 5;
        NLEBBE3D::NELEM = 2;
        NLEBBE3D::NDOF = 6;
        NLEBBE3D::NLS = 10;
        NLEBBE3D::NEN = 3;
        NLEBBE3D::NDM = 3;

        NLEBBE3D::NODE = Eigen::MatrixXd::Zero(NLEBBE3D::NNODE, NLEBBE3D::NDM);
        NLEBBE3D::ELEM = Eigen::MatrixXd::Zero(NLEBBE3D::NELEM, NLEBBE3D::NEN + 1);

        //Nodal Information
        NLEBBE3D::NODE(0, 0) = -2;
        NLEBBE3D::NODE(0, 1) = 1.5;
        NLEBBE3D::NODE(0, 2) = 0;
        NLEBBE3D::NODE(1, 0) = -1;
        NLEBBE3D::NODE(1, 1) = 1.5;
        NLEBBE3D::NODE(1, 2) = 0;
        NLEBBE3D::NODE(2, 0) = 0;
        NLEBBE3D::NODE(2, 1) = 1.5;
        NLEBBE3D::NODE(2, 2) = 0;
        NLEBBE3D::NODE(3, 0) = 1;
        NLEBBE3D::NODE(3, 1) = 1.5;
        NLEBBE3D::NODE(3, 2) = 0;
        NLEBBE3D::NODE(4, 0) = 2;
        NLEBBE3D::NODE(4, 1) = 1.5;
        NLEBBE3D::NODE(4, 2) = 0;

        //Element connectivity
        NLEBBE3D::ELEM(0, 0) = 1;
        NLEBBE3D::ELEM(0, 1) = 1;
        NLEBBE3D::ELEM(0, 2) = 2;
        NLEBBE3D::ELEM(0, 3) = 3;

        NLEBBE3D::ELEM(1, 0) = 1;
        NLEBBE3D::ELEM(1, 1) = 3;
        NLEBBE3D::ELEM(1, 2) = 4;
        NLEBBE3D::ELEM(1, 3) = 5;

        NLEBBE3D::E = 2.4e9;
        NLEBBE3D::nu = 0.39;
        NLEBBE3D::Bp = 1;
        NLEBBE3D::Hp = 1;
        NLEBBE3D::Zx = 0;
        NLEBBE3D::Zy = 0;
        NLEBBE3D::Zz = 1;

        NLEBBE3D::DIA = 0.1;
    }
    else if (choice == 2)
    {
        NLEBBE3D::MAT = str;

        NLEBBE3D::NNODE = 5;
        NLEBBE3D::NELEM = 2;
        NLEBBE3D::NDOF = 6;
        NLEBBE3D::NLS = 10;
        NLEBBE3D::NEN = 3;
        NLEBBE3D::NDM = 3;

        NLEBBE3D::NODE = Eigen::MatrixXd::Zero(NLEBBE3D::NNODE, NLEBBE3D::NDM);
        NLEBBE3D::ELEM = Eigen::MatrixXd::Zero(NLEBBE3D::NELEM, NLEBBE3D::NEN + 1);

        //Nodal Information
        //Nodal Information
        NLEBBE3D::NODE(0, 0) = -1.2;
        NLEBBE3D::NODE(0, 1) = -2;
        NLEBBE3D::NODE(0, 2) = 0;
        NLEBBE3D::NODE(1, 0) = -1.2;
        NLEBBE3D::NODE(1, 1) = -1;
        NLEBBE3D::NODE(1, 2) = 0;
        NLEBBE3D::NODE(2, 0) = -1.2;
        NLEBBE3D::NODE(2, 1) = 0;
        NLEBBE3D::NODE(2, 2) = 0;
        NLEBBE3D::NODE(3, 0) = -1.2;
        NLEBBE3D::NODE(3, 1) = 1;
        NLEBBE3D::NODE(3, 2) = 0;
        NLEBBE3D::NODE(4, 0) = -1.2;
        NLEBBE3D::NODE(4, 1) = 2;
        NLEBBE3D::NODE(4, 2) = 0;

        //Element connectivity
        NLEBBE3D::ELEM(0, 0) = 1;
        NLEBBE3D::ELEM(0, 1) = 1;
        NLEBBE3D::ELEM(0, 2) = 2;
        NLEBBE3D::ELEM(0, 3) = 3;

        NLEBBE3D::ELEM(1, 0) = 1;
        NLEBBE3D::ELEM(1, 1) = 3;
        NLEBBE3D::ELEM(1, 2) = 4;
        NLEBBE3D::ELEM(1, 3) = 5;

        NLEBBE3D::E = 2.4e9;
        NLEBBE3D::nu = 0.39;
        NLEBBE3D::Bp = 1;
        NLEBBE3D::Hp = 1;
        NLEBBE3D::Zx = 0;
        NLEBBE3D::Zy = 0;
        NLEBBE3D::Zz = 1;

        NLEBBE3D::DIA = 0.1;
    }
    else if (choice == 4)
    {
        NLEBBE3D::MAT = str;

        NLEBBE3D::NNODE = 5;
        NLEBBE3D::NELEM = 2;
        NLEBBE3D::NDOF = 6;
        NLEBBE3D::NLS = 10;
        NLEBBE3D::NEN = 3;
        NLEBBE3D::NDM = 3;

        NLEBBE3D::NODE = Eigen::MatrixXd::Zero(NLEBBE3D::NNODE, NLEBBE3D::NDM);
        NLEBBE3D::ELEM = Eigen::MatrixXd::Zero(NLEBBE3D::NELEM, NLEBBE3D::NEN + 1);

        //Nodal Information
        NLEBBE3D::NODE(0, 0) = 1.2;
        NLEBBE3D::NODE(0, 1) = -2;
        NLEBBE3D::NODE(0, 2) = 0;
        NLEBBE3D::NODE(1, 0) = 1.2;
        NLEBBE3D::NODE(1, 1) = -1;
        NLEBBE3D::NODE(1, 2) = 0;
        NLEBBE3D::NODE(2, 0) = 1.2;
        NLEBBE3D::NODE(2, 1) = 0;
        NLEBBE3D::NODE(2, 2) = 0;
        NLEBBE3D::NODE(3, 0) = 1.2;
        NLEBBE3D::NODE(3, 1) = 1;
        NLEBBE3D::NODE(3, 2) = 0;
        NLEBBE3D::NODE(4, 0) = 1.2;
        NLEBBE3D::NODE(4, 1) = 2;
        NLEBBE3D::NODE(4, 2) = 0;

        //Element connectivity
        NLEBBE3D::ELEM(0, 0) = 1;
        NLEBBE3D::ELEM(0, 1) = 1;
        NLEBBE3D::ELEM(0, 2) = 2;
        NLEBBE3D::ELEM(0, 3) = 3;

        NLEBBE3D::ELEM(1, 0) = 1;
        NLEBBE3D::ELEM(1, 1) = 3;
        NLEBBE3D::ELEM(1, 2) = 4;
        NLEBBE3D::ELEM(1, 3) = 5;

        NLEBBE3D::E = 2.4e9;
        NLEBBE3D::nu = 0.39;
        NLEBBE3D::Bp = 1;
        NLEBBE3D::Hp = 1;
        NLEBBE3D::Zx = 0;
        NLEBBE3D::Zy = 0;
        NLEBBE3D::Zz = 1;

        NLEBBE3D::DIA = 0.1;
        }
}*/


std::string NLEBBE3D::get_matmodel()
{
    return this->MAT;
}

int NLEBBE3D::get_nen()
{
    return this->NEN;
}

int NLEBBE3D::get_ndof()
{
    return this->NDOF;
}

int NLEBBE3D::get_nnode()
{
    return this->NNODE;
}

int NLEBBE3D::get_nelem()
{
    return this->NELEM;
}

int NLEBBE3D::get_nls()
{
    return this->NLS;
}

double NLEBBE3D::get_coordinates(int i, int j)
{
    return this->NODE(i, j);
}

double NLEBBE3D::get_matprop(std::string str)
{
    if (str == "E")
        return this->E;
    else if (str == "nu")
        return this->nu;
    else if (str == "Bp")
        return this->Bp;
    else if (str == "Hp")
        return this->Hp;
    else if (str == "Zx")
        return this->Zx;
    else if (str == "Zy")
        return this->Zy;
    else if (str == "Zz")
        return this->Zz;
    else if (str == "C10")
        return this->C10;
    else
        return -99;
}

int NLEBBE3D::get_connectivity(int i, int j)
{
    return this->ELEM(i, j);
}

int NLEBBE3D::get_loadnode()
{
    return this->loadnode;
}

int NLEBBE3D::get_nbeams()
{
    return this->NBEAMS;
}

double NLEBBE3D::get_diameter()
{
    return this->DIA;
}

int NLEBBE3D::get_ndim()
{
    return this->NDM;
}


double GaussIntegrationPoints(int i, int j)
{
    if (i == 1)
    {
        if (j == 1)
            return -0.577350269;
        else if (j == 2)
            return -0.577350269;
        else if (j == 3)
            return 0.577350269;
        else
            return 1.00;
    }
    else if (i == 2)
    {
        if (j == 1)
            return 0.577350269;
        else if (j == 2)
            return -0.577350269;
        else if (j == 3)
            return  0.577350269;
        else
            return 1.00;
    }
    else if (i == 3)
    {
        if (j == 1)
            return 0.577350269;
        else if (j == 2)
            return 0.577350269;
        else if (j == 3)
            return 0.577350269;
        else
            return 1.00;
    }
    else if (i == 4)
    {
        if (j == 1)
            return -0.577350269;
        else if (j == 2)
            return 0.577350269;
        else if (j == 3)
            return 0.577350269;
        else
            return 1.00;
    }
    else if (i == 5)
    {
        if (j == 1)
            return -0.577350269;
        else if (j == 2)
            return -0.577350269;
        else if (j == 3)
            return -0.577350269;
        else
            return 1.00;
    }
    else if (i == 6)
    {
        if (j == 1)
            return 0.577350269;
        else if (j == 2)
            return -0.577350269;
        else if (j == 3)
            return -0.577350269;
        else
            return 1.00;
    }
    else if (i == 7)
    {
        if (j == 1)
            return 0.577350269;
        else if (j == 2)
            return 0.577350269;
        else if (j == 3)
            return -0.577350269;
        else
            return 1.00;
    }
    else if (i == 8)
    {
        if (j == 1)
            return -0.577350269;
        else if (j == 2)
            return 0.577350269;
        else if (j == 3)
            return -0.577350269;
        else
            return 1.00;
    }
}

void NLEBBE3D::UpdateNodes(NLEBBE3D* EBBE3D, Eigen::VectorXd* GU, int m)
{
    int temp = 0;
    for (int i = m; i > 0; i--)
        temp += EBBE3D[i].get_nnode() * EBBE3D[i].get_ndof();
    for (int i = 0; i < EBBE3D[m].get_nnode(); i++)
    {
        for (int k = 0; k < 3; k++)
        {
            EBBE3D[m].NODE(i, k) = EBBE3D[m].get_coordinates(i, k) + (*GU)(EBBE3D[m].get_ndof() * i + k + temp);
            //std::cout << EBBE3D[m].NODE(i, k);
        }
    }
}

void NLEBBE3D::RKt_Orthotropic(double D[14], double X[3][3], double U[3][6]
    , double** T, double* R)
{
    double v[1525];
    int i43, i227, ii386, b105, b106, b278, b444, b704;
    /* 1 = E_1 */
    v[1] = D[0];
    /* 2 = E_2 */
    v[2] = D[1];
    /* 3 = E_3 */
    v[3] = D[2];
    /* 4 = \[Nu]_12 */
    v[4] = D[3];
    /* 221 = \[Nu]_21 */
    v[221] = (v[2] * v[4]) / v[1];
    /* 5 = \[Nu]_13 */
    v[5] = D[4];
    /* 222 = \[Nu]_31 */
    v[222] = (v[3] * v[5]) / v[1];
    /* 6 = \[Nu]_23 */
    v[6] = D[5];
    /* 223 = \[Nu]_32 */
    v[223] = (v[3] * v[6]) / v[2];
    v[225] = 1e0 - v[221] * v[4] - v[222] * v[5] - v[223] * v[6] - 2e0 * v[222] * v[4] * v[6];
    /* 224 = \[CapitalDelta] */
    v[224] = v[225] / (v[1] * v[2] * v[3]);
    /* 7 = G_12 */
    v[7] = D[6];
    /* 8 = G_13 */
    v[8] = D[7];
    /* 9 = G_23 */
    v[9] = D[8];
    /* 10 = Bp */
    v[10] = D[9];
    /* 11 = Hp */
    v[11] = D[10];
    /* 12 = Zx */
    v[12] = D[11];
    /* 13 = Zy */
    v[13] = D[12];
    /* 14 = Zz */
    v[14] = D[13];
    /* 15 = XIO_1|1 */
    v[15] = X[0][0];
    /* 16 = XIO_1|2 */
    v[16] = X[0][1];
    /* 17 = XIO_1|3 */
    v[17] = X[0][2];
    /* 18 = XIO_2|1 */
    v[18] = X[1][0];
    /* 19 = XIO_2|2 */
    v[19] = X[1][1];
    /* 20 = XIO_2|3 */
    v[20] = X[1][2];
    /* 21 = XIO_3|1 */
    v[21] = X[2][0];
    /* 75 = X0_1;\[Xi];\[Xi] */
    v[75] = v[15] - 2e0 * v[18] + v[21];
    /* 22 = XIO_3|2 */
    v[22] = X[2][1];
    /* 74 = X0_2;\[Xi];\[Xi] */
    v[74] = v[16] - 2e0 * v[19] + v[22];
    /* 23 = XIO_3|3 */
    v[23] = X[2][2];
    /* 73 = X0_3;\[Xi];\[Xi] */
    v[73] = v[17] - 2e0 * v[20] + v[23];
    /* 24 = peIO_1|1 */
    v[24] = U[0][0];
    /* 25 = peIO_1|2 */
    v[25] = U[0][1];
    /* 26 = peIO_1|3 */
    v[26] = U[0][2];
    /* 27 = peIO_1|4 */
    v[27] = U[0][3];
    /* 28 = peIO_1|5 */
    v[28] = U[0][4];
    /* 29 = peIO_1|6 */
    v[29] = U[0][5];
    /* 30 = peIO_2|1 */
    v[30] = U[1][0];
    /* 31 = peIO_2|2 */
    v[31] = U[1][1];
    /* 32 = peIO_2|3 */
    v[32] = U[1][2];
    /* 33 = peIO_2|4 */
    v[33] = U[1][3];
    /* 34 = peIO_2|5 */
    v[34] = U[1][4];
    /* 35 = peIO_2|6 */
    v[35] = U[1][5];
    /* 36 = peIO_3|1 */
    v[36] = U[2][0];
    /* 37 = peIO_3|2 */
    v[37] = U[2][1];
    /* 38 = peIO_3|3 */
    v[38] = U[2][2];
    /* 39 = peIO_3|4 */
    v[39] = U[2][3];
    /* 40 = peIO_3|5 */
    v[40] = U[2][4];
    /* 41 = peIO_3|6 */
    v[41] = U[2][5];
    v[1149] = v[24];
    v[1150] = v[25];
    v[1151] = v[26];
    v[1152] = v[27];
    v[1153] = v[28];
    v[1154] = v[29];
    v[1155] = v[30];
    v[1156] = v[31];
    v[1157] = v[32];
    v[1158] = v[33];
    v[1159] = v[34];
    v[1160] = v[35];
    v[1161] = v[36];
    v[1162] = v[37];
    v[1163] = v[38];
    v[1164] = v[39];
    v[1165] = v[40];
    v[1166] = v[41];
    v[387] = 0e0;/*debug*/
    /* 42 = Le */
    v[42] = sqrt(Power(v[15] - v[18], 2) + Power(v[16] - v[19], 2) + Power(v[17] - v[20], 2));
    for (i43 = 1; i43 <= 8; i43++) {
        v[43] = i43;
        /* 44 = \[Xi] */
        v[44] = GaussIntegrationPoints(i43, 1);
        v[58] = v[44] / 2e0;
        v[51] = 1e0 + v[44];
        /* 59 = Nh_3;\[Xi] */
        v[59] = v[51] / 2e0 + v[58];
        v[49] = -1e0 + v[44];
        /* 57 = Nh_2;\[Xi] */
        v[57] = -v[49] - v[51];
        /* 56 = Nh_1;\[Xi] */
        v[56] = v[49] / 2e0 + v[58];
        /* 150 = u\[Phi]_6;\[Xi] */
        v[150] = v[29] * v[56] + v[35] * v[57] + v[41] * v[59];
        /* 147 = u\[Phi]_5;\[Xi] */
        v[147] = v[28] * v[56] + v[34] * v[57] + v[40] * v[59];
        /* 145 = u\[Phi]_4;\[Xi] */
        v[145] = v[27] * v[56] + v[33] * v[57] + v[39] * v[59];
        /* 144 = u\[Phi]_3;\[Xi] */
        v[144] = v[26] * v[56] + v[32] * v[57] + v[38] * v[59];
        /* 143 = u\[Phi]_2;\[Xi] */
        v[143] = v[25] * v[56] + v[31] * v[57] + v[37] * v[59];
        /* 142 = u\[Phi]_1;\[Xi] */
        v[142] = v[24] * v[56] + v[30] * v[57] + v[36] * v[59];
        /* 62 = X0_3;\[Xi] */
        v[62] = v[17] * v[56] + v[20] * v[57] + v[23] * v[59];
        /* 61 = X0_2;\[Xi] */
        v[61] = v[16] * v[56] + v[19] * v[57] + v[22] * v[59];
        /* 60 = X0_1;\[Xi] */
        v[60] = v[15] * v[56] + v[18] * v[57] + v[21] * v[59];
        v[140] = 2e0 * v[62] * v[73] + 2e0 * v[61] * v[74] + 2e0 * v[60] * v[75];
        v[77] = (v[60] * v[60]) + (v[61] * v[61]) + (v[62] * v[62]);
        v[141] = -0.5e0 * v[140] / (v[77] * sqrt(v[77]));
        v[76] = 1e0 / sqrt(v[77]);
        v[78] = v[141];
        v[64] = v[76];
        /* 81 = t1_3;\[Xi] */
        v[81] = v[64] * v[73] + v[62] * v[78];
        /* 80 = t1_2;\[Xi] */
        v[80] = v[64] * v[74] + v[61] * v[78];
        /* 82 = t2_1;\[Xi] */
        v[82] = -(v[14] * v[80]) + v[13] * v[81];
        /* 79 = t1_1;\[Xi] */
        v[79] = v[64] * v[75] + v[60] * v[78];
        /* 84 = t2_3;\[Xi] */
        v[84] = -(v[13] * v[79]) + v[12] * v[80];
        /* 83 = t2_2;\[Xi] */
        v[83] = v[14] * v[79] - v[12] * v[81];
        /* 45 = \[Eta] */
        v[45] = GaussIntegrationPoints(i43, 2);
        v[266] = (v[10] * v[45] * v[84]) / 2e0;
        v[263] = (v[10] * v[45] * v[83]) / 2e0;
        v[260] = (v[10] * v[45] * v[82]) / 2e0;
        /* 46 = \[Zeta] */
        v[46] = GaussIntegrationPoints(i43, 3);
        /* 47 = wgp */
        v[47] = GaussIntegrationPoints(i43, 4);
        /* 48 = Nh_1 */
        v[48] = (v[44] * v[49]) / 2e0;
        /* 50 = Nh_2 */
        v[50] = -(v[49] * v[51]);
        /* 52 = Nh_3 */
        v[52] = (v[44] * v[51]) / 2e0;
        /* 53 = X0_1 */
        v[53] = v[15] * v[48] + v[18] * v[50] + v[21] * v[52];
        /* 54 = X0_2 */
        v[54] = v[16] * v[48] + v[19] * v[50] + v[22] * v[52];
        /* 55 = X0_3 */
        v[55] = v[17] * v[48] + v[20] * v[50] + v[23] * v[52];
        /* 63 = t1_1 */
        v[63] = v[60] * v[64];
        /* 65 = t1_2 */
        v[65] = v[61] * v[64];
        /* 66 = t1_3 */
        v[66] = v[62] * v[64];
        /* 67 = t2_1 */
        v[67] = -(v[14] * v[65]) + v[13] * v[66];
        /* 68 = t2_2 */
        v[68] = v[14] * v[63] - v[12] * v[66];
        /* 87 = t3_3;\[Xi] */
        v[87] = v[68] * v[79] - v[67] * v[80] - v[65] * v[82] + v[63] * v[83];
        v[267] = (v[11] * v[46] * v[87]) / 2e0;
        v[273] = v[266] + v[267];
        /* 69 = t2_3 */
        v[69] = -(v[13] * v[63]) + v[12] * v[65];
        /* 86 = t3_2;\[Xi] */
        v[86] = -(v[69] * v[79]) + v[67] * v[81] + v[66] * v[82] - v[63] * v[84];
        v[264] = (v[11] * v[46] * v[86]) / 2e0;
        v[271] = v[263] + v[264];
        /* 85 = t3_1;\[Xi] */
        v[85] = v[69] * v[80] - v[68] * v[81] - v[66] * v[83] + v[65] * v[84];
        v[261] = (v[11] * v[46] * v[85]) / 2e0;
        v[269] = v[260] + v[261];
        /* 70 = t3_1 */
        v[70] = -(v[66] * v[68]) + v[65] * v[69];
        v[253] = (v[10] * v[45] * v[67]) / 2e0 + (v[11] * v[46] * v[70]) / 2e0;
        /* 71 = t3_2 */
        v[71] = v[66] * v[67] - v[63] * v[69];
        v[251] = (v[10] * v[45] * v[68]) / 2e0 + (v[11] * v[46] * v[71]) / 2e0;
        /* 72 = t3_3 */
        v[72] = -(v[65] * v[67]) + v[63] * v[68];
        v[249] = (v[10] * v[45] * v[69]) / 2e0 + (v[11] * v[46] * v[72]) / 2e0;
        /* 88 = Je_1|1 */
        v[88] = v[260] + v[261] + v[60];
        /* 89 = Je_1|2 */
        v[89] = (v[10] * v[67]) / 2e0;
        /* 90 = Je_1|3 */
        v[90] = (v[11] * v[70]) / 2e0;
        /* 91 = Je_2|1 */
        v[91] = v[263] + v[264] + v[61];
        /* 92 = Je_2|2 */
        v[92] = (v[10] * v[68]) / 2e0;
        v[202] = -(v[89] * v[91]) + v[88] * v[92];
        /* 93 = Je_2|3 */
        v[93] = (v[11] * v[71]) / 2e0;
        v[201] = v[90] * v[91] - v[88] * v[93];
        v[200] = -(v[90] * v[92]) + v[89] * v[93];
        /* 94 = Je_3|1 */
        v[94] = v[266] + v[267] + v[62];
        /* 95 = Je_3|2 */
        v[95] = (v[10] * v[69]) / 2e0;
        v[208] = -(v[92] * v[94]) + v[91] * v[95];
        v[205] = v[89] * v[94] - v[88] * v[95];
        /* 96 = Je_3|3 */
        v[96] = (v[11] * v[72]) / 2e0;
        v[207] = v[93] * v[94] - v[91] * v[96];
        v[206] = -(v[93] * v[95]) + v[92] * v[96];
        v[204] = -(v[90] * v[94]) + v[88] * v[96];
        v[203] = v[90] * v[95] - v[89] * v[96];
        /* 97 = Jed */
        v[97] = -(v[90] * v[92] * v[94]) + v[89] * v[93] * v[94] + v[90] * v[91] * v[95] - v[88] * v[93] * v[95] - v[89] * v[91] * v[96]
            + v[88] * v[92] * v[96];
        v[247] = (v[207] * v[63]) / v[97] + (v[204] * v[65]) / v[97] + (v[201] * v[66]) / v[97];
        v[245] = (v[208] * v[63]) / v[97] + (v[205] * v[65]) / v[97] + (v[202] * v[66]) / v[97];
        v[236] = (v[206] * v[63]) / v[97] + (v[203] * v[65]) / v[97] + (v[200] * v[66]) / v[97];
        /* 98 = u\[Phi]_1 */
        v[98] = v[24] * v[48] + v[30] * v[50] + v[36] * v[52];
        /* 99 = u\[Phi]_2 */
        v[99] = v[25] * v[48] + v[31] * v[50] + v[37] * v[52];
        /* 100 = u\[Phi]_3 */
        v[100] = v[26] * v[48] + v[32] * v[50] + v[38] * v[52];
        /* 101 = u\[Phi]_4 */
        v[101] = v[27] * v[48] + v[33] * v[50] + v[39] * v[52];
        v[146] = 2e0 * v[101] * v[145];
        v[117] = (v[101] * v[101]);
        /* 102 = u\[Phi]_5 */
        v[102] = v[28] * v[48] + v[34] * v[50] + v[40] * v[52];
        v[488] = v[102] * v[145] + v[101] * v[147];
        v[148] = 2e0 * v[102] * v[147];
        v[149] = v[146] + v[148];
        v[107] = (v[102] * v[102]);
        v[139] = v[107] + v[117];
        /* 103 = u\[Phi]_6 */
        v[103] = v[29] * v[48] + v[35] * v[50] + v[41] * v[52];
        v[494] = v[103] * v[147] + v[102] * v[150];
        v[491] = v[103] * v[145] + v[101] * v[150];
        v[151] = 2e0 * v[103] * v[150];
        /* 154 = ff_;\[Xi] */
        v[154] = v[146] + v[148] + v[151];
        v[153] = v[148] + v[151];
        v[152] = v[146] + v[151];
        v[108] = (v[103] * v[103]);
        v[134] = v[108] + v[117];
        v[127] = v[107] + v[108];
        /* 104 = ff */
        v[104] = v[107] + v[108] + v[117];
        b105 = v[104] < 0.1e-7;
        v[105] = b105;
        b106 = b105;
        v[106] = b106;
        if (b106) {
            v[161] = v[104] * v[154];
            v[160] = -30e0 + v[104];
            v[162] = v[154] * v[160] + v[161];
            v[155] = -20e0 + v[104];
            v[156] = v[154] * v[155] + v[161];
            v[112] = 120e0 + v[104] * v[155];
            v[159] = -6e0 * v[112] * v[150] - 6e0 * v[103] * v[156];
            v[158] = 6e0 * v[112] * v[147] + 6e0 * v[102] * v[156];
            v[157] = -6e0 * v[112] * v[145] - 6e0 * v[101] * v[156];
            v[124] = -6e0 * v[101] * v[112];
            v[121] = 6e0 * v[102] * v[112];
            v[115] = -6e0 * v[103] * v[112];
            v[110] = 360e0 + v[104] * v[160];
            v[165] = v[102] * v[110] * v[145] + v[101] * v[110] * v[147] + v[101] * v[102] * v[162];
            v[164] = v[103] * v[110] * v[145] + v[101] * v[110] * v[150] + v[101] * v[103] * v[162];
            v[163] = v[103] * v[110] * v[147] + v[102] * v[110] * v[150] + v[102] * v[103] * v[162];
            v[123] = v[102] * v[103] * v[110];
            v[120] = v[101] * v[103] * v[110];
            v[114] = v[101] * v[102] * v[110];
            /* 166 = v_1|1;\[Xi] */
            v[166] = (-1e0 / 720e0) * (v[110] * v[153]) - (v[127] * v[162]) / 720e0;
            /* 109 = v_1|1 */
            v[109] = 1e0 - (v[110] * v[127]) / 720e0;
            /* 167 = v_1|2;\[Xi] */
            v[167] = (v[159] + v[165]) / 720e0;
            /* 111 = v_1|2 */
            v[111] = (v[114] + v[115]) / 720e0;
            /* 168 = v_1|3;\[Xi] */
            v[168] = (v[158] + v[164]) / 720e0;
            /* 113 = v_1|3 */
            v[113] = (v[120] + v[121]) / 720e0;
            /* 169 = v_2|1;\[Xi] */
            v[169] = (-v[159] + v[165]) / 720e0;
            /* 116 = v_2|1 */
            v[116] = (v[114] - v[115]) / 720e0;
            /* 170 = v_2|2;\[Xi] */
            v[170] = (-1e0 / 720e0) * (v[110] * v[152]) - (v[134] * v[162]) / 720e0;
            /* 118 = v_2|2 */
            v[118] = 1e0 - (v[110] * v[134]) / 720e0;
            /* 171 = v_2|3;\[Xi] */
            v[171] = (v[157] + v[163]) / 720e0;
            /* 119 = v_2|3 */
            v[119] = (v[123] + v[124]) / 720e0;
            /* 172 = v_3|1;\[Xi] */
            v[172] = (-v[158] + v[164]) / 720e0;
            /* 122 = v_3|1 */
            v[122] = (v[120] - v[121]) / 720e0;
            /* 173 = v_3|2;\[Xi] */
            v[173] = (-v[157] + v[163]) / 720e0;
            /* 125 = v_3|2 */
            v[125] = (v[123] - v[124]) / 720e0;
            /* 174 = v_3|3;\[Xi] */
            v[174] = (-1e0 / 720e0) * (v[110] * v[149]) - (v[139] * v[162]) / 720e0;
            /* 126 = v_3|3 */
            v[126] = 1e0 - (v[110] * v[139]) / 720e0;
        }
        else {
            v[185] = 1e0 / (v[104] * v[104]);
            v[175] = 1e0 / sqrt(v[104]);
            v[176] = (v[154] * v[175]) / 2e0;
            v[128] = sqrt(v[104]);
            v[345] = 1e0 / (v[128] * v[128]);
            v[178] = cos(v[128]);
            v[179] = v[176] * v[178];
            v[177] = -(v[176] * v[345]);
            v[759] = -(v[150] * v[175]) - v[103] * v[177];
            v[758] = v[147] * v[175] + v[102] * v[177];
            v[757] = -(v[145] * v[175]) - v[101] * v[177];
            v[131] = v[175];
            v[130] = sin(v[128]);
            v[770] = -(v[130] * v[150]) - v[103] * v[179];
            v[769] = v[130] * v[147] + v[102] * v[179];
            v[768] = -(v[130] * v[145]) - v[101] * v[179];
            v[183] = -(v[130] * v[176]);
            v[767] = v[101] * v[102] * v[183] * v[185];
            v[764] = v[101] * v[103] * v[183] * v[185];
            v[761] = v[102] * v[103] * v[183] * v[185];
            v[182] = -(v[130] * v[131] * v[150]) - v[103] * v[130] * v[177] - v[103] * v[131] * v[179];
            v[181] = v[130] * v[131] * v[147] + v[102] * v[130] * v[177] + v[102] * v[131] * v[179];
            v[180] = -(v[130] * v[131] * v[145]) - v[101] * v[130] * v[177] - v[101] * v[131] * v[179];
            v[137] = -(v[101] * v[130] * v[131]);
            v[135] = v[102] * v[130] * v[131];
            v[132] = -(v[103] * v[130] * v[131]);
            v[129] = -1e0 + v[178];
            v[766] = v[101] * v[129] * v[147] * v[185];
            v[765] = v[102] * v[129] * v[145] * v[185];
            v[763] = v[101] * v[129] * v[150] * v[185];
            v[762] = v[103] * v[129] * v[145] * v[185];
            v[760] = v[102] * v[129] * v[150] * v[185];
            v[581] = v[129] * v[145] * v[185] + v[101] * v[183] * v[185];
            v[579] = v[129] * v[147] * v[185] + v[102] * v[183] * v[185];
            v[559] = v[129] * v[150] * v[185] + v[103] * v[183] * v[185];
            v[187] = -((v[102] * v[129] * v[145]) / v[104]) - (v[101] * v[129] * v[147]) / v[104] - (v[101] * v[102] * v[183])
                / v[104] + v[101] * v[102] * v[129] * v[154] * v[185];
            v[186] = -((v[103] * v[129] * v[145]) / v[104]) - (v[101] * v[129] * v[150]) / v[104] - (v[101] * v[103] * v[183])
                / v[104] + v[101] * v[103] * v[129] * v[154] * v[185];
            v[184] = -((v[103] * v[129] * v[147]) / v[104]) - (v[102] * v[129] * v[150]) / v[104] - (v[102] * v[103] * v[183])
                / v[104] + v[102] * v[103] * v[129] * v[154] * v[185];
            v[138] = -((v[102] * v[103] * v[129]) / v[104]);
            v[136] = -((v[101] * v[103] * v[129]) / v[104]);
            v[133] = -((v[101] * v[102] * v[129]) / v[104]);
            /* 166 = v_1|1;\[Xi] */
            v[166] = (v[129] * v[153]) / v[104] + (v[127] * v[183]) / v[104] - v[127] * v[129] * v[154] * v[185];
            /* 109 = v_1|1 */
            v[109] = 1e0 + (v[127] * v[129]) / v[104];
            /* 167 = v_1|2;\[Xi] */
            v[167] = v[182] + v[187];
            /* 111 = v_1|2 */
            v[111] = v[132] + v[133];
            /* 168 = v_1|3;\[Xi] */
            v[168] = v[181] + v[186];
            /* 113 = v_1|3 */
            v[113] = v[135] + v[136];
            /* 169 = v_2|1;\[Xi] */
            v[169] = -v[182] + v[187];
            /* 116 = v_2|1 */
            v[116] = -v[132] + v[133];
            /* 170 = v_2|2;\[Xi] */
            v[170] = (v[129] * v[152]) / v[104] + (v[134] * v[183]) / v[104] - v[129] * v[134] * v[154] * v[185];
            /* 118 = v_2|2 */
            v[118] = 1e0 + (v[129] * v[134]) / v[104];
            /* 171 = v_2|3;\[Xi] */
            v[171] = v[180] + v[184];
            /* 119 = v_2|3 */
            v[119] = v[137] + v[138];
            /* 172 = v_3|1;\[Xi] */
            v[172] = -v[181] + v[186];
            /* 122 = v_3|1 */
            v[122] = -v[135] + v[136];
            /* 173 = v_3|2;\[Xi] */
            v[173] = -v[180] + v[184];
            /* 125 = v_3|2 */
            v[125] = -v[137] + v[138];
            /* 174 = v_3|3;\[Xi] */
            v[174] = (v[129] * v[149]) / v[104] + (v[139] * v[183]) / v[104] - v[129] * v[139] * v[154] * v[185];
            /* 126 = v_3|3 */
            v[126] = 1e0 + (v[129] * v[139]) / v[104];
        };
        /* 214 = Ft_3|3 */
        v[214] = v[122] * v[70] + v[125] * v[71] + v[126] * v[72];
        /* 213 = Ft_3|2 */
        v[213] = v[122] * v[67] + v[125] * v[68] + v[126] * v[69];
        /* 211 = Ft_2|3 */
        v[211] = v[116] * v[70] + v[118] * v[71] + v[119] * v[72];
        /* 210 = Ft_2|2 */
        v[210] = v[116] * v[67] + v[118] * v[68] + v[119] * v[69];
        /* 199 = Ft_1|3 */
        v[199] = v[109] * v[70] + v[111] * v[71] + v[113] * v[72];
        /* 198 = Ft_1|2 */
        v[198] = v[109] * v[67] + v[111] * v[68] + v[113] * v[69];
        v[239] = v[198] * v[199] + v[210] * v[211] + v[213] * v[214];
        /* 188 = Fref_1|1 */
        v[188] = v[142] + v[60] + (v[10] * v[45] * (v[166] * v[67] + v[167] * v[68] + v[168] * v[69] + v[109] * v[82] + v[111] * v[83]
            + v[113] * v[84])) / 2e0 + (v[11] * v[46] * (v[166] * v[70] + v[167] * v[71] + v[168] * v[72] + v[109] * v[85] + v[111] * v[86]
                + v[113] * v[87])) / 2e0;
        /* 189 = Fref_1|2 */
        v[189] = (v[10] * v[198]) / 2e0;
        /* 190 = Fref_1|3 */
        v[190] = (v[11] * v[199]) / 2e0;
        /* 191 = Fref_2|1 */
        v[191] = v[143] + v[61] + (v[10] * v[45] * (v[169] * v[67] + v[170] * v[68] + v[171] * v[69] + v[116] * v[82] + v[118] * v[83]
            + v[119] * v[84])) / 2e0 + (v[11] * v[46] * (v[169] * v[70] + v[170] * v[71] + v[171] * v[72] + v[116] * v[85] + v[118] * v[86]
                + v[119] * v[87])) / 2e0;
        /* 192 = Fref_2|2 */
        v[192] = (v[10] * v[210]) / 2e0;
        /* 193 = Fref_2|3 */
        v[193] = (v[11] * v[211]) / 2e0;
        /* 194 = Fref_3|1 */
        v[194] = v[144] + v[62] + (v[10] * v[45] * (v[172] * v[67] + v[173] * v[68] + v[174] * v[69] + v[122] * v[82] + v[125] * v[83]
            + v[126] * v[84])) / 2e0 + (v[11] * v[46] * (v[172] * v[70] + v[173] * v[71] + v[174] * v[72] + v[122] * v[85] + v[125] * v[86]
                + v[126] * v[87])) / 2e0;
        /* 195 = Fref_3|2 */
        v[195] = (v[10] * v[213]) / 2e0;
        /* 196 = Fref_3|3 */
        v[196] = (v[11] * v[214]) / 2e0;
        /* 197 = Ft_1|1 */
        v[197] = v[66] * ((v[188] * v[200]) / v[97] + (v[189] * v[201]) / v[97] + (v[190] * v[202]) / v[97]) + v[65] * (
            (v[188] * v[203]) / v[97] + (v[189] * v[204]) / v[97] + (v[190] * v[205]) / v[97]) + v[63] * ((v[188] * v[206]) / v[97] +
                (v[189] * v[207]) / v[97] + (v[190] * v[208]) / v[97]);
        /* 209 = Ft_2|1 */
        v[209] = v[66] * ((v[191] * v[200]) / v[97] + (v[192] * v[201]) / v[97] + (v[193] * v[202]) / v[97]) + v[65] * (
            (v[191] * v[203]) / v[97] + (v[192] * v[204]) / v[97] + (v[193] * v[205]) / v[97]) + v[63] * ((v[191] * v[206]) / v[97] +
                (v[192] * v[207]) / v[97] + (v[193] * v[208]) / v[97]);
        /* 212 = Ft_3|1 */
        v[212] = v[66] * ((v[194] * v[200]) / v[97] + (v[195] * v[201]) / v[97] + (v[196] * v[202]) / v[97]) + v[65] * (
            (v[194] * v[203]) / v[97] + (v[195] * v[204]) / v[97] + (v[196] * v[205]) / v[97]) + v[63] * ((v[194] * v[206]) / v[97] +
                (v[195] * v[207]) / v[97] + (v[196] * v[208]) / v[97]);
        v[231] = -1e0 + (v[197] * v[197]) + (v[209] * v[209]) + (v[212] * v[212]);
        v[230] = v[197] * v[199] + v[209] * v[211] + v[212] * v[214];
        v[229] = v[197] * v[198] + v[209] * v[210] + v[212] * v[213];
        /* 215 = Et_1|1 */
        v[215] = v[231] / 2e0;
        /* 216 = [Et_2|1][Et_1|2] */
        v[216] = v[229] / 2e0;
        /* 217 = [Et_3|1][Et_1|3] */
        v[217] = v[230] / 2e0;
        /* 218 = Et_2|2 */
        v[218] = (-1e0 + (v[198] * v[198]) + (v[210] * v[210]) + (v[213] * v[213])) / 2e0;
        /* 219 = [Et_3|2][Et_2|3] */
        v[219] = v[239] / 2e0;
        /* 220 = Et_3|3 */
        v[220] = (-1e0 + (v[199] * v[199]) + (v[211] * v[211]) + (v[214] * v[214])) / 2e0;
        /* 226 = W */
        v[226] = ((v[215] * v[215]) * v[225]) / (2e0 * v[2] * v[224] * v[3]) + 2e0 * ((v[216] * v[216]) * v[7] + (v[217] * v[217]
            ) * v[8] + (v[219] * v[219]) * v[9]);
        /* 232 = \[OverBracket]_Ft_3|1(W|W) */
        v[232] = (v[212] * v[225] * v[231]) / (2e0 * v[2] * v[224] * v[3]) + 2e0 * ((v[213] * v[229] * v[7]) / 2e0 +
            (v[214] * v[230] * v[8]) / 2e0);
        /* 233 = \[OverBracket]_Ft_2|1(W|W) */
        v[233] = (v[209] * v[225] * v[231]) / (2e0 * v[2] * v[224] * v[3]) + 2e0 * ((v[210] * v[229] * v[7]) / 2e0 +
            (v[211] * v[230] * v[8]) / 2e0);
        /* 234 = \[OverBracket]_Ft_1|1(W|W) */
        v[234] = (v[197] * v[225] * v[231]) / (2e0 * v[2] * v[224] * v[3]) + 2e0 * ((v[198] * v[229] * v[7]) / 2e0 +
            (v[199] * v[230] * v[8]) / 2e0);
        /* 235 = \[OverBracket]_Fref_3|1(Ft|W)_3|1 */
        v[235] = v[232] * v[236];
        /* 237 = \[OverBracket]_Fref_2|1(Ft|W)_2|1 */
        v[237] = v[233] * v[236];
        /* 238 = \[OverBracket]_Fref_1|1(Ft|W)_1|1 */
        v[238] = v[234] * v[236];
        /* 240 = \[OverBracket]_Ft_3|3(Fref|W)_3|3 */
        v[240] = (v[11] * v[232] * v[245]) / 2e0 + 2e0 * ((v[212] * v[230] * v[8]) / 2e0 + (v[213] * v[239] * v[9]) / 2e0);
        /* 241 = \[OverBracket]_Ft_3|2(Fref|W)_3|2 */
        v[241] = (v[10] * v[232] * v[247]) / 2e0 + 2e0 * ((v[212] * v[229] * v[7]) / 2e0 + (v[214] * v[239] * v[9]) / 2e0);
        /* 242 = \[OverBracket]_R_3|3;\[Xi](Fref|W)_3|1 */
        v[242] = v[235] * v[249];
        /* 243 = \[OverBracket]_R_3|2;\[Xi](Fref|W)_3|1 */
        v[243] = v[235] * v[251];
        /* 244 = \[OverBracket]_R_3|1;\[Xi](Fref|W)_3|1 */
        v[244] = v[235] * v[253];
        /* 246 = \[OverBracket]_Ft_2|3(Fref|W)_2|3 */
        v[246] = (v[11] * v[233] * v[245]) / 2e0 + 2e0 * ((v[209] * v[230] * v[8]) / 2e0 + (v[210] * v[239] * v[9]) / 2e0);
        /* 248 = \[OverBracket]_Ft_2|2(Fref|W)_2|2 */
        v[248] = (v[10] * v[233] * v[247]) / 2e0 + 2e0 * ((v[209] * v[229] * v[7]) / 2e0 + (v[211] * v[239] * v[9]) / 2e0);
        /* 250 = \[OverBracket]_R_2|3;\[Xi](Fref|W)_2|1 */
        v[250] = v[237] * v[249];
        /* 252 = \[OverBracket]_R_2|2;\[Xi](Fref|W)_2|1 */
        v[252] = v[237] * v[251];
        /* 254 = \[OverBracket]_R_2|1;\[Xi](Fref|W)_2|1 */
        v[254] = v[237] * v[253];
        /* 255 = \[OverBracket]_Ft_1|3(Fref|W)_1|3 */
        v[255] = (v[11] * v[234] * v[245]) / 2e0 + 2e0 * ((v[197] * v[230] * v[8]) / 2e0 + (v[198] * v[239] * v[9]) / 2e0);
        /* 256 = \[OverBracket]_Ft_1|2(Fref|W)_1|2 */
        v[256] = (v[10] * v[234] * v[247]) / 2e0 + 2e0 * ((v[197] * v[229] * v[7]) / 2e0 + (v[199] * v[239] * v[9]) / 2e0);
        /* 257 = \[OverBracket]_R_1|3;\[Xi](Fref|W)_1|1 */
        v[257] = v[238] * v[249];
        /* 258 = \[OverBracket]_R_1|2;\[Xi](Fref|W)_1|1 */
        v[258] = v[238] * v[251];
        /* 259 = \[OverBracket]_R_1|1;\[Xi](Fref|W)_1|1 */
        v[259] = v[238] * v[253];
        /* 262 = \[OverBracket]_R_1|1(Ft|W)_1|3 */
        v[262] = v[238] * v[269] + v[256] * v[67] + v[255] * v[70];
        /* 265 = \[OverBracket]_R_1|2(Ft|W)_1|3 */
        v[265] = v[238] * v[271] + v[256] * v[68] + v[255] * v[71];
        /* 268 = \[OverBracket]_R_1|3(Ft|W)_1|3 */
        v[268] = v[238] * v[273] + v[256] * v[69] + v[255] * v[72];
        /* 270 = \[OverBracket]_R_2|1(Ft|W)_2|3 */
        v[270] = v[237] * v[269] + v[248] * v[67] + v[246] * v[70];
        /* 272 = \[OverBracket]_R_2|2(Ft|W)_2|3 */
        v[272] = v[237] * v[271] + v[248] * v[68] + v[246] * v[71];
        /* 274 = \[OverBracket]_R_2|3(Ft|W)_2|3 */
        v[274] = v[237] * v[273] + v[248] * v[69] + v[246] * v[72];
        /* 275 = \[OverBracket]_R_3|1(Ft|W)_3|3 */
        v[275] = v[235] * v[269] + v[241] * v[67] + v[240] * v[70];
        /* 276 = \[OverBracket]_R_3|2(Ft|W)_3|3 */
        v[276] = v[235] * v[271] + v[241] * v[68] + v[240] * v[71];
        /* 277 = \[OverBracket]_R_3|3(Ft|W)_3|3 */
        v[277] = v[235] * v[273] + v[241] * v[69] + v[240] * v[72];
        b278 = b105;
        v[278] = b278;
        if (b278) {
            v[319] = v[110] * v[145] + v[101] * v[162];
            v[318] = v[110] * v[147] + v[102] * v[162];
            v[313] = v[110] * v[150] + v[103] * v[162];
            v[305] = v[258] / 720e0;
            v[304] = v[254] / 720e0;
            v[301] = v[270] / 720e0;
            v[300] = v[265] / 720e0;
            v[297] = v[257] / 720e0;
            v[296] = v[244] / 720e0;
            v[293] = v[275] / 720e0;
            v[292] = v[268] / 720e0;
            v[287] = v[250] / 720e0;
            v[286] = v[243] / 720e0;
            v[283] = v[276] / 720e0;
            v[282] = v[274] / 720e0;
            /* 279 = \[OverBracket]_\[Yen]_139(v|W)_3|3;\[Xi] */
            v[279] = (-1e0 / 720e0) * (v[162] * v[242]) - (v[110] * v[277]) / 720e0;
            /* 280 = \[OverBracket]_\[Yen]_149(v|W)_3|3;\[Xi] */
            v[280] = (-1e0 / 720e0) * (v[110] * v[242]);
            /* 281 = \[OverBracket]_\[Yen]_123(v|W)_2|3 */
            v[281] = v[282] + v[283];
            v[524] = v[103] * v[281];
            v[521] = v[102] * v[281];
            /* 284 = \[OverBracket]_\[Yen]_124(v|W)_2|3 */
            v[284] = v[282] - v[283];
            /* 285 = \[OverBracket]_\[Yen]_163(v|W)_2|3;\[Xi] */
            v[285] = v[286] + v[287];
            /* 288 = \[OverBracket]_\[Yen]_157(v|W)_2|3;\[Xi] */
            v[288] = -v[286] + v[287];
            /* 289 = \[OverBracket]_\[Yen]_134(v|W)_2|2;\[Xi] */
            v[289] = (-1e0 / 720e0) * (v[162] * v[252]) - (v[110] * v[272]) / 720e0;
            /* 290 = \[OverBracket]_\[Yen]_152(v|W)_2|2;\[Xi] */
            v[290] = (-1e0 / 720e0) * (v[110] * v[252]);
            /* 291 = \[OverBracket]_\[Yen]_120(v|W)_1|3 */
            v[291] = v[292] + v[293];
            v[527] = v[103] * v[291];
            v[522] = v[101] * v[291];
            /* 294 = \[OverBracket]_\[Yen]_121(v|W)_1|3 */
            v[294] = v[292] - v[293];
            /* 295 = \[OverBracket]_\[Yen]_164(v|W)_1|3;\[Xi] */
            v[295] = v[296] + v[297];
            v[523] = v[102] * v[285] + v[101] * v[295];
            /* 298 = \[OverBracket]_\[Yen]_158(v|W)_1|3;\[Xi] */
            v[298] = -v[296] + v[297];
            /* 299 = \[OverBracket]_\[Yen]_114(v|W)_1|2 */
            v[299] = v[300] + v[301];
            v[528] = v[102] * v[299];
            v[525] = v[101] * v[299];
            /* 302 = \[OverBracket]_\[Yen]_115(v|W)_1|2 */
            v[302] = v[300] - v[301];
            /* 303 = \[OverBracket]_\[Yen]_165(v|W)_1|2;\[Xi] */
            v[303] = v[304] + v[305];
            v[529] = v[103] * v[295] + v[102] * v[303];
            v[526] = v[103] * v[285] + v[101] * v[303];
            /* 306 = \[OverBracket]_\[Yen]_159(v|W)_1|2;\[Xi] */
            v[306] = -v[304] + v[305];
            /* 307 = \[OverBracket]_\[Yen]_127(v|W)_1|1;\[Xi] */
            v[307] = (-1e0 / 720e0) * (v[162] * v[259]) - (v[110] * v[262]) / 720e0;
            /* 308 = \[OverBracket]_\[Yen]_153(v|W)_1|1;\[Xi] */
            v[308] = (-1e0 / 720e0) * (v[110] * v[259]);
            /* 309 = \[OverBracket]_\[Yen]_110(\[Yen]|W)_165 */
            v[309] = (-1e0 / 720e0) * (v[149] * v[242]) - (v[152] * v[252]) / 720e0 - (v[153] * v[259]) / 720e0 - (v[127] * v[262])
                / 720e0 - (v[134] * v[272]) / 720e0 - (v[139] * v[277]) / 720e0 + v[102] * v[103] * v[281] + v[101] * v[103] * v[291]
                + v[101] * v[102] * v[299] + v[303] * v[488] + v[295] * v[491] + v[285] * v[494];
            /* 310 = \[OverBracket]_\[Yen]_162(\[Yen]|W)_165 */
            v[310] = (-1e0 / 720e0) * (v[139] * v[242]) - (v[134] * v[252]) / 720e0 - (v[127] * v[259]) / 720e0
                + v[102] * v[103] * v[285] + v[101] * v[103] * v[295] + v[101] * v[102] * v[303];
            /* 311 = \[OverBracket]_u\[Phi]_4(\[Yen]|W)_157 */
            v[311] = -6e0 * v[112] * v[284] - 6e0 * v[156] * v[288] + v[103] * v[110] * v[291] + v[102] * v[110] * v[299]
                + v[295] * v[313] + v[303] * v[318];
            /* 312 = \[OverBracket]_u\[Phi]_4;\[Xi](\[Yen]|W)_157 */
            v[312] = -6e0 * v[112] * v[288] + v[103] * v[110] * v[295] + v[102] * v[110] * v[303];
            /* 314 = \[OverBracket]_u\[Phi]_5(\[Yen]|W)_158 */
            v[314] = v[103] * v[110] * v[281] + 6e0 * v[112] * v[294] + 6e0 * v[156] * v[298] + v[101] * v[110] * v[299]
                + v[285] * v[313] + v[303] * v[319];
            /* 315 = \[OverBracket]_u\[Phi]_5;\[Xi](\[Yen]|W)_158 */
            v[315] = v[103] * v[110] * v[285] + 6e0 * v[112] * v[298] + v[101] * v[110] * v[303];
            /* 316 = \[OverBracket]_\[Yen]_112(\[Yen]|W)_159 */
            v[316] = -6e0 * v[101] * v[284] - 6e0 * v[145] * v[288] + 6e0 * v[102] * v[294] + 6e0 * v[147] * v[298] - 6e0 * v[103] * v[302]
                - 6e0 * v[150] * v[306];
            /* 317 = \[OverBracket]_\[Yen]_156(\[Yen]|W)_159 */
            v[317] = -6e0 * v[101] * v[288] + 6e0 * v[102] * v[298] - 6e0 * v[103] * v[306];
            /* 320 = \[OverBracket]_u\[Phi]_6(\[Yen]|W)_159 */
            v[320] = v[102] * v[110] * v[281] + v[101] * v[110] * v[291] - 6e0 * v[112] * v[302] - 6e0 * v[156] * v[306]
                + v[285] * v[318] + v[295] * v[319];
            /* 321 = \[OverBracket]_u\[Phi]_6;\[Xi](\[Yen]|W)_159 */
            v[321] = v[102] * v[110] * v[285] + v[101] * v[110] * v[295] - 6e0 * v[112] * v[306];
            /* 322 = \[OverBracket]_\[Yen]_161(\[Yen]|W)_162 */
            v[322] = v[310] + v[317];
            /* 323 = \[OverBracket]_ff_(\[Yen]|W)_161 */
            v[323] = v[104] * v[309] + v[160] * v[309] + v[154] * v[310] + v[104] * v[316] + v[155] * v[316] + v[154] * v[317]
                + v[154] * v[322];
            /* 324 = \[OverBracket]_ff_;\[Xi](\[Yen]|W)_161 */
            v[324] = v[160] * v[310] + v[155] * v[317] + v[104] * v[322];
        }
        else {
            v[344] = -((v[129] * v[145]) / v[104]) - (v[101] * v[183]) / v[104] + v[101] * v[129] * v[154] * v[185];
            v[343] = -((v[129] * v[147]) / v[104]) - (v[102] * v[183]) / v[104] + v[102] * v[129] * v[154] * v[185];
            v[340] = -(v[130] * v[177]) - v[175] * v[179];
            v[339] = -((v[129] * v[150]) / v[104]) - (v[103] * v[183]) / v[104] + v[103] * v[129] * v[154] * v[185];
            v[329] = v[183] / v[104] - v[129] * v[154] * v[185];
            /* 279 = \[OverBracket]_\[Yen]_139(v|W)_3|3;\[Xi] */
            v[279] = (v[129] * v[277]) / v[104] + v[242] * v[329];
            /* 280 = \[OverBracket]_\[Yen]_149(v|W)_3|3;\[Xi] */
            v[280] = (v[129] * v[242]) / v[104];
            /* 325 = \[OverBracket]_\[Yen]_138(v|W)_2|3 */
            v[325] = v[274] + v[276];
            v[578] = v[102] * v[129] * v[185] * v[325];
            /* 326 = \[OverBracket]_\[Yen]_137(v|W)_2|3 */
            v[326] = v[274] - v[276];
            /* 327 = \[OverBracket]_\[Yen]_184(v|W)_2|3;\[Xi] */
            v[327] = v[243] + v[250];
            v[625] = v[103] * v[129] * v[185] * v[327];
            v[622] = v[102] * v[129] * v[185] * v[327];
            v[553] = (v[103] * v[147] * v[185] + v[102] * v[150] * v[185]) * v[327];
            /* 328 = \[OverBracket]_\[Yen]_180(v|W)_2|3;\[Xi] */
            v[328] = -v[243] + v[250];
            /* 289 = \[OverBracket]_\[Yen]_134(v|W)_2|2;\[Xi] */
            v[289] = (v[129] * v[272]) / v[104] + v[252] * v[329];
            /* 290 = \[OverBracket]_\[Yen]_152(v|W)_2|2;\[Xi] */
            v[290] = (v[129] * v[252]) / v[104];
            /* 330 = \[OverBracket]_\[Yen]_136(v|W)_1|3 */
            v[330] = v[268] + v[275];
            v[580] = v[101] * v[129] * v[185] * v[330];
            v[554] = v[101] * v[103] * v[185] * v[330];
            /* 331 = \[OverBracket]_\[Yen]_135(v|W)_1|3 */
            v[331] = v[268] - v[275];
            /* 332 = \[OverBracket]_\[Yen]_186(v|W)_1|3;\[Xi] */
            v[332] = v[244] + v[257];
            v[630] = v[103] * v[129] * v[185] * v[332];
            v[623] = v[101] * v[129] * v[185] * v[332];
            v[577] = v[622] + v[623];
            v[555] = (v[103] * v[145] * v[185] + v[101] * v[150] * v[185]) * v[332];
            /* 333 = \[OverBracket]_\[Yen]_181(v|W)_1|3;\[Xi] */
            v[333] = -v[244] + v[257];
            /* 334 = \[OverBracket]_\[Yen]_133(v|W)_1|2 */
            v[334] = v[265] + v[270];
            v[613] = v[102] * v[129] * v[185] * v[334];
            v[597] = v[101] * v[129] * v[185] * v[334];
            v[556] = v[101] * v[102] * v[185] * v[334];
            /* 335 = \[OverBracket]_\[Yen]_132(v|W)_1|2 */
            v[335] = v[265] - v[270];
            /* 336 = \[OverBracket]_\[Yen]_187(v|W)_1|2;\[Xi] */
            v[336] = v[254] + v[258];
            v[628] = -(v[129] * v[139] * v[242]) - v[129] * v[134] * v[252] - v[127] * v[129] * v[259]
                + v[102] * v[103] * v[129] * v[327] + v[101] * v[103] * v[129] * v[332] + v[101] * v[102] * v[129] * v[336];
            v[627] = -(v[139] * v[154] * v[242]) - v[134] * v[154] * v[252] - v[127] * v[154] * v[259]
                + v[102] * v[103] * v[154] * v[327] + v[101] * v[103] * v[154] * v[332] + v[101] * v[102] * v[154] * v[336];
            v[626] = -(v[139] * v[185] * v[242]) - v[134] * v[185] * v[252] - v[127] * v[185] * v[259]
                + v[102] * v[103] * v[185] * v[327] + v[101] * v[103] * v[185] * v[332] + v[101] * v[102] * v[185] * v[336];
            v[601] = v[102] * v[129] * v[185] * v[336];
            v[576] = v[101] * v[129] * v[185] * v[336];
            v[558] = -(v[129] * v[139] * v[154] * v[242]) - v[129] * v[134] * v[154] * v[252] - v[127] * v[129] * v[154] * v[259]
                + v[102] * v[103] * v[129] * v[154] * v[327] + v[101] * v[103] * v[129] * v[154] * v[332]
                + v[101] * v[102] * v[129] * v[154] * v[336];
            v[557] = (v[102] * v[145] * v[185] + v[101] * v[147] * v[185]) * v[336];
            /* 337 = \[OverBracket]_\[Yen]_182(v|W)_1|2;\[Xi] */
            v[337] = -v[254] + v[258];
            v[575] = -(v[101] * v[328]) + v[102] * v[333] - v[103] * v[337];
            v[549] = -(v[101] * v[326]) - v[145] * v[328] + v[102] * v[331] + v[147] * v[333] - v[103] * v[335] - v[150] * v[337];
            /* 307 = \[OverBracket]_\[Yen]_127(v|W)_1|1;\[Xi] */
            v[307] = (v[129] * v[262]) / v[104] + v[259] * v[329];
            /* 308 = \[OverBracket]_\[Yen]_153(v|W)_1|1;\[Xi] */
            v[308] = (v[129] * v[259]) / v[104];
            /* 338 = \[OverBracket]_\[Yen]_183(\[Yen]|W)_187 */
            v[338] = (v[139] * v[242]) / v[104] + (v[134] * v[252]) / v[104] + (v[127] * v[259]) / v[104] - (v[102] * v[103] * v[327]
                ) / v[104] - (v[101] * v[103] * v[332]) / v[104] - (v[101] * v[102] * v[336]) / v[104];
            v[541] = -(v[101] * v[175] * v[326]) + (-(v[145] * v[175]) - v[101] * v[177]) * v[328] + v[102] * v[175] * v[331] +
                (v[147] * v[175] + v[102] * v[177]) * v[333] - v[103] * v[175] * v[335] + (-(v[150] * v[175]) - v[103] * v[177]) * v[337]
                - v[176] * v[338];
            v[550] = v[178] * v[541];
            /* 311 = \[OverBracket]_u\[Phi]_4(\[Yen]|W)_180 */
            v[311] = -(v[130] * v[175] * v[326]) - (v[103] * v[129] * v[330]) / v[104] - (v[102] * v[129] * v[334]) / v[104]
                + v[332] * v[339] + v[328] * v[340] + v[336] * v[343];
            /* 312 = \[OverBracket]_u\[Phi]_4;\[Xi](\[Yen]|W)_180 */
            v[312] = -(v[130] * v[175] * v[328]) - (v[103] * v[129] * v[332]) / v[104] - (v[102] * v[129] * v[336]) / v[104];
            /* 314 = \[OverBracket]_u\[Phi]_5(\[Yen]|W)_181 */
            v[314] = -((v[103] * v[129] * v[325]) / v[104]) + v[130] * v[175] * v[331] - (v[101] * v[129] * v[334]) / v[104]
                + v[327] * v[339] - v[333] * v[340] + v[336] * v[344];
            /* 315 = \[OverBracket]_u\[Phi]_5;\[Xi](\[Yen]|W)_181 */
            v[315] = -((v[103] * v[129] * v[327]) / v[104]) + v[130] * v[175] * v[333] - (v[101] * v[129] * v[336]) / v[104];
            /* 341 = \[OverBracket]_\[Yen]_177(\[Yen]|W)_182 */
            v[341] = -(v[101] * v[130] * v[328]) + v[102] * v[130] * v[333] - v[103] * v[130] * v[337];
            /* 342 = \[OverBracket]_\[Yen]_179(\[Yen]|W)_182 */
            v[342] = -(v[101] * v[175] * v[328]) + v[102] * v[175] * v[333] - v[103] * v[175] * v[337];
            v[548] = v[176] * v[342];
            /* 320 = \[OverBracket]_u\[Phi]_6(\[Yen]|W)_182 */
            v[320] = -((v[102] * v[129] * v[325]) / v[104]) - (v[101] * v[129] * v[330]) / v[104] - v[130] * v[175] * v[335]
                + v[337] * v[340] + v[327] * v[343] + v[332] * v[344];
            /* 321 = \[OverBracket]_u\[Phi]_6;\[Xi](\[Yen]|W)_182 */
            v[321] = -((v[102] * v[129] * v[327]) / v[104]) - (v[101] * v[129] * v[332]) / v[104] - v[130] * v[175] * v[337];
            /* 346 = \[OverBracket]_\[Yen]_176(\[Yen]|W)_179 */
            v[346] = -(v[130] * v[338]) + v[178] * v[342] - v[341] * v[345];
            v[552] = -(v[101] * v[130] * v[326]) + (-(v[130] * v[145]) - v[101] * v[179]) * v[328] + v[102] * v[130] * v[331] +
                (v[130] * v[147] + v[102] * v[179]) * v[333] - v[103] * v[130] * v[335] + (-(v[130] * v[150]) - v[103] * v[179]) * v[337]
                + (v[154] * v[346]) / 2e0;
            /* 324 = \[OverBracket]_ff_;\[Xi](\[Yen]|W)_176 */
            v[324] = -(v[129] * v[139] * v[185] * v[242]) - v[129] * v[134] * v[185] * v[252] - v[127] * v[129] * v[185] * v[259]
                + v[102] * v[103] * v[129] * v[185] * v[327] + v[101] * v[103] * v[129] * v[185] * v[332]
                + v[101] * v[102] * v[129] * v[185] * v[336] + (v[175] * v[346]) / 2e0;
            /* 323 = \[OverBracket]_ff_(\[Yen]|W)_185 */
            v[323] = (-(v[129] * v[149] * v[185]) - v[139] * v[183] * v[185]) * v[242] + (-(v[129] * v[152] * v[185])
                - v[134] * v[183] * v[185]) * v[252] + (-(v[129] * v[153] * v[185]) - v[127] * v[183] * v[185]) * v[259]
                - v[127] * v[129] * v[185] * v[262] - v[129] * v[134] * v[185] * v[272] - v[129] * v[139] * v[185] * v[277]
                + v[102] * v[103] * v[129] * v[185] * v[325] + (v[103] * v[129] * v[147] * v[185] + v[102] * v[129] * v[150] * v[185]
                    + v[102] * v[103] * v[183] * v[185]) * v[327] + v[101] * v[103] * v[129] * v[185] * v[330] +
                (v[103] * v[129] * v[145] * v[185] + v[101] * v[129] * v[150] * v[185] + v[101] * v[103] * v[183] * v[185]) * v[332]
                + v[101] * v[102] * v[129] * v[185] * v[334] + (v[102] * v[129] * v[145] * v[185] + v[101] * v[129] * v[147] * v[185]
                    + v[101] * v[102] * v[183] * v[185]) * v[336] + (v[175] * ((2e0 * v[176] * v[341]) / Power(v[128], 3) - v[130] * ((v[149]
                        / v[104] - v[139] * v[154] * v[185]) * v[242] + (v[152] / v[104] - v[134] * v[154] * v[185]) * v[252] + (v[153] / v[104]
                            - v[127] * v[154] * v[185]) * v[259] + (v[127] * v[262]) / v[104] + (v[134] * v[272]) / v[104] + (v[139] * v[277])
                        / v[104] - (v[102] * v[103] * v[325]) / v[104] + (-((v[103] * v[147]) / v[104]) - (v[102] * v[150]) / v[104]
                            + v[102] * v[103] * v[154] * v[185]) * v[327] - (v[101] * v[103] * v[330]) / v[104] + (-((v[103] * v[145]) / v[104]) -
                                (v[101] * v[150]) / v[104] + v[101] * v[103] * v[154] * v[185]) * v[332] - (v[101] * v[102] * v[334]) / v[104] + (-(
                                    (v[102] * v[145]) / v[104]) - (v[101] * v[147]) / v[104] + v[101] * v[102] * v[154] * v[185]) * v[336] + v[548]) + v[550]
                        )) / 2e0 - (v[175] * v[552]) / (2e0 * v[104]) - (2e0 * v[558]) / Power(v[104], 3);
        };
        /* 347 = \[OverBracket]_\[Yen]_108(ff|W) */
        v[347] = v[323];
        /* 348 = \[OverBracket]_\[Yen]_107(ff|W) */
        v[348] = v[323];
        /* 349 = \[OverBracket]_\[Yen]_117(ff|W) */
        v[349] = v[323];
        /* 350 = \[OverBracket]_\[Yen]_108(\[Yen]|W)_127 */
        v[350] = v[307] + v[347];
        /* 351 = \[OverBracket]_\[Yen]_107(\[Yen]|W)_127 */
        v[351] = v[350];
        /* 352 = \[OverBracket]_\[Yen]_108(\[Yen]|W)_134 */
        v[352] = v[289] + v[350];
        /* 353 = \[OverBracket]_\[Yen]_117(\[Yen]|W)_134 */
        v[353] = v[289] + v[349];
        /* 354 = \[OverBracket]_\[Yen]_151(\[Yen]|W)_152 */
        v[354] = v[290];
        /* 355 = \[OverBracket]_\[Yen]_146(\[Yen]|W)_152 */
        v[355] = v[290];
        /* 356 = \[OverBracket]_\[Yen]_151(\[Yen]|W)_153 */
        v[356] = v[308] + v[354];
        /* 357 = \[OverBracket]_\[Yen]_148(\[Yen]|W)_153 */
        v[357] = v[308];
        /* 358 = \[OverBracket]_\[Yen]_151(ff|W)_;\[Xi] */
        v[358] = v[324] + v[356];
        /* 359 = \[OverBracket]_\[Yen]_148(ff|W)_;\[Xi] */
        v[359] = v[324] + v[357];
        /* 360 = \[OverBracket]_\[Yen]_146(ff|W)_;\[Xi] */
        v[360] = v[324] + v[355];
        /* 320 = \[OverBracket]_u\[Phi]_6(\[Yen]|W)_151 */
        v[320] = v[320] + 2e0 * v[103] * v[352] + 2e0 * v[150] * v[358];
        /* 321 = \[OverBracket]_u\[Phi]_6;\[Xi](\[Yen]|W)_151 */
        v[321] = v[321] + 2e0 * v[103] * v[358];
        /* 361 = \[OverBracket]_peIO_3|6(u\[Phi]|W)_6 */
        v[361] = v[320] * v[52];
        /* 362 = \[OverBracket]_peIO_2|6(u\[Phi]|W)_6 */
        v[362] = v[320] * v[50];
        /* 363 = \[OverBracket]_peIO_1|6(u\[Phi]|W)_6 */
        v[363] = v[320] * v[48];
        /* 364 = \[OverBracket]_\[Yen]_107(\[Yen]|W)_139 */
        v[364] = v[279] + v[351];
        /* 365 = \[OverBracket]_\[Yen]_117(\[Yen]|W)_139 */
        v[365] = v[279] + v[353];
        /* 366 = \[OverBracket]_\[Yen]_148(\[Yen]|W)_149 */
        v[366] = v[280] + v[359];
        /* 367 = \[OverBracket]_\[Yen]_146(\[Yen]|W)_149 */
        v[367] = v[280] + v[360];
        /* 314 = \[OverBracket]_u\[Phi]_5(\[Yen]|W)_148 */
        v[314] = v[314] + 2e0 * v[102] * v[364] + 2e0 * v[147] * v[366];
        /* 315 = \[OverBracket]_u\[Phi]_5;\[Xi](\[Yen]|W)_148 */
        v[315] = v[315] + 2e0 * v[102] * v[366];
        /* 368 = \[OverBracket]_peIO_3|5(u\[Phi]|W)_5 */
        v[368] = v[314] * v[52];
        /* 369 = \[OverBracket]_peIO_2|5(u\[Phi]|W)_5 */
        v[369] = v[314] * v[50];
        /* 370 = \[OverBracket]_peIO_1|5(u\[Phi]|W)_5 */
        v[370] = v[314] * v[48];
        /* 311 = \[OverBracket]_u\[Phi]_4(\[Yen]|W)_146 */
        v[311] = v[311] + 2e0 * v[101] * v[365] + 2e0 * v[145] * v[367];
        /* 312 = \[OverBracket]_u\[Phi]_4;\[Xi](\[Yen]|W)_146 */
        v[312] = v[312] + 2e0 * v[101] * v[367];
        /* 371 = \[OverBracket]_peIO_3|4(u\[Phi]|W)_4 */
        v[371] = v[311] * v[52];
        /* 372 = \[OverBracket]_peIO_2|4(u\[Phi]|W)_4 */
        v[372] = v[311] * v[50];
        /* 373 = \[OverBracket]_peIO_1|4(u\[Phi]|W)_4 */
        v[373] = v[311] * v[48];
        /* 374 = \[OverBracket]_peIO_3|4(u\[Phi]|W)_4;\[Xi] */
        v[374] = v[371] + v[312] * v[59];
        /* 375 = \[OverBracket]_peIO_2|4(u\[Phi]|W)_4;\[Xi] */
        v[375] = v[372] + v[312] * v[57];
        /* 376 = \[OverBracket]_peIO_1|4(u\[Phi]|W)_4;\[Xi] */
        v[376] = v[373] + v[312] * v[56];
        /* 377 = \[OverBracket]_peIO_3|5(u\[Phi]|W)_5;\[Xi] */
        v[377] = v[368] + v[315] * v[59];
        /* 378 = \[OverBracket]_peIO_2|5(u\[Phi]|W)_5;\[Xi] */
        v[378] = v[369] + v[315] * v[57];
        /* 379 = \[OverBracket]_peIO_1|5(u\[Phi]|W)_5;\[Xi] */
        v[379] = v[370] + v[315] * v[56];
        /* 380 = \[OverBracket]_peIO_3|6(u\[Phi]|W)_6;\[Xi] */
        v[380] = v[361] + v[321] * v[59];
        /* 381 = \[OverBracket]_peIO_2|6(u\[Phi]|W)_6;\[Xi] */
        v[381] = v[362] + v[321] * v[57];
        /* 382 = \[OverBracket]_peIO_1|6(u\[Phi]|W)_6;\[Xi] */
        v[382] = v[363] + v[321] * v[56];
        v[1167] = v[238] * v[56];
        v[1168] = v[237] * v[56];
        v[1169] = v[235] * v[56];
        v[1170] = v[376];
        v[1171] = v[379];
        v[1172] = v[382];
        v[1173] = v[238] * v[57];
        v[1174] = v[237] * v[57];
        v[1175] = v[235] * v[57];
        v[1176] = v[375];
        v[1177] = v[378];
        v[1178] = v[381];
        v[1179] = v[238] * v[59];
        v[1180] = v[237] * v[59];
        v[1181] = v[235] * v[59];
        v[1182] = v[374];
        v[1183] = v[377];
        v[1184] = v[380];
        v[383] = 0e0;/*debug*/
        for (i227 = 1; i227 <= 18; i227++) {
            v[227] = i227;
            /* 228 = \[DoubleStruckCapitalG]_m */
            v[228] = v[1148 + i227];
            /* 384 = Rgm */
            v[384] = v[1166 + i227] * v[97];
            /* 389 = \[OverBracket]_\[OverBracket]_peIO_1|6(u\[Phi]|W)(Rgm|Rgm)_6;\[Xi] */
            v[389] = (6 == i227 ? 1 : 0) * v[97];
            /* 390 = \[OverBracket]_\[OverBracket]_peIO_2|6(u\[Phi]|W)(Rgm|Rgm)_6;\[Xi] */
            v[390] = (12 == i227 ? 1 : 0) * v[97];
            /* 391 = \[OverBracket]_\[OverBracket]_peIO_3|6(u\[Phi]|W)(Rgm|Rgm)_6;\[Xi] */
            v[391] = (18 == i227 ? 1 : 0) * v[97];
            v[416] = v[389] * v[48] + v[390] * v[50] + v[391] * v[52];
            /* 392 = \[OverBracket]_\[OverBracket]_peIO_1|5(u\[Phi]|W)(Rgm|Rgm)_5;\[Xi] */
            v[392] = (5 == i227 ? 1 : 0) * v[97];
            /* 393 = \[OverBracket]_\[OverBracket]_peIO_2|5(u\[Phi]|W)(Rgm|Rgm)_5;\[Xi] */
            v[393] = (11 == i227 ? 1 : 0) * v[97];
            /* 394 = \[OverBracket]_\[OverBracket]_peIO_3|5(u\[Phi]|W)(Rgm|Rgm)_5;\[Xi] */
            v[394] = (17 == i227 ? 1 : 0) * v[97];
            v[407] = v[392] * v[48] + v[393] * v[50] + v[394] * v[52];
            /* 395 = \[OverBracket]_\[OverBracket]_peIO_1|4(u\[Phi]|W)(Rgm|Rgm)_4;\[Xi] */
            v[395] = (4 == i227 ? 1 : 0) * v[97];
            /* 396 = \[OverBracket]_\[OverBracket]_peIO_2|4(u\[Phi]|W)(Rgm|Rgm)_4;\[Xi] */
            v[396] = (10 == i227 ? 1 : 0) * v[97];
            /* 397 = \[OverBracket]_\[OverBracket]_peIO_3|4(u\[Phi]|W)(Rgm|Rgm)_4;\[Xi] */
            v[397] = (16 == i227 ? 1 : 0) * v[97];
            v[400] = v[395] * v[48] + v[396] * v[50] + v[397] * v[52];
            /* 398 = \[OverBracket]_\[Yen]_312|3(\[Yen]|Rgm)_312|3 */
            v[398] = v[395] * v[56] + v[396] * v[57] + v[397] * v[59];
            /* 399 = [\[OverBracket]_\[OverBracket]_\[Yen]_146(\[Yen]|W)\[OverBracket]_149(u\[Phi]|Rgm)_4(\[Yen]|W)_146][\[OverBracket]_\[OverBracket]_\[Yen]_146(ff|W)\[OverBracket]_;\[Xi](\[Yen]|Rgm)_146(\[Yen]|W)_149] */
            v[399] = 2e0 * v[101] * v[398] + 2e0 * v[145] * v[400];
            /* 401 = [\[OverBracket]_\[OverBracket]_\[Yen]_117(\[Yen]|W)\[OverBracket]_139(u\[Phi]|Rgm)_4(\[Yen]|W)_146][\[OverBracket]_\[OverBracket]_\[Yen]_117(\[Yen]|W)\[OverBracket]_134(\[Yen]|Rgm)_117(\[Yen]|W)_139] */
            v[401] = 2e0 * v[101] * v[400];
            /* 402 = \[OverBracket]_u\[Phi]_4\[OverBracket]_(u\[Phi]|Rgm)_4(\[Yen]|W)_146 */
            v[402] = 2e0 * v[367] * v[398] + 2e0 * v[365] * v[400];
            v[532] = v[402];
            /* 403 = \[OverBracket]_u\[Phi]_4;\[Xi]\[OverBracket]_(u\[Phi]|Rgm)_4(\[Yen]|W)_146 */
            v[403] = 2e0 * v[367] * v[400];
            v[533] = v[403];
            /* 404 = \[OverBracket]_\[Yen]_311|3(\[Yen]|Rgm)_311|3 */
            v[404] = v[400];
            /* 405 = \[OverBracket]_\[Yen]_315|3(\[Yen]|Rgm)_315|3 */
            v[405] = v[392] * v[56] + v[393] * v[57] + v[394] * v[59];
            /* 406 = [\[OverBracket]_\[OverBracket]_\[Yen]_148(\[Yen]|W)\[OverBracket]_149(u\[Phi]|Rgm)_5(\[Yen]|W)_148][\[OverBracket]_\[OverBracket]_\[Yen]_148(ff|W)\[OverBracket]_;\[Xi](\[Yen]|Rgm)_148(\[Yen]|W)_149] */
            v[406] = 2e0 * v[102] * v[405] + 2e0 * v[147] * v[407];
            /* 408 = \[OverBracket]_\[OverBracket]_\[Yen]_107(\[Yen]|W)\[OverBracket]_139(u\[Phi]|Rgm)_5(\[Yen]|W)_148 */
            v[408] = 2e0 * v[102] * v[407];
            /* 409 = \[OverBracket]_u\[Phi]_5\[OverBracket]_(u\[Phi]|Rgm)_5(\[Yen]|W)_148 */
            v[409] = 2e0 * v[366] * v[405] + 2e0 * v[364] * v[407];
            v[531] = v[409];
            /* 410 = \[OverBracket]_u\[Phi]_5;\[Xi]\[OverBracket]_(u\[Phi]|Rgm)_5(\[Yen]|W)_148 */
            v[410] = 2e0 * v[366] * v[407];
            v[534] = v[410];
            /* 411 = \[OverBracket]_\[Yen]_314|3(\[Yen]|Rgm)_314|3 */
            v[411] = v[407];
            /* 412 = \[OverBracket]_\[Yen]_280|3\[OverBracket]_(\[Yen]|Rgm)_148(\[Yen]|W)_149 */
            v[412] = v[399] + v[406];
            /* 413 = \[OverBracket]_\[Yen]_279|3\[OverBracket]_(\[Yen]|Rgm)_107(\[Yen]|W)_139 */
            v[413] = v[401] + v[408];
            /* 414 = \[OverBracket]_\[Yen]_321|3(\[Yen]|Rgm)_321|3 */
            v[414] = v[389] * v[56] + v[390] * v[57] + v[391] * v[59];
            /* 415 = [\[OverBracket]_\[OverBracket]_\[Yen]_151(ff|W)\[OverBracket]_;\[Xi](u\[Phi]|Rgm)_6(\[Yen]|W)_151][\[OverBracket]_\[OverBracket]_\[Yen]_151(\[Yen]|W)\[OverBracket]_153(\[Yen]|Rgm)_151(ff|W)_;\[Xi]] */
            v[415] = 2e0 * v[103] * v[414] + 2e0 * v[150] * v[416];
            /* 417 = \[OverBracket]_\[OverBracket]_\[Yen]_108(\[Yen]|W)\[OverBracket]_134(u\[Phi]|Rgm)_6(\[Yen]|W)_151 */
            v[417] = 2e0 * v[103] * v[416];
            /* 418 = \[OverBracket]_u\[Phi]_6\[OverBracket]_(u\[Phi]|Rgm)_6(\[Yen]|W)_151 */
            v[418] = 2e0 * v[358] * v[414] + 2e0 * v[352] * v[416];
            v[530] = v[418];
            /* 419 = \[OverBracket]_u\[Phi]_6;\[Xi]\[OverBracket]_(u\[Phi]|Rgm)_6(\[Yen]|W)_151 */
            v[419] = 2e0 * v[358] * v[416];
            v[535] = v[419];
            /* 420 = \[OverBracket]_\[Yen]_320|3(\[Yen]|Rgm)_320|3 */
            v[420] = v[416];
            /* 421 = \[OverBracket]_\[Yen]_324|3\[OverBracket]_(\[Yen]|Rgm)_151(ff|W)_;\[Xi] */
            v[421] = v[399] + v[406] + v[415];
            /* 422 = \[OverBracket]_\[Yen]_308|3\[OverBracket]_(\[Yen]|Rgm)_151(\[Yen]|W)_153 */
            v[422] = v[406] + v[415];
            /* 423 = \[OverBracket]_\[Yen]_290|3\[OverBracket]_(\[Yen]|Rgm)_151(\[Yen]|W)_152 */
            v[423] = v[399] + v[415];
            /* 424 = [\[OverBracket]_\[OverBracket]_\[Yen]_108(\[Yen]|W)\[OverBracket]_127(\[Yen]|Rgm)_108(\[Yen]|W)_134][\[OverBracket]_\[Yen]_307|3\[OverBracket]_(\[Yen]|Rgm)_108(\[Yen]|W)_127] */
            v[424] = v[408] + v[417];
            /* 425 = \[OverBracket]_\[Yen]_289|3\[OverBracket]_(\[Yen]|Rgm)_108(\[Yen]|W)_134 */
            v[425] = v[401] + v[417];
            /* 426 = \[OverBracket]_\[Yen]_323|3\[OverBracket]_(\[Yen]|Rgm)_108(ff|W) */
            v[426] = v[401] + v[424];
            /* 427 = \[OverBracket]_\[Yen]_110(\[Yen]|Rgm)_110 */
            v[427] = 0e0;
            /* 428 = \[OverBracket]_\[Yen]_112(\[Yen]|Rgm)_112 */
            v[428] = 0e0;
            /* 429 = \[OverBracket]_\[Yen]_128(\[Yen]|Rgm)_128 */
            v[429] = 0e0;
            /* 430 = \[OverBracket]_\[Yen]_129(\[Yen]|Rgm)_129 */
            v[430] = 0e0;
            /* 431 = \[OverBracket]_\[Yen]_130(\[Yen]|Rgm)_130 */
            v[431] = 0e0;
            /* 432 = \[OverBracket]_\[Yen]_155(\[Yen]|Rgm)_155 */
            v[432] = 0e0;
            /* 433 = \[OverBracket]_\[Yen]_156(\[Yen]|Rgm)_156 */
            v[433] = 0e0;
            /* 434 = \[OverBracket]_\[Yen]_160(\[Yen]|Rgm)_160 */
            v[434] = 0e0;
            /* 435 = \[OverBracket]_\[Yen]_162(\[Yen]|Rgm)_162 */
            v[435] = 0e0;
            /* 436 = \[OverBracket]_\[Yen]_175(\[Yen]|Rgm)_175 */
            v[436] = 0e0;
            /* 437 = \[OverBracket]_\[Yen]_176(\[Yen]|Rgm)_176 */
            v[437] = 0e0;
            /* 438 = \[OverBracket]_\[Yen]_177(\[Yen]|Rgm)_177 */
            v[438] = 0e0;
            /* 439 = \[OverBracket]_\[Yen]_178(\[Yen]|Rgm)_178 */
            v[439] = 0e0;
            /* 440 = \[OverBracket]_\[Yen]_179(\[Yen]|Rgm)_179 */
            v[440] = 0e0;
            /* 441 = \[OverBracket]_\[Yen]_183(\[Yen]|Rgm)_183 */
            v[441] = 0e0;
            /* 442 = \[OverBracket]_\[Yen]_185(\[Yen]|Rgm)_185 */
            v[442] = 0e0;
            /* 443 = \[OverBracket]_\[Yen]_345(\[Yen]|Rgm)_345 */
            v[443] = 0e0;
            b444 = b105;
            v[444] = b444;
            if (b444) {
                v[450] = v[154] * v[426];
                /* 445 = \[OverBracket]_\[OverBracket]_\[Yen]_161(\[Yen]|W)\[OverBracket]_162(ff|Rgm)_(\[Yen]|W)_161 */
                v[445] = v[104] * v[421] + v[450];
                /* 446 = \[OverBracket]_\[OverBracket]_\[Yen]_112(\[Yen]|W)\[OverBracket]_159(ff|Rgm)_(\[Yen]|W)_161 */
                v[446] = (v[104] + v[155]) * v[426];
                /* 447 = \[OverBracket]_\[OverBracket]_\[Yen]_110(\[Yen]|W)\[OverBracket]_165(ff|Rgm)_(\[Yen]|W)_161 */
                v[447] = (v[104] + v[160]) * v[426];
                /* 432 = \[OverBracket]_\[Yen]_155\[OverBracket]_(ff|Rgm)_(\[Yen]|W)_161 */
                v[432] = v[317] * v[421] + v[316] * v[426];
                /* 434 = \[OverBracket]_\[Yen]_160\[OverBracket]_(ff|Rgm)_(\[Yen]|W)_161 */
                v[434] = v[310] * v[421] + v[309] * v[426];
                /* 448 = \[OverBracket]_ff_\[OverBracket]_(ff|Rgm)_(\[Yen]|W)_161 */
                v[448] = v[322] * v[421] + (v[309] + v[316]) * v[426];
                /* 449 = \[OverBracket]_ff_;\[Xi]\[OverBracket]_(ff|Rgm)_(\[Yen]|W)_161 */
                v[449] = (v[310] + v[317] + v[322]) * v[426];
                /* 451 = \[OverBracket]_\[OverBracket]_\[Yen]_156(\[Yen]|W)\[OverBracket]_159(\[Yen]|Rgm)_161(\[Yen]|W)_162 */
                v[451] = v[155] * v[421] + v[445] + v[450];
                /* 452 = \[OverBracket]_\[OverBracket]_\[Yen]_162(\[Yen]|W)\[OverBracket]_165(\[Yen]|Rgm)_161(\[Yen]|W)_162 */
                v[452] = v[160] * v[421] + v[445] + v[450];
                /* 453 = \[OverBracket]_\[OverBracket]_\[Yen]_159(v|W)\[OverBracket]_1|2;\[Xi](u\[Phi]|Rgm)_6;\[Xi](\[Yen]|W)_159 */
                v[453] = -6e0 * v[112] * v[414];
                /* 454 = \[OverBracket]_\[OverBracket]_\[Yen]_164(v|W)\[OverBracket]_1|3;\[Xi](u\[Phi]|Rgm)_6;\[Xi](\[Yen]|W)_159 */
                v[454] = v[101] * v[110] * v[414];
                /* 455 = \[OverBracket]_\[OverBracket]_\[Yen]_163(v|W)\[OverBracket]_2|3;\[Xi](u\[Phi]|Rgm)_6;\[Xi](\[Yen]|W)_159 */
                v[455] = v[102] * v[110] * v[414];
                /* 427 = \[OverBracket]_\[Yen]_110\[OverBracket]_(u\[Phi]|Rgm)_6;\[Xi](\[Yen]|W)_159 */
                v[427] = v[414] * v[523];
                /* 428 = \[OverBracket]_\[Yen]_112\[OverBracket]_(u\[Phi]|Rgm)_6;\[Xi](\[Yen]|W)_159 */
                v[428] = -6e0 * v[306] * v[414];
                /* 409 = \[OverBracket]_u\[Phi]_5\[OverBracket]_(u\[Phi]|Rgm)_6;\[Xi](\[Yen]|W)_159 */
                v[409] = v[409] + v[110] * v[285] * v[414];
                /* 402 = \[OverBracket]_u\[Phi]_4\[OverBracket]_(u\[Phi]|Rgm)_6;\[Xi](\[Yen]|W)_159 */
                v[402] = v[402] + v[110] * v[295] * v[414];
                /* 456 = \[OverBracket]_\[OverBracket]_\[Yen]_159(v|W)\[OverBracket]_1|2;\[Xi](u\[Phi]|Rgm)_6(\[Yen]|W)_159 */
                v[456] = -6e0 * v[156] * v[420] + v[453];
                /* 457 = \[OverBracket]_\[OverBracket]_\[Yen]_115(v|W)\[OverBracket]_1|2(u\[Phi]|Rgm)_6(\[Yen]|W)_159 */
                v[457] = -6e0 * v[112] * v[420];
                /* 458 = \[OverBracket]_\[OverBracket]_\[Yen]_164(v|W)\[OverBracket]_1|3;\[Xi](u\[Phi]|Rgm)_6(\[Yen]|W)_159 */
                v[458] = v[319] * v[420] + v[454];
                /* 459 = \[OverBracket]_\[OverBracket]_\[Yen]_120(v|W)\[OverBracket]_1|3(u\[Phi]|Rgm)_6(\[Yen]|W)_159 */
                v[459] = v[101] * v[110] * v[420];
                /* 460 = \[OverBracket]_\[OverBracket]_\[Yen]_163(v|W)\[OverBracket]_2|3;\[Xi](u\[Phi]|Rgm)_6(\[Yen]|W)_159 */
                v[460] = v[318] * v[420] + v[455];
                /* 461 = \[OverBracket]_\[OverBracket]_\[Yen]_123(v|W)\[OverBracket]_2|3(u\[Phi]|Rgm)_6(\[Yen]|W)_159 */
                v[461] = v[102] * v[110] * v[420];
                /* 462 = \[OverBracket]_\[Yen]_318\[OverBracket]_(u\[Phi]|Rgm)_6(\[Yen]|W)_159 */
                v[462] = v[285] * v[420];
                /* 463 = \[OverBracket]_\[Yen]_319\[OverBracket]_(u\[Phi]|Rgm)_6(\[Yen]|W)_159 */
                v[463] = v[295] * v[420];
                /* 427 = \[OverBracket]_\[Yen]_110\[OverBracket]_(u\[Phi]|Rgm)_6(\[Yen]|W)_159 */
                v[427] = v[427] + v[420] * (v[521] + v[522]);
                /* 428 = \[OverBracket]_\[Yen]_112\[OverBracket]_(u\[Phi]|Rgm)_6(\[Yen]|W)_159 */
                v[428] = -6e0 * v[302] * v[420] + v[428];
                /* 433 = \[OverBracket]_\[Yen]_156\[OverBracket]_(u\[Phi]|Rgm)_6(\[Yen]|W)_159 */
                v[433] = -6e0 * v[306] * v[420];
                /* 409 = \[OverBracket]_u\[Phi]_5\[OverBracket]_(u\[Phi]|Rgm)_6(\[Yen]|W)_159 */
                v[409] = v[409] + v[110] * v[281] * v[420];
                /* 402 = \[OverBracket]_u\[Phi]_4\[OverBracket]_(u\[Phi]|Rgm)_6(\[Yen]|W)_159 */
                v[402] = v[402] + v[110] * v[291] * v[420];
                /* 464 = \[OverBracket]_\[OverBracket]_\[Yen]_159(v|W)\[OverBracket]_1|2;\[Xi](\[Yen]|Rgm)_112(\[Yen]|W)_159 */
                v[464] = -6e0 * v[150] * v[446] - 6e0 * v[103] * v[451] + v[456];
                /* 465 = \[OverBracket]_\[OverBracket]_\[Yen]_115(v|W)\[OverBracket]_1|2(\[Yen]|Rgm)_112(\[Yen]|W)_159 */
                v[465] = -6e0 * v[103] * v[446] + v[457];
                /* 466 = \[OverBracket]_\[OverBracket]_\[Yen]_165(v|W)\[OverBracket]_1|2;\[Xi](u\[Phi]|Rgm)_5;\[Xi](\[Yen]|W)_158 */
                v[466] = v[101] * v[110] * v[405];
                /* 467 = \[OverBracket]_\[OverBracket]_\[Yen]_158(v|W)\[OverBracket]_1|3;\[Xi](u\[Phi]|Rgm)_5;\[Xi](\[Yen]|W)_158 */
                v[467] = 6e0 * v[112] * v[405] + 6e0 * v[147] * v[446] + 6e0 * v[102] * v[451];
                /* 468 = \[OverBracket]_\[OverBracket]_\[Yen]_163(v|W)\[OverBracket]_2|3;\[Xi](u\[Phi]|Rgm)_5;\[Xi](\[Yen]|W)_158 */
                v[468] = v[103] * v[110] * v[405] + v[460];
                /* 427 = \[OverBracket]_\[Yen]_110\[OverBracket]_(u\[Phi]|Rgm)_5;\[Xi](\[Yen]|W)_158 */
                v[427] = v[427] + v[405] * v[526];
                /* 428 = \[OverBracket]_\[Yen]_112\[OverBracket]_(u\[Phi]|Rgm)_5;\[Xi](\[Yen]|W)_158 */
                v[428] = 6e0 * v[298] * v[405] + v[428];
                /* 418 = \[OverBracket]_u\[Phi]_6\[OverBracket]_(u\[Phi]|Rgm)_5;\[Xi](\[Yen]|W)_158 */
                v[418] = v[110] * v[285] * v[405] + v[418] - 6e0 * v[302] * v[446] - 6e0 * v[306] * v[451];
                /* 402 = \[OverBracket]_u\[Phi]_4\[OverBracket]_(u\[Phi]|Rgm)_5;\[Xi](\[Yen]|W)_158 */
                v[402] = v[402] + v[110] * v[303] * v[405] - 6e0 * v[284] * v[446] - 6e0 * v[288] * v[451];
                /* 469 = \[OverBracket]_\[OverBracket]_\[Yen]_165(v|W)\[OverBracket]_1|2;\[Xi](u\[Phi]|Rgm)_5(\[Yen]|W)_158 */
                v[469] = v[319] * v[411] + v[466];
                /* 470 = \[OverBracket]_\[OverBracket]_\[Yen]_114(v|W)\[OverBracket]_1|2(u\[Phi]|Rgm)_5(\[Yen]|W)_158 */
                v[470] = v[101] * v[110] * v[411];
                /* 471 = \[OverBracket]_\[OverBracket]_\[Yen]_158(v|W)\[OverBracket]_1|3;\[Xi](u\[Phi]|Rgm)_5(\[Yen]|W)_158 */
                v[471] = 6e0 * v[156] * v[411] + v[467];
                /* 472 = \[OverBracket]_\[OverBracket]_\[Yen]_121(v|W)\[OverBracket]_1|3(u\[Phi]|Rgm)_5(\[Yen]|W)_158 */
                v[472] = 6e0 * v[112] * v[411] + 6e0 * v[102] * v[446];
                /* 473 = \[OverBracket]_\[OverBracket]_\[Yen]_163(v|W)\[OverBracket]_2|3;\[Xi](u\[Phi]|Rgm)_5(\[Yen]|W)_158 */
                v[473] = v[313] * v[411] + v[468];
                /* 474 = \[OverBracket]_\[OverBracket]_\[Yen]_123(v|W)\[OverBracket]_2|3(u\[Phi]|Rgm)_5(\[Yen]|W)_158 */
                v[474] = v[103] * v[110] * v[411] + v[461];
                /* 475 = \[OverBracket]_\[Yen]_313\[OverBracket]_(u\[Phi]|Rgm)_5(\[Yen]|W)_158 */
                v[475] = v[285] * v[411];
                /* 476 = \[OverBracket]_\[Yen]_319\[OverBracket]_(u\[Phi]|Rgm)_5(\[Yen]|W)_158 */
                v[476] = v[303] * v[411] + v[463];
                /* 427 = \[OverBracket]_\[Yen]_110\[OverBracket]_(u\[Phi]|Rgm)_5(\[Yen]|W)_158 */
                v[427] = v[427] + v[411] * (v[524] + v[525]);
                /* 428 = \[OverBracket]_\[Yen]_112\[OverBracket]_(u\[Phi]|Rgm)_5(\[Yen]|W)_158 */
                v[428] = 6e0 * v[294] * v[411] + v[428];
                /* 433 = \[OverBracket]_\[Yen]_156\[OverBracket]_(u\[Phi]|Rgm)_5(\[Yen]|W)_158 */
                v[433] = 6e0 * v[298] * v[411] + v[433];
                /* 418 = \[OverBracket]_u\[Phi]_6\[OverBracket]_(u\[Phi]|Rgm)_5(\[Yen]|W)_158 */
                v[418] = v[110] * v[281] * v[411] + v[418];
                /* 402 = \[OverBracket]_u\[Phi]_4\[OverBracket]_(u\[Phi]|Rgm)_5(\[Yen]|W)_158 */
                v[402] = v[402] + v[110] * v[299] * v[411];
                /* 477 = \[OverBracket]_\[OverBracket]_\[Yen]_165(v|W)\[OverBracket]_1|2;\[Xi](u\[Phi]|Rgm)_4;\[Xi](\[Yen]|W)_157 */
                v[477] = v[102] * v[110] * v[398] + v[469];
                /* 478 = \[OverBracket]_\[OverBracket]_\[Yen]_164(v|W)\[OverBracket]_1|3;\[Xi](u\[Phi]|Rgm)_4;\[Xi](\[Yen]|W)_157 */
                v[478] = v[103] * v[110] * v[398] + v[458];
                /* 479 = \[OverBracket]_\[OverBracket]_\[Yen]_157(v|W)\[OverBracket]_2|3;\[Xi](u\[Phi]|Rgm)_4;\[Xi](\[Yen]|W)_157 */
                v[479] = -6e0 * v[112] * v[398] - 6e0 * v[145] * v[446] - 6e0 * v[101] * v[451];
                /* 427 = \[OverBracket]_\[Yen]_110\[OverBracket]_(u\[Phi]|Rgm)_4;\[Xi](\[Yen]|W)_157 */
                v[427] = v[427] + v[398] * v[529];
                /* 428 = \[OverBracket]_\[Yen]_112\[OverBracket]_(u\[Phi]|Rgm)_4;\[Xi](\[Yen]|W)_157 */
                v[428] = -6e0 * v[288] * v[398] + v[428];
                /* 418 = \[OverBracket]_u\[Phi]_6\[OverBracket]_(u\[Phi]|Rgm)_4;\[Xi](\[Yen]|W)_157 */
                v[418] = v[110] * v[295] * v[398] + v[418];
                /* 409 = \[OverBracket]_u\[Phi]_5\[OverBracket]_(u\[Phi]|Rgm)_4;\[Xi](\[Yen]|W)_157 */
                v[409] = v[110] * v[303] * v[398] + v[409] + 6e0 * v[294] * v[446] + 6e0 * v[298] * v[451];
                /* 480 = \[OverBracket]_\[OverBracket]_\[Yen]_165(v|W)\[OverBracket]_1|2;\[Xi](u\[Phi]|Rgm)_4(\[Yen]|W)_157 */
                v[480] = v[318] * v[404] + v[477];
                /* 481 = \[OverBracket]_\[OverBracket]_\[Yen]_114(v|W)\[OverBracket]_1|2(u\[Phi]|Rgm)_4(\[Yen]|W)_157 */
                v[481] = v[102] * v[110] * v[404] + v[470];
                /* 482 = \[OverBracket]_\[OverBracket]_\[Yen]_164(v|W)\[OverBracket]_1|3;\[Xi](u\[Phi]|Rgm)_4(\[Yen]|W)_157 */
                v[482] = v[313] * v[404] + v[478];
                /* 483 = \[OverBracket]_\[OverBracket]_\[Yen]_120(v|W)\[OverBracket]_1|3(u\[Phi]|Rgm)_4(\[Yen]|W)_157 */
                v[483] = v[103] * v[110] * v[404] + v[459];
                /* 484 = \[OverBracket]_\[OverBracket]_\[Yen]_157(v|W)\[OverBracket]_2|3;\[Xi](u\[Phi]|Rgm)_4(\[Yen]|W)_157 */
                v[484] = -6e0 * v[156] * v[404] + v[479];
                /* 485 = \[OverBracket]_\[OverBracket]_\[Yen]_124(v|W)\[OverBracket]_2|3(u\[Phi]|Rgm)_4(\[Yen]|W)_157 */
                v[485] = -6e0 * v[112] * v[404] - 6e0 * v[101] * v[446];
                /* 486 = \[OverBracket]_\[Yen]_313\[OverBracket]_(u\[Phi]|Rgm)_4(\[Yen]|W)_157 */
                v[486] = v[295] * v[404] + v[475];
                /* 487 = \[OverBracket]_\[Yen]_318\[OverBracket]_(u\[Phi]|Rgm)_4(\[Yen]|W)_157 */
                v[487] = v[303] * v[404] + v[462];
                /* 427 = \[OverBracket]_\[Yen]_110\[OverBracket]_(u\[Phi]|Rgm)_4(\[Yen]|W)_157 */
                v[427] = v[427] + v[404] * (v[527] + v[528]);
                /* 428 = \[OverBracket]_\[Yen]_112\[OverBracket]_(u\[Phi]|Rgm)_4(\[Yen]|W)_157 */
                v[428] = -6e0 * v[284] * v[404] + v[428];
                /* 433 = \[OverBracket]_\[Yen]_156\[OverBracket]_(u\[Phi]|Rgm)_4(\[Yen]|W)_157 */
                v[433] = -6e0 * v[288] * v[404] + v[433];
                /* 418 = \[OverBracket]_u\[Phi]_6\[OverBracket]_(u\[Phi]|Rgm)_4(\[Yen]|W)_157 */
                v[418] = v[110] * v[291] * v[404] + v[418];
                /* 409 = \[OverBracket]_u\[Phi]_5\[OverBracket]_(u\[Phi]|Rgm)_4(\[Yen]|W)_157 */
                v[409] = v[110] * v[299] * v[404] + v[409];
                /* 489 = \[OverBracket]_\[OverBracket]_\[Yen]_165(v|W)\[OverBracket]_1|2;\[Xi](\[Yen]|Rgm)_110(\[Yen]|W)_165 */
                v[489] = v[101] * v[102] * v[452] + v[480] + v[447] * v[488];
                /* 490 = \[OverBracket]_\[OverBracket]_\[Yen]_114(v|W)\[OverBracket]_1|2(\[Yen]|Rgm)_110(\[Yen]|W)_165 */
                v[490] = v[101] * v[102] * v[447] + v[481];
                /* 492 = \[OverBracket]_\[OverBracket]_\[Yen]_164(v|W)\[OverBracket]_1|3;\[Xi](\[Yen]|Rgm)_110(\[Yen]|W)_165 */
                v[492] = v[101] * v[103] * v[452] + v[482] + v[447] * v[491];
                /* 493 = \[OverBracket]_\[OverBracket]_\[Yen]_120(v|W)\[OverBracket]_1|3(\[Yen]|Rgm)_110(\[Yen]|W)_165 */
                v[493] = v[101] * v[103] * v[447] + v[483];
                /* 495 = \[OverBracket]_\[OverBracket]_\[Yen]_163(v|W)\[OverBracket]_2|3;\[Xi](\[Yen]|Rgm)_110(\[Yen]|W)_165 */
                v[495] = v[102] * v[103] * v[452] + v[473] + v[447] * v[494];
                /* 496 = \[OverBracket]_\[OverBracket]_\[Yen]_123(v|W)\[OverBracket]_2|3(\[Yen]|Rgm)_110(\[Yen]|W)_165 */
                v[496] = v[102] * v[103] * v[447] + v[474];
                /* 497 = \[OverBracket]_\[Yen]_127\[OverBracket]_(\[Yen]|Rgm)_110(\[Yen]|W)_165 */
                v[497] = (-1e0 / 720e0) * (v[262] * v[447]) - (v[259] * v[452]) / 720e0;
                /* 498 = \[OverBracket]_\[Yen]_134\[OverBracket]_(\[Yen]|Rgm)_110(\[Yen]|W)_165 */
                v[498] = (-1e0 / 720e0) * (v[272] * v[447]) - (v[252] * v[452]) / 720e0;
                /* 499 = \[OverBracket]_\[Yen]_152\[OverBracket]_(\[Yen]|Rgm)_110(\[Yen]|W)_165 */
                v[499] = (-1e0 / 720e0) * (v[252] * v[447]);
                /* 500 = \[OverBracket]_\[Yen]_153\[OverBracket]_(\[Yen]|Rgm)_110(\[Yen]|W)_165 */
                v[500] = (-1e0 / 720e0) * (v[259] * v[447]);
                /* 501 = \[OverBracket]_\[Yen]_139\[OverBracket]_(\[Yen]|Rgm)_110(\[Yen]|W)_165 */
                v[501] = (-1e0 / 720e0) * (v[277] * v[447]) - (v[242] * v[452]) / 720e0;
                /* 502 = \[OverBracket]_\[Yen]_149\[OverBracket]_(\[Yen]|Rgm)_110(\[Yen]|W)_165 */
                v[502] = (-1e0 / 720e0) * (v[242] * v[447]);
                /* 503 = \[OverBracket]_\[OverBracket]_R_1|1(Ft|W)\[OverBracket]_1|3(\[Yen]|Rgm)_127(v|W)_1|1;\[Xi] */
                v[503] = (-1e0 / 720e0) * (v[110] * v[424]) - (v[127] * v[447]) / 720e0;
                /* 504 = \[OverBracket]_\[OverBracket]_R_1|1;\[Xi](Fref|W)\[OverBracket]_1|1(\[Yen]|Rgm)_127(v|W)_1|1;\[Xi] */
                v[504] = (-1e0 / 720e0) * (v[110] * v[422]) - (v[162] * v[424]) / 720e0 - (v[153] * v[447]) / 720e0 - (v[127] * v[452])
                    / 720e0;
                /* 505 = \[OverBracket]_\[OverBracket]_R_2|2(Ft|W)\[OverBracket]_2|3(\[Yen]|Rgm)_134(v|W)_2|2;\[Xi] */
                v[505] = (-1e0 / 720e0) * (v[110] * v[425]) - (v[134] * v[447]) / 720e0;
                /* 506 = \[OverBracket]_\[OverBracket]_R_2|2;\[Xi](Fref|W)\[OverBracket]_2|1(\[Yen]|Rgm)_134(v|W)_2|2;\[Xi] */
                v[506] = (-1e0 / 720e0) * (v[110] * v[423]) - (v[162] * v[425]) / 720e0 - (v[152] * v[447]) / 720e0 - (v[134] * v[452])
                    / 720e0;
                /* 507 = \[OverBracket]_\[OverBracket]_R_3|3(Ft|W)\[OverBracket]_3|3(\[Yen]|Rgm)_139(v|W)_3|3;\[Xi] */
                v[507] = (-1e0 / 720e0) * (v[110] * v[413]) - (v[139] * v[447]) / 720e0;
                /* 508 = \[OverBracket]_\[OverBracket]_R_3|3;\[Xi](Fref|W)\[OverBracket]_3|1(\[Yen]|Rgm)_139(v|W)_3|3;\[Xi] */
                v[508] = (-1e0 / 720e0) * (v[110] * v[412]) - (v[162] * v[413]) / 720e0 - (v[149] * v[447]) / 720e0 - (v[139] * v[452])
                    / 720e0;
                /* 509 = \[OverBracket]_\[OverBracket]_R_2|3(Ft|W)(\[Yen]|Rgm)_2|3282 */
                v[509] = (v[485] + v[496]) / 720e0;
                /* 510 = \[OverBracket]_\[OverBracket]_R_3|2(Ft|W)(\[Yen]|Rgm)_3|3283 */
                v[510] = (-v[485] + v[496]) / 720e0;
                /* 511 = \[OverBracket]_\[OverBracket]_R_3|2;\[Xi](Fref|W)(\[Yen]|Rgm)_3|1286 */
                v[511] = (-v[484] + v[495]) / 720e0;
                /* 512 = \[OverBracket]_\[OverBracket]_R_2|3;\[Xi](Fref|W)(\[Yen]|Rgm)_2|1287 */
                v[512] = (v[484] + v[495]) / 720e0;
                /* 513 = \[OverBracket]_\[OverBracket]_R_1|3(Ft|W)(\[Yen]|Rgm)_1|3292 */
                v[513] = (v[472] + v[493]) / 720e0;
                /* 514 = \[OverBracket]_\[OverBracket]_R_3|1(Ft|W)(\[Yen]|Rgm)_3|3293 */
                v[514] = (-v[472] + v[493]) / 720e0;
                /* 515 = \[OverBracket]_\[OverBracket]_R_3|1;\[Xi](Fref|W)(\[Yen]|Rgm)_3|1296 */
                v[515] = (-v[471] + v[492]) / 720e0;
                /* 516 = \[OverBracket]_\[OverBracket]_R_1|3;\[Xi](Fref|W)(\[Yen]|Rgm)_1|1297 */
                v[516] = (v[471] + v[492]) / 720e0;
                /* 517 = \[OverBracket]_\[OverBracket]_R_1|2(Ft|W)(\[Yen]|Rgm)_1|3300 */
                v[517] = (v[465] + v[490]) / 720e0;
                /* 518 = \[OverBracket]_\[OverBracket]_R_2|1(Ft|W)(\[Yen]|Rgm)_2|3301 */
                v[518] = (-v[465] + v[490]) / 720e0;
                /* 519 = \[OverBracket]_\[OverBracket]_R_2|1;\[Xi](Fref|W)(\[Yen]|Rgm)_2|1304 */
                v[519] = (-v[464] + v[489]) / 720e0;
                /* 520 = \[OverBracket]_\[OverBracket]_R_1|2;\[Xi](Fref|W)(\[Yen]|Rgm)_1|1305 */
                v[520] = (v[464] + v[489]) / 720e0;
                /* 418 = \[OverBracket]_u\[Phi]_6(\[Yen]|Rgm)_313 */
                v[418] = v[418] + v[162] * v[486] + v[447] * (v[147] * v[285] + v[145] * v[295] + v[521] + v[522]) + v[452] * v[523];
                /* 419 = \[OverBracket]_u\[Phi]_6;\[Xi](\[Yen]|Rgm)_313 */
                v[419] = v[419] - 6e0 * v[306] * v[446] + v[110] * v[486] + v[447] * v[523];
                /* 409 = \[OverBracket]_u\[Phi]_5(\[Yen]|Rgm)_318 */
                v[409] = v[409] + v[162] * v[487] + v[447] * (v[150] * v[285] + v[145] * v[303] + v[524] + v[525]) + v[452] * v[526];
                /* 410 = \[OverBracket]_u\[Phi]_5;\[Xi](\[Yen]|Rgm)_318 */
                v[410] = v[410] + 6e0 * v[298] * v[446] + v[110] * v[487] + v[447] * v[526];
                /* 427 = \[OverBracket]_\[Yen]_110(\[Yen]|Rgm)_319 */
                v[427] = (-1e0 / 720e0) * (v[242] * v[412]) - (v[277] * v[413]) / 720e0 - (v[259] * v[422]) / 720e0 - (v[252] * v[423])
                    / 720e0 - (v[262] * v[424]) / 720e0 - (v[272] * v[425]) / 720e0 + v[427] + v[145] * v[476] + v[150] * v[486]
                    + v[147] * v[487];
                /* 435 = \[OverBracket]_\[Yen]_162(\[Yen]|Rgm)_319 */
                v[435] = (-1e0 / 720e0) * (v[242] * v[413]) - (v[259] * v[424]) / 720e0 - (v[252] * v[425]) / 720e0 + v[101] * v[476]
                    + v[103] * v[486] + v[102] * v[487];
                /* 402 = \[OverBracket]_u\[Phi]_4(\[Yen]|Rgm)_319 */
                v[402] = v[402] + v[162] * v[476] + v[447] * (v[150] * v[295] + v[147] * v[303] + v[527] + v[528]) + v[452] * v[529];
                /* 403 = \[OverBracket]_u\[Phi]_4;\[Xi](\[Yen]|Rgm)_319 */
                v[403] = v[403] - 6e0 * v[288] * v[446] + v[110] * v[476] + v[447] * v[529];
            }
            else {
                v[629] = -((v[103] * v[332]) / v[104]) - (v[102] * v[336]) / v[104];
                v[624] = -((v[103] * v[327]) / v[104]) - (v[101] * v[336]) / v[104];
                v[621] = -((v[102] * v[327]) / v[104]) - (v[101] * v[332]) / v[104];
                v[620] = v[149] / v[104] - v[139] * v[154] * v[185];
                v[618] = v[152] / v[104] - v[134] * v[154] * v[185];
                v[617] = v[153] / v[104] - v[127] * v[154] * v[185];
                v[612] = -((v[102] * v[334]) / v[104]);
                v[611] = -((v[103] * v[330]) / v[104]);
                v[596] = -((v[101] * v[334]) / v[104]);
                v[595] = -((v[103] * v[325]) / v[104]);
                v[585] = -(v[145] / v[104]) + v[101] * v[154] * v[185];
                v[584] = -((v[101] * v[330]) / v[104]);
                v[583] = -(v[147] / v[104]) + v[102] * v[154] * v[185];
                v[582] = -((v[102] * v[325]) / v[104]);
                v[571] = -((v[102] * v[145]) / v[104]) - (v[101] * v[147]) / v[104] + v[101] * v[102] * v[154] * v[185];
                v[560] = -(v[150] / v[104]) + v[103] * v[154] * v[185];
                v[547] = -((v[103] * v[145]) / v[104]) - (v[101] * v[150]) / v[104] + v[101] * v[103] * v[154] * v[185];
                v[546] = -((v[103] * v[147]) / v[104]) - (v[102] * v[150]) / v[104] + v[102] * v[103] * v[154] * v[185];
                v[551] = -((v[127] * v[262]) / v[104]) - (v[134] * v[272]) / v[104] - (v[139] * v[277]) / v[104] +
                    (v[102] * v[103] * v[325]) / v[104] + (v[101] * v[103] * v[330]) / v[104] + (v[101] * v[102] * v[334]) / v[104]
                    - v[327] * v[546] - v[332] * v[547] - v[548] - v[336] * v[571] - v[259] * v[617] - v[252] * v[618] - v[242] * v[620];
                v[544] = 1e0 / Power(v[104], 3);
                v[538] = 1e0 / Power(v[128], 3);
                /* 429 = \[OverBracket]_\[Yen]_128\[OverBracket]_(ff|Rgm)_(\[Yen]|W)_185 */
                v[429] = (-3e0 * v[175] * v[176] * v[341] * v[426]) / Power(v[128], 4);
                /* 437 = \[OverBracket]_\[Yen]_176\[OverBracket]_(ff|Rgm)_(\[Yen]|W)_185 */
                v[437] = (v[175] * v[426] * (-(v[178] * v[338]) - v[130] * v[342] + 2e0 * v[341] * v[538])) / 2e0;
                /* 499 = \[OverBracket]_\[Yen]_152\[OverBracket]_(ff|Rgm)_(\[Yen]|W)_185 */
                v[499] = (-0.5e0 * (v[130] * v[175] * v[252]) / v[104] - v[129] * v[185] * v[252]) * v[426];
                /* 500 = \[OverBracket]_\[Yen]_153\[OverBracket]_(ff|Rgm)_(\[Yen]|W)_185 */
                v[500] = (-0.5e0 * (v[130] * v[175] * v[259]) / v[104] - v[129] * v[185] * v[259]) * v[426];
                /* 502 = \[OverBracket]_\[Yen]_149\[OverBracket]_(ff|Rgm)_(\[Yen]|W)_185 */
                v[502] = (-0.5e0 * (v[130] * v[175] * v[242]) / v[104] - v[129] * v[185] * v[242]) * v[426];
                /* 536 = \[OverBracket]_\[OverBracket]_\[Yen]_176(\[Yen]|W)\[OverBracket]_179(ff|Rgm)_;\[Xi](\[Yen]|W)_176 */
                v[536] = (v[175] * v[421]) / 2e0 - (v[154] * v[175] * v[426]) / (4e0 * v[104]);
                /* 537 = \[OverBracket]_\[OverBracket]_\[Yen]_179(\[Yen]|W)\[OverBracket]_182(\[Yen]|Rgm)_176(\[Yen]|W)_179 */
                v[537] = -0.5e0 * (v[130] * v[175] * v[176] * v[426]) + v[178] * v[536];
                /* 539 = \[OverBracket]_\[OverBracket]_\[Yen]_177(\[Yen]|W)\[OverBracket]_182(\[Yen]|Rgm)_176(\[Yen]|W)_179 */
                v[539] = -(v[345] * v[536]) + v[175] * v[176] * v[426] * v[538];
                /* 540 = \[OverBracket]_\[OverBracket]_\[Yen]_183(\[Yen]|W)\[OverBracket]_187(\[Yen]|Rgm)_176(\[Yen]|W)_179 */
                v[540] = -0.5e0 * (v[175] * v[176] * v[178] * v[426]) - v[130] * v[536];
                /* 439 = \[OverBracket]_\[Yen]_178\[OverBracket]_(\[Yen]|Rgm)_176(\[Yen]|W)_179 */
                v[439] = v[342] * v[536] + (v[175] * v[426] * v[541]) / 2e0;
                /* 443 = \[OverBracket]_\[Yen]_345\[OverBracket]_(\[Yen]|Rgm)_176(\[Yen]|W)_179 */
                v[443] = -(v[341] * v[536]);
                /* 542 = \[OverBracket]_\[OverBracket]_\[Yen]_182(v|W)\[OverBracket]_1|2;\[Xi](u\[Phi]|Rgm)_6;\[Xi](\[Yen]|W)_182 */
                v[542] = -(v[130] * v[175] * v[414]) + v[426] * ((v[175] * v[178] * v[759]) / 2e0 - (v[175] * v[770]) / (2e0 * v[104]));
                /* 543 = \[OverBracket]_\[OverBracket]_\[Yen]_186(v|W)\[OverBracket]_1|3;\[Xi](u\[Phi]|Rgm)_6;\[Xi](\[Yen]|W)_182 */
                v[543] = -((v[101] * v[129] * v[414]) / v[104]) + v[101] * v[103] * v[129] * v[185] * v[421] + v[426] * (
                    -2e0 * v[101] * v[103] * v[129] * v[154] * v[544] - (v[130] * v[175] * v[547]) / 2e0 + v[762] + v[763] + v[764]);
                /* 545 = \[OverBracket]_\[OverBracket]_\[Yen]_184(v|W)\[OverBracket]_2|3;\[Xi](u\[Phi]|Rgm)_6;\[Xi](\[Yen]|W)_182 */
                v[545] = -((v[102] * v[129] * v[414]) / v[104]) + v[102] * v[103] * v[129] * v[185] * v[421] + v[426] *
                    (v[103] * v[129] * v[147] * v[185] - 2e0 * v[102] * v[103] * v[129] * v[154] * v[544] - (v[130] * v[175] * v[546]) / 2e0
                        + v[760] + v[761]);
                /* 430 = \[OverBracket]_\[Yen]_129\[OverBracket]_(u\[Phi]|Rgm)_6;\[Xi](\[Yen]|W)_182 */
                v[430] = v[414] * v[621] + v[421] * v[626] + v[426] * (-(v[149] * v[185] * v[242]) - v[152] * v[185] * v[252]
                    - v[153] * v[185] * v[259] - v[127] * v[185] * v[262] - v[134] * v[185] * v[272] - v[139] * v[185] * v[277]
                    + v[102] * v[103] * v[185] * v[325] + v[553] + v[554] + v[555] + v[556] + v[557] - 2e0 * v[544] * v[627]);
                /* 431 = \[OverBracket]_\[Yen]_130\[OverBracket]_(u\[Phi]|Rgm)_6;\[Xi](\[Yen]|W)_182 */
                v[431] = -(v[175] * v[337] * v[414]) - v[338] * v[536] + v[426] * (-0.5e0 * (v[175] * v[549]) / v[104] +
                    (v[175] * v[551]) / 2e0);
                /* 436 = \[OverBracket]_\[Yen]_175\[OverBracket]_(u\[Phi]|Rgm)_6;\[Xi](\[Yen]|W)_182 */
                v[436] = -(v[130] * v[337] * v[414]) + (v[346] * v[421]) / 2e0 + v[426] * ((v[175] * v[178] * v[549]) / 2e0 +
                    (2e0 * v[176] * v[341] * v[538] + v[550] + v[130] * v[551]) / 2e0 - v[552] / (2e0 * v[104]));
                /* 448 = \[OverBracket]_ff_\[OverBracket]_(u\[Phi]|Rgm)_6;\[Xi](\[Yen]|W)_182 */
                v[448] = v[426] * ((v[175] * v[185] * v[552]) / 2e0 - (v[130] * v[175] * (-(v[149] * v[185] * v[242])
                    - v[152] * v[185] * v[252] - v[153] * v[185] * v[259] - v[127] * v[185] * v[262] - v[134] * v[185] * v[272]
                    - v[139] * v[185] * v[277] + v[102] * v[103] * v[185] * v[325] + v[553] + v[554] + v[555] + v[556] + v[557])) / 2e0 +
                    (6e0 * v[558]) / Power(v[104], 4)) + v[414] * v[577];
                /* 409 = \[OverBracket]_u\[Phi]_5\[OverBracket]_(u\[Phi]|Rgm)_6;\[Xi](\[Yen]|W)_182 */
                v[409] = -((v[129] * v[327] * v[414]) / v[104]) + v[531] + v[421] * (v[103] * v[129] * v[185] * v[327] + v[576])
                    + v[426] * (v[103] * v[129] * v[185] * v[325] - (v[175] * (v[130] * v[331] + v[179] * v[333])) / (2e0 * v[104]) - 2e0 *
                        (v[103] * v[129] * v[154] * v[327] + v[101] * v[129] * v[154] * v[336]) * v[544] + v[327] * v[559] + v[336] * v[581] +
                        (v[175] * (v[178] * (v[175] * v[331] + v[177] * v[333]) - v[130] * (v[327] * v[560] + v[336] * v[585] + v[595] + v[596])
                            )) / 2e0 + v[597]);
                /* 402 = \[OverBracket]_u\[Phi]_4\[OverBracket]_(u\[Phi]|Rgm)_6;\[Xi](\[Yen]|W)_182 */
                v[402] = -((v[129] * v[332] * v[414]) / v[104]) + v[532] + v[421] * (v[103] * v[129] * v[185] * v[332] + v[601])
                    + v[426] * (-0.5e0 * (v[175] * (-(v[130] * v[326]) - v[179] * v[328])) / v[104] + v[103] * v[129] * v[185] * v[330]
                        - 2e0 * (v[103] * v[129] * v[154] * v[332] + v[102] * v[129] * v[154] * v[336]) * v[544] + v[332] * v[559]
                        + v[336] * v[579] + (v[175] * (v[178] * (-(v[175] * v[326]) - v[177] * v[328]) - v[130] * (v[332] * v[560]
                            + v[336] * v[583] + v[611] + v[612]))) / 2e0 + v[613]);
                /* 561 = \[OverBracket]_\[OverBracket]_\[Yen]_182(v|W)\[OverBracket]_1|2;\[Xi](u\[Phi]|Rgm)_6(\[Yen]|W)_182 */
                v[561] = v[340] * v[420] + v[542];
                /* 562 = \[OverBracket]_\[OverBracket]_\[Yen]_132(v|W)\[OverBracket]_1|2(u\[Phi]|Rgm)_6(\[Yen]|W)_182 */
                v[562] = -(v[130] * v[175] * v[420]) + ((v[103] * v[130] * v[175]) / (2e0 * v[104]) - (v[103] * v[178] * v[345]) / 2e0
                    ) * v[426];
                /* 563 = \[OverBracket]_\[OverBracket]_\[Yen]_186(v|W)\[OverBracket]_1|3;\[Xi](u\[Phi]|Rgm)_6(\[Yen]|W)_182 */
                v[563] = v[344] * v[420] + v[543];
                /* 564 = \[OverBracket]_\[OverBracket]_\[Yen]_136(v|W)\[OverBracket]_1|3(u\[Phi]|Rgm)_6(\[Yen]|W)_182 */
                v[564] = -((v[101] * v[129] * v[420]) / v[104]) + ((v[101] * v[103] * v[130] * v[175]) / (2e0 * v[104])
                    + v[101] * v[103] * v[129] * v[185]) * v[426];
                /* 565 = \[OverBracket]_\[OverBracket]_\[Yen]_184(v|W)\[OverBracket]_2|3;\[Xi](u\[Phi]|Rgm)_6(\[Yen]|W)_182 */
                v[565] = v[343] * v[420] + v[545];
                /* 566 = \[OverBracket]_\[OverBracket]_\[Yen]_138(v|W)\[OverBracket]_2|3(u\[Phi]|Rgm)_6(\[Yen]|W)_182 */
                v[566] = -((v[102] * v[129] * v[420]) / v[104]) + ((v[102] * v[103] * v[130] * v[175]) / (2e0 * v[104])
                    + v[102] * v[103] * v[129] * v[185]) * v[426];
                /* 567 = \[OverBracket]_\[Yen]_340\[OverBracket]_(u\[Phi]|Rgm)_6(\[Yen]|W)_182 */
                v[567] = v[337] * v[420];
                /* 568 = \[OverBracket]_\[Yen]_343\[OverBracket]_(u\[Phi]|Rgm)_6(\[Yen]|W)_182 */
                v[568] = v[327] * v[420];
                /* 569 = \[OverBracket]_\[Yen]_344\[OverBracket]_(u\[Phi]|Rgm)_6(\[Yen]|W)_182 */
                v[569] = v[332] * v[420];
                /* 430 = \[OverBracket]_\[Yen]_129\[OverBracket]_(u\[Phi]|Rgm)_6(\[Yen]|W)_182 */
                v[430] = v[430] + v[420] * (v[582] + v[584]);
                /* 431 = \[OverBracket]_\[Yen]_130\[OverBracket]_(u\[Phi]|Rgm)_6(\[Yen]|W)_182 */
                v[431] = -(v[175] * v[335] * v[420]) + v[431];
                /* 436 = \[OverBracket]_\[Yen]_175\[OverBracket]_(u\[Phi]|Rgm)_6(\[Yen]|W)_182 */
                v[436] = -(v[130] * v[335] * v[420]) + v[436];
                /* 448 = \[OverBracket]_ff_\[OverBracket]_(u\[Phi]|Rgm)_6(\[Yen]|W)_182 */
                v[448] = v[448] + v[420] * (v[578] + v[580]);
                /* 409 = \[OverBracket]_u\[Phi]_5\[OverBracket]_(u\[Phi]|Rgm)_6(\[Yen]|W)_182 */
                v[409] = v[409] - (v[129] * v[325] * v[420]) / v[104];
                /* 402 = \[OverBracket]_u\[Phi]_4\[OverBracket]_(u\[Phi]|Rgm)_6(\[Yen]|W)_182 */
                v[402] = v[402] - (v[129] * v[330] * v[420]) / v[104];
                /* 570 = \[OverBracket]_\[OverBracket]_\[Yen]_182(v|W)\[OverBracket]_1|2;\[Xi](\[Yen]|Rgm)_177(\[Yen]|W)_182 */
                v[570] = -(v[103] * v[175] * v[537]) - v[103] * v[130] * v[539] + v[561];
                /* 572 = \[OverBracket]_\[OverBracket]_\[Yen]_187(v|W)\[OverBracket]_1|2;\[Xi](u\[Phi]|Rgm)_5;\[Xi](\[Yen]|W)_181 */
                v[572] = -((v[101] * v[129] * v[405]) / v[104]) + v[101] * v[102] * v[129] * v[185] * v[421] + v[426] * (
                    -2e0 * v[101] * v[102] * v[129] * v[154] * v[544] - (v[130] * v[175] * v[571]) / 2e0 + v[765] + v[766] + v[767]);
                /* 573 = \[OverBracket]_\[OverBracket]_\[Yen]_181(v|W)\[OverBracket]_1|3;\[Xi](u\[Phi]|Rgm)_5;\[Xi](\[Yen]|W)_181 */
                v[573] = v[130] * v[175] * v[405] + v[102] * v[175] * v[537] + v[102] * v[130] * v[539] + v[426] * (
                    (v[175] * v[178] * v[758]) / 2e0 - (v[175] * v[769]) / (2e0 * v[104]));
                /* 574 = \[OverBracket]_\[OverBracket]_\[Yen]_184(v|W)\[OverBracket]_2|3;\[Xi](u\[Phi]|Rgm)_5;\[Xi](\[Yen]|W)_181 */
                v[574] = -((v[103] * v[129] * v[405]) / v[104]) + v[565];
                /* 430 = \[OverBracket]_\[Yen]_129\[OverBracket]_(u\[Phi]|Rgm)_5;\[Xi](\[Yen]|W)_181 */
                v[430] = v[430] + v[405] * v[624];
                /* 431 = \[OverBracket]_\[Yen]_130\[OverBracket]_(u\[Phi]|Rgm)_5;\[Xi](\[Yen]|W)_181 */
                v[431] = v[175] * v[333] * v[405] + v[431] + v[539] * v[575];
                /* 436 = \[OverBracket]_\[Yen]_175\[OverBracket]_(u\[Phi]|Rgm)_5;\[Xi](\[Yen]|W)_181 */
                v[436] = v[130] * v[333] * v[405] + v[436] + v[537] * v[575];
                /* 448 = \[OverBracket]_ff_\[OverBracket]_(u\[Phi]|Rgm)_5;\[Xi](\[Yen]|W)_181 */
                v[448] = v[448] + v[405] * (v[576] + v[625]);
                /* 418 = \[OverBracket]_u\[Phi]_6\[OverBracket]_(u\[Phi]|Rgm)_5;\[Xi](\[Yen]|W)_181 */
                v[418] = -((v[129] * v[327] * v[405]) / v[104]) + v[530] - v[175] * v[337] * v[537] - v[130] * v[337] * v[539]
                    + v[421] * v[577] + v[426] * (-0.5e0 * (v[175] * (-(v[130] * v[335]) - v[179] * v[337])) / v[104] - 2e0 *
                        (v[102] * v[129] * v[154] * v[327] + v[101] * v[129] * v[154] * v[332]) * v[544] + v[578] + v[327] * v[579] + v[580]
                        + v[332] * v[581] + (v[175] * (v[178] * (-(v[175] * v[335]) - v[177] * v[337]) - v[130] * (v[582] + v[327] * v[583]
                            + v[584] + v[332] * v[585]))) / 2e0);
                /* 402 = \[OverBracket]_u\[Phi]_4\[OverBracket]_(u\[Phi]|Rgm)_5;\[Xi](\[Yen]|W)_181 */
                v[402] = v[402] - (v[129] * v[336] * v[405]) / v[104] - v[175] * v[328] * v[537] - v[130] * v[328] * v[539];
                /* 586 = \[OverBracket]_\[OverBracket]_\[Yen]_187(v|W)\[OverBracket]_1|2;\[Xi](u\[Phi]|Rgm)_5(\[Yen]|W)_181 */
                v[586] = v[344] * v[411] + v[572];
                /* 587 = \[OverBracket]_\[OverBracket]_\[Yen]_133(v|W)\[OverBracket]_1|2(u\[Phi]|Rgm)_5(\[Yen]|W)_181 */
                v[587] = -((v[101] * v[129] * v[411]) / v[104]) + ((v[101] * v[102] * v[130] * v[175]) / (2e0 * v[104])
                    + v[101] * v[102] * v[129] * v[185]) * v[426];
                /* 588 = \[OverBracket]_\[OverBracket]_\[Yen]_181(v|W)\[OverBracket]_1|3;\[Xi](u\[Phi]|Rgm)_5(\[Yen]|W)_181 */
                v[588] = -(v[340] * v[411]) + v[573];
                /* 589 = \[OverBracket]_\[OverBracket]_\[Yen]_135(v|W)\[OverBracket]_1|3(u\[Phi]|Rgm)_5(\[Yen]|W)_181 */
                v[589] = v[130] * v[175] * v[411] + (-0.5e0 * (v[102] * v[130] * v[175]) / v[104] + (v[102] * v[178] * v[345]) / 2e0
                    ) * v[426];
                /* 590 = \[OverBracket]_\[OverBracket]_\[Yen]_184(v|W)\[OverBracket]_2|3;\[Xi](u\[Phi]|Rgm)_5(\[Yen]|W)_181 */
                v[590] = v[339] * v[411] + v[574];
                /* 591 = \[OverBracket]_\[OverBracket]_\[Yen]_138(v|W)\[OverBracket]_2|3(u\[Phi]|Rgm)_5(\[Yen]|W)_181 */
                v[591] = -((v[103] * v[129] * v[411]) / v[104]) + v[566];
                /* 592 = \[OverBracket]_\[Yen]_339\[OverBracket]_(u\[Phi]|Rgm)_5(\[Yen]|W)_181 */
                v[592] = v[327] * v[411];
                /* 593 = \[OverBracket]_\[Yen]_340\[OverBracket]_(u\[Phi]|Rgm)_5(\[Yen]|W)_181 */
                v[593] = -(v[333] * v[411]) + v[567];
                /* 594 = \[OverBracket]_\[Yen]_344\[OverBracket]_(u\[Phi]|Rgm)_5(\[Yen]|W)_181 */
                v[594] = v[336] * v[411] + v[569];
                /* 430 = \[OverBracket]_\[Yen]_129\[OverBracket]_(u\[Phi]|Rgm)_5(\[Yen]|W)_181 */
                v[430] = v[430] + v[411] * (v[595] + v[596]);
                /* 431 = \[OverBracket]_\[Yen]_130\[OverBracket]_(u\[Phi]|Rgm)_5(\[Yen]|W)_181 */
                v[431] = v[175] * v[331] * v[411] + v[431];
                /* 436 = \[OverBracket]_\[Yen]_175\[OverBracket]_(u\[Phi]|Rgm)_5(\[Yen]|W)_181 */
                v[436] = v[130] * v[331] * v[411] + v[436];
                /* 448 = \[OverBracket]_ff_\[OverBracket]_(u\[Phi]|Rgm)_5(\[Yen]|W)_181 */
                v[448] = v[448] + v[411] * (v[103] * v[129] * v[185] * v[325] + v[597]);
                /* 418 = \[OverBracket]_u\[Phi]_6\[OverBracket]_(u\[Phi]|Rgm)_5(\[Yen]|W)_181 */
                v[418] = -((v[129] * v[325] * v[411]) / v[104]) + v[418];
                /* 402 = \[OverBracket]_u\[Phi]_4\[OverBracket]_(u\[Phi]|Rgm)_5(\[Yen]|W)_181 */
                v[402] = v[402] - (v[129] * v[334] * v[411]) / v[104];
                /* 598 = \[OverBracket]_\[OverBracket]_\[Yen]_187(v|W)\[OverBracket]_1|2;\[Xi](u\[Phi]|Rgm)_4;\[Xi](\[Yen]|W)_180 */
                v[598] = -((v[102] * v[129] * v[398]) / v[104]) + v[586];
                /* 599 = \[OverBracket]_\[OverBracket]_\[Yen]_186(v|W)\[OverBracket]_1|3;\[Xi](u\[Phi]|Rgm)_4;\[Xi](\[Yen]|W)_180 */
                v[599] = -((v[103] * v[129] * v[398]) / v[104]) + v[563];
                /* 600 = \[OverBracket]_\[OverBracket]_\[Yen]_180(v|W)\[OverBracket]_2|3;\[Xi](u\[Phi]|Rgm)_4;\[Xi](\[Yen]|W)_180 */
                v[600] = -(v[130] * v[175] * v[398]) - v[101] * v[175] * v[537] - v[101] * v[130] * v[539] + v[426] * (
                    (v[175] * v[178] * v[757]) / 2e0 - (v[175] * v[768]) / (2e0 * v[104]));
                /* 430 = \[OverBracket]_\[Yen]_129\[OverBracket]_(u\[Phi]|Rgm)_4;\[Xi](\[Yen]|W)_180 */
                v[430] = v[430] + v[398] * v[629];
                /* 431 = \[OverBracket]_\[Yen]_130\[OverBracket]_(u\[Phi]|Rgm)_4;\[Xi](\[Yen]|W)_180 */
                v[431] = -(v[175] * v[328] * v[398]) + v[431];
                /* 436 = \[OverBracket]_\[Yen]_175\[OverBracket]_(u\[Phi]|Rgm)_4;\[Xi](\[Yen]|W)_180 */
                v[436] = -(v[130] * v[328] * v[398]) + v[436];
                /* 448 = \[OverBracket]_ff_\[OverBracket]_(u\[Phi]|Rgm)_4;\[Xi](\[Yen]|W)_180 */
                v[448] = v[448] + v[398] * (v[601] + v[630]);
                /* 418 = \[OverBracket]_u\[Phi]_6\[OverBracket]_(u\[Phi]|Rgm)_4;\[Xi](\[Yen]|W)_180 */
                v[418] = -((v[129] * v[332] * v[398]) / v[104]) + v[418];
                /* 409 = \[OverBracket]_u\[Phi]_5\[OverBracket]_(u\[Phi]|Rgm)_4;\[Xi](\[Yen]|W)_180 */
                v[409] = -((v[129] * v[336] * v[398]) / v[104]) + v[409] + v[175] * v[333] * v[537] + v[130] * v[333] * v[539];
                /* 602 = \[OverBracket]_\[OverBracket]_\[Yen]_187(v|W)\[OverBracket]_1|2;\[Xi](u\[Phi]|Rgm)_4(\[Yen]|W)_180 */
                v[602] = v[343] * v[404] + v[598];
                /* 603 = \[OverBracket]_\[OverBracket]_\[Yen]_133(v|W)\[OverBracket]_1|2(u\[Phi]|Rgm)_4(\[Yen]|W)_180 */
                v[603] = -((v[102] * v[129] * v[404]) / v[104]) + v[587];
                /* 604 = \[OverBracket]_\[OverBracket]_\[Yen]_186(v|W)\[OverBracket]_1|3;\[Xi](u\[Phi]|Rgm)_4(\[Yen]|W)_180 */
                v[604] = v[339] * v[404] + v[599];
                /* 605 = \[OverBracket]_\[OverBracket]_\[Yen]_136(v|W)\[OverBracket]_1|3(u\[Phi]|Rgm)_4(\[Yen]|W)_180 */
                v[605] = -((v[103] * v[129] * v[404]) / v[104]) + v[564];
                /* 606 = \[OverBracket]_\[OverBracket]_\[Yen]_180(v|W)\[OverBracket]_2|3;\[Xi](u\[Phi]|Rgm)_4(\[Yen]|W)_180 */
                v[606] = v[340] * v[404] + v[600];
                /* 607 = \[OverBracket]_\[OverBracket]_\[Yen]_137(v|W)\[OverBracket]_2|3(u\[Phi]|Rgm)_4(\[Yen]|W)_180 */
                v[607] = -(v[130] * v[175] * v[404]) + ((v[101] * v[130] * v[175]) / (2e0 * v[104]) - (v[101] * v[178] * v[345]) / 2e0
                    ) * v[426];
                /* 608 = \[OverBracket]_\[Yen]_339\[OverBracket]_(u\[Phi]|Rgm)_4(\[Yen]|W)_180 */
                v[608] = v[332] * v[404] + v[592];
                /* 609 = \[OverBracket]_\[Yen]_340\[OverBracket]_(u\[Phi]|Rgm)_4(\[Yen]|W)_180 */
                v[609] = v[328] * v[404] + v[593];
                /* 610 = \[OverBracket]_\[Yen]_343\[OverBracket]_(u\[Phi]|Rgm)_4(\[Yen]|W)_180 */
                v[610] = v[336] * v[404] + v[568];
                /* 430 = \[OverBracket]_\[Yen]_129\[OverBracket]_(u\[Phi]|Rgm)_4(\[Yen]|W)_180 */
                v[430] = v[430] + v[404] * (v[611] + v[612]);
                /* 431 = \[OverBracket]_\[Yen]_130\[OverBracket]_(u\[Phi]|Rgm)_4(\[Yen]|W)_180 */
                v[431] = -(v[175] * v[326] * v[404]) + v[431];
                /* 436 = \[OverBracket]_\[Yen]_175\[OverBracket]_(u\[Phi]|Rgm)_4(\[Yen]|W)_180 */
                v[436] = -(v[130] * v[326] * v[404]) + v[436];
                /* 448 = \[OverBracket]_ff_\[OverBracket]_(u\[Phi]|Rgm)_4(\[Yen]|W)_180 */
                v[448] = v[448] + v[404] * (v[103] * v[129] * v[185] * v[330] + v[613]);
                /* 418 = \[OverBracket]_u\[Phi]_6\[OverBracket]_(u\[Phi]|Rgm)_4(\[Yen]|W)_180 */
                v[418] = -((v[129] * v[330] * v[404]) / v[104]) + v[418];
                /* 409 = \[OverBracket]_u\[Phi]_5\[OverBracket]_(u\[Phi]|Rgm)_4(\[Yen]|W)_180 */
                v[409] = -((v[129] * v[334] * v[404]) / v[104]) + v[409];
                /* 614 = \[OverBracket]_\[OverBracket]_\[Yen]_187(v|W)\[OverBracket]_1|2;\[Xi](\[Yen]|Rgm)_183(\[Yen]|W)_187 */
                v[614] = -((v[101] * v[102] * v[540]) / v[104]) + v[602];
                /* 615 = \[OverBracket]_\[OverBracket]_\[Yen]_186(v|W)\[OverBracket]_1|3;\[Xi](\[Yen]|Rgm)_183(\[Yen]|W)_187 */
                v[615] = -((v[101] * v[103] * v[540]) / v[104]) + v[604];
                /* 616 = \[OverBracket]_\[OverBracket]_\[Yen]_184(v|W)\[OverBracket]_2|3;\[Xi](\[Yen]|Rgm)_183(\[Yen]|W)_187 */
                v[616] = -((v[102] * v[103] * v[540]) / v[104]) + v[590];
                /* 497 = \[OverBracket]_\[Yen]_127\[OverBracket]_(\[Yen]|Rgm)_183(\[Yen]|W)_187 */
                v[497] = -(v[129] * v[185] * v[259] * v[421]) + (v[259] * v[540]) / v[104] + v[426] * (-(v[183] * v[185] * v[259])
                    - v[129] * v[185] * v[262] - (v[130] * v[175] * (-(v[154] * v[185] * v[259]) + v[262] / v[104])) / 2e0
                    + 2e0 * v[129] * v[154] * v[259] * v[544]);
                /* 498 = \[OverBracket]_\[Yen]_134\[OverBracket]_(\[Yen]|Rgm)_183(\[Yen]|W)_187 */
                v[498] = -(v[129] * v[185] * v[252] * v[421]) + (v[252] * v[540]) / v[104] + v[426] * (-(v[183] * v[185] * v[252])
                    - v[129] * v[185] * v[272] - (v[130] * v[175] * (-(v[154] * v[185] * v[252]) + v[272] / v[104])) / 2e0
                    + 2e0 * v[129] * v[154] * v[252] * v[544]);
                /* 501 = \[OverBracket]_\[Yen]_139\[OverBracket]_(\[Yen]|Rgm)_183(\[Yen]|W)_187 */
                v[501] = -(v[129] * v[185] * v[242] * v[421]) + (v[242] * v[540]) / v[104] + v[426] * (-(v[183] * v[185] * v[242])
                    - v[129] * v[185] * v[277] - (v[130] * v[175] * (-(v[154] * v[185] * v[242]) + v[277] / v[104])) / 2e0
                    + 2e0 * v[129] * v[154] * v[242] * v[544]);
                /* 503 = \[OverBracket]_\[OverBracket]_R_1|1(Ft|W)\[OverBracket]_1|3(\[Yen]|Rgm)_127(v|W)_1|1;\[Xi] */
                v[503] = (v[129] * v[424]) / v[104] + (-0.5e0 * (v[127] * v[130] * v[175]) / v[104] - v[127] * v[129] * v[185]
                    ) * v[426];
                /* 504 = \[OverBracket]_\[OverBracket]_R_1|1;\[Xi](Fref|W)\[OverBracket]_1|1(\[Yen]|Rgm)_127(v|W)_1|1;\[Xi] */
                v[504] = -(v[127] * v[129] * v[185] * v[421]) + (v[129] * v[422]) / v[104] + v[329] * v[424] + (v[127] * v[540])
                    / v[104] + v[426] * (-(v[129] * v[153] * v[185]) - v[127] * v[183] * v[185] + 2e0 * v[127] * v[129] * v[154] * v[544] -
                        (v[130] * v[175] * v[617]) / 2e0);
                /* 520 = \[OverBracket]_\[OverBracket]_R_1|2;\[Xi](Fref|W)\[OverBracket]_1|1(\[Yen]|Rgm)_187(v|W)_1|2;\[Xi] */
                v[520] = v[570] + v[614];
                /* 519 = \[OverBracket]_\[OverBracket]_R_2|1;\[Xi](Fref|W)\[OverBracket]_2|1(\[Yen]|Rgm)_187(v|W)_1|2;\[Xi] */
                v[519] = -v[570] + v[614];
                /* 518 = \[OverBracket]_\[OverBracket]_R_2|1(Ft|W)\[OverBracket]_2|3(\[Yen]|Rgm)_133(v|W)_1|2 */
                v[518] = -v[562] + v[603];
                /* 517 = \[OverBracket]_\[OverBracket]_R_1|2(Ft|W)\[OverBracket]_1|3(\[Yen]|Rgm)_133(v|W)_1|2 */
                v[517] = v[562] + v[603];
                /* 516 = \[OverBracket]_\[OverBracket]_R_1|3;\[Xi](Fref|W)\[OverBracket]_1|1(\[Yen]|Rgm)_186(v|W)_1|3;\[Xi] */
                v[516] = v[588] + v[615];
                /* 515 = \[OverBracket]_\[OverBracket]_R_3|1;\[Xi](Fref|W)\[OverBracket]_3|1(\[Yen]|Rgm)_186(v|W)_1|3;\[Xi] */
                v[515] = -v[588] + v[615];
                /* 514 = \[OverBracket]_\[OverBracket]_R_3|1(Ft|W)\[OverBracket]_3|3(\[Yen]|Rgm)_136(v|W)_1|3 */
                v[514] = -v[589] + v[605];
                /* 513 = \[OverBracket]_\[OverBracket]_R_1|3(Ft|W)\[OverBracket]_1|3(\[Yen]|Rgm)_136(v|W)_1|3 */
                v[513] = v[589] + v[605];
                /* 505 = \[OverBracket]_\[OverBracket]_R_2|2(Ft|W)\[OverBracket]_2|3(\[Yen]|Rgm)_134(v|W)_2|2;\[Xi] */
                v[505] = (v[129] * v[425]) / v[104] + (-0.5e0 * (v[130] * v[134] * v[175]) / v[104] - v[129] * v[134] * v[185]
                    ) * v[426];
                /* 506 = \[OverBracket]_\[OverBracket]_R_2|2;\[Xi](Fref|W)\[OverBracket]_2|1(\[Yen]|Rgm)_134(v|W)_2|2;\[Xi] */
                v[506] = -(v[129] * v[134] * v[185] * v[421]) + (v[129] * v[423]) / v[104] + v[329] * v[425] + (v[134] * v[540])
                    / v[104] + v[426] * (-(v[129] * v[152] * v[185]) - v[134] * v[183] * v[185] + 2e0 * v[129] * v[134] * v[154] * v[544] -
                        (v[130] * v[175] * v[618]) / 2e0);
                /* 512 = \[OverBracket]_\[OverBracket]_R_2|3;\[Xi](Fref|W)\[OverBracket]_2|1(\[Yen]|Rgm)_184(v|W)_2|3;\[Xi] */
                v[512] = v[606] + v[616];
                /* 511 = \[OverBracket]_\[OverBracket]_R_3|2;\[Xi](Fref|W)\[OverBracket]_3|1(\[Yen]|Rgm)_184(v|W)_2|3;\[Xi] */
                v[511] = -v[606] + v[616];
                /* 510 = \[OverBracket]_\[OverBracket]_R_3|2(Ft|W)\[OverBracket]_3|3(\[Yen]|Rgm)_138(v|W)_2|3 */
                v[510] = v[591] - v[607];
                /* 509 = \[OverBracket]_\[OverBracket]_R_2|3(Ft|W)\[OverBracket]_2|3(\[Yen]|Rgm)_138(v|W)_2|3 */
                v[509] = v[591] + v[607];
                /* 619 = \[OverBracket]_\[Yen]_329\[OverBracket]_(\[Yen]|Rgm)_139(v|W)_3|3;\[Xi] */
                v[619] = v[242] * v[413] + v[259] * v[424] + v[252] * v[425];
                /* 507 = \[OverBracket]_\[OverBracket]_R_3|3(Ft|W)\[OverBracket]_3|3(\[Yen]|Rgm)_139(v|W)_3|3;\[Xi] */
                v[507] = (v[129] * v[413]) / v[104] + (-0.5e0 * (v[130] * v[139] * v[175]) / v[104] - v[129] * v[139] * v[185]
                    ) * v[426];
                /* 508 = \[OverBracket]_\[OverBracket]_R_3|3;\[Xi](Fref|W)\[OverBracket]_3|1(\[Yen]|Rgm)_139(v|W)_3|3;\[Xi] */
                v[508] = (v[129] * v[412]) / v[104] + v[329] * v[413] - v[129] * v[139] * v[185] * v[421] + (v[139] * v[540]) / v[104]
                    + v[426] * (-(v[129] * v[149] * v[185]) - v[139] * v[183] * v[185] + 2e0 * v[129] * v[139] * v[154] * v[544] -
                        (v[130] * v[175] * v[620]) / 2e0);
                /* 418 = \[OverBracket]_u\[Phi]_6(\[Yen]|Rgm)_339 */
                v[418] = v[418] - v[329] * v[608] + v[540] * v[621];
                /* 419 = \[OverBracket]_u\[Phi]_6;\[Xi](\[Yen]|Rgm)_339 */
                v[419] = v[535] - (v[129] * v[608]) / v[104] + v[426] * ((v[130] * v[175] * v[337]) / (2e0 * v[104]) + (v[175] * (-
                    (v[175] * v[178] * v[337]) - v[130] * v[621])) / 2e0 + v[622] + v[623]);
                /* 431 = \[OverBracket]_\[Yen]_130(\[Yen]|Rgm)_340 */
                v[431] = v[431] - v[177] * v[609];
                /* 438 = \[OverBracket]_\[Yen]_177(\[Yen]|Rgm)_340 */
                v[438] = (v[175] * v[178] * v[426] * v[575]) / 2e0 - v[130] * v[609];
                /* 440 = \[OverBracket]_\[Yen]_179(\[Yen]|Rgm)_340 */
                v[440] = -0.5e0 * (v[175] * v[426] * v[575]) / v[104] - v[175] * v[609];
                /* 436 = \[OverBracket]_\[Yen]_175(\[Yen]|Rgm)_340 */
                v[436] = v[436] - v[179] * v[609];
                /* 409 = \[OverBracket]_u\[Phi]_5(\[Yen]|Rgm)_343 */
                v[409] = v[409] - v[329] * v[610] + v[540] * v[624];
                /* 410 = \[OverBracket]_u\[Phi]_5;\[Xi](\[Yen]|Rgm)_343 */
                v[410] = v[534] - (v[129] * v[610]) / v[104] + v[426] * (-0.5e0 * (v[130] * v[175] * v[333]) / v[104] + v[576] +
                    (v[175] * (v[175] * v[178] * v[333] - v[130] * v[624])) / 2e0 + v[625]);
                /* 430 = \[OverBracket]_\[Yen]_129(\[Yen]|Rgm)_344 */
                v[430] = (v[242] * v[412]) / v[104] + (v[277] * v[413]) / v[104] + (v[259] * v[422]) / v[104] + (v[252] * v[423])
                    / v[104] + (v[262] * v[424]) / v[104] + (v[272] * v[425]) / v[104] + v[430] + v[585] * v[594] + v[560] * v[608]
                    + v[583] * v[610] - v[154] * v[185] * v[619];
                /* 441 = \[OverBracket]_\[Yen]_183(\[Yen]|Rgm)_344 */
                v[441] = -((v[101] * v[594]) / v[104]) - (v[103] * v[608]) / v[104] - (v[102] * v[610]) / v[104] + v[619] / v[104]
                    + v[426] * v[626];
                /* 442 = \[OverBracket]_\[Yen]_185(\[Yen]|Rgm)_344 */
                v[442] = v[101] * v[129] * v[154] * v[594] + v[103] * v[129] * v[154] * v[608] + v[102] * v[129] * v[154] * v[610]
                    - v[129] * v[154] * v[619] + v[426] * ((-(v[129] * v[149]) - v[139] * v[183]) * v[242] + (-(v[129] * v[152])
                        - v[134] * v[183]) * v[252] + (-(v[129] * v[153]) - v[127] * v[183]) * v[259] - v[127] * v[129] * v[262]
                        - v[129] * v[134] * v[272] - v[129] * v[139] * v[277] + v[102] * v[103] * v[129] * v[325] + (v[103] * v[129] * v[147]
                            + v[102] * v[129] * v[150] + v[102] * v[103] * v[183]) * v[327] + v[101] * v[103] * v[129] * v[330] +
                        (v[103] * v[129] * v[145] + v[101] * v[129] * v[150] + v[101] * v[103] * v[183]) * v[332]
                        + v[101] * v[102] * v[129] * v[334] + (v[102] * v[129] * v[145] + v[101] * v[129] * v[147] + v[101] * v[102] * v[183]
                            ) * v[336] - (v[130] * v[175] * v[627]) / 2e0) + v[421] * v[628];
                /* 448 = \[OverBracket]_ff_(\[Yen]|Rgm)_344 */
                v[448] = -(v[129] * v[185] * v[242] * v[412]) - v[129] * v[185] * v[277] * v[413] - v[129] * v[185] * v[259] * v[422]
                    - v[129] * v[185] * v[252] * v[423] - v[129] * v[185] * v[262] * v[424] - v[129] * v[185] * v[272] * v[425] + v[448]
                    + v[581] * v[594] + v[559] * v[608] + v[579] * v[610] - v[183] * v[185] * v[619] + v[540] * v[626];
                /* 449 = \[OverBracket]_ff_;\[Xi](\[Yen]|Rgm)_344 */
                v[449] = v[101] * v[129] * v[185] * v[594] + v[103] * v[129] * v[185] * v[608] + v[102] * v[129] * v[185] * v[610]
                    - v[129] * v[185] * v[619] + v[426] * (-0.25e0 * (v[175] * v[346]) / v[104] - (v[130] * v[175] * v[626]) / 2e0
                        - 2e0 * v[544] * v[628]);
                /* 402 = \[OverBracket]_u\[Phi]_4(\[Yen]|Rgm)_344 */
                v[402] = v[402] - v[329] * v[594] + v[540] * v[629];
                /* 403 = \[OverBracket]_u\[Phi]_4;\[Xi](\[Yen]|Rgm)_344 */
                v[403] = v[533] - (v[129] * v[594]) / v[104] + v[426] * ((v[130] * v[175] * v[328]) / (2e0 * v[104]) + v[601] + (v[175] *
                    (-(v[175] * v[178] * v[328]) - v[130] * v[629])) / 2e0 + v[630]);
            };
            v[733] = v[449];
            v[731] = v[448];
            v[756] = v[403];
            v[755] = v[419];
            v[754] = v[410];
            v[753] = v[418];
            v[752] = v[402];
            v[751] = v[409];
            v[750] = v[500];
            v[749] = v[497];
            v[740] = v[499];
            v[739] = v[498];
            v[734] = v[502];
            v[732] = v[501];
            /* 631 = \[OverBracket]_\[OverBracket]_Ft_3|2(Fref|W)\[OverBracket]_3|2(R|Rgm)_3|3(Ft|W)_3|3 */
            v[631] = v[507] * v[69];
            /* 632 = \[OverBracket]_\[OverBracket]_Ft_3|3(Fref|W)\[OverBracket]_3|3(R|Rgm)_3|3(Ft|W)_3|3 */
            v[632] = v[507] * v[72];
            /* 633 = \[OverBracket]_\[OverBracket]_Fref_3|1(Ft|W)\[OverBracket]_3|1(R|Rgm)_3|3(Ft|W)_3|3 */
            v[1390] = 0e0;
            v[1391] = 0e0;
            v[1392] = v[56];
            v[1393] = 0e0;
            v[1394] = 0e0;
            v[1395] = 0e0;
            v[1396] = 0e0;
            v[1397] = 0e0;
            v[1398] = v[57];
            v[1399] = 0e0;
            v[1400] = 0e0;
            v[1401] = 0e0;
            v[1402] = 0e0;
            v[1403] = 0e0;
            v[1404] = v[59];
            v[1405] = 0e0;
            v[1406] = 0e0;
            v[1407] = 0e0;
            v[633] = v[273] * v[507] + v[1389 + i227] * v[97];
            /* 634 = \[OverBracket]_\[OverBracket]_Ft_3|2(Fref|W)\[OverBracket]_3|2(R|Rgm)_3|2(Ft|W)_3|3 */
            v[634] = v[631] + v[510] * v[68];
            /* 635 = \[OverBracket]_\[OverBracket]_Ft_3|3(Fref|W)\[OverBracket]_3|3(R|Rgm)_3|2(Ft|W)_3|3 */
            v[635] = v[632] + v[510] * v[71];
            /* 636 = \[OverBracket]_\[OverBracket]_Fref_3|1(Ft|W)\[OverBracket]_3|1(R|Rgm)_3|2(Ft|W)_3|3 */
            v[636] = v[271] * v[510] + v[633];
            /* 637 = \[OverBracket]_\[OverBracket]_Ft_3|2(Fref|W)\[OverBracket]_3|2(R|Rgm)_3|1(Ft|W)_3|3 */
            v[637] = v[634] + v[514] * v[67];
            /* 638 = \[OverBracket]_\[OverBracket]_Ft_3|3(Fref|W)\[OverBracket]_3|3(R|Rgm)_3|1(Ft|W)_3|3 */
            v[638] = v[635] + v[514] * v[70];
            /* 639 = \[OverBracket]_\[OverBracket]_Fref_3|1(Ft|W)\[OverBracket]_3|1(R|Rgm)_3|1(Ft|W)_3|3 */
            v[639] = v[269] * v[514] + v[636];
            /* 640 = \[OverBracket]_\[OverBracket]_Ft_2|2(Fref|W)\[OverBracket]_2|2(R|Rgm)_2|3(Ft|W)_2|3 */
            v[640] = v[509] * v[69];
            /* 641 = \[OverBracket]_\[OverBracket]_Ft_2|3(Fref|W)\[OverBracket]_2|3(R|Rgm)_2|3(Ft|W)_2|3 */
            v[641] = v[509] * v[72];
            /* 642 = \[OverBracket]_\[OverBracket]_Fref_2|1(Ft|W)\[OverBracket]_2|1(R|Rgm)_2|3(Ft|W)_2|3 */
            v[1372] = 0e0;
            v[1373] = v[56];
            v[1374] = 0e0;
            v[1375] = 0e0;
            v[1376] = 0e0;
            v[1377] = 0e0;
            v[1378] = 0e0;
            v[1379] = v[57];
            v[1380] = 0e0;
            v[1381] = 0e0;
            v[1382] = 0e0;
            v[1383] = 0e0;
            v[1384] = 0e0;
            v[1385] = v[59];
            v[1386] = 0e0;
            v[1387] = 0e0;
            v[1388] = 0e0;
            v[1389] = 0e0;
            v[642] = v[273] * v[509] + v[1371 + i227] * v[97];
            /* 643 = \[OverBracket]_\[OverBracket]_Ft_2|2(Fref|W)\[OverBracket]_2|2(R|Rgm)_2|2(Ft|W)_2|3 */
            v[643] = v[640] + v[505] * v[68];
            /* 644 = \[OverBracket]_\[OverBracket]_Ft_2|3(Fref|W)\[OverBracket]_2|3(R|Rgm)_2|2(Ft|W)_2|3 */
            v[644] = v[641] + v[505] * v[71];
            /* 645 = \[OverBracket]_\[OverBracket]_Fref_2|1(Ft|W)\[OverBracket]_2|1(R|Rgm)_2|2(Ft|W)_2|3 */
            v[645] = v[271] * v[505] + v[642];
            /* 646 = \[OverBracket]_\[OverBracket]_Ft_2|2(Fref|W)\[OverBracket]_2|2(R|Rgm)_2|1(Ft|W)_2|3 */
            v[646] = v[643] + v[518] * v[67];
            /* 647 = \[OverBracket]_\[OverBracket]_Ft_2|3(Fref|W)\[OverBracket]_2|3(R|Rgm)_2|1(Ft|W)_2|3 */
            v[647] = v[644] + v[518] * v[70];
            /* 648 = \[OverBracket]_\[OverBracket]_Fref_2|1(Ft|W)\[OverBracket]_2|1(R|Rgm)_2|1(Ft|W)_2|3 */
            v[648] = v[269] * v[518] + v[645];
            /* 649 = \[OverBracket]_\[OverBracket]_Ft_1|2(Fref|W)\[OverBracket]_1|2(R|Rgm)_1|3(Ft|W)_1|3 */
            v[649] = v[513] * v[69];
            /* 650 = \[OverBracket]_\[OverBracket]_Ft_1|3(Fref|W)\[OverBracket]_1|3(R|Rgm)_1|3(Ft|W)_1|3 */
            v[650] = v[513] * v[72];
            /* 651 = \[OverBracket]_\[OverBracket]_Fref_1|1(Ft|W)\[OverBracket]_1|1(R|Rgm)_1|3(Ft|W)_1|3 */
            v[1354] = v[56];
            v[1355] = 0e0;
            v[1356] = 0e0;
            v[1357] = 0e0;
            v[1358] = 0e0;
            v[1359] = 0e0;
            v[1360] = v[57];
            v[1361] = 0e0;
            v[1362] = 0e0;
            v[1363] = 0e0;
            v[1364] = 0e0;
            v[1365] = 0e0;
            v[1366] = v[59];
            v[1367] = 0e0;
            v[1368] = 0e0;
            v[1369] = 0e0;
            v[1370] = 0e0;
            v[1371] = 0e0;
            v[651] = v[273] * v[513] + v[1353 + i227] * v[97];
            /* 652 = \[OverBracket]_\[OverBracket]_Ft_1|2(Fref|W)\[OverBracket]_1|2(R|Rgm)_1|2(Ft|W)_1|3 */
            v[652] = v[649] + v[517] * v[68];
            /* 653 = \[OverBracket]_\[OverBracket]_Ft_1|3(Fref|W)\[OverBracket]_1|3(R|Rgm)_1|2(Ft|W)_1|3 */
            v[653] = v[650] + v[517] * v[71];
            /* 654 = \[OverBracket]_\[OverBracket]_Fref_1|1(Ft|W)\[OverBracket]_1|1(R|Rgm)_1|2(Ft|W)_1|3 */
            v[654] = v[271] * v[517] + v[651];
            /* 655 = \[OverBracket]_\[OverBracket]_Ft_1|2(Fref|W)\[OverBracket]_1|2(R|Rgm)_1|1(Ft|W)_1|3 */
            v[655] = v[652] + v[503] * v[67];
            /* 656 = \[OverBracket]_\[OverBracket]_Ft_1|3(Fref|W)\[OverBracket]_1|3(R|Rgm)_1|1(Ft|W)_1|3 */
            v[656] = v[653] + v[503] * v[70];
            /* 657 = \[OverBracket]_\[OverBracket]_Fref_1|1(Ft|W)\[OverBracket]_1|1(R|Rgm)_1|1(Ft|W)_1|3 */
            v[657] = v[269] * v[503] + v[654];
            /* 658 = \[OverBracket]_\[OverBracket]_Fref_1|1(Ft|W)\[OverBracket]_1|1(R|Rgm)_1|1;\[Xi](Fref|W)_1|1 */
            v[658] = v[253] * v[504] + v[657];
            /* 659 = \[OverBracket]_\[OverBracket]_Fref_1|1(Ft|W)\[OverBracket]_1|1(R|Rgm)_1|2;\[Xi](Fref|W)_1|1 */
            v[659] = v[251] * v[520] + v[658];
            /* 660 = \[OverBracket]_\[OverBracket]_Fref_1|1(Ft|W)\[OverBracket]_1|1(R|Rgm)_1|3;\[Xi](Fref|W)_1|1 */
            v[660] = v[249] * v[516] + v[659];
            /* 661 = \[OverBracket]_\[OverBracket]_Fref_2|1(Ft|W)\[OverBracket]_2|1(R|Rgm)_2|1;\[Xi](Fref|W)_2|1 */
            v[661] = v[253] * v[519] + v[648];
            /* 662 = \[OverBracket]_\[OverBracket]_Fref_2|1(Ft|W)\[OverBracket]_2|1(R|Rgm)_2|2;\[Xi](Fref|W)_2|1 */
            v[662] = v[251] * v[506] + v[661];
            /* 663 = \[OverBracket]_\[OverBracket]_Fref_2|1(Ft|W)\[OverBracket]_2|1(R|Rgm)_2|3;\[Xi](Fref|W)_2|1 */
            v[663] = v[249] * v[512] + v[662];
            /* 664 = \[OverBracket]_\[OverBracket]_Fref_3|1(Ft|W)\[OverBracket]_3|1(R|Rgm)_3|1;\[Xi](Fref|W)_3|1 */
            v[664] = v[253] * v[515] + v[639];
            /* 665 = \[OverBracket]_\[OverBracket]_Fref_3|1(Ft|W)\[OverBracket]_3|1(R|Rgm)_3|2;\[Xi](Fref|W)_3|1 */
            v[665] = v[251] * v[511] + v[664];
            /* 666 = \[OverBracket]_\[OverBracket]_Fref_3|1(Ft|W)\[OverBracket]_3|1(R|Rgm)_3|3;\[Xi](Fref|W)_3|1 */
            v[666] = v[249] * v[508] + v[665];
            /* 667 = \[OverBracket]_\[Yen]_239\[OverBracket]_(Ft|Rgm)_3|3(Fref|W)_3|3 */
            v[667] = v[214] * v[637] * v[9] + v[213] * v[638] * v[9] + v[211] * v[646] * v[9] + v[210] * v[647] * v[9]
                + v[199] * v[655] * v[9] + v[198] * v[656] * v[9];
            /* 668 = \[OverBracket]_\[OverBracket]_Ft_1|1(W|W)\[OverBracket]_(Fref|Rgm)_1|1(Ft|W)_1|1 */
            v[668] = (v[10] * v[247] * v[655]) / 2e0 + (v[11] * v[245] * v[656]) / 2e0 + v[236] * v[660];
            /* 669 = \[OverBracket]_\[OverBracket]_Ft_2|1(W|W)\[OverBracket]_(Fref|Rgm)_2|1(Ft|W)_2|1 */
            v[669] = (v[10] * v[247] * v[646]) / 2e0 + (v[11] * v[245] * v[647]) / 2e0 + v[236] * v[663];
            /* 670 = \[OverBracket]_\[OverBracket]_Ft_3|1(W|W)\[OverBracket]_(Fref|Rgm)_3|1(Ft|W)_3|1 */
            v[670] = (v[10] * v[247] * v[637]) / 2e0 + (v[11] * v[245] * v[638]) / 2e0 + v[236] * v[666];
            /* 671 = \[OverBracket]_\[Yen]_229\[OverBracket]_(Ft|Rgm)_3|1(W|W) */
            v[671] = v[212] * v[637] * v[7] + v[209] * v[646] * v[7] + v[197] * v[655] * v[7] + v[198] * v[668] * v[7]
                + v[210] * v[669] * v[7] + v[213] * v[670] * v[7];
            /* 672 = \[OverBracket]_\[Yen]_230\[OverBracket]_(Ft|Rgm)_3|1(W|W) */
            v[672] = v[212] * v[638] * v[8] + v[209] * v[647] * v[8] + v[197] * v[656] * v[8] + v[199] * v[668] * v[8]
                + v[211] * v[669] * v[8] + v[214] * v[670] * v[8];
            /* 673 = \[OverBracket]_\[Yen]_231\[OverBracket]_(Ft|Rgm)_3|1(W|W) */
            v[673] = (v[197] * v[225] * v[668]) / (2e0 * v[2] * v[224] * v[3]) + (v[209] * v[225] * v[669]) / (2e0 * v[2] * v[224] * v[3]
                ) + (v[212] * v[225] * v[670]) / (2e0 * v[2] * v[224] * v[3]);
            /* 674 = \[OverBracket]_Ft_3|1(\[Yen]|Rgm)_231 */
            v[674] = (v[225] * v[231] * v[670]) / (2e0 * v[2] * v[224] * v[3]) + v[213] * v[671] + v[214] * v[672]
                + 2e0 * v[212] * v[673] + v[229] * v[637] * v[7] + v[230] * v[638] * v[8];
            /* 675 = \[OverBracket]_Ft_2|1(\[Yen]|Rgm)_231 */
            v[675] = (v[225] * v[231] * v[669]) / (2e0 * v[2] * v[224] * v[3]) + v[210] * v[671] + v[211] * v[672]
                + 2e0 * v[209] * v[673] + v[229] * v[646] * v[7] + v[230] * v[647] * v[8];
            /* 676 = \[OverBracket]_Ft_1|1(\[Yen]|Rgm)_231 */
            v[676] = (v[225] * v[231] * v[668]) / (2e0 * v[2] * v[224] * v[3]) + v[198] * v[671] + v[199] * v[672]
                + 2e0 * v[197] * v[673] + v[229] * v[655] * v[7] + v[230] * v[656] * v[8];
            /* 677 = \[OverBracket]_Fref_3|1(Ft|Rgm)_3|1 */
            v[677] = v[236] * v[674];
            /* 678 = \[OverBracket]_Fref_2|1(Ft|Rgm)_2|1 */
            v[678] = v[236] * v[675];
            /* 679 = \[OverBracket]_Fref_1|1(Ft|Rgm)_1|1 */
            v[679] = v[236] * v[676];
            /* 680 = \[OverBracket]_R_3|3;\[Xi](Fref|Rgm)_3|1 */
            v[680] = v[249] * v[677];
            /* 681 = \[OverBracket]_R_3|2;\[Xi](Fref|Rgm)_3|1 */
            v[681] = v[251] * v[677];
            /* 682 = \[OverBracket]_R_3|1;\[Xi](Fref|Rgm)_3|1 */
            v[682] = v[253] * v[677];
            /* 683 = \[OverBracket]_R_2|3;\[Xi](Fref|Rgm)_2|1 */
            v[683] = v[249] * v[678];
            /* 684 = \[OverBracket]_R_2|2;\[Xi](Fref|Rgm)_2|1 */
            v[684] = v[251] * v[678];
            /* 685 = \[OverBracket]_R_2|1;\[Xi](Fref|Rgm)_2|1 */
            v[685] = v[253] * v[678];
            /* 686 = \[OverBracket]_R_1|3;\[Xi](Fref|Rgm)_1|1 */
            v[686] = v[249] * v[679];
            /* 687 = \[OverBracket]_R_1|2;\[Xi](Fref|Rgm)_1|1 */
            v[687] = v[251] * v[679];
            /* 688 = \[OverBracket]_R_1|1;\[Xi](Fref|Rgm)_1|1 */
            v[688] = v[253] * v[679];
            /* 689 = \[OverBracket]_Ft_1|2(\[Yen]|Rgm)_239 */
            v[689] = v[199] * v[667] + v[197] * v[671] + (v[10] * v[247] * v[676]) / 2e0 + v[229] * v[668] * v[7]
                + v[239] * v[656] * v[9];
            /* 690 = \[OverBracket]_Ft_1|3(\[Yen]|Rgm)_239 */
            v[690] = v[198] * v[667] + v[197] * v[672] + (v[11] * v[245] * v[676]) / 2e0 + v[230] * v[668] * v[8]
                + v[239] * v[655] * v[9];
            /* 691 = \[OverBracket]_Ft_2|2(\[Yen]|Rgm)_239 */
            v[691] = v[211] * v[667] + v[209] * v[671] + (v[10] * v[247] * v[675]) / 2e0 + v[229] * v[669] * v[7]
                + v[239] * v[647] * v[9];
            /* 692 = \[OverBracket]_Ft_2|3(\[Yen]|Rgm)_239 */
            v[692] = v[210] * v[667] + v[209] * v[672] + (v[11] * v[245] * v[675]) / 2e0 + v[230] * v[669] * v[8]
                + v[239] * v[646] * v[9];
            /* 693 = \[OverBracket]_Ft_3|2(\[Yen]|Rgm)_239 */
            v[693] = v[214] * v[667] + v[212] * v[671] + (v[10] * v[247] * v[674]) / 2e0 + v[229] * v[670] * v[7]
                + v[239] * v[638] * v[9];
            /* 694 = \[OverBracket]_Ft_3|3(\[Yen]|Rgm)_239 */
            v[694] = v[213] * v[667] + v[212] * v[672] + (v[11] * v[245] * v[674]) / 2e0 + v[230] * v[670] * v[8]
                + v[239] * v[637] * v[9];
            /* 695 = \[OverBracket]_R_1|1(Ft|Rgm)_1|3 */
            v[695] = v[269] * v[679] + v[67] * v[689] + v[690] * v[70];
            /* 696 = \[OverBracket]_R_1|2(Ft|Rgm)_1|3 */
            v[696] = v[271] * v[679] + v[68] * v[689] + v[690] * v[71];
            /* 697 = \[OverBracket]_R_1|3(Ft|Rgm)_1|3 */
            v[697] = v[273] * v[679] + v[689] * v[69] + v[690] * v[72];
            /* 698 = \[OverBracket]_R_2|1(Ft|Rgm)_2|3 */
            v[698] = v[269] * v[678] + v[67] * v[691] + v[692] * v[70];
            /* 699 = \[OverBracket]_R_2|2(Ft|Rgm)_2|3 */
            v[699] = v[271] * v[678] + v[68] * v[691] + v[692] * v[71];
            /* 700 = \[OverBracket]_R_2|3(Ft|Rgm)_2|3 */
            v[700] = v[273] * v[678] + v[69] * v[691] + v[692] * v[72];
            /* 701 = \[OverBracket]_R_3|1(Ft|Rgm)_3|3 */
            v[701] = v[269] * v[677] + v[67] * v[693] + v[694] * v[70];
            /* 702 = \[OverBracket]_R_3|2(Ft|Rgm)_3|3 */
            v[702] = v[271] * v[677] + v[68] * v[693] + v[694] * v[71];
            /* 703 = \[OverBracket]_R_3|3(Ft|Rgm)_3|3 */
            v[703] = v[273] * v[677] + v[69] * v[693] + v[694] * v[72];
            b704 = b105;
            v[704] = b704;
            if (b704) {
                v[727] = v[687] / 720e0;
                v[726] = v[685] / 720e0;
                v[723] = v[698] / 720e0;
                v[722] = v[696] / 720e0;
                v[719] = v[686] / 720e0;
                v[718] = v[682] / 720e0;
                v[715] = v[701] / 720e0;
                v[714] = v[697] / 720e0;
                v[711] = v[683] / 720e0;
                v[710] = v[681] / 720e0;
                v[707] = v[702] / 720e0;
                v[706] = v[700] / 720e0;
                /* 501 = \[OverBracket]_\[Yen]_139(v|Rgm)_3|3;\[Xi] */
                v[501] = v[501] - (v[162] * v[680]) / 720e0 - (v[110] * v[703]) / 720e0;
                /* 502 = \[OverBracket]_\[Yen]_149(v|Rgm)_3|3;\[Xi] */
                v[502] = v[502] - (v[110] * v[680]) / 720e0;
                /* 705 = \[OverBracket]_\[Yen]_123(v|Rgm)_2|3 */
                v[705] = v[706] + v[707];
                /* 708 = \[OverBracket]_\[Yen]_124(v|Rgm)_2|3 */
                v[708] = v[706] - v[707];
                /* 709 = \[OverBracket]_\[Yen]_163(v|Rgm)_2|3;\[Xi] */
                v[709] = v[710] + v[711];
                /* 712 = \[OverBracket]_\[Yen]_157(v|Rgm)_2|3;\[Xi] */
                v[712] = -v[710] + v[711];
                /* 498 = \[OverBracket]_\[Yen]_134(v|Rgm)_2|2;\[Xi] */
                v[498] = v[498] - (v[162] * v[684]) / 720e0 - (v[110] * v[699]) / 720e0;
                /* 499 = \[OverBracket]_\[Yen]_152(v|Rgm)_2|2;\[Xi] */
                v[499] = v[499] - (v[110] * v[684]) / 720e0;
                /* 713 = \[OverBracket]_\[Yen]_120(v|Rgm)_1|3 */
                v[713] = v[714] + v[715];
                /* 716 = \[OverBracket]_\[Yen]_121(v|Rgm)_1|3 */
                v[716] = v[714] - v[715];
                /* 717 = \[OverBracket]_\[Yen]_164(v|Rgm)_1|3;\[Xi] */
                v[717] = v[718] + v[719];
                /* 720 = \[OverBracket]_\[Yen]_158(v|Rgm)_1|3;\[Xi] */
                v[720] = -v[718] + v[719];
                /* 721 = \[OverBracket]_\[Yen]_114(v|Rgm)_1|2 */
                v[721] = v[722] + v[723];
                /* 724 = \[OverBracket]_\[Yen]_115(v|Rgm)_1|2 */
                v[724] = v[722] - v[723];
                /* 725 = \[OverBracket]_\[Yen]_165(v|Rgm)_1|2;\[Xi] */
                v[725] = v[726] + v[727];
                /* 728 = \[OverBracket]_\[Yen]_159(v|Rgm)_1|2;\[Xi] */
                v[728] = -v[726] + v[727];
                /* 497 = \[OverBracket]_\[Yen]_127(v|Rgm)_1|1;\[Xi] */
                v[497] = v[497] - (v[162] * v[688]) / 720e0 - (v[110] * v[695]) / 720e0;
                /* 500 = \[OverBracket]_\[Yen]_153(v|Rgm)_1|1;\[Xi] */
                v[500] = v[500] - (v[110] * v[688]) / 720e0;
                /* 427 = \[OverBracket]_\[Yen]_110(\[Yen]|Rgm)_165 */
                v[427] = v[427] - (v[149] * v[680]) / 720e0 - (v[152] * v[684]) / 720e0 - (v[153] * v[688]) / 720e0 - (v[127] * v[695])
                    / 720e0 - (v[134] * v[699]) / 720e0 - (v[139] * v[703]) / 720e0 + v[102] * v[103] * v[705] + v[494] * v[709]
                    + v[101] * v[103] * v[713] + v[491] * v[717] + v[101] * v[102] * v[721] + v[488] * v[725];
                /* 435 = \[OverBracket]_\[Yen]_162(\[Yen]|Rgm)_165 */
                v[435] = v[435] - (v[139] * v[680]) / 720e0 - (v[134] * v[684]) / 720e0 - (v[127] * v[688]) / 720e0
                    + v[102] * v[103] * v[709] + v[101] * v[103] * v[717] + v[101] * v[102] * v[725];
                /* 434 = \[OverBracket]_\[Yen]_160(\[Yen]|Rgm)_110 */
                v[434] = v[104] * v[427] + v[434];
                /* 448 = \[OverBracket]_ff_(\[Yen]|Rgm)_110 */
                v[448] = v[160] * v[427] + v[448];
                /* 402 = \[OverBracket]_u\[Phi]_4(\[Yen]|Rgm)_157 */
                v[402] = v[402] - 6e0 * v[112] * v[708] - 6e0 * v[156] * v[712] + v[103] * v[110] * v[713] + v[313] * v[717]
                    + v[102] * v[110] * v[721] + v[318] * v[725];
                /* 403 = \[OverBracket]_u\[Phi]_4;\[Xi](\[Yen]|Rgm)_157 */
                v[403] = v[403] - 6e0 * v[112] * v[712] + v[103] * v[110] * v[717] + v[102] * v[110] * v[725];
                /* 409 = \[OverBracket]_u\[Phi]_5(\[Yen]|Rgm)_158 */
                v[409] = v[409] + v[103] * v[110] * v[705] + v[313] * v[709] + 6e0 * v[112] * v[716] + 6e0 * v[156] * v[720]
                    + v[101] * v[110] * v[721] + v[319] * v[725];
                /* 410 = \[OverBracket]_u\[Phi]_5;\[Xi](\[Yen]|Rgm)_158 */
                v[410] = v[410] + v[103] * v[110] * v[709] + 6e0 * v[112] * v[720] + v[101] * v[110] * v[725];
                /* 428 = \[OverBracket]_\[Yen]_112(\[Yen]|Rgm)_159 */
                v[428] = v[428] - 6e0 * v[101] * v[708] - 6e0 * v[145] * v[712] + 6e0 * v[102] * v[716] + 6e0 * v[147] * v[720]
                    - 6e0 * v[103] * v[724] - 6e0 * v[150] * v[728];
                /* 433 = \[OverBracket]_\[Yen]_156(\[Yen]|Rgm)_159 */
                v[433] = v[433] - 6e0 * v[101] * v[712] + 6e0 * v[102] * v[720] - 6e0 * v[103] * v[728];
                /* 418 = \[OverBracket]_u\[Phi]_6(\[Yen]|Rgm)_159 */
                v[418] = v[418] + v[102] * v[110] * v[705] + v[318] * v[709] + v[101] * v[110] * v[713] + v[319] * v[717]
                    - 6e0 * v[112] * v[724] - 6e0 * v[156] * v[728];
                /* 419 = \[OverBracket]_u\[Phi]_6;\[Xi](\[Yen]|Rgm)_159 */
                v[419] = v[419] + v[102] * v[110] * v[709] + v[101] * v[110] * v[717] - 6e0 * v[112] * v[728];
                /* 432 = \[OverBracket]_\[Yen]_155(\[Yen]|Rgm)_112 */
                v[432] = v[104] * v[428] + v[432];
                /* 448 = \[OverBracket]_ff_(\[Yen]|Rgm)_112 */
                v[448] = v[155] * v[428] + v[448];
                /* 432 = \[OverBracket]_\[Yen]_155(\[Yen]|Rgm)_156 */
                v[432] = v[432] + v[154] * v[433];
                /* 729 = \[OverBracket]_\[Yen]_161(\[Yen]|Rgm)_156 */
                v[729] = v[433];
                /* 449 = \[OverBracket]_ff_;\[Xi](\[Yen]|Rgm)_156 */
                v[449] = v[155] * v[433] + v[449];
                /* 448 = \[OverBracket]_ff_(\[Yen]|Rgm)_155 */
                v[448] = v[432] + v[448];
                /* 434 = \[OverBracket]_\[Yen]_160(\[Yen]|Rgm)_162 */
                v[434] = v[434] + v[154] * v[435];
                /* 730 = \[OverBracket]_\[Yen]_161(\[Yen]|Rgm)_162 */
                v[730] = v[435] + v[729];
                /* 449 = \[OverBracket]_ff_;\[Xi](\[Yen]|Rgm)_162 */
                v[449] = v[160] * v[435] + v[449];
                /* 448 = \[OverBracket]_ff_(\[Yen]|Rgm)_160 */
                v[448] = v[434] + v[448];
                /* 448 = \[OverBracket]_ff_(\[Yen]|Rgm)_161 */
                v[448] = v[448] + v[154] * v[730];
                /* 449 = \[OverBracket]_ff_;\[Xi](\[Yen]|Rgm)_161 */
                v[449] = v[449] + v[104] * v[730];
            }
            else {
                /* 501 = \[OverBracket]_\[Yen]_139(v|Rgm)_3|3;\[Xi] */
                v[501] = v[329] * v[680] + (v[129] * v[703]) / v[104] + v[732];
                /* 502 = \[OverBracket]_\[Yen]_149(v|Rgm)_3|3;\[Xi] */
                v[502] = (v[129] * v[680]) / v[104] + v[734];
                /* 735 = \[OverBracket]_\[Yen]_138(v|Rgm)_2|3 */
                v[735] = v[700] + v[702];
                /* 736 = \[OverBracket]_\[Yen]_137(v|Rgm)_2|3 */
                v[736] = v[700] - v[702];
                /* 737 = \[OverBracket]_\[Yen]_184(v|Rgm)_2|3;\[Xi] */
                v[737] = v[681] + v[683];
                /* 738 = \[OverBracket]_\[Yen]_180(v|Rgm)_2|3;\[Xi] */
                v[738] = -v[681] + v[683];
                /* 498 = \[OverBracket]_\[Yen]_134(v|Rgm)_2|2;\[Xi] */
                v[498] = v[329] * v[684] + (v[129] * v[699]) / v[104] + v[739];
                /* 499 = \[OverBracket]_\[Yen]_152(v|Rgm)_2|2;\[Xi] */
                v[499] = (v[129] * v[684]) / v[104] + v[740];
                /* 741 = \[OverBracket]_\[Yen]_136(v|Rgm)_1|3 */
                v[741] = v[697] + v[701];
                /* 742 = \[OverBracket]_\[Yen]_135(v|Rgm)_1|3 */
                v[742] = v[697] - v[701];
                /* 743 = \[OverBracket]_\[Yen]_186(v|Rgm)_1|3;\[Xi] */
                v[743] = v[682] + v[686];
                /* 744 = \[OverBracket]_\[Yen]_181(v|Rgm)_1|3;\[Xi] */
                v[744] = -v[682] + v[686];
                /* 745 = \[OverBracket]_\[Yen]_133(v|Rgm)_1|2 */
                v[745] = v[696] + v[698];
                /* 746 = \[OverBracket]_\[Yen]_132(v|Rgm)_1|2 */
                v[746] = v[696] - v[698];
                /* 747 = \[OverBracket]_\[Yen]_187(v|Rgm)_1|2;\[Xi] */
                v[747] = v[685] + v[687];
                /* 748 = \[OverBracket]_\[Yen]_182(v|Rgm)_1|2;\[Xi] */
                v[748] = -v[685] + v[687];
                /* 497 = \[OverBracket]_\[Yen]_127(v|Rgm)_1|1;\[Xi] */
                v[497] = v[329] * v[688] + (v[129] * v[695]) / v[104] + v[749];
                /* 500 = \[OverBracket]_\[Yen]_153(v|Rgm)_1|1;\[Xi] */
                v[500] = (v[129] * v[688]) / v[104] + v[750];
                /* 430 = \[OverBracket]_\[Yen]_129(\[Yen]|Rgm)_187 */
                v[430] = v[430] + v[620] * v[680] + v[618] * v[684] + v[617] * v[688] + (v[127] * v[695]) / v[104] + (v[134] * v[699])
                    / v[104] + (v[139] * v[703]) / v[104] - (v[102] * v[103] * v[735]) / v[104] + v[546] * v[737] - (v[101] * v[103] * v[741]
                        ) / v[104] + v[547] * v[743] - (v[101] * v[102] * v[745]) / v[104] + v[571] * v[747];
                /* 441 = \[OverBracket]_\[Yen]_183(\[Yen]|Rgm)_187 */
                v[441] = v[441] + (v[139] * v[680]) / v[104] + (v[134] * v[684]) / v[104] + (v[127] * v[688]) / v[104] -
                    (v[102] * v[103] * v[737]) / v[104] - (v[101] * v[103] * v[743]) / v[104] - (v[101] * v[102] * v[747]) / v[104];
                /* 442 = \[OverBracket]_\[Yen]_185(\[Yen]|Rgm)_187 */
                v[442] = v[442] - v[129] * v[139] * v[154] * v[680] - v[129] * v[134] * v[154] * v[684]
                    - v[127] * v[129] * v[154] * v[688] + v[102] * v[103] * v[129] * v[154] * v[737]
                    + v[101] * v[103] * v[129] * v[154] * v[743] + v[101] * v[102] * v[129] * v[154] * v[747];
                /* 439 = \[OverBracket]_\[Yen]_178(\[Yen]|Rgm)_129 */
                v[439] = v[430] + v[439];
                /* 402 = \[OverBracket]_u\[Phi]_4(\[Yen]|Rgm)_180 */
                v[402] = -(v[130] * v[175] * v[736]) + v[340] * v[738] - (v[103] * v[129] * v[741]) / v[104] + v[339] * v[743] -
                    (v[102] * v[129] * v[745]) / v[104] + v[343] * v[747] + v[752];
                /* 403 = \[OverBracket]_u\[Phi]_4;\[Xi](\[Yen]|Rgm)_180 */
                v[403] = -(v[130] * v[175] * v[738]) - (v[103] * v[129] * v[743]) / v[104] - (v[102] * v[129] * v[747]) / v[104]
                    + v[756];
                /* 409 = \[OverBracket]_u\[Phi]_5(\[Yen]|Rgm)_181 */
                v[409] = -((v[103] * v[129] * v[735]) / v[104]) + v[339] * v[737] + v[130] * v[175] * v[742] - v[340] * v[744] -
                    (v[101] * v[129] * v[745]) / v[104] + v[344] * v[747] + v[751];
                /* 410 = \[OverBracket]_u\[Phi]_5;\[Xi](\[Yen]|Rgm)_181 */
                v[410] = -((v[103] * v[129] * v[737]) / v[104]) + v[130] * v[175] * v[744] - (v[101] * v[129] * v[747]) / v[104]
                    + v[754];
                /* 438 = \[OverBracket]_\[Yen]_177(\[Yen]|Rgm)_182 */
                v[438] = v[438] - v[101] * v[130] * v[738] + v[102] * v[130] * v[744] - v[103] * v[130] * v[748];
                /* 440 = \[OverBracket]_\[Yen]_179(\[Yen]|Rgm)_182 */
                v[440] = v[440] - v[101] * v[175] * v[738] + v[102] * v[175] * v[744] - v[103] * v[175] * v[748];
                /* 418 = \[OverBracket]_u\[Phi]_6(\[Yen]|Rgm)_182 */
                v[418] = -((v[102] * v[129] * v[735]) / v[104]) + v[343] * v[737] - (v[101] * v[129] * v[741]) / v[104]
                    + v[344] * v[743] - v[130] * v[175] * v[746] + v[340] * v[748] + v[753];
                /* 419 = \[OverBracket]_u\[Phi]_6;\[Xi](\[Yen]|Rgm)_182 */
                v[419] = -((v[102] * v[129] * v[737]) / v[104]) - (v[101] * v[129] * v[743]) / v[104] - v[130] * v[175] * v[748]
                    + v[755];
                /* 431 = \[OverBracket]_\[Yen]_130(\[Yen]|Rgm)_183 */
                v[431] = v[431] - v[176] * v[441] - v[101] * v[175] * v[736] + v[102] * v[175] * v[742] - v[103] * v[175] * v[746]
                    + v[738] * v[757] + v[744] * v[758] + v[748] * v[759];
                /* 437 = \[OverBracket]_\[Yen]_176(\[Yen]|Rgm)_183 */
                v[437] = v[437] - v[130] * v[441];
                /* 429 = \[OverBracket]_\[Yen]_128(\[Yen]|Rgm)_130 */
                v[429] = v[429] + v[178] * v[431];
                /* 443 = \[OverBracket]_\[Yen]_345(\[Yen]|Rgm)_177 */
                v[443] = -(v[176] * v[438]) + v[443];
                /* 437 = \[OverBracket]_\[Yen]_176(\[Yen]|Rgm)_177 */
                v[437] = v[437] - v[345] * v[438];
                /* 439 = \[OverBracket]_\[Yen]_178(\[Yen]|Rgm)_179 */
                v[439] = v[439] + v[176] * v[440];
                /* 437 = \[OverBracket]_\[Yen]_176(\[Yen]|Rgm)_179 */
                v[437] = v[437] + v[178] * v[440];
                /* 429 = \[OverBracket]_\[Yen]_128(\[Yen]|Rgm)_178 */
                v[429] = v[429] - v[130] * v[439];
                /* 429 = \[OverBracket]_\[Yen]_128(\[Yen]|Rgm)_345 */
                v[429] = v[429] - 2e0 * v[443] * v[538];
                /* 448 = \[OverBracket]_ff_(\[Yen]|Rgm)_128 */
                v[448] = (v[175] * v[429]) / 2e0 + (-(v[129] * v[149] * v[185]) - v[139] * v[183] * v[185]) * v[680] + (-
                    (v[129] * v[152] * v[185]) - v[134] * v[183] * v[185]) * v[684] + (-(v[129] * v[153] * v[185])
                        - v[127] * v[183] * v[185]) * v[688] - v[127] * v[129] * v[185] * v[695] - v[129] * v[134] * v[185] * v[699]
                    - v[129] * v[139] * v[185] * v[703] + v[731] + v[102] * v[103] * v[129] * v[185] * v[735]
                    + v[101] * v[103] * v[129] * v[185] * v[741] + v[101] * v[102] * v[129] * v[185] * v[745] + v[737] *
                    (v[103] * v[129] * v[147] * v[185] + v[760] + v[761]) + v[743] * (v[762] + v[763] + v[764]) + v[747] * (v[765] + v[766]
                        + v[767]);
                /* 436 = \[OverBracket]_\[Yen]_175(\[Yen]|Rgm)_176 */
                v[436] = v[436] + (v[154] * v[437]) / 2e0 - v[101] * v[130] * v[736] + v[102] * v[130] * v[742] - v[103] * v[130] * v[746]
                    + v[738] * v[768] + v[744] * v[769] + v[748] * v[770];
                /* 449 = \[OverBracket]_ff_;\[Xi](\[Yen]|Rgm)_176 */
                v[449] = (v[175] * v[437]) / 2e0 - v[129] * v[139] * v[185] * v[680] - v[129] * v[134] * v[185] * v[684]
                    - v[127] * v[129] * v[185] * v[688] + v[733] + v[102] * v[103] * v[129] * v[185] * v[737]
                    + v[101] * v[103] * v[129] * v[185] * v[743] + v[101] * v[102] * v[129] * v[185] * v[747];
                /* 448 = \[OverBracket]_ff_(\[Yen]|Rgm)_175 */
                v[448] = -0.5e0 * (v[175] * v[436]) / v[104] + v[448];
                /* 448 = \[OverBracket]_ff_(\[Yen]|Rgm)_185 */
                v[448] = v[448] - 2e0 * v[442] * v[544];
            };
            /* 771 = \[OverBracket]_\[Yen]_108(ff|Rgm) */
            v[771] = v[448];
            /* 772 = \[OverBracket]_\[Yen]_107(ff|Rgm) */
            v[772] = v[448];
            /* 773 = \[OverBracket]_\[Yen]_117(ff|Rgm) */
            v[773] = v[448];
            /* 774 = \[OverBracket]_\[Yen]_108(\[Yen]|Rgm)_127 */
            v[774] = v[497] + v[771];
            /* 775 = \[OverBracket]_\[Yen]_107(\[Yen]|Rgm)_127 */
            v[775] = v[774];
            /* 776 = \[OverBracket]_\[Yen]_108(\[Yen]|Rgm)_134 */
            v[776] = v[498] + v[774];
            /* 777 = \[OverBracket]_\[Yen]_117(\[Yen]|Rgm)_134 */
            v[777] = v[498] + v[773];
            /* 778 = \[OverBracket]_\[Yen]_151(\[Yen]|Rgm)_152 */
            v[778] = v[499];
            /* 779 = \[OverBracket]_\[Yen]_146(\[Yen]|Rgm)_152 */
            v[779] = v[499];
            /* 780 = \[OverBracket]_\[Yen]_151(\[Yen]|Rgm)_153 */
            v[780] = v[500] + v[778];
            /* 781 = \[OverBracket]_\[Yen]_148(\[Yen]|Rgm)_153 */
            v[781] = v[500];
            /* 782 = \[OverBracket]_\[Yen]_151(ff|Rgm)_;\[Xi] */
            v[782] = v[449] + v[780];
            /* 783 = \[OverBracket]_\[Yen]_148(ff|Rgm)_;\[Xi] */
            v[783] = v[449] + v[781];
            /* 784 = \[OverBracket]_\[Yen]_146(ff|Rgm)_;\[Xi] */
            v[784] = v[449] + v[779];
            /* 418 = \[OverBracket]_u\[Phi]_6(\[Yen]|Rgm)_151 */
            v[418] = v[418] + 2e0 * v[103] * v[776] + 2e0 * v[150] * v[782];
            /* 419 = \[OverBracket]_u\[Phi]_6;\[Xi](\[Yen]|Rgm)_151 */
            v[419] = v[419] + 2e0 * v[103] * v[782];
            /* 785 = \[OverBracket]_peIO_3|6(u\[Phi]|Rgm)_6 */
            v[785] = v[418] * v[52];
            /* 786 = \[OverBracket]_peIO_2|6(u\[Phi]|Rgm)_6 */
            v[786] = v[418] * v[50];
            /* 787 = \[OverBracket]_peIO_1|6(u\[Phi]|Rgm)_6 */
            v[787] = v[418] * v[48];
            /* 788 = \[OverBracket]_\[Yen]_107(\[Yen]|Rgm)_139 */
            v[788] = v[501] + v[775];
            /* 789 = \[OverBracket]_\[Yen]_117(\[Yen]|Rgm)_139 */
            v[789] = v[501] + v[777];
            /* 790 = \[OverBracket]_\[Yen]_148(\[Yen]|Rgm)_149 */
            v[790] = v[502] + v[783];
            /* 791 = \[OverBracket]_\[Yen]_146(\[Yen]|Rgm)_149 */
            v[791] = v[502] + v[784];
            /* 409 = \[OverBracket]_u\[Phi]_5(\[Yen]|Rgm)_148 */
            v[409] = v[409] + 2e0 * v[102] * v[788] + 2e0 * v[147] * v[790];
            /* 410 = \[OverBracket]_u\[Phi]_5;\[Xi](\[Yen]|Rgm)_148 */
            v[410] = v[410] + 2e0 * v[102] * v[790];
            /* 792 = \[OverBracket]_peIO_3|5(u\[Phi]|Rgm)_5 */
            v[792] = v[409] * v[52];
            /* 793 = \[OverBracket]_peIO_2|5(u\[Phi]|Rgm)_5 */
            v[793] = v[409] * v[50];
            /* 794 = \[OverBracket]_peIO_1|5(u\[Phi]|Rgm)_5 */
            v[794] = v[409] * v[48];
            /* 402 = \[OverBracket]_u\[Phi]_4(\[Yen]|Rgm)_146 */
            v[402] = v[402] + 2e0 * v[101] * v[789] + 2e0 * v[145] * v[791];
            /* 403 = \[OverBracket]_u\[Phi]_4;\[Xi](\[Yen]|Rgm)_146 */
            v[403] = v[403] + 2e0 * v[101] * v[791];
            /* 795 = \[OverBracket]_peIO_3|4(u\[Phi]|Rgm)_4 */
            v[795] = v[402] * v[52];
            /* 796 = \[OverBracket]_peIO_2|4(u\[Phi]|Rgm)_4 */
            v[796] = v[402] * v[50];
            /* 797 = \[OverBracket]_peIO_1|4(u\[Phi]|Rgm)_4 */
            v[797] = v[402] * v[48];
            /* 798 = \[OverBracket]_peIO_3|4(u\[Phi]|Rgm)_4;\[Xi] */
            v[798] = v[403] * v[59] + v[795];
            /* 799 = \[OverBracket]_peIO_2|4(u\[Phi]|Rgm)_4;\[Xi] */
            v[799] = v[403] * v[57] + v[796];
            /* 800 = \[OverBracket]_peIO_1|4(u\[Phi]|Rgm)_4;\[Xi] */
            v[800] = v[403] * v[56] + v[797];
            /* 801 = \[OverBracket]_peIO_3|5(u\[Phi]|Rgm)_5;\[Xi] */
            v[801] = v[410] * v[59] + v[792];
            /* 802 = \[OverBracket]_peIO_2|5(u\[Phi]|Rgm)_5;\[Xi] */
            v[802] = v[410] * v[57] + v[793];
            /* 803 = \[OverBracket]_peIO_1|5(u\[Phi]|Rgm)_5;\[Xi] */
            v[803] = v[410] * v[56] + v[794];
            /* 804 = \[OverBracket]_peIO_3|6(u\[Phi]|Rgm)_6;\[Xi] */
            v[804] = v[419] * v[59] + v[785];
            /* 805 = \[OverBracket]_peIO_2|6(u\[Phi]|Rgm)_6;\[Xi] */
            v[805] = v[419] * v[57] + v[786];
            /* 806 = \[OverBracket]_peIO_1|6(u\[Phi]|Rgm)_6;\[Xi] */
            v[806] = v[419] * v[56] + v[787];
            v[1408] = v[56] * v[679];
            v[1409] = v[56] * v[678];
            v[1410] = v[56] * v[677];
            v[1411] = v[800];
            v[1412] = v[803];
            v[1413] = v[806];
            v[1414] = v[57] * v[679];
            v[1415] = v[57] * v[678];
            v[1416] = v[57] * v[677];
            v[1417] = v[799];
            v[1418] = v[802];
            v[1419] = v[805];
            v[1420] = v[59] * v[679];
            v[1421] = v[59] * v[678];
            v[1422] = v[59] * v[677];
            v[1423] = v[798];
            v[1424] = v[801];
            v[1425] = v[804];
            v[807] = 0e0;/*debug*/
            R[i227 - 1] += v[384] * v[47];
            v[385] = 0e0;/*debug*/
            for (ii386 = i227; ii386 <= 18; ii386++) {
                v[386] = ii386;
                /* 388 = \[DoubleStruckCapitalG]_n */
                v[388] = v[1148 + ii386];
                /* 808 = Kgmn */
                v[808] = v[1407 + ii386];
                T[i227 - 1][ii386 - 1] += v[47] * v[808];
                v[809] = 0e0;/*debug*/
            };/* end for */
        };/* end for */
    };/* end for */

    for (int i = 0; i < 18; i++)
    {
        for (int j = 0; j < 18; j++)
            if (j != i)
                T[j][i] = T[i][j];
    }
};


/******************* S U B R O U T I N E *********************/
void NLEBBE3D::RKt(double D[7], double X[3][3], double U[3][6]
    , double **T, double* R)
{
    double v[1518];

    int i36, i218, i379, b98, b99, b271, b437, b697;
    /* 1 = Em */
    v[1] = D[0];
    /* 2 = \[Nu] */
    v[2] = D[1];
    v[215] = 1e0 / (1e0 + v[2]);
    /* 216 = \[Mu] */
    v[216] = (v[1] * v[215]) / 2e0;
    /* 214 = \[Lambda] */
    v[214] = (v[1] * v[2] * v[215]) / (1e0 - 2e0 * v[2]);
    v[224] = 3e0 * v[214] + 2e0 * v[216];
    v[223] = 1e0 / (v[214] + v[216]);
    /* 3 = Bp */
    v[3] = D[2];
    /* 4 = Hp */
    v[4] = D[3];
    /* 5 = Zx */
    v[5] = D[4];
    /* 6 = Zy */
    v[6] = D[5];
    /* 7 = Zz */
    v[7] = D[6];
    /* 8 = XIO_1|1 */
    v[8] = X[0][0];
    /* 9 = XIO_1|2 */
    v[9] = X[0][1];
    /* 10 = XIO_1|3 */
    v[10] = X[0][2];
    /* 11 = XIO_2|1 */
    v[11] = X[1][0];
    /* 12 = XIO_2|2 */
    v[12] = X[1][1];
    /* 13 = XIO_2|3 */
    v[13] = X[1][2];
    /* 14 = XIO_3|1 */
    v[14] = X[2][0];
    /* 68 = X0_1;\[Xi];\[Xi] */
    v[68] = -2e0 * v[11] + v[14] + v[8];
    /* 15 = XIO_3|2 */
    v[15] = X[2][1];
    /* 67 = X0_2;\[Xi];\[Xi] */
    v[67] = -2e0 * v[12] + v[15] + v[9];
    /* 16 = XIO_3|3 */
    v[16] = X[2][2];
    /* 66 = X0_3;\[Xi];\[Xi] */
    v[66] = v[10] - 2e0 * v[13] + v[16];
    /* 17 = peIO_1|1 */
    v[17] = U[0][0];
    /* 18 = peIO_1|2 */
    v[18] = U[0][1];
    /* 19 = peIO_1|3 */
    v[19] = U[0][2];
    /* 20 = peIO_1|4 */
    v[20] = U[0][3];
    /* 21 = peIO_1|5 */
    v[21] = U[0][4];
    /* 22 = peIO_1|6 */
    v[22] = U[0][5];
    /* 23 = peIO_2|1 */
    v[23] = U[1][0];
    /* 24 = peIO_2|2 */
    v[24] = U[1][1];
    /* 25 = peIO_2|3 */
    v[25] = U[1][2];
    /* 26 = peIO_2|4 */
    v[26] = U[1][3];
    /* 27 = peIO_2|5 */
    v[27] = U[1][4];
    /* 28 = peIO_2|6 */
    v[28] = U[1][5];
    /* 29 = peIO_3|1 */
    v[29] = U[2][0];
    /* 30 = peIO_3|2 */
    v[30] = U[2][1];
    /* 31 = peIO_3|3 */
    v[31] = U[2][2];
    /* 32 = peIO_3|4 */
    v[32] = U[2][3];
    /* 33 = peIO_3|5 */
    v[33] = U[2][4];
    /* 34 = peIO_3|6 */
    v[34] = U[2][5];
    v[1142] = v[17];
    v[1143] = v[18];
    v[1144] = v[19];
    v[1145] = v[20];
    v[1146] = v[21];
    v[1147] = v[22];
    v[1148] = v[23];
    v[1149] = v[24];
    v[1150] = v[25];
    v[1151] = v[26];
    v[1152] = v[27];
    v[1153] = v[28];
    v[1154] = v[29];
    v[1155] = v[30];
    v[1156] = v[31];
    v[1157] = v[32];
    v[1158] = v[33];
    v[1159] = v[34];
    v[380] = 0e0;/*debug*/
    /* 35 = Le */
    v[35] = sqrt(Power(v[10] - v[13], 2) + Power(-v[11] + v[8], 2) + Power(-v[12] + v[9], 2));
    for (i36 = 1; i36 <= 8; i36++) {
        v[36] = i36;
        /* 37 = \[Xi] */
        v[37] = GaussIntegrationPoints(i36,1);
        v[51] = v[37] / 2e0;
        v[44] = 1e0 + v[37];
        /* 52 = Nh_3;\[Xi] */
        v[52] = v[44] / 2e0 + v[51];
        v[42] = -1e0 + v[37];
        /* 50 = Nh_2;\[Xi] */
        v[50] = -v[42] - v[44];
        /* 49 = Nh_1;\[Xi] */
        v[49] = v[42] / 2e0 + v[51];
        /* 143 = u\[Phi]_6;\[Xi] */
        v[143] = v[22] * v[49] + v[28] * v[50] + v[34] * v[52];
        /* 140 = u\[Phi]_5;\[Xi] */
        v[140] = v[21] * v[49] + v[27] * v[50] + v[33] * v[52];
        /* 138 = u\[Phi]_4;\[Xi] */
        v[138] = v[20] * v[49] + v[26] * v[50] + v[32] * v[52];
        /* 137 = u\[Phi]_3;\[Xi] */
        v[137] = v[19] * v[49] + v[25] * v[50] + v[31] * v[52];
        /* 136 = u\[Phi]_2;\[Xi] */
        v[136] = v[18] * v[49] + v[24] * v[50] + v[30] * v[52];
        /* 135 = u\[Phi]_1;\[Xi] */
        v[135] = v[17] * v[49] + v[23] * v[50] + v[29] * v[52];
        /* 55 = X0_3;\[Xi] */
        v[55] = v[10] * v[49] + v[13] * v[50] + v[16] * v[52];
        /* 54 = X0_2;\[Xi] */
        v[54] = v[12] * v[50] + v[15] * v[52] + v[49] * v[9];
        /* 53 = X0_1;\[Xi] */
        v[53] = v[11] * v[50] + v[14] * v[52] + v[49] * v[8];
        v[133] = 2e0 * v[55] * v[66] + 2e0 * v[54] * v[67] + 2e0 * v[53] * v[68];
        v[70] = (v[53] * v[53]) + (v[54] * v[54]) + (v[55] * v[55]);
        v[134] = -0.5e0 * v[133] / (v[70] * sqrt(v[70]));
        v[69] = 1e0 / sqrt(v[70]);
        v[71] = v[134];
        v[57] = v[69];
        /* 74 = t1_3;\[Xi] */
        v[74] = v[57] * v[66] + v[55] * v[71];
        /* 73 = t1_2;\[Xi] */
        v[73] = v[57] * v[67] + v[54] * v[71];
        /* 75 = t2_1;\[Xi] */
        v[75] = -(v[7] * v[73]) + v[6] * v[74];
        /* 72 = t1_1;\[Xi] */
        v[72] = v[57] * v[68] + v[53] * v[71];
        /* 77 = t2_3;\[Xi] */
        v[77] = -(v[6] * v[72]) + v[5] * v[73];
        /* 76 = t2_2;\[Xi] */
        v[76] = v[7] * v[72] - v[5] * v[74];
        /* 38 = \[Eta] */
        v[38] = GaussIntegrationPoints(i36, 2);
        v[259] = (v[3] * v[38] * v[77]) / 2e0;
        v[256] = (v[3] * v[38] * v[76]) / 2e0;
        v[253] = (v[3] * v[38] * v[75]) / 2e0;
        /* 39 = \[Zeta] */
        v[39] = GaussIntegrationPoints(i36, 3);
        /* 40 = wgp */
        v[40] = GaussIntegrationPoints(i36, 4);
        /* 41 = Nh_1 */
        v[41] = (v[37] * v[42]) / 2e0;
        /* 43 = Nh_2 */
        v[43] = -(v[42] * v[44]);
        /* 45 = Nh_3 */
        v[45] = (v[37] * v[44]) / 2e0;
        /* 46 = X0_1 */
        v[46] = v[11] * v[43] + v[14] * v[45] + v[41] * v[8];
        /* 47 = X0_2 */
        v[47] = v[12] * v[43] + v[15] * v[45] + v[41] * v[9];
        /* 48 = X0_3 */
        v[48] = v[10] * v[41] + v[13] * v[43] + v[16] * v[45];
        /* 56 = t1_1 */
        v[56] = v[53] * v[57];
        /* 58 = t1_2 */
        v[58] = v[54] * v[57];
        /* 59 = t1_3 */
        v[59] = v[55] * v[57];
        /* 60 = t2_1 */
        v[60] = v[59] * v[6] - v[58] * v[7];
        /* 61 = t2_2 */
        v[61] = -(v[5] * v[59]) + v[56] * v[7];
        /* 80 = t3_3;\[Xi] */
        v[80] = v[61] * v[72] - v[60] * v[73] - v[58] * v[75] + v[56] * v[76];
        v[260] = (v[39] * v[4] * v[80]) / 2e0;
        v[266] = v[259] + v[260];
        /* 62 = t2_3 */
        v[62] = v[5] * v[58] - v[56] * v[6];
        /* 79 = t3_2;\[Xi] */
        v[79] = -(v[62] * v[72]) + v[60] * v[74] + v[59] * v[75] - v[56] * v[77];
        v[257] = (v[39] * v[4] * v[79]) / 2e0;
        v[264] = v[256] + v[257];
        /* 78 = t3_1;\[Xi] */
        v[78] = v[62] * v[73] - v[61] * v[74] - v[59] * v[76] + v[58] * v[77];
        v[254] = (v[39] * v[4] * v[78]) / 2e0;
        v[262] = v[253] + v[254];
        /* 63 = t3_1 */
        v[63] = -(v[59] * v[61]) + v[58] * v[62];
        v[246] = (v[3] * v[38] * v[60]) / 2e0 + (v[39] * v[4] * v[63]) / 2e0;
        /* 64 = t3_2 */
        v[64] = v[59] * v[60] - v[56] * v[62];
        v[244] = (v[3] * v[38] * v[61]) / 2e0 + (v[39] * v[4] * v[64]) / 2e0;
        /* 65 = t3_3 */
        v[65] = -(v[58] * v[60]) + v[56] * v[61];
        v[242] = (v[3] * v[38] * v[62]) / 2e0 + (v[39] * v[4] * v[65]) / 2e0;
        /* 81 = Je_1|1 */
        v[81] = v[253] + v[254] + v[53];
        /* 82 = Je_1|2 */
        v[82] = (v[3] * v[60]) / 2e0;
        /* 83 = Je_1|3 */
        v[83] = (v[4] * v[63]) / 2e0;
        /* 84 = Je_2|1 */
        v[84] = v[256] + v[257] + v[54];
        /* 85 = Je_2|2 */
        v[85] = (v[3] * v[61]) / 2e0;
        v[195] = -(v[82] * v[84]) + v[81] * v[85];
        /* 86 = Je_2|3 */
        v[86] = (v[4] * v[64]) / 2e0;
        v[194] = v[83] * v[84] - v[81] * v[86];
        v[193] = -(v[83] * v[85]) + v[82] * v[86];
        /* 87 = Je_3|1 */
        v[87] = v[259] + v[260] + v[55];
        /* 88 = Je_3|2 */
        v[88] = (v[3] * v[62]) / 2e0;
        v[201] = -(v[85] * v[87]) + v[84] * v[88];
        v[198] = v[82] * v[87] - v[81] * v[88];
        /* 89 = Je_3|3 */
        v[89] = (v[4] * v[65]) / 2e0;
        v[200] = v[86] * v[87] - v[84] * v[89];
        v[199] = -(v[86] * v[88]) + v[85] * v[89];
        v[197] = -(v[83] * v[87]) + v[81] * v[89];
        v[196] = v[83] * v[88] - v[82] * v[89];
        /* 90 = Jed */
        v[90] = -(v[83] * v[85] * v[87]) + v[82] * v[86] * v[87] + v[83] * v[84] * v[88] - v[81] * v[86] * v[88] - v[82] * v[84] * v[89]
            + v[81] * v[85] * v[89];
        v[240] = (v[200] * v[56]) / v[90] + (v[197] * v[58]) / v[90] + (v[194] * v[59]) / v[90];
        v[238] = (v[201] * v[56]) / v[90] + (v[198] * v[58]) / v[90] + (v[195] * v[59]) / v[90];
        v[229] = (v[199] * v[56]) / v[90] + (v[196] * v[58]) / v[90] + (v[193] * v[59]) / v[90];
        /* 91 = u\[Phi]_1 */
        v[91] = v[17] * v[41] + v[23] * v[43] + v[29] * v[45];
        /* 92 = u\[Phi]_2 */
        v[92] = v[18] * v[41] + v[24] * v[43] + v[30] * v[45];
        /* 93 = u\[Phi]_3 */
        v[93] = v[19] * v[41] + v[25] * v[43] + v[31] * v[45];
        /* 94 = u\[Phi]_4 */
        v[94] = v[20] * v[41] + v[26] * v[43] + v[32] * v[45];
        v[139] = 2e0 * v[138] * v[94];
        v[110] = (v[94] * v[94]);
        /* 95 = u\[Phi]_5 */
        v[95] = v[21] * v[41] + v[27] * v[43] + v[33] * v[45];
        v[481] = v[140] * v[94] + v[138] * v[95];
        v[141] = 2e0 * v[140] * v[95];
        v[142] = v[139] + v[141];
        v[100] = (v[95] * v[95]);
        v[132] = v[100] + v[110];
        /* 96 = u\[Phi]_6 */
        v[96] = v[22] * v[41] + v[28] * v[43] + v[34] * v[45];
        v[487] = v[143] * v[95] + v[140] * v[96];
        v[484] = v[143] * v[94] + v[138] * v[96];
        v[144] = 2e0 * v[143] * v[96];
        /* 147 = ff_;\[Xi] */
        v[147] = v[139] + v[141] + v[144];
        v[146] = v[141] + v[144];
        v[145] = v[139] + v[144];
        v[101] = (v[96] * v[96]);
        v[127] = v[101] + v[110];
        v[120] = v[100] + v[101];
        /* 97 = ff */
        v[97] = v[100] + v[101] + v[110];
        b98 = v[97] < 0.1e-7;
        v[98] = b98;
        b99 = b98;
        v[99] = b99;
        if (b99) {
            v[154] = v[147] * v[97];
            v[153] = -30e0 + v[97];
            v[155] = v[147] * v[153] + v[154];
            v[148] = -20e0 + v[97];
            v[149] = v[147] * v[148] + v[154];
            v[105] = 120e0 + v[148] * v[97];
            v[152] = -6e0 * v[105] * v[143] - 6e0 * v[149] * v[96];
            v[151] = 6e0 * v[105] * v[140] + 6e0 * v[149] * v[95];
            v[150] = -6e0 * v[105] * v[138] - 6e0 * v[149] * v[94];
            v[117] = -6e0 * v[105] * v[94];
            v[114] = 6e0 * v[105] * v[95];
            v[108] = -6e0 * v[105] * v[96];
            v[103] = 360e0 + v[153] * v[97];
            v[158] = v[103] * v[140] * v[94] + v[103] * v[138] * v[95] + v[155] * v[94] * v[95];
            v[157] = v[103] * v[143] * v[94] + v[103] * v[138] * v[96] + v[155] * v[94] * v[96];
            v[156] = v[103] * v[143] * v[95] + v[103] * v[140] * v[96] + v[155] * v[95] * v[96];
            v[116] = v[103] * v[95] * v[96];
            v[113] = v[103] * v[94] * v[96];
            v[107] = v[103] * v[94] * v[95];
            /* 159 = v_1|1;\[Xi] */
            v[159] = (-1e0 / 720e0) * (v[103] * v[146]) - (v[120] * v[155]) / 720e0;
            /* 102 = v_1|1 */
            v[102] = 1e0 - (v[103] * v[120]) / 720e0;
            /* 160 = v_1|2;\[Xi] */
            v[160] = (v[152] + v[158]) / 720e0;
            /* 104 = v_1|2 */
            v[104] = (v[107] + v[108]) / 720e0;
            /* 161 = v_1|3;\[Xi] */
            v[161] = (v[151] + v[157]) / 720e0;
            /* 106 = v_1|3 */
            v[106] = (v[113] + v[114]) / 720e0;
            /* 162 = v_2|1;\[Xi] */
            v[162] = (-v[152] + v[158]) / 720e0;
            /* 109 = v_2|1 */
            v[109] = (v[107] - v[108]) / 720e0;
            /* 163 = v_2|2;\[Xi] */
            v[163] = (-1e0 / 720e0) * (v[103] * v[145]) - (v[127] * v[155]) / 720e0;
            /* 111 = v_2|2 */
            v[111] = 1e0 - (v[103] * v[127]) / 720e0;
            /* 164 = v_2|3;\[Xi] */
            v[164] = (v[150] + v[156]) / 720e0;
            /* 112 = v_2|3 */
            v[112] = (v[116] + v[117]) / 720e0;
            /* 165 = v_3|1;\[Xi] */
            v[165] = (-v[151] + v[157]) / 720e0;
            /* 115 = v_3|1 */
            v[115] = (v[113] - v[114]) / 720e0;
            /* 166 = v_3|2;\[Xi] */
            v[166] = (-v[150] + v[156]) / 720e0;
            /* 118 = v_3|2 */
            v[118] = (v[116] - v[117]) / 720e0;
            /* 167 = v_3|3;\[Xi] */
            v[167] = (-1e0 / 720e0) * (v[103] * v[142]) - (v[132] * v[155]) / 720e0;
            /* 119 = v_3|3 */
            v[119] = 1e0 - (v[103] * v[132]) / 720e0;
        }
        else {
            v[178] = 1e0 / (v[97] * v[97]);
            v[168] = 1e0 / sqrt(v[97]);
            v[169] = (v[147] * v[168]) / 2e0;
            v[121] = sqrt(v[97]);
            v[338] = 1e0 / (v[121] * v[121]);
            v[171] = cos(v[121]);
            v[172] = v[169] * v[171];
            v[170] = -(v[169] * v[338]);
            v[752] = -(v[143] * v[168]) - v[170] * v[96];
            v[751] = v[140] * v[168] + v[170] * v[95];
            v[750] = -(v[138] * v[168]) - v[170] * v[94];
            v[124] = v[168];
            v[123] = sin(v[121]);
            v[763] = -(v[123] * v[143]) - v[172] * v[96];
            v[762] = v[123] * v[140] + v[172] * v[95];
            v[761] = -(v[123] * v[138]) - v[172] * v[94];
            v[176] = -(v[123] * v[169]);
            v[760] = v[176] * v[178] * v[94] * v[95];
            v[757] = v[176] * v[178] * v[94] * v[96];
            v[754] = v[176] * v[178] * v[95] * v[96];
            v[175] = -(v[123] * v[124] * v[143]) - v[123] * v[170] * v[96] - v[124] * v[172] * v[96];
            v[174] = v[123] * v[124] * v[140] + v[123] * v[170] * v[95] + v[124] * v[172] * v[95];
            v[173] = -(v[123] * v[124] * v[138]) - v[123] * v[170] * v[94] - v[124] * v[172] * v[94];
            v[130] = -(v[123] * v[124] * v[94]);
            v[128] = v[123] * v[124] * v[95];
            v[125] = -(v[123] * v[124] * v[96]);
            v[122] = -1e0 + v[171];
            v[759] = v[122] * v[140] * v[178] * v[94];
            v[758] = v[122] * v[138] * v[178] * v[95];
            v[756] = v[122] * v[143] * v[178] * v[94];
            v[755] = v[122] * v[138] * v[178] * v[96];
            v[753] = v[122] * v[143] * v[178] * v[95];
            v[574] = v[122] * v[138] * v[178] + v[176] * v[178] * v[94];
            v[572] = v[122] * v[140] * v[178] + v[176] * v[178] * v[95];
            v[552] = v[122] * v[143] * v[178] + v[176] * v[178] * v[96];
            v[180] = v[122] * v[147] * v[178] * v[94] * v[95] - (v[122] * v[140] * v[94]) / v[97] - (v[122] * v[138] * v[95]) / v[97] -
                (v[176] * v[94] * v[95]) / v[97];
            v[179] = v[122] * v[147] * v[178] * v[94] * v[96] - (v[122] * v[143] * v[94]) / v[97] - (v[122] * v[138] * v[96]) / v[97] -
                (v[176] * v[94] * v[96]) / v[97];
            v[177] = v[122] * v[147] * v[178] * v[95] * v[96] - (v[122] * v[143] * v[95]) / v[97] - (v[122] * v[140] * v[96]) / v[97] -
                (v[176] * v[95] * v[96]) / v[97];
            v[131] = -((v[122] * v[95] * v[96]) / v[97]);
            v[129] = -((v[122] * v[94] * v[96]) / v[97]);
            v[126] = -((v[122] * v[94] * v[95]) / v[97]);
            /* 159 = v_1|1;\[Xi] */
            v[159] = -(v[120] * v[122] * v[147] * v[178]) + (v[122] * v[146]) / v[97] + (v[120] * v[176]) / v[97];
            /* 102 = v_1|1 */
            v[102] = 1e0 + (v[120] * v[122]) / v[97];
            /* 160 = v_1|2;\[Xi] */
            v[160] = v[175] + v[180];
            /* 104 = v_1|2 */
            v[104] = v[125] + v[126];
            /* 161 = v_1|3;\[Xi] */
            v[161] = v[174] + v[179];
            /* 106 = v_1|3 */
            v[106] = v[128] + v[129];
            /* 162 = v_2|1;\[Xi] */
            v[162] = -v[175] + v[180];
            /* 109 = v_2|1 */
            v[109] = -v[125] + v[126];
            /* 163 = v_2|2;\[Xi] */
            v[163] = -(v[122] * v[127] * v[147] * v[178]) + (v[122] * v[145]) / v[97] + (v[127] * v[176]) / v[97];
            /* 111 = v_2|2 */
            v[111] = 1e0 + (v[122] * v[127]) / v[97];
            /* 164 = v_2|3;\[Xi] */
            v[164] = v[173] + v[177];
            /* 112 = v_2|3 */
            v[112] = v[130] + v[131];
            /* 165 = v_3|1;\[Xi] */
            v[165] = -v[174] + v[179];
            /* 115 = v_3|1 */
            v[115] = -v[128] + v[129];
            /* 166 = v_3|2;\[Xi] */
            v[166] = -v[173] + v[177];
            /* 118 = v_3|2 */
            v[118] = -v[130] + v[131];
            /* 167 = v_3|3;\[Xi] */
            v[167] = -(v[122] * v[132] * v[147] * v[178]) + (v[122] * v[142]) / v[97] + (v[132] * v[176]) / v[97];
            /* 119 = v_3|3 */
            v[119] = 1e0 + (v[122] * v[132]) / v[97];
        };
        /* 207 = Ft_3|3 */
        v[207] = v[115] * v[63] + v[118] * v[64] + v[119] * v[65];
        /* 206 = Ft_3|2 */
        v[206] = v[115] * v[60] + v[118] * v[61] + v[119] * v[62];
        /* 204 = Ft_2|3 */
        v[204] = v[109] * v[63] + v[111] * v[64] + v[112] * v[65];
        /* 203 = Ft_2|2 */
        v[203] = v[109] * v[60] + v[111] * v[61] + v[112] * v[62];
        /* 192 = Ft_1|3 */
        v[192] = v[102] * v[63] + v[104] * v[64] + v[106] * v[65];
        /* 191 = Ft_1|2 */
        v[191] = v[102] * v[60] + v[104] * v[61] + v[106] * v[62];
        v[232] = v[191] * v[192] + v[203] * v[204] + v[206] * v[207];
        /* 181 = Fref_1|1 */
        v[181] = v[135] + v[53] + (v[3] * v[38] * (v[159] * v[60] + v[160] * v[61] + v[161] * v[62] + v[102] * v[75] + v[104] * v[76]
            + v[106] * v[77])) / 2e0 + (v[39] * v[4] * (v[159] * v[63] + v[160] * v[64] + v[161] * v[65] + v[102] * v[78] + v[104] * v[79]
                + v[106] * v[80])) / 2e0;
        /* 182 = Fref_1|2 */
        v[182] = (v[191] * v[3]) / 2e0;
        /* 183 = Fref_1|3 */
        v[183] = (v[192] * v[4]) / 2e0;
        /* 184 = Fref_2|1 */
        v[184] = v[136] + v[54] + (v[3] * v[38] * (v[162] * v[60] + v[163] * v[61] + v[164] * v[62] + v[109] * v[75] + v[111] * v[76]
            + v[112] * v[77])) / 2e0 + (v[39] * v[4] * (v[162] * v[63] + v[163] * v[64] + v[164] * v[65] + v[109] * v[78] + v[111] * v[79]
                + v[112] * v[80])) / 2e0;
        /* 185 = Fref_2|2 */
        v[185] = (v[203] * v[3]) / 2e0;
        /* 186 = Fref_2|3 */
        v[186] = (v[204] * v[4]) / 2e0;
        /* 187 = Fref_3|1 */
        v[187] = v[137] + v[55] + (v[3] * v[38] * (v[165] * v[60] + v[166] * v[61] + v[167] * v[62] + v[115] * v[75] + v[118] * v[76]
            + v[119] * v[77])) / 2e0 + (v[39] * v[4] * (v[165] * v[63] + v[166] * v[64] + v[167] * v[65] + v[115] * v[78] + v[118] * v[79]
                + v[119] * v[80])) / 2e0;
        /* 188 = Fref_3|2 */
        v[188] = (v[206] * v[3]) / 2e0;
        /* 189 = Fref_3|3 */
        v[189] = (v[207] * v[4]) / 2e0;
        /* 190 = Ft_1|1 */
        v[190] = v[59] * ((v[181] * v[193]) / v[90] + (v[182] * v[194]) / v[90] + (v[183] * v[195]) / v[90]) + v[58] * (
            (v[181] * v[196]) / v[90] + (v[182] * v[197]) / v[90] + (v[183] * v[198]) / v[90]) + v[56] * ((v[181] * v[199]) / v[90] +
                (v[182] * v[200]) / v[90] + (v[183] * v[201]) / v[90]);
        /* 202 = Ft_2|1 */
        v[202] = v[59] * ((v[184] * v[193]) / v[90] + (v[185] * v[194]) / v[90] + (v[186] * v[195]) / v[90]) + v[58] * (
            (v[184] * v[196]) / v[90] + (v[185] * v[197]) / v[90] + (v[186] * v[198]) / v[90]) + v[56] * ((v[184] * v[199]) / v[90] +
                (v[185] * v[200]) / v[90] + (v[186] * v[201]) / v[90]);
        /* 205 = Ft_3|1 */
        v[205] = v[59] * ((v[187] * v[193]) / v[90] + (v[188] * v[194]) / v[90] + (v[189] * v[195]) / v[90]) + v[58] * (
            (v[187] * v[196]) / v[90] + (v[188] * v[197]) / v[90] + (v[189] * v[198]) / v[90]) + v[56] * ((v[187] * v[199]) / v[90] +
                (v[188] * v[200]) / v[90] + (v[189] * v[201]) / v[90]);
        v[222] = -1e0 + (v[190] * v[190]) + (v[202] * v[202]) + (v[205] * v[205]);
        v[221] = v[190] * v[192] + v[202] * v[204] + v[205] * v[207];
        v[220] = v[190] * v[191] + v[202] * v[203] + v[205] * v[206];
        /* 208 = Et_1|1 */
        v[208] = v[222] / 2e0;
        /* 209 = [Et_2|1][Et_1|2] */
        v[209] = v[220] / 2e0;
        /* 210 = [Et_3|1][Et_1|3] */
        v[210] = v[221] / 2e0;
        /* 211 = Et_2|2 */
        v[211] = (-1e0 + (v[191] * v[191]) + (v[203] * v[203]) + (v[206] * v[206])) / 2e0;
        /* 212 = [Et_3|2][Et_2|3] */
        v[212] = v[232] / 2e0;
        /* 213 = Et_3|3 */
        v[213] = (-1e0 + (v[192] * v[192]) + (v[204] * v[204]) + (v[207] * v[207])) / 2e0;
        /* 217 = W */
        v[217] = (v[216] * (4e0 * ((v[209] * v[209]) + (v[210] * v[210]) + (v[212] * v[212])) + (v[208] * v[208]
            ) * v[223] * v[224])) / 2e0;
        /* 225 = \[OverBracket]_Ft_3|1(W|W) */
        v[225] = (v[216] * (4e0 * ((v[206] * v[220]) / 2e0 + (v[207] * v[221]) / 2e0) + v[205] * v[222] * v[223] * v[224])) / 2e0;
        /* 226 = \[OverBracket]_Ft_2|1(W|W) */
        v[226] = (v[216] * (4e0 * ((v[203] * v[220]) / 2e0 + (v[204] * v[221]) / 2e0) + v[202] * v[222] * v[223] * v[224])) / 2e0;
        /* 227 = \[OverBracket]_Ft_1|1(W|W) */
        v[227] = (v[216] * (4e0 * ((v[191] * v[220]) / 2e0 + (v[192] * v[221]) / 2e0) + v[190] * v[222] * v[223] * v[224])) / 2e0;
        /* 228 = \[OverBracket]_Fref_3|1(Ft|W)_3|1 */
        v[228] = v[225] * v[229];
        /* 230 = \[OverBracket]_Fref_2|1(Ft|W)_2|1 */
        v[230] = v[226] * v[229];
        /* 231 = \[OverBracket]_Fref_1|1(Ft|W)_1|1 */
        v[231] = v[227] * v[229];
        /* 233 = \[OverBracket]_Ft_3|3(Fref|W)_3|3 */
        v[233] = 2e0 * v[216] * ((v[205] * v[221]) / 2e0 + (v[206] * v[232]) / 2e0) + (v[225] * v[238] * v[4]) / 2e0;
        /* 234 = \[OverBracket]_Ft_3|2(Fref|W)_3|2 */
        v[234] = 2e0 * v[216] * ((v[205] * v[220]) / 2e0 + (v[207] * v[232]) / 2e0) + (v[225] * v[240] * v[3]) / 2e0;
        /* 235 = \[OverBracket]_R_3|3;\[Xi](Fref|W)_3|1 */
        v[235] = v[228] * v[242];
        /* 236 = \[OverBracket]_R_3|2;\[Xi](Fref|W)_3|1 */
        v[236] = v[228] * v[244];
        /* 237 = \[OverBracket]_R_3|1;\[Xi](Fref|W)_3|1 */
        v[237] = v[228] * v[246];
        /* 239 = \[OverBracket]_Ft_2|3(Fref|W)_2|3 */
        v[239] = 2e0 * v[216] * ((v[202] * v[221]) / 2e0 + (v[203] * v[232]) / 2e0) + (v[226] * v[238] * v[4]) / 2e0;
        /* 241 = \[OverBracket]_Ft_2|2(Fref|W)_2|2 */
        v[241] = 2e0 * v[216] * ((v[202] * v[220]) / 2e0 + (v[204] * v[232]) / 2e0) + (v[226] * v[240] * v[3]) / 2e0;
        /* 243 = \[OverBracket]_R_2|3;\[Xi](Fref|W)_2|1 */
        v[243] = v[230] * v[242];
        /* 245 = \[OverBracket]_R_2|2;\[Xi](Fref|W)_2|1 */
        v[245] = v[230] * v[244];
        /* 247 = \[OverBracket]_R_2|1;\[Xi](Fref|W)_2|1 */
        v[247] = v[230] * v[246];
        /* 248 = \[OverBracket]_Ft_1|3(Fref|W)_1|3 */
        v[248] = 2e0 * v[216] * ((v[190] * v[221]) / 2e0 + (v[191] * v[232]) / 2e0) + (v[227] * v[238] * v[4]) / 2e0;
        /* 249 = \[OverBracket]_Ft_1|2(Fref|W)_1|2 */
        v[249] = 2e0 * v[216] * ((v[190] * v[220]) / 2e0 + (v[192] * v[232]) / 2e0) + (v[227] * v[240] * v[3]) / 2e0;
        /* 250 = \[OverBracket]_R_1|3;\[Xi](Fref|W)_1|1 */
        v[250] = v[231] * v[242];
        /* 251 = \[OverBracket]_R_1|2;\[Xi](Fref|W)_1|1 */
        v[251] = v[231] * v[244];
        /* 252 = \[OverBracket]_R_1|1;\[Xi](Fref|W)_1|1 */
        v[252] = v[231] * v[246];
        /* 255 = \[OverBracket]_R_1|1(Ft|W)_1|3 */
        v[255] = v[231] * v[262] + v[249] * v[60] + v[248] * v[63];
        /* 258 = \[OverBracket]_R_1|2(Ft|W)_1|3 */
        v[258] = v[231] * v[264] + v[249] * v[61] + v[248] * v[64];
        /* 261 = \[OverBracket]_R_1|3(Ft|W)_1|3 */
        v[261] = v[231] * v[266] + v[249] * v[62] + v[248] * v[65];
        /* 263 = \[OverBracket]_R_2|1(Ft|W)_2|3 */
        v[263] = v[230] * v[262] + v[241] * v[60] + v[239] * v[63];
        /* 265 = \[OverBracket]_R_2|2(Ft|W)_2|3 */
        v[265] = v[230] * v[264] + v[241] * v[61] + v[239] * v[64];
        /* 267 = \[OverBracket]_R_2|3(Ft|W)_2|3 */
        v[267] = v[230] * v[266] + v[241] * v[62] + v[239] * v[65];
        /* 268 = \[OverBracket]_R_3|1(Ft|W)_3|3 */
        v[268] = v[228] * v[262] + v[234] * v[60] + v[233] * v[63];
        /* 269 = \[OverBracket]_R_3|2(Ft|W)_3|3 */
        v[269] = v[228] * v[264] + v[234] * v[61] + v[233] * v[64];
        /* 270 = \[OverBracket]_R_3|3(Ft|W)_3|3 */
        v[270] = v[228] * v[266] + v[234] * v[62] + v[233] * v[65];
        b271 = b98;
        v[271] = b271;
        if (b271) {
            v[312] = v[103] * v[138] + v[155] * v[94];
            v[311] = v[103] * v[140] + v[155] * v[95];
            v[306] = v[103] * v[143] + v[155] * v[96];
            v[298] = v[251] / 720e0;
            v[297] = v[247] / 720e0;
            v[294] = v[263] / 720e0;
            v[293] = v[258] / 720e0;
            v[290] = v[250] / 720e0;
            v[289] = v[237] / 720e0;
            v[286] = v[268] / 720e0;
            v[285] = v[261] / 720e0;
            v[280] = v[243] / 720e0;
            v[279] = v[236] / 720e0;
            v[276] = v[269] / 720e0;
            v[275] = v[267] / 720e0;
            /* 272 = \[OverBracket]_\[Yen]_132(v|W)_3|3;\[Xi] */
            v[272] = (-1e0 / 720e0) * (v[155] * v[235]) - (v[103] * v[270]) / 720e0;
            /* 273 = \[OverBracket]_\[Yen]_142(v|W)_3|3;\[Xi] */
            v[273] = (-1e0 / 720e0) * (v[103] * v[235]);
            /* 274 = \[OverBracket]_\[Yen]_116(v|W)_2|3 */
            v[274] = v[275] + v[276];
            v[517] = v[274] * v[96];
            v[514] = v[274] * v[95];
            /* 277 = \[OverBracket]_\[Yen]_117(v|W)_2|3 */
            v[277] = v[275] - v[276];
            /* 278 = \[OverBracket]_\[Yen]_156(v|W)_2|3;\[Xi] */
            v[278] = v[279] + v[280];
            /* 281 = \[OverBracket]_\[Yen]_150(v|W)_2|3;\[Xi] */
            v[281] = -v[279] + v[280];
            /* 282 = \[OverBracket]_\[Yen]_127(v|W)_2|2;\[Xi] */
            v[282] = (-1e0 / 720e0) * (v[155] * v[245]) - (v[103] * v[265]) / 720e0;
            /* 283 = \[OverBracket]_\[Yen]_145(v|W)_2|2;\[Xi] */
            v[283] = (-1e0 / 720e0) * (v[103] * v[245]);
            /* 284 = \[OverBracket]_\[Yen]_113(v|W)_1|3 */
            v[284] = v[285] + v[286];
            v[520] = v[284] * v[96];
            v[515] = v[284] * v[94];
            /* 287 = \[OverBracket]_\[Yen]_114(v|W)_1|3 */
            v[287] = v[285] - v[286];
            /* 288 = \[OverBracket]_\[Yen]_157(v|W)_1|3;\[Xi] */
            v[288] = v[289] + v[290];
            v[516] = v[288] * v[94] + v[278] * v[95];
            /* 291 = \[OverBracket]_\[Yen]_151(v|W)_1|3;\[Xi] */
            v[291] = -v[289] + v[290];
            /* 292 = \[OverBracket]_\[Yen]_107(v|W)_1|2 */
            v[292] = v[293] + v[294];
            v[521] = v[292] * v[95];
            v[518] = v[292] * v[94];
            /* 295 = \[OverBracket]_\[Yen]_108(v|W)_1|2 */
            v[295] = v[293] - v[294];
            /* 296 = \[OverBracket]_\[Yen]_158(v|W)_1|2;\[Xi] */
            v[296] = v[297] + v[298];
            v[522] = v[296] * v[95] + v[288] * v[96];
            v[519] = v[296] * v[94] + v[278] * v[96];
            /* 299 = \[OverBracket]_\[Yen]_152(v|W)_1|2;\[Xi] */
            v[299] = -v[297] + v[298];
            /* 300 = \[OverBracket]_\[Yen]_120(v|W)_1|1;\[Xi] */
            v[300] = (-1e0 / 720e0) * (v[155] * v[252]) - (v[103] * v[255]) / 720e0;
            /* 301 = \[OverBracket]_\[Yen]_146(v|W)_1|1;\[Xi] */
            v[301] = (-1e0 / 720e0) * (v[103] * v[252]);
            /* 302 = \[OverBracket]_\[Yen]_103(\[Yen]|W)_158 */
            v[302] = (-1e0 / 720e0) * (v[142] * v[235]) - (v[145] * v[245]) / 720e0 - (v[146] * v[252]) / 720e0 - (v[120] * v[255])
                / 720e0 - (v[127] * v[265]) / 720e0 - (v[132] * v[270]) / 720e0 + v[296] * v[481] + v[288] * v[484] + v[278] * v[487]
                + v[292] * v[94] * v[95] + v[284] * v[94] * v[96] + v[274] * v[95] * v[96];
            /* 303 = \[OverBracket]_\[Yen]_155(\[Yen]|W)_158 */
            v[303] = (-1e0 / 720e0) * (v[132] * v[235]) - (v[127] * v[245]) / 720e0 - (v[120] * v[252]) / 720e0
                + v[296] * v[94] * v[95] + v[288] * v[94] * v[96] + v[278] * v[95] * v[96];
            /* 304 = \[OverBracket]_u\[Phi]_4(\[Yen]|W)_150 */
            v[304] = -6e0 * v[105] * v[277] - 6e0 * v[149] * v[281] + v[288] * v[306] + v[296] * v[311] + v[103] * v[292] * v[95]
                + v[103] * v[284] * v[96];
            /* 305 = \[OverBracket]_u\[Phi]_4;\[Xi](\[Yen]|W)_150 */
            v[305] = -6e0 * v[105] * v[281] + v[103] * v[296] * v[95] + v[103] * v[288] * v[96];
            /* 307 = \[OverBracket]_u\[Phi]_5(\[Yen]|W)_151 */
            v[307] = 6e0 * v[105] * v[287] + 6e0 * v[149] * v[291] + v[278] * v[306] + v[296] * v[312] + v[103] * v[292] * v[94]
                + v[103] * v[274] * v[96];
            /* 308 = \[OverBracket]_u\[Phi]_5;\[Xi](\[Yen]|W)_151 */
            v[308] = 6e0 * v[105] * v[291] + v[103] * v[296] * v[94] + v[103] * v[278] * v[96];
            /* 309 = \[OverBracket]_\[Yen]_105(\[Yen]|W)_152 */
            v[309] = -6e0 * v[138] * v[281] + 6e0 * v[140] * v[291] - 6e0 * v[143] * v[299] - 6e0 * v[277] * v[94] + 6e0 * v[287] * v[95]
                - 6e0 * v[295] * v[96];
            /* 310 = \[OverBracket]_\[Yen]_149(\[Yen]|W)_152 */
            v[310] = -6e0 * v[281] * v[94] + 6e0 * v[291] * v[95] - 6e0 * v[299] * v[96];
            /* 313 = \[OverBracket]_u\[Phi]_6(\[Yen]|W)_152 */
            v[313] = -6e0 * v[105] * v[295] - 6e0 * v[149] * v[299] + v[278] * v[311] + v[288] * v[312] + v[103] * v[284] * v[94]
                + v[103] * v[274] * v[95];
            /* 314 = \[OverBracket]_u\[Phi]_6;\[Xi](\[Yen]|W)_152 */
            v[314] = -6e0 * v[105] * v[299] + v[103] * v[288] * v[94] + v[103] * v[278] * v[95];
            /* 315 = \[OverBracket]_\[Yen]_154(\[Yen]|W)_155 */
            v[315] = v[303] + v[310];
            /* 316 = \[OverBracket]_ff_(\[Yen]|W)_154 */
            v[316] = v[153] * v[302] + v[147] * v[303] + v[148] * v[309] + v[147] * v[310] + v[147] * v[315] + v[302] * v[97]
                + v[309] * v[97];
            /* 317 = \[OverBracket]_ff_;\[Xi](\[Yen]|W)_154 */
            v[317] = v[153] * v[303] + v[148] * v[310] + v[315] * v[97];
        }
        else {
            v[337] = v[122] * v[147] * v[178] * v[94] - (v[122] * v[138]) / v[97] - (v[176] * v[94]) / v[97];
            v[336] = v[122] * v[147] * v[178] * v[95] - (v[122] * v[140]) / v[97] - (v[176] * v[95]) / v[97];
            v[333] = -(v[123] * v[170]) - v[168] * v[172];
            v[332] = v[122] * v[147] * v[178] * v[96] - (v[122] * v[143]) / v[97] - (v[176] * v[96]) / v[97];
            v[322] = -(v[122] * v[147] * v[178]) + v[176] / v[97];
            /* 272 = \[OverBracket]_\[Yen]_132(v|W)_3|3;\[Xi] */
            v[272] = v[235] * v[322] + (v[122] * v[270]) / v[97];
            /* 273 = \[OverBracket]_\[Yen]_142(v|W)_3|3;\[Xi] */
            v[273] = (v[122] * v[235]) / v[97];
            /* 318 = \[OverBracket]_\[Yen]_131(v|W)_2|3 */
            v[318] = v[267] + v[269];
            v[571] = v[122] * v[178] * v[318] * v[95];
            /* 319 = \[OverBracket]_\[Yen]_130(v|W)_2|3 */
            v[319] = v[267] - v[269];
            /* 320 = \[OverBracket]_\[Yen]_177(v|W)_2|3;\[Xi] */
            v[320] = v[236] + v[243];
            v[618] = v[122] * v[178] * v[320] * v[96];
            v[615] = v[122] * v[178] * v[320] * v[95];
            v[546] = v[320] * (v[143] * v[178] * v[95] + v[140] * v[178] * v[96]);
            /* 321 = \[OverBracket]_\[Yen]_173(v|W)_2|3;\[Xi] */
            v[321] = -v[236] + v[243];
            /* 282 = \[OverBracket]_\[Yen]_127(v|W)_2|2;\[Xi] */
            v[282] = v[245] * v[322] + (v[122] * v[265]) / v[97];
            /* 283 = \[OverBracket]_\[Yen]_145(v|W)_2|2;\[Xi] */
            v[283] = (v[122] * v[245]) / v[97];
            /* 323 = \[OverBracket]_\[Yen]_129(v|W)_1|3 */
            v[323] = v[261] + v[268];
            v[573] = v[122] * v[178] * v[323] * v[94];
            v[547] = v[178] * v[323] * v[94] * v[96];
            /* 324 = \[OverBracket]_\[Yen]_128(v|W)_1|3 */
            v[324] = v[261] - v[268];
            /* 325 = \[OverBracket]_\[Yen]_179(v|W)_1|3;\[Xi] */
            v[325] = v[237] + v[250];
            v[623] = v[122] * v[178] * v[325] * v[96];
            v[616] = v[122] * v[178] * v[325] * v[94];
            v[570] = v[615] + v[616];
            v[548] = v[325] * (v[143] * v[178] * v[94] + v[138] * v[178] * v[96]);
            /* 326 = \[OverBracket]_\[Yen]_174(v|W)_1|3;\[Xi] */
            v[326] = -v[237] + v[250];
            /* 327 = \[OverBracket]_\[Yen]_126(v|W)_1|2 */
            v[327] = v[258] + v[263];
            v[606] = v[122] * v[178] * v[327] * v[95];
            v[590] = v[122] * v[178] * v[327] * v[94];
            v[549] = v[178] * v[327] * v[94] * v[95];
            /* 328 = \[OverBracket]_\[Yen]_125(v|W)_1|2 */
            v[328] = v[258] - v[263];
            /* 329 = \[OverBracket]_\[Yen]_180(v|W)_1|2;\[Xi] */
            v[329] = v[247] + v[251];
            v[621] = -(v[122] * v[132] * v[235]) - v[122] * v[127] * v[245] - v[120] * v[122] * v[252]
                + v[122] * v[329] * v[94] * v[95] + v[122] * v[325] * v[94] * v[96] + v[122] * v[320] * v[95] * v[96];
            v[620] = -(v[132] * v[147] * v[235]) - v[127] * v[147] * v[245] - v[120] * v[147] * v[252]
                + v[147] * v[329] * v[94] * v[95] + v[147] * v[325] * v[94] * v[96] + v[147] * v[320] * v[95] * v[96];
            v[619] = -(v[132] * v[178] * v[235]) - v[127] * v[178] * v[245] - v[120] * v[178] * v[252]
                + v[178] * v[329] * v[94] * v[95] + v[178] * v[325] * v[94] * v[96] + v[178] * v[320] * v[95] * v[96];
            v[594] = v[122] * v[178] * v[329] * v[95];
            v[569] = v[122] * v[178] * v[329] * v[94];
            v[551] = -(v[122] * v[132] * v[147] * v[235]) - v[122] * v[127] * v[147] * v[245] - v[120] * v[122] * v[147] * v[252]
                + v[122] * v[147] * v[329] * v[94] * v[95] + v[122] * v[147] * v[325] * v[94] * v[96]
                + v[122] * v[147] * v[320] * v[95] * v[96];
            v[550] = v[329] * (v[140] * v[178] * v[94] + v[138] * v[178] * v[95]);
            /* 330 = \[OverBracket]_\[Yen]_175(v|W)_1|2;\[Xi] */
            v[330] = -v[247] + v[251];
            v[568] = -(v[321] * v[94]) + v[326] * v[95] - v[330] * v[96];
            v[542] = -(v[138] * v[321]) + v[140] * v[326] - v[143] * v[330] - v[319] * v[94] + v[324] * v[95] - v[328] * v[96];
            /* 300 = \[OverBracket]_\[Yen]_120(v|W)_1|1;\[Xi] */
            v[300] = v[252] * v[322] + (v[122] * v[255]) / v[97];
            /* 301 = \[OverBracket]_\[Yen]_146(v|W)_1|1;\[Xi] */
            v[301] = (v[122] * v[252]) / v[97];
            /* 331 = \[OverBracket]_\[Yen]_176(\[Yen]|W)_180 */
            v[331] = (v[132] * v[235]) / v[97] + (v[127] * v[245]) / v[97] + (v[120] * v[252]) / v[97] - (v[329] * v[94] * v[95])
                / v[97] - (v[325] * v[94] * v[96]) / v[97] - (v[320] * v[95] * v[96]) / v[97];
            v[534] = -(v[169] * v[331]) - v[168] * v[319] * v[94] + v[321] * (-(v[138] * v[168]) - v[170] * v[94])
                + v[168] * v[324] * v[95] + v[326] * (v[140] * v[168] + v[170] * v[95]) - v[168] * v[328] * v[96] + v[330] * (-
                    (v[143] * v[168]) - v[170] * v[96]);
            v[543] = v[171] * v[534];
            /* 304 = \[OverBracket]_u\[Phi]_4(\[Yen]|W)_173 */
            v[304] = -(v[123] * v[168] * v[319]) + v[325] * v[332] + v[321] * v[333] + v[329] * v[336] - (v[122] * v[327] * v[95])
                / v[97] - (v[122] * v[323] * v[96]) / v[97];
            /* 305 = \[OverBracket]_u\[Phi]_4;\[Xi](\[Yen]|W)_173 */
            v[305] = -(v[123] * v[168] * v[321]) - (v[122] * v[329] * v[95]) / v[97] - (v[122] * v[325] * v[96]) / v[97];
            /* 307 = \[OverBracket]_u\[Phi]_5(\[Yen]|W)_174 */
            v[307] = v[123] * v[168] * v[324] + v[320] * v[332] - v[326] * v[333] + v[329] * v[337] - (v[122] * v[327] * v[94]) / v[97]
                - (v[122] * v[318] * v[96]) / v[97];
            /* 308 = \[OverBracket]_u\[Phi]_5;\[Xi](\[Yen]|W)_174 */
            v[308] = v[123] * v[168] * v[326] - (v[122] * v[329] * v[94]) / v[97] - (v[122] * v[320] * v[96]) / v[97];
            /* 334 = \[OverBracket]_\[Yen]_170(\[Yen]|W)_175 */
            v[334] = -(v[123] * v[321] * v[94]) + v[123] * v[326] * v[95] - v[123] * v[330] * v[96];
            /* 335 = \[OverBracket]_\[Yen]_172(\[Yen]|W)_175 */
            v[335] = -(v[168] * v[321] * v[94]) + v[168] * v[326] * v[95] - v[168] * v[330] * v[96];
            v[541] = v[169] * v[335];
            /* 313 = \[OverBracket]_u\[Phi]_6(\[Yen]|W)_175 */
            v[313] = -(v[123] * v[168] * v[328]) + v[330] * v[333] + v[320] * v[336] + v[325] * v[337] - (v[122] * v[323] * v[94])
                / v[97] - (v[122] * v[318] * v[95]) / v[97];
            /* 314 = \[OverBracket]_u\[Phi]_6;\[Xi](\[Yen]|W)_175 */
            v[314] = -(v[123] * v[168] * v[330]) - (v[122] * v[325] * v[94]) / v[97] - (v[122] * v[320] * v[95]) / v[97];
            /* 339 = \[OverBracket]_\[Yen]_169(\[Yen]|W)_172 */
            v[339] = -(v[123] * v[331]) + v[171] * v[335] - v[334] * v[338];
            v[545] = (v[147] * v[339]) / 2e0 - v[123] * v[319] * v[94] + v[321] * (-(v[123] * v[138]) - v[172] * v[94])
                + v[123] * v[324] * v[95] + v[326] * (v[123] * v[140] + v[172] * v[95]) - v[123] * v[328] * v[96] + v[330] * (-
                    (v[123] * v[143]) - v[172] * v[96]);
            /* 317 = \[OverBracket]_ff_;\[Xi](\[Yen]|W)_169 */
            v[317] = -(v[122] * v[132] * v[178] * v[235]) - v[122] * v[127] * v[178] * v[245] - v[120] * v[122] * v[178] * v[252] +
                (v[168] * v[339]) / 2e0 + v[122] * v[178] * v[329] * v[94] * v[95] + v[122] * v[178] * v[325] * v[94] * v[96]
                + v[122] * v[178] * v[320] * v[95] * v[96];
            /* 316 = \[OverBracket]_ff_(\[Yen]|W)_178 */
            v[316] = (-(v[122] * v[142] * v[178]) - v[132] * v[176] * v[178]) * v[235] + (-(v[122] * v[145] * v[178])
                - v[127] * v[176] * v[178]) * v[245] + (-(v[122] * v[146] * v[178]) - v[120] * v[176] * v[178]) * v[252]
                - v[120] * v[122] * v[178] * v[255] - v[122] * v[127] * v[178] * v[265] - v[122] * v[132] * v[178] * v[270]
                + v[122] * v[178] * v[327] * v[94] * v[95] + v[329] * (v[122] * v[140] * v[178] * v[94] + v[122] * v[138] * v[178] * v[95]
                    + v[176] * v[178] * v[94] * v[95]) + v[122] * v[178] * v[323] * v[94] * v[96] + v[122] * v[178] * v[318] * v[95] * v[96]
                + v[325] * (v[122] * v[143] * v[178] * v[94] + v[122] * v[138] * v[178] * v[96] + v[176] * v[178] * v[94] * v[96]) + v[320] *
                (v[122] * v[143] * v[178] * v[95] + v[122] * v[140] * v[178] * v[96] + v[176] * v[178] * v[95] * v[96]) + (v[168] * (
                    (2e0 * v[169] * v[334]) / Power(v[121], 3) + v[543] - v[123] * (v[541] + v[235] * (-(v[132] * v[147] * v[178]) + v[142]
                        / v[97]) + v[245] * (-(v[127] * v[147] * v[178]) + v[145] / v[97]) + v[252] * (-(v[120] * v[147] * v[178]) + v[146]
                            / v[97]) + v[329] * (v[147] * v[178] * v[94] * v[95] - (v[140] * v[94]) / v[97] - (v[138] * v[95]) / v[97]) + v[325] *
                        (v[147] * v[178] * v[94] * v[96] - (v[143] * v[94]) / v[97] - (v[138] * v[96]) / v[97]) + v[320] *
                        (v[147] * v[178] * v[95] * v[96] - (v[143] * v[95]) / v[97] - (v[140] * v[96]) / v[97]) + (v[120] * v[255]) / v[97] +
                        (v[127] * v[265]) / v[97] + (v[132] * v[270]) / v[97] - (v[327] * v[94] * v[95]) / v[97] - (v[323] * v[94] * v[96]) / v[97]
                        - (v[318] * v[95] * v[96]) / v[97]))) / 2e0 - (2e0 * v[551]) / Power(v[97], 3) - (v[168] * v[545]) / (2e0 * v[97]);
        };
        /* 340 = \[OverBracket]_\[Yen]_101(ff|W) */
        v[340] = v[316];
        /* 341 = \[OverBracket]_\[Yen]_100(ff|W) */
        v[341] = v[316];
        /* 342 = \[OverBracket]_\[Yen]_110(ff|W) */
        v[342] = v[316];
        /* 343 = \[OverBracket]_\[Yen]_101(\[Yen]|W)_120 */
        v[343] = v[300] + v[340];
        /* 344 = \[OverBracket]_\[Yen]_100(\[Yen]|W)_120 */
        v[344] = v[343];
        /* 345 = \[OverBracket]_\[Yen]_101(\[Yen]|W)_127 */
        v[345] = v[282] + v[343];
        /* 346 = \[OverBracket]_\[Yen]_110(\[Yen]|W)_127 */
        v[346] = v[282] + v[342];
        /* 347 = \[OverBracket]_\[Yen]_144(\[Yen]|W)_145 */
        v[347] = v[283];
        /* 348 = \[OverBracket]_\[Yen]_139(\[Yen]|W)_145 */
        v[348] = v[283];
        /* 349 = \[OverBracket]_\[Yen]_144(\[Yen]|W)_146 */
        v[349] = v[301] + v[347];
        /* 350 = \[OverBracket]_\[Yen]_141(\[Yen]|W)_146 */
        v[350] = v[301];
        /* 351 = \[OverBracket]_\[Yen]_144(ff|W)_;\[Xi] */
        v[351] = v[317] + v[349];
        /* 352 = \[OverBracket]_\[Yen]_141(ff|W)_;\[Xi] */
        v[352] = v[317] + v[350];
        /* 353 = \[OverBracket]_\[Yen]_139(ff|W)_;\[Xi] */
        v[353] = v[317] + v[348];
        /* 313 = \[OverBracket]_u\[Phi]_6(\[Yen]|W)_144 */
        v[313] = v[313] + 2e0 * v[143] * v[351] + 2e0 * v[345] * v[96];
        /* 314 = \[OverBracket]_u\[Phi]_6;\[Xi](\[Yen]|W)_144 */
        v[314] = v[314] + 2e0 * v[351] * v[96];
        /* 354 = \[OverBracket]_peIO_3|6(u\[Phi]|W)_6 */
        v[354] = v[313] * v[45];
        /* 355 = \[OverBracket]_peIO_2|6(u\[Phi]|W)_6 */
        v[355] = v[313] * v[43];
        /* 356 = \[OverBracket]_peIO_1|6(u\[Phi]|W)_6 */
        v[356] = v[313] * v[41];
        /* 357 = \[OverBracket]_\[Yen]_100(\[Yen]|W)_132 */
        v[357] = v[272] + v[344];
        /* 358 = \[OverBracket]_\[Yen]_110(\[Yen]|W)_132 */
        v[358] = v[272] + v[346];
        /* 359 = \[OverBracket]_\[Yen]_141(\[Yen]|W)_142 */
        v[359] = v[273] + v[352];
        /* 360 = \[OverBracket]_\[Yen]_139(\[Yen]|W)_142 */
        v[360] = v[273] + v[353];
        /* 307 = \[OverBracket]_u\[Phi]_5(\[Yen]|W)_141 */
        v[307] = v[307] + 2e0 * v[140] * v[359] + 2e0 * v[357] * v[95];
        /* 308 = \[OverBracket]_u\[Phi]_5;\[Xi](\[Yen]|W)_141 */
        v[308] = v[308] + 2e0 * v[359] * v[95];
        /* 361 = \[OverBracket]_peIO_3|5(u\[Phi]|W)_5 */
        v[361] = v[307] * v[45];
        /* 362 = \[OverBracket]_peIO_2|5(u\[Phi]|W)_5 */
        v[362] = v[307] * v[43];
        /* 363 = \[OverBracket]_peIO_1|5(u\[Phi]|W)_5 */
        v[363] = v[307] * v[41];
        /* 304 = \[OverBracket]_u\[Phi]_4(\[Yen]|W)_139 */
        v[304] = v[304] + 2e0 * v[138] * v[360] + 2e0 * v[358] * v[94];
        /* 305 = \[OverBracket]_u\[Phi]_4;\[Xi](\[Yen]|W)_139 */
        v[305] = v[305] + 2e0 * v[360] * v[94];
        /* 364 = \[OverBracket]_peIO_3|4(u\[Phi]|W)_4 */
        v[364] = v[304] * v[45];
        /* 365 = \[OverBracket]_peIO_2|4(u\[Phi]|W)_4 */
        v[365] = v[304] * v[43];
        /* 366 = \[OverBracket]_peIO_1|4(u\[Phi]|W)_4 */
        v[366] = v[304] * v[41];
        /* 367 = \[OverBracket]_peIO_3|4(u\[Phi]|W)_4;\[Xi] */
        v[367] = v[364] + v[305] * v[52];
        /* 368 = \[OverBracket]_peIO_2|4(u\[Phi]|W)_4;\[Xi] */
        v[368] = v[365] + v[305] * v[50];
        /* 369 = \[OverBracket]_peIO_1|4(u\[Phi]|W)_4;\[Xi] */
        v[369] = v[366] + v[305] * v[49];
        /* 370 = \[OverBracket]_peIO_3|5(u\[Phi]|W)_5;\[Xi] */
        v[370] = v[361] + v[308] * v[52];
        /* 371 = \[OverBracket]_peIO_2|5(u\[Phi]|W)_5;\[Xi] */
        v[371] = v[362] + v[308] * v[50];
        /* 372 = \[OverBracket]_peIO_1|5(u\[Phi]|W)_5;\[Xi] */
        v[372] = v[363] + v[308] * v[49];
        /* 373 = \[OverBracket]_peIO_3|6(u\[Phi]|W)_6;\[Xi] */
        v[373] = v[354] + v[314] * v[52];
        /* 374 = \[OverBracket]_peIO_2|6(u\[Phi]|W)_6;\[Xi] */
        v[374] = v[355] + v[314] * v[50];
        /* 375 = \[OverBracket]_peIO_1|6(u\[Phi]|W)_6;\[Xi] */
        v[375] = v[356] + v[314] * v[49];
        v[1160] = v[231] * v[49];
        v[1161] = v[230] * v[49];
        v[1162] = v[228] * v[49];
        v[1163] = v[369];
        v[1164] = v[372];
        v[1165] = v[375];
        v[1166] = v[231] * v[50];
        v[1167] = v[230] * v[50];
        v[1168] = v[228] * v[50];
        v[1169] = v[368];
        v[1170] = v[371];
        v[1171] = v[374];
        v[1172] = v[231] * v[52];
        v[1173] = v[230] * v[52];
        v[1174] = v[228] * v[52];
        v[1175] = v[367];
        v[1176] = v[370];
        v[1177] = v[373];
        v[376] = 0e0;/*debug*/
        for (i218 = 1; i218 <= 18; i218++) {
            v[218] = i218;
            /* 219 = \[DoubleStruckCapitalG]_m */
            v[219] = v[1141 + i218];
            /* 377 = Rgm */
            v[377] = v[1159 + i218] * v[90];
            /* 382 = \[OverBracket]_\[OverBracket]_peIO_1|6(u\[Phi]|W)(Rgm|Rgm)_6;\[Xi] */
            v[382] = (6 == i218 ? 1 : 0) * v[90];
            /* 383 = \[OverBracket]_\[OverBracket]_peIO_2|6(u\[Phi]|W)(Rgm|Rgm)_6;\[Xi] */
            v[383] = (12 == i218 ? 1 : 0) * v[90];
            /* 384 = \[OverBracket]_\[OverBracket]_peIO_3|6(u\[Phi]|W)(Rgm|Rgm)_6;\[Xi] */
            v[384] = (18 == i218 ? 1 : 0) * v[90];
            v[409] = v[382] * v[41] + v[383] * v[43] + v[384] * v[45];
            /* 385 = \[OverBracket]_\[OverBracket]_peIO_1|5(u\[Phi]|W)(Rgm|Rgm)_5;\[Xi] */
            v[385] = (5 == i218 ? 1 : 0) * v[90];
            /* 386 = \[OverBracket]_\[OverBracket]_peIO_2|5(u\[Phi]|W)(Rgm|Rgm)_5;\[Xi] */
            v[386] = (11 == i218 ? 1 : 0) * v[90];
            /* 387 = \[OverBracket]_\[OverBracket]_peIO_3|5(u\[Phi]|W)(Rgm|Rgm)_5;\[Xi] */
            v[387] = (17 == i218 ? 1 : 0) * v[90];
            v[400] = v[385] * v[41] + v[386] * v[43] + v[387] * v[45];
            /* 388 = \[OverBracket]_\[OverBracket]_peIO_1|4(u\[Phi]|W)(Rgm|Rgm)_4;\[Xi] */
            v[388] = (4 == i218 ? 1 : 0) * v[90];
            /* 389 = \[OverBracket]_\[OverBracket]_peIO_2|4(u\[Phi]|W)(Rgm|Rgm)_4;\[Xi] */
            v[389] = (10 == i218 ? 1 : 0) * v[90];
            /* 390 = \[OverBracket]_\[OverBracket]_peIO_3|4(u\[Phi]|W)(Rgm|Rgm)_4;\[Xi] */
            v[390] = (16 == i218 ? 1 : 0) * v[90];
            v[393] = v[388] * v[41] + v[389] * v[43] + v[390] * v[45];
            /* 391 = \[OverBracket]_\[Yen]_305|3(\[Yen]|Rgm)_305|3 */
            v[391] = v[388] * v[49] + v[389] * v[50] + v[390] * v[52];
            /* 392 = [\[OverBracket]_\[OverBracket]_\[Yen]_139(\[Yen]|W)\[OverBracket]_142(u\[Phi]|Rgm)_4(\[Yen]|W)_139][\[OverBracket]_\[OverBracket]_\[Yen]_139(ff|W)\[OverBracket]_;\[Xi](\[Yen]|Rgm)_139(\[Yen]|W)_142] */
            v[392] = 2e0 * v[138] * v[393] + 2e0 * v[391] * v[94];
            /* 394 = [\[OverBracket]_\[OverBracket]_\[Yen]_110(\[Yen]|W)\[OverBracket]_132(u\[Phi]|Rgm)_4(\[Yen]|W)_139][\[OverBracket]_\[OverBracket]_\[Yen]_110(\[Yen]|W)\[OverBracket]_127(\[Yen]|Rgm)_110(\[Yen]|W)_132] */
            v[394] = 2e0 * v[393] * v[94];
            /* 395 = \[OverBracket]_u\[Phi]_4\[OverBracket]_(u\[Phi]|Rgm)_4(\[Yen]|W)_139 */
            v[395] = 2e0 * v[360] * v[391] + 2e0 * v[358] * v[393];
            v[525] = v[395];
            /* 396 = \[OverBracket]_u\[Phi]_4;\[Xi]\[OverBracket]_(u\[Phi]|Rgm)_4(\[Yen]|W)_139 */
            v[396] = 2e0 * v[360] * v[393];
            v[526] = v[396];
            /* 397 = \[OverBracket]_\[Yen]_304|3(\[Yen]|Rgm)_304|3 */
            v[397] = v[393];
            /* 398 = \[OverBracket]_\[Yen]_308|3(\[Yen]|Rgm)_308|3 */
            v[398] = v[385] * v[49] + v[386] * v[50] + v[387] * v[52];
            /* 399 = [\[OverBracket]_\[OverBracket]_\[Yen]_141(\[Yen]|W)\[OverBracket]_142(u\[Phi]|Rgm)_5(\[Yen]|W)_141][\[OverBracket]_\[OverBracket]_\[Yen]_141(ff|W)\[OverBracket]_;\[Xi](\[Yen]|Rgm)_141(\[Yen]|W)_142] */
            v[399] = 2e0 * v[140] * v[400] + 2e0 * v[398] * v[95];
            /* 401 = \[OverBracket]_\[OverBracket]_\[Yen]_100(\[Yen]|W)\[OverBracket]_132(u\[Phi]|Rgm)_5(\[Yen]|W)_141 */
            v[401] = 2e0 * v[400] * v[95];
            /* 402 = \[OverBracket]_u\[Phi]_5\[OverBracket]_(u\[Phi]|Rgm)_5(\[Yen]|W)_141 */
            v[402] = 2e0 * v[359] * v[398] + 2e0 * v[357] * v[400];
            v[524] = v[402];
            /* 403 = \[OverBracket]_u\[Phi]_5;\[Xi]\[OverBracket]_(u\[Phi]|Rgm)_5(\[Yen]|W)_141 */
            v[403] = 2e0 * v[359] * v[400];
            v[527] = v[403];
            /* 404 = \[OverBracket]_\[Yen]_307|3(\[Yen]|Rgm)_307|3 */
            v[404] = v[400];
            /* 405 = \[OverBracket]_\[Yen]_273|3\[OverBracket]_(\[Yen]|Rgm)_141(\[Yen]|W)_142 */
            v[405] = v[392] + v[399];
            /* 406 = \[OverBracket]_\[Yen]_272|3\[OverBracket]_(\[Yen]|Rgm)_100(\[Yen]|W)_132 */
            v[406] = v[394] + v[401];
            /* 407 = \[OverBracket]_\[Yen]_314|3(\[Yen]|Rgm)_314|3 */
            v[407] = v[382] * v[49] + v[383] * v[50] + v[384] * v[52];
            /* 408 = [\[OverBracket]_\[OverBracket]_\[Yen]_144(ff|W)\[OverBracket]_;\[Xi](u\[Phi]|Rgm)_6(\[Yen]|W)_144][\[OverBracket]_\[OverBracket]_\[Yen]_144(\[Yen]|W)\[OverBracket]_146(\[Yen]|Rgm)_144(ff|W)_;\[Xi]] */
            v[408] = 2e0 * v[143] * v[409] + 2e0 * v[407] * v[96];
            /* 410 = \[OverBracket]_\[OverBracket]_\[Yen]_101(\[Yen]|W)\[OverBracket]_127(u\[Phi]|Rgm)_6(\[Yen]|W)_144 */
            v[410] = 2e0 * v[409] * v[96];
            /* 411 = \[OverBracket]_u\[Phi]_6\[OverBracket]_(u\[Phi]|Rgm)_6(\[Yen]|W)_144 */
            v[411] = 2e0 * v[351] * v[407] + 2e0 * v[345] * v[409];
            v[523] = v[411];
            /* 412 = \[OverBracket]_u\[Phi]_6;\[Xi]\[OverBracket]_(u\[Phi]|Rgm)_6(\[Yen]|W)_144 */
            v[412] = 2e0 * v[351] * v[409];
            v[528] = v[412];
            /* 413 = \[OverBracket]_\[Yen]_313|3(\[Yen]|Rgm)_313|3 */
            v[413] = v[409];
            /* 414 = \[OverBracket]_\[Yen]_317|3\[OverBracket]_(\[Yen]|Rgm)_144(ff|W)_;\[Xi] */
            v[414] = v[392] + v[399] + v[408];
            /* 415 = \[OverBracket]_\[Yen]_301|3\[OverBracket]_(\[Yen]|Rgm)_144(\[Yen]|W)_146 */
            v[415] = v[399] + v[408];
            /* 416 = \[OverBracket]_\[Yen]_283|3\[OverBracket]_(\[Yen]|Rgm)_144(\[Yen]|W)_145 */
            v[416] = v[392] + v[408];
            /* 417 = [\[OverBracket]_\[OverBracket]_\[Yen]_101(\[Yen]|W)\[OverBracket]_120(\[Yen]|Rgm)_101(\[Yen]|W)_127][\[OverBracket]_\[Yen]_300|3\[OverBracket]_(\[Yen]|Rgm)_101(\[Yen]|W)_120] */
            v[417] = v[401] + v[410];
            /* 418 = \[OverBracket]_\[Yen]_282|3\[OverBracket]_(\[Yen]|Rgm)_101(\[Yen]|W)_127 */
            v[418] = v[394] + v[410];
            /* 419 = \[OverBracket]_\[Yen]_316|3\[OverBracket]_(\[Yen]|Rgm)_101(ff|W) */
            v[419] = v[394] + v[417];
            /* 420 = \[OverBracket]_\[Yen]_103(\[Yen]|Rgm)_103 */
            v[420] = 0e0;
            /* 421 = \[OverBracket]_\[Yen]_105(\[Yen]|Rgm)_105 */
            v[421] = 0e0;
            /* 422 = \[OverBracket]_\[Yen]_121(\[Yen]|Rgm)_121 */
            v[422] = 0e0;
            /* 423 = \[OverBracket]_\[Yen]_122(\[Yen]|Rgm)_122 */
            v[423] = 0e0;
            /* 424 = \[OverBracket]_\[Yen]_123(\[Yen]|Rgm)_123 */
            v[424] = 0e0;
            /* 425 = \[OverBracket]_\[Yen]_148(\[Yen]|Rgm)_148 */
            v[425] = 0e0;
            /* 426 = \[OverBracket]_\[Yen]_149(\[Yen]|Rgm)_149 */
            v[426] = 0e0;
            /* 427 = \[OverBracket]_\[Yen]_153(\[Yen]|Rgm)_153 */
            v[427] = 0e0;
            /* 428 = \[OverBracket]_\[Yen]_155(\[Yen]|Rgm)_155 */
            v[428] = 0e0;
            /* 429 = \[OverBracket]_\[Yen]_168(\[Yen]|Rgm)_168 */
            v[429] = 0e0;
            /* 430 = \[OverBracket]_\[Yen]_169(\[Yen]|Rgm)_169 */
            v[430] = 0e0;
            /* 431 = \[OverBracket]_\[Yen]_170(\[Yen]|Rgm)_170 */
            v[431] = 0e0;
            /* 432 = \[OverBracket]_\[Yen]_171(\[Yen]|Rgm)_171 */
            v[432] = 0e0;
            /* 433 = \[OverBracket]_\[Yen]_172(\[Yen]|Rgm)_172 */
            v[433] = 0e0;
            /* 434 = \[OverBracket]_\[Yen]_176(\[Yen]|Rgm)_176 */
            v[434] = 0e0;
            /* 435 = \[OverBracket]_\[Yen]_178(\[Yen]|Rgm)_178 */
            v[435] = 0e0;
            /* 436 = \[OverBracket]_\[Yen]_338(\[Yen]|Rgm)_338 */
            v[436] = 0e0;
            b437 = b98;
            v[437] = b437;
            if (b437) {
                v[443] = v[147] * v[419];
                /* 438 = \[OverBracket]_\[OverBracket]_\[Yen]_154(\[Yen]|W)\[OverBracket]_155(ff|Rgm)_(\[Yen]|W)_154 */
                v[438] = v[443] + v[414] * v[97];
                /* 439 = \[OverBracket]_\[OverBracket]_\[Yen]_105(\[Yen]|W)\[OverBracket]_152(ff|Rgm)_(\[Yen]|W)_154 */
                v[439] = v[419] * (v[148] + v[97]);
                /* 440 = \[OverBracket]_\[OverBracket]_\[Yen]_103(\[Yen]|W)\[OverBracket]_158(ff|Rgm)_(\[Yen]|W)_154 */
                v[440] = v[419] * (v[153] + v[97]);
                /* 425 = \[OverBracket]_\[Yen]_148\[OverBracket]_(ff|Rgm)_(\[Yen]|W)_154 */
                v[425] = v[310] * v[414] + v[309] * v[419];
                /* 427 = \[OverBracket]_\[Yen]_153\[OverBracket]_(ff|Rgm)_(\[Yen]|W)_154 */
                v[427] = v[303] * v[414] + v[302] * v[419];
                /* 441 = \[OverBracket]_ff_\[OverBracket]_(ff|Rgm)_(\[Yen]|W)_154 */
                v[441] = v[315] * v[414] + (v[302] + v[309]) * v[419];
                /* 442 = \[OverBracket]_ff_;\[Xi]\[OverBracket]_(ff|Rgm)_(\[Yen]|W)_154 */
                v[442] = (v[303] + v[310] + v[315]) * v[419];
                /* 444 = \[OverBracket]_\[OverBracket]_\[Yen]_149(\[Yen]|W)\[OverBracket]_152(\[Yen]|Rgm)_154(\[Yen]|W)_155 */
                v[444] = v[148] * v[414] + v[438] + v[443];
                /* 445 = \[OverBracket]_\[OverBracket]_\[Yen]_155(\[Yen]|W)\[OverBracket]_158(\[Yen]|Rgm)_154(\[Yen]|W)_155 */
                v[445] = v[153] * v[414] + v[438] + v[443];
                /* 446 = \[OverBracket]_\[OverBracket]_\[Yen]_152(v|W)\[OverBracket]_1|2;\[Xi](u\[Phi]|Rgm)_6;\[Xi](\[Yen]|W)_152 */
                v[446] = -6e0 * v[105] * v[407];
                /* 447 = \[OverBracket]_\[OverBracket]_\[Yen]_157(v|W)\[OverBracket]_1|3;\[Xi](u\[Phi]|Rgm)_6;\[Xi](\[Yen]|W)_152 */
                v[447] = v[103] * v[407] * v[94];
                /* 448 = \[OverBracket]_\[OverBracket]_\[Yen]_156(v|W)\[OverBracket]_2|3;\[Xi](u\[Phi]|Rgm)_6;\[Xi](\[Yen]|W)_152 */
                v[448] = v[103] * v[407] * v[95];
                /* 420 = \[OverBracket]_\[Yen]_103\[OverBracket]_(u\[Phi]|Rgm)_6;\[Xi](\[Yen]|W)_152 */
                v[420] = v[407] * v[516];
                /* 421 = \[OverBracket]_\[Yen]_105\[OverBracket]_(u\[Phi]|Rgm)_6;\[Xi](\[Yen]|W)_152 */
                v[421] = -6e0 * v[299] * v[407];
                /* 402 = \[OverBracket]_u\[Phi]_5\[OverBracket]_(u\[Phi]|Rgm)_6;\[Xi](\[Yen]|W)_152 */
                v[402] = v[402] + v[103] * v[278] * v[407];
                /* 395 = \[OverBracket]_u\[Phi]_4\[OverBracket]_(u\[Phi]|Rgm)_6;\[Xi](\[Yen]|W)_152 */
                v[395] = v[395] + v[103] * v[288] * v[407];
                /* 449 = \[OverBracket]_\[OverBracket]_\[Yen]_152(v|W)\[OverBracket]_1|2;\[Xi](u\[Phi]|Rgm)_6(\[Yen]|W)_152 */
                v[449] = -6e0 * v[149] * v[413] + v[446];
                /* 450 = \[OverBracket]_\[OverBracket]_\[Yen]_108(v|W)\[OverBracket]_1|2(u\[Phi]|Rgm)_6(\[Yen]|W)_152 */
                v[450] = -6e0 * v[105] * v[413];
                /* 451 = \[OverBracket]_\[OverBracket]_\[Yen]_157(v|W)\[OverBracket]_1|3;\[Xi](u\[Phi]|Rgm)_6(\[Yen]|W)_152 */
                v[451] = v[312] * v[413] + v[447];
                /* 452 = \[OverBracket]_\[OverBracket]_\[Yen]_113(v|W)\[OverBracket]_1|3(u\[Phi]|Rgm)_6(\[Yen]|W)_152 */
                v[452] = v[103] * v[413] * v[94];
                /* 453 = \[OverBracket]_\[OverBracket]_\[Yen]_156(v|W)\[OverBracket]_2|3;\[Xi](u\[Phi]|Rgm)_6(\[Yen]|W)_152 */
                v[453] = v[311] * v[413] + v[448];
                /* 454 = \[OverBracket]_\[OverBracket]_\[Yen]_116(v|W)\[OverBracket]_2|3(u\[Phi]|Rgm)_6(\[Yen]|W)_152 */
                v[454] = v[103] * v[413] * v[95];
                /* 455 = \[OverBracket]_\[Yen]_311\[OverBracket]_(u\[Phi]|Rgm)_6(\[Yen]|W)_152 */
                v[455] = v[278] * v[413];
                /* 456 = \[OverBracket]_\[Yen]_312\[OverBracket]_(u\[Phi]|Rgm)_6(\[Yen]|W)_152 */
                v[456] = v[288] * v[413];
                /* 420 = \[OverBracket]_\[Yen]_103\[OverBracket]_(u\[Phi]|Rgm)_6(\[Yen]|W)_152 */
                v[420] = v[420] + v[413] * (v[514] + v[515]);
                /* 421 = \[OverBracket]_\[Yen]_105\[OverBracket]_(u\[Phi]|Rgm)_6(\[Yen]|W)_152 */
                v[421] = -6e0 * v[295] * v[413] + v[421];
                /* 426 = \[OverBracket]_\[Yen]_149\[OverBracket]_(u\[Phi]|Rgm)_6(\[Yen]|W)_152 */
                v[426] = -6e0 * v[299] * v[413];
                /* 402 = \[OverBracket]_u\[Phi]_5\[OverBracket]_(u\[Phi]|Rgm)_6(\[Yen]|W)_152 */
                v[402] = v[402] + v[103] * v[274] * v[413];
                /* 395 = \[OverBracket]_u\[Phi]_4\[OverBracket]_(u\[Phi]|Rgm)_6(\[Yen]|W)_152 */
                v[395] = v[395] + v[103] * v[284] * v[413];
                /* 457 = \[OverBracket]_\[OverBracket]_\[Yen]_152(v|W)\[OverBracket]_1|2;\[Xi](\[Yen]|Rgm)_105(\[Yen]|W)_152 */
                v[457] = -6e0 * v[143] * v[439] + v[449] - 6e0 * v[444] * v[96];
                /* 458 = \[OverBracket]_\[OverBracket]_\[Yen]_108(v|W)\[OverBracket]_1|2(\[Yen]|Rgm)_105(\[Yen]|W)_152 */
                v[458] = v[450] - 6e0 * v[439] * v[96];
                /* 459 = \[OverBracket]_\[OverBracket]_\[Yen]_158(v|W)\[OverBracket]_1|2;\[Xi](u\[Phi]|Rgm)_5;\[Xi](\[Yen]|W)_151 */
                v[459] = v[103] * v[398] * v[94];
                /* 460 = \[OverBracket]_\[OverBracket]_\[Yen]_151(v|W)\[OverBracket]_1|3;\[Xi](u\[Phi]|Rgm)_5;\[Xi](\[Yen]|W)_151 */
                v[460] = 6e0 * v[105] * v[398] + 6e0 * v[140] * v[439] + 6e0 * v[444] * v[95];
                /* 461 = \[OverBracket]_\[OverBracket]_\[Yen]_156(v|W)\[OverBracket]_2|3;\[Xi](u\[Phi]|Rgm)_5;\[Xi](\[Yen]|W)_151 */
                v[461] = v[453] + v[103] * v[398] * v[96];
                /* 420 = \[OverBracket]_\[Yen]_103\[OverBracket]_(u\[Phi]|Rgm)_5;\[Xi](\[Yen]|W)_151 */
                v[420] = v[420] + v[398] * v[519];
                /* 421 = \[OverBracket]_\[Yen]_105\[OverBracket]_(u\[Phi]|Rgm)_5;\[Xi](\[Yen]|W)_151 */
                v[421] = 6e0 * v[291] * v[398] + v[421];
                /* 411 = \[OverBracket]_u\[Phi]_6\[OverBracket]_(u\[Phi]|Rgm)_5;\[Xi](\[Yen]|W)_151 */
                v[411] = v[103] * v[278] * v[398] + v[411] - 6e0 * v[295] * v[439] - 6e0 * v[299] * v[444];
                /* 395 = \[OverBracket]_u\[Phi]_4\[OverBracket]_(u\[Phi]|Rgm)_5;\[Xi](\[Yen]|W)_151 */
                v[395] = v[395] + v[103] * v[296] * v[398] - 6e0 * v[277] * v[439] - 6e0 * v[281] * v[444];
                /* 462 = \[OverBracket]_\[OverBracket]_\[Yen]_158(v|W)\[OverBracket]_1|2;\[Xi](u\[Phi]|Rgm)_5(\[Yen]|W)_151 */
                v[462] = v[312] * v[404] + v[459];
                /* 463 = \[OverBracket]_\[OverBracket]_\[Yen]_107(v|W)\[OverBracket]_1|2(u\[Phi]|Rgm)_5(\[Yen]|W)_151 */
                v[463] = v[103] * v[404] * v[94];
                /* 464 = \[OverBracket]_\[OverBracket]_\[Yen]_151(v|W)\[OverBracket]_1|3;\[Xi](u\[Phi]|Rgm)_5(\[Yen]|W)_151 */
                v[464] = 6e0 * v[149] * v[404] + v[460];
                /* 465 = \[OverBracket]_\[OverBracket]_\[Yen]_114(v|W)\[OverBracket]_1|3(u\[Phi]|Rgm)_5(\[Yen]|W)_151 */
                v[465] = 6e0 * v[105] * v[404] + 6e0 * v[439] * v[95];
                /* 466 = \[OverBracket]_\[OverBracket]_\[Yen]_156(v|W)\[OverBracket]_2|3;\[Xi](u\[Phi]|Rgm)_5(\[Yen]|W)_151 */
                v[466] = v[306] * v[404] + v[461];
                /* 467 = \[OverBracket]_\[OverBracket]_\[Yen]_116(v|W)\[OverBracket]_2|3(u\[Phi]|Rgm)_5(\[Yen]|W)_151 */
                v[467] = v[454] + v[103] * v[404] * v[96];
                /* 468 = \[OverBracket]_\[Yen]_306\[OverBracket]_(u\[Phi]|Rgm)_5(\[Yen]|W)_151 */
                v[468] = v[278] * v[404];
                /* 469 = \[OverBracket]_\[Yen]_312\[OverBracket]_(u\[Phi]|Rgm)_5(\[Yen]|W)_151 */
                v[469] = v[296] * v[404] + v[456];
                /* 420 = \[OverBracket]_\[Yen]_103\[OverBracket]_(u\[Phi]|Rgm)_5(\[Yen]|W)_151 */
                v[420] = v[420] + v[404] * (v[517] + v[518]);
                /* 421 = \[OverBracket]_\[Yen]_105\[OverBracket]_(u\[Phi]|Rgm)_5(\[Yen]|W)_151 */
                v[421] = 6e0 * v[287] * v[404] + v[421];
                /* 426 = \[OverBracket]_\[Yen]_149\[OverBracket]_(u\[Phi]|Rgm)_5(\[Yen]|W)_151 */
                v[426] = 6e0 * v[291] * v[404] + v[426];
                /* 411 = \[OverBracket]_u\[Phi]_6\[OverBracket]_(u\[Phi]|Rgm)_5(\[Yen]|W)_151 */
                v[411] = v[103] * v[274] * v[404] + v[411];
                /* 395 = \[OverBracket]_u\[Phi]_4\[OverBracket]_(u\[Phi]|Rgm)_5(\[Yen]|W)_151 */
                v[395] = v[395] + v[103] * v[292] * v[404];
                /* 470 = \[OverBracket]_\[OverBracket]_\[Yen]_158(v|W)\[OverBracket]_1|2;\[Xi](u\[Phi]|Rgm)_4;\[Xi](\[Yen]|W)_150 */
                v[470] = v[462] + v[103] * v[391] * v[95];
                /* 471 = \[OverBracket]_\[OverBracket]_\[Yen]_157(v|W)\[OverBracket]_1|3;\[Xi](u\[Phi]|Rgm)_4;\[Xi](\[Yen]|W)_150 */
                v[471] = v[451] + v[103] * v[391] * v[96];
                /* 472 = \[OverBracket]_\[OverBracket]_\[Yen]_150(v|W)\[OverBracket]_2|3;\[Xi](u\[Phi]|Rgm)_4;\[Xi](\[Yen]|W)_150 */
                v[472] = -6e0 * v[105] * v[391] - 6e0 * v[138] * v[439] - 6e0 * v[444] * v[94];
                /* 420 = \[OverBracket]_\[Yen]_103\[OverBracket]_(u\[Phi]|Rgm)_4;\[Xi](\[Yen]|W)_150 */
                v[420] = v[420] + v[391] * v[522];
                /* 421 = \[OverBracket]_\[Yen]_105\[OverBracket]_(u\[Phi]|Rgm)_4;\[Xi](\[Yen]|W)_150 */
                v[421] = -6e0 * v[281] * v[391] + v[421];
                /* 411 = \[OverBracket]_u\[Phi]_6\[OverBracket]_(u\[Phi]|Rgm)_4;\[Xi](\[Yen]|W)_150 */
                v[411] = v[103] * v[288] * v[391] + v[411];
                /* 402 = \[OverBracket]_u\[Phi]_5\[OverBracket]_(u\[Phi]|Rgm)_4;\[Xi](\[Yen]|W)_150 */
                v[402] = v[103] * v[296] * v[391] + v[402] + 6e0 * v[287] * v[439] + 6e0 * v[291] * v[444];
                /* 473 = \[OverBracket]_\[OverBracket]_\[Yen]_158(v|W)\[OverBracket]_1|2;\[Xi](u\[Phi]|Rgm)_4(\[Yen]|W)_150 */
                v[473] = v[311] * v[397] + v[470];
                /* 474 = \[OverBracket]_\[OverBracket]_\[Yen]_107(v|W)\[OverBracket]_1|2(u\[Phi]|Rgm)_4(\[Yen]|W)_150 */
                v[474] = v[463] + v[103] * v[397] * v[95];
                /* 475 = \[OverBracket]_\[OverBracket]_\[Yen]_157(v|W)\[OverBracket]_1|3;\[Xi](u\[Phi]|Rgm)_4(\[Yen]|W)_150 */
                v[475] = v[306] * v[397] + v[471];
                /* 476 = \[OverBracket]_\[OverBracket]_\[Yen]_113(v|W)\[OverBracket]_1|3(u\[Phi]|Rgm)_4(\[Yen]|W)_150 */
                v[476] = v[452] + v[103] * v[397] * v[96];
                /* 477 = \[OverBracket]_\[OverBracket]_\[Yen]_150(v|W)\[OverBracket]_2|3;\[Xi](u\[Phi]|Rgm)_4(\[Yen]|W)_150 */
                v[477] = -6e0 * v[149] * v[397] + v[472];
                /* 478 = \[OverBracket]_\[OverBracket]_\[Yen]_117(v|W)\[OverBracket]_2|3(u\[Phi]|Rgm)_4(\[Yen]|W)_150 */
                v[478] = -6e0 * v[105] * v[397] - 6e0 * v[439] * v[94];
                /* 479 = \[OverBracket]_\[Yen]_306\[OverBracket]_(u\[Phi]|Rgm)_4(\[Yen]|W)_150 */
                v[479] = v[288] * v[397] + v[468];
                /* 480 = \[OverBracket]_\[Yen]_311\[OverBracket]_(u\[Phi]|Rgm)_4(\[Yen]|W)_150 */
                v[480] = v[296] * v[397] + v[455];
                /* 420 = \[OverBracket]_\[Yen]_103\[OverBracket]_(u\[Phi]|Rgm)_4(\[Yen]|W)_150 */
                v[420] = v[420] + v[397] * (v[520] + v[521]);
                /* 421 = \[OverBracket]_\[Yen]_105\[OverBracket]_(u\[Phi]|Rgm)_4(\[Yen]|W)_150 */
                v[421] = -6e0 * v[277] * v[397] + v[421];
                /* 426 = \[OverBracket]_\[Yen]_149\[OverBracket]_(u\[Phi]|Rgm)_4(\[Yen]|W)_150 */
                v[426] = -6e0 * v[281] * v[397] + v[426];
                /* 411 = \[OverBracket]_u\[Phi]_6\[OverBracket]_(u\[Phi]|Rgm)_4(\[Yen]|W)_150 */
                v[411] = v[103] * v[284] * v[397] + v[411];
                /* 402 = \[OverBracket]_u\[Phi]_5\[OverBracket]_(u\[Phi]|Rgm)_4(\[Yen]|W)_150 */
                v[402] = v[103] * v[292] * v[397] + v[402];
                /* 482 = \[OverBracket]_\[OverBracket]_\[Yen]_158(v|W)\[OverBracket]_1|2;\[Xi](\[Yen]|Rgm)_103(\[Yen]|W)_158 */
                v[482] = v[473] + v[440] * v[481] + v[445] * v[94] * v[95];
                /* 483 = \[OverBracket]_\[OverBracket]_\[Yen]_107(v|W)\[OverBracket]_1|2(\[Yen]|Rgm)_103(\[Yen]|W)_158 */
                v[483] = v[474] + v[440] * v[94] * v[95];
                /* 485 = \[OverBracket]_\[OverBracket]_\[Yen]_157(v|W)\[OverBracket]_1|3;\[Xi](\[Yen]|Rgm)_103(\[Yen]|W)_158 */
                v[485] = v[475] + v[440] * v[484] + v[445] * v[94] * v[96];
                /* 486 = \[OverBracket]_\[OverBracket]_\[Yen]_113(v|W)\[OverBracket]_1|3(\[Yen]|Rgm)_103(\[Yen]|W)_158 */
                v[486] = v[476] + v[440] * v[94] * v[96];
                /* 488 = \[OverBracket]_\[OverBracket]_\[Yen]_156(v|W)\[OverBracket]_2|3;\[Xi](\[Yen]|Rgm)_103(\[Yen]|W)_158 */
                v[488] = v[466] + v[440] * v[487] + v[445] * v[95] * v[96];
                /* 489 = \[OverBracket]_\[OverBracket]_\[Yen]_116(v|W)\[OverBracket]_2|3(\[Yen]|Rgm)_103(\[Yen]|W)_158 */
                v[489] = v[467] + v[440] * v[95] * v[96];
                /* 490 = \[OverBracket]_\[Yen]_120\[OverBracket]_(\[Yen]|Rgm)_103(\[Yen]|W)_158 */
                v[490] = (-1e0 / 720e0) * (v[255] * v[440]) - (v[252] * v[445]) / 720e0;
                /* 491 = \[OverBracket]_\[Yen]_127\[OverBracket]_(\[Yen]|Rgm)_103(\[Yen]|W)_158 */
                v[491] = (-1e0 / 720e0) * (v[265] * v[440]) - (v[245] * v[445]) / 720e0;
                /* 492 = \[OverBracket]_\[Yen]_145\[OverBracket]_(\[Yen]|Rgm)_103(\[Yen]|W)_158 */
                v[492] = (-1e0 / 720e0) * (v[245] * v[440]);
                /* 493 = \[OverBracket]_\[Yen]_146\[OverBracket]_(\[Yen]|Rgm)_103(\[Yen]|W)_158 */
                v[493] = (-1e0 / 720e0) * (v[252] * v[440]);
                /* 494 = \[OverBracket]_\[Yen]_132\[OverBracket]_(\[Yen]|Rgm)_103(\[Yen]|W)_158 */
                v[494] = (-1e0 / 720e0) * (v[270] * v[440]) - (v[235] * v[445]) / 720e0;
                /* 495 = \[OverBracket]_\[Yen]_142\[OverBracket]_(\[Yen]|Rgm)_103(\[Yen]|W)_158 */
                v[495] = (-1e0 / 720e0) * (v[235] * v[440]);
                /* 496 = \[OverBracket]_\[OverBracket]_R_1|1(Ft|W)\[OverBracket]_1|3(\[Yen]|Rgm)_120(v|W)_1|1;\[Xi] */
                v[496] = (-1e0 / 720e0) * (v[103] * v[417]) - (v[120] * v[440]) / 720e0;
                /* 497 = \[OverBracket]_\[OverBracket]_R_1|1;\[Xi](Fref|W)\[OverBracket]_1|1(\[Yen]|Rgm)_120(v|W)_1|1;\[Xi] */
                v[497] = (-1e0 / 720e0) * (v[103] * v[415]) - (v[155] * v[417]) / 720e0 - (v[146] * v[440]) / 720e0 - (v[120] * v[445])
                    / 720e0;
                /* 498 = \[OverBracket]_\[OverBracket]_R_2|2(Ft|W)\[OverBracket]_2|3(\[Yen]|Rgm)_127(v|W)_2|2;\[Xi] */
                v[498] = (-1e0 / 720e0) * (v[103] * v[418]) - (v[127] * v[440]) / 720e0;
                /* 499 = \[OverBracket]_\[OverBracket]_R_2|2;\[Xi](Fref|W)\[OverBracket]_2|1(\[Yen]|Rgm)_127(v|W)_2|2;\[Xi] */
                v[499] = (-1e0 / 720e0) * (v[103] * v[416]) - (v[155] * v[418]) / 720e0 - (v[145] * v[440]) / 720e0 - (v[127] * v[445])
                    / 720e0;
                /* 500 = \[OverBracket]_\[OverBracket]_R_3|3(Ft|W)\[OverBracket]_3|3(\[Yen]|Rgm)_132(v|W)_3|3;\[Xi] */
                v[500] = (-1e0 / 720e0) * (v[103] * v[406]) - (v[132] * v[440]) / 720e0;
                /* 501 = \[OverBracket]_\[OverBracket]_R_3|3;\[Xi](Fref|W)\[OverBracket]_3|1(\[Yen]|Rgm)_132(v|W)_3|3;\[Xi] */
                v[501] = (-1e0 / 720e0) * (v[103] * v[405]) - (v[155] * v[406]) / 720e0 - (v[142] * v[440]) / 720e0 - (v[132] * v[445])
                    / 720e0;
                /* 502 = \[OverBracket]_\[OverBracket]_R_2|3(Ft|W)(\[Yen]|Rgm)_2|3275 */
                v[502] = (v[478] + v[489]) / 720e0;
                /* 503 = \[OverBracket]_\[OverBracket]_R_3|2(Ft|W)(\[Yen]|Rgm)_3|3276 */
                v[503] = (-v[478] + v[489]) / 720e0;
                /* 504 = \[OverBracket]_\[OverBracket]_R_3|2;\[Xi](Fref|W)(\[Yen]|Rgm)_3|1279 */
                v[504] = (-v[477] + v[488]) / 720e0;
                /* 505 = \[OverBracket]_\[OverBracket]_R_2|3;\[Xi](Fref|W)(\[Yen]|Rgm)_2|1280 */
                v[505] = (v[477] + v[488]) / 720e0;
                /* 506 = \[OverBracket]_\[OverBracket]_R_1|3(Ft|W)(\[Yen]|Rgm)_1|3285 */
                v[506] = (v[465] + v[486]) / 720e0;
                /* 507 = \[OverBracket]_\[OverBracket]_R_3|1(Ft|W)(\[Yen]|Rgm)_3|3286 */
                v[507] = (-v[465] + v[486]) / 720e0;
                /* 508 = \[OverBracket]_\[OverBracket]_R_3|1;\[Xi](Fref|W)(\[Yen]|Rgm)_3|1289 */
                v[508] = (-v[464] + v[485]) / 720e0;
                /* 509 = \[OverBracket]_\[OverBracket]_R_1|3;\[Xi](Fref|W)(\[Yen]|Rgm)_1|1290 */
                v[509] = (v[464] + v[485]) / 720e0;
                /* 510 = \[OverBracket]_\[OverBracket]_R_1|2(Ft|W)(\[Yen]|Rgm)_1|3293 */
                v[510] = (v[458] + v[483]) / 720e0;
                /* 511 = \[OverBracket]_\[OverBracket]_R_2|1(Ft|W)(\[Yen]|Rgm)_2|3294 */
                v[511] = (-v[458] + v[483]) / 720e0;
                /* 512 = \[OverBracket]_\[OverBracket]_R_2|1;\[Xi](Fref|W)(\[Yen]|Rgm)_2|1297 */
                v[512] = (-v[457] + v[482]) / 720e0;
                /* 513 = \[OverBracket]_\[OverBracket]_R_1|2;\[Xi](Fref|W)(\[Yen]|Rgm)_1|1298 */
                v[513] = (v[457] + v[482]) / 720e0;
                /* 411 = \[OverBracket]_u\[Phi]_6(\[Yen]|Rgm)_306 */
                v[411] = v[411] + v[155] * v[479] + v[440] * (v[140] * v[278] + v[138] * v[288] + v[514] + v[515]) + v[445] * v[516];
                /* 412 = \[OverBracket]_u\[Phi]_6;\[Xi](\[Yen]|Rgm)_306 */
                v[412] = v[412] - 6e0 * v[299] * v[439] + v[103] * v[479] + v[440] * v[516];
                /* 402 = \[OverBracket]_u\[Phi]_5(\[Yen]|Rgm)_311 */
                v[402] = v[402] + v[155] * v[480] + v[440] * (v[143] * v[278] + v[138] * v[296] + v[517] + v[518]) + v[445] * v[519];
                /* 403 = \[OverBracket]_u\[Phi]_5;\[Xi](\[Yen]|Rgm)_311 */
                v[403] = v[403] + 6e0 * v[291] * v[439] + v[103] * v[480] + v[440] * v[519];
                /* 420 = \[OverBracket]_\[Yen]_103(\[Yen]|Rgm)_312 */
                v[420] = (-1e0 / 720e0) * (v[235] * v[405]) - (v[270] * v[406]) / 720e0 - (v[252] * v[415]) / 720e0 - (v[245] * v[416])
                    / 720e0 - (v[255] * v[417]) / 720e0 - (v[265] * v[418]) / 720e0 + v[420] + v[138] * v[469] + v[143] * v[479]
                    + v[140] * v[480];
                /* 428 = \[OverBracket]_\[Yen]_155(\[Yen]|Rgm)_312 */
                v[428] = (-1e0 / 720e0) * (v[235] * v[406]) - (v[252] * v[417]) / 720e0 - (v[245] * v[418]) / 720e0 + v[469] * v[94]
                    + v[480] * v[95] + v[479] * v[96];
                /* 395 = \[OverBracket]_u\[Phi]_4(\[Yen]|Rgm)_312 */
                v[395] = v[395] + v[155] * v[469] + v[440] * (v[143] * v[288] + v[140] * v[296] + v[520] + v[521]) + v[445] * v[522];
                /* 396 = \[OverBracket]_u\[Phi]_4;\[Xi](\[Yen]|Rgm)_312 */
                v[396] = v[396] - 6e0 * v[281] * v[439] + v[103] * v[469] + v[440] * v[522];
            }
            else {
                v[622] = -((v[329] * v[95]) / v[97]) - (v[325] * v[96]) / v[97];
                v[617] = -((v[329] * v[94]) / v[97]) - (v[320] * v[96]) / v[97];
                v[614] = -((v[325] * v[94]) / v[97]) - (v[320] * v[95]) / v[97];
                v[613] = -(v[132] * v[147] * v[178]) + v[142] / v[97];
                v[611] = -(v[127] * v[147] * v[178]) + v[145] / v[97];
                v[610] = -(v[120] * v[147] * v[178]) + v[146] / v[97];
                v[605] = -((v[327] * v[95]) / v[97]);
                v[604] = -((v[323] * v[96]) / v[97]);
                v[589] = -((v[327] * v[94]) / v[97]);
                v[588] = -((v[318] * v[96]) / v[97]);
                v[578] = v[147] * v[178] * v[94] - v[138] / v[97];
                v[577] = -((v[323] * v[94]) / v[97]);
                v[576] = v[147] * v[178] * v[95] - v[140] / v[97];
                v[575] = -((v[318] * v[95]) / v[97]);
                v[564] = v[147] * v[178] * v[94] * v[95] - (v[140] * v[94]) / v[97] - (v[138] * v[95]) / v[97];
                v[553] = v[147] * v[178] * v[96] - v[143] / v[97];
                v[540] = v[147] * v[178] * v[94] * v[96] - (v[143] * v[94]) / v[97] - (v[138] * v[96]) / v[97];
                v[539] = v[147] * v[178] * v[95] * v[96] - (v[143] * v[95]) / v[97] - (v[140] * v[96]) / v[97];
                v[544] = -(v[320] * v[539]) - v[325] * v[540] - v[541] - v[329] * v[564] - v[252] * v[610] - v[245] * v[611]
                    - v[235] * v[613] - (v[120] * v[255]) / v[97] - (v[127] * v[265]) / v[97] - (v[132] * v[270]) / v[97] +
                    (v[327] * v[94] * v[95]) / v[97] + (v[323] * v[94] * v[96]) / v[97] + (v[318] * v[95] * v[96]) / v[97];
                v[537] = 1e0 / Power(v[97], 3);
                v[531] = 1e0 / Power(v[121], 3);
                /* 422 = \[OverBracket]_\[Yen]_121\[OverBracket]_(ff|Rgm)_(\[Yen]|W)_178 */
                v[422] = (-3e0 * v[168] * v[169] * v[334] * v[419]) / Power(v[121], 4);
                /* 430 = \[OverBracket]_\[Yen]_169\[OverBracket]_(ff|Rgm)_(\[Yen]|W)_178 */
                v[430] = (v[168] * v[419] * (-(v[171] * v[331]) - v[123] * v[335] + 2e0 * v[334] * v[531])) / 2e0;
                /* 492 = \[OverBracket]_\[Yen]_145\[OverBracket]_(ff|Rgm)_(\[Yen]|W)_178 */
                v[492] = v[419] * (-(v[122] * v[178] * v[245]) - (v[123] * v[168] * v[245]) / (2e0 * v[97]));
                /* 493 = \[OverBracket]_\[Yen]_146\[OverBracket]_(ff|Rgm)_(\[Yen]|W)_178 */
                v[493] = v[419] * (-(v[122] * v[178] * v[252]) - (v[123] * v[168] * v[252]) / (2e0 * v[97]));
                /* 495 = \[OverBracket]_\[Yen]_142\[OverBracket]_(ff|Rgm)_(\[Yen]|W)_178 */
                v[495] = v[419] * (-(v[122] * v[178] * v[235]) - (v[123] * v[168] * v[235]) / (2e0 * v[97]));
                /* 529 = \[OverBracket]_\[OverBracket]_\[Yen]_169(\[Yen]|W)\[OverBracket]_172(ff|Rgm)_;\[Xi](\[Yen]|W)_169 */
                v[529] = (v[168] * v[414]) / 2e0 - (v[147] * v[168] * v[419]) / (4e0 * v[97]);
                /* 530 = \[OverBracket]_\[OverBracket]_\[Yen]_172(\[Yen]|W)\[OverBracket]_175(\[Yen]|Rgm)_169(\[Yen]|W)_172 */
                v[530] = -0.5e0 * (v[123] * v[168] * v[169] * v[419]) + v[171] * v[529];
                /* 532 = \[OverBracket]_\[OverBracket]_\[Yen]_170(\[Yen]|W)\[OverBracket]_175(\[Yen]|Rgm)_169(\[Yen]|W)_172 */
                v[532] = -(v[338] * v[529]) + v[168] * v[169] * v[419] * v[531];
                /* 533 = \[OverBracket]_\[OverBracket]_\[Yen]_176(\[Yen]|W)\[OverBracket]_180(\[Yen]|Rgm)_169(\[Yen]|W)_172 */
                v[533] = -0.5e0 * (v[168] * v[169] * v[171] * v[419]) - v[123] * v[529];
                /* 432 = \[OverBracket]_\[Yen]_171\[OverBracket]_(\[Yen]|Rgm)_169(\[Yen]|W)_172 */
                v[432] = v[335] * v[529] + (v[168] * v[419] * v[534]) / 2e0;
                /* 436 = \[OverBracket]_\[Yen]_338\[OverBracket]_(\[Yen]|Rgm)_169(\[Yen]|W)_172 */
                v[436] = -(v[334] * v[529]);
                /* 535 = \[OverBracket]_\[OverBracket]_\[Yen]_175(v|W)\[OverBracket]_1|2;\[Xi](u\[Phi]|Rgm)_6;\[Xi](\[Yen]|W)_175 */
                v[535] = -(v[123] * v[168] * v[407]) + v[419] * ((v[168] * v[171] * v[752]) / 2e0 - (v[168] * v[763]) / (2e0 * v[97]));
                /* 536 = \[OverBracket]_\[OverBracket]_\[Yen]_179(v|W)\[OverBracket]_1|3;\[Xi](u\[Phi]|Rgm)_6;\[Xi](\[Yen]|W)_175 */
                v[536] = v[122] * v[178] * v[414] * v[94] * v[96] + v[419] * (-0.5e0 * (v[123] * v[168] * v[540]) + v[755] + v[756]
                    + v[757] - 2e0 * v[122] * v[147] * v[537] * v[94] * v[96]) - (v[122] * v[407] * v[94]) / v[97];
                /* 538 = \[OverBracket]_\[OverBracket]_\[Yen]_177(v|W)\[OverBracket]_2|3;\[Xi](u\[Phi]|Rgm)_6;\[Xi](\[Yen]|W)_175 */
                v[538] = v[122] * v[178] * v[414] * v[95] * v[96] + v[419] * (-0.5e0 * (v[123] * v[168] * v[539]) + v[753] + v[754]
                    + v[122] * v[140] * v[178] * v[96] - 2e0 * v[122] * v[147] * v[537] * v[95] * v[96]) - (v[122] * v[407] * v[95]) / v[97];
                /* 423 = \[OverBracket]_\[Yen]_122\[OverBracket]_(u\[Phi]|Rgm)_6;\[Xi](\[Yen]|W)_175 */
                v[423] = v[407] * v[614] + v[414] * v[619] + v[419] * (-(v[142] * v[178] * v[235]) - v[145] * v[178] * v[245]
                    - v[146] * v[178] * v[252] - v[120] * v[178] * v[255] - v[127] * v[178] * v[265] - v[132] * v[178] * v[270] + v[546]
                    + v[547] + v[548] + v[549] + v[550] - 2e0 * v[537] * v[620] + v[178] * v[318] * v[95] * v[96]);
                /* 424 = \[OverBracket]_\[Yen]_123\[OverBracket]_(u\[Phi]|Rgm)_6;\[Xi](\[Yen]|W)_175 */
                v[424] = -(v[168] * v[330] * v[407]) - v[331] * v[529] + v[419] * ((v[168] * v[544]) / 2e0 - (v[168] * v[542]) /
                    (2e0 * v[97]));
                /* 429 = \[OverBracket]_\[Yen]_168\[OverBracket]_(u\[Phi]|Rgm)_6;\[Xi](\[Yen]|W)_175 */
                v[429] = -(v[123] * v[330] * v[407]) + (v[339] * v[414]) / 2e0 + v[419] * ((v[168] * v[171] * v[542]) / 2e0 +
                    (2e0 * v[169] * v[334] * v[531] + v[543] + v[123] * v[544]) / 2e0 - v[545] / (2e0 * v[97]));
                /* 441 = \[OverBracket]_ff_\[OverBracket]_(u\[Phi]|Rgm)_6;\[Xi](\[Yen]|W)_175 */
                v[441] = v[407] * v[570] + v[419] * ((v[168] * v[178] * v[545]) / 2e0 - (v[123] * v[168] * (-(v[142] * v[178] * v[235])
                    - v[145] * v[178] * v[245] - v[146] * v[178] * v[252] - v[120] * v[178] * v[255] - v[127] * v[178] * v[265]
                    - v[132] * v[178] * v[270] + v[546] + v[547] + v[548] + v[549] + v[550] + v[178] * v[318] * v[95] * v[96])) / 2e0 +
                    (6e0 * v[551]) / Power(v[97], 4));
                /* 402 = \[OverBracket]_u\[Phi]_5\[OverBracket]_(u\[Phi]|Rgm)_6;\[Xi](\[Yen]|W)_175 */
                v[402] = v[524] + v[414] * (v[569] + v[122] * v[178] * v[320] * v[96]) + v[419] * (v[320] * v[552] + v[329] * v[574] +
                    (v[168] * (v[171] * (v[168] * v[324] + v[170] * v[326]) - v[123] * (v[320] * v[553] + v[329] * v[578] + v[588] + v[589])
                        )) / 2e0 + v[590] + v[122] * v[178] * v[318] * v[96] - 2e0 * v[537] * (v[122] * v[147] * v[329] * v[94]
                            + v[122] * v[147] * v[320] * v[96]) - (v[168] * (v[123] * v[324] + v[172] * v[326])) / (2e0 * v[97])) -
                    (v[122] * v[320] * v[407]) / v[97];
                /* 395 = \[OverBracket]_u\[Phi]_4\[OverBracket]_(u\[Phi]|Rgm)_6;\[Xi](\[Yen]|W)_175 */
                v[395] = v[525] + v[414] * (v[594] + v[122] * v[178] * v[325] * v[96]) + v[419] * (v[325] * v[552] + v[329] * v[572] +
                    (v[168] * (v[171] * (-(v[168] * v[319]) - v[170] * v[321]) - v[123] * (v[325] * v[553] + v[329] * v[576] + v[604]
                        + v[605]))) / 2e0 + v[606] + v[122] * v[178] * v[323] * v[96] - 2e0 * v[537] * (v[122] * v[147] * v[329] * v[95]
                            + v[122] * v[147] * v[325] * v[96]) - (v[168] * (-(v[123] * v[319]) - v[172] * v[321])) / (2e0 * v[97])) -
                    (v[122] * v[325] * v[407]) / v[97];
                /* 554 = \[OverBracket]_\[OverBracket]_\[Yen]_175(v|W)\[OverBracket]_1|2;\[Xi](u\[Phi]|Rgm)_6(\[Yen]|W)_175 */
                v[554] = v[333] * v[413] + v[535];
                /* 555 = \[OverBracket]_\[OverBracket]_\[Yen]_125(v|W)\[OverBracket]_1|2(u\[Phi]|Rgm)_6(\[Yen]|W)_175 */
                v[555] = -(v[123] * v[168] * v[413]) + v[419] * (-0.5e0 * (v[171] * v[338] * v[96]) + (v[123] * v[168] * v[96]) /
                    (2e0 * v[97]));
                /* 556 = \[OverBracket]_\[OverBracket]_\[Yen]_179(v|W)\[OverBracket]_1|3;\[Xi](u\[Phi]|Rgm)_6(\[Yen]|W)_175 */
                v[556] = v[337] * v[413] + v[536];
                /* 557 = \[OverBracket]_\[OverBracket]_\[Yen]_129(v|W)\[OverBracket]_1|3(u\[Phi]|Rgm)_6(\[Yen]|W)_175 */
                v[557] = v[419] * (v[122] * v[178] * v[94] * v[96] + (v[123] * v[168] * v[94] * v[96]) / (2e0 * v[97])) -
                    (v[122] * v[413] * v[94]) / v[97];
                /* 558 = \[OverBracket]_\[OverBracket]_\[Yen]_177(v|W)\[OverBracket]_2|3;\[Xi](u\[Phi]|Rgm)_6(\[Yen]|W)_175 */
                v[558] = v[336] * v[413] + v[538];
                /* 559 = \[OverBracket]_\[OverBracket]_\[Yen]_131(v|W)\[OverBracket]_2|3(u\[Phi]|Rgm)_6(\[Yen]|W)_175 */
                v[559] = v[419] * (v[122] * v[178] * v[95] * v[96] + (v[123] * v[168] * v[95] * v[96]) / (2e0 * v[97])) -
                    (v[122] * v[413] * v[95]) / v[97];
                /* 560 = \[OverBracket]_\[Yen]_333\[OverBracket]_(u\[Phi]|Rgm)_6(\[Yen]|W)_175 */
                v[560] = v[330] * v[413];
                /* 561 = \[OverBracket]_\[Yen]_336\[OverBracket]_(u\[Phi]|Rgm)_6(\[Yen]|W)_175 */
                v[561] = v[320] * v[413];
                /* 562 = \[OverBracket]_\[Yen]_337\[OverBracket]_(u\[Phi]|Rgm)_6(\[Yen]|W)_175 */
                v[562] = v[325] * v[413];
                /* 423 = \[OverBracket]_\[Yen]_122\[OverBracket]_(u\[Phi]|Rgm)_6(\[Yen]|W)_175 */
                v[423] = v[423] + v[413] * (v[575] + v[577]);
                /* 424 = \[OverBracket]_\[Yen]_123\[OverBracket]_(u\[Phi]|Rgm)_6(\[Yen]|W)_175 */
                v[424] = -(v[168] * v[328] * v[413]) + v[424];
                /* 429 = \[OverBracket]_\[Yen]_168\[OverBracket]_(u\[Phi]|Rgm)_6(\[Yen]|W)_175 */
                v[429] = -(v[123] * v[328] * v[413]) + v[429];
                /* 441 = \[OverBracket]_ff_\[OverBracket]_(u\[Phi]|Rgm)_6(\[Yen]|W)_175 */
                v[441] = v[441] + v[413] * (v[571] + v[573]);
                /* 402 = \[OverBracket]_u\[Phi]_5\[OverBracket]_(u\[Phi]|Rgm)_6(\[Yen]|W)_175 */
                v[402] = v[402] - (v[122] * v[318] * v[413]) / v[97];
                /* 395 = \[OverBracket]_u\[Phi]_4\[OverBracket]_(u\[Phi]|Rgm)_6(\[Yen]|W)_175 */
                v[395] = v[395] - (v[122] * v[323] * v[413]) / v[97];
                /* 563 = \[OverBracket]_\[OverBracket]_\[Yen]_175(v|W)\[OverBracket]_1|2;\[Xi](\[Yen]|Rgm)_170(\[Yen]|W)_175 */
                v[563] = v[554] - v[168] * v[530] * v[96] - v[123] * v[532] * v[96];
                /* 565 = \[OverBracket]_\[OverBracket]_\[Yen]_180(v|W)\[OverBracket]_1|2;\[Xi](u\[Phi]|Rgm)_5;\[Xi](\[Yen]|W)_174 */
                v[565] = v[122] * v[178] * v[414] * v[94] * v[95] + v[419] * (-0.5e0 * (v[123] * v[168] * v[564]) + v[758] + v[759]
                    + v[760] - 2e0 * v[122] * v[147] * v[537] * v[94] * v[95]) - (v[122] * v[398] * v[94]) / v[97];
                /* 566 = \[OverBracket]_\[OverBracket]_\[Yen]_174(v|W)\[OverBracket]_1|3;\[Xi](u\[Phi]|Rgm)_5;\[Xi](\[Yen]|W)_174 */
                v[566] = v[123] * v[168] * v[398] + v[168] * v[530] * v[95] + v[123] * v[532] * v[95] + v[419] * (
                    (v[168] * v[171] * v[751]) / 2e0 - (v[168] * v[762]) / (2e0 * v[97]));
                /* 567 = \[OverBracket]_\[OverBracket]_\[Yen]_177(v|W)\[OverBracket]_2|3;\[Xi](u\[Phi]|Rgm)_5;\[Xi](\[Yen]|W)_174 */
                v[567] = v[558] - (v[122] * v[398] * v[96]) / v[97];
                /* 423 = \[OverBracket]_\[Yen]_122\[OverBracket]_(u\[Phi]|Rgm)_5;\[Xi](\[Yen]|W)_174 */
                v[423] = v[423] + v[398] * v[617];
                /* 424 = \[OverBracket]_\[Yen]_123\[OverBracket]_(u\[Phi]|Rgm)_5;\[Xi](\[Yen]|W)_174 */
                v[424] = v[168] * v[326] * v[398] + v[424] + v[532] * v[568];
                /* 429 = \[OverBracket]_\[Yen]_168\[OverBracket]_(u\[Phi]|Rgm)_5;\[Xi](\[Yen]|W)_174 */
                v[429] = v[123] * v[326] * v[398] + v[429] + v[530] * v[568];
                /* 441 = \[OverBracket]_ff_\[OverBracket]_(u\[Phi]|Rgm)_5;\[Xi](\[Yen]|W)_174 */
                v[441] = v[441] + v[398] * (v[569] + v[618]);
                /* 411 = \[OverBracket]_u\[Phi]_6\[OverBracket]_(u\[Phi]|Rgm)_5;\[Xi](\[Yen]|W)_174 */
                v[411] = v[523] - v[168] * v[330] * v[530] - v[123] * v[330] * v[532] + v[414] * v[570] + v[419] * (v[571]
                    + v[320] * v[572] + v[573] + v[325] * v[574] + (v[168] * (v[171] * (-(v[168] * v[328]) - v[170] * v[330]) - v[123] *
                        (v[575] + v[320] * v[576] + v[577] + v[325] * v[578]))) / 2e0 - 2e0 * v[537] * (v[122] * v[147] * v[325] * v[94]
                            + v[122] * v[147] * v[320] * v[95]) - (v[168] * (-(v[123] * v[328]) - v[172] * v[330])) / (2e0 * v[97])) -
                    (v[122] * v[320] * v[398]) / v[97];
                /* 395 = \[OverBracket]_u\[Phi]_4\[OverBracket]_(u\[Phi]|Rgm)_5;\[Xi](\[Yen]|W)_174 */
                v[395] = v[395] - v[168] * v[321] * v[530] - v[123] * v[321] * v[532] - (v[122] * v[329] * v[398]) / v[97];
                /* 579 = \[OverBracket]_\[OverBracket]_\[Yen]_180(v|W)\[OverBracket]_1|2;\[Xi](u\[Phi]|Rgm)_5(\[Yen]|W)_174 */
                v[579] = v[337] * v[404] + v[565];
                /* 580 = \[OverBracket]_\[OverBracket]_\[Yen]_126(v|W)\[OverBracket]_1|2(u\[Phi]|Rgm)_5(\[Yen]|W)_174 */
                v[580] = v[419] * (v[122] * v[178] * v[94] * v[95] + (v[123] * v[168] * v[94] * v[95]) / (2e0 * v[97])) -
                    (v[122] * v[404] * v[94]) / v[97];
                /* 581 = \[OverBracket]_\[OverBracket]_\[Yen]_174(v|W)\[OverBracket]_1|3;\[Xi](u\[Phi]|Rgm)_5(\[Yen]|W)_174 */
                v[581] = -(v[333] * v[404]) + v[566];
                /* 582 = \[OverBracket]_\[OverBracket]_\[Yen]_128(v|W)\[OverBracket]_1|3(u\[Phi]|Rgm)_5(\[Yen]|W)_174 */
                v[582] = v[123] * v[168] * v[404] + v[419] * ((v[171] * v[338] * v[95]) / 2e0 - (v[123] * v[168] * v[95]) / (2e0 * v[97])
                    );
                /* 583 = \[OverBracket]_\[OverBracket]_\[Yen]_177(v|W)\[OverBracket]_2|3;\[Xi](u\[Phi]|Rgm)_5(\[Yen]|W)_174 */
                v[583] = v[332] * v[404] + v[567];
                /* 584 = \[OverBracket]_\[OverBracket]_\[Yen]_131(v|W)\[OverBracket]_2|3(u\[Phi]|Rgm)_5(\[Yen]|W)_174 */
                v[584] = v[559] - (v[122] * v[404] * v[96]) / v[97];
                /* 585 = \[OverBracket]_\[Yen]_332\[OverBracket]_(u\[Phi]|Rgm)_5(\[Yen]|W)_174 */
                v[585] = v[320] * v[404];
                /* 586 = \[OverBracket]_\[Yen]_333\[OverBracket]_(u\[Phi]|Rgm)_5(\[Yen]|W)_174 */
                v[586] = -(v[326] * v[404]) + v[560];
                /* 587 = \[OverBracket]_\[Yen]_337\[OverBracket]_(u\[Phi]|Rgm)_5(\[Yen]|W)_174 */
                v[587] = v[329] * v[404] + v[562];
                /* 423 = \[OverBracket]_\[Yen]_122\[OverBracket]_(u\[Phi]|Rgm)_5(\[Yen]|W)_174 */
                v[423] = v[423] + v[404] * (v[588] + v[589]);
                /* 424 = \[OverBracket]_\[Yen]_123\[OverBracket]_(u\[Phi]|Rgm)_5(\[Yen]|W)_174 */
                v[424] = v[168] * v[324] * v[404] + v[424];
                /* 429 = \[OverBracket]_\[Yen]_168\[OverBracket]_(u\[Phi]|Rgm)_5(\[Yen]|W)_174 */
                v[429] = v[123] * v[324] * v[404] + v[429];
                /* 441 = \[OverBracket]_ff_\[OverBracket]_(u\[Phi]|Rgm)_5(\[Yen]|W)_174 */
                v[441] = v[441] + v[404] * (v[590] + v[122] * v[178] * v[318] * v[96]);
                /* 411 = \[OverBracket]_u\[Phi]_6\[OverBracket]_(u\[Phi]|Rgm)_5(\[Yen]|W)_174 */
                v[411] = v[411] - (v[122] * v[318] * v[404]) / v[97];
                /* 395 = \[OverBracket]_u\[Phi]_4\[OverBracket]_(u\[Phi]|Rgm)_5(\[Yen]|W)_174 */
                v[395] = v[395] - (v[122] * v[327] * v[404]) / v[97];
                /* 591 = \[OverBracket]_\[OverBracket]_\[Yen]_180(v|W)\[OverBracket]_1|2;\[Xi](u\[Phi]|Rgm)_4;\[Xi](\[Yen]|W)_173 */
                v[591] = v[579] - (v[122] * v[391] * v[95]) / v[97];
                /* 592 = \[OverBracket]_\[OverBracket]_\[Yen]_179(v|W)\[OverBracket]_1|3;\[Xi](u\[Phi]|Rgm)_4;\[Xi](\[Yen]|W)_173 */
                v[592] = v[556] - (v[122] * v[391] * v[96]) / v[97];
                /* 593 = \[OverBracket]_\[OverBracket]_\[Yen]_173(v|W)\[OverBracket]_2|3;\[Xi](u\[Phi]|Rgm)_4;\[Xi](\[Yen]|W)_173 */
                v[593] = -(v[123] * v[168] * v[391]) - v[168] * v[530] * v[94] - v[123] * v[532] * v[94] + v[419] * (
                    (v[168] * v[171] * v[750]) / 2e0 - (v[168] * v[761]) / (2e0 * v[97]));
                /* 423 = \[OverBracket]_\[Yen]_122\[OverBracket]_(u\[Phi]|Rgm)_4;\[Xi](\[Yen]|W)_173 */
                v[423] = v[423] + v[391] * v[622];
                /* 424 = \[OverBracket]_\[Yen]_123\[OverBracket]_(u\[Phi]|Rgm)_4;\[Xi](\[Yen]|W)_173 */
                v[424] = -(v[168] * v[321] * v[391]) + v[424];
                /* 429 = \[OverBracket]_\[Yen]_168\[OverBracket]_(u\[Phi]|Rgm)_4;\[Xi](\[Yen]|W)_173 */
                v[429] = -(v[123] * v[321] * v[391]) + v[429];
                /* 441 = \[OverBracket]_ff_\[OverBracket]_(u\[Phi]|Rgm)_4;\[Xi](\[Yen]|W)_173 */
                v[441] = v[441] + v[391] * (v[594] + v[623]);
                /* 411 = \[OverBracket]_u\[Phi]_6\[OverBracket]_(u\[Phi]|Rgm)_4;\[Xi](\[Yen]|W)_173 */
                v[411] = v[411] - (v[122] * v[325] * v[391]) / v[97];
                /* 402 = \[OverBracket]_u\[Phi]_5\[OverBracket]_(u\[Phi]|Rgm)_4;\[Xi](\[Yen]|W)_173 */
                v[402] = v[402] + v[168] * v[326] * v[530] + v[123] * v[326] * v[532] - (v[122] * v[329] * v[391]) / v[97];
                /* 595 = \[OverBracket]_\[OverBracket]_\[Yen]_180(v|W)\[OverBracket]_1|2;\[Xi](u\[Phi]|Rgm)_4(\[Yen]|W)_173 */
                v[595] = v[336] * v[397] + v[591];
                /* 596 = \[OverBracket]_\[OverBracket]_\[Yen]_126(v|W)\[OverBracket]_1|2(u\[Phi]|Rgm)_4(\[Yen]|W)_173 */
                v[596] = v[580] - (v[122] * v[397] * v[95]) / v[97];
                /* 597 = \[OverBracket]_\[OverBracket]_\[Yen]_179(v|W)\[OverBracket]_1|3;\[Xi](u\[Phi]|Rgm)_4(\[Yen]|W)_173 */
                v[597] = v[332] * v[397] + v[592];
                /* 598 = \[OverBracket]_\[OverBracket]_\[Yen]_129(v|W)\[OverBracket]_1|3(u\[Phi]|Rgm)_4(\[Yen]|W)_173 */
                v[598] = v[557] - (v[122] * v[397] * v[96]) / v[97];
                /* 599 = \[OverBracket]_\[OverBracket]_\[Yen]_173(v|W)\[OverBracket]_2|3;\[Xi](u\[Phi]|Rgm)_4(\[Yen]|W)_173 */
                v[599] = v[333] * v[397] + v[593];
                /* 600 = \[OverBracket]_\[OverBracket]_\[Yen]_130(v|W)\[OverBracket]_2|3(u\[Phi]|Rgm)_4(\[Yen]|W)_173 */
                v[600] = -(v[123] * v[168] * v[397]) + v[419] * (-0.5e0 * (v[171] * v[338] * v[94]) + (v[123] * v[168] * v[94]) /
                    (2e0 * v[97]));
                /* 601 = \[OverBracket]_\[Yen]_332\[OverBracket]_(u\[Phi]|Rgm)_4(\[Yen]|W)_173 */
                v[601] = v[325] * v[397] + v[585];
                /* 602 = \[OverBracket]_\[Yen]_333\[OverBracket]_(u\[Phi]|Rgm)_4(\[Yen]|W)_173 */
                v[602] = v[321] * v[397] + v[586];
                /* 603 = \[OverBracket]_\[Yen]_336\[OverBracket]_(u\[Phi]|Rgm)_4(\[Yen]|W)_173 */
                v[603] = v[329] * v[397] + v[561];
                /* 423 = \[OverBracket]_\[Yen]_122\[OverBracket]_(u\[Phi]|Rgm)_4(\[Yen]|W)_173 */
                v[423] = v[423] + v[397] * (v[604] + v[605]);
                /* 424 = \[OverBracket]_\[Yen]_123\[OverBracket]_(u\[Phi]|Rgm)_4(\[Yen]|W)_173 */
                v[424] = -(v[168] * v[319] * v[397]) + v[424];
                /* 429 = \[OverBracket]_\[Yen]_168\[OverBracket]_(u\[Phi]|Rgm)_4(\[Yen]|W)_173 */
                v[429] = -(v[123] * v[319] * v[397]) + v[429];
                /* 441 = \[OverBracket]_ff_\[OverBracket]_(u\[Phi]|Rgm)_4(\[Yen]|W)_173 */
                v[441] = v[441] + v[397] * (v[606] + v[122] * v[178] * v[323] * v[96]);
                /* 411 = \[OverBracket]_u\[Phi]_6\[OverBracket]_(u\[Phi]|Rgm)_4(\[Yen]|W)_173 */
                v[411] = v[411] - (v[122] * v[323] * v[397]) / v[97];
                /* 402 = \[OverBracket]_u\[Phi]_5\[OverBracket]_(u\[Phi]|Rgm)_4(\[Yen]|W)_173 */
                v[402] = v[402] - (v[122] * v[327] * v[397]) / v[97];
                /* 607 = \[OverBracket]_\[OverBracket]_\[Yen]_180(v|W)\[OverBracket]_1|2;\[Xi](\[Yen]|Rgm)_176(\[Yen]|W)_180 */
                v[607] = v[595] - (v[533] * v[94] * v[95]) / v[97];
                /* 608 = \[OverBracket]_\[OverBracket]_\[Yen]_179(v|W)\[OverBracket]_1|3;\[Xi](\[Yen]|Rgm)_176(\[Yen]|W)_180 */
                v[608] = v[597] - (v[533] * v[94] * v[96]) / v[97];
                /* 609 = \[OverBracket]_\[OverBracket]_\[Yen]_177(v|W)\[OverBracket]_2|3;\[Xi](\[Yen]|Rgm)_176(\[Yen]|W)_180 */
                v[609] = v[583] - (v[533] * v[95] * v[96]) / v[97];
                /* 490 = \[OverBracket]_\[Yen]_120\[OverBracket]_(\[Yen]|Rgm)_176(\[Yen]|W)_180 */
                v[490] = -(v[122] * v[178] * v[252] * v[414]) + v[419] * (-(v[176] * v[178] * v[252]) - v[122] * v[178] * v[255]
                    + 2e0 * v[122] * v[147] * v[252] * v[537] - (v[123] * v[168] * (-(v[147] * v[178] * v[252]) + v[255] / v[97])) / 2e0) +
                    (v[252] * v[533]) / v[97];
                /* 491 = \[OverBracket]_\[Yen]_127\[OverBracket]_(\[Yen]|Rgm)_176(\[Yen]|W)_180 */
                v[491] = -(v[122] * v[178] * v[245] * v[414]) + v[419] * (-(v[176] * v[178] * v[245]) - v[122] * v[178] * v[265]
                    + 2e0 * v[122] * v[147] * v[245] * v[537] - (v[123] * v[168] * (-(v[147] * v[178] * v[245]) + v[265] / v[97])) / 2e0) +
                    (v[245] * v[533]) / v[97];
                /* 494 = \[OverBracket]_\[Yen]_132\[OverBracket]_(\[Yen]|Rgm)_176(\[Yen]|W)_180 */
                v[494] = -(v[122] * v[178] * v[235] * v[414]) + v[419] * (-(v[176] * v[178] * v[235]) - v[122] * v[178] * v[270]
                    + 2e0 * v[122] * v[147] * v[235] * v[537] - (v[123] * v[168] * (-(v[147] * v[178] * v[235]) + v[270] / v[97])) / 2e0) +
                    (v[235] * v[533]) / v[97];
                /* 496 = \[OverBracket]_\[OverBracket]_R_1|1(Ft|W)\[OverBracket]_1|3(\[Yen]|Rgm)_120(v|W)_1|1;\[Xi] */
                v[496] = v[419] * (-(v[120] * v[122] * v[178]) - (v[120] * v[123] * v[168]) / (2e0 * v[97])) + (v[122] * v[417])
                    / v[97];
                /* 497 = \[OverBracket]_\[OverBracket]_R_1|1;\[Xi](Fref|W)\[OverBracket]_1|1(\[Yen]|Rgm)_120(v|W)_1|1;\[Xi] */
                v[497] = -(v[120] * v[122] * v[178] * v[414]) + v[322] * v[417] + v[419] * (-(v[122] * v[146] * v[178])
                    - v[120] * v[176] * v[178] + 2e0 * v[120] * v[122] * v[147] * v[537] - (v[123] * v[168] * v[610]) / 2e0) + (v[122] * v[415]
                        ) / v[97] + (v[120] * v[533]) / v[97];
                /* 513 = \[OverBracket]_\[OverBracket]_R_1|2;\[Xi](Fref|W)\[OverBracket]_1|1(\[Yen]|Rgm)_180(v|W)_1|2;\[Xi] */
                v[513] = v[563] + v[607];
                /* 512 = \[OverBracket]_\[OverBracket]_R_2|1;\[Xi](Fref|W)\[OverBracket]_2|1(\[Yen]|Rgm)_180(v|W)_1|2;\[Xi] */
                v[512] = -v[563] + v[607];
                /* 511 = \[OverBracket]_\[OverBracket]_R_2|1(Ft|W)\[OverBracket]_2|3(\[Yen]|Rgm)_126(v|W)_1|2 */
                v[511] = -v[555] + v[596];
                /* 510 = \[OverBracket]_\[OverBracket]_R_1|2(Ft|W)\[OverBracket]_1|3(\[Yen]|Rgm)_126(v|W)_1|2 */
                v[510] = v[555] + v[596];
                /* 509 = \[OverBracket]_\[OverBracket]_R_1|3;\[Xi](Fref|W)\[OverBracket]_1|1(\[Yen]|Rgm)_179(v|W)_1|3;\[Xi] */
                v[509] = v[581] + v[608];
                /* 508 = \[OverBracket]_\[OverBracket]_R_3|1;\[Xi](Fref|W)\[OverBracket]_3|1(\[Yen]|Rgm)_179(v|W)_1|3;\[Xi] */
                v[508] = -v[581] + v[608];
                /* 507 = \[OverBracket]_\[OverBracket]_R_3|1(Ft|W)\[OverBracket]_3|3(\[Yen]|Rgm)_129(v|W)_1|3 */
                v[507] = -v[582] + v[598];
                /* 506 = \[OverBracket]_\[OverBracket]_R_1|3(Ft|W)\[OverBracket]_1|3(\[Yen]|Rgm)_129(v|W)_1|3 */
                v[506] = v[582] + v[598];
                /* 498 = \[OverBracket]_\[OverBracket]_R_2|2(Ft|W)\[OverBracket]_2|3(\[Yen]|Rgm)_127(v|W)_2|2;\[Xi] */
                v[498] = v[419] * (-(v[122] * v[127] * v[178]) - (v[123] * v[127] * v[168]) / (2e0 * v[97])) + (v[122] * v[418])
                    / v[97];
                /* 499 = \[OverBracket]_\[OverBracket]_R_2|2;\[Xi](Fref|W)\[OverBracket]_2|1(\[Yen]|Rgm)_127(v|W)_2|2;\[Xi] */
                v[499] = -(v[122] * v[127] * v[178] * v[414]) + v[322] * v[418] + v[419] * (-(v[122] * v[145] * v[178])
                    - v[127] * v[176] * v[178] + 2e0 * v[122] * v[127] * v[147] * v[537] - (v[123] * v[168] * v[611]) / 2e0) + (v[122] * v[416]
                        ) / v[97] + (v[127] * v[533]) / v[97];
                /* 505 = \[OverBracket]_\[OverBracket]_R_2|3;\[Xi](Fref|W)\[OverBracket]_2|1(\[Yen]|Rgm)_177(v|W)_2|3;\[Xi] */
                v[505] = v[599] + v[609];
                /* 504 = \[OverBracket]_\[OverBracket]_R_3|2;\[Xi](Fref|W)\[OverBracket]_3|1(\[Yen]|Rgm)_177(v|W)_2|3;\[Xi] */
                v[504] = -v[599] + v[609];
                /* 503 = \[OverBracket]_\[OverBracket]_R_3|2(Ft|W)\[OverBracket]_3|3(\[Yen]|Rgm)_131(v|W)_2|3 */
                v[503] = v[584] - v[600];
                /* 502 = \[OverBracket]_\[OverBracket]_R_2|3(Ft|W)\[OverBracket]_2|3(\[Yen]|Rgm)_131(v|W)_2|3 */
                v[502] = v[584] + v[600];
                /* 612 = \[OverBracket]_\[Yen]_322\[OverBracket]_(\[Yen]|Rgm)_132(v|W)_3|3;\[Xi] */
                v[612] = v[235] * v[406] + v[252] * v[417] + v[245] * v[418];
                /* 500 = \[OverBracket]_\[OverBracket]_R_3|3(Ft|W)\[OverBracket]_3|3(\[Yen]|Rgm)_132(v|W)_3|3;\[Xi] */
                v[500] = v[419] * (-(v[122] * v[132] * v[178]) - (v[123] * v[132] * v[168]) / (2e0 * v[97])) + (v[122] * v[406])
                    / v[97];
                /* 501 = \[OverBracket]_\[OverBracket]_R_3|3;\[Xi](Fref|W)\[OverBracket]_3|1(\[Yen]|Rgm)_132(v|W)_3|3;\[Xi] */
                v[501] = v[322] * v[406] - v[122] * v[132] * v[178] * v[414] + v[419] * (-(v[122] * v[142] * v[178])
                    - v[132] * v[176] * v[178] + 2e0 * v[122] * v[132] * v[147] * v[537] - (v[123] * v[168] * v[613]) / 2e0) + (v[122] * v[405]
                        ) / v[97] + (v[132] * v[533]) / v[97];
                /* 411 = \[OverBracket]_u\[Phi]_6(\[Yen]|Rgm)_332 */
                v[411] = v[411] - v[322] * v[601] + v[533] * v[614];
                /* 412 = \[OverBracket]_u\[Phi]_6;\[Xi](\[Yen]|Rgm)_332 */
                v[412] = v[528] + v[419] * ((v[168] * (-(v[168] * v[171] * v[330]) - v[123] * v[614])) / 2e0 + v[615] + v[616] +
                    (v[123] * v[168] * v[330]) / (2e0 * v[97])) - (v[122] * v[601]) / v[97];
                /* 424 = \[OverBracket]_\[Yen]_123(\[Yen]|Rgm)_333 */
                v[424] = v[424] - v[170] * v[602];
                /* 431 = \[OverBracket]_\[Yen]_170(\[Yen]|Rgm)_333 */
                v[431] = (v[168] * v[171] * v[419] * v[568]) / 2e0 - v[123] * v[602];
                /* 433 = \[OverBracket]_\[Yen]_172(\[Yen]|Rgm)_333 */
                v[433] = -(v[168] * v[602]) - (v[168] * v[419] * v[568]) / (2e0 * v[97]);
                /* 429 = \[OverBracket]_\[Yen]_168(\[Yen]|Rgm)_333 */
                v[429] = v[429] - v[172] * v[602];
                /* 402 = \[OverBracket]_u\[Phi]_5(\[Yen]|Rgm)_336 */
                v[402] = v[402] - v[322] * v[603] + v[533] * v[617];
                /* 403 = \[OverBracket]_u\[Phi]_5;\[Xi](\[Yen]|Rgm)_336 */
                v[403] = v[527] + v[419] * (v[569] + (v[168] * (v[168] * v[171] * v[326] - v[123] * v[617])) / 2e0 + v[618] -
                    (v[123] * v[168] * v[326]) / (2e0 * v[97])) - (v[122] * v[603]) / v[97];
                /* 423 = \[OverBracket]_\[Yen]_122(\[Yen]|Rgm)_337 */
                v[423] = v[423] + v[578] * v[587] + v[553] * v[601] + v[576] * v[603] - v[147] * v[178] * v[612] + (v[235] * v[405])
                    / v[97] + (v[270] * v[406]) / v[97] + (v[252] * v[415]) / v[97] + (v[245] * v[416]) / v[97] + (v[255] * v[417]) / v[97] +
                    (v[265] * v[418]) / v[97];
                /* 434 = \[OverBracket]_\[Yen]_176(\[Yen]|Rgm)_337 */
                v[434] = v[419] * v[619] + v[612] / v[97] - (v[587] * v[94]) / v[97] - (v[603] * v[95]) / v[97] - (v[601] * v[96])
                    / v[97];
                /* 435 = \[OverBracket]_\[Yen]_178(\[Yen]|Rgm)_337 */
                v[435] = -(v[122] * v[147] * v[612]) + v[414] * v[621] + v[122] * v[147] * v[587] * v[94]
                    + v[122] * v[147] * v[603] * v[95] + v[122] * v[147] * v[601] * v[96] + v[419] * ((-(v[122] * v[142]) - v[132] * v[176]
                        ) * v[235] + (-(v[122] * v[145]) - v[127] * v[176]) * v[245] + (-(v[122] * v[146]) - v[120] * v[176]) * v[252]
                        - v[120] * v[122] * v[255] - v[122] * v[127] * v[265] - v[122] * v[132] * v[270] - (v[123] * v[168] * v[620]) / 2e0
                        + v[122] * v[327] * v[94] * v[95] + v[329] * (v[122] * v[140] * v[94] + v[122] * v[138] * v[95] + v[176] * v[94] * v[95])
                        + v[122] * v[323] * v[94] * v[96] + v[122] * v[318] * v[95] * v[96] + v[325] * (v[122] * v[143] * v[94]
                            + v[122] * v[138] * v[96] + v[176] * v[94] * v[96]) + v[320] * (v[122] * v[143] * v[95] + v[122] * v[140] * v[96]
                                + v[176] * v[95] * v[96]));
                /* 441 = \[OverBracket]_ff_(\[Yen]|Rgm)_337 */
                v[441] = -(v[122] * v[178] * v[235] * v[405]) - v[122] * v[178] * v[270] * v[406] - v[122] * v[178] * v[252] * v[415]
                    - v[122] * v[178] * v[245] * v[416] - v[122] * v[178] * v[255] * v[417] - v[122] * v[178] * v[265] * v[418] + v[441]
                    + v[574] * v[587] + v[552] * v[601] + v[572] * v[603] - v[176] * v[178] * v[612] + v[533] * v[619];
                /* 442 = \[OverBracket]_ff_;\[Xi](\[Yen]|Rgm)_337 */
                v[442] = -(v[122] * v[178] * v[612]) + v[122] * v[178] * v[587] * v[94] + v[122] * v[178] * v[603] * v[95]
                    + v[122] * v[178] * v[601] * v[96] + v[419] * (-0.5e0 * (v[123] * v[168] * v[619]) - 2e0 * v[537] * v[621] -
                        (v[168] * v[339]) / (4e0 * v[97]));
                /* 395 = \[OverBracket]_u\[Phi]_4(\[Yen]|Rgm)_337 */
                v[395] = v[395] - v[322] * v[587] + v[533] * v[622];
                /* 396 = \[OverBracket]_u\[Phi]_4;\[Xi](\[Yen]|Rgm)_337 */
                v[396] = v[526] + v[419] * (v[594] + (v[168] * (-(v[168] * v[171] * v[321]) - v[123] * v[622])) / 2e0 + v[623] +
                    (v[123] * v[168] * v[321]) / (2e0 * v[97])) - (v[122] * v[587]) / v[97];
            };
            v[726] = v[442];
            v[724] = v[441];
            v[749] = v[396];
            v[748] = v[412];
            v[747] = v[403];
            v[746] = v[411];
            v[745] = v[395];
            v[744] = v[402];
            v[743] = v[493];
            v[742] = v[490];
            v[733] = v[492];
            v[732] = v[491];
            v[727] = v[495];
            v[725] = v[494];
            /* 624 = \[OverBracket]_\[OverBracket]_Ft_3|2(Fref|W)\[OverBracket]_3|2(R|Rgm)_3|3(Ft|W)_3|3 */
            v[624] = v[500] * v[62];
            /* 625 = \[OverBracket]_\[OverBracket]_Ft_3|3(Fref|W)\[OverBracket]_3|3(R|Rgm)_3|3(Ft|W)_3|3 */
            v[625] = v[500] * v[65];
            /* 626 = \[OverBracket]_\[OverBracket]_Fref_3|1(Ft|W)\[OverBracket]_3|1(R|Rgm)_3|3(Ft|W)_3|3 */
            v[1383] = 0e0;
            v[1384] = 0e0;
            v[1385] = v[49];
            v[1386] = 0e0;
            v[1387] = 0e0;
            v[1388] = 0e0;
            v[1389] = 0e0;
            v[1390] = 0e0;
            v[1391] = v[50];
            v[1392] = 0e0;
            v[1393] = 0e0;
            v[1394] = 0e0;
            v[1395] = 0e0;
            v[1396] = 0e0;
            v[1397] = v[52];
            v[1398] = 0e0;
            v[1399] = 0e0;
            v[1400] = 0e0;
            v[626] = v[266] * v[500] + v[1382 + i218] * v[90];
            /* 627 = \[OverBracket]_\[OverBracket]_Ft_3|2(Fref|W)\[OverBracket]_3|2(R|Rgm)_3|2(Ft|W)_3|3 */
            v[627] = v[503] * v[61] + v[624];
            /* 628 = \[OverBracket]_\[OverBracket]_Ft_3|3(Fref|W)\[OverBracket]_3|3(R|Rgm)_3|2(Ft|W)_3|3 */
            v[628] = v[625] + v[503] * v[64];
            /* 629 = \[OverBracket]_\[OverBracket]_Fref_3|1(Ft|W)\[OverBracket]_3|1(R|Rgm)_3|2(Ft|W)_3|3 */
            v[629] = v[264] * v[503] + v[626];
            /* 630 = \[OverBracket]_\[OverBracket]_Ft_3|2(Fref|W)\[OverBracket]_3|2(R|Rgm)_3|1(Ft|W)_3|3 */
            v[630] = v[507] * v[60] + v[627];
            /* 631 = \[OverBracket]_\[OverBracket]_Ft_3|3(Fref|W)\[OverBracket]_3|3(R|Rgm)_3|1(Ft|W)_3|3 */
            v[631] = v[628] + v[507] * v[63];
            /* 632 = \[OverBracket]_\[OverBracket]_Fref_3|1(Ft|W)\[OverBracket]_3|1(R|Rgm)_3|1(Ft|W)_3|3 */
            v[632] = v[262] * v[507] + v[629];
            /* 633 = \[OverBracket]_\[OverBracket]_Ft_2|2(Fref|W)\[OverBracket]_2|2(R|Rgm)_2|3(Ft|W)_2|3 */
            v[633] = v[502] * v[62];
            /* 634 = \[OverBracket]_\[OverBracket]_Ft_2|3(Fref|W)\[OverBracket]_2|3(R|Rgm)_2|3(Ft|W)_2|3 */
            v[634] = v[502] * v[65];
            /* 635 = \[OverBracket]_\[OverBracket]_Fref_2|1(Ft|W)\[OverBracket]_2|1(R|Rgm)_2|3(Ft|W)_2|3 */
            v[1365] = 0e0;
            v[1366] = v[49];
            v[1367] = 0e0;
            v[1368] = 0e0;
            v[1369] = 0e0;
            v[1370] = 0e0;
            v[1371] = 0e0;
            v[1372] = v[50];
            v[1373] = 0e0;
            v[1374] = 0e0;
            v[1375] = 0e0;
            v[1376] = 0e0;
            v[1377] = 0e0;
            v[1378] = v[52];
            v[1379] = 0e0;
            v[1380] = 0e0;
            v[1381] = 0e0;
            v[1382] = 0e0;
            v[635] = v[266] * v[502] + v[1364 + i218] * v[90];
            /* 636 = \[OverBracket]_\[OverBracket]_Ft_2|2(Fref|W)\[OverBracket]_2|2(R|Rgm)_2|2(Ft|W)_2|3 */
            v[636] = v[498] * v[61] + v[633];
            /* 637 = \[OverBracket]_\[OverBracket]_Ft_2|3(Fref|W)\[OverBracket]_2|3(R|Rgm)_2|2(Ft|W)_2|3 */
            v[637] = v[634] + v[498] * v[64];
            /* 638 = \[OverBracket]_\[OverBracket]_Fref_2|1(Ft|W)\[OverBracket]_2|1(R|Rgm)_2|2(Ft|W)_2|3 */
            v[638] = v[264] * v[498] + v[635];
            /* 639 = \[OverBracket]_\[OverBracket]_Ft_2|2(Fref|W)\[OverBracket]_2|2(R|Rgm)_2|1(Ft|W)_2|3 */
            v[639] = v[511] * v[60] + v[636];
            /* 640 = \[OverBracket]_\[OverBracket]_Ft_2|3(Fref|W)\[OverBracket]_2|3(R|Rgm)_2|1(Ft|W)_2|3 */
            v[640] = v[511] * v[63] + v[637];
            /* 641 = \[OverBracket]_\[OverBracket]_Fref_2|1(Ft|W)\[OverBracket]_2|1(R|Rgm)_2|1(Ft|W)_2|3 */
            v[641] = v[262] * v[511] + v[638];
            /* 642 = \[OverBracket]_\[OverBracket]_Ft_1|2(Fref|W)\[OverBracket]_1|2(R|Rgm)_1|3(Ft|W)_1|3 */
            v[642] = v[506] * v[62];
            /* 643 = \[OverBracket]_\[OverBracket]_Ft_1|3(Fref|W)\[OverBracket]_1|3(R|Rgm)_1|3(Ft|W)_1|3 */
            v[643] = v[506] * v[65];
            /* 644 = \[OverBracket]_\[OverBracket]_Fref_1|1(Ft|W)\[OverBracket]_1|1(R|Rgm)_1|3(Ft|W)_1|3 */
            v[1347] = v[49];
            v[1348] = 0e0;
            v[1349] = 0e0;
            v[1350] = 0e0;
            v[1351] = 0e0;
            v[1352] = 0e0;
            v[1353] = v[50];
            v[1354] = 0e0;
            v[1355] = 0e0;
            v[1356] = 0e0;
            v[1357] = 0e0;
            v[1358] = 0e0;
            v[1359] = v[52];
            v[1360] = 0e0;
            v[1361] = 0e0;
            v[1362] = 0e0;
            v[1363] = 0e0;
            v[1364] = 0e0;
            v[644] = v[266] * v[506] + v[1346 + i218] * v[90];
            /* 645 = \[OverBracket]_\[OverBracket]_Ft_1|2(Fref|W)\[OverBracket]_1|2(R|Rgm)_1|2(Ft|W)_1|3 */
            v[645] = v[510] * v[61] + v[642];
            /* 646 = \[OverBracket]_\[OverBracket]_Ft_1|3(Fref|W)\[OverBracket]_1|3(R|Rgm)_1|2(Ft|W)_1|3 */
            v[646] = v[510] * v[64] + v[643];
            /* 647 = \[OverBracket]_\[OverBracket]_Fref_1|1(Ft|W)\[OverBracket]_1|1(R|Rgm)_1|2(Ft|W)_1|3 */
            v[647] = v[264] * v[510] + v[644];
            /* 648 = \[OverBracket]_\[OverBracket]_Ft_1|2(Fref|W)\[OverBracket]_1|2(R|Rgm)_1|1(Ft|W)_1|3 */
            v[648] = v[496] * v[60] + v[645];
            /* 649 = \[OverBracket]_\[OverBracket]_Ft_1|3(Fref|W)\[OverBracket]_1|3(R|Rgm)_1|1(Ft|W)_1|3 */
            v[649] = v[496] * v[63] + v[646];
            /* 650 = \[OverBracket]_\[OverBracket]_Fref_1|1(Ft|W)\[OverBracket]_1|1(R|Rgm)_1|1(Ft|W)_1|3 */
            v[650] = v[262] * v[496] + v[647];
            /* 651 = \[OverBracket]_\[OverBracket]_Fref_1|1(Ft|W)\[OverBracket]_1|1(R|Rgm)_1|1;\[Xi](Fref|W)_1|1 */
            v[651] = v[246] * v[497] + v[650];
            /* 652 = \[OverBracket]_\[OverBracket]_Fref_1|1(Ft|W)\[OverBracket]_1|1(R|Rgm)_1|2;\[Xi](Fref|W)_1|1 */
            v[652] = v[244] * v[513] + v[651];
            /* 653 = \[OverBracket]_\[OverBracket]_Fref_1|1(Ft|W)\[OverBracket]_1|1(R|Rgm)_1|3;\[Xi](Fref|W)_1|1 */
            v[653] = v[242] * v[509] + v[652];
            /* 654 = \[OverBracket]_\[OverBracket]_Fref_2|1(Ft|W)\[OverBracket]_2|1(R|Rgm)_2|1;\[Xi](Fref|W)_2|1 */
            v[654] = v[246] * v[512] + v[641];
            /* 655 = \[OverBracket]_\[OverBracket]_Fref_2|1(Ft|W)\[OverBracket]_2|1(R|Rgm)_2|2;\[Xi](Fref|W)_2|1 */
            v[655] = v[244] * v[499] + v[654];
            /* 656 = \[OverBracket]_\[OverBracket]_Fref_2|1(Ft|W)\[OverBracket]_2|1(R|Rgm)_2|3;\[Xi](Fref|W)_2|1 */
            v[656] = v[242] * v[505] + v[655];
            /* 657 = \[OverBracket]_\[OverBracket]_Fref_3|1(Ft|W)\[OverBracket]_3|1(R|Rgm)_3|1;\[Xi](Fref|W)_3|1 */
            v[657] = v[246] * v[508] + v[632];
            /* 658 = \[OverBracket]_\[OverBracket]_Fref_3|1(Ft|W)\[OverBracket]_3|1(R|Rgm)_3|2;\[Xi](Fref|W)_3|1 */
            v[658] = v[244] * v[504] + v[657];
            /* 659 = \[OverBracket]_\[OverBracket]_Fref_3|1(Ft|W)\[OverBracket]_3|1(R|Rgm)_3|3;\[Xi](Fref|W)_3|1 */
            v[659] = v[242] * v[501] + v[658];
            /* 660 = \[OverBracket]_\[Yen]_232\[OverBracket]_(Ft|Rgm)_3|3(Fref|W)_3|3 */
            v[660] = v[207] * v[216] * v[630] + v[206] * v[216] * v[631] + v[204] * v[216] * v[639] + v[203] * v[216] * v[640]
                + v[192] * v[216] * v[648] + v[191] * v[216] * v[649];
            /* 661 = \[OverBracket]_\[OverBracket]_Ft_1|1(W|W)\[OverBracket]_(Fref|Rgm)_1|1(Ft|W)_1|1 */
            v[661] = (v[240] * v[3] * v[648]) / 2e0 + (v[238] * v[4] * v[649]) / 2e0 + v[229] * v[653];
            /* 662 = \[OverBracket]_\[OverBracket]_Ft_2|1(W|W)\[OverBracket]_(Fref|Rgm)_2|1(Ft|W)_2|1 */
            v[662] = (v[240] * v[3] * v[639]) / 2e0 + (v[238] * v[4] * v[640]) / 2e0 + v[229] * v[656];
            /* 663 = \[OverBracket]_\[OverBracket]_Ft_3|1(W|W)\[OverBracket]_(Fref|Rgm)_3|1(Ft|W)_3|1 */
            v[663] = (v[240] * v[3] * v[630]) / 2e0 + (v[238] * v[4] * v[631]) / 2e0 + v[229] * v[659];
            /* 664 = \[OverBracket]_\[Yen]_220\[OverBracket]_(Ft|Rgm)_3|1(W|W) */
            v[664] = v[205] * v[216] * v[630] + v[202] * v[216] * v[639] + v[190] * v[216] * v[648] + v[191] * v[216] * v[661]
                + v[203] * v[216] * v[662] + v[206] * v[216] * v[663];
            /* 665 = \[OverBracket]_\[Yen]_221\[OverBracket]_(Ft|Rgm)_3|1(W|W) */
            v[665] = v[205] * v[216] * v[631] + v[202] * v[216] * v[640] + v[190] * v[216] * v[649] + v[192] * v[216] * v[661]
                + v[204] * v[216] * v[662] + v[207] * v[216] * v[663];
            /* 666 = \[OverBracket]_\[Yen]_222\[OverBracket]_(Ft|Rgm)_3|1(W|W) */
            v[666] = (v[190] * v[216] * v[223] * v[224] * v[661]) / 2e0 + (v[202] * v[216] * v[223] * v[224] * v[662]) / 2e0 +
                (v[205] * v[216] * v[223] * v[224] * v[663]) / 2e0;
            /* 667 = \[OverBracket]_Ft_3|1(\[Yen]|Rgm)_222 */
            v[667] = v[216] * v[220] * v[630] + v[216] * v[221] * v[631] + (v[216] * v[222] * v[223] * v[224] * v[663]) / 2e0
                + v[206] * v[664] + v[207] * v[665] + 2e0 * v[205] * v[666];
            /* 668 = \[OverBracket]_Ft_2|1(\[Yen]|Rgm)_222 */
            v[668] = v[216] * v[220] * v[639] + v[216] * v[221] * v[640] + (v[216] * v[222] * v[223] * v[224] * v[662]) / 2e0
                + v[203] * v[664] + v[204] * v[665] + 2e0 * v[202] * v[666];
            /* 669 = \[OverBracket]_Ft_1|1(\[Yen]|Rgm)_222 */
            v[669] = v[216] * v[220] * v[648] + v[216] * v[221] * v[649] + (v[216] * v[222] * v[223] * v[224] * v[661]) / 2e0
                + v[191] * v[664] + v[192] * v[665] + 2e0 * v[190] * v[666];
            /* 670 = \[OverBracket]_Fref_3|1(Ft|Rgm)_3|1 */
            v[670] = v[229] * v[667];
            /* 671 = \[OverBracket]_Fref_2|1(Ft|Rgm)_2|1 */
            v[671] = v[229] * v[668];
            /* 672 = \[OverBracket]_Fref_1|1(Ft|Rgm)_1|1 */
            v[672] = v[229] * v[669];
            /* 673 = \[OverBracket]_R_3|3;\[Xi](Fref|Rgm)_3|1 */
            v[673] = v[242] * v[670];
            /* 674 = \[OverBracket]_R_3|2;\[Xi](Fref|Rgm)_3|1 */
            v[674] = v[244] * v[670];
            /* 675 = \[OverBracket]_R_3|1;\[Xi](Fref|Rgm)_3|1 */
            v[675] = v[246] * v[670];
            /* 676 = \[OverBracket]_R_2|3;\[Xi](Fref|Rgm)_2|1 */
            v[676] = v[242] * v[671];
            /* 677 = \[OverBracket]_R_2|2;\[Xi](Fref|Rgm)_2|1 */
            v[677] = v[244] * v[671];
            /* 678 = \[OverBracket]_R_2|1;\[Xi](Fref|Rgm)_2|1 */
            v[678] = v[246] * v[671];
            /* 679 = \[OverBracket]_R_1|3;\[Xi](Fref|Rgm)_1|1 */
            v[679] = v[242] * v[672];
            /* 680 = \[OverBracket]_R_1|2;\[Xi](Fref|Rgm)_1|1 */
            v[680] = v[244] * v[672];
            /* 681 = \[OverBracket]_R_1|1;\[Xi](Fref|Rgm)_1|1 */
            v[681] = v[246] * v[672];
            /* 682 = \[OverBracket]_Ft_1|2(\[Yen]|Rgm)_232 */
            v[682] = v[216] * v[232] * v[649] + v[192] * v[660] + v[216] * v[220] * v[661] + v[190] * v[664] + (v[240] * v[3] * v[669])
                / 2e0;
            /* 683 = \[OverBracket]_Ft_1|3(\[Yen]|Rgm)_232 */
            v[683] = v[216] * v[232] * v[648] + v[191] * v[660] + v[216] * v[221] * v[661] + v[190] * v[665] + (v[238] * v[4] * v[669])
                / 2e0;
            /* 684 = \[OverBracket]_Ft_2|2(\[Yen]|Rgm)_232 */
            v[684] = v[216] * v[232] * v[640] + v[204] * v[660] + v[216] * v[220] * v[662] + v[202] * v[664] + (v[240] * v[3] * v[668])
                / 2e0;
            /* 685 = \[OverBracket]_Ft_2|3(\[Yen]|Rgm)_232 */
            v[685] = v[216] * v[232] * v[639] + v[203] * v[660] + v[216] * v[221] * v[662] + v[202] * v[665] + (v[238] * v[4] * v[668])
                / 2e0;
            /* 686 = \[OverBracket]_Ft_3|2(\[Yen]|Rgm)_232 */
            v[686] = v[216] * v[232] * v[631] + v[207] * v[660] + v[216] * v[220] * v[663] + v[205] * v[664] + (v[240] * v[3] * v[667])
                / 2e0;
            /* 687 = \[OverBracket]_Ft_3|3(\[Yen]|Rgm)_232 */
            v[687] = v[216] * v[232] * v[630] + v[206] * v[660] + v[216] * v[221] * v[663] + v[205] * v[665] + (v[238] * v[4] * v[667])
                / 2e0;
            /* 688 = \[OverBracket]_R_1|1(Ft|Rgm)_1|3 */
            v[688] = v[262] * v[672] + v[60] * v[682] + v[63] * v[683];
            /* 689 = \[OverBracket]_R_1|2(Ft|Rgm)_1|3 */
            v[689] = v[264] * v[672] + v[61] * v[682] + v[64] * v[683];
            /* 690 = \[OverBracket]_R_1|3(Ft|Rgm)_1|3 */
            v[690] = v[266] * v[672] + v[62] * v[682] + v[65] * v[683];
            /* 691 = \[OverBracket]_R_2|1(Ft|Rgm)_2|3 */
            v[691] = v[262] * v[671] + v[60] * v[684] + v[63] * v[685];
            /* 692 = \[OverBracket]_R_2|2(Ft|Rgm)_2|3 */
            v[692] = v[264] * v[671] + v[61] * v[684] + v[64] * v[685];
            /* 693 = \[OverBracket]_R_2|3(Ft|Rgm)_2|3 */
            v[693] = v[266] * v[671] + v[62] * v[684] + v[65] * v[685];
            /* 694 = \[OverBracket]_R_3|1(Ft|Rgm)_3|3 */
            v[694] = v[262] * v[670] + v[60] * v[686] + v[63] * v[687];
            /* 695 = \[OverBracket]_R_3|2(Ft|Rgm)_3|3 */
            v[695] = v[264] * v[670] + v[61] * v[686] + v[64] * v[687];
            /* 696 = \[OverBracket]_R_3|3(Ft|Rgm)_3|3 */
            v[696] = v[266] * v[670] + v[62] * v[686] + v[65] * v[687];
            b697 = b98;
            v[697] = b697;
            if (b697) {
                v[720] = v[680] / 720e0;
                v[719] = v[678] / 720e0;
                v[716] = v[691] / 720e0;
                v[715] = v[689] / 720e0;
                v[712] = v[679] / 720e0;
                v[711] = v[675] / 720e0;
                v[708] = v[694] / 720e0;
                v[707] = v[690] / 720e0;
                v[704] = v[676] / 720e0;
                v[703] = v[674] / 720e0;
                v[700] = v[695] / 720e0;
                v[699] = v[693] / 720e0;
                /* 494 = \[OverBracket]_\[Yen]_132(v|Rgm)_3|3;\[Xi] */
                v[494] = v[494] - (v[155] * v[673]) / 720e0 - (v[103] * v[696]) / 720e0;
                /* 495 = \[OverBracket]_\[Yen]_142(v|Rgm)_3|3;\[Xi] */
                v[495] = v[495] - (v[103] * v[673]) / 720e0;
                /* 698 = \[OverBracket]_\[Yen]_116(v|Rgm)_2|3 */
                v[698] = v[699] + v[700];
                /* 701 = \[OverBracket]_\[Yen]_117(v|Rgm)_2|3 */
                v[701] = v[699] - v[700];
                /* 702 = \[OverBracket]_\[Yen]_156(v|Rgm)_2|3;\[Xi] */
                v[702] = v[703] + v[704];
                /* 705 = \[OverBracket]_\[Yen]_150(v|Rgm)_2|3;\[Xi] */
                v[705] = -v[703] + v[704];
                /* 491 = \[OverBracket]_\[Yen]_127(v|Rgm)_2|2;\[Xi] */
                v[491] = v[491] - (v[155] * v[677]) / 720e0 - (v[103] * v[692]) / 720e0;
                /* 492 = \[OverBracket]_\[Yen]_145(v|Rgm)_2|2;\[Xi] */
                v[492] = v[492] - (v[103] * v[677]) / 720e0;
                /* 706 = \[OverBracket]_\[Yen]_113(v|Rgm)_1|3 */
                v[706] = v[707] + v[708];
                /* 709 = \[OverBracket]_\[Yen]_114(v|Rgm)_1|3 */
                v[709] = v[707] - v[708];
                /* 710 = \[OverBracket]_\[Yen]_157(v|Rgm)_1|3;\[Xi] */
                v[710] = v[711] + v[712];
                /* 713 = \[OverBracket]_\[Yen]_151(v|Rgm)_1|3;\[Xi] */
                v[713] = -v[711] + v[712];
                /* 714 = \[OverBracket]_\[Yen]_107(v|Rgm)_1|2 */
                v[714] = v[715] + v[716];
                /* 717 = \[OverBracket]_\[Yen]_108(v|Rgm)_1|2 */
                v[717] = v[715] - v[716];
                /* 718 = \[OverBracket]_\[Yen]_158(v|Rgm)_1|2;\[Xi] */
                v[718] = v[719] + v[720];
                /* 721 = \[OverBracket]_\[Yen]_152(v|Rgm)_1|2;\[Xi] */
                v[721] = -v[719] + v[720];
                /* 490 = \[OverBracket]_\[Yen]_120(v|Rgm)_1|1;\[Xi] */
                v[490] = v[490] - (v[155] * v[681]) / 720e0 - (v[103] * v[688]) / 720e0;
                /* 493 = \[OverBracket]_\[Yen]_146(v|Rgm)_1|1;\[Xi] */
                v[493] = v[493] - (v[103] * v[681]) / 720e0;
                /* 420 = \[OverBracket]_\[Yen]_103(\[Yen]|Rgm)_158 */
                v[420] = v[420] - (v[142] * v[673]) / 720e0 - (v[145] * v[677]) / 720e0 - (v[146] * v[681]) / 720e0 - (v[120] * v[688])
                    / 720e0 - (v[127] * v[692]) / 720e0 - (v[132] * v[696]) / 720e0 + v[487] * v[702] + v[484] * v[710] + v[481] * v[718]
                    + v[714] * v[94] * v[95] + v[706] * v[94] * v[96] + v[698] * v[95] * v[96];
                /* 428 = \[OverBracket]_\[Yen]_155(\[Yen]|Rgm)_158 */
                v[428] = v[428] - (v[132] * v[673]) / 720e0 - (v[127] * v[677]) / 720e0 - (v[120] * v[681]) / 720e0
                    + v[718] * v[94] * v[95] + v[710] * v[94] * v[96] + v[702] * v[95] * v[96];
                /* 427 = \[OverBracket]_\[Yen]_153(\[Yen]|Rgm)_103 */
                v[427] = v[427] + v[420] * v[97];
                /* 441 = \[OverBracket]_ff_(\[Yen]|Rgm)_103 */
                v[441] = v[153] * v[420] + v[441];
                /* 395 = \[OverBracket]_u\[Phi]_4(\[Yen]|Rgm)_150 */
                v[395] = v[395] - 6e0 * v[105] * v[701] - 6e0 * v[149] * v[705] + v[306] * v[710] + v[311] * v[718]
                    + v[103] * v[714] * v[95] + v[103] * v[706] * v[96];
                /* 396 = \[OverBracket]_u\[Phi]_4;\[Xi](\[Yen]|Rgm)_150 */
                v[396] = v[396] - 6e0 * v[105] * v[705] + v[103] * v[718] * v[95] + v[103] * v[710] * v[96];
                /* 402 = \[OverBracket]_u\[Phi]_5(\[Yen]|Rgm)_151 */
                v[402] = v[402] + v[306] * v[702] + 6e0 * v[105] * v[709] + 6e0 * v[149] * v[713] + v[312] * v[718]
                    + v[103] * v[714] * v[94] + v[103] * v[698] * v[96];
                /* 403 = \[OverBracket]_u\[Phi]_5;\[Xi](\[Yen]|Rgm)_151 */
                v[403] = v[403] + 6e0 * v[105] * v[713] + v[103] * v[718] * v[94] + v[103] * v[702] * v[96];
                /* 421 = \[OverBracket]_\[Yen]_105(\[Yen]|Rgm)_152 */
                v[421] = v[421] - 6e0 * v[138] * v[705] + 6e0 * v[140] * v[713] - 6e0 * v[143] * v[721] - 6e0 * v[701] * v[94]
                    + 6e0 * v[709] * v[95] - 6e0 * v[717] * v[96];
                /* 426 = \[OverBracket]_\[Yen]_149(\[Yen]|Rgm)_152 */
                v[426] = v[426] - 6e0 * v[705] * v[94] + 6e0 * v[713] * v[95] - 6e0 * v[721] * v[96];
                /* 411 = \[OverBracket]_u\[Phi]_6(\[Yen]|Rgm)_152 */
                v[411] = v[411] + v[311] * v[702] + v[312] * v[710] - 6e0 * v[105] * v[717] - 6e0 * v[149] * v[721]
                    + v[103] * v[706] * v[94] + v[103] * v[698] * v[95];
                /* 412 = \[OverBracket]_u\[Phi]_6;\[Xi](\[Yen]|Rgm)_152 */
                v[412] = v[412] - 6e0 * v[105] * v[721] + v[103] * v[710] * v[94] + v[103] * v[702] * v[95];
                /* 425 = \[OverBracket]_\[Yen]_148(\[Yen]|Rgm)_105 */
                v[425] = v[425] + v[421] * v[97];
                /* 441 = \[OverBracket]_ff_(\[Yen]|Rgm)_105 */
                v[441] = v[148] * v[421] + v[441];
                /* 425 = \[OverBracket]_\[Yen]_148(\[Yen]|Rgm)_149 */
                v[425] = v[425] + v[147] * v[426];
                /* 722 = \[OverBracket]_\[Yen]_154(\[Yen]|Rgm)_149 */
                v[722] = v[426];
                /* 442 = \[OverBracket]_ff_;\[Xi](\[Yen]|Rgm)_149 */
                v[442] = v[148] * v[426] + v[442];
                /* 441 = \[OverBracket]_ff_(\[Yen]|Rgm)_148 */
                v[441] = v[425] + v[441];
                /* 427 = \[OverBracket]_\[Yen]_153(\[Yen]|Rgm)_155 */
                v[427] = v[427] + v[147] * v[428];
                /* 723 = \[OverBracket]_\[Yen]_154(\[Yen]|Rgm)_155 */
                v[723] = v[428] + v[722];
                /* 442 = \[OverBracket]_ff_;\[Xi](\[Yen]|Rgm)_155 */
                v[442] = v[153] * v[428] + v[442];
                /* 441 = \[OverBracket]_ff_(\[Yen]|Rgm)_153 */
                v[441] = v[427] + v[441];
                /* 441 = \[OverBracket]_ff_(\[Yen]|Rgm)_154 */
                v[441] = v[441] + v[147] * v[723];
                /* 442 = \[OverBracket]_ff_;\[Xi](\[Yen]|Rgm)_154 */
                v[442] = v[442] + v[723] * v[97];
            }
            else {
                /* 494 = \[OverBracket]_\[Yen]_132(v|Rgm)_3|3;\[Xi] */
                v[494] = v[322] * v[673] + v[725] + (v[122] * v[696]) / v[97];
                /* 495 = \[OverBracket]_\[Yen]_142(v|Rgm)_3|3;\[Xi] */
                v[495] = v[727] + (v[122] * v[673]) / v[97];
                /* 728 = \[OverBracket]_\[Yen]_131(v|Rgm)_2|3 */
                v[728] = v[693] + v[695];
                /* 729 = \[OverBracket]_\[Yen]_130(v|Rgm)_2|3 */
                v[729] = v[693] - v[695];
                /* 730 = \[OverBracket]_\[Yen]_177(v|Rgm)_2|3;\[Xi] */
                v[730] = v[674] + v[676];
                /* 731 = \[OverBracket]_\[Yen]_173(v|Rgm)_2|3;\[Xi] */
                v[731] = -v[674] + v[676];
                /* 491 = \[OverBracket]_\[Yen]_127(v|Rgm)_2|2;\[Xi] */
                v[491] = v[322] * v[677] + v[732] + (v[122] * v[692]) / v[97];
                /* 492 = \[OverBracket]_\[Yen]_145(v|Rgm)_2|2;\[Xi] */
                v[492] = v[733] + (v[122] * v[677]) / v[97];
                /* 734 = \[OverBracket]_\[Yen]_129(v|Rgm)_1|3 */
                v[734] = v[690] + v[694];
                /* 735 = \[OverBracket]_\[Yen]_128(v|Rgm)_1|3 */
                v[735] = v[690] - v[694];
                /* 736 = \[OverBracket]_\[Yen]_179(v|Rgm)_1|3;\[Xi] */
                v[736] = v[675] + v[679];
                /* 737 = \[OverBracket]_\[Yen]_174(v|Rgm)_1|3;\[Xi] */
                v[737] = -v[675] + v[679];
                /* 738 = \[OverBracket]_\[Yen]_126(v|Rgm)_1|2 */
                v[738] = v[689] + v[691];
                /* 739 = \[OverBracket]_\[Yen]_125(v|Rgm)_1|2 */
                v[739] = v[689] - v[691];
                /* 740 = \[OverBracket]_\[Yen]_180(v|Rgm)_1|2;\[Xi] */
                v[740] = v[678] + v[680];
                /* 741 = \[OverBracket]_\[Yen]_175(v|Rgm)_1|2;\[Xi] */
                v[741] = -v[678] + v[680];
                /* 490 = \[OverBracket]_\[Yen]_120(v|Rgm)_1|1;\[Xi] */
                v[490] = v[322] * v[681] + v[742] + (v[122] * v[688]) / v[97];
                /* 493 = \[OverBracket]_\[Yen]_146(v|Rgm)_1|1;\[Xi] */
                v[493] = v[743] + (v[122] * v[681]) / v[97];
                /* 423 = \[OverBracket]_\[Yen]_122(\[Yen]|Rgm)_180 */
                v[423] = v[423] + v[613] * v[673] + v[611] * v[677] + v[610] * v[681] + v[539] * v[730] + v[540] * v[736]
                    + v[564] * v[740] + (v[120] * v[688]) / v[97] + (v[127] * v[692]) / v[97] + (v[132] * v[696]) / v[97] -
                    (v[738] * v[94] * v[95]) / v[97] - (v[734] * v[94] * v[96]) / v[97] - (v[728] * v[95] * v[96]) / v[97];
                /* 434 = \[OverBracket]_\[Yen]_176(\[Yen]|Rgm)_180 */
                v[434] = v[434] + (v[132] * v[673]) / v[97] + (v[127] * v[677]) / v[97] + (v[120] * v[681]) / v[97] -
                    (v[740] * v[94] * v[95]) / v[97] - (v[736] * v[94] * v[96]) / v[97] - (v[730] * v[95] * v[96]) / v[97];
                /* 435 = \[OverBracket]_\[Yen]_178(\[Yen]|Rgm)_180 */
                v[435] = v[435] - v[122] * v[132] * v[147] * v[673] - v[122] * v[127] * v[147] * v[677]
                    - v[120] * v[122] * v[147] * v[681] + v[122] * v[147] * v[740] * v[94] * v[95] + v[122] * v[147] * v[736] * v[94] * v[96]
                    + v[122] * v[147] * v[730] * v[95] * v[96];
                /* 432 = \[OverBracket]_\[Yen]_171(\[Yen]|Rgm)_122 */
                v[432] = v[423] + v[432];
                /* 395 = \[OverBracket]_u\[Phi]_4(\[Yen]|Rgm)_173 */
                v[395] = -(v[123] * v[168] * v[729]) + v[333] * v[731] + v[332] * v[736] + v[336] * v[740] + v[745] -
                    (v[122] * v[738] * v[95]) / v[97] - (v[122] * v[734] * v[96]) / v[97];
                /* 396 = \[OverBracket]_u\[Phi]_4;\[Xi](\[Yen]|Rgm)_173 */
                v[396] = -(v[123] * v[168] * v[731]) + v[749] - (v[122] * v[740] * v[95]) / v[97] - (v[122] * v[736] * v[96]) / v[97];
                /* 402 = \[OverBracket]_u\[Phi]_5(\[Yen]|Rgm)_174 */
                v[402] = v[332] * v[730] + v[123] * v[168] * v[735] - v[333] * v[737] + v[337] * v[740] + v[744] -
                    (v[122] * v[738] * v[94]) / v[97] - (v[122] * v[728] * v[96]) / v[97];
                /* 403 = \[OverBracket]_u\[Phi]_5;\[Xi](\[Yen]|Rgm)_174 */
                v[403] = v[123] * v[168] * v[737] + v[747] - (v[122] * v[740] * v[94]) / v[97] - (v[122] * v[730] * v[96]) / v[97];
                /* 431 = \[OverBracket]_\[Yen]_170(\[Yen]|Rgm)_175 */
                v[431] = v[431] - v[123] * v[731] * v[94] + v[123] * v[737] * v[95] - v[123] * v[741] * v[96];
                /* 433 = \[OverBracket]_\[Yen]_172(\[Yen]|Rgm)_175 */
                v[433] = v[433] - v[168] * v[731] * v[94] + v[168] * v[737] * v[95] - v[168] * v[741] * v[96];
                /* 411 = \[OverBracket]_u\[Phi]_6(\[Yen]|Rgm)_175 */
                v[411] = v[336] * v[730] + v[337] * v[736] - v[123] * v[168] * v[739] + v[333] * v[741] + v[746] -
                    (v[122] * v[734] * v[94]) / v[97] - (v[122] * v[728] * v[95]) / v[97];
                /* 412 = \[OverBracket]_u\[Phi]_6;\[Xi](\[Yen]|Rgm)_175 */
                v[412] = -(v[123] * v[168] * v[741]) + v[748] - (v[122] * v[736] * v[94]) / v[97] - (v[122] * v[730] * v[95]) / v[97];
                /* 424 = \[OverBracket]_\[Yen]_123(\[Yen]|Rgm)_176 */
                v[424] = v[424] - v[169] * v[434] + v[731] * v[750] + v[737] * v[751] + v[741] * v[752] - v[168] * v[729] * v[94]
                    + v[168] * v[735] * v[95] - v[168] * v[739] * v[96];
                /* 430 = \[OverBracket]_\[Yen]_169(\[Yen]|Rgm)_176 */
                v[430] = v[430] - v[123] * v[434];
                /* 422 = \[OverBracket]_\[Yen]_121(\[Yen]|Rgm)_123 */
                v[422] = v[422] + v[171] * v[424];
                /* 436 = \[OverBracket]_\[Yen]_338(\[Yen]|Rgm)_170 */
                v[436] = -(v[169] * v[431]) + v[436];
                /* 430 = \[OverBracket]_\[Yen]_169(\[Yen]|Rgm)_170 */
                v[430] = v[430] - v[338] * v[431];
                /* 432 = \[OverBracket]_\[Yen]_171(\[Yen]|Rgm)_172 */
                v[432] = v[432] + v[169] * v[433];
                /* 430 = \[OverBracket]_\[Yen]_169(\[Yen]|Rgm)_172 */
                v[430] = v[430] + v[171] * v[433];
                /* 422 = \[OverBracket]_\[Yen]_121(\[Yen]|Rgm)_171 */
                v[422] = v[422] - v[123] * v[432];
                /* 422 = \[OverBracket]_\[Yen]_121(\[Yen]|Rgm)_338 */
                v[422] = v[422] - 2e0 * v[436] * v[531];
                /* 441 = \[OverBracket]_ff_(\[Yen]|Rgm)_121 */
                v[441] = (v[168] * v[422]) / 2e0 + (-(v[122] * v[142] * v[178]) - v[132] * v[176] * v[178]) * v[673] + (-
                    (v[122] * v[145] * v[178]) - v[127] * v[176] * v[178]) * v[677] + (-(v[122] * v[146] * v[178])
                        - v[120] * v[176] * v[178]) * v[681] - v[120] * v[122] * v[178] * v[688] - v[122] * v[127] * v[178] * v[692]
                    - v[122] * v[132] * v[178] * v[696] + v[724] + v[736] * (v[755] + v[756] + v[757]) + v[740] * (v[758] + v[759] + v[760])
                    + v[122] * v[178] * v[738] * v[94] * v[95] + v[122] * v[178] * v[734] * v[94] * v[96]
                    + v[122] * v[178] * v[728] * v[95] * v[96] + v[730] * (v[753] + v[754] + v[122] * v[140] * v[178] * v[96]);
                /* 429 = \[OverBracket]_\[Yen]_168(\[Yen]|Rgm)_169 */
                v[429] = v[429] + (v[147] * v[430]) / 2e0 + v[731] * v[761] + v[737] * v[762] + v[741] * v[763] - v[123] * v[729] * v[94]
                    + v[123] * v[735] * v[95] - v[123] * v[739] * v[96];
                /* 442 = \[OverBracket]_ff_;\[Xi](\[Yen]|Rgm)_169 */
                v[442] = (v[168] * v[430]) / 2e0 - v[122] * v[132] * v[178] * v[673] - v[122] * v[127] * v[178] * v[677]
                    - v[120] * v[122] * v[178] * v[681] + v[726] + v[122] * v[178] * v[740] * v[94] * v[95]
                    + v[122] * v[178] * v[736] * v[94] * v[96] + v[122] * v[178] * v[730] * v[95] * v[96];
                /* 441 = \[OverBracket]_ff_(\[Yen]|Rgm)_168 */
                v[441] = v[441] - (v[168] * v[429]) / (2e0 * v[97]);
                /* 441 = \[OverBracket]_ff_(\[Yen]|Rgm)_178 */
                v[441] = v[441] - 2e0 * v[435] * v[537];
            };
            /* 764 = \[OverBracket]_\[Yen]_101(ff|Rgm) */
            v[764] = v[441];
            /* 765 = \[OverBracket]_\[Yen]_100(ff|Rgm) */
            v[765] = v[441];
            /* 766 = \[OverBracket]_\[Yen]_110(ff|Rgm) */
            v[766] = v[441];
            /* 767 = \[OverBracket]_\[Yen]_101(\[Yen]|Rgm)_120 */
            v[767] = v[490] + v[764];
            /* 768 = \[OverBracket]_\[Yen]_100(\[Yen]|Rgm)_120 */
            v[768] = v[767];
            /* 769 = \[OverBracket]_\[Yen]_101(\[Yen]|Rgm)_127 */
            v[769] = v[491] + v[767];
            /* 770 = \[OverBracket]_\[Yen]_110(\[Yen]|Rgm)_127 */
            v[770] = v[491] + v[766];
            /* 771 = \[OverBracket]_\[Yen]_144(\[Yen]|Rgm)_145 */
            v[771] = v[492];
            /* 772 = \[OverBracket]_\[Yen]_139(\[Yen]|Rgm)_145 */
            v[772] = v[492];
            /* 773 = \[OverBracket]_\[Yen]_144(\[Yen]|Rgm)_146 */
            v[773] = v[493] + v[771];
            /* 774 = \[OverBracket]_\[Yen]_141(\[Yen]|Rgm)_146 */
            v[774] = v[493];
            /* 775 = \[OverBracket]_\[Yen]_144(ff|Rgm)_;\[Xi] */
            v[775] = v[442] + v[773];
            /* 776 = \[OverBracket]_\[Yen]_141(ff|Rgm)_;\[Xi] */
            v[776] = v[442] + v[774];
            /* 777 = \[OverBracket]_\[Yen]_139(ff|Rgm)_;\[Xi] */
            v[777] = v[442] + v[772];
            /* 411 = \[OverBracket]_u\[Phi]_6(\[Yen]|Rgm)_144 */
            v[411] = v[411] + 2e0 * v[143] * v[775] + 2e0 * v[769] * v[96];
            /* 412 = \[OverBracket]_u\[Phi]_6;\[Xi](\[Yen]|Rgm)_144 */
            v[412] = v[412] + 2e0 * v[775] * v[96];
            /* 778 = \[OverBracket]_peIO_3|6(u\[Phi]|Rgm)_6 */
            v[778] = v[411] * v[45];
            /* 779 = \[OverBracket]_peIO_2|6(u\[Phi]|Rgm)_6 */
            v[779] = v[411] * v[43];
            /* 780 = \[OverBracket]_peIO_1|6(u\[Phi]|Rgm)_6 */
            v[780] = v[41] * v[411];
            /* 781 = \[OverBracket]_\[Yen]_100(\[Yen]|Rgm)_132 */
            v[781] = v[494] + v[768];
            /* 782 = \[OverBracket]_\[Yen]_110(\[Yen]|Rgm)_132 */
            v[782] = v[494] + v[770];
            /* 783 = \[OverBracket]_\[Yen]_141(\[Yen]|Rgm)_142 */
            v[783] = v[495] + v[776];
            /* 784 = \[OverBracket]_\[Yen]_139(\[Yen]|Rgm)_142 */
            v[784] = v[495] + v[777];
            /* 402 = \[OverBracket]_u\[Phi]_5(\[Yen]|Rgm)_141 */
            v[402] = v[402] + 2e0 * v[140] * v[783] + 2e0 * v[781] * v[95];
            /* 403 = \[OverBracket]_u\[Phi]_5;\[Xi](\[Yen]|Rgm)_141 */
            v[403] = v[403] + 2e0 * v[783] * v[95];
            /* 785 = \[OverBracket]_peIO_3|5(u\[Phi]|Rgm)_5 */
            v[785] = v[402] * v[45];
            /* 786 = \[OverBracket]_peIO_2|5(u\[Phi]|Rgm)_5 */
            v[786] = v[402] * v[43];
            /* 787 = \[OverBracket]_peIO_1|5(u\[Phi]|Rgm)_5 */
            v[787] = v[402] * v[41];
            /* 395 = \[OverBracket]_u\[Phi]_4(\[Yen]|Rgm)_139 */
            v[395] = v[395] + 2e0 * v[138] * v[784] + 2e0 * v[782] * v[94];
            /* 396 = \[OverBracket]_u\[Phi]_4;\[Xi](\[Yen]|Rgm)_139 */
            v[396] = v[396] + 2e0 * v[784] * v[94];
            /* 788 = \[OverBracket]_peIO_3|4(u\[Phi]|Rgm)_4 */
            v[788] = v[395] * v[45];
            /* 789 = \[OverBracket]_peIO_2|4(u\[Phi]|Rgm)_4 */
            v[789] = v[395] * v[43];
            /* 790 = \[OverBracket]_peIO_1|4(u\[Phi]|Rgm)_4 */
            v[790] = v[395] * v[41];
            /* 791 = \[OverBracket]_peIO_3|4(u\[Phi]|Rgm)_4;\[Xi] */
            v[791] = v[396] * v[52] + v[788];
            /* 792 = \[OverBracket]_peIO_2|4(u\[Phi]|Rgm)_4;\[Xi] */
            v[792] = v[396] * v[50] + v[789];
            /* 793 = \[OverBracket]_peIO_1|4(u\[Phi]|Rgm)_4;\[Xi] */
            v[793] = v[396] * v[49] + v[790];
            /* 794 = \[OverBracket]_peIO_3|5(u\[Phi]|Rgm)_5;\[Xi] */
            v[794] = v[403] * v[52] + v[785];
            /* 795 = \[OverBracket]_peIO_2|5(u\[Phi]|Rgm)_5;\[Xi] */
            v[795] = v[403] * v[50] + v[786];
            /* 796 = \[OverBracket]_peIO_1|5(u\[Phi]|Rgm)_5;\[Xi] */
            v[796] = v[403] * v[49] + v[787];
            /* 797 = \[OverBracket]_peIO_3|6(u\[Phi]|Rgm)_6;\[Xi] */
            v[797] = v[412] * v[52] + v[778];
            /* 798 = \[OverBracket]_peIO_2|6(u\[Phi]|Rgm)_6;\[Xi] */
            v[798] = v[412] * v[50] + v[779];
            /* 799 = \[OverBracket]_peIO_1|6(u\[Phi]|Rgm)_6;\[Xi] */
            v[799] = v[412] * v[49] + v[780];
            v[1401] = v[49] * v[672];
            v[1402] = v[49] * v[671];
            v[1403] = v[49] * v[670];
            v[1404] = v[793];
            v[1405] = v[796];
            v[1406] = v[799];
            v[1407] = v[50] * v[672];
            v[1408] = v[50] * v[671];
            v[1409] = v[50] * v[670];
            v[1410] = v[792];
            v[1411] = v[795];
            v[1412] = v[798];
            v[1413] = v[52] * v[672];
            v[1414] = v[52] * v[671];
            v[1415] = v[52] * v[670];
            v[1416] = v[791];
            v[1417] = v[794];
            v[1418] = v[797];
            v[800] = 0e0;/*debug*/
            R[i218 - 1] += v[377] * v[40];
            v[378] = 0e0;/*debug*/
            for (i379 = i218; i379 <= 18; i379++) {
                v[379] = i379;
                /* 381 = \[DoubleStruckCapitalG]_n */
                v[381] = v[1141 + i379];
                /* 801 = Kgmn */
                v[801] = v[1400 + i379];
                T[i218 - 1][i379 - 1] += v[40] * v[801];
                v[802] = 0e0;/*debug*/
            };/* end for */
        };/* end for */
    };/* end for */

    for (int i = 0; i < 18; i++)
    {
        for (int j = 0; j < 18; j++)
            if (j != i)
                T[j][i] = T[i][j];
    }
};

#endif
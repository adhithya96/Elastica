#ifndef VARIABLES_H
#define VARIABLES_H
#include "Variables.h"
#include<iomanip>

/*************************************************************
* AceGen    7.505 Windows (16 Aug 22)                        *
*           Co. J. Korelc  2020           19 May 23 14:24:13 *
**************************************************************
User     : Full professional version
Notebook : ContactElement
Evaluation time                 : 1 s     Mode  : Debug
Number of formulae              : 52      Method: Automatic
Subroutine                      : Beam Contact size: 554
Total size of Mathematica  code : 554 subexpressions
Total size of C code            : 3181 bytes */

/*NLEBBE3D::NLEBBE3D(int choice)
{
    if (choice == 1)
    {
        this->NNODE = 5;
        this->NELEM = 2;
        this->NDOF = 6;
        this->NLS = 30;
        this->NEN = 3;
        this->NDM = 3;

        this->NODE = Eigen::MatrixXd::Zero(this->NNODE, this->NDM);
        this->ELEM = Eigen::MatrixXd::Zero(this->NELEM, this->NEN + 1);

        this->NODE(0, 0) = 0;
        this->NODE(0, 1) = 0;
        this->NODE(0, 2) = 0;
        this->NODE(1, 0) = 2.5;
        this->NODE(1, 1) = 0;
        this->NODE(1, 2) = 0;
        this->NODE(2, 0) = 5;
        this->NODE(2, 1) = 0;
        this->NODE(2, 2) = 0;
        this->NODE(3, 0) = 7.5;
        this->NODE(3, 1) = 0;
        this->NODE(3, 2) = 0;
        this->NODE(4, 0) = 10;
        this->NODE(4, 1) = 0;
        this->NODE(4, 2) = 0;

        this->ELEM(0, 0) = 1;
        this->ELEM(0, 1) = 1;
        this->ELEM(0, 2) = 2;
        this->ELEM(0, 3) = 3;
        this->ELEM(1, 0) = 1;
        this->ELEM(1, 1) = 3;
        this->ELEM(1, 2) = 4;
        this->ELEM(1, 3) = 5;

        this->E = 1e7;
        this->nu = 0.0001;
        this->Bp = 1;
        this->Hp = 1;
        this->Zx = 0;
        this->Zy = 0;
        this->Zz = 1;

        this->DIA = 0.5;
    }
    else if (choice == 2)
    {
        this->NNODE = 5;
        this->NELEM = 2;
        this->NDOF = 6;
        this->NLS = 30;
        this->NEN = 3;
        this->NDM = 3;

        this->NODE = Eigen::MatrixXd::Zero(this->NNODE, this->NDM);
        this->ELEM = Eigen::MatrixXd::Zero(this->NELEM, this->NEN + 1);

        this->NODE(0, 0) = 0;
        this->NODE(0, 1) = 2;
        this->NODE(0, 2) = 0;
        this->NODE(1, 0) = 2.5;
        this->NODE(1, 1) = 2;
        this->NODE(1, 2) = 0;
        this->NODE(2, 0) = 5;
        this->NODE(2, 1) = 2;
        this->NODE(2, 2) = 0;
        this->NODE(3, 0) = 7.5;
        this->NODE(3, 1) = 2;
        this->NODE(3, 2) = 0;
        this->NODE(4, 0) = 10;
        this->NODE(4, 1) = 2;
        this->NODE(4, 2) = 0;

        this->ELEM(0, 0) = 1;
        this->ELEM(0, 1) = 1;
        this->ELEM(0, 2) = 2;
        this->ELEM(0, 3) = 3;
        this->ELEM(1, 0) = 1;
        this->ELEM(1, 1) = 3;
        this->ELEM(1, 2) = 4;
        this->ELEM(1, 3) = 5;

        this->E = 1e7;
        this->nu = 0.0001;
        this->Bp = 1;
        this->Hp = 1;
        this->Zx = 0;
        this->Zy = 0;
        this->Zz = 1;

        this->DIA = 0.5;
    }
    else if (choice == 3)
    {
        this->NNODE = 5;
        this->NELEM = 2;
        this->NDOF = 6;
        this->NLS = 30;
        this->NEN = 3;
        this->NDM = 3;

        this->NODE = Eigen::MatrixXd::Zero(this->NNODE, this->NDM);
        this->ELEM = Eigen::MatrixXd::Zero(this->NELEM, this->NEN + 1);

        this->NODE(0, 0) = 0;
        this->NODE(0, 1) = 0;
        this->NODE(0, 2) = 2;
        this->NODE(1, 0) = 2.5;
        this->NODE(1, 1) = 0;
        this->NODE(1, 2) = 2;
        this->NODE(2, 0) = 5;
        this->NODE(2, 1) = 0;
        this->NODE(2, 2) = 2;
        this->NODE(3, 0) = 7.5;
        this->NODE(3, 1) = 0;
        this->NODE(3, 2) = 2;
        this->NODE(4, 0) = 10;
        this->NODE(4, 1) = 0;
        this->NODE(4, 2) = 2;

        this->ELEM(0, 0) = 1;
        this->ELEM(0, 1) = 1;
        this->ELEM(0, 2) = 2;
        this->ELEM(0, 3) = 3;
        this->ELEM(1, 0) = 1;
        this->ELEM(1, 1) = 3;
        this->ELEM(1, 2) = 4;
        this->ELEM(1, 3) = 5;

        this->E = 1e7;
        this->nu = 0.0001;
        this->Bp = 1;
        this->Hp = 1;
        this->Zx = 0;
        this->Zy = 0;
        this->Zz = 1;

        this->DIA = 0.5;
    }
    else if (choice == 4)
    {
        this->NNODE = 5;
        this->NELEM = 2;
        this->NDOF = 6;
        this->NLS = 30;
        this->NEN = 3;
        this->NDM = 3;

        this->NODE = Eigen::MatrixXd::Zero(this->NNODE, this->NDM);
        this->ELEM = Eigen::MatrixXd::Zero(this->NELEM, this->NEN + 1);

        this->NODE(0, 0) = 0;
        this->NODE(0, 1) = 2;
        this->NODE(0, 2) = 2;
        this->NODE(1, 0) = 2.5;
        this->NODE(1, 1) = 2;
        this->NODE(1, 2) = 2;
        this->NODE(2, 0) = 5;
        this->NODE(2, 1) = 2;
        this->NODE(2, 2) = 2;
        this->NODE(3, 0) = 7.5;
        this->NODE(3, 1) = 2;
        this->NODE(3, 2) = 2;
        this->NODE(4, 0) = 10;
        this->NODE(4, 1) = 2;
        this->NODE(4, 2) = 2;

        this->ELEM(0, 0) = 1;
        this->ELEM(0, 1) = 1;
        this->ELEM(0, 2) = 2;
        this->ELEM(0, 3) = 3;
        this->ELEM(1, 0) = 1;
        this->ELEM(1, 1) = 3;
        this->ELEM(1, 2) = 4;
        this->ELEM(1, 3) = 5;

        this->E = 1e7;
        this->nu = 0.0001;
        this->Bp = 1;
        this->Hp = 1;
        this->Zx = 0;
        this->Zy = 0;
        this->Zz = 1;

        this->DIA = 0.5;
    }
}*/

/*NLEBBE3D::NLEBBE3D(int choice, std::string str)
{
    if (choice == 1)
    {
        NLEBBE3D::MAT = str;

        this->NNODE = 5;
        this->NELEM = 2;
        this->NDOF = 6;
        this->NLS = 11;
        this->NEN = 3;
        this->NDM = 3;

        this->NODE = Eigen::MatrixXd::Zero(this->NNODE, this->NDM);
        this->ELEM = Eigen::MatrixXd::Zero(this->NELEM, this->NEN + 1);

        this->NODE(0, 0) = 0;
        this->NODE(0, 1) = 0;
        this->NODE(0, 2) = 0;
        this->NODE(1, 0) = 2.5;
        this->NODE(1, 1) = 0;
        this->NODE(1, 2) = 0;
        this->NODE(2, 0) = 5;
        this->NODE(2, 1) = 0;
        this->NODE(2, 2) = 0;
        this->NODE(3, 0) = 7.5;
        this->NODE(3, 1) = 0;
        this->NODE(3, 2) = 0;
        this->NODE(4, 0) = 10;
        this->NODE(4, 1) = 0;
        this->NODE(4, 2) = 0;

        this->ELEM(0, 0) = 1;
        this->ELEM(0, 1) = 1;
        this->ELEM(0, 2) = 2;
        this->ELEM(0, 3) = 3;
        this->ELEM(1, 0) = 1;
        this->ELEM(1, 1) = 3;
        this->ELEM(1, 2) = 4;
        this->ELEM(1, 3) = 5;

        this->E = 1e7;
        this->nu = 0.0001;
        this->Bp = 1;
        this->Hp = 1;
        this->Zx = 0;
        this->Zy = 0;
        this->Zz = 1;

        this->DIA = 0.5;
    }
    else if (choice == 2)
    {
        NLEBBE3D::MAT = str;

        this->NNODE = 5;
        this->NELEM = 2;
        this->NDOF = 6;
        this->NLS = 11;
        this->NEN = 3;
        this->NDM = 3;

        this->NODE = Eigen::MatrixXd::Zero(this->NNODE, this->NDM);
        this->ELEM = Eigen::MatrixXd::Zero(this->NELEM, this->NEN + 1);

        this->NODE(0, 0) = 0;
        this->NODE(0, 1) = 0.50;
        this->NODE(0, 2) = 0;
        this->NODE(1, 0) = 2.5;
        this->NODE(1, 1) = 0.50;
        this->NODE(1, 2) = 0;
        this->NODE(2, 0) = 5;
        this->NODE(2, 1) = 0.50;
        this->NODE(2, 2) = 0;
        this->NODE(3, 0) = 7.5;
        this->NODE(3, 1) = 0.50;
        this->NODE(3, 2) = 0;
        this->NODE(4, 0) = 10;
        this->NODE(4, 1) = 0.50;
        this->NODE(4, 2) = 0;

        this->ELEM(0, 0) = 1;
        this->ELEM(0, 1) = 1;
        this->ELEM(0, 2) = 2;
        this->ELEM(0, 3) = 3;
        this->ELEM(1, 0) = 1;
        this->ELEM(1, 1) = 3;
        this->ELEM(1, 2) = 4;
        this->ELEM(1, 3) = 5;

        this->E = 1e7;
        this->nu = 0.0001;
        this->Bp = 1;
        this->Hp = 1;
        this->Zx = 0;
        this->Zy = 0;
        this->Zz = 1;

        this->DIA = 0.5;
    }
}*/




#endif

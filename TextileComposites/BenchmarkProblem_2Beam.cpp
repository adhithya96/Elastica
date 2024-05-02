#pragma once

#include "Variables.h"

NLEBBE3D::NLEBBE3D(int choice, std::string str)
{
    if (choice == 1)
    {
        NLEBBE3D::MAT = str;

        NLEBBE3D::NNODE = 21;
        NLEBBE3D::NELEM = 10;
        NLEBBE3D::NDOF = 6;
        NLEBBE3D::NLS = 40;
        NLEBBE3D::NEN = 3;
        NLEBBE3D::NDM = 3;

        NLEBBE3D::NODE = Eigen::MatrixXd::Zero(NLEBBE3D::NNODE, NLEBBE3D::NDM);
        NLEBBE3D::ELEM = Eigen::MatrixXd::Zero(NLEBBE3D::NELEM, NLEBBE3D::NEN + 1);

        //Nodal Information
        //Nodal Information
        NLEBBE3D::NODE(0, 0) = -50;
        NLEBBE3D::NODE(0, 1) = 0;
        NLEBBE3D::NODE(0, 2) = 0;
        NLEBBE3D::NODE(1, 0) = -45;
        NLEBBE3D::NODE(1, 1) = 0;
        NLEBBE3D::NODE(1, 2) = 0;
        NLEBBE3D::NODE(2, 0) = -40;
        NLEBBE3D::NODE(2, 1) = 0;
        NLEBBE3D::NODE(2, 2) = 0;
        NLEBBE3D::NODE(3, 0) = -35;
        NLEBBE3D::NODE(3, 1) = 0;
        NLEBBE3D::NODE(3, 2) = 0;
        NLEBBE3D::NODE(4, 0) = -30;
        NLEBBE3D::NODE(4, 1) = 0;
        NLEBBE3D::NODE(4, 2) = 0;
        NLEBBE3D::NODE(5, 0) = -25;
        NLEBBE3D::NODE(5, 1) = 0;
        NLEBBE3D::NODE(5, 2) = 0;
        NLEBBE3D::NODE(6, 0) = -20;
        NLEBBE3D::NODE(6, 1) = 0;
        NLEBBE3D::NODE(6, 2) = 0;
        NLEBBE3D::NODE(7, 0) = -15;
        NLEBBE3D::NODE(7, 1) = 0;
        NLEBBE3D::NODE(7, 2) = 0;
        NLEBBE3D::NODE(8, 0) = -10;
        NLEBBE3D::NODE(8, 1) = 0;
        NLEBBE3D::NODE(8, 2) = 0;
        NLEBBE3D::NODE(9, 0) = -5;
        NLEBBE3D::NODE(9, 1) = 0;
        NLEBBE3D::NODE(9, 2) = 0;
        NLEBBE3D::NODE(10, 0) = 0;
        NLEBBE3D::NODE(10, 1) = 0;
        NLEBBE3D::NODE(10, 2) = 0;
        NLEBBE3D::NODE(11, 0) = 5;
        NLEBBE3D::NODE(11, 1) = 0;
        NLEBBE3D::NODE(11, 2) = 0;
        NLEBBE3D::NODE(12, 0) = 10;
        NLEBBE3D::NODE(12, 1) = 0;
        NLEBBE3D::NODE(12, 2) = 0;
        NLEBBE3D::NODE(13, 0) = 15;
        NLEBBE3D::NODE(13, 1) = 0;
        NLEBBE3D::NODE(13, 2) = 0;
        NLEBBE3D::NODE(14, 0) = 20;
        NLEBBE3D::NODE(14, 1) = 0;
        NLEBBE3D::NODE(14, 2) = 0;
        NLEBBE3D::NODE(15, 0) = 25;
        NLEBBE3D::NODE(15, 1) = 0;
        NLEBBE3D::NODE(15, 2) = 0;
        NLEBBE3D::NODE(16, 0) = 30;
        NLEBBE3D::NODE(16, 1) = 0;
        NLEBBE3D::NODE(16, 2) = 0;
        NLEBBE3D::NODE(17, 0) = 35;
        NLEBBE3D::NODE(17, 1) = 0;
        NLEBBE3D::NODE(17, 2) = 0;
        NLEBBE3D::NODE(18, 0) = 40;
        NLEBBE3D::NODE(18, 1) = 0;
        NLEBBE3D::NODE(18, 2) = 0;
        NLEBBE3D::NODE(19, 0) = 45;
        NLEBBE3D::NODE(19, 1) = 0;
        NLEBBE3D::NODE(19, 2) = 0;
        NLEBBE3D::NODE(20, 0) = 50;
        NLEBBE3D::NODE(20, 1) = 0;
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
        
        NLEBBE3D::E = 20,000;
        NLEBBE3D::nu = 0.3;
        NLEBBE3D::Bp = 5;
        NLEBBE3D::Hp = 5;
        NLEBBE3D::Zx = 0;
        NLEBBE3D::Zy = 0;
        NLEBBE3D::Zz = 1;

        NLEBBE3D::DIA = 1;
    }
    else if (choice == 2)
    {
        NLEBBE3D::MAT = str;

        NLEBBE3D::NNODE = 21;
        NLEBBE3D::NELEM = 10;
        NLEBBE3D::NDOF = 6;
        NLEBBE3D::NLS = 40;
        NLEBBE3D::NEN = 3;
        NLEBBE3D::NDM = 3;

        NLEBBE3D::NODE = Eigen::MatrixXd::Zero(NLEBBE3D::NNODE, NLEBBE3D::NDM);
        NLEBBE3D::ELEM = Eigen::MatrixXd::Zero(NLEBBE3D::NELEM, NLEBBE3D::NEN + 1);

        //Nodal Information
        NLEBBE3D::NODE(0, 0) = 0;
        NLEBBE3D::NODE(0, 1) = -50;
        NLEBBE3D::NODE(0, 2) = -2.393;
        NLEBBE3D::NODE(1, 0) = 0;
        NLEBBE3D::NODE(1, 1) = -45;
        NLEBBE3D::NODE(1, 2) = -2.393;
        NLEBBE3D::NODE(2, 0) = 0;
        NLEBBE3D::NODE(2, 1) = -40;
        NLEBBE3D::NODE(2, 2) = -2.393;
        NLEBBE3D::NODE(3, 0) = 0;
        NLEBBE3D::NODE(3, 1) = -35;
        NLEBBE3D::NODE(3, 2) = -2.393;
        NLEBBE3D::NODE(4, 0) = 0;
        NLEBBE3D::NODE(4, 1) = -30;
        NLEBBE3D::NODE(4, 2) = -2.393;
        NLEBBE3D::NODE(5, 0) = 0;
        NLEBBE3D::NODE(5, 1) = -25;
        NLEBBE3D::NODE(5, 2) = -2.393;
        NLEBBE3D::NODE(6, 0) = 0;
        NLEBBE3D::NODE(6, 1) = -20;
        NLEBBE3D::NODE(6, 2) = -2.393;
        NLEBBE3D::NODE(7, 0) = 0;
        NLEBBE3D::NODE(7, 1) = -15;
        NLEBBE3D::NODE(7, 2) = -2.393;
        NLEBBE3D::NODE(8, 0) = 0;
        NLEBBE3D::NODE(8, 1) = -10;
        NLEBBE3D::NODE(8, 2) = -2.393;
        NLEBBE3D::NODE(9, 0) = 0;
        NLEBBE3D::NODE(9, 1) = -5;
        NLEBBE3D::NODE(9, 2) = -2.393;
        NLEBBE3D::NODE(10, 0) = 0;
        NLEBBE3D::NODE(10, 1) = 0;
        NLEBBE3D::NODE(10, 2) = -2.393;
        NLEBBE3D::NODE(11, 0) = 0;
        NLEBBE3D::NODE(11, 1) = 5;
        NLEBBE3D::NODE(11, 2) = -2.393;
        NLEBBE3D::NODE(12, 0) = 0;
        NLEBBE3D::NODE(12, 1) = 10;
        NLEBBE3D::NODE(12, 2) = -2.393;
        NLEBBE3D::NODE(13, 0) = 0;
        NLEBBE3D::NODE(13, 1) = 15;
        NLEBBE3D::NODE(13, 2) = -2.393;
        NLEBBE3D::NODE(14, 0) = 0;
        NLEBBE3D::NODE(14, 1) = 20;
        NLEBBE3D::NODE(14, 2) = -2.393;
        NLEBBE3D::NODE(15, 0) = 0;
        NLEBBE3D::NODE(15, 1) = 25;
        NLEBBE3D::NODE(15, 2) = -2.393;
        NLEBBE3D::NODE(16, 0) = 0;
        NLEBBE3D::NODE(16, 1) = 30;
        NLEBBE3D::NODE(16, 2) = -2.393;
        NLEBBE3D::NODE(17, 0) = 0;
        NLEBBE3D::NODE(17, 1) = 35;
        NLEBBE3D::NODE(17, 2) = -2.393;
        NLEBBE3D::NODE(18, 0) = 0;
        NLEBBE3D::NODE(18, 1) = 40;
        NLEBBE3D::NODE(18, 2) = -2.393;
        NLEBBE3D::NODE(19, 0) = 0;
        NLEBBE3D::NODE(19, 1) = 45;
        NLEBBE3D::NODE(19, 2) = -2.393;
        NLEBBE3D::NODE(20, 0) = 0;
        NLEBBE3D::NODE(20, 1) = 50;
        NLEBBE3D::NODE(20, 2) = -2.393;

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

        NLEBBE3D::E = 30,000;
        NLEBBE3D::nu = 0.17;
        NLEBBE3D::Bp = 5;
        NLEBBE3D::Hp = 5;
        NLEBBE3D::Zx = 0;
        NLEBBE3D::Zy = 0;
        NLEBBE3D::Zz = 1;

        NLEBBE3D::DIA = 1;
    }
}

//
//Loading::Loading(int beamnum)
//{
//    if (beamnum == 1)
//    {
//        this->ELEM = Eigen::MatrixXd::Zero(0, 0);
//        //this->ELEM = Eigen::MatrixXd::Zero(1, 4);
//        //node
//        //this->ELEM(0, 0) = 21;
//        ////type
//        //this->ELEM(0, 1) = -2;
//        ////direction
//        //this->ELEM(0, 2) = -3;
//        ////value
//        //this->ELEM(0, 3) = 1;
//    }
//    else if (beamnum == 2)
//    {
//        //this->ELEM = Eigen::MatrixXd::Zero(0, 0);
//        this->ELEM = Eigen::MatrixXd::Zero(0, 0);
//        ////node
//        //this->ELEM(0, 0) = 5;
//        ////type
//        //this->ELEM(0, 1) = -1;
//        ////direction
//        //this->ELEM(0, 2) = -2;
//        ////value
//        //this->ELEM(0, 3) = 0.3;
//    }
//}
//
//Boundary::Boundary(int beamnum)
//{
//    if (beamnum == 2)
//    {
//        this->ELEM = Eigen::MatrixXd::Zero(1, 2);
//        this->ELEM(0, 0) = 1;
//        this->ELEM(0, 1) = -1;
//    }
//    else if (beamnum == 1)
//    {
//        this->ELEM = Eigen::MatrixXd::Zero(2, 4);
//        this->ELEM(0, 0) = 1;
//        this->ELEM(0, 1) = -1;
//        this->ELEM(1, 0) = 21;
//        this->ELEM(1, 1) = -3;
//        this->ELEM(1, 2) = -3;
//        this->ELEM(1, 3) = 30;
//    }
//}

//BeamContact::BeamContact(const int nbeams, std::string str)
//{
//    this->NBEAMS = nbeams;
//
//    this->epsilon = 50000;
//
//    if (str == "STS")
//        this->NEN = 4;
//    else if (str == "NTN")
//        this->NEN = 2;
//
//    this->GlobalELEM = Eigen::MatrixXd(1, 2);
//    for (int i = 0; i < this->GlobalELEM.rows(); i++)
//    {
//        for (int j = 0; j < this->GlobalELEM.cols(); j++)
//            this->GlobalELEM(i, j) = -1;
//    }
//
//    this->n = Eigen::VectorXd(3);
//    this->n(0) = 0;
//    this->n(1) = 0;
//    this->n(2) = -1;
//
//    this->start = 1;
//
//    this->precontactiter = 0;
//}


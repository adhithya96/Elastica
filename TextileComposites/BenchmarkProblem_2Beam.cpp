#pragma once

#include "Variables.h"

NLEBBE3D::NLEBBE3D(int choice, std::string str)
{
    if (choice == 1)
    {
        NLEBBE3D::MAT = str;

        NLEBBE3D::NELEM = 10;
        NLEBBE3D::NDOF = 6;
        NLEBBE3D::NLS = 80;
        NLEBBE3D::NEN = 3;
        NLEBBE3D::NDM = 3;

        //Nodal Information
        //Nodal Information
        NNODE = NELEM * 2 + 1;

        NLEBBE3D::NODE = Eigen::MatrixXd::Zero(NLEBBE3D::NNODE, NLEBBE3D::NDM);
        NLEBBE3D::ELEM = Eigen::MatrixXd::Zero(NLEBBE3D::NELEM, NLEBBE3D::NEN + 1);

        for (int i = 0; i < NNODE; i++)
        {
            NODE(i, 0) = i * 100 / (NNODE - 1);
            NODE(i, 1) = 0;
            NODE(i, 2) = 0;
        }

        for (int i = 0; i < NELEM; i++)
        {
            ELEM(i, 0) = 1;
            ELEM(i, 1) = 2 * i + 1;
            ELEM(i, 2) = 2 * i + 2;
            ELEM(i, 3) = 2 * i + 3;
        }
        
        NLEBBE3D::E = 20000;
        NLEBBE3D::nu = 0.3;
        NLEBBE3D::Bp = 5;
        NLEBBE3D::Hp = 5;
        NLEBBE3D::Zx = 0;
        NLEBBE3D::Zy = 0;
        NLEBBE3D::Zz = 1;

        NLEBBE3D::DIA = 5;
    }
    else if (choice == 2)
    {
        NLEBBE3D::MAT = str;

        NLEBBE3D::NELEM = 10;
        NLEBBE3D::NDOF = 6;
        NLEBBE3D::NLS = 80;
        NLEBBE3D::NEN = 3;
        NLEBBE3D::NDM = 3;

        NNODE = NELEM * 2 + 1;

        NLEBBE3D::NODE = Eigen::MatrixXd::Zero(NLEBBE3D::NNODE, NLEBBE3D::NDM);
        NLEBBE3D::ELEM = Eigen::MatrixXd::Zero(NLEBBE3D::NELEM, NLEBBE3D::NEN + 1);

        for (int i = 0; i < NNODE; i++)
        {
            NODE(i, 0) = i * 100 / (NNODE - 1);
            NODE(i, 1) = 10;
            NODE(i, 2) = 0;
        }

        for (int i = 0; i < NELEM; i++)
        {
            ELEM(i, 0) = 1;
            ELEM(i, 1) = 2 * i + 1;
            ELEM(i, 2) = 2 * i + 2;
            ELEM(i, 3) = 2 * i + 3;
        }


        NLEBBE3D::E = 20000;
        NLEBBE3D::nu = 0.3;
        NLEBBE3D::Bp = 5;
        NLEBBE3D::Hp = 5;
        NLEBBE3D::Zx = 0;
        NLEBBE3D::Zy = 0;
        NLEBBE3D::Zz = 1;

        NLEBBE3D::DIA = 5;
    }
}
//
//
//Loading::Loading(int beamnum)
//{
//    if (beamnum == 1)
//    {
//        //this->ELEM = Eigen::MatrixXd::Zero(0, 0);
//        this->ELEM = Eigen::MatrixXd::Zero(1, 4);
//        //node
//        this->ELEM(0, 0) = 1;
//        //type
//        this->ELEM(0, 1) = -1;
//        //direction
//        this->ELEM(0, 2) = 2;
//        //value
//        this->ELEM(0, 3) = 100;
//    }
//    else if (beamnum == 2)
//    {
//        //this->ELEM = Eigen::MatrixXd::Zero(0, 0);
//        this->ELEM = Eigen::MatrixXd::Zero(1, 4);
//        //node
//        this->ELEM(0, 0) = 1;
//        //type
//        this->ELEM(0, 1) = -1;
//        //direction
//        this->ELEM(0, 2) = -2;
//        //value
//        this->ELEM(0, 3) = 100;
//    }
//}
//
//Boundary::Boundary(int beamnum)
//{
//    if (beamnum == 1)
//    {
//        this->ELEM = Eigen::MatrixXd::Zero(1, 2);
//        this->ELEM(0, 0) = 1;
//        this->ELEM(0, 1) = -1;
//    }
//    else if (beamnum == 2)
//    {
//        this->ELEM = Eigen::MatrixXd::Zero(1, 2);
//        this->ELEM(0, 0) = 1;
//        this->ELEM(0, 1) = -1;
//    }
//}
//
//BeamContact::BeamContact(const int nbeams, std::string str)
//{
//    this->NBEAMS = nbeams;
//
//    this->epsilon = 1e10;
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
//    this->n(1) = 1;
//    this->n(2) = 0;
//
//    this->start = 1;
//
//    this->precontactiter = 0;
//}
//

#pragma once

#include "Variables.h"

NLEBBE3D::NLEBBE3D(int nele, double bd, double bl, int nls, int nen, int ndm, double E, double nu, int beamnum)
{
    int nnode = nele * 2 + 1;

    this->NNODE = nnode;
    this->NELEM = nele;
    this->NDOF = 6;
    this->NLS = nls;
    this->NEN = nen;
    this->NDM = ndm;

    this->NODE = Eigen::MatrixXd::Zero(nnode, ndm);
    this->ELEM = Eigen::MatrixXd::Zero(this->NELEM, nen + 1);

    this->E = E;
    this->nu = nu;
    this->Bp = 1;
    this->Hp = 1;
    this->Zx = 0;
    this->Zy = 0;
    this->Zz = 1;

    this->DIA = bd;
    
    int sign;
    if (beamnum == 0)
        sign = -1;
    else
        sign = 1;

    for (int i = 0; i < nnode; i++)
    {
        this->NODE(i, 0) = i * bl / (nnode - 1);
        this->NODE(i, 1) = sign * bd / 2 + sign * 0.001;
        this->NODE(i, 2) = 0;
    }
    for (int i = 0; i < this->NELEM; i++)
    {
        this->ELEM(i, 0) = 1;
        this->ELEM(i, 1) = 2 * i + 1;
        this->ELEM(i, 2) = 2 * i + 2;
        this->ELEM(i, 3) = 2 * i + 3;
    }
}

//BeamContact::BeamContact(const int nbeams, std::string str)
//{
//    this->NBEAMS = nbeams;
//
//    this->epsilon = pow(10, 10);
//
//    if (str == "STS")
//        this->NEN = 4;
//    else if (str == "NTN")
//        this->NEN = 2;
//
//    this->n = Eigen::VectorXd(3);
//    this->n(0) = 0;
//    this->n(1) = -1;
//    this->n(2) = 0;
//
//    this->start = 1;
//
//    this->precontactiter = 0;
//}

//Loading::Loading(int beamnum)
//{
//    if (beamnum == 1)
//    {
//        this->ELEM = Eigen::MatrixXd::Zero(1, 4);
//        //node
//        this->ELEM(0, 0) = 13;
//        //type
//        this->ELEM(0, 1) = -1;
//        //direction
//        this->ELEM(0, 2) = -2;
//        //value
//        this->ELEM(0, 3) = pow(10, 4);
//    }
//    else if (beamnum == 2)
//    {
//
//        this->ELEM = Eigen::MatrixXd::Zero(1, 4);
//        //node
//        this->ELEM(0, 0) = 13;
//        //type
//        this->ELEM(0, 1) = -1;
//        //direction
//        this->ELEM(0, 2) = 2;
//        //value
//        this->ELEM(0, 3) = pow(10, 4);
//    }
//}
//
//Boundary::Boundary(int beamnum)
//{
//    if (beamnum == 1)
//    {
//        this->ELEM = Eigen::MatrixXd::Zero(2, 2);
//        this->ELEM(0, 0) = 1;
//        //type
//        this->ELEM(0, 1) = -1;
//
//        this->ELEM(1, 0) = 81;
//        //type
//        this->ELEM(1, 1) = -1;
//        //this->ELEM(0, 0) = 1;
//        //type
//        //this->ELEM(0, 1) = -1;
//
//        //this->ELEM(1, 0) = 13;
//        //this->ELEM(1, 1) = -3;
//        //this->ELEM(1, 2) = 2;
//        //this->ELEM(1, 3) = 0.15;
//    }
//    else if (beamnum == 2)
//    {
//        this->ELEM = Eigen::MatrixXd::Zero(2, 2);
//        this->ELEM(0, 0) = 1;
//        //type
//        this->ELEM(0, 1) = -1;
//
//        this->ELEM(1, 0) = 81;
//        //type
//        this->ELEM(1, 1) = -1;
//        //this->ELEM(1, 0) = 13;
//        //this->ELEM(1, 1) = -3;
//        //this->ELEM(1, 2) = 2;
//        //this->ELEM(1, 3) = 0.15;
//    }
//}

//NLEBBE3D::NLEBBE3D(int choice, std::string str)
//{
//    if (choice == 1)
//    {
//        NLEBBE3D::MAT = str;
//
//        NLEBBE3D::NNODE = 5;
//        NLEBBE3D::NELEM = 2;
//        NLEBBE3D::NDOF = 6;
//        NLEBBE3D::NLS = 20;
//        NLEBBE3D::NEN = 3;
//        NLEBBE3D::NDM = 3;
//
//        NLEBBE3D::NODE = Eigen::MatrixXd::Zero(NLEBBE3D::NNODE, NLEBBE3D::NDM);
//        NLEBBE3D::ELEM = Eigen::MatrixXd::Zero(NLEBBE3D::NELEM, NLEBBE3D::NEN + 1);
//
//        //Nodal Information
//        //Nodal Information
//        NLEBBE3D::NODE(0, 0) = -5;
//        NLEBBE3D::NODE(0, 1) = -0.1;
//        NLEBBE3D::NODE(0, 2) = 0;
//        NLEBBE3D::NODE(1, 0) = -4;
//        NLEBBE3D::NODE(1, 1) = -0.5;
//        NLEBBE3D::NODE(1, 2) = 0;
//        NLEBBE3D::NODE(2, 0) = -3;
//        NLEBBE3D::NODE(2, 1) = -0.5;
//        NLEBBE3D::NODE(2, 2) = 0;
//        NLEBBE3D::NODE(3, 0) = -2;
//        NLEBBE3D::NODE(3, 1) = -0.1;
//        NLEBBE3D::NODE(3, 2) = 0;
//        NLEBBE3D::NODE(4, 0) = -1;
//        NLEBBE3D::NODE(4, 1) = -0.5;
//        NLEBBE3D::NODE(4, 2) = 0;
//        NLEBBE3D::NODE(5, 0) = 0;
//        NLEBBE3D::NODE(5, 1) = -0.5;
//        NLEBBE3D::NODE(5, 2) = 0;
//        NLEBBE3D::NODE(6, 0) = -5;
//        NLEBBE3D::NODE(6, 1) = -0.1;
//        NLEBBE3D::NODE(6, 2) = 0;
//        NLEBBE3D::NODE(7, 0) = -4;
//        NLEBBE3D::NODE(7, 1) = -0.5;
//        NLEBBE3D::NODE(7, 2) = 0;
//        NLEBBE3D::NODE(8, 0) = -3;
//        NLEBBE3D::NODE(8, 1) = -0.5;
//        NLEBBE3D::NODE(8, 2) = 0;
//
//
//
//        //Element connectivity
//        NLEBBE3D::ELEM(0, 0) = 1;
//        NLEBBE3D::ELEM(0, 1) = 1;
//        NLEBBE3D::ELEM(0, 2) = 2;
//        NLEBBE3D::ELEM(0, 3) = 3;
//
//        NLEBBE3D::ELEM(1, 0) = 1;
//        NLEBBE3D::ELEM(1, 1) = 3;
//        NLEBBE3D::ELEM(1, 2) = 4;
//        NLEBBE3D::ELEM(1, 3) = 5;
//
//        NLEBBE3D::E = 2.4e3;
//        NLEBBE3D::nu = 0.39;
//        NLEBBE3D::Bp = 1;
//        NLEBBE3D::Hp = 1;
//        NLEBBE3D::Zx = 0;
//        NLEBBE3D::Zy = 0;
//        NLEBBE3D::Zz = 1;
//
//        NLEBBE3D::DIA = 0.1;
//    }
//    else if (choice == 2)
//    {
//        NLEBBE3D::MAT = str;
//
//        NLEBBE3D::NNODE = 5;
//        NLEBBE3D::NELEM = 2;
//        NLEBBE3D::NDOF = 6;
//        NLEBBE3D::NLS = 20;
//        NLEBBE3D::NEN = 3;
//        NLEBBE3D::NDM = 3;
//
//        NLEBBE3D::NODE = Eigen::MatrixXd::Zero(NLEBBE3D::NNODE, NLEBBE3D::NDM);
//        NLEBBE3D::ELEM = Eigen::MatrixXd::Zero(NLEBBE3D::NELEM, NLEBBE3D::NEN + 1);
//
//        //Nodal Information
//        NLEBBE3D::NODE(0, 0) = 0;
//        NLEBBE3D::NODE(0, 1) = 0.5;
//        NLEBBE3D::NODE(0, 2) = -2;
//        NLEBBE3D::NODE(1, 0) = 0;
//        NLEBBE3D::NODE(1, 1) = 0.5;
//        NLEBBE3D::NODE(1, 2) = -1;
//        NLEBBE3D::NODE(2, 0) = 0;
//        NLEBBE3D::NODE(2, 1) = 0.5;
//        NLEBBE3D::NODE(2, 2) = 0;
//        NLEBBE3D::NODE(3, 0) = 0;
//        NLEBBE3D::NODE(3, 1) = 0.5;
//        NLEBBE3D::NODE(3, 2) = 1;
//        NLEBBE3D::NODE(4, 0) = 0;
//        NLEBBE3D::NODE(4, 1) = 0.5;
//        NLEBBE3D::NODE(4, 2) = 2;
//
//        //Element connectivity
//        NLEBBE3D::ELEM(0, 0) = 1;
//        NLEBBE3D::ELEM(0, 1) = 1;
//        NLEBBE3D::ELEM(0, 2) = 2;
//        NLEBBE3D::ELEM(0, 3) = 3;
//
//        NLEBBE3D::ELEM(1, 0) = 1;
//        NLEBBE3D::ELEM(1, 1) = 3;
//        NLEBBE3D::ELEM(1, 2) = 4;
//        NLEBBE3D::ELEM(1, 3) = 5;
//
//        NLEBBE3D::E = 2.4e3;
//        NLEBBE3D::nu = 0.39;
//        NLEBBE3D::Bp = 1;
//        NLEBBE3D::Hp = 1;
//        NLEBBE3D::Zx = 0;
//        NLEBBE3D::Zy = 1;
//        NLEBBE3D::Zz = 0;
//
//        NLEBBE3D::DIA = 0.1;
//    }
//}

/*NLEBBE3D::NLEBBE3D(int choice, std::string str)
{
    if (choice == 1)
    {
        NLEBBE3D::MAT = str;

        NLEBBE3D::NNODE = 21;
        NLEBBE3D::NELEM = 10;
        NLEBBE3D::NDOF = 6;
        NLEBBE3D::NLS = 20;
        NLEBBE3D::NEN = 3;
        NLEBBE3D::NDM = 3;

        NLEBBE3D::NODE = Eigen::MatrixXd::Zero(NLEBBE3D::NNODE, NLEBBE3D::NDM);
        NLEBBE3D::ELEM = Eigen::MatrixXd::Zero(NLEBBE3D::NELEM, NLEBBE3D::NEN + 1);

        //Nodal Information
        NLEBBE3D::NODE(0, 0) = -250;
        NLEBBE3D::NODE(0, 1) = -83.33;
        NLEBBE3D::NODE(0, 2) = 0;
        NLEBBE3D::NODE(1, 0) = -224.785;
        NLEBBE3D::NODE(1, 1) = -83.33;
        NLEBBE3D::NODE(1, 2) = -2.3709;
        NLEBBE3D::NODE(2, 0) = -199.652;
        NLEBBE3D::NODE(2, 1) = -83.33;
        NLEBBE3D::NODE(2, 2) = -5.49695;
        NLEBBE3D::NODE(3, 0) = -174.452;
        NLEBBE3D::NODE(3, 1) = -83.33;
        NLEBBE3D::NODE(3, 2) = -8.02484;
        NLEBBE3D::NODE(4, 0) = -149.309;
        NLEBBE3D::NODE(4, 1) = -83.33;
        NLEBBE3D::NODE(4, 2) = -11.0731;
        NLEBBE3D::NODE(5, 0) = -124.113;
        NLEBBE3D::NODE(5, 1) = -83.33;
        NLEBBE3D::NODE(5, 2) = -13.6469;
        NLEBBE3D::NODE(6, 0) = -98.9726;
        NLEBBE3D::NODE(6, 1) = -83.33;
        NLEBBE3D::NODE(6, 2) = -16.7119;
        NLEBBE3D::NODE(7, 0) = -73.7445;
        NLEBBE3D::NODE(7, 1) = -83.33;
        NLEBBE3D::NODE(7, 2) = -18;
        NLEBBE3D::NODE(8, 0) = -48.9798;
        NLEBBE3D::NODE(8, 1) = -83.33;
        NLEBBE3D::NODE(8, 2) = -12.9843;
        NLEBBE3D::NODE(9, 0) = -24.5143;
        NLEBBE3D::NODE(9, 1) = -83.33;
        NLEBBE3D::NODE(9, 2) = -6.4;
        NLEBBE3D::NODE(10, 0) = 0;
        NLEBBE3D::NODE(10, 1) = -83.33;
        NLEBBE3D::NODE(10, 2) = 0;
        NLEBBE3D::NODE(11, 0) = 24.5143;
        NLEBBE3D::NODE(11, 1) = -83.33;
        NLEBBE3D::NODE(11, 2) = 6.4;
        NLEBBE3D::NODE(12, 0) = 48.9798;
        NLEBBE3D::NODE(12, 1) = -83.33;
        NLEBBE3D::NODE(12, 2) = 12.9843;
        NLEBBE3D::NODE(13, 0) = 73.7445;
        NLEBBE3D::NODE(13, 1) = -83.33;
        NLEBBE3D::NODE(13, 2) = 18;
        NLEBBE3D::NODE(14, 0) = 98.9726;
        NLEBBE3D::NODE(14, 1) = -83.33;
        NLEBBE3D::NODE(14, 2) = 16.7119;
        NLEBBE3D::NODE(15, 0) = 124.113;
        NLEBBE3D::NODE(15, 1) = -83.33;
        NLEBBE3D::NODE(15, 2) = 13.6469;
        NLEBBE3D::NODE(16, 0) = 149.309;
        NLEBBE3D::NODE(16, 1) = -83.33;
        NLEBBE3D::NODE(16, 2) = 11.0731;
        NLEBBE3D::NODE(17, 0) = 174.452;
        NLEBBE3D::NODE(17, 1) = -83.33;
        NLEBBE3D::NODE(17, 2) = 8.02484;
        NLEBBE3D::NODE(18, 0) = 199.652;
        NLEBBE3D::NODE(18, 1) = -83.33;
        NLEBBE3D::NODE(18, 2) = 5.49695;
        NLEBBE3D::NODE(19, 0) = 224.785;
        NLEBBE3D::NODE(19, 1) = -83.33;
        NLEBBE3D::NODE(19, 2) = 2.37094;
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
        NLEBBE3D::Bp = 30;
        NLEBBE3D::Hp = 30;
        NLEBBE3D::Zx = 0;
        NLEBBE3D::Zy = 0;
        NLEBBE3D::Zz = 1;

        NLEBBE3D::DIA = 30;
    }
    else if (choice == 2)
    {
        NLEBBE3D::MAT = str;

        NLEBBE3D::NNODE = 21;
        NLEBBE3D::NELEM = 10;
        NLEBBE3D::NDOF = 6;
        NLEBBE3D::NLS = 20;
        NLEBBE3D::NEN = 3;
        NLEBBE3D::NDM = 3;

        NLEBBE3D::NODE = Eigen::MatrixXd::Zero(NLEBBE3D::NNODE, NLEBBE3D::NDM);
        NLEBBE3D::ELEM = Eigen::MatrixXd::Zero(NLEBBE3D::NELEM, NLEBBE3D::NEN + 1);

        //Nodal Information
        NLEBBE3D::NODE(0, 0) = -250;
        NLEBBE3D::NODE(0, 1) = 83.33;
        NLEBBE3D::NODE(0, 2) = 0;
        NLEBBE3D::NODE(1, 0) = -224.785;
        NLEBBE3D::NODE(1, 1) = 83.33;
        NLEBBE3D::NODE(1, 2) = 2.3709;
        NLEBBE3D::NODE(2, 0) = -199.652;
        NLEBBE3D::NODE(2, 1) = 83.33;
        NLEBBE3D::NODE(2, 2) = 5.49695;
        NLEBBE3D::NODE(3, 0) = -174.452;
        NLEBBE3D::NODE(3, 1) = 83.33;
        NLEBBE3D::NODE(3, 2) = 8.02484;
        NLEBBE3D::NODE(4, 0) = -149.309;
        NLEBBE3D::NODE(4, 1) = 83.33;
        NLEBBE3D::NODE(4, 2) = 11.0731;
        NLEBBE3D::NODE(5, 0) = -124.113;
        NLEBBE3D::NODE(5, 1) = 83.33;
        NLEBBE3D::NODE(5, 2) = 13.6469;
        NLEBBE3D::NODE(6, 0) = -98.9726;
        NLEBBE3D::NODE(6, 1) = 83.33;
        NLEBBE3D::NODE(6, 2) = 16.7119;
        NLEBBE3D::NODE(7, 0) = -73.7445;
        NLEBBE3D::NODE(7, 1) = 83.33;
        NLEBBE3D::NODE(7, 2) = 18;
        NLEBBE3D::NODE(8, 0) = -48.9798;
        NLEBBE3D::NODE(8, 1) = 83.33;
        NLEBBE3D::NODE(8, 2) = 12.9843;
        NLEBBE3D::NODE(9, 0) = -24.5143;
        NLEBBE3D::NODE(9, 1) = 83.33;
        NLEBBE3D::NODE(9, 2) = 6.4;
        NLEBBE3D::NODE(10, 0) = 0;
        NLEBBE3D::NODE(10, 1) = 83.33;
        NLEBBE3D::NODE(10, 2) = 0;
        NLEBBE3D::NODE(11, 0) = 24.5143;
        NLEBBE3D::NODE(11, 1) = 83.33;
        NLEBBE3D::NODE(11, 2) = -6.4;
        NLEBBE3D::NODE(12, 0) = 48.9798;
        NLEBBE3D::NODE(12, 1) = 83.33;
        NLEBBE3D::NODE(12, 2) = -12.9843;
        NLEBBE3D::NODE(13, 0) = 73.7445;
        NLEBBE3D::NODE(13, 1) = 83.33;
        NLEBBE3D::NODE(13, 2) = -18;
        NLEBBE3D::NODE(14, 0) = 98.9726;
        NLEBBE3D::NODE(14, 1) = 83.33;
        NLEBBE3D::NODE(14, 2) = -16.7119;
        NLEBBE3D::NODE(15, 0) = 124.113;
        NLEBBE3D::NODE(15, 1) = 83.33;
        NLEBBE3D::NODE(15, 2) = -13.6469;
        NLEBBE3D::NODE(16, 0) = 149.309;
        NLEBBE3D::NODE(16, 1) = 83.33;
        NLEBBE3D::NODE(16, 2) = -11.0731;
        NLEBBE3D::NODE(17, 0) = 174.452;
        NLEBBE3D::NODE(17, 1) = 83.33;
        NLEBBE3D::NODE(17, 2) = -8.02484;
        NLEBBE3D::NODE(18, 0) = 199.652;
        NLEBBE3D::NODE(18, 1) = 83.33;
        NLEBBE3D::NODE(18, 2) = -5.49695;
        NLEBBE3D::NODE(19, 0) = 224.785;
        NLEBBE3D::NODE(19, 1) = 83.33;
        NLEBBE3D::NODE(19, 2) = -2.37094;
        NLEBBE3D::NODE(20, 0) = 250;
        NLEBBE3D::NODE(20, 1) = 83.33;
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
        NLEBBE3D::Bp = 30;
        NLEBBE3D::Hp = 30;
        NLEBBE3D::Zx = 0;
        NLEBBE3D::Zy = 0;
        NLEBBE3D::Zz = 1;

        NLEBBE3D::DIA = 30;
    }
    else if (choice == 3)
    {
        NLEBBE3D::MAT = str;

        NLEBBE3D::NNODE = 21;
        NLEBBE3D::NELEM = 10;
        NLEBBE3D::NDOF = 6;
        NLEBBE3D::NLS = 20;
        NLEBBE3D::NEN = 3;
        NLEBBE3D::NDM = 3;

        NLEBBE3D::NODE = Eigen::MatrixXd::Zero(NLEBBE3D::NNODE, NLEBBE3D::NDM);
        NLEBBE3D::ELEM = Eigen::MatrixXd::Zero(NLEBBE3D::NELEM, NLEBBE3D::NEN + 1);

        //Nodal Information
        NLEBBE3D::NODE(0, 0) = -83.33;
        NLEBBE3D::NODE(0, 1) = -250;
        NLEBBE3D::NODE(0, 2) = 0;
        NLEBBE3D::NODE(1, 0) = -83.33;
        NLEBBE3D::NODE(1, 1) = -224.785;
        NLEBBE3D::NODE(1, 2) = 2.37094;
        NLEBBE3D::NODE(2, 0) = -83.33;
        NLEBBE3D::NODE(2, 1) = -199.652;
        NLEBBE3D::NODE(2, 2) = 5.49695;
        NLEBBE3D::NODE(3, 0) = -83.33;
        NLEBBE3D::NODE(3, 1) = -174.452;
        NLEBBE3D::NODE(3, 2) = 8.02484;
        NLEBBE3D::NODE(4, 0) = -83.33;
        NLEBBE3D::NODE(4, 1) = -149.309;
        NLEBBE3D::NODE(4, 2) = 11.0731;
        NLEBBE3D::NODE(5, 0) = -83.33;
        NLEBBE3D::NODE(5, 1) = -124.113;
        NLEBBE3D::NODE(5, 2) = 13.6469;
        NLEBBE3D::NODE(6, 0) = -83.33;
        NLEBBE3D::NODE(6, 1) = -98.973;
        NLEBBE3D::NODE(6, 2) = 16.7119;
        NLEBBE3D::NODE(7, 0) = -83.33;
        NLEBBE3D::NODE(7, 1) = -73.744;
        NLEBBE3D::NODE(7, 2) = 18;
        NLEBBE3D::NODE(8, 0) = -83.33;
        NLEBBE3D::NODE(8, 1) = -48.9798;
        NLEBBE3D::NODE(8, 2) = 12.9843;
        NLEBBE3D::NODE(9, 0) = -83.33;
        NLEBBE3D::NODE(9, 1) = -24.5143;
        NLEBBE3D::NODE(9, 2) = 6.4;
        NLEBBE3D::NODE(10, 0) = -83.33;
        NLEBBE3D::NODE(10, 1) = 0;
        NLEBBE3D::NODE(10, 2) = 0;
        NLEBBE3D::NODE(11, 0) = -83.33;
        NLEBBE3D::NODE(11, 1) = 24.5143;
        NLEBBE3D::NODE(11, 2) = -6.4;
        NLEBBE3D::NODE(12, 0) = -83.33;
        NLEBBE3D::NODE(12, 1) = 48.9798;
        NLEBBE3D::NODE(12, 2) = -12.9843;
        NLEBBE3D::NODE(13, 0) = -83.33;
        NLEBBE3D::NODE(13, 1) = 73.7445;
        NLEBBE3D::NODE(13, 2) = -18;
        NLEBBE3D::NODE(14, 0) = -83.33;
        NLEBBE3D::NODE(14, 1) = 98.9726;
        NLEBBE3D::NODE(14, 2) = -16.7119;
        NLEBBE3D::NODE(15, 0) = -83.33;
        NLEBBE3D::NODE(15, 1) = 124.113;
        NLEBBE3D::NODE(15, 2) = -13.6469;
        NLEBBE3D::NODE(16, 0) = -83.33;
        NLEBBE3D::NODE(16, 1) = 149.309;
        NLEBBE3D::NODE(16, 2) = -11.0731;
        NLEBBE3D::NODE(17, 0) = -83.33;
        NLEBBE3D::NODE(17, 1) = 174.452;
        NLEBBE3D::NODE(17, 2) = -8.02484;
        NLEBBE3D::NODE(18, 0) = -83.33;
        NLEBBE3D::NODE(18, 1) = 199.652;
        NLEBBE3D::NODE(18, 2) = -5.49695;
        NLEBBE3D::NODE(19, 0) = -83.33;
        NLEBBE3D::NODE(19, 1) = 224.785;
        NLEBBE3D::NODE(19, 2) = -2.37094;
        NLEBBE3D::NODE(20, 0) = -83.33;
        NLEBBE3D::NODE(20, 1) = 250;
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
        NLEBBE3D::Bp = 30;
        NLEBBE3D::Hp = 30;
        NLEBBE3D::Zx = 0;
        NLEBBE3D::Zy = 0;
        NLEBBE3D::Zz = 1;

        NLEBBE3D::DIA = 30;
    }

}

Loading::Loading(int beamnum, int contactswitch)
{
    if (contactswitch == 1)
    {
        if (beamnum == 1)
        {
            this->ELEM = Eigen::MatrixXd::Zero(0, 0);
        }
        else if (beamnum == 1)
        {
            this->ELEM = Eigen::MatrixXd::Zero(0, 0);
        }
        else if (beamnum == 3)
        {
            //this->ELEM = Eigen::MatrixXd::Zero(0, 0);
            this->ELEM = Eigen::MatrixXd::Zero(1, 4);
            //node
            this->ELEM(0, 0) = 1;
            //type
            this->ELEM(0, 1) = -2;
            //direction
            this->ELEM(0, 2) = -2;
            //value
            this->ELEM(0, 3) = 100000;
        }
    }
}

Boundary::Boundary(int beamnum, int contactswitch)
{
    if (contactswitch == 1)
    {
        if (beamnum == 1)
        {
            this->ELEM = Eigen::MatrixXd::Zero(2, 2);
            this->ELEM(0, 0) = 21;
            this->ELEM(0, 1) = -1;
            this->ELEM(0, 0) = 1;
            this->ELEM(0, 1) = -1;
        }
        else if (beamnum == 2)
        {
            this->ELEM = Eigen::MatrixXd::Zero(2, 2);
            this->ELEM(0, 0) = 21;
            this->ELEM(0, 1) = -1;
            this->ELEM(0, 0) = 1;
            this->ELEM(0, 1) = -1;
        }
        else if (beamnum == 3)
        {
            this->ELEM = Eigen::MatrixXd::Zero(1, 2);
            this->ELEM(0, 0) = 21;
            this->ELEM(0, 1) = -1;
        }
    }
    else
    {
        if (beamnum == 1)
        {
            this->ELEM = Eigen::MatrixXd::Zero(2, 2);
            this->ELEM(0, 0) = 1;
            this->ELEM(0, 1) = -1;
            this->ELEM(1, 0) = 21;
            this->ELEM(1, 1) = -1;
        }
        else if (beamnum == 2)
        {
            this->ELEM = Eigen::MatrixXd::Zero(2, 2);
            this->ELEM(0, 0) = 1;
            this->ELEM(0, 1) = -1;
            this->ELEM(1, 0) = 21;
            this->ELEM(1, 1) = -1;
        }
    }
}

BeamContact::BeamContact(const int nbeams, std::string str)
{
    this->NBEAMS = nbeams;

    this->epsilon = 1000;

    if (str == "STS")
        this->NEN = 4;
    else if (str == "NTN")
        this->NEN = 2;

    this->GlobalELEM = Eigen::MatrixXd(1, 2);
    for (int i = 0; i < this->GlobalELEM.rows(); i++)
    {
        for (int j = 0; j < this->GlobalELEM.cols(); j++)
            this->GlobalELEM(i, j) = -1;
    }

    this->n = Eigen::VectorXd(3);
    this->n(0) = 0;
    this->n(1) = 0;
    this->n(2) = 1;

    this->start = 1;

    this->precontactiter = 0;
}*/
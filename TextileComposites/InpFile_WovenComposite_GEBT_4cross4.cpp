#pragma once

#include "Variables.h"


//BeamContact::BeamContact(const int nbeams, std::string str)
//{
//    this->NBEAMS = nbeams;
//
//    this->epsilon = 1000000;
//
//    if (str == "STS")
//        this->NEN = 4;
//    else if (str == "NTN")
//        this->NEN = 2;
//
//    this->n = Eigen::VectorXd(3);
//    this->n(0) = 0;
//    this->n(1) = 0;
//    this->n(2) = 1;
//
//    this->start = 1;
//
//    this->precontactiter = 0;
//}
//
//// -1 distributed load
//// -2 point load
////2 * 2 microstructure simulation
////Loading conditions
//Loading::Loading(int beamnum, int contactswitch)
//{
//    if (contactswitch == 1)
//    {
//        if (beamnum == 1)
//        {
//            this->ELEM = Eigen::MatrixXd::Zero(1, 4);
//            this->ELEM(0, 0) = 21;
//            this->ELEM(0, 1) = -2;
//            this->ELEM(0, 2) = 2;
//            this->ELEM(0, 3) = 100;
//        }
//        else if (beamnum == 2)
//        {
//            this->ELEM = Eigen::MatrixXd::Zero(1, 4);
//            this->ELEM(0, 0) = 21;
//            this->ELEM(0, 1) = -2;
//            this->ELEM(0, 2) = 2;
//            this->ELEM(0, 3) = 100;
//        }
//        else if (beamnum == 3)
//        {
//            this->ELEM = Eigen::MatrixXd::Zero(1, 4);
//            this->ELEM(0, 0) = 21;
//            this->ELEM(0, 1) = -2;
//            this->ELEM(0, 2) = 2;
//            this->ELEM(0, 3) = 100;
//        }
//        else if (beamnum == 4)
//        {
//            this->ELEM = Eigen::MatrixXd::Zero(1, 4);
//            this->ELEM(0, 0) = 21;
//            this->ELEM(0, 1) = -2;
//            this->ELEM(0, 2) = 2;
//            this->ELEM(0, 3) = 100;
//        }
//        else if (beamnum == 5)
//        {
//            this->ELEM = Eigen::MatrixXd::Zero(1, 4);
//            this->ELEM(0, 0) = 21;
//            this->ELEM(0, 1) = -2;
//            this->ELEM(0, 2) = 1;
//            this->ELEM(0, 3) = 100;
//        }
//        else if (beamnum == 6)
//        {
//            this->ELEM = Eigen::MatrixXd::Zero(1, 4);
//            this->ELEM(0, 0) = 21;
//            this->ELEM(0, 1) = -2;
//            this->ELEM(0, 2) = 1;
//            this->ELEM(0, 3) = 100;
//        }
//        else if (beamnum == 7)
//        {
//            this->ELEM = Eigen::MatrixXd::Zero(1, 4);
//            this->ELEM(0, 0) = 21;
//            this->ELEM(0, 1) = -2;
//            this->ELEM(0, 2) = 1;
//            this->ELEM(0, 3) = 100;
//        }
//        else if (beamnum == 8)
//        {
//            this->ELEM = Eigen::MatrixXd::Zero(1, 4);
//            this->ELEM(0, 0) = 21;
//            this->ELEM(0, 1) = -2;
//            this->ELEM(0, 2) = 1;
//            this->ELEM(0, 3) = 100;
//        }
//    }
//    else
//    {
//        if (beamnum == 1)
//        {
//            this->ELEM = Eigen::MatrixXd::Zero(0, 0);
//        }
//        else if (beamnum == 2)
//        {
//            this->ELEM = Eigen::MatrixXd::Zero(0, 0);
//        }
//    }
//
//}
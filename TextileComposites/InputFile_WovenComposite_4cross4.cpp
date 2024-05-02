#pragma once
#include "Variables.h"

//Boundary::Boundary(int yn, int nnode, int beamnum, double val)
//{
//    for (int i = 0; i < yn / 2; i++)
//    {
//        this->ELEM(i, 0) = 
//    }
//}

//post contact loading
Loading::Loading(int beamnum)
{
    if (beamnum == 1)
    {
        this->ELEM = Eigen::MatrixXd::Zero(1, 4);
        this->ELEM(0,0) = 21;
        this->ELEM(0,1) = -2;
        this->ELEM(0, 2) = 2;
        this->ELEM(0, 3) = 200;
    }
    else if (beamnum == 2)
    {
        this->ELEM = Eigen::MatrixXd::Zero(1, 4);
        this->ELEM(0, 0) = 21;
        this->ELEM(0, 1) = -2;
        this->ELEM(0, 2) = 2;
        this->ELEM(0, 3) = 200;
    }
    if (beamnum == 3)
    {
        this->ELEM = Eigen::MatrixXd::Zero(1, 4);
        this->ELEM(0, 0) = 21;
        this->ELEM(0, 1) = -2;
        this->ELEM(0, 2) = 2;
        this->ELEM(0, 3) = 200;
    }
    else if (beamnum == 4)
    {
        this->ELEM = Eigen::MatrixXd::Zero(1, 4);
        this->ELEM(0, 0) = 21;
        this->ELEM(0, 1) = -2;
        this->ELEM(0, 2) = 2;
        this->ELEM(0, 3) = 200;
    }
    if (beamnum == 5)
    {
        this->ELEM = Eigen::MatrixXd::Zero(1, 4);
        this->ELEM(0, 0) = 21;
        this->ELEM(0, 1) = -2;
        this->ELEM(0, 2) = 1;
        this->ELEM(0, 3) = 200;
    }
    else if (beamnum == 6)
    {
        this->ELEM = Eigen::MatrixXd::Zero(1, 4);
        this->ELEM(0, 0) = 21;
        this->ELEM(0, 1) = -2;
        this->ELEM(0, 2) = 1;
        this->ELEM(0, 3) = 200;
    }
    if (beamnum == 7)
    {
        this->ELEM = Eigen::MatrixXd::Zero(1, 4);
        this->ELEM(0, 0) = 21;
        this->ELEM(0, 1) = -2;
        this->ELEM(0, 2) = 1;
        this->ELEM(0, 3) = 200;
    }
    else if (beamnum == 8)
    {
        this->ELEM = Eigen::MatrixXd::Zero(1, 4);
        this->ELEM(0, 0) = 21;
        this->ELEM(0, 1) = -2;
        this->ELEM(0, 2) = 1;
        this->ELEM(0, 3) = 200;
    }
}
//
////precontact 
Loading::Loading(int beamnum, int contactswitch)
{
    if (beamnum == 1)
    {
        this->ELEM = Eigen::MatrixXd::Zero(0, 0);
    }
    else if (beamnum == 2)
    {
        this->ELEM = Eigen::MatrixXd::Zero(0, 0);
    }
    if (beamnum == 3)
    {
        this->ELEM = Eigen::MatrixXd::Zero(0, 0);
    }
    else if (beamnum == 4)
    {
        this->ELEM = Eigen::MatrixXd::Zero(0, 0);
    }
    if (beamnum == 5)
    {
        this->ELEM = Eigen::MatrixXd::Zero(0, 0);
    }
    else if (beamnum == 6)
    {
        this->ELEM = Eigen::MatrixXd::Zero(0, 0);
    }
    if (beamnum == 7)
    {
        this->ELEM = Eigen::MatrixXd::Zero(0, 0);
    }
    else if (beamnum == 8)
    {
        this->ELEM = Eigen::MatrixXd::Zero(0, 0);
    }
}
//
////after contact switch is on
Boundary::Boundary(int beamnum)
{
    if (beamnum == 1)
    {
        this->ELEM = Eigen::MatrixXd::Zero(1, 2);
        this->ELEM(0, 0) = 1;
        //type
        this->ELEM(0, 1) = -1;

        //this->ELEM(1, 0) = 11;
        //this->ELEM(1, 1) = -3;
        //this->ELEM(1, 2) = 2;
        //this->ELEM(1, 3) = 0.2;
    }
    else if (beamnum == 2)
    {
        this->ELEM = Eigen::MatrixXd::Zero(1, 2);
        this->ELEM(0, 0) = 1;
        //type
        this->ELEM(0, 1) = -1;

        //this->ELEM(1, 0) = 11;
        //this->ELEM(1, 1) = -3;
        //this->ELEM(1, 2) = 2;
        //this->ELEM(1, 3) = 0.2;
    }
    else if (beamnum == 3)
    {
        this->ELEM = Eigen::MatrixXd::Zero(1, 2);
        this->ELEM(0, 0) = 1;
        //type
        this->ELEM(0, 1) = -1;

        //this->ELEM(1, 0) = 11;
        //this->ELEM(1, 1) = -3;
        //this->ELEM(1, 2) = 2;
        //this->ELEM(1, 3) = 0.2;
    }
    else if (beamnum == 4)
    {
        this->ELEM = Eigen::MatrixXd::Zero(1, 2);
        this->ELEM(0, 0) = 1;
        //type
        this->ELEM(0, 1) = -1;

        //this->ELEM(1, 0) = 11;
        //this->ELEM(1, 1) = -3;
        //this->ELEM(1, 2) = 2;
        //this->ELEM(1, 3) = 0.2;
    }
    else if (beamnum == 5)
    {
        this->ELEM = Eigen::MatrixXd::Zero(1, 2);
        this->ELEM(0, 0) = 1;
        //type
        this->ELEM(0, 1) = -1;
    }
    else if (beamnum == 6)
    {
        this->ELEM = Eigen::MatrixXd::Zero(1, 2);
        this->ELEM(0, 0) = 1;
        //type
        this->ELEM(0, 1) = -1;
    }
    else if (beamnum == 7)
    {
        this->ELEM = Eigen::MatrixXd::Zero(2, 2);
        this->ELEM(0, 0) = 1;
        //type
        this->ELEM(0, 1) = -1;
    }
    else if (beamnum == 8)
    {
        this->ELEM = Eigen::MatrixXd::Zero(1, 2);
        this->ELEM(0, 0) = 1;
        //type
        this->ELEM(0, 1) = -1;
    }
}
//
//
////precontact scenario
Boundary::Boundary(int beamnum, double val)
{
    if (beamnum == 1)
    {
        this->ELEM = Eigen::MatrixXd::Zero(6, 4);
        this->ELEM(0, 0) = 5;
        //type
        this->ELEM(0, 1) = -3;
        this->ELEM(0, 2) = -3;
        this->ELEM(0, 3) = val;

        this->ELEM(1, 0) = 9;
        //type
        this->ELEM(1, 1) = -3;
        this->ELEM(1, 2) = 3;
        this->ELEM(1, 3) = val;


        this->ELEM(2, 0) = 13;
        //type
        this->ELEM(2, 1) = -3;
        this->ELEM(2, 2) = -3;
        this->ELEM(2, 3) = val;


        this->ELEM(3, 0) = 17;
        //type
        this->ELEM(3, 1) = -3;
        this->ELEM(3, 2) = 3;
        this->ELEM(3, 3) = val;

        this->ELEM(4, 0) = 1;
        //type
        this->ELEM(4, 1) = -1;

        this->ELEM(5, 0) = 21;
        //type
        this->ELEM(5, 1) = -1;
    }
    else if (beamnum == 2)
    {
        this->ELEM = Eigen::MatrixXd::Zero(6, 4);

        this->ELEM(0, 0) = 5;
        //type
        this->ELEM(0, 1) = -3;
        this->ELEM(0, 2) = 3;
        this->ELEM(0, 3) = val;

        this->ELEM(1, 0) = 9;
        //type
        this->ELEM(1, 1) = -3;
        this->ELEM(1, 2) = -3;
        this->ELEM(1, 3) = val;


        this->ELEM(2, 0) = 13;
        //type
        this->ELEM(2, 1) = -3;
        this->ELEM(2, 2) = 3;
        this->ELEM(2, 3) = val;


        this->ELEM(3, 0) = 17;
        //type
        this->ELEM(3, 1) = -3;
        this->ELEM(3, 2) = -3;
        this->ELEM(3, 3) = val;

        this->ELEM(4, 0) = 1;
        //type
        this->ELEM(4, 1) = -1;

        this->ELEM(5, 0) = 21;
        //type
        this->ELEM(5, 1) = -1;
    }
    else if (beamnum == 3)
    {
        this->ELEM = Eigen::MatrixXd::Zero(6, 4);
        this->ELEM(0, 0) = 5;
        //type
        this->ELEM(0, 1) = -3;
        this->ELEM(0, 2) = -3;
        this->ELEM(0, 3) = val;

        this->ELEM(1, 0) = 9;
        //type
        this->ELEM(1, 1) = -3;
        this->ELEM(1, 2) = 3;
        this->ELEM(1, 3) = val;


        this->ELEM(2, 0) = 13;
        //type
        this->ELEM(2, 1) = -3;
        this->ELEM(2, 2) = -3;
        this->ELEM(2, 3) = val;


        this->ELEM(3, 0) = 17;
        //type
        this->ELEM(3, 1) = -3;
        this->ELEM(3, 2) = 3;
        this->ELEM(3, 3) = val;

        this->ELEM(4, 0) = 1;
        //type
        this->ELEM(4, 1) = -1;

        this->ELEM(5, 0) = 21;
        //type
        this->ELEM(5, 1) = -1;
    }
    else if (beamnum == 4)
    {
        this->ELEM = Eigen::MatrixXd::Zero(6, 4);
        this->ELEM(0, 0) = 5;
        //type
        this->ELEM(0, 1) = -3;
        this->ELEM(0, 2) = 3;
        this->ELEM(0, 3) = val;

        this->ELEM(1, 0) = 9;
        //type
        this->ELEM(1, 1) = -3;
        this->ELEM(1, 2) = -3;
        this->ELEM(1, 3) = val;


        this->ELEM(2, 0) = 13;
        //type
        this->ELEM(2, 1) = -3;
        this->ELEM(2, 2) = 3;
        this->ELEM(2, 3) = val;


        this->ELEM(3, 0) = 17;
        //type
        this->ELEM(3, 1) = -3;
        this->ELEM(3, 2) = -3;
        this->ELEM(3, 3) = val;

        this->ELEM(4, 0) = 1;
        //type
        this->ELEM(4, 1) = -1;

        this->ELEM(5, 0) = 21;
        //type
        this->ELEM(5, 1) = -1;
    }
    else if (beamnum == 5)
    {
        this->ELEM = Eigen::MatrixXd::Zero(6, 4);
        this->ELEM(0, 0) = 5;
        //type
        this->ELEM(0, 1) = -3;
        this->ELEM(0, 2) = 3;
        this->ELEM(0, 3) = val;

        this->ELEM(1, 0) = 9;
        //type
        this->ELEM(1, 1) = -3;
        this->ELEM(1, 2) = -3;
        this->ELEM(1, 3) = val;


        this->ELEM(2, 0) = 13;
        //type
        this->ELEM(2, 1) = -3;
        this->ELEM(2, 2) = 3;
        this->ELEM(2, 3) = val;


        this->ELEM(3, 0) = 17;
        //type
        this->ELEM(3, 1) = -3;
        this->ELEM(3, 2) = -3;
        this->ELEM(3, 3) = val;

        this->ELEM(4, 0) = 1;
        //type
        this->ELEM(4, 1) = -1;

        this->ELEM(5, 0) = 21;
        //type
        this->ELEM(5, 1) = -1;
    }
    else if (beamnum == 6)
    {
        this->ELEM = Eigen::MatrixXd::Zero(6, 4);
        this->ELEM(0, 0) = 5;
        //type
        this->ELEM(0, 1) = -3;
        this->ELEM(0, 2) = -3;
        this->ELEM(0, 3) = val;

        this->ELEM(1, 0) = 9;
        //type
        this->ELEM(1, 1) = -3;
        this->ELEM(1, 2) = 3;
        this->ELEM(1, 3) = val;


        this->ELEM(2, 0) = 13;
        //type
        this->ELEM(2, 1) = -3;
        this->ELEM(2, 2) = -3;
        this->ELEM(2, 3) = val;


        this->ELEM(3, 0) = 17;
        //type
        this->ELEM(3, 1) = -3;
        this->ELEM(3, 2) = 3;
        this->ELEM(3, 3) = val;

        this->ELEM(4, 0) = 1;
        //type
        this->ELEM(4, 1) = -1;

        this->ELEM(5, 0) = 21;
        //type
        this->ELEM(5, 1) = -1;
    }
    else if (beamnum == 7)
    {
        this->ELEM = Eigen::MatrixXd::Zero(6, 4);
        this->ELEM(0, 0) = 5;
        //type
        this->ELEM(0, 1) = -3;
        this->ELEM(0, 2) = 3;
        this->ELEM(0, 3) = val;

        this->ELEM(1, 0) = 9;
        //type
        this->ELEM(1, 1) = -3;
        this->ELEM(1, 2) = -3;
        this->ELEM(1, 3) = val;


        this->ELEM(2, 0) = 13;
        //type
        this->ELEM(2, 1) = -3;
        this->ELEM(2, 2) = 3;
        this->ELEM(2, 3) = val;


        this->ELEM(3, 0) = 17;
        //type
        this->ELEM(3, 1) = -3;
        this->ELEM(3, 2) = -3;
        this->ELEM(3, 3) = val;

        this->ELEM(4, 0) = 1;
        //type
        this->ELEM(4, 1) = -1;

        this->ELEM(5, 0) = 21;
        //type
        this->ELEM(5, 1) = -1;
    }
    else if (beamnum == 8)
    {
        this->ELEM = Eigen::MatrixXd::Zero(6, 4);
        this->ELEM(0, 0) = 5;
        //type
        this->ELEM(0, 1) = -3;
        this->ELEM(0, 2) = -3;
        this->ELEM(0, 3) = val;

        this->ELEM(1, 0) = 9;
        //type
        this->ELEM(1, 1) = -3;
        this->ELEM(1, 2) = 3;
        this->ELEM(1, 3) = val;


        this->ELEM(2, 0) = 13;
        //type
        this->ELEM(2, 1) = -3;
        this->ELEM(2, 2) = -3;
        this->ELEM(2, 3) = val;


        this->ELEM(3, 0) = 17;
        //type
        this->ELEM(3, 1) = -3;
        this->ELEM(3, 2) = 3;
        this->ELEM(3, 3) = val;

        this->ELEM(4, 0) = 1;
        //type
        this->ELEM(4, 1) = -1;

        this->ELEM(5, 0) = 21;
        //type
        this->ELEM(5, 1) = -1;
    }
}

BeamContact::BeamContact(const int nbeams, std::string str)
{
    this->NBEAMS = nbeams;

    this->epsilon = 10000;

    if (str == "STS")
        this->NEN = 4;
    else if (str == "NTN")
        this->NEN = 2;

    this->n = Eigen::VectorXd(3);
    this->n(0) = 0;
    this->n(1) = 0;
    this->n(2) = 1;

    this->start = 0;

    this->precontactiter = 10;
}
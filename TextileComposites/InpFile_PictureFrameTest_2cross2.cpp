#pragma once

#include "Variables.h"

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

////precontact
Boundary::Boundary(int beamnum, double val)
{
    if (beamnum == 1)
    {
        this->ELEM = Eigen::MatrixXd::Zero(4, 4);
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

        this->ELEM(2, 0) = 1;
        //type
        this->ELEM(2, 1) = -1;

        this->ELEM(3, 0) = 13;
        //type
        this->ELEM(3, 1) = -1;
    }
    else if (beamnum == 2)
    {
        this->ELEM = Eigen::MatrixXd::Zero(4, 4);

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

        this->ELEM(2, 0) = 1;
        //type
        this->ELEM(2, 1) = -1;

        this->ELEM(3, 0) = 13;
        //type
        this->ELEM(3, 1) = -1;


    }
    else if (beamnum == 3)
    {
        this->ELEM = Eigen::MatrixXd::Zero(4, 4);
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

        this->ELEM(2, 0) = 1;
        //type
        this->ELEM(2, 1) = -1;

        this->ELEM(3, 0) = 13;
        //type
        this->ELEM(3, 1) = -1;
    }
    else if (beamnum == 4)
    {
        this->ELEM = Eigen::MatrixXd::Zero(4, 4);
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

        this->ELEM(2, 0) = 1;
        //type
        this->ELEM(2, 1) = -1;

        this->ELEM(3, 0) = 13;
        //type
        this->ELEM(3, 1) = -1;
    }
}

//after contact is established
Boundary::Boundary(int beamnum)
{
    if (beamnum == 1)
    {
        this->ELEM = Eigen::MatrixXd::Zero(0, 0);
        //this->ELEM(0, 0) = 1;
        //type
        //this->ELEM(0, 1) = -1;

        //this->ELEM(1, 0) = 13;
        //this->ELEM(1, 1) = -3;
        //this->ELEM(1, 2) = 2;
        //this->ELEM(1, 3) = 0.15;
    }
    else if (beamnum == 2)
    {
        this->ELEM = Eigen::MatrixXd::Zero(0, 0);
        //this->ELEM(0, 0) = 1;
        //type
        //this->ELEM(0, 1) = -1;

        //this->ELEM(1, 0) = 13;
        //this->ELEM(1, 1) = -3;
        //this->ELEM(1, 2) = 2;
        //this->ELEM(1, 3) = 0.15;
    }
    else if (beamnum == 3)
    {
        this->ELEM = Eigen::MatrixXd::Zero(0, 0);
        //this->ELEM(0, 0) = 1;
        //type
        //this->ELEM(0, 1) = -1;

        //this->ELEM(1, 0) = 13;
        ////type
        //this->ELEM(1, 1) = -1;
        //this->ELEM(1, 0) = 7;
        //this->ELEM(1, 1) = -3;
        //this->ELEM(1, 2) = 2;
        //this->ELEM(1, 3) = 0.2;
    }
    else if (beamnum == 4)
    {
        this->ELEM = Eigen::MatrixXd::Zero(0, 0);
        //this->ELEM(0, 0) = 1;
        //type
        //this->ELEM(0, 1) = -1;

        //this->ELEM(1, 0) = 13;
        ////type
        //this->ELEM(1, 1) = -1;
        //this->ELEM(1, 0) = 7;
        //this->ELEM(1, 1) = -3;
        //this->ELEM(1, 2) = 2;
        //this->ELEM(1, 3) = 0.2;
    }
}
//
////postcontact
Loading::Loading(int beamnum)
{
    if (beamnum == 1)
    {
        this->ELEM = Eigen::MatrixXd::Zero(4, 4);
        this->ELEM(0, 0) = 13;
        this->ELEM(0, 1) = -2;
        this->ELEM(0, 2) = 2;
        this->ELEM(0, 3) = 200;

        this->ELEM(1, 0) = 13;
        this->ELEM(1, 1) = -2;
        this->ELEM(1, 2) = 1;
        this->ELEM(1, 3) = 200;

        this->ELEM(2, 0) = 1;
        this->ELEM(2, 1) = -2;
        this->ELEM(2, 2) = -2;
        this->ELEM(2, 3) = 200;

        this->ELEM(3, 0) = 1;
        this->ELEM(3, 1) = -2;
        this->ELEM(3, 2) = -1;
        this->ELEM(3, 3) = 200;
    }
    else if (beamnum == 2)
    {
        this->ELEM = Eigen::MatrixXd::Zero(4, 4);
        this->ELEM(0, 0) = 13;
        this->ELEM(0, 1) = -2;
        this->ELEM(0, 2) = 2;
        this->ELEM(0, 3) = 200;

        this->ELEM(1, 0) = 13;
        this->ELEM(1, 1) = -2;
        this->ELEM(1, 2) = 1;
        this->ELEM(1, 3) = 200;

        this->ELEM(2, 0) = 1;
        this->ELEM(2, 1) = -2;
        this->ELEM(2, 2) = -2;
        this->ELEM(2, 3) = 200;

        this->ELEM(3, 0) = 1;
        this->ELEM(3, 1) = -2;
        this->ELEM(3, 2) = -1;
        this->ELEM(3, 3) = 200;
    }
    if (beamnum == 3)
    {
        this->ELEM = Eigen::MatrixXd::Zero(4, 4);
        this->ELEM(0, 0) = 13;
        this->ELEM(0, 1) = -2;
        this->ELEM(0, 2) = 2;
        this->ELEM(0, 3) = 200;

        this->ELEM(1, 0) = 13;
        this->ELEM(1, 1) = -2;
        this->ELEM(1, 2) = 1;
        this->ELEM(1, 3) = 200;

        this->ELEM(2, 0) = 1;
        this->ELEM(2, 1) = -2;
        this->ELEM(2, 2) = -2;
        this->ELEM(2, 3) = 200;

        this->ELEM(3, 0) = 1;
        this->ELEM(3, 1) = -2;
        this->ELEM(3, 2) = -1;
        this->ELEM(3, 3) = 200;
    }
    else if (beamnum == 4)
    {
        this->ELEM = Eigen::MatrixXd::Zero(4, 4);
        this->ELEM(0, 0) = 13;
        this->ELEM(0, 1) = -2;
        this->ELEM(0, 2) = 2;
        this->ELEM(0, 3) = 200;

        this->ELEM(1, 0) = 13;
        this->ELEM(1, 1) = -2;
        this->ELEM(1, 2) = 1;
        this->ELEM(1, 3) = 200;

        this->ELEM(2, 0) = 1;
        this->ELEM(2, 1) = -2;
        this->ELEM(2, 2) = -2;
        this->ELEM(2, 3) = 200;

        this->ELEM(3, 0) = 1;
        this->ELEM(3, 1) = -2;
        this->ELEM(3, 2) = -1;
        this->ELEM(3, 3) = 200;
    }
}
//
////precontact
Loading::Loading(int beamnum, int choice)
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
}

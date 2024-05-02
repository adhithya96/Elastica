#pragma once

#include "Variables.h"


//BeamContact::BeamContact(const int nbeams, std::string str)
//{
//    this->NBEAMS = nbeams;
//
//    this->epsilon = 100;
//
//    if (str == "STS")
//        this->NEN = 4;
//    else if (str == "NTN")
//        this->NEN = 2;
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




// -1 distributed load
// -2 point load
//2 * 2 microstructure simulation
//Loading conditions
//Loading::Loading(int beamnum, int contactswitch)
//{
//    if (contactswitch == 1)
//    {
//        if (beamnum == 1)
//        {
//            this->ELEM = Eigen::MatrixXd::Zero(1, 4);
//            //node
//            this->ELEM(0, 0) = 3;
//            //type
//            this->ELEM(0, 1) = -2;
//            //direction
//            this->ELEM(0, 2) = 2;
//            //value
//            this->ELEM(0, 3) = 400;
//        }
//        else if (beamnum == 2)
//        {
//            this->ELEM = Eigen::MatrixXd::Zero(1, 4);
//            //node
//            this->ELEM(0, 0) = 3;
//            //type
//            this->ELEM(0, 1) = -2;
//            //direction
//            this->ELEM(0, 2) = -2;
//            //value
//            this->ELEM(0, 3) = 400;
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

VAMBeamElement::VAMBeamElement(int beamnum)
{
    if (beamnum == 1)
    {
        this->NNODE = 6;
        this->NELEM = 5;
        this->NDOF = 12;
        this->NLS = 11;
        this->NDM = 3;

        this->NODE = Eigen::MatrixXd::Zero(this->NNODE, 3);
        this->ELEM = Eigen::MatrixXd::Zero(this->NELEM, 3);

        //    Nodal Information
        this->NODE(0, 0) = 0;
        this->NODE(0, 1) = 0;
        this->NODE(0, 2) = 0;
        this->NODE(1, 0) = 20;
        this->NODE(1, 1) = 0;
        this->NODE(1, 2) = 0;
        this->NODE(2, 0) = 40;
        this->NODE(2, 1) = 0;
        this->NODE(2, 2) = 0;
        this->NODE(3, 0) = 60;
        this->NODE(3, 1) = 0;
        this->NODE(3, 2) = 0;
        this->NODE(4, 0) = 80;
        this->NODE(4, 1) = 0;
        this->NODE(4, 2) = 0;
        this->NODE(5, 0) = 100;
        this->NODE(5, 1) = 0;
        this->NODE(5, 2) = 0;

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

        this->D = 8;
    }
    else if (beamnum == 2)
    {
        this->NNODE = 6;
        this->NELEM = 5;
        this->NDOF = 12;
        this->NLS = 11;
        this->NDM = 3;

        this->NODE = Eigen::MatrixXd::Zero(this->NNODE, 3);
        this->ELEM = Eigen::MatrixXd::Zero(this->NELEM, 3);

        //    Nodal Information
        this->NODE(0, 0) = 0;
        this->NODE(0, 1) = 25;
        this->NODE(0, 2) = 0;
        this->NODE(1, 0) = 20;
        this->NODE(1, 1) = 25;
        this->NODE(1, 2) = 0;
        this->NODE(2, 0) = 40;
        this->NODE(2, 1) = 25;
        this->NODE(2, 2) = 0;
        this->NODE(3, 0) = 60;
        this->NODE(3, 1) = 25;
        this->NODE(3, 2) = 0;
        this->NODE(4, 0) = 80;
        this->NODE(4, 1) = 25;
        this->NODE(4, 2) = 0;
        this->NODE(5, 0) = 100;
        this->NODE(5, 1) = 25;
        this->NODE(5, 2) = 0;

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

        this->D = 8;
    }
    else
    {
        std::cout << "Check the no. of beams" << std::endl;
    }
}

//Inp file
//Single beam case
//Isotropic validation with elastica
//VAMBeamElement::VAMBeamElement(int beamnum)
//{
//    if (beamnum == 1)
//    {
//        this->NNODE = 17;
//        this->NELEM = 16;
//        this->NDOF = 12;
//        this->NLS = 11;
//
//        this->NODE = Eigen::MatrixXd::Zero(this->NNODE, 3);
//        this->ELEM = Eigen::MatrixXd::Zero(this->NELEM, 3);
//
//        //    Nodal Information
//        this->NODE(0, 0) = 0;
//        this->NODE(0, 1) = 0;
//        this->NODE(0, 2) = 0;
//        this->NODE(1, 0) = 6.25;
//        this->NODE(1, 1) = 0;
//        this->NODE(1, 2) = 0;
//        this->NODE(2, 0) = 12.5;
//        this->NODE(2, 1) = 0;
//        this->NODE(2, 2) = 0;
//        this->NODE(3, 0) = 18.75;
//        this->NODE(3, 1) = 0;
//        this->NODE(3, 2) = 0;
//        this->NODE(4, 0) = 25.0;
//        this->NODE(4, 1) = 0;
//        this->NODE(4, 2) = 0;
//        this->NODE(5, 0) = 31.25;
//        this->NODE(5, 1) = 0;
//        this->NODE(5, 2) = 0;
//        this->NODE(6, 0) = 37.5;
//        this->NODE(6, 1) = 0;
//        this->NODE(6, 2) = 0;
//        this->NODE(7, 0) = 43.75;
//        this->NODE(7, 1) = 0;
//        this->NODE(7, 2) = 0;
//        this->NODE(8, 0) = 50.0;
//        this->NODE(8, 1) = 0;
//        this->NODE(8, 2) = 0;
//        this->NODE(9, 0) = 56.25;
//        this->NODE(9, 1) = 0;
//        this->NODE(9, 2) = 0;
//        this->NODE(10, 0) = 62.5;
//        this->NODE(10, 1) = 0;
//        this->NODE(10, 2) = 0;
//        this->NODE(11, 0) = 68.75;
//        this->NODE(11, 1) = 0;
//        this->NODE(11, 2) = 0;
//        this->NODE(12, 0) = 75.0;
//        this->NODE(12, 1) = 0;
//        this->NODE(12, 2) = 0;
//        this->NODE(13, 0) = 81.25;
//        this->NODE(13, 1) = 0;
//        this->NODE(13, 2) = 0;
//        this->NODE(14, 0) = 87.5;
//        this->NODE(14, 1) = 0;
//        this->NODE(14, 2) = 0;
//        this->NODE(15, 0) = 93.75;
//        this->NODE(15, 1) = 0;
//        this->NODE(15, 2) = 0;
//        this->NODE(16, 0) = 100;
//        this->NODE(16, 1) = 0;
//        this->NODE(16, 2) = 0;
//
//        //Connectivity information
//        this->ELEM(0, 0) = 1;
//        this->ELEM(0, 1) = 1;
//        this->ELEM(0, 2) = 2;
//        this->ELEM(1, 0) = 1;
//        this->ELEM(1, 1) = 2;
//        this->ELEM(1, 2) = 3;
//        this->ELEM(2, 0) = 1;
//        this->ELEM(2, 1) = 3;
//        this->ELEM(2, 2) = 4;
//        this->ELEM(3, 0) = 1;
//        this->ELEM(3, 1) = 4;
//        this->ELEM(3, 2) = 5;
//        this->ELEM(4, 0) = 1;
//        this->ELEM(4, 1) = 5;
//        this->ELEM(4, 2) = 6;
//        this->ELEM(5, 0) = 1;
//        this->ELEM(5, 1) = 6;
//        this->ELEM(5, 2) = 7;
//        this->ELEM(6, 0) = 1;
//        this->ELEM(6, 1) = 7;
//        this->ELEM(6, 2) = 8;
//        this->ELEM(7, 0) = 1;
//        this->ELEM(7, 1) = 8;
//        this->ELEM(7, 2) = 9;
//        this->ELEM(8, 0) = 1;
//        this->ELEM(8, 1) = 9;
//        this->ELEM(8, 2) = 10;
//        this->ELEM(9, 0) = 1;
//        this->ELEM(9, 1) = 10;
//        this->ELEM(9, 2) = 11;
//        this->ELEM(10, 0) = 1;
//        this->ELEM(10, 1) = 11;
//        this->ELEM(10, 2) = 12;
//        this->ELEM(11, 0) = 1;
//        this->ELEM(11, 1) = 12;
//        this->ELEM(11, 2) = 13;
//        this->ELEM(12, 0) = 1;
//        this->ELEM(12, 1) = 13;
//        this->ELEM(12, 2) = 14;
//        this->ELEM(13, 0) = 1;
//        this->ELEM(13, 1) = 14;
//        this->ELEM(13, 2) = 15;
//        this->ELEM(14, 0) = 1;
//        this->ELEM(14, 1) = 15;
//        this->ELEM(14, 2) = 16;
//        this->ELEM(15, 0) = 1;
//        this->ELEM(15, 1) = 16;
//        this->ELEM(15, 2) = 17;
//
//    }
//    else if(beamnum == 2)
//    {
//        this->NNODE = 17;
//        this->NELEM = 16;
//        this->NDOF = 12;
//        this->NLS = 11;
//
//        this->NODE = Eigen::MatrixXd::Zero(this->NNODE, 3);
//        this->ELEM = Eigen::MatrixXd::Zero(this->NELEM, 3);
//
//        //    Nodal Information
//        this->NODE(0, 0) = 0;
//        this->NODE(0, 1) = 10;
//        this->NODE(0, 2) = 0;
//        this->NODE(1, 0) = 6.25;
//        this->NODE(1, 1) = 10;
//        this->NODE(1, 2) = 0;
//        this->NODE(2, 0) = 12.5;
//        this->NODE(2, 1) = 10;
//        this->NODE(2, 2) = 0;
//        this->NODE(3, 0) = 18.75;
//        this->NODE(3, 1) = 10;
//        this->NODE(3, 2) = 0;
//        this->NODE(4, 0) = 25.0;
//        this->NODE(4, 1) = 10;
//        this->NODE(4, 2) = 0;
//        this->NODE(5, 0) = 31.25;
//        this->NODE(5, 1) = 10;
//        this->NODE(5, 2) = 0;
//        this->NODE(6, 0) = 37.5;
//        this->NODE(6, 1) = 10;
//        this->NODE(6, 2) = 0;
//        this->NODE(7, 0) = 43.75;
//        this->NODE(7, 1) = 10;
//        this->NODE(7, 2) = 0;
//        this->NODE(8, 0) = 50.0;
//        this->NODE(8, 1) = 10;
//        this->NODE(8, 2) = 0;
//        this->NODE(9, 0) = 56.25;
//        this->NODE(9, 1) = 10;
//        this->NODE(9, 2) = 0;
//        this->NODE(10, 0) = 62.5;
//        this->NODE(10, 1) = 10;
//        this->NODE(10, 2) = 0;
//        this->NODE(11, 0) = 68.75;
//        this->NODE(11, 1) = 10;
//        this->NODE(11, 2) = 0;
//        this->NODE(12, 0) = 75.0;
//        this->NODE(12, 1) = 10;
//        this->NODE(12, 2) = 0;
//        this->NODE(13, 0) = 81.25;
//        this->NODE(13, 1) = 10;
//        this->NODE(13, 2) = 0;
//        this->NODE(14, 0) = 87.5;
//        this->NODE(14, 1) = 10;
//        this->NODE(14, 2) = 0;
//        this->NODE(15, 0) = 93.75;
//        this->NODE(15, 1) = 10;
//        this->NODE(15, 2) = 0;
//        this->NODE(16, 0) = 100;
//        this->NODE(16, 1) = 10;
//        this->NODE(16, 2) = 0;
//
//        //Connectivity information
//        this->ELEM(0, 0) = 1;
//        this->ELEM(0, 1) = 1;
//        this->ELEM(0, 2) = 2;
//        this->ELEM(1, 0) = 1;
//        this->ELEM(1, 1) = 2;
//        this->ELEM(1, 2) = 3;
//        this->ELEM(2, 0) = 1;
//        this->ELEM(2, 1) = 3;
//        this->ELEM(2, 2) = 4;
//        this->ELEM(3, 0) = 1;
//        this->ELEM(3, 1) = 4;
//        this->ELEM(3, 2) = 5;
//        this->ELEM(4, 0) = 1;
//        this->ELEM(4, 1) = 5;
//        this->ELEM(4, 2) = 6;
//        this->ELEM(5, 0) = 1;
//        this->ELEM(5, 1) = 6;
//        this->ELEM(5, 2) = 7;
//        this->ELEM(6, 0) = 1;
//        this->ELEM(6, 1) = 7;
//        this->ELEM(6, 2) = 8;
//        this->ELEM(7, 0) = 1;
//        this->ELEM(7, 1) = 8;
//        this->ELEM(7, 2) = 9;
//        this->ELEM(8, 0) = 1;
//        this->ELEM(8, 1) = 9;
//        this->ELEM(8, 2) = 10;
//        this->ELEM(9, 0) = 1;
//        this->ELEM(9, 1) = 10;
//        this->ELEM(9, 2) = 11;
//        this->ELEM(10, 0) = 1;
//        this->ELEM(10, 1) = 11;
//        this->ELEM(10, 2) = 12;
//        this->ELEM(11, 0) = 1;
//        this->ELEM(11, 1) = 12;
//        this->ELEM(11, 2) = 13;
//        this->ELEM(12, 0) = 1;
//        this->ELEM(12, 1) = 13;
//        this->ELEM(12, 2) = 14;
//        this->ELEM(13, 0) = 1;
//        this->ELEM(13, 1) = 14;
//        this->ELEM(13, 2) = 15;
//        this->ELEM(14, 0) = 1;
//        this->ELEM(14, 1) = 15;
//        this->ELEM(14, 2) = 16;
//        this->ELEM(15, 0) = 1;
//        this->ELEM(15, 1) = 16;
//        this->ELEM(15, 2) = 17;
//
//    }
//    
//
//}
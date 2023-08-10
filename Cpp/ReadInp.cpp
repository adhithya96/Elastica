#ifndef VARIABLES_H
#define VARIABLES_H
#include "Variables.h"

LinearBarElement ReadLBEFile()
{
    struct LinearBarElement LBE;

    LBE.NBODIES = 1;
    LBE.NNODE = 5;
    LBE.NELEM = 4;
    LBE.NMAT = 2;
    LBE.NDOF = 1;

    LBE.NODE = Eigen::MatrixXd::Zero(LBE.NNODE, LBE.NDOF);
    LBE.ELEM = Eigen::MatrixXd::Zero(LBE.NELEM, 3);
    LBE.MAT = Eigen::MatrixXd::Zero(LBE.NMAT, 2);

    LBE.MAT(0, 0) = 30e6;
    LBE.MAT(0, 1) = 0.33;
    LBE.MAT(1, 0) = 30e6;
    LBE.MAT(1, 1) = 0.33;

    LBE.CS.choice = "RECT";
    LBE.CS.Rect.height = 1;
    LBE.CS.Rect.height = 1;

    LBE.NODE(0, 0) = 0;
    LBE.NODE(1, 0) = 0.25;
    LBE.NODE(2, 0) = 0.5;
    LBE.NODE(3, 0) = 0.75;
    LBE.NODE(4, 0) = 1;

    LBE.ELEM(0, 0) = 1;
    LBE.ELEM(0, 1) = 1;
    LBE.ELEM(0, 2) = 2;
    LBE.ELEM(1, 0) = 1;
    LBE.ELEM(1, 1) = 2;
    LBE.ELEM(1, 2) = 3;
    LBE.ELEM(2, 0) = 1;
    LBE.ELEM(2, 1) = 3;
    LBE.ELEM(2, 2) = 4;
    LBE.ELEM(3, 0) = 2;
    LBE.ELEM(3, 1) = 4;
    LBE.ELEM(3, 2) = 5;

    LBE.CNODE = Eigen::MatrixXd::Zero(1, 2);
    LBE.CNODE(0, 0) = 1;
    LBE.CNODE(0, 1) = 0;


    LBE.LOAD = Eigen::MatrixXd::Zero(1, 2);
    LBE.LOAD(0, 0) = 5;
    LBE.LOAD(0, 1) = 1;

    return LBE;
}

NonLinearBarElement ReadNLBEFile()
{
    struct NonLinearBarElement NLBE;

    NLBE.NBODIES = 1;
    NLBE.NNODE = 3;
    NLBE.NELEM = 2;
    NLBE.NMAT = 1;
    NLBE.NDOF = 1;

    NLBE.NODE = Eigen::MatrixXd::Zero(NLBE.NNODE, 2);
    NLBE.ELEM = Eigen::MatrixXd::Zero(NLBE.NELEM, 2);

    NLBE.CS.choice = "RECT";
    NLBE.CS.Rect.height = 1;
    NLBE.CS.Rect.width = 1;

    NLBE.NODE(0, 0) = 0;
    NLBE.NODE(0, 1) = 0;
    NLBE.NODE(1, 0) = 0.5;
    NLBE.NODE(1, 1) = 0;
    NLBE.NODE(2, 0) = 1;
    NLBE.NODE(2, 1) = 0;

    NLBE.ELEM(0, 0) = 1;
    NLBE.ELEM(0, 1) = 2;
    NLBE.ELEM(1, 0) = 2;
    NLBE.ELEM(1, 1) = 3;

    NLBE.CNODE = Eigen::MatrixXd::Zero(1, 2);
    NLBE.CNODE(0, 0) = 3;
    NLBE.CNODE(0, 1) = sqrt(2);

    return NLBE;
}

NLEBBE2D::NLEBBE2D()
{
    NLEBBE2D::NNODE = 9;
    NLEBBE2D::NELEM = 8;
    NLEBBE2D::NMAT = 2;
    NLEBBE2D::NDOF = 3;
    NLEBBE2D::NLS = 12;

    NLEBBE2D::NODE = Eigen::MatrixXd::Zero(NLEBBE2D::NNODE, 3);
    NLEBBE2D::ELEM = Eigen::MatrixXd::Zero(NLEBBE2D::NELEM, 3);

    NLEBBE2D::E = 30e6;
    NLEBBE2D::nu = 0.33;

    NLEBBE2D::h = 1;
    NLEBBE2D::b = 1;

    NLEBBE2D::NODE(0, 0) = 0;
    NLEBBE2D::NODE(0, 1) = 0;
    NLEBBE2D::NODE(0, 2) = 0;
    NLEBBE2D::NODE(1, 0) = 12.5;
    NLEBBE2D::NODE(1, 1) = 0;
    NLEBBE2D::NODE(1, 2) = 0;
    NLEBBE2D::NODE(2, 0) = 25;
    NLEBBE2D::NODE(2, 1) = 0;
    NLEBBE2D::NODE(2, 2) = 0;
    NLEBBE2D::NODE(3, 0) = 37.5;
    NLEBBE2D::NODE(3, 1) = 0;
    NLEBBE2D::NODE(3, 2) = 0;
    NLEBBE2D::NODE(4, 0) = 50;
    NLEBBE2D::NODE(4, 1) = 0;
    NLEBBE2D::NODE(4, 2) = 0;
    NLEBBE2D::NODE(5, 0) = 62.5;
    NLEBBE2D::NODE(5, 1) = 0;
    NLEBBE2D::NODE(5, 2) = 0;
    NLEBBE2D::NODE(6, 0) = 75;
    NLEBBE2D::NODE(6, 1) = 0;
    NLEBBE2D::NODE(6, 2) = 0;
    NLEBBE2D::NODE(7, 0) = 87.5;
    NLEBBE2D::NODE(7, 1) = 0;
    NLEBBE2D::NODE(7, 2) = 0;
    NLEBBE2D::NODE(8, 0) = 100;
    NLEBBE2D::NODE(8, 1) = 0;
    NLEBBE2D::NODE(8, 2) = 0;

    NLEBBE2D::ELEM(0, 0) = 1;
    NLEBBE2D::ELEM(0, 1) = 1;
    NLEBBE2D::ELEM(0, 2) = 2;
    NLEBBE2D::ELEM(1, 0) = 1;
    NLEBBE2D::ELEM(1, 1) = 2;
    NLEBBE2D::ELEM(1, 2) = 3;
    NLEBBE2D::ELEM(2, 0) = 1;
    NLEBBE2D::ELEM(2, 1) = 3;
    NLEBBE2D::ELEM(2, 2) = 4;
    NLEBBE2D::ELEM(3, 0) = 1;
    NLEBBE2D::ELEM(3, 1) = 4;
    NLEBBE2D::ELEM(3, 2) = 5;
    NLEBBE2D::ELEM(4, 0) = 1;
    NLEBBE2D::ELEM(4, 1) = 5;
    NLEBBE2D::ELEM(4, 2) = 6;
    NLEBBE2D::ELEM(5, 0) = 1;
    NLEBBE2D::ELEM(5, 1) = 6;
    NLEBBE2D::ELEM(5, 2) = 7;
    NLEBBE2D::ELEM(6, 0) = 1;
    NLEBBE2D::ELEM(6, 1) = 7;
    NLEBBE2D::ELEM(6, 2) = 8;
    NLEBBE2D::ELEM(7, 0) = 1;
    NLEBBE2D::ELEM(7, 1) = 8;
    NLEBBE2D::ELEM(7, 2) = 9;

    NLEBBE2D::CNODE = Eigen::MatrixXd::Zero(4, 2);
    NLEBBE2D::CNODE(0, 0) = 1;
    NLEBBE2D::CNODE(0, 1) = 2;
    NLEBBE2D::CNODE(1, 0) = 5;
    NLEBBE2D::CNODE(1, 1) = 1;
    NLEBBE2D::CNODE(2, 0) = 5;
    NLEBBE2D::CNODE(2, 1) = 3;
    NLEBBE2D::CNODE(3, 0) = 9;
    NLEBBE2D::CNODE(3, 1) = 2;

    NLEBBE2D::vf = 0;
    NLEBBE2D::af = 0;
}

//Validation case for GEBT (12 dofs per node) axial load 0/0 ply.
VAMBeamElement ReadVAMBEFile()
{
    struct VAMBeamElement M;
    M.NBODIES = 1;
    M.NNODE = 3;
    M.NELEM = 2;
    M.NMAT = 1;
    M.NDOF = 12;
    M.NLS = 10;

    M.NODE = Eigen::MatrixXd::Zero(M.NNODE, 3);
    M.ELEM = Eigen::MatrixXd::Zero(M.NELEM, 3);

    //    Nodal Information
    M.NODE(0, 0) = 0;
    M.NODE(0, 1) = 0;
    M.NODE(0, 2) = 0;
    M.NODE(1, 0) = 0.127;
    M.NODE(1, 1) = 0;
    M.NODE(1, 2) = 0;
    M.NODE(2, 0) = 0.254;
    M.NODE(2, 1) = 0;
    M.NODE(2, 2) = 0;

    //Connectivity information
    M.ELEM(0, 0) = 1;
    M.ELEM(0, 1) = 1;
    M.ELEM(0, 2) = 2;
    M.ELEM(1, 0) = 1;
    M.ELEM(1, 1) = 2;
    M.ELEM(1, 2) = 3;

    M.inittwist = Eigen::VectorXd::Zero(3);
    M.inittwist(0) = 0;
    M.inittwist(1) = 0;
    M.inittwist(2) = 0;

    M.B.a1 = Eigen::VectorXd::Zero(3);
    M.B.b1 = Eigen::VectorXd::Zero(3);
    M.B.aN = Eigen::VectorXd::Zero(3);
    M.B.bN = Eigen::VectorXd::Zero(3);

    M.B.x1 = 1;
    M.B.y1 = 2;
    M.B.xN = 3;
    M.B.yN = 4;

    M.B.aN(0) = 600;

    return M;
}

//Validation case for Euler Bernouli Beam element (6 dofs per node)
NLEBBE3D::NLEBBE3D(int choice)
{
    if (choice == 1)
    {
        NLEBBE3D::NNODE = 17;
        NLEBBE3D::NELEM = 8;
        NLEBBE3D::NDOF = 6;
        NLEBBE3D::NLS = 12;
        NLEBBE3D::NEN = 3;

        NLEBBE3D::NODE = Eigen::MatrixXd::Zero(NLEBBE3D::NNODE, 3);
        NLEBBE3D::ELEM = Eigen::MatrixXd::Zero(NLEBBE3D::NELEM, 4);

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
    }
    else
        std::cout << "Wrong input" << std::endl;
}

/*NonLinearEulerBernouliBeamElement3D ReadEBBE3DElement()
{
    struct NonLinearEulerBernouliBeamElement3D M;
    M.NNODE = 5;
    M.NELEM = 2;
    M.NDOF = 6;
    M.NLS = 24;
    M.NEN = 3;

    M.NODE = Eigen::MatrixXd::Zero(M.NNODE, 3);
    M.ELEM = Eigen::MatrixXd::Zero(M.NELEM, 4);

    //Nodal Information
    M.NODE(0, 0) = 0;
    M.NODE(0, 1) = 0;
    M.NODE(0, 2) = 0;
    M.NODE(1, 0) = 2.5;
    M.NODE(1, 1) = 0;
    M.NODE(1, 2) = 0;
    M.NODE(2, 0) = 5;
    M.NODE(2, 1) = 0;
    M.NODE(2, 2) = 0;
    M.NODE(3, 0) = 7.5;
    M.NODE(3, 1) = 0;
    M.NODE(3, 2) = 0;
    M.NODE(4, 0) = 10;
    M.NODE(4, 1) = 0;
    M.NODE(4, 2) = 0;

    //Element connectivity
    M.ELEM(0, 0) = 1;
    M.ELEM(0, 1) = 1;
    M.ELEM(0, 2) = 2;
    M.ELEM(0, 3) = 3;

    M.ELEM(1, 0) = 1;
    M.ELEM(1, 1) = 3;
    M.ELEM(1, 2) = 4;
    M.ELEM(1, 3) = 5;

    M.E = 1e7;
    M.nu = 0.0001;
    M.Bp = 1;
    M.Hp = 1;
    M.Zx = 0;
    M.Zy = 0;
    M.Zz = 1;

    return M;
}*/

//Validation case for node-to-node contact with 3-D EB beams
/*NLEBBE3D::NLEBBE3D(int choice)
{
    //Master
    if (choice == 1)
    {
        NLEBBE3D::NNODE = 11;
        NLEBBE3D::NELEM = 5;
        NLEBBE3D::NDOF = 6;
        NLEBBE3D::NLS = 24;
        NLEBBE3D::NEN = 3;

        NLEBBE3D::NODE = Eigen::MatrixXd::Zero(NLEBBE3D::NNODE, 3);
        NLEBBE3D::ELEM = Eigen::MatrixXd::Zero(NLEBBE3D::NELEM, 4);

        //Nodal Information
        NLEBBE3D::NODE(0, 0) = 0;
        NLEBBE3D::NODE(0, 1) = 0;
        NLEBBE3D::NODE(0, 2) = 0;
        NLEBBE3D::NODE(1, 0) = 1;
        NLEBBE3D::NODE(1, 1) = 0;
        NLEBBE3D::NODE(1, 2) = 0;
        NLEBBE3D::NODE(2, 0) = 2;
        NLEBBE3D::NODE(2, 1) = 0;
        NLEBBE3D::NODE(2, 2) = 0;
        NLEBBE3D::NODE(3, 0) = 3;
        NLEBBE3D::NODE(3, 1) = 0;
        NLEBBE3D::NODE(3, 2) = 0;
        NLEBBE3D::NODE(4, 0) = 4;
        NLEBBE3D::NODE(4, 1) = 0;
        NLEBBE3D::NODE(4, 2) = 0;
        NLEBBE3D::NODE(5, 0) = 5;
        NLEBBE3D::NODE(5, 1) = 0;
        NLEBBE3D::NODE(5, 2) = 0;
        NLEBBE3D::NODE(6, 0) = 6;
        NLEBBE3D::NODE(6, 1) = 0;
        NLEBBE3D::NODE(6, 2) = 0;
        NLEBBE3D::NODE(7, 0) = 7;
        NLEBBE3D::NODE(7, 1) = 0;
        NLEBBE3D::NODE(7, 2) = 0;
        NLEBBE3D::NODE(8, 0) = 8;
        NLEBBE3D::NODE(8, 1) = 0;
        NLEBBE3D::NODE(8, 2) = 0;
        NLEBBE3D::NODE(9, 0) = 9;
        NLEBBE3D::NODE(9, 1) = 0;
        NLEBBE3D::NODE(9, 2) = 0;
        NLEBBE3D::NODE(10, 0) = 10;
        NLEBBE3D::NODE(10, 1) = 0;
        NLEBBE3D::NODE(10, 2) = 0;

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

        NLEBBE3D::E = 1e7;
        NLEBBE3D::nu = 0.0001;
        NLEBBE3D::Bp = 1;
        NLEBBE3D::Hp = 1;
        NLEBBE3D::Zx = 0;
        NLEBBE3D::Zy = 0;
        NLEBBE3D::Zz = 1;

        NLEBBE3D::loadnode = 6;

        NLEBBE3D::DIA = 0.0002;
    }
    //Slave
    else if (choice == 2)
    {
        NLEBBE3D::NNODE = 11;
        NLEBBE3D::NELEM = 5;
        NLEBBE3D::NDOF = 6;
        NLEBBE3D::NLS = 24;
        NLEBBE3D::NEN = 3;

        NLEBBE3D::NODE = Eigen::MatrixXd::Zero(NLEBBE3D::NNODE, 3);
        NLEBBE3D::ELEM = Eigen::MatrixXd::Zero(NLEBBE3D::NELEM, 4);

        //Nodal Information
        NLEBBE3D::NODE(0, 0) = 0;
        NLEBBE3D::NODE(0, 1) = 0.001;
        NLEBBE3D::NODE(0, 2) = 0;
        NLEBBE3D::NODE(1, 0) = 1;
        NLEBBE3D::NODE(1, 1) = 0.001;
        NLEBBE3D::NODE(1, 2) = 0;
        NLEBBE3D::NODE(2, 0) = 2;
        NLEBBE3D::NODE(2, 1) = 0.001;
        NLEBBE3D::NODE(2, 2) = 0;
        NLEBBE3D::NODE(3, 0) = 3;
        NLEBBE3D::NODE(3, 1) = 0.001;
        NLEBBE3D::NODE(3, 2) = 0;
        NLEBBE3D::NODE(4, 0) = 4;
        NLEBBE3D::NODE(4, 1) = 0.001;
        NLEBBE3D::NODE(4, 2) = 0;
        NLEBBE3D::NODE(5, 0) = 5;
        NLEBBE3D::NODE(5, 1) = 0.001;
        NLEBBE3D::NODE(5, 2) = 0;
        NLEBBE3D::NODE(6, 0) = 6;
        NLEBBE3D::NODE(6, 1) = 0.001;
        NLEBBE3D::NODE(6, 2) = 0;
        NLEBBE3D::NODE(7, 0) = 7;
        NLEBBE3D::NODE(7, 1) = 0.001;
        NLEBBE3D::NODE(7, 2) = 0;
        NLEBBE3D::NODE(8, 0) = 8;
        NLEBBE3D::NODE(8, 1) = 0.001;
        NLEBBE3D::NODE(8, 2) = 0;
        NLEBBE3D::NODE(9, 0) = 9;
        NLEBBE3D::NODE(9, 1) = 0.001;
        NLEBBE3D::NODE(9, 2) = 0;
        NLEBBE3D::NODE(10, 0) = 10;
        NLEBBE3D::NODE(10, 1) = 0.001;
        NLEBBE3D::NODE(10, 2) = 0;

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

        NLEBBE3D::E = 1e7;
        NLEBBE3D::nu = 0.0001;
        NLEBBE3D::Bp = 1;
        NLEBBE3D::Hp = 1;
        NLEBBE3D::Zx = 0;
        NLEBBE3D::Zy = 0;
        NLEBBE3D::Zz = 1;

        NLEBBE3D::loadnode = 6;

        NLEBBE3D::DIA = 0.0002;

    }
}

BeamContact::BeamContact(const int nbeams, std::string str)
{
    BeamContact::NBEAMS = nbeams;

    BeamContact::epsilon = 1e10;

    if (str == "STS")
        BeamContact::NEN = 4;
    else if (str == "NTN")
        BeamContact::NEN = 2;

    BeamContact::GlobalELEM = { 0, 1 };
    
    

}

//Validation case for segment-to-segment contact with 3-D EB beams
/*NonLinearEulerBernouliBeamElement3D ReadEBBE3DElement(double ms)
{
    struct NonLinearEulerBernouliBeamElement3D M;
    //Master
    if (ms == 1)
    {
        M.NNODE = 3;
        M.NELEM = 1;
        M.NDOF = 6;
        M.NLS = 40;
        M.NEN = 3;

        M.NODE = Eigen::MatrixXd::Zero(M.NNODE, 3);
        M.ELEM = Eigen::MatrixXd::Zero(M.NELEM, 4);

        //Nodal Information
        M.NODE(0, 0) = 0;
        M.NODE(0, 1) = 0;
        M.NODE(0, 2) = 0;
        M.NODE(1, 0) = 5;
        M.NODE(1, 1) = 0;
        M.NODE(1, 2) = 0;
        M.NODE(2, 0) = 10;
        M.NODE(2, 1) = 0;
        M.NODE(2, 2) = 0;

        //Element connectivity
        M.ELEM(0, 0) = 1;
        M.ELEM(0, 1) = 1;
        M.ELEM(0, 2) = 2;
        M.ELEM(0, 3) = 3;

        M.E = 1e7;
        M.nu = 0.0001;
        M.Bp = 1;
        M.Hp = 1;
        M.Zx = 0;
        M.Zy = 0;
        M.Zz = 1;

        //node where point load is acting
        M.load = 2;

        M.D = 0.0001;

        M.epsilon = 1e10;
    }
    //Slave
    else if (ms == 2)
    {
        M.NNODE = 3;
        M.NELEM = 1;
        M.NDOF = 6;
        M.NLS = 40;
        M.NEN = 3;

        M.NODE = Eigen::MatrixXd::Zero(M.NNODE, 3);
        M.ELEM = Eigen::MatrixXd::Zero(M.NELEM, 4);

        //Nodal Information
        M.NODE(0, 0) = 0;
        M.NODE(0, 1) = 0.001;
        M.NODE(0, 2) = 0;
        M.NODE(1, 0) = 5;
        M.NODE(1, 1) = 0.001;
        M.NODE(1, 2) = 0;
        M.NODE(2, 0) = 10;
        M.NODE(2, 1) = 0.001;
        M.NODE(2, 2) = 0;

        //Element connectivity
        M.ELEM(0, 0) = 1;
        M.ELEM(0, 1) = 1;
        M.ELEM(0, 2) = 2;
        M.ELEM(0, 3) = 3;

        M.E = 1e7;
        M.nu = 0.0001;
        M.Bp = 1;
        M.Hp = 1;
        M.Zx = 0;
        M.Zy = 0;
        M.Zz = 1;

        M.load = 2;

        M.D = 0.0001;

        M.epsilon = 1e10;
    }
    return M;
}*/






/*
Mesh ReadInpFile(std::string filename)
{
    struct Mesh M;
    std::fstream file;
    file.open(filename.c_str());
    if(file.is_open())
    {
        //Read the file
        while(!file.eof())
        {
            std::string str;
            getline(file,str);
            std::stringstream s(str);
            s>>str;
            //Type of analysis
            if(std::strcmp(str.c_str(),"*STATIC"))
            {
                getline(file,str);
                M.choice = stoi(str);
            }
            else
            {
                std::cout<<"You seem to have forgotten to enter the kind of analysis that you want me to do."<<std::endl;
                std::cout<<"Use the numbers listed below to tell me about the same"<<std::endl;
                std::cout<<"1. Linear Static Analysis using Bar Elements"<<std::endl;
                std::cout<<"2. Non Linear Static Anlaysis using Bar Elements"<<std::endl;
                std::cout<<"3. Non Linear Static Analysis using Euler Bernoulli Beam Elements"<<std::endl;
                std::cout<<"4. Non Linear Static Analysis using VAM Beam Elements"<<std::endl;
                std::cout<<"Check the sample .inp file to know more."<<std::endl;

            }
            //Number of Bodies
            if(std::strcmp(str.c_str(),"*NBODIES"))
            {
                getline(file, str);
                M.NBODIES = stoi(str);
            }
            else
            {
                std::cout<<"You seem to have forgotten to enter the number of bodies to taken for me to solve."<<std::endl;
                std::cout<<"I'm assuming there's only one body"<<std::endl;
                M.NBODIES = 1;
            }
            //Get nodal data.
            if(std::strcmp(str.c_str(),"*NODE"))
            {
                std::cout<<"Accessing Nodal Data"<<std::endl;
                int i=0;
                do
                {
                    getline(file, str);
                    std::stringstream s(str);
                    s>>M.NODE(i,0)>>M.NODE(i,1)>>M.NODE(i,3);
                    i++;
                } while (!std::strcmp(str.c_str()," "));
            }
            else
            {
                std::cout<<"I dont' see any nodal data"<<std::endl;
                std::cout<<"Check the sample .inp file to know more"<<std::endl;
            }
            //Get Elemental Data
            if(std::strcmp(str.c_str(), "*ELEMENT"))
            {
                std::cout<<"Accessing Connectivity Data"<<std::endl;
                int i=0;
                do
                {
                    getline(file, str);
                    std::stringstream s(str);
                    s>>M.ELEM(i,0)>>M.ELEM(i,1)>>M.ELEM(i,2);
                    i++;
                } while (!std::strcmp(str.c_str()," "));

            }
            else
            {
                std::cout<<"I dont' see any connectivity data"<<std::endl;
                std::cout<<"Check the sample .inp file to know more"<<std::endl;
            }
            //Get Material Parameters
            if(std::strcmp(str.c_str(),"*MATERIAL"))
            {
                do
                {
                    if(std::strcmp(str.c_str(),"NMAT"))
                    {
                        getline(file, str);
                        std::stringstream s(str);
                        s>>M.NMAT;
                    }
                    std::string temp;
                    getline(file, str);
                    std::stringstream s(str);
                    s>>temp;
                    if(strcmp(temp.c_str(),"*ELASTIC"))
                    {
                        s>>temp>>temp>>temp;
                        if(strcmp(temp.c_str(),"ISOTROPIC"))
                        {
                            int i=0;
                            do
                            {
                                getline(file, str);
                                std::stringstream s(str);
                                s>>M.Mat[i].Iso.E>>M.Mat[i].Iso.nu;
                                i++;
                            }
                            while (!std::strcmp(str.c_str()," "));

                        }
                        else if(strcmp(temp.c_str(),"ORTHOTROPIC"))
                        {
                            getline(file, str);
                            std::stringstream s(str);
                        }

                    }


                }while (!std::strcmp(str.c_str()," "));
            }
            else
            {

            }
            //Get Loading and Boundary conditions
            //Element Sets

        }
    }
    else
    {
        std::cout<<"File not found"<<std::endl;
        std::cout<<"Maybe there's a spelling error"<<std::endl;
        std::cout<<"Maybe the slash is reversed"<<std::endl;
    }

}*/

#endif

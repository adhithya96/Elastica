#ifndef VARIABLES_H
#define VARIABLES_H
#include "Variables.h"

//Mesh ReadInpFile()
//{
//    struct Mesh M;
//    M.NBODIES = 2;
//    M.NNODE = 12;
//    M.NELEM = 10;
//    M.NMAT = 2;
//    M.NDOF = 3;
//    M.NLS = 20;

//    M.NODE = Eigen::MatrixXd::Zero(M.NNODE,3);
//    M.ELEM = Eigen::MatrixXd::Zero(M.NELEM, 3);
//    M.MAT = Eigen::MatrixXd::Zero(M.NMAT, M.NBODIES);
//    M.CNODE = Eigen::MatrixXd::Zero(6,2);
//    M.LOAD = Eigen::MatrixXd::Zero(2,3);
////    Nodal Information
//    M.NODE(0,0) = 0;
//    M.NODE(0,1) = 0;
//    M.NODE(0,2) = 0;
//    M.NODE(1,0) = 0.508/2;
//    M.NODE(1,1) = 0;
//    M.NODE(1,2) = 0;
//    M.NODE(2,0) = 1.016/2;
//    M.NODE(2,1) = 0;
//    M.NODE(2,2) = 0;
//    M.NODE(3,0) = 1.524/2;
//    M.NODE(3,1) = 0;
//    M.NODE(3,2) = 0;
//    M.NODE(4,0) = 2.032/2;
//    M.NODE(4,1) = 0;
//    M.NODE(4,2) = 0;
//    M.NODE(5,0) = 2.54/2;
//    M.NODE(5,1) = 0;
//    M.NODE(5,2) = 0;
//    M.NODE(6,0) = 0;
//    M.NODE(6,1) = 1;
//    M.NODE(6,2) = 0;
//    M.NODE(7,0) = 0.508/2;
//    M.NODE(7,1) = 1;
//    M.NODE(7,2) = 0;
//    M.NODE(8,0) = 1.016/2;
//    M.NODE(8,1) = 1;
//    M.NODE(8,2) = 0;
//    M.NODE(9,0) = 1.524/2;
//    M.NODE(9,1) = 1;
//    M.NODE(9,2) = 0;
//    M.NODE(10,0) = 2.032/2;
//    M.NODE(10,1) = 1;
//    M.NODE(10,2) = 0;
//    M.NODE(11,0) = 2.54/2;
//    M.NODE(11,1) = 1;
//    M.NODE(11,2) = 0;
////    Connectivity information
//    M.ELEM(0,0) = 1;
//    M.ELEM(0,1) = 1;
//    M.ELEM(0,2) = 2;
//    M.ELEM(1,0) = 1;
//    M.ELEM(1,1) = 2;
//    M.ELEM(1,2) = 3;
//    M.ELEM(2,0) = 1;
//    M.ELEM(2,1) = 3;
//    M.ELEM(2,2) = 4;
//    M.ELEM(3,0) = 1;
//    M.ELEM(3,1) = 4;
//    M.ELEM(3,2) = 5;
//    M.ELEM(4,0) = 1;
//    M.ELEM(4,1) = 5;
//    M.ELEM(4,2) = 6;
//    M.ELEM(5,0) = 2;
//    M.ELEM(5,1) = 7;
//    M.ELEM(5,2) = 8;
//    M.ELEM(6,0) = 2;
//    M.ELEM(6,1) = 8;
//    M.ELEM(6,2) = 9;
//    M.ELEM(7,0) = 2;
//    M.ELEM(7,1) = 9;
//    M.ELEM(7,2) = 10;
//    M.ELEM(8,0) = 2;
//    M.ELEM(8,1) = 10;
//    M.ELEM(8,2) = 11;
//    M.ELEM(9,0) = 2;
//    M.ELEM(9,1) = 11;
//    M.ELEM(9,2) = 12;
////    Beam material information
//    M.MAT(0,0) = 2.068423e+11;
//    M.MAT(0,1) = 0.33;
//    M.MAT(1,0) = 2.068423e+11;
//    M.MAT(1,1) = 0.33;
////    Beam cross section information
//    M.b1 = 0.0254;
//    M.h1 = 0.0254;
//    M.b2 = 0.0254;
//    M.h2 = 0.0254;
////    Boundary conditions and constrainst
//    //There are 3 dof for each node. Numbering followed is as follows:
//    //Dispalcement
//        //0. Horizontal deflection
//        //1. Vertical Deflection
//        //2. Angle of rotation
//    //Force/Neumann
//        //1. Axial Force
//        //2. Vertical Force
//        //3. Moment
//    //Displacement Boundary
//    M.CNODE(0,0) = 1;
//    M.CNODE(0,1) = 1;
//    M.CNODE(1,0) = 6;
//    M.CNODE(1,1) = 0;
//    M.CNODE(2,0) = 6;
//    M.CNODE(2,1) = 2;
//    M.CNODE(3,0) = 7;
//    M.CNODE(3,1) = 1;
//    M.CNODE(4,0) = 12;
//    M.CNODE(4,1) = 0;
//    M.CNODE(5,0) = 12;
//    M.CNODE(5,1) = 2;
////Loading conditions
//    //Element sets have to be written in inp file for distributed load
////    Load can act in 3 ways depending on the no. of degrees of freedom
////    0. Axial Force
////    1. Vertical Force
////    2. Moment
////    Element sets are described for each load indicating the nodes on which load act
//  /*  M.LOAD(0,0) = 0;
//    M.LOAD(0,1) = -175.126835;
//    M.LOAD(0,2) = 0;
//    M.LOAD(1,0) = 0;
//    M.LOAD(1,1) = -175.126835;
//    M.LOAD(1,2) = 0;
//    M.choice = 3;
////Element Sets
//    M.ELSET(0,0) = 1;
//    M.ELSET(0,1) = 2;
//    M.ELSET(0,2) = 3;
//    M.ELSET(0,3) = 4;
//    M.ELSET(0,4) = 5;
//    M.ELSET(0,5) = 6;
//    M.ELSET(1,0) = 7;
//    M.ELSET(1,1) = 8;
//    M.ELSET(1,2) = 9;
//    M.ELSET(1,3) = 10;
//    M.ELSET(1,4) = 11;
//    M.ELSET(1,5) = 12;*/
//    M.choice = 3;
//    return M;
//}

/*void LoadVector(Eigen::MatrixXd force, Eigen::MatrixXd& ELSET, Eigen::MatrixXd &LOAD, int NNODE)
{
    //Initialie load vector to zero
    for(int i=0;i<NNODE;i++)
    {
        LOAD(i,0) = 0;
        LOAD(i,1) = 0;
        LOAD(i,2) = 0;
    }
//    initialize load vector with the values of load given by the user
    for(int j=0;j<ELSET.rows();j++)
        for(int i=0;i<ELSET.cols();i++)
        {
            LOAD((int)ELSET(j,i)-1,0) = force(j,0);
            LOAD((int)ELSET(j,i)-1,1) = force(j,1);
            LOAD((int)ELSET(j,i)-1,2) = force(j,2);
        }
}
*/
/*
GaussPoints GaussQuadraturePoints()
{
    struct GaussPoints GP;
    //Three Point
    GP.w1 = 5/9;
    GP.w2 = 8/9;
    GP.w3 = 5/9;
    GP.x1 = -sqrt(3/5);
    GP.x2 = 0;
    GP.x3 = sqrt(3/5);

    return GP;
}*/


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

//Single Body
NonLinearEulerBernouliBeamElement ReadNLEBBEFile()
{
    struct NonLinearEulerBernouliBeamElement NLEBBE;

    NLEBBE.NBODIES = 1;
    NLEBBE.NNODE = 5;
    NLEBBE.NELEM = 4;
    NLEBBE.NMAT = 2;
    NLEBBE.NDOF = 3;
    NLEBBE.NLS = 10;

    NLEBBE.NODE = Eigen::MatrixXd::Zero(NLEBBE.NNODE,3);
    NLEBBE.ELEM = Eigen::MatrixXd::Zero(NLEBBE.NELEM, 3);
    NLEBBE.MAT = Eigen::MatrixXd::Zero(NLEBBE.NMAT, 2);

    NLEBBE.MAT(0,0) = 30e6;
    NLEBBE.MAT(0,1) = 0.33;
    NLEBBE.MAT(1,0) = 30e6;
    NLEBBE.MAT(1,1) = 0.33;

    NLEBBE.CS.choice = "RECT";
    NLEBBE.CS.Rect.height = 1;
    NLEBBE.CS.Rect.width = 1;
    NLEBBE.CS.choice = "CIRCLE";
    NLEBBE.CS.Cir.radius = 1;


    NLEBBE.NODE(0,0) = 0;
    NLEBBE.NODE(0,1) = 0;
    NLEBBE.NODE(0,2) = 0;
    NLEBBE.NODE(1,0) = 12.5;
    NLEBBE.NODE(1,1) = 0;
    NLEBBE.NODE(1,2) = 0;
    NLEBBE.NODE(2,0) = 25;
    NLEBBE.NODE(2,1) = 0;
    NLEBBE.NODE(2,2) = 0;
    NLEBBE.NODE(3,0) = 37.5;
    NLEBBE.NODE(3,1) = 0;
    NLEBBE.NODE(3,2) = 0;
    NLEBBE.NODE(4,0) = 50;
    NLEBBE.NODE(4,1) = 0;
    NLEBBE.NODE(4,2) = 0;

    NLEBBE.ELEM(0,0) = 1;
    NLEBBE.ELEM(0,1) = 1;
    NLEBBE.ELEM(0,2) = 2;
    NLEBBE.ELEM(1,0) = 1;
    NLEBBE.ELEM(1,1) = 2;
    NLEBBE.ELEM(1,2) = 3;
    NLEBBE.ELEM(2,0) = 1;
    NLEBBE.ELEM(2,1) = 3;
    NLEBBE.ELEM(2,2) = 4;
    NLEBBE.ELEM(3,0) = 1;
    NLEBBE.ELEM(3,1) = 4;
    NLEBBE.ELEM(3,2) = 5;

    NLEBBE.CNODE = Eigen::MatrixXd::Zero(3, 2);
    NLEBBE.CNODE(0, 0) = 1;
    NLEBBE.CNODE(0, 1) = 2;
    NLEBBE.CNODE(1, 0) = 5;
    NLEBBE.CNODE(1, 1) = 1;
    NLEBBE.CNODE(2, 0) = 5;
    NLEBBE.CNODE(2, 1) = 3;

    NLEBBE.vf = 1;
    NLEBBE.af = 0;

    return NLEBBE;
}

//Multiple Bodies
/*NonLinearEulerBernouliBeamElement ReadNLEBBEFile()
{
    struct NonLinearEulerBernouliBeamElement NLEBBE;

    NLEBBE.NBODIES = 2;
    NLEBBE.NNODE = 8;
    NLEBBE.NELEM = 6;
    NLEBBE.NMAT = 2;
    NLEBBE.NDOF = 3;
    NLEBBE.NLS = 1;
    NLEBBE.NDM = 2;

    NLEBBE.NODE = Eigen::MatrixXd::Zero(NLEBBE.NNODE, 2);
    NLEBBE.ELEM = Eigen::MatrixXd::Zero(NLEBBE.NELEM, 3);
    NLEBBE.MAT = Eigen::MatrixXd::Zero(NLEBBE.NMAT, 2);

    NLEBBE.MAT(0, 0) = 30e6;
    NLEBBE.MAT(0, 1) = 0.33;
    NLEBBE.MAT(1, 0) = 30e6;
    NLEBBE.MAT(1, 1) = 0.33;

    NLEBBE.CS.choice = "CIRCLE";
    NLEBBE.CS.Cir.radius = 1;

    //Beam 1
    NLEBBE.NODE(0, 0) = 0;
    NLEBBE.NODE(0, 1) = 0;
    NLEBBE.NODE(1, 0) = 0.25;
    NLEBBE.NODE(1, 1) = 0.25;
    NLEBBE.NODE(2, 0) = 0.75;
    NLEBBE.NODE(2, 1) = 0.75;
    NLEBBE.NODE(3, 0) = 1;
    NLEBBE.NODE(3, 1) = 1;

    //Beam 2
    NLEBBE.NODE(4, 0) = 0;
    NLEBBE.NODE(4, 1) = 1;
    NLEBBE.NODE(5, 0) = 0.25;
    NLEBBE.NODE(5, 1) = 0.75;
    NLEBBE.NODE(6, 0) = 0.75;
    NLEBBE.NODE(6, 1) = 0.25;
    NLEBBE.NODE(7, 0) = 1;
    NLEBBE.NODE(7, 1) = 0;

    //Connectivity 1
    NLEBBE.ELEM(0, 0) = 1;
    NLEBBE.ELEM(0, 1) = 1;
    NLEBBE.ELEM(0, 2) = 2;
    NLEBBE.ELEM(1, 0) = 1;
    NLEBBE.ELEM(1, 1) = 2;
    NLEBBE.ELEM(1, 2) = 3;
    NLEBBE.ELEM(2, 0) = 1;
    NLEBBE.ELEM(2, 1) = 3;
    NLEBBE.ELEM(2, 2) = 4;

    //Connectivity 2
    NLEBBE.ELEM(3, 0) = 2;
    NLEBBE.ELEM(3, 1) = 5;
    NLEBBE.ELEM(3, 2) = 6;
    NLEBBE.ELEM(4, 0) = 2;
    NLEBBE.ELEM(4, 1) = 6;
    NLEBBE.ELEM(4, 2) = 7;
    NLEBBE.ELEM(5, 0) = 2;
    NLEBBE.ELEM(5, 1) = 7;
    NLEBBE.ELEM(5, 2) = 8;

    NLEBBE.CNODE = Eigen::MatrixXd::Zero(3, 2);
    NLEBBE.CNODE(0, 0) = 1;
    NLEBBE.CNODE(0, 1) = 1;
    NLEBBE.CNODE(1, 0) = 1;
    NLEBBE.CNODE(1, 1) = 2;
    NLEBBE.CNODE(2, 0) = 5;
    NLEBBE.CNODE(2, 1) = 3;

    NLEBBE.vf = 1;
    NLEBBE.af = 0;

    return NLEBBE;
}*/


VAMBeamElement ReadVAMBEFile()
{
    struct VAMBeamElement M;
    M.NBODIES = 1;
    M.NNODE = 12;
    M.NELEM = 11;
    M.NMAT = 1;
    M.NDOF = 12;
    M.NLS = 20;
    //M.CMP.np = 16;
    M.NCS = 1;

    M.NODE = Eigen::MatrixXd::Zero(M.NNODE, 3);
    M.ELEM = Eigen::MatrixXd::Zero(M.NELEM, 3);
    //M.CMP.Orient = Eigen::VectorXd::Zero(M.CMP.np);
    M.CMP.inittwist = Eigen::VectorXd::Zero(3);
    //M.LOAD = Eigen::MatrixXd::Zero(1,3);
//    Nodal Information
    M.NODE(0, 0) = 0;
    M.NODE(0, 1) = 0;
    M.NODE(0, 2) = 0;
    M.NODE(1, 0) = 1.97;
    M.NODE(1, 1) = 0;
    M.NODE(1, 2) = 0;
    M.NODE(2, 0) = 3.94;
    M.NODE(2, 1) = 0;
    M.NODE(2, 2) = 0;
    M.NODE(3, 0) = 5.91;
    M.NODE(3, 1) = 0;
    M.NODE(3, 2) = 0;
    M.NODE(4, 0) = 7.88;
    M.NODE(4, 1) = 0;
    M.NODE(4, 2) = 0;
    M.NODE(5, 0) = 9.85;
    M.NODE(5, 1) = 0;
    M.NODE(5, 2) = 0;
    M.NODE(6, 0) = 11.82;
    M.NODE(6, 1) = 0;
    M.NODE(6, 2) = 0;
    M.NODE(7, 0) = 13.79;
    M.NODE(7, 1) = 0;
    M.NODE(7, 2) = 0;
    M.NODE(8, 0) = 15.76;
    M.NODE(8, 1) = 0;
    M.NODE(8, 2) = 0;
    M.NODE(9, 0) = 17.73;
    M.NODE(9, 1) = 0;
    M.NODE(9, 2) = 0;
    M.NODE(10, 0) = 19.7;
    M.NODE(10, 1) = 0;
    M.NODE(10, 2) = 0;
    M.NODE(11, 0) = 21.67;
    M.NODE(11, 1) = 0;
    M.NODE(11, 2) = 0;

    M.ELEM(0, 0) = 1;
    M.ELEM(0, 1) = 1;
    M.ELEM(0, 2) = 2;
    M.ELEM(1, 0) = 1;
    M.ELEM(1, 1) = 2;
    M.ELEM(1, 2) = 3;
    M.ELEM(2, 0) = 1;
    M.ELEM(2, 1) = 3;
    M.ELEM(2, 2) = 4;
    M.ELEM(3, 0) = 1;
    M.ELEM(3, 1) = 4;
    M.ELEM(3, 2) = 5;
    M.ELEM(4, 0) = 1;
    M.ELEM(4, 1) = 5;
    M.ELEM(4, 2) = 6;
    M.ELEM(5, 0) = 1;
    M.ELEM(5, 1) = 6;
    M.ELEM(5, 2) = 7;
    M.ELEM(6, 0) = 1;
    M.ELEM(6, 1) = 7;
    M.ELEM(6, 2) = 8;
    M.ELEM(7, 0) = 1;
    M.ELEM(7, 1) = 8;
    M.ELEM(7, 2) = 9;
    M.ELEM(8, 0) = 1;
    M.ELEM(8, 1) = 9;
    M.ELEM(8, 2) = 10;
    M.ELEM(9, 0) = 1;
    M.ELEM(9, 1) = 10;
    M.ELEM(9, 2) = 11;
    M.ELEM(10, 0) = 1;
    M.ELEM(10, 1) = 11;
    M.ELEM(10, 2) = 12;

    /*M.NODE(0, 0) = 0;
    M.NODE(0, 1) = 0;
    M.NODE(0, 2) = 0;
    M.NODE(1, 0) = 84.67e-3;
    M.NODE(1, 1) = 0;
    M.NODE(1, 2) = 0;
    M.NODE(2, 0) = 169.33e-3;
    M.NODE(2, 1) = 0;
    M.NODE(2, 2) = 0;
    M.NODE(3, 0) = 254e-3;
    M.NODE(3, 1) = 0;
    M.NODE(3, 2) = 0;


    M.ELEM(0, 0) = 1;
    M.ELEM(0, 1) = 1;
    M.ELEM(0, 2) = 2;
    M.ELEM(1, 0) = 1;
    M.ELEM(1, 1) = 2;
    M.ELEM(1, 2) = 3;
    M.ELEM(2, 0) = 1;
    M.ELEM(2, 1) = 3;
    M.ELEM(2, 2) = 4;*/
    //    Beam material information
        //M.MAT(0,0) = 30e6;
        //M.MAT(0,1) = 0.33;
    /*M.CMP.E11 = 135.6e9;
    M.CMP.E22 = 9.9e9;
    M.CMP.E33 = 9.9e9;
    M.CMP.G12 = 4.2e9;
    M.CMP.G13 = 4.2e9;
    M.CMP.G23 = 3.3e9;
    M.CMP.nu12 = 0.3;
    M.CMP.nu13 = 0.3;
    M.CMP.nu23 = 0.5;

    M.CMP.Orient(0) = 20;
    M.CMP.Orient(1) = 20;
    M.CMP.Orient(2) = -70;
    M.CMP.Orient(3) = -70;
    M.CMP.Orient(4) = -70;
    M.CMP.Orient(5) = -70;
    M.CMP.Orient(6) = 20;
    M.CMP.Orient(7) = 20;
    M.CMP.Orient(8) = -20;
    M.CMP.Orient(9) = -20;
    M.CMP.Orient(10) = 70;
    M.CMP.Orient(11) = 70;
    M.CMP.Orient(12) = 70;
    M.CMP.Orient(13) = 70;
    M.CMP.Orient(14) = -20;
    M.CMP.Orient(15) = -20;*/

    M.CMP.inittwist(0) = 0;
    M.CMP.inittwist(1) = 0;
    M.CMP.inittwist(2) = 0;

    //    Beam cross section information
    /*M.CS.Rect.width = 0.0254;
    M.CS.Rect.height = 1.168e-3;*/

    //Loading
    M.B.a1 = Eigen::VectorXd::Zero(3);
    M.B.b1 = Eigen::VectorXd::Zero(3);
    M.B.aN = Eigen::VectorXd::Zero(3);
    M.B.bN = Eigen::VectorXd::Zero(3);

    M.B.x1 = 1;
    M.B.y1 = 2;
    M.B.xN = 3;
    M.B.yN = 4;

    M.B.aN(1) = 0;

    return M;
}

NonLinearEulerBernouliBeamElement3D ReadEBBE3DElement()
{
    struct NonLinearEulerBernouliBeamElement3D M;
    M.NNODE = 17;
    M.NELEM = 16;
    M.NDOF = 6;
    M.NLS = 12;
    //M.CMP.np = 16;

    M.NODE = Eigen::MatrixXd::Zero(M.NNODE, 3);
    M.ELEM = Eigen::MatrixXd::Zero(M.NELEM, 3);
    //M.CMP.Orient = Eigen::VectorXd::Zero(M.CMP.np);
    //M.LOAD = Eigen::MatrixXd::Zero(1,3);

    //Nodal Information
    M.NODE(0, 0) = 0;
    M.NODE(0, 1) = 0;
    M.NODE(0, 2) = 0;
    M.NODE(1, 0) = 0.1204;
    M.NODE(1, 1) = 4.9067;
    M.NODE(1, 2) = 0;
    M.NODE(2, 0) = 0.4815;
    M.NODE(2, 1) = 9.8017;
    M.NODE(2, 2) = 0;
    M.NODE(3, 0) = 1.082;
    M.NODE(3, 1) = 14.673;
    M.NODE(3, 2) = 0;
    M.NODE(4, 0) = 1.9214;
    M.NODE(4, 1) = 19.509;
    M.NODE(4, 2) = 0;
    M.NODE(5, 0) = 2.9968;
    M.NODE(5, 1) = 24.298;
    M.NODE(5, 2) = 0;
    M.NODE(6, 0) = 4.306;
    M.NODE(6, 1) = 29.03;
    M.NODE(6, 2) = 0;
    M.NODE(7, 0) = 5.845;
    M.NODE(7, 1) = 33.69;
    M.NODE(7, 2) = 0;
    M.NODE(8, 0) = 7.612;
    M.NODE(8, 1) = 38.27;
    M.NODE(8, 2) = 0;
    M.NODE(9, 0) = 9.601;
    M.NODE(9, 1) = 42.75;
    M.NODE(9, 2) = 0;
    M.NODE(10, 0) = 11.807;
    M.NODE(10, 1) = 47.14;
    M.NODE(10, 2) = 0;
    M.NODE(11, 0) = 14.23;
    M.NODE(11, 1) = 51.41;
    M.NODE(11, 2) = 0;
    M.NODE(12, 0) = 16.853;
    M.NODE(12, 1) = 55.56;
    M.NODE(12, 2) = 0;
    M.NODE(13, 0) = 19.68;
    M.NODE(13, 1) = 59.57;
    M.NODE(13, 2) = 0;
    M.NODE(14, 0) = 22.698;
    M.NODE(14, 1) = 63.44;
    M.NODE(14, 2) = 0;
    M.NODE(15, 0) = 25.9;
    M.NODE(15, 1) = 67.15;
    M.NODE(15, 2) = 0;
    M.NODE(16, 0) = 29.3;
    M.NODE(16, 1) = 70.71;
    M.NODE(16, 2) = 0;

    //Element connectivity
    M.ELEM(0, 0) = 1;
    M.ELEM(0, 1) = 1;
    M.ELEM(0, 2) = 2;
    M.ELEM(1, 0) = 1;
    M.ELEM(1, 1) = 2;
    M.ELEM(1, 2) = 3;
    M.ELEM(2, 0) = 1;
    M.ELEM(2, 1) = 3;
    M.ELEM(2, 2) = 4;
    M.ELEM(3, 0) = 1;
    M.ELEM(3, 1) = 4;
    M.ELEM(3, 2) = 5;
    M.ELEM(4, 0) = 1;
    M.ELEM(4, 1) = 5;
    M.ELEM(4, 2) = 6;
    M.ELEM(5, 0) = 1;
    M.ELEM(5, 1) = 6;
    M.ELEM(5, 2) = 7;
    M.ELEM(6, 0) = 1;
    M.ELEM(6, 1) = 7;
    M.ELEM(6, 2) = 8;
    M.ELEM(7, 0) = 1;
    M.ELEM(7, 1) = 8;
    M.ELEM(7, 2) = 9;
    M.ELEM(8, 0) = 1;
    M.ELEM(8, 1) = 9;
    M.ELEM(8, 2) = 10;
    M.ELEM(9, 0) = 1;
    M.ELEM(9, 1) = 10;
    M.ELEM(9, 2) = 11;
    M.ELEM(10, 0) = 1;
    M.ELEM(10, 1) = 11;
    M.ELEM(10, 2) = 12;
    M.ELEM(11, 0) = 1;
    M.ELEM(11, 1) = 12;
    M.ELEM(11, 2) = 13;
    M.ELEM(12, 0) = 1;
    M.ELEM(12, 1) = 13;
    M.ELEM(12, 2) = 14;
    M.ELEM(13, 0) = 1;
    M.ELEM(13, 1) = 14;
    M.ELEM(13, 2) = 15;
    M.ELEM(14, 0) = 1;
    M.ELEM(14, 1) = 15;
    M.ELEM(14, 2) = 16;
    M.ELEM(15, 0) = 1;
    M.ELEM(15, 1) = 16;
    M.ELEM(15, 2) = 17;

    M.E = 1e7;
    M.nu = 0.0001;
    M.Bp = 1;
    M.Hp = 1;
    M.Zx = 0;
    M.Zy = 0;
    M.Zz = 1;

    return M;
}

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

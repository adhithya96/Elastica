#include "variables.h"
#include<iostream>

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

Mesh ReadInpFile()
{
    struct Mesh M;
    M.NBODIES = 1;
    M.NNODE = 5;
    M.NELEM = 4;
    M.NMAT = 1;
    M.NDOF = 3;
    M.NLS = 20;

    M.NODE = Eigen::MatrixXd::Zero(M.NNODE,3);
    M.ELEM = Eigen::MatrixXd::Zero(M.NELEM, 3);
    M.MAT = Eigen::MatrixXd::Zero(M.NMAT, 2);
    M.CNODE = Eigen::MatrixXd::Zero(3,2);
    //M.LOAD = Eigen::MatrixXd::Zero(1,3);
//    Nodal Information
    M.NODE(0,0) = 0;
    M.NODE(0,1) = 0;
    M.NODE(0,2) = 0;
    M.NODE(1,0) = 12.5;
    M.NODE(1,1) = 0;
    M.NODE(1,2) = 0;
    M.NODE(2,0) = 25;
    M.NODE(2,1) = 0;
    M.NODE(2,2) = 0;
    M.NODE(3,0) = 37.5;
    M.NODE(3,1) = 0;
    M.NODE(3,2) = 0;
    M.NODE(4,0) = 50;
    M.NODE(4,1) = 0;
    M.NODE(4,2) = 0;
//    Connectivity information
    M.ELEM(0,0) = 1;
    M.ELEM(0,1) = 1;
    M.ELEM(0,2) = 2;
    M.ELEM(1,0) = 1;
    M.ELEM(1,1) = 2;
    M.ELEM(1,2) = 3;
    M.ELEM(2,0) = 1;
    M.ELEM(2,1) = 3;
    M.ELEM(2,2) = 4;
    M.ELEM(3,0) = 1;
    M.ELEM(3,1) = 4;
    M.ELEM(3,2) = 5;
//    Beam material information
    M.MAT(0,0) = 30e6;
    M.MAT(0,1) = 0.33;
//    Beam cross section information
    M.b1 = 1;
    M.h1 = 1;
//    Boundary conditions and constrainst
    //There are 3 dof for each node. Numbering followed is as follows:
    //Dispalcement
        //0. Horizontal deflection
        //1. Vertical Deflection
        //2. Angle of rotation
    //Force/Neumann
        //1. Axial Force
        //2. Vertical Force
        //3. Moment
    //Displacement Boundary
    M.CNODE(0,0) = 1;
    M.CNODE(0,1) = 1;
    M.CNODE(1,0) = 5;
    M.CNODE(1,1) = 0;
    M.CNODE(2,0) = 5;
    M.CNODE(2,1) = 2;
//Loading conditions
    //Element sets have to be written in inp file for distributed load
//    Load can act in 3 ways depending on the no. of degrees of freedom
//    0. Axial Force
//    1. Vertical Force
//    2. Moment
//    Element sets are described for each load indicating the nodes on which load act
    //M.LOAD(0,0) = -1;
    //M.LOAD(0,1) = "DLOAD";
//    M.LOAD(1,0) = 0;
//    M.LOAD(1,1) = -175.126835;
//    M.LOAD(1,2) = 0;
    M.choice = 3;
//Element Sets
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
//    M.ELSET(1,5) = 12;
    M.vf = 1;
    M.af = 0;
    return M;
}

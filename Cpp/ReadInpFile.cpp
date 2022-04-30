#include "variables.h"

Mesh ReadInpFile()
{
    struct Mesh M;
    M.NBODIES = 2;
    M.NNODE = 12;
    M.NELEM = 10;
    M.NMAT = 2;
    M.NDOF = 3;
    M.NLS = 20;


    M.NODE = Eigen::MatrixXd::Zero(M.NNODE,3);
    M.ELEM = Eigen::MatrixXd::Zero(M.NELEM, 3);
    M.MAT = Eigen::MatrixXd::Zero(M.NMAT, M.NBODIES);

    M.CNODE = Eigen::MatrixXd(M.NNODE, M.NDOF);
    M.FORCE = Eigen::MatrixXd(M.NNODE, M.NDOF);
    for(int i=0;i<M.NNODE;i++)
        for(int j=0;j<M.NDOF;j++)
        {
            M.CNODE(i,j) = -1;
            M.FORCE(i,j) = 0;
        }
//    Nodal Information
    M.NODE(0,0) = 0;
    M.NODE(0,1) = 0;
    M.NODE(0,2) = 0;
    M.NODE(1,0) = 0.2;
    M.NODE(1,1) = 0;
    M.NODE(1,2) = 0;
    M.NODE(2,0) = 0.4;
    M.NODE(2,1) = 0;
    M.NODE(2,2) = 0;
    M.NODE(3,0) = 0.6;
    M.NODE(3,1) = 0;
    M.NODE(3,2) = 0;
    M.NODE(4,0) = 0.8;
    M.NODE(4,1) = 0;
    M.NODE(4,2) = 0;
    M.NODE(5,0) = 1;
    M.NODE(5,1) = 0;
    M.NODE(5,2) = 0;
    M.NODE(6,0) = 0;
    M.NODE(6,1) = 1;
    M.NODE(6,2) = 0;
    M.NODE(7,0) = 0.2;
    M.NODE(7,1) = 1;
    M.NODE(7,2) = 0;
    M.NODE(8,0) = 0.4;
    M.NODE(8,1) = 1;
    M.NODE(8,2) = 0;
    M.NODE(9,0) = 0.6;
    M.NODE(9,1) = 1;
    M.NODE(9,2) = 0;
    M.NODE(10,0) = 0.8;
    M.NODE(10,1) = 1;
    M.NODE(10,2) = 0;
    M.NODE(11,0) = 1;
    M.NODE(11,1) = 1;
    M.NODE(11,2) = 0;
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
    M.ELEM(4,0) = 1;
    M.ELEM(4,1) = 5;
    M.ELEM(4,2) = 6;
    M.ELEM(5,0) = 2;
    M.ELEM(5,1) = 7;
    M.ELEM(5,2) = 8;
    M.ELEM(6,0) = 2;
    M.ELEM(6,1) = 8;
    M.ELEM(6,2) = 9;
    M.ELEM(7,0) = 2;
    M.ELEM(7,1) = 9;
    M.ELEM(7,2) = 10;
    M.ELEM(8,0) = 2;
    M.ELEM(8,1) = 10;
    M.ELEM(8,2) = 11;
    M.ELEM(9,0) = 2;
    M.ELEM(9,1) = 11;
    M.ELEM(9,2) = 12;
//    Beam material information
    M.MAT(0,0) = 2.0e9;
    M.MAT(0,1) = 0.33;
    M.MAT(1,0) = 2.0e9;
    M.MAT(1,1) = 0.33;
//    Beam cross section information
    M.b1 = 1;
    M.h1 = 1;
    M.b2 = 1;
    M.h2 = 1;
//    Boundary conditions and constrainst
    //Displacement Boundary
    M.CNODE(0,0) = 0;
    M.CNODE(0,1) = 0;
    M.FORCE(11,1) = -1;
    M.CNODE(6,0) = 0;
    M.CNODE(6,1) = 0;
    M.FORCE(12,1) = -1;

    M.maxiter = 100;
    M.fiter = 10;
    M.choice = 3;

    return M;
}

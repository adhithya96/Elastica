/*
#ifndef VARIABLES_H
#define VARIABLES_H
#include "Variables.h"


void StiffnessMatrix_NLBE(Eigen::MatrixXd& k,Eigen::VectorXd U, int StartNode, int EndNode, double h)
{
    //Symmetric part
    k(0,0) = (U(StartNode)+U(EndNode))/(2*h);
    k(0,1) = -(U(StartNode)+U(EndNode))/(2*h);
    k(1,0) = -(U(StartNode)+U(EndNode))/(2*h);
    k(1,1) = (U(StartNode)+U(EndNode))/(2*h);
}

void TangentStiffnessMatrix_NLBE(Eigen::MatrixXd& t, Eigen::VectorXd U, int StartNode, int EndNode, double h)
{
    t(0,0) = (U(StartNode)+U(EndNode))/(2*h)-(U(EndNode)-U(StartNode))/(2*h);
    t(0,1) = -(U(StartNode)+U(EndNode))/(2*h)-(U(EndNode)-U(StartNode))/(2*h);
    t(1,0) = -(U(StartNode)+U(EndNode))/(2*h)+(U(EndNode)-U(StartNode))/(2*h);
    t(1,1) = (U(StartNode)+U(EndNode))/(2*h)+(U(EndNode)-U(StartNode))/(2*h);
}

void LocalForceVec_NLBE(Eigen::VectorXd& f, int StartNode, int EndNode, double h)
{
    //f Vector
    f(0) =(-h/2);
    f(1) =(-h/2);
}

//Boundary conditions
void ApplyConstraints_NLBE(Eigen::SparseMatrix<double, Eigen::RowMajor>& T, Eigen::VectorXd &U, Eigen::VectorXd &R, int NNODE, int cNODE, double value)
{
    U(cNODE) = value;
    for(int i=0;i<NNODE;i++)
        if(i!=cNODE)
            T.coeffRef(cNODE,i)=0;
    T.coeffRef(cNODE,cNODE)=1;
    R.coeffRef(cNODE) = 0;
}
#endif
*/

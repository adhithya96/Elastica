#ifndef VARIABLES_H
#define VARIABLES_H
#include"Variables.h"


void StiffnessMatrix_NLBE(Eigen::MatrixXd& k, Eigen::VectorXd U, int StartNode, int EndNode, double h)
{
    //Symmetric part
    k(0, 0) = (U(StartNode) + U(EndNode)) / (2 * h);
    k(0, 1) = -(U(StartNode) + U(EndNode)) / (2 * h);
    k(1, 0) = -(U(StartNode) + U(EndNode)) / (2 * h);
    k(1, 1) = (U(StartNode) + U(EndNode)) / (2 * h);
}

void TangentStiffnessMatrix_NLBE(Eigen::MatrixXd& t, Eigen::VectorXd U, int StartNode, int EndNode, double h)
{
    t(0, 0) = (U(StartNode) + U(EndNode)) / (2 * h) - (U(EndNode) - U(StartNode)) / (2 * h);
    t(0, 1) = -(U(StartNode) + U(EndNode)) / (2 * h) - (U(EndNode) - U(StartNode)) / (2 * h);
    t(1, 0) = -(U(StartNode) + U(EndNode)) / (2 * h) + (U(EndNode) - U(StartNode)) / (2 * h);
    t(1, 1) = (U(StartNode) + U(EndNode)) / (2 * h) + (U(EndNode) - U(StartNode)) / (2 * h);
}

void LocalForceVec_NLBE(Eigen::VectorXd& f, int StartNode, int EndNode, double h)
{
    //f Vector
    f(0) = (-h / 2);
    f(1) = (-h / 2);
}

//Boundary conditions
void ApplyConstraints_NLBE(Eigen::SparseMatrix<double, Eigen::RowMajor>& T, Eigen::VectorXd& U, Eigen::VectorXd& R, int NNODE, Eigen::MatrixXd CNODE)
{
    for (int i = 0; i < CNODE.rows(); i++)
    {
        int nodenum = CNODE(i, 0);
        double value = CNODE(i, 1);
        for (int j = 0; j < NNODE; j++)
            T.coeffRef(nodenum - 1, i) = 0;
        T.coeffRef(nodenum - 1, nodenum - 1) = 1;
        R.coeffRef(nodenum - 1) = 0;
        U(nodenum - 1) = value;
    }
}

#endif
#ifndef VARIABLES_H
#define VARIABLES_H
#include"Variables.h"
//How to input loads and boundary conditions

void StiffnessMatrix_LBE(double A, double E, double xa, double xb, Eigen::MatrixXd& k)
{
    double h = xb - xa;
    k(0, 0) = A * E / h;
    k(0, 1) = -A * E / h;
    k(1, 0) = -A * E / h;
    k(1, 1) = A * E / h;
}

void LocalForceVec_LBE(double f1, double f2, Eigen::VectorXd& f)
{
    f(0) = f1;
    f(1) = f2;
}

void ApplyConstraints_LBE(Eigen::MatrixXd CNODE, int NNODE, Eigen::VectorXd F,
    Eigen::SparseMatrix<double, Eigen::ColMajor>& K)
{
    for (int j = 0; j < CNODE.rows(); j++)
    {
        int nodenum = CNODE(j, 0);
        double value = CNODE(j, 1);
        for (int i = 0; i < NNODE; i++)
            K.coeffRef(nodenum - 1, i) = 0;
        F(nodenum - 1) = 0;
        K.coeffRef(nodenum - 1, nodenum - 1) = 1;
    }

}
#endif
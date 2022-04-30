#include <QCoreApplication>
#include<iostream>
#include<cassert>
#include "variables.h"

int main(int argc, char *argv[])
{

    //Read variables from input file
    struct Mesh M = ReadInpFile();
    Eigen::VectorXd dU = Eigen::VectorXd::Zero(M.NDOF*M.NNODE);
    Eigen::VectorXd U = Eigen::VectorXd::Zero(M.NDOF*M.NNODE);
    Eigen::VectorXd error = Eigen::VectorXd::Zero(M.NDOF*M.NNODE);
    Eigen::MatrixXd k = Eigen::MatrixXd::Zero(6,6);
    Eigen::MatrixXd t = Eigen::MatrixXd::Zero(6,6);
    Eigen::VectorXd f = Eigen::VectorXd::Zero(6);
    Eigen::VectorXd F = Eigen::VectorXd::Zero(M.NDOF*M.NNODE);
    Eigen::VectorXd U_new = Eigen::VectorXd::Zero(M.NDOF*M.NNODE);
    Eigen::VectorXd R = Eigen::VectorXd::Zero(M.NDOF*M.NNODE);
    Eigen::SparseMatrix<double, Eigen::RowMajor> K, T;
    double max;
    int iter=0;
    int fiter = 1;


    //Loop for load increment
    do
    {
        for(int i = 0; i<M.NNODE;i++)
        {
            if(M.FORCE(i,0)!=0)
                M.FORCE(i,0) = M.FORCE(i,0)*fiter/M.NLS;
            else if(M.NDOF==2 && M.FORCE(i,1)!=0)
                M.FORCE(i,1) = M.FORCE(i,1)*fiter/M.NLS;
            else if(M.NDOF==3 && M.FORCE(i,2)!=0)
                M.FORCE(i,2) = M.FORCE(i,2)*fiter/M.NLS;
            else
                continue;
        }
        if(fiter<M.NLS)
            fiter++;
    //Newton Raphson loop
        do
        {
            //Initialize all matrices to zero
            U_new.setZero();
            F.setZero();
            R.setZero();
            K.setZero();
            T.setZero();

            //Construct Global Stiffness Matrix

            for(int i=0;i<M.NELEM;i++)
            {
                int xa = M.ELEM(i,0);
                int xb = M.ELEM(i,1);
                //if(choice==1)
                //   BarElement(A,E,l,K,StartNode,EndNode);
    //            Bar Element or Poisson Equation Given in Reddy
                if(M.choice==2)
                {
                    NonLinearStiffnessMatrix(k,U,xa,xb,M.h1);
                    K.coeffRef(xa,xa) += k(0,0);
                    K.coeffRef(xa,xb) += k(0,1);
                    K.coeffRef(xb,xa) += k(1,0);
                    K.coeffRef(xb,xb) += k(1,1);

                    TangentStiffnessMatrix(t,U,xa,xb,M.h1);
                    T.coeffRef(xa,xa) += t(0,0);
                    T.coeffRef(xa,xb) += t(0,1);
                    T.coeffRef(xb,xa) += t(1,0);
                    T.coeffRef(xb,xb) += t(1,1);

                    ForceVector(f,xa,xb,M.h1);
                    F.coeffRef(xa) += f(0);
                    F.coeffRef(xb) += f(1);
                }
    //            Euler Bernoulli Beam in 2D
                if(M.choice==3)
                {
                    int a = M.ELEM(i,1)-1;
                    int b = M.ELEM(i,2)-1;
                    double xa = M.NODE((int)a,0);
                    double xb = M.NODE((int)b,0);
                    double ya = M.NODE((int)a,1);
                    double yb = M.NODE((int)b,1);
                    double E = M.MAT((int)M.ELEM(i,0)-1,0);
                    double nu = M.MAT((int)M.ELEM(i,0)-1,1);
                    //Local stiffness matrix
                    StiffnessEulerBernoulli(k,xa,xb,E,nu,M.b1,M.h1,U,a,b);
                    //Each node has 3 dof. The global node is expanded into 3 nodes to store the values of

                    //Local Force Vector

                    //Tangent Stiffness matrix
                    TangentStiffnessEulerBernoulli(t,k,xa,xb,E,nu,M.b1,M.h1,U,a,b);
                    RearrangeElementStiffness(k);


                    //Assembly stiffness matrix
                    for(int i=0;i<3;i++)
                        for(int j=0;j<3;j++)
                            K.coeffRef(3*a+i,3*a+j) = k(i,j);
                    for(int i=0;i<3;i++)
                        for(int j=0;j<3;j++)
                            K.coeffRef(3*b+i,3*b+j) = k(i+2,j+2);

                    //Assembly Tangent stiffness matrix
                    for(int i=0;i<3;i++)
                        for(int j=0;j<3;j++)
                            T.coeffRef(3*a+i,3*a+j) = t(i,j);
                    for(int i=0;i<3;i++)
                        for(int j=0;j<3;j++)
                            T.coeffRef(3*b+i,3*b+j) = t(i+2,j+2);

                }
            }

            //Assemble global stiffness and tangent matrix
            //std::cout<<"Stiffness matrix"<<K<<std::endl;

            //std::cout<<"Displacements"<<U<<std::endl;

            //std::cout<<"Force Vector"<<F<<std::endl;

            //Residual Matrix
            R = K*U-F;


            //std::cout<<"Residual Matrix"<<"  "<<R<<std::endl;

           //Apply Constraints
            //ApplyConstraints(T,U,R,NNODE,cNODE1,sqrt(2));
            //ApplyConstraints(T,U,R,NNODE,cNODE2,sqrt(2));
            K.makeCompressed();
            T.makeCompressed();

            std::cout<<"Tangent Stiffness Matrix"<<"  "<<T<<std::endl;

            std::cout<<"Residual Matrix"<<R<<std::endl;

            //Matrix Solution
            Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::RowMajor>> solver;
            solver.analyzePattern(T);
            solver.factorize(T); //LU decomposition

            assert(solver.info()==Eigen::Success);

            dU = solver.solve(-R);

            std::cout<<"Increment in displacement"<<"  "<<dU<<std::endl;

            //Next iteration
            for(int i=0;i<M.NNODE;i++)
                U_new(i) = U(i) + dU(i);

            //Error calculation
            for(int i=0;i<M.NNODE;i++)
                error(i) = abs((U_new(i)-U(i))/U(i));
            max = error(0);
            for(int i=0;i<M.NNODE;i++)
                if(max<error(i))
                    max = error(i);

            //Assignment for next iteration
            for(int i=0;i<M.NNODE;i++)
                U(i) = U_new(i);
            iter++;

            //std::cout<<"Displacements"<<U<<std::endl;
        }while(max>pow(10,-8)&&iter<M.maxiter);
    //Load step increment
    }while(M.fiter<M.NLS);
    std::cout << "Displacements  "<< std::endl << U << std::endl;

    if(iter==M.maxiter)
        std::cout<<"Maximum iteration limit reached"<<std::endl<<"Don't kid yourself"<<std::endl<<"Solution did not converge"<<std::endl;
}

/*void BarElement(double A, double E, double l, Eigen::SparseMatrix<double, Eigen::RowMajor> K,std::vector<Eigen::Triplet<double>>matTriplets, int StartNode, int EndNode)
{
    double temp = A*E/l;
    //std::cout<<temp<<std::endl;
    K(StartNode,StartNode) += temp;
    //std::cout<<K(StartNode,StartNode)<<std::endl;
    K(StartNode,StartNode+1) += -temp;
    K(EndNode,StartNode) += -temp;
    K(EndNode,StartNode+1) += temp;
}*/


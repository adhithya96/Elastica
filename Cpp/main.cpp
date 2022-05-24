//Input force vector and then divide it based on increments
//Based on the degree of integrand, no. of integration points are selected appropriately.
//How to input distributed loads?
//What happens in post processing after finding displacements?
//How to implement reduced integration

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
    Eigen::SparseMatrix<double, Eigen::ColMajor> K(M.NNODE*M.NDOF, M.NNODE*M.NDOF), T(M.NNODE*M.NDOF, M.NNODE*M.NDOF);
    //Eigen::MatrixXd LOAD;
    double max;
    int iter = 1;
    int fiter = 1;
    int maxiter = 30;
    //LoadVector(M.LOAD,M.ELSET,LOAD,M.NNODE);

    //Loop for load increment
    do
    {
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


                //if(choice==1)
                //   BarElement(A,E,l,K,StartNode,EndNode);
    //            Bar Element or Poisson Equation Given in Reddy
                if(M.choice==2)
                {
                    for(int i=0;i<(int)M.NELEM;i++)
                    {
                        int xa = M.ELEM(i,0);
                        int xb = M.ELEM(i,1);
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
                }
    //            Euler Bernoulli Beam in 2D
                if(M.choice==3)
                {
                    for(int i=0;i<(int)M.NELEM;i++)
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

                        //std::vector<Eigen::Triplet<double>>kTriplets;
                        //Each node has 3 dof. The global node is expanded into 3 nodes to store the values of

                        //Tangent Stiffness matrix
                        TangentStiffnessEulerBernoulli(t,k,xa,xb,E,nu,M.b1,M.h1,U,a,b);

                        //Local Force Vector
                        ForceVec(f,xa,xb,a,b, M.vf, M.af);
//                        std::cout<<"Element Stiffness Matrix before rearranging"<<std::endl;
//                        std::cout<<k<<std::endl;

//                        std::cout<<"Element Tangent Matrix before rerarranging"<<std::endl;
//                        std::cout<<t<<std::endl;

//                        std::cout<<"Force Vector before rearranging"<<std::endl;
//                        std::cout<<f<<std::endl;

                        RearrangeElementMatrices(k,t,f);

//                        std::cout<<"Element Stiffness Matrix"<<std::endl;
//                        std::cout<<k<<std::endl;

//                        std::cout<<"Element Tangent Matrix"<<std::endl;
//                        std::cout<<t<<std::endl;

//                        std::cout<<"Force Vector"<<std::endl;
//                        std::cout<<f<<std::endl;

                        //Assembly stiffness matrix
                        for(int i=0;i<3;i++)
                            for(int j=0;j<3;j++)
                            {
                                K.coeffRef(3*a+i,3*a+j) += k(i,j);
                                K.coeffRef(3*b+i,3*b+j) += k(i+3,j+3);
                                K.coeffRef(3*a+i,3*b+j) += k(i,j+3);
                                K.coeffRef(3*b+i,3*a+j) += k(i+3,j);
                            }

                        //Assembly Tangent stiffness matrix
                        for(int i=0;i<3;i++)
                            for(int j=0;j<3;j++)
                            {
                                T.coeffRef(3*a+i,3*a+j) += t(i,j);
                                T.coeffRef(3*b+i,3*b+j) += t(i+3,j+3);
                                T.coeffRef(3*a+i,3*b+j) += t(i,j+3);
                                T.coeffRef(3*b+i,3*a+j) += t(i+3,j);
                            }

                        //Assembly Force Vector
                        for(int i=0;i<3;i++)
                        {
                            F(3*a+i) += f(i);
                            F(3*b+i) += f(i+3);
                        }
                    }
                }
                //VAM - Variational Asymptotic method
                if(M.choice==4)
                {

                }
            //Assemble global stiffness and tangent matrix
//            std::cout<<"Stiffness matrix"<<K<<std::endl;

//            std::cout<<"Displacements"<<U<<std::endl;

//            std::cout<<"Force Vector"<<F<<std::endl;

            //Residual Matrix
            R = K*U-F;


//            std::cout<<"Residual Matrix before applying boundary conditions"<<std::endl<<R<<std::endl;

//            std::cout<<"Stiffness matrix before applying bcs"<<std::endl<<K<<std::endl;

           //Apply Constraints
            //ApplyConstraints(T,U,R,NNODE,cNODE1,sqrt(2));
            //ApplyConstraints(T,U,R,NNODE,cNODE2,sqrt(2));
            ApplyConstraints2DEB(T,U,M.CNODE,M.NNODE,R);
            K.makeCompressed();
            T.makeCompressed();

//            std::cout<<"Tangent Stiffness Matrix"<<std::endl<<T<<std::endl;

//            std::cout<<"Residual Matrix"<<std::endl<<R<<std::endl;

            //Matrix Solution
            Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>> solver;
            solver.analyzePattern(T);
            solver.factorize(T); //LU decomposition

            assert(solver.info()==Eigen::Success);

            dU = solver.solve(-R);

//            std::cout<<"Increment in displacement"<<std::endl<<dU<<std::endl;

            //Next iteration
            for(int i=0;i<M.NNODE*M.NDOF;i++)
                U_new(i) = U(i) + dU(i);

            //Error calculation
            for(int i=0;i<M.NNODE*M.NDOF;i++)
                error(i) = abs((U_new(i)-U(i))/U(i));
            max = error(0);
            for(int i=0;i<M.NNODE*M.NDOF;i++)
                if(max<error(i))
                    max = error(i);

            //Assignment for next iteration
            for(int i=0;i<M.NNODE*M.NDOF;i++)
                U(i) = U_new(i);
            iter++;

            //std::cout<<"Displacements"<<U<<std::endl;
        }while(max>pow(10,-3)&&iter<maxiter);
        if(fiter<M.NLS)
            fiter++;
        if(iter==maxiter)
        {
            M.vf = M.vf-0.5;
            std::cout<<"Solution not converging for the given load"<<std::endl;
            std::cout<<"Trying with a load with less increment"<<std::endl;
        }
        std::cout<<M.vf<<std::endl;
        std::cout << "Displacements  "<< std::endl << U << std::endl;
//        std::cout<<"Force"<<std::endl<<F<<std::endl;
//        std::cout<<"Element Stiffness Matrix"<<std::endl<<k<<std::endl;
//        std::cout<<"Element Tangent Matrix"<<std::endl<<t<<std::endl;
        M.vf += 1;
    //Load step increment
    }while(fiter<M.NLS);

    if(iter==maxiter)
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



#ifndef VARIABLES_H
#define VARIABLES_H
#include <QCoreApplication>

#include "Variables.h"


int main()
{
    //Read variables from input file
    //std::string filename="Beam.inp";
//    std::cout<<"Enter the name of the input file"<<std::endl;
//    std::cin>>filename;
   // struct Mesh M = ReadInpFile(filename);
    int choice = 3;
    //LoadVector(M.LOAD,M.ELSET,LOAD,M.NNODE);
    // -------------------------------------
    //--------------------------------------
    //----------BAR ELEMENT LINEAR----------
    //--------------------------------------
    //--------------------------------------
    //This code works only for linear bar element with a concentrated load at the end
    /*if(choice==1)
    {
        LinearBarElement LBE = ReadLBEFile();

        Eigen::MatrixXd k = Eigen::MatrixXd::Zero(2,2);
        Eigen::VectorXd f = Eigen::VectorXd::Zero(2);
        Eigen::SparseMatrix<double, Eigen::ColMajor> K(LBE.NNODE, LBE.NNODE);
        Eigen::VectorXd F = Eigen::VectorXd::Zero(LBE.NNODE);
        Eigen::VectorXd U = Eigen::VectorXd::Zero(LBE.NNODE);
        double height = LBE.CS.height;
        double width = LBE.CS.width;
        double Area = height*width;
        for(int i=0;i<(int)LBE.NELEM;i++)
        {
            double E = LBE.E;
            int a = LBE.ELEM(i,1)-1;
            int b = LBE.ELEM(i,2)-1;
            int xa = LBE.NODE((int)a,0);
            int xb = LBE.NODE((int)b,0);
            StiffnessMatrix_LBE(Area,E,xa,xb,k);

            double f1 = LBE.LOAD.cload[a];
            double f2 = LBE.LOAD.cload[b];
            LocalForceVec_LBE(f1, f2, f);

            //Assembly of Global Stiffness Matrix
            K.coeffRef(a,a) += k(0,0);
            K.coeffRef(a,b) += k(0,1);
            K.coeffRef(b,a) += k(1,0);
            K.coeffRef(b,b) += k(1,1);

            //Assembly of Global Force Vector
            F(a) += f(0);
            F(b) += f(1);
        }

        //Apply constraints
        ApplyConstraints_LBE(LBE.CNODE,LBE.NNODE,F,K);

        //Solve for displacements
        Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>> solver;
        solver.analyzePattern(K);
        solver.factorize(K); //LU decomposition

        assert(solver.info()==Eigen::Success);

        U = solver.solve(F);
        std::cout<<"Displacements"<<std::endl;
        std::cout<<U<<std::endl;
    }
    //--------------------------------------
    //--------------------------------------
    //-------BAR ELEMENT NONLINEAR 1D------
    //--------------------------------------
    //--------------------------------------
    else if(choice==2)
    {
        NonLinearBarElement NLBE = ReadNLBEFile();

        Eigen::VectorXd dU = Eigen::VectorXd::Zero(NLBE.NDOF*NLBE.NNODE);
        Eigen::VectorXd U = Eigen::VectorXd::Zero(NLBE.NDOF*NLBE.NNODE);
        Eigen::VectorXd error = Eigen::VectorXd::Zero(NLBE.NDOF*NLBE.NNODE);
        Eigen::MatrixXd k = Eigen::MatrixXd::Zero(6,6);
        Eigen::MatrixXd t = Eigen::MatrixXd::Zero(6,6);
        Eigen::VectorXd f = Eigen::VectorXd::Zero(6);
        Eigen::VectorXd F = Eigen::VectorXd::Zero(NLBE.NDOF*NLBE.NNODE);
        Eigen::VectorXd U_new = Eigen::VectorXd::Zero(NLBE.NDOF*NLBE.NNODE);
        Eigen::VectorXd R = Eigen::VectorXd::Zero(NLBE.NDOF*NLBE.NNODE);
        Eigen::SparseMatrix<double, Eigen::ColMajor> K(NLBE.NNODE*NLBE.NDOF, NLBE.NNODE*NLBE.NDOF),
                T(NLBE.NNODE*NLBE.NDOF, NLBE.NNODE*NLBE.NDOF);
        //Eigen::MatrixXd LOAD;
        double max;
        int iter = 1;
        int maxiter = 30;
        do
        {
            //Initialize all matrices to zero
            U_new.setZero();
            F.setZero();
            R.setZero();
            K.setZero();
            T.setZero();

            //Construct Global Stiffness Matrix

            for(int i=0;i<(int)NLBE.NELEM;i++)
            {
                int xa = NLBE.ELEM(i,0);
                int xb = NLBE.ELEM(i,1);
                StiffnessMatrix_NLBE(k,U,xa,xb,NLBE.CS.height);
                K.coeffRef(xa,xa) += k(0,0);
                K.coeffRef(xa,xb) += k(0,1);
                K.coeffRef(xb,xa) += k(1,0);
                K.coeffRef(xb,xb) += k(1,1);

                TangentStiffnessMatrix_NLBE(t,U,xa,xb,NLBE.CS.height);
                T.coeffRef(xa,xa) += t(0,0);
                T.coeffRef(xa,xb) += t(0,1);
                T.coeffRef(xb,xa) += t(1,0);
                T.coeffRef(xb,xb) += t(1,1);

                LocalForceVec_NLBE(f,xa,xb,NLBE.CS.height);
                F.coeffRef(xa) += f(0);
                F.coeffRef(xb) += f(1);
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
            //ApplyConstraints_NLBE(T,U,R,NLBE.NNODE,NLBE.CNODE,)
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
            for(int i=0;i<NLBE.NNODE*NLBE.NDOF;i++)
                U_new(i) = U(i) + dU(i);

            //Error calculation
            for(int i=0;i<NLBE.NNODE*NLBE.NDOF;i++)
                error(i) = abs((U_new(i)-U(i))/U(i));
            max = error(0);
            for(int i=0;i<NLBE.NNODE*NLBE.NDOF;i++)
                if(max<error(i))
                    max = error(i);

            //Assignment for next iteration
            for(int i=0;i<NLBE.NNODE*NLBE.NDOF;i++)
                U(i) = U_new(i);
            iter++;

            //std::cout<<"Displacements"<<U<<std::endl;
        }while(max>pow(10,-3)&&iter<maxiter);
    }
    //----------------------------------------------
    //----------------------------------------------
    //---------EULER BERNOULLI NONLINEAR 2D---------
    //----------------------------------------------
    //----------------------------------------------
    //The problem description for which the Euler Bernouli beam code works is a unit distributed load on a cantilever beam
    //To know more look at Introduction to Nonlinear FEM by JN Reddy
    else if(choice==3)
    {
        NonLinearEulerBernouliBeamElement NLEBBE = ReadNLEBBEFile();
        Eigen::VectorXd dU = Eigen::VectorXd::Zero(NLEBBE.NDOF*NLEBBE.NNODE);
        Eigen::VectorXd U = Eigen::VectorXd::Zero(NLEBBE.NDOF*NLEBBE.NNODE);
        Eigen::VectorXd error = Eigen::VectorXd::Zero(NLEBBE.NDOF*NLEBBE.NNODE);
        Eigen::MatrixXd k = Eigen::MatrixXd::Zero(6,6);
        Eigen::MatrixXd t = Eigen::MatrixXd::Zero(6,6);
        Eigen::VectorXd f = Eigen::VectorXd::Zero(6);
        Eigen::VectorXd F = Eigen::VectorXd::Zero(NLEBBE.NDOF*NLEBBE.NNODE);
        Eigen::VectorXd U_new = Eigen::VectorXd::Zero(NLEBBE.NDOF*NLEBBE.NNODE);
        Eigen::VectorXd R = Eigen::VectorXd::Zero(NLEBBE.NDOF*NLEBBE.NNODE);
        Eigen::SparseMatrix<double, Eigen::ColMajor> K(NLEBBE.NNODE*NLEBBE.NDOF, NLEBBE.NNODE*NLEBBE.NDOF),
                T(NLEBBE.NNODE*NLEBBE.NDOF, NLEBBE.NNODE*NLEBBE.NDOF);
        //Eigen::MatrixXd LOAD;
        double max;
        int iter = 1;
        int fiter = 1;
        int maxiter = 30;
        //Loop for force iterations.
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
                for(int i=0;i<(int)NLEBBE.NELEM;i++)
                {
                    int a = NLEBBE.ELEM(i,1)-1;
                    int b = NLEBBE.ELEM(i,2)-1;
                    double xa = NLEBBE.NODE((int)a,0);
                    double xb = NLEBBE.NODE((int)b,0);
                    double E = NLEBBE.MAT[(int)NLEBBE.ELEM(i,0)-1,0].E;
                    double nu = NLEBBE.MAT[(int)NLEBBE.ELEM(i,0)-1,0].nu;
                    //Local stiffness matrix
                    StiffnessMatrix_NLEBBE(k,xa,xb,E,nu,NLEBBE.CS.width,NLEBBE.CS.height,U,a,b);

                    //std::vector<Eigen::Triplet<double>>kTriplets;
                    //Each node has 3 dof. The global node is expanded into 3 nodes to store the values of

                    //Tangent Stiffness matrix
                    TangentStiffnessMatrix_NLEBBE(t,k,xa,xb,E,nu,NLEBBE.CS.width,NLEBBE.CS.height,U,a,b);

                    //Local Force Vector
                    LocalFoceVec_NLEBBE(f,xa,xb,a,b, NLEBBE.vf, NLEBBE.af);
        //                        std::cout<<"Element Stiffness Matrix before rearranging"<<std::endl;
        //                        std::cout<<k<<std::endl;

        //                        std::cout<<"Element Tangent Matrix before rerarranging"<<std::endl;
        //                        std::cout<<t<<std::endl;

        //                        std::cout<<"Force Vector before rearranging"<<std::endl;
        //                        std::cout<<f<<std::endl;

                    RearrangeElementStiffness_NLEBBE(k,t,f);

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
                ApplyConstraints_NLEBBE(T,U,NLEBBE.CNODE,NLEBBE.NNODE,R);
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
                for(int i=0;i<NLEBBE.NNODE*NLEBBE.NDOF;i++)
                    U_new(i) = U(i) + dU(i);

                //Error calculation
                for(int i=0;i<NLEBBE.NNODE*NLEBBE.NDOF;i++)
                    error(i) = abs((U_new(i)-U(i))/U(i));
                max = error(0);
                for(int i=0;i<NLEBBE.NNODE*NLEBBE.NDOF;i++)
                    if(max<error(i))
                        max = error(i);

                //Assignment for next iteration
                for(int i=0;i<NLEBBE.NNODE*NLEBBE.NDOF;i++)
                    U(i) = U_new(i);
                iter++;

                //std::cout<<"Displacements"<<U<<std::endl;
            }while(max>pow(10,-3)&&iter<maxiter);
            if(fiter<NLEBBE.NLS)
                fiter++;
            if(iter==maxiter)
            {
                NLEBBE.vf = NLEBBE.vf-0.5;
                std::cout<<"Solution not converging for the given load"<<std::endl;
                std::cout<<"Trying with a load with less increment"<<std::endl;
            }
            std::cout<<NLEBBE.vf<<std::endl;
            std::cout << "Displacements  "<< std::endl << U << std::endl;
    //        std::cout<<"Force"<<std::endl<<F<<std::endl;
    //        std::cout<<"Element Stiffness Matrix"<<std::endl<<k<<std::endl;
    //        std::cout<<"Element Tangent Matrix"<<std::endl<<t<<std::endl;
            NLEBBE.vf += 1;
        //Load step increment
        }while(fiter<NLEBBE.NLS);

        if(iter==maxiter)
            std::cout<<"Maximum iteration limit reached"<<std::endl<<"Don't kid yourself"<<std::endl<<"Solution did not converge"<<std::endl;
    }*/
    //----------------------------------------------
    //----------------------------------------------
    //-------------VAM Beam Element 3D--------------
    //----------------------------------------------
    //----------------------------------------------
    //GEBT
    if(choice==4)
    {
        VAMBeamElement VAMBE = ReadVAMBEFile();
        Eigen::VectorXd dU = Eigen::VectorXd::Zero(VAMBE.NDOF*VAMBE.NNODE);
        Eigen::VectorXd U = Eigen::VectorXd::Zero(VAMBE.NDOF*VAMBE.NNODE);
        Eigen::VectorXd error = Eigen::VectorXd::Zero(VAMBE.NDOF*VAMBE.NNODE);
        Eigen::MatrixXd k = Eigen::MatrixXd::Zero(6,6);
        Eigen::MatrixXd t = Eigen::MatrixXd::Zero(6,6);
        Eigen::VectorXd f = Eigen::VectorXd::Zero(6);
        Eigen::VectorXd F = Eigen::VectorXd::Zero(VAMBE.NDOF*VAMBE.NNODE);
        Eigen::VectorXd U_new = Eigen::VectorXd::Zero(VAMBE.NDOF*VAMBE.NNODE);
        Eigen::VectorXd R = Eigen::VectorXd::Zero(VAMBE.NDOF*VAMBE.NNODE);
        Eigen::SparseMatrix<double, Eigen::ColMajor> K(VAMBE.NNODE*VAMBE.NDOF, VAMBE.NNODE*VAMBE.NDOF),
                J(VAMBE.NNODE*VAMBE.NDOF, VAMBE.NNODE*VAMBE.NDOF);
        Eigen::VectorXd Tau = Eigen::VectorXd::Zero(VAMBE.NDOF*VAMBE.NNODE);
        //Eigen::MatrixXd LOAD;
        double max;
        int iter = 1;
        //int fiter = 1;
        int maxiter = 30;
        //Loop for force iterations.
    //    do
    //    {
            do
            {

                //Initialize all matrices to zero
                U_new.setZero();
                F.setZero();
                R.setZero();
                K.setZero();
                J.setZero();

                Eigen::MatrixXd S1 = Eigen::MatrixXd::Zero(9,9);

                S1 = StiffnessMatrix_VAM(VAMBE.CMP.E11,VAMBE.CMP.E22,VAMBE.CMP.nu12,VAMBE.CMP.G12,
                                        VAMBE.CMP.E22*VAMBE.CMP.nu12/VAMBE.CMP.E11,VAMBE.CS.width,VAMBE.CS.height,
                                        VAMBE.CMP.inittwist(0),VAMBE.CMP.Orient,VAMBE.CMP.np);

                Eigen::MatrixXd Seq_1 = Eigen::MatrixXd::Identity(4,4);
                Eigen::MatrixXd Seq_2 = Eigen::MatrixXd::Identity(4,4);
                Eigen::VectorXd U1(6),F1(6),U2(6),F2(6);

                for(int i=0; i<VAMBE.NELEM-1; i++)
                {
                    int a = VAMBE.ELEM(i,1)-1;
                    int b = VAMBE.ELEM(i,2)-1;
                    double xa = VAMBE.NODE((int)a,0);
                    double xb = VAMBE.NODE((int)b,0);
                    double h = xb-xa;

                    //Dont' know what Theta means but it is used in C_ab formula
                    double Theta = 30;

                    for(int j=0; j<5 ;j++)
                    {
                        U1(j) = U(12*a+j);
                        F1(j) = U(12*a+6+j);
                        U2(j) = U(12*b+j);
                        F2(j) = U(12*b+6+j);
                    }
                    //2D Cross-Section Analysis
                    //S is the cross sectional stiffness matrix obtained from the cross sectional analysis
                    //9*9 stiffness matrix is borrowed from "Non-classical e⁄ects in non-linear analysis of pretwisted
                    //anisotropic strips" by Harursampath et. al
                    //It is converted to 4*4 by taking derivative of one dimensional strain energy with respect to the strains.
                    //The matrix is in general converted to 6*6 to generalize the code such that there's 3 dofs at each node accounting
                    //for a total of 6 dofs for each element.
                    //This can be extended to Timoshenko formulation as well.
                    //Formula for 1D strain energy can be found in Eq. (31) in the paper.
                    Seq_1 = ClassicalStiffnessModel_VAM(S1,U1);
                    Seq_2 = ClassicalStiffnessModel_VAM(S1,U2);

                    //1D Non Linear Beam Analysis
                    //Implemented from GEBT: A general-purpose nonlinear analysis tool for composite beams
                    //by Wenbin Yu et. al
                    //The Geometrically Exact Beam Theory equation can be found in Eq. 5.50 of
                    //"Non linear Composite beam theory by Hodges"
                    //This equation is simplified using Finite elements in the paper by Wenbin Yu et. al
                    //The equations required for implementation are (45), (46) and (47) from the paper.
                    //They are solved using Newton Raphson method.
                    Eigen::VectorXd Element_R(12);
                    Eigen::MatrixXd Element_J(12,24);

                    //Element Residual
                    Element_R = Element_Residual(U1, U2, F1, F2, h, Seq_1, Seq_2, Theta);
                    //Element Jacobian
                    Element_J = Element_Jacobian(U1, U2, F1, F2, h, Seq_1, Seq_2, Theta);

                    //Assembly
                    for(int j=0;j<12;j++)
                        for(int k=0;k<12;k++)
                        {
                            J.coeffRef(12*a+j,12*a+k) = Element_J(j,k);
                            J.coeffRef(12*b+j,12*b+k) = Element_J(j,k+12);
                        }
                    for(int j=0;j<12;j++)
                    {
                        R.coeffRef(6*a+ j) = Element_R(j);
                        R.coeffRef(6*b+ j) = Element_R(j+6);
                    }
                }

                //Add Jacobian and Residual for equations in first and last node
                Eigen::VectorXd Element_R(6);
                Eigen::MatrixXd Element_J(12,12);
                //Node 0
                for(int j=0; j<5 ;j++)
                {
                    U1(j) = U(12*0+j);
                    F1(j) = U(12*0+6+j);
                }
                double xa = VAMBE.NODE(0,0);
                double xb = VAMBE.NODE(1,0);
                double h = xb-xa;

                Seq_1 = ClassicalStiffnessModel_VAM(S1,U1);

                double Theta = 30;
                Element_Residual(U1,F1,h,Seq_1,Theta,0,VAMBE.B);
                Element_Jacobian(U1,F1,h,Seq_1,Theta,0);

                //Assembly of Node 0 into Global Stiffness matrix
                for(int j=0;j<12;j++)
                    for(int k=0;k<12;k++)
                        J.coeffRef(12*0+j,12*0+k) = Element_J(j,k);
                for(int j=0;j<12;j++)
                    R.coeffRef(6*0+j) = Element_R(j);

                //Node N
                int lastnode = VAMBE.NNODE-1;
                for(int j=0; j<5 ;j++)
                {
                    U1(j) = U(12*lastnode+j);
                    F1(j) = U(12*lastnode+6+j);
                }
                xa = VAMBE.NODE(VAMBE.NNODE-2,0);
                xb = VAMBE.NODE(VAMBE.NNODE-1,0);
                h = xb-xa;
                Seq_1 = ClassicalStiffnessModel_VAM(S1,U1);
                Element_Residual(U1,F1,h,Seq_1,Theta,VAMBE.NNODE-1,VAMBE.B);
                Element_Jacobian(U1,F1,h,Seq_1,Theta,VAMBE.NNODE-1);

                //Assembly of Node n into global stiffness matrix
                for(int j=0;j<12;j++)
                    for(int k=0;k<12;k++)
                        J.coeffRef(12*lastnode+j,12*lastnode+k) = Element_J(j,k);
                for(int j=0;j<12;j++)
                    R.coeffRef(6*lastnode+j) = Element_R(j);

                //Solve the equation
                K.makeCompressed();
                J.makeCompressed();

                Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>> solver;
                solver.analyzePattern(J);
                solver.factorize(J); //LU decomposition

                //std::cout<<J<<std::endl;

                assert(solver.info()==Eigen::Success);

                dU = solver.solve(-R);

                //Calculate errors and update for next iteration
                for(int i=0;i<VAMBE.NNODE*VAMBE.NDOF;i++)
                    U_new(i) = U(i) + dU(i);

                //Error calculation
                for(int i=0;i<VAMBE.NNODE*VAMBE.NDOF;i++)
                    error(i) = abs((U_new(i)-U(i))/U(i));
                max = error(0);
                for(int i=0;i<VAMBE.NNODE*VAMBE.NDOF;i++)
                    if(max<error(i))
                        max = error(i);

                //Assignment for next iteration
                for(int i=0;i<VAMBE.NNODE*VAMBE.NDOF;i++)
                    U(i) = U_new(i);
                iter++;

                //After calculating the internal forces through 1D nonlinear beam analysis,
                //the strains are updated using the constitutive relation.
                Update_Strains(Tau,VAMBE.NNODE,S1,U);

             }while(max>pow(10,-3)&&iter<maxiter);
   //     }while(max>fiter);
    }
    else
        std::cout<<"Wrong option"<<std::endl;


}

#endif




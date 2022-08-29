
#ifndef VARIABLES_H
#define VARIABLES_H
#include "Variables.h"
#include<fstream>


int main()
{
    //Read variables from input file
    //std::string filename="Beam.inp";
//    std::cout<<"Enter the name of the input file"<<std::endl;
//    std::cin>>filename;
   // struct Mesh M = ReadInpFile(filename);
    int choice = 4;
    //LoadVector(M.LOAD,M.ELSET,LOAD,M.NNODE);
    // -------------------------------------
    //--------------------------------------
    //----------BAR ELEMENT LINEAR----------
    //--------------------------------------
    //--------------------------------------
    //This code works only for linear bar element with a concentrated loads.
    if(choice==1)
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
        int count = 0;
        for(int i=0;i<(int)LBE.NELEM;i++)
        {
            double E = LBE.MAT((int)LBE.ELEM(i, 0)-1, 0);
            int a = (int)LBE.ELEM(i,1)-1;
            int b = (int)LBE.ELEM(i,2)-1;
            double xa = LBE.NODE((int)a,0);
            double xb = LBE.NODE((int)b,0);
            StiffnessMatrix_LBE(Area,E,xa,xb,k);

            int nodenum = (int)LBE.LOAD(count, 0) - 1;
            double value = (int)LBE.LOAD(count, 1);
            if (a == nodenum)
            {
                LocalForceVec_LBE(value / 2, 0, f);
                count++;
            }
            else if (b == nodenum)
            {
                if (nodenum == LBE.NNODE - 1)
                    LocalForceVec_LBE(0, value, f);
                else 
                    LocalForceVec_LBE(0, value / 2, f);
                count++;
            }
            else
                LocalForceVec_LBE(0, 0, f);

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

        std::cout << "Stiffness Matrix" << K << std::endl;
        std::cout << "Force Vector" << F << std::endl;

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
    //This part of the code solves the non linear poison's equation given by
    //d(udu/dx) = f_0, 0<x<1
    //Check Chapter 3 from "Introduction to Non linear Finite Element Analysis" by JN Reddy
    //Examples 3.3.1, 3.3.2, and 3.4.1
    //Example 3.4.1 contains the results which can be compared
    else if(choice==2)
    {
        NonLinearBarElement NLBE = ReadNLBEFile();

        Eigen::VectorXd dU = Eigen::VectorXd::Zero(NLBE.NDOF*NLBE.NNODE);
        Eigen::VectorXd U = Eigen::VectorXd::Zero(NLBE.NDOF*NLBE.NNODE);
        Eigen::VectorXd error = Eigen::VectorXd::Zero(NLBE.NDOF*NLBE.NNODE);
        Eigen::MatrixXd k = Eigen::MatrixXd::Zero(2,2);
        Eigen::MatrixXd t = Eigen::MatrixXd::Zero(2,2);
        Eigen::VectorXd f = Eigen::VectorXd::Zero(2);
        Eigen::VectorXd F = Eigen::VectorXd::Zero(NLBE.NDOF*NLBE.NNODE);
        Eigen::VectorXd U_new = Eigen::VectorXd::Zero(NLBE.NDOF*NLBE.NNODE);
        Eigen::VectorXd R = Eigen::VectorXd::Zero(NLBE.NDOF*NLBE.NNODE);
        Eigen::SparseMatrix<double, Eigen::RowMajor> K(NLBE.NNODE*NLBE.NDOF, NLBE.NNODE*NLBE.NDOF),
                T(NLBE.NNODE*NLBE.NDOF, NLBE.NNODE*NLBE.NDOF);
        double max;
        int iter = 1;
        int maxiter = 30;
        
        //Initial guess 
        for (int i = 0; i < NLBE.NNODE; i++)
            U(i) = 1;

        do
        {
            //Initialize all matrices to zero
            U_new.setZero();
            F.setZero();
            R.setZero();
            K.setZero();
            T.setZero();

            //Construct Global Stiffness Matrix
            //Assemble global stiffness, tangent matrix and residual vector
            for(int i=0;i<(int)NLBE.NELEM;i++)
            {
                int a = (int)NLBE.ELEM(i,0)-1;
                int b = (int)NLBE.ELEM(i,1)-1;
                double xa = NLBE.NODE(a, 0);
                double xb = NLBE.NODE(b, 0);
                double h = xb - xa;

                StiffnessMatrix_NLBE(k, U, a, b, h);
                K.coeffRef(a,a) += k(0,0);
                K.coeffRef(a,b) += k(0,1);
                K.coeffRef(b,a) += k(1,0);
                K.coeffRef(b,b) += k(1,1);

                TangentStiffnessMatrix_NLBE(t, U, a, b, h);
                T.coeffRef(a,a) += t(0,0);
                T.coeffRef(a,b) += t(0,1);
                T.coeffRef(b,a) += t(1,0);
                T.coeffRef(b,b) += t(1,1);

                LocalForceVec_NLBE(f, a, b, h);
                F.coeffRef(a) += f(0);
                F.coeffRef(b) += f(1);
            }


            //Residual Matrix
            R = K*U-F;

           //Apply Constraints
            ApplyConstraints_NLBE(T, U, R, NLBE.NNODE, NLBE.CNODE);
            K.makeCompressed();
            T.makeCompressed();

            std::cout<<"Tangent Stiffness Matrix"<<std::endl<<T<<std::endl;

            std::cout << "Stiffness Matrix" << std::endl << K << std::endl;

            std::cout << "Residual" << std::endl << R << std::endl;

            //Matrix Solution
            Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::RowMajor>> solver;
            solver.analyzePattern(T);
            solver.factorize(T); //LU decomposition

            assert(solver.info()==Eigen::Success);

            dU = solver.solve(-R);

//            std::cout<<"Increment in displacement"<<std::endl<<dU<<std::endl;

            //Next iteration
            for(int i=0;i<NLBE.NNODE*NLBE.NDOF;i++)
                U_new(i) = U(i) + dU(i);

            std::cout << "Displacements" << U_new << std::endl;
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
    //To know more look at Chapter 4 of "Introduction to Nonlinear FEM by JN Reddy"
    //Default parameters are taken to solve Example 4.2.1. 
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
                    int a = (int)NLEBBE.ELEM(i,1)-1;
                    int b = (int)NLEBBE.ELEM(i,2)-1;
                    double xa = NLEBBE.NODE((int)a,0);
                    double xb = NLEBBE.NODE((int)b,0);
                    double E = NLEBBE.MAT((int)NLEBBE.ELEM(i,0)-1,0);
                    double nu = NLEBBE.MAT((int)NLEBBE.ELEM(i,0)-1,1);
                    //Local stiffness matrix
                    StiffnessMatrix_NLEBBE(k,xa,xb,E,nu,NLEBBE.CS.width,NLEBBE.CS.height,U,a,b);

                    //Tangent Stiffness matrix
                    TangentStiffnessMatrix_NLEBBE(t,k,xa,xb,E,nu,NLEBBE.CS.width,NLEBBE.CS.height,U,a,b);

                    //Local Force Vector
                    LocalFoceVec_NLEBBE(f,xa,xb,a,b, NLEBBE.vf, NLEBBE.af);

                    RearrangeElementStiffness_NLEBBE(k,t,f);

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

                //Residual Matrix
                R = K*U-F;

               //Apply Constraints
                ApplyConstraints_NLEBBE(T,U,NLEBBE.CNODE,NLEBBE.NNODE,R);
                K.makeCompressed();
                T.makeCompressed();

                //Matrix Solution
                Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>> solver;

                solver.analyzePattern(T);
                solver.factorize(T); //LU decomposition

                assert(solver.info()==Eigen::Success);

                dU = solver.solve(-R);

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

            }while(max>pow(10,-3)&&iter<maxiter);
            if(fiter<NLEBBE.NLS)
                fiter++;
            if(iter==maxiter)
            {
                NLEBBE.vf = NLEBBE.vf-0.5;
                std::cout<<"Solution not converging for the given load"<<std::endl;
                std::cout<<"Trying with a load with less increment"<<std::endl;
            }
            std::cout<<"Load  "<<NLEBBE.vf << std::endl;
            std::cout << "Displacements  "<< std::endl << U << std::endl;

            NLEBBE.vf += 1;
        //Load step increment
        }while(fiter<NLEBBE.NLS);

        if(iter==maxiter)
            std::cout<<"Maximum iteration limit reached"<<std::endl<<"Don't kid yourself"<<std::endl<<"Solution did not converge"<<std::endl;
    }
    //----------------------------------------------
    //----------------------------------------------
    //-------------VAM Beam Element 3D--------------
    //----------------------------------------------
    //----------------------------------------------
    //GEBT
    else if (choice == 4)
    {
        VAMBeamElement VAMBE = ReadVAMBEFile();
        Eigen::VectorXd dU = Eigen::VectorXd::Zero(VAMBE.NDOF * (VAMBE.NNODE + 1));
        Eigen::VectorXd U = Eigen::VectorXd::Zero(VAMBE.NDOF * (VAMBE.NNODE + 1));
        Eigen::VectorXd error = Eigen::VectorXd::Zero(VAMBE.NDOF * (VAMBE.NNODE + 1));
        Eigen::VectorXd F = Eigen::VectorXd::Zero(VAMBE.NDOF * (VAMBE.NNODE + 1));
        Eigen::VectorXd U_new = Eigen::VectorXd::Zero(VAMBE.NDOF * (VAMBE.NNODE + 1));
        Eigen::VectorXd R = Eigen::VectorXd::Zero(VAMBE.NDOF * (VAMBE.NNODE + 1));
        Eigen::SparseMatrix<double, Eigen::ColMajor> J((VAMBE.NNODE + 1)* VAMBE.NDOF, (VAMBE.NNODE + 1)* VAMBE.NDOF);
        Eigen::VectorXd Tau = Eigen::VectorXd::Zero(6 * (VAMBE.NNODE + 1));
        double max;
        int iter = 1;
        int fiter = 1;
        int maxiter = 30;
        std::fstream file1, file2;
        file1.open("G:/MTech_Aerospace/MTech_Res_Thesis/Cpp/TextileComposite_MicroMechanics/TextileComposite_MicroMechanics/TextileComposite_MicroMechanics/Result_Log.txt", std::fstream::in | std::fstream::out);
        file2.open("G:/MTech_Aerospace/MTech_Res_Thesis/Cpp/TextileComposite_MicroMechanics/TextileComposite_MicroMechanics/TextileComposite_MicroMechanics/Results.txt", std::fstream::in | std::fstream::out);
      /*
        file1 << "        Printing input values needed to calculate the 9*9 stiffness matrix" << std::endl;
        file1 << "            E11: " << VAMBE.CMP.E11 << std::endl;
        file1 << "            E22: " << VAMBE.CMP.E22 << std::endl;
        file1 << "            nu12: " << VAMBE.CMP.nu12 << std::endl;
        file1 << "            G12: " << VAMBE.CMP.G12 << std::endl;
        file1 << "            nu21: " << VAMBE.CMP.E22 * VAMBE.CMP.nu12 / VAMBE.CMP.E11 << std::endl;
        file1 << "            width of the laminate: " << VAMBE.CS.width << std::endl;
        file1 << "            height of laminate: " << VAMBE.CS.height << std::endl;
        file1 << "            Initial twist vector: " << VAMBE.CMP.inittwist(0) << std::endl;
        for (int j = 0; j < 5; j++)
        {
            file1 << "            " << VAMBE.CMP.Orient(j) << std::endl;
        }
        file1 << "            " << VAMBE.CMP.np << std::endl;
        */
        Stiffness S1;
        S1 = StiffnessMatrix_VAM(VAMBE.CMP.E11, VAMBE.CMP.E22, VAMBE.CMP.nu12, VAMBE.CMP.G12,
            VAMBE.CMP.E22 * VAMBE.CMP.nu12 / VAMBE.CMP.E11, VAMBE.CS.width, VAMBE.CS.height,
            VAMBE.CMP.inittwist(0), VAMBE.CMP.Orient, (int)VAMBE.CMP.np, file1);
        /*
        file1 << "                Printing linear stiffness matrix Sl" << std::endl;
        for (int i = 0; i < 4; i++)
        {
            file1 << "                    ";
            for (int j = 0; j < 4; j++)
            {
                file1 << S1.Sl(i, j) << "  ";
            }
            file1 << std::endl;
        }
        file1 << "                Printing non-linear stiffness matrix Sln" << std::endl;
        for (int i = 0; i < 4; i++)
        {
            file1 << "                    ";
            for (int j = 0; j < 5; j++)
            {
                file1 << S1.Sln(i, j) << "  ";
            }
            file1 << std::endl;
        }
        file1 << "                Printing non-linear stiffness matrix Sn" << std::endl;
        for (int i = 0; i < 5; i++)
        {
            file1 << "                    ";
            for (int j = 0; j < 5; j++)
            {
                file1 << S1.Sn(i, j) << "  ";
            }
            file1 << std::endl;
        }*/
        //Loop for force iterations.
        if (file1.is_open() && file2.is_open())
        {
            file1 << "Evaluating for force " << VAMBE.B.FN(0) << std::endl;
            do
            {
                file1 << "    Evaluating primary unkowns" << std::endl;
                file1 << "    Displacements, Rotations, Internal forces and moments" << std::endl;
                do
                {

                    //Initialize all matrices to zero
                    U_new.setZero();
                    F.setZero();
                    R.setZero();
                    J.setZero();

                    Eigen::VectorXd U1(6), F1(6), U2(6), F2(6);
                    Eigen::MatrixXd Seq_1 = Eigen::MatrixXd::Identity(6, 6);
                    Eigen::MatrixXd Seq_2 = Eigen::MatrixXd::Identity(6, 6);

                    for (int i = 0; i < VAMBE.NELEM; i++)
                    {
                        file1 << "INTERIOR                "<< std::endl;
                        int a = (int)VAMBE.ELEM(i, 1) - 1;
                        int b = (int)VAMBE.ELEM(i, 2) - 1;
                        double xa = VAMBE.NODE(a, 0);
                        double xb = VAMBE.NODE(b, 0);
                        double h = xb - xa;

                        //Dont' know what Theta means but it is used in C_ab formula
                        double Theta = 30;

                        for (int j = 0; j < 6; j++)
                        {
                            U1(j) = U(12 * a + j);
                            F1(j) = U(12 * a + 6 + j);
                            U2(j) = U(12 * b + j);
                            F2(j) = U(12 * b + 6 + j);
                        }
                        /*
                        file1 << "        Printing nodal values for node" << i << std::endl;
                        file1 << "            Nodal displacements and rotations" << std::endl;
                        for (int j = 0; j < 6; j++)
                            file1 << "            " << U1(j) << std::endl;
                        file1 << "            Internal forces and moments" << std::endl;
                        for (int j = 0; j < 6; j++)
                            file1 << "            " << F1(j) << std::endl;
                        file1 << "        Printing nodal values for node" << i + 1 << std::endl;
                        file1 << "            Nodal displacements and rotations" << std::endl;
                        for (int j = 0; j < 6; j++)
                            file1 << "            " << U2(j) << std::endl;
                        file1 << "            Internal forces and moments" << std::endl;
                        for (int j = 0; j < 6; j++)
                            file1 << "            " << F2(j) << std::endl;*/

                        //2D Cross-Section Analysis
                        //S is the cross sectional stiffness matrix obtained from the cross sectional analysis
                        //9*9 stiffness matrix is borrowed from "Non-classical effects in non-linear analysis of pretwisted
                        //anisotropic strips" by Harursampath et. al
                        //It is converted to 4*4 by taking derivative of one dimensional strain energy with respect to the strains.
                        //The matrix is in general converted to 6*6 to generalize the code such that there's 3 dofs at each node accounting
                        //for a total of 6 dofs for each element.
                        //This can be extended to Timoshenko formulation as well.
                        //Formula for 1D strain energy can be found in Eq. (31) in the paper.
                        Eigen::VectorXd Strain1(6), Strain2(6);
                        for (int k = 0; k < 6; k++)
                        {
                            Strain1(k) = Tau(6 * a + k);
                            Strain2(k) = Tau(6 * b + k);
                        }

                        Seq_1 = Equivalent_ClassicalStiffnessModel_VAM(S1, Strain1, file1);
                        Seq_2 = Equivalent_ClassicalStiffnessModel_VAM(S1, Strain2, file1);

                        file1 << "        Printing equivalent cross-sectional stiffness matrix for node " << i << std::endl;
                        for (int j = 0; j < 6; j++)
                        {
                            file1 << "            ";
                            for (int k = 0; k < 6; k++)
                            {
                                file1 <<Seq_1(j, k) << "  ";
                            }
                            file1 << std::endl;
                        }
                        file1 << "        Printing equivalent cross-sectional stiffness matrix for node " << i + 1 << std::endl;
                        for (int j = 0; j < 6; j++)
                        {
                            file1 << "            ";
                            for (int k = 0; k < 6; k++)
                            {
                                file1 << Seq_2(j, k) << "  ";
                            }
                            file1 << std::endl;
                        }
                        //1D Non Linear Beam Analysis
                        //Implemented from "GEBT: A general-purpose nonlinear analysis tool for composite beams"
                        //by Wenbin Yu et. al
                        //The Geometrically Exact Beam Theory equation can be found in Eq. 5.50 of
                        //"Non linear Composite beam theory by Hodges"
                        //This equation is simplified using Finite elements in the paper by Wenbin Yu et. al
                        //The equations required for implementation are (45), (46) and (47) from the paper.
                        //They are solved using Newton Raphson method.
                        Eigen::VectorXd Element_R(12);
                        Eigen::MatrixXd Element_J(12, 24);

                        //Element Residual
                        Element_R = Element_Residual(U1, U2, F1, F2, h, Seq_1, Seq_2, Theta);
                        //Element Jacobian
                        Element_J = Element_Jacobian(U1, U2, F1, F2, h, Seq_1, Seq_2, Theta, file1);

                        file1 << "        Element Residual Vector for node " << i + 1 << std::endl;
                        for (int j = 0; j < 12; j++)
                        {
                            file1 << "            " << Element_R(j) << std::endl;
                        }
                        file1 << "        Element Jacobian Matrix for node " << i + 1 << std::endl;
                        for (int j = 0; j < 12; j++)
                        {
                            file1 << "            ";
                            for (int k = 0; k < 24; k++)
                            {
                                file1 << Element_J(j, k) << "  ";
                            }
                            file1 << std::endl;
                        }

                        //Assembly
                        for (int j = 0; j < 12; j++)
                            for (int k = 0; k < 24; k++)
                                J.coeffRef(12 * b + j, 12 * a + k) += Element_J(j, k);
                        for (int j = 0; j < 12; j++)
                            R.coeffRef(12 * b + j) += Element_R(j);
                    }

                    //Boundary Conditions
                    //Add Jacobian and Residual for equations in first and last node
                    file1 << "BOUNDARY           " << std::endl;
                    Eigen::VectorXd Element_R(12);
                    Eigen::MatrixXd Element_J(12, 12);

                    //Node 0
                    for (int j = 0; j < 6; j++)
                    {
                        U1(j) = U(12 * 0 + j);
                        F1(j) = U(12 * 0 + 6 + j);
                    }
                    double xa = VAMBE.NODE(0, 0);
                    double xb = VAMBE.NODE(1, 0);
                    double h = xb - xa;

                    Seq_1 = Equivalent_ClassicalStiffnessModel_VAM(S1, U1, file1);

                    double Theta = 30;
                    Element_R = Element_Residual(U1, F1, h, Seq_1, Theta, 0, VAMBE.B);
                    Element_J = Element_Jacobian(U1, F1, h, Seq_1, Theta, 0, file1);

                    file1 << "        Element Residual Vector for node " << 0 + 1 << std::endl;
                    for (int j = 0; j < 12; j++)
                    {
                        file1 << "            " << Element_R(j) << std::endl;
                    }
                    file1 << "        Element Jacobian Matrix for node " << 0 + 1 << std::endl;
                    for (int j = 0; j < 12; j++)
                    {
                        file1 << "            ";
                        for (int k = 0; k < 12; k++)
                        {
                            file1 << Element_J(j, k) << "  ";
                        }
                        file1 << std::endl;
                    }

                    //Assembly of Node 0 into Global Stiffness matrix
                    for (int j = 0; j < 12; j++)
                        for (int k = 0; k < 12; k++)
                            J.coeffRef(12 * 0 + j, 12 * 0 + k) += Element_J(j, k);
                    for (int j = 0; j < 12; j++)
                        R.coeffRef(12 * 0 + j) += Element_R(j);

                    //Node N
                    //The idea is to add a dof N+1 and after solving assign the values of displacements
                    //of node N+1 to N
                    int lastnode = VAMBE.NNODE - 1;
                    for (int j = 0; j < 6; j++)
                    {
                        U1(j) = U(12 * lastnode + j);
                        F1(j) = U(12 * lastnode + 6 + j);
                    }
                    xa = VAMBE.NODE(VAMBE.NNODE - 2, 0);
                    xb = VAMBE.NODE(VAMBE.NNODE - 1, 0);
                    h = xb - xa;

                    Seq_1 = Equivalent_ClassicalStiffnessModel_VAM(S1, U1, file1);

                    Element_R = Element_Residual(U1, F1, h, Seq_1, Theta, VAMBE.NNODE - 1, VAMBE.B);
                    Element_J = Element_Jacobian(U1, F1, h, Seq_1, Theta, VAMBE.NNODE - 1, file1);

                    file1 << "        Element Residual Vector for node " << VAMBE.NNODE << std::endl;
                    for (int j = 0; j < 12; j++)
                    {
                        file1 << "            " << Element_R(j) << std::endl;
                    }
                    file1 << "        Element Jacobian Matrix for node " << VAMBE.NNODE << std::endl;
                    for (int j = 0; j < 12; j++)
                    {
                        file1 << "            ";
                        for (int k = 0; k < 12; k++)
                        {
                            file1 << Element_J(j, k) << "  ";
                        }
                        file1 << std::endl;
                    }

                    //Assembly of Node n into global stiffness matrix
                    for (int j = 0; j < 12; j++)
                        for (int k = 0; k < 12; k++)
                            J.coeffRef(12 * VAMBE.NNODE + j, 12 * VAMBE.NNODE + k) += Element_J(j, k);
                    for (int j = 0; j < 12; j++)
                        R.coeffRef(12 * VAMBE.NNODE + j) += Element_R(j);

                    file1 << "        Global Jacobian/Tangent Matrix" << std::endl;
                    file1 << J << std::endl;
                    file1 << "        Global Residual Vector" << std::endl;
                    file1 << R << std::endl;

                    //Solve the equation
                    J.makeCompressed();

                    Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>> solver;
                    solver.analyzePattern(J);
                    solver.factorize(J); //LU decomposition

                    //std::cout<<J<<std::endl;

                    assert(solver.info() == Eigen::Success);

                    dU = solver.solve(-R);

                    file1 << "            Change in displacements after " << iter << " iteration" << std::endl;
                    for (int j = 0; j < (VAMBE.NNODE + 1)  * VAMBE.NDOF; j++)
                    {
                        file1 << "                " << dU(j) << std::endl;
                    }

                    //Calculate errors and update for next iteration
                    for (int j = 0; j < (VAMBE.NNODE + 1) * VAMBE.NDOF; j++)
                        U_new(j) = U(j) + dU(j);

                    for (int j = 0; j < 12; j++)
                    {
                        U_new(12 * (VAMBE.NNODE - 1) + j) = U_new(12 * VAMBE.NNODE + j);
                    }

                    file1 << "            Corrected Displacements After " << iter << " iteration" << std::endl;
                    for (int j = 0; j < (VAMBE.NNODE + 1) * VAMBE.NDOF; j++)
                    {
                        file1 << "                " << U_new(j) << std::endl;
                    }
                   
                    //Error calculation
                    for (int j = 0; j < (VAMBE.NNODE + 1) * VAMBE.NDOF; j++)
                        error(j) = abs(U_new(j) - U(j));

                    file1 << "            Error " << iter << " iteration" << std::endl;
                    for (int j = 0; j < (VAMBE.NNODE + 1) * VAMBE.NDOF; j++)
                    {
                        file1 << "                " << error(j) << std::endl;
                    }

                    max = error(0);
                    for (int j = 0; j < (VAMBE.NNODE + 1) * VAMBE.NDOF; j++)
                        if (max < error(j))
                            max = error(j);

                    file1 << "            Maximum Error " << std::endl;
                    file1 << "                " << max << std::endl;
                    //Assignment for next iteration
                    for (int j = 0; j < (VAMBE.NNODE + 1) * VAMBE.NDOF; j++)
                        U(j) = U_new(j);
                    iter++;

                    //After calculating the internal forces through 1D nonlinear beam analysis,
                    //the strains are updated using the constitutive relation.
                    Update_Strains(Tau, VAMBE.NNODE, S1, U, file1);

                    file1 << "            Updates Strains After "<< iter << " iteration" << std::endl;
                    for (int j = 0; j < VAMBE.NNODE * 6; j++)
                    {
                        file1 << "                " << Tau(j) << std::endl;
                    }

                } while (max > pow(10, -6) && iter < maxiter);
                if (iter < maxiter)
                {
                    fiter++;
                    file2 << VAMBE.B.FN(0) << "   " << U(12 * (VAMBE.NNODE - 1) + 0);
                    file2 << "Converged Displacements for load " << VAMBE.B.FN(0) << std::endl;
                    for (int j = 0; j < (VAMBE.NNODE + 1) * VAMBE.NDOF; j++)
                        file2 << U(j) << std::endl;
                    VAMBE.B.FN(0) = VAMBE.B.FN(0) + 200.0;
                }
                else 
                    file2 << "Code didn't converge for force   " << VAMBE.B.FN(0) << std::endl;
            } while (VAMBE.NLS > fiter);
        }
        else
        {
            if(!file1.is_open())
                std::cout << "Couldn't open file Result_Log.txt" << std::endl;
            if(!file2.is_open())
                std::cout << "Couldn't open file Results.txt" << std::endl;
        }
        file1.close();
        file2.close();
    }
    //--------------------------------------------------------//
    //--------------------------------------------------------//
    //----------3D Corotational Formulation------------------//
    //-----------------Euler Bernoulli-----------------------//
    //-------------------------------------------------------//
    //-------------------------------------------------------//
    else if (choice == 5)
    {
        

    }
    else
        std::cout << "Wrong option" << std::endl;
}

#endif




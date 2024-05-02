
#pragma once
#include "Variables.h"
#include<fstream>
#include<limits>
#include<chrono>

int main()
{
    auto start = std::chrono::high_resolution_clock::now();
    // 1 - Linear Bar Element
    // 2 - Nonlinear Bar Element
    // 3 - Nonlinear Euler Bernouli Beam Element (3 dofs per node)
    // 4 - GEBT(Geometrically Exact Beam theory) based beam element
    // 5 - Nonlinear Euler Bernouli Beam Element (6 dofs per node)(quadratic)
    // 6 - Nonlinear Euler Bernouli Beam Element with contact using Node to Node (n bodies)
    // 7 - Nonlinear Euler Bernouli Beam Element with contact using Segment to Segment (2 bodies)
    // 8 - GEBT(Geometrically Exact Beam theory) based beam element with contact using Node to Node (n bodies)
    int choice = 6;
    //LoadVector(M.LOAD,M.ELSET,LOAD,M.NNODE);
    // -------------------------------------//
    //--------------------------------------//
    //----------BAR ELEMENT LINEAR----------//
    //--------------------------------------//
    //--------------------------------------//
    //This code works only for linear bar element with a concentrated loads.
    if (choice == 1)
    {
        LinearBarElement LBE = ReadLBEFile();

        Eigen::MatrixXd k = Eigen::MatrixXd::Zero(2, 2);
        Eigen::VectorXd f = Eigen::VectorXd::Zero(2);
        Eigen::SparseMatrix<double, Eigen::ColMajor> K(LBE.NNODE, LBE.NNODE);
        Eigen::VectorXd F = Eigen::VectorXd::Zero(LBE.NNODE);
        Eigen::VectorXd U = Eigen::VectorXd::Zero(LBE.NNODE);
        double height = LBE.CS.Rect.height;
        double width = LBE.CS.Rect.width;
        double Area = height * width;
        int count = 0;
        for (int i = 0; i < (int)LBE.NELEM; i++)
        {
            double E = LBE.MAT((int)LBE.ELEM(i, 0) - 1, 0);
            int a = (int)LBE.ELEM(i, 1) - 1;
            int b = (int)LBE.ELEM(i, 2) - 1;
            double xa = LBE.NODE((int)a, 0);
            double xb = LBE.NODE((int)b, 0);
            StiffnessMatrix_LBE(Area, E, xa, xb, k);

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
            K.coeffRef(a, a) += k(0, 0);
            K.coeffRef(a, b) += k(0, 1);
            K.coeffRef(b, a) += k(1, 0);
            K.coeffRef(b, b) += k(1, 1);

            //Assembly of Global Force Vector
            F(a) += f(0);
            F(b) += f(1);
        }

        //Apply constraints
        ApplyConstraints_LBE(LBE.CNODE, LBE.NNODE, F, K);

        std::cout << "Stiffness Matrix" << K << std::endl;
        std::cout << "Force Vector" << F << std::endl;

        //Solve for displacements
        Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>> solver;
        solver.analyzePattern(K);
        solver.factorize(K); //LU decomposition

        assert(solver.info() == Eigen::Success);

        U = solver.solve(F);
        std::cout << "Displacements" << std::endl;
        std::cout << U << std::endl;
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
    else if (choice == 2)
    {
        NonLinearBarElement NLBE = ReadNLBEFile();

        Eigen::VectorXd dU = Eigen::VectorXd::Zero(NLBE.NDOF * NLBE.NNODE);
        Eigen::VectorXd U = Eigen::VectorXd::Zero(NLBE.NDOF * NLBE.NNODE);
        Eigen::VectorXd error = Eigen::VectorXd::Zero(NLBE.NDOF * NLBE.NNODE);
        Eigen::MatrixXd k = Eigen::MatrixXd::Zero(2, 2);
        Eigen::MatrixXd t = Eigen::MatrixXd::Zero(2, 2);
        Eigen::VectorXd f = Eigen::VectorXd::Zero(2);
        Eigen::VectorXd F = Eigen::VectorXd::Zero(NLBE.NDOF * NLBE.NNODE);
        Eigen::VectorXd U_new = Eigen::VectorXd::Zero(NLBE.NDOF * NLBE.NNODE);
        Eigen::VectorXd R = Eigen::VectorXd::Zero(NLBE.NDOF * NLBE.NNODE);
        Eigen::SparseMatrix<double, Eigen::RowMajor> K(NLBE.NNODE * NLBE.NDOF, NLBE.NNODE * NLBE.NDOF),
            T(NLBE.NNODE * NLBE.NDOF, NLBE.NNODE * NLBE.NDOF);
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
            for (int i = 0; i < (int)NLBE.NELEM; i++)
            {
                int a = (int)NLBE.ELEM(i, 0) - 1;
                int b = (int)NLBE.ELEM(i, 1) - 1;
                double xa = NLBE.NODE(a, 0);
                double xb = NLBE.NODE(b, 0);
                double h = xb - xa;

                StiffnessMatrix_NLBE(k, U, a, b, h);
                K.coeffRef(a, a) += k(0, 0);
                K.coeffRef(a, b) += k(0, 1);
                K.coeffRef(b, a) += k(1, 0);
                K.coeffRef(b, b) += k(1, 1);

                TangentStiffnessMatrix_NLBE(t, U, a, b, h);
                T.coeffRef(a, a) += t(0, 0);
                T.coeffRef(a, b) += t(0, 1);
                T.coeffRef(b, a) += t(1, 0);
                T.coeffRef(b, b) += t(1, 1);

                LocalForceVec_NLBE(f, a, b, h);
                F.coeffRef(a) += f(0);
                F.coeffRef(b) += f(1);
            }


            //Residual Matrix
            R = K * U - F;

            //Apply Constraints
            ApplyConstraints_NLBE(T, U, R, NLBE.NNODE, NLBE.CNODE);
            K.makeCompressed();
            T.makeCompressed();

            std::cout << "Tangent Stiffness Matrix" << std::endl << T << std::endl;

            std::cout << "Stiffness Matrix" << std::endl << K << std::endl;

            std::cout << "Residual" << std::endl << R << std::endl;

            //Matrix Solution
            Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::RowMajor>> solver;
            solver.analyzePattern(T);
            solver.factorize(T); //LU decomposition

            assert(solver.info() == Eigen::Success);

            dU = solver.solve(-R);

            //            std::cout<<"Increment in displacement"<<std::endl<<dU<<std::endl;

            //Next iteration
            for (int i = 0; i < NLBE.NNODE * NLBE.NDOF; i++)
                U_new(i) = U(i) + dU(i);

            std::cout << "Displacements" << U_new << std::endl;
            //Error calculation
            for (int i = 0; i < NLBE.NNODE * NLBE.NDOF; i++)
                error(i) = abs((U_new(i) - U(i)) / U(i));
            max = error(0);
            for (int i = 0; i < NLBE.NNODE * NLBE.NDOF; i++)
                if (max < error(i))
                    max = error(i);

            //Assignment for next iteration
            for (int i = 0; i < NLBE.NNODE * NLBE.NDOF; i++)
                U(i) = U_new(i);

            iter++;

            //std::cout<<"Displacements"<<U<<std::endl;
        } while (max > pow(10, -3) && iter < maxiter);
    }
    //----------------------------------------------
    //----------------------------------------------
    //---------EULER BERNOULLI NONLINEAR 2D---------
    //----------------------------------------------
    //----------------------------------------------
    //The problem description for which the Euler Bernouli beam code works is a unit distributed load on a cantilever beam
    //To know more look at Chapter 4 of "Introduction to Nonlinear FEM by JN Reddy"
    //Default parameters are taken to solve Example 4.2.1. 
    else if (choice == 3)
    {
        const int nbeams = 1;

        NLEBBE2D EBBE2D = NLEBBE2D::NLEBBE2D();

        Eigen::VectorXd dU = Eigen::VectorXd::Zero(EBBE2D.get_ndof() * EBBE2D.get_nnode());
        Eigen::VectorXd U = Eigen::VectorXd::Zero(EBBE2D.get_ndof() * EBBE2D.get_nnode());
        Eigen::VectorXd error = Eigen::VectorXd::Zero(EBBE2D.get_ndof() * EBBE2D.get_nnode());
        Eigen::MatrixXd k = Eigen::MatrixXd::Zero(6, 6);
        Eigen::MatrixXd t = Eigen::MatrixXd::Zero(6, 6);
        Eigen::VectorXd f = Eigen::VectorXd::Zero(6);
        Eigen::VectorXd F = Eigen::VectorXd::Zero(EBBE2D.get_ndof() * EBBE2D.get_nnode());
        Eigen::VectorXd U_new = Eigen::VectorXd::Zero(EBBE2D.get_ndof() * EBBE2D.get_nnode());
        Eigen::VectorXd R = Eigen::VectorXd::Zero(EBBE2D.get_ndof() * EBBE2D.get_nnode());
        Eigen::SparseMatrix<double, Eigen::ColMajor> K(EBBE2D.get_nnode() * EBBE2D.get_ndof(), EBBE2D.get_nnode() * EBBE2D.get_ndof()),
            T(EBBE2D.get_nnode() * EBBE2D.get_ndof(), EBBE2D.get_nnode() * EBBE2D.get_ndof());
        //Eigen::VectorXd X_ref = Eigen::VectorXd::Zero(NLEBBE.NDM * EBBE2D.get_nnode());
        //Eigen::VectorXd X_curr = Eigen::VectorXd::Zero(NLEBBE.NDM * EBBE2D.get_nnode());
     
        //Eigen::MatrixXd LOAD;
        double max = 1;
        int iter = 1;
        int fiter = 1;
        int maxiter = 50;
        //Initialize Reference and current configuration
        //for (int i = 0; i < (int)EBBE2D.get_nnode(); i++)
        //{
        //    X_ref(NLEBBE.NDM * i) = EBBE2D.get_coordinates(i, 0);
        //    X_curr(NLEBBE.NDM* i) = EBBE2D.get_coordinates(i, 0);
        //    X_ref(NLEBBE.NDM * i + 1) = EBBE2D.get_coordinates(i, 1);
        //    X_curr(NLEBBE.NDM* i + 1) = EBBE2D.get_coordinates(i, 1);
        //}
        //Loop for force iterations.
        do
        {
            //Newton Raphson loop
            //Update Reference configuration
            //for (int i = 0; i < (int)EBBE2D.get_nnode(); i++)
            //{
            //    X_ref(NLEBBE.NDM * i) += U(EBBE2D.get_ndof() * i);
            //    X_ref(NLEBBE.NDM * i + 1) += U(EBBE2D.get_ndof() * i + 1);
            //}
            do
            {
                //Initialize all matrices to zero
                U_new.setZero();
                F.setZero();
                R.setZero();
                K.setZero();
                T.setZero();
                //f.setZero();
                

                for (int i = 0; i < (int)EBBE2D.get_nelem(); i++)
                {
                    int a = (int)EBBE2D.get_connectivity(i, 1) - 1;
                    int b = (int)EBBE2D.get_connectivity(i, 2) - 1;
                    double xa = EBBE2D.get_coordinates((int)a, 0);
                    double xb = EBBE2D.get_coordinates((int)b, 0);
                    double E = EBBE2D.get_modelprop("E");
                    double nu = EBBE2D.get_modelprop("nu");

                    double A = EBBE2D.get_modelprop("b") * EBBE2D.get_modelprop("h");

                    double I = EBBE2D.get_modelprop("b") * pow(EBBE2D.get_modelprop("h"), 3) / 12;

                    double vf = EBBE2D.get_loadprop("vf");
                    double af = EBBE2D.get_loadprop("af");

                    //Local stiffness matrix
                    k = EBBE2D.StiffnessMatrix_NLEBBE(xa, xb, E, nu, A, I, U, a, b);

                    //Tangent Stiffness matrix
                    t = EBBE2D.TangentStiffnessMatrix_NLEBBE(k, xa, xb, E, nu, A, I, U, a, b);

                    //Local Force Vector
                    f = EBBE2D.LocalFoceVec_NLEBBE(xa, xb, a, b, vf, af, U, E, nu);

                    EBBE2D.RearrangeElementStiffness_NLEBBE(k, t, f);

                    //std::cout << k << std::endl;

                    //std::cout << t << std::endl;

                    //std::cout << f << std::endl;

                    //Assembly stiffness matrix
                    for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                        {
                            K.coeffRef(3 * a + i, 3 * a + j) += k(i, j);
                            K.coeffRef(3 * b + i, 3 * b + j) += k(i + 3, j + 3);
                            K.coeffRef(3 * a + i, 3 * b + j) += k(i, j + 3);
                            K.coeffRef(3 * b + i, 3 * a + j) += k(i + 3, j);
                        }

                    //Assembly Tangent stiffness matrix
                    for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                        {
                            T.coeffRef(3 * a + i, 3 * a + j) += t(i, j);
                            T.coeffRef(3 * b + i, 3 * b + j) += t(i + 3, j + 3);
                            T.coeffRef(3 * a + i, 3 * b + j) += t(i, j + 3);
                            T.coeffRef(3 * b + i, 3 * a + j) += t(i + 3, j);
                        }

                    //Assembly Force Vector
                    for (int i = 0; i < 3; i++)
                    {
                        F(3 * a + i) += f(i);
                        F(3 * b + i) += f(i + 3);
                    }
                }

                //std::cout << K << std::endl;

                //Residual Matrix
                R = K * U - F;

               //Apply Constraints
                EBBE2D.ApplyConstraints_NLEBBE(T, U, EBBE2D.get_nnode(), R);

                K.makeCompressed();
                T.makeCompressed();

                //Matrix Solution
                Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>> solver;

                solver.analyzePattern(T);
                solver.factorize(T); //LU decomposition

                assert(solver.info()==Eigen::Success);

                dU = solver.solve(-R);

                //std::cout << dU << std::endl;

                //Next iteration
                for(int i=0;i<EBBE2D.get_nnode()*EBBE2D.get_ndof();i++)
                    U_new(i) = U(i) + dU(i);

                //Error calculation
                for (int i = 0; i < EBBE2D.get_nnode() * EBBE2D.get_ndof(); i++)
                    error(i) = abs(U_new(i)-U(i));
                max = error(0);
                for (int i = 0; i < EBBE2D.get_nnode() * EBBE2D.get_ndof(); i++)
                    if (max < error(i))
                        max = error(i);

                //Assignment for next iteration
                for (int i = 0; i < EBBE2D.get_nnode() * EBBE2D.get_ndof(); i++)
                    U(i) = U_new(i);
                iter++;
                
            } while (max > pow(10, -3) && iter < maxiter);
            if (fiter < EBBE2D.get_nls())
            {
                VTKGrid VTK1 = VTKGrid(EBBE2D, EBBE2D.get_loadprop("vf"), "EBBE2D", U, 2, 9, EBBE2D.get_ndof(), 4);
                fiter++;
                std::cout << EBBE2D.get_loadprop("vf") << "  " << U(((int)EBBE2D.get_nnode() / 2) * EBBE2D.get_ndof() + 1) << std::endl;
                //std::cout << ((int)EBBE2D.get_nnode() / 2) * EBBE2D.get_ndof() + 2;
            }
            if (iter == maxiter)
            {
                EBBE2D.set_loadprop("vf", EBBE2D.get_loadprop("vf") - 0.5);
                std::cout << "Solution not converging for the given load" << std::endl;
                std::cout << "Trying with a load with less increment" << std::endl;
            }
            //std::cout << "Load  " << EBBE2D.get_loadprop("vf") << std::endl;
            //std::cout << "Displacements  " << std::endl << U << std::endl;

            EBBE2D.set_loadprop("vf", EBBE2D.get_loadprop("vf") + 1);
            //Load step increment
        } while (fiter < EBBE2D.get_nls());

        if (iter == maxiter)
            std::cout << "Maximum iteration limit reached" << std::endl << "Don't kid yourself" << std::endl << "Solution did not converge" << std::endl;
    }

    //----------------------------------------------
    //----------------------------------------------
    //-------------VAM Beam Element 3D--------------
    //----------------------------------------------
    //----------------------------------------------
    else if (choice == 4)
    {
        VAMBeamElement VAMBE = VAMBeamElement();

        Eigen::VectorXd dU = Eigen::VectorXd::Zero(VAMBE.get_ndof() * VAMBE.get_nnode());
        Eigen::VectorXd U = Eigen::VectorXd::Zero(VAMBE.get_ndof() * VAMBE.get_nnode());
        Eigen::VectorXd error = Eigen::VectorXd::Zero(VAMBE.get_ndof() * VAMBE.get_nnode());
        Eigen::VectorXd U_new = Eigen::VectorXd::Zero(VAMBE.get_ndof() * VAMBE.get_nnode());
        Eigen::VectorXd R = Eigen::VectorXd::Zero(VAMBE.get_ndof() * VAMBE.get_nnode());
        Eigen::SparseMatrix<double, Eigen::ColMajor> J(VAMBE.get_nnode() * VAMBE.get_ndof(), VAMBE.get_nnode() * VAMBE.get_ndof());
        Eigen::VectorXd Tau = Eigen::VectorXd::Zero(6 * VAMBE.get_nnode());
        double max;
        int iter = 1;
        int fiter = 0;
        int maxiter = 1000;
        std::fstream file1, file2;
        file1.open("E:/Adhithya/MTech_Res_Thesis/Cpp/ThesisCode/TextileComposites/TextileComposites/Result_Log.txt", std::fstream::in | std::fstream::out);
        file2.open("E:/Adhithya/MTech_Res_Thesis/Cpp/ThesisCode/TextileComposites/TextileComposites/Results.txt", std::fstream::in | std::fstream::out);

        //Don't know what Theta means but it is used in C_ab formula
        double Theta = 0 * M_PI / 180;
        Eigen::VectorXd aN = VAMBE.get_load(3);
        double load = aN(2);
        //Loop for force iterations.
        if (file1.is_open() && file2.is_open())
        {
            do
            {
                iter = 0;
                VAMBE.set_load(3, 3, fiter * 400/ VAMBE.get_nls());
                //std::cout << VAMBE.get_load(3) << std::endl;
                do
                {
                    //Initialize all matrices to zero
                    U_new.setZero();
                    R.setZero();
                    J.setZero();
                    Tau.setZero();

                    Eigen::VectorXd U1(6), F1(6), U2(6), F2(6);
                    Eigen::MatrixXd Seq_1 = Eigen::MatrixXd::Zero(6, 6);
                    Eigen::MatrixXd Seq_2 = Eigen::MatrixXd::Zero(6, 6);

                    for (int i = 0; i < VAMBE.get_nelem(); i++)
                    {
                        int a = (int)VAMBE.get_connectivity(i, 1) - 1;
                        int b = (int)VAMBE.get_connectivity(i, 2) - 1;
                        double xa = VAMBE.get_coordinates(a, 0);
                        double xb = VAMBE.get_coordinates(b, 0);
                        double h = xb - xa;
                        if (i < VAMBE.get_nelem() - 1)
                        {
                            if (i == 0)
                            {
                                Eigen::VectorXd Element_R = Eigen::VectorXd::Zero(12);
                                Eigen::MatrixXd Element_J = Eigen::MatrixXd::Zero(12, 18);

                                Eigen::VectorXd U0(6);
                                for (int j = 0; j < 6; j++)
                                {
                                    U1(j) = U(a * 12 + 6 + j);
                                    F1(j) = U(a * 12 + 12 + j);
                                    //Fhat
                                    //Mhat
                                    U0(j) = U(j);
                                }

                                Eigen::VectorXd Strain(6);
                                for (int j = 0; j < 6; j++)
                                {
                                    Strain(j) = Tau(6 * a + j);
                                }
                                Seq_1 = VAMBE.Equivalent_StiffnessMatrix_FirstOrder(file1);

                                Element_R = VAMBE.Element_Residual(U1, F1, h, Seq_1, Theta, 0, VAMBE.get_load(1), VAMBE.get_load(2), VAMBE.get_load(3), VAMBE.get_load(4), U0);
                                Element_J = VAMBE.Element_Jacobian(U1, F1, h, Seq_1, Theta, 0, file1);

                                //Assembly of Node 0 into Global Stiffness matrix
                                //fu1 - F1 = 0
                                //fpsi1 - M1 = 0
                                //fF1 - u1 = 0
                                //fM1 - theta1 = 0
                                for (int j = 0; j < 12; j++)
                                {
                                    R.coeffRef(j) += Element_R(j);
                                    for (int k = 0; k < 18; k++)
                                    {
                                        J.coeffRef(j, k) += Element_J(j, k);
                                    }
                                }
                            }

                            for (int j = 0; j < 6; j++)
                            {
                                U1(j) = U(a * 12 + 6 + j);
                                F1(j) = U(a * 12 + 12 + j);
                                U2(j) = U(b * 12 + 6 + j);
                                F2(j) = U(b * 12 + 12 + j);
                            }

                            Eigen::VectorXd Strain1(6), Strain2(6);
                            for (int j = 0; j < 6; j++)
                            {
                                Strain1(j) = Tau(6 * a + j);
                                Strain2(j) = Tau(6 * b + j);
                            }

                            Seq_1 = VAMBE.Equivalent_StiffnessMatrix_FirstOrder(file1);
                            Seq_2 = VAMBE.Equivalent_StiffnessMatrix_FirstOrder(file1);

                            Eigen::VectorXd Element_R = Eigen::VectorXd::Zero(12);
                            Eigen::MatrixXd Element_J = Eigen::MatrixXd::Zero(12, 24);

                            //Element Residual
                            Element_R = VAMBE.Element_Residual(U1, U2, F1, F2, h, Seq_1, Seq_2, Theta);
                            //Element Jacobian
                            Element_J = VAMBE.Element_Jacobian(U1, U2, F1, F2, h, Seq_1, Seq_2, Theta, file1);

                            //Assembly
                            //Solve for intermediate nodes
                            //fui + fui+1
                            //fpsii + fpsii+1
                            //fFi + fFi+1
                            //fMi + fMi+1
                            for (int j = 0; j < 12; j++)
                            {
                                R.coeffRef(12 * b + j) += Element_R(j);
                                for (int k = 0; k < 24; k++)
                                    J.coeffRef(12 * b + j, 12 * a + 6 + k) += Element_J(j, k);
                            }
                        }
                        else if (i == VAMBE.get_nelem() - 1)
                        {
                            Eigen::VectorXd Element_R = Eigen::VectorXd::Zero(12);
                            Eigen::MatrixXd Element_J = Eigen::MatrixXd::Zero(12, 18);

                            Eigen::VectorXd U0(6);
                            for (int j = 0; j < 6; j++)
                            {
                                U1(j) = U(12 * a + 6 + j);
                                F1(j) = U(12 * a + 12 + j);
                                //FN+1
                                //MN+1
                                U0(j) = U(12 * b + 6 + j);
                            }

                            Eigen::VectorXd Strain = Eigen::VectorXd::Zero(6);
                            for (int j = 0; j < 6; j++)
                            {
                                Strain(j) = Tau(6 * a + j);
                            }

                            Seq_1 = VAMBE.Equivalent_StiffnessMatrix_FirstOrder(file1);

                            Element_R = VAMBE.Element_Residual(U1, F1, h, Seq_1, Theta, VAMBE.get_nnode() - 1, VAMBE.get_load(1), VAMBE.get_load(2), VAMBE.get_load(3), VAMBE.get_load(4), U0);
                            Element_J = VAMBE.Element_Jacobian(U1, F1, h, Seq_1, Theta, VAMBE.get_nnode() - 1, file1);

                            //fuN - FN+1 = 0
                            //fpsiN - MN+1 = 0
                            //fFN + uN+1 = 0
                            //fMN + thetaN+1 = 0
                            for (int j = 0; j < 12; j++)
                            {
                                R.coeffRef(12 * b + j) += Element_R(j);
                                for (int k = 0; k < 18; k++)
                                    J.coeffRef(12 * b + j, 12 * a + 6 + k) += Element_J(j, k);
                            }
                        }
                    }

                    //Solve the equation
                    J.makeCompressed();

                    Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>> solver;
                    solver.analyzePattern(J);
                    solver.factorize(J); //LU decomposition

                    assert(solver.info() == Eigen::Success);

                    dU = solver.solve(-R);

                    //std::cout << dU << std::endl;
                    //double alpha = 0;
                    //Calculate errors and update for next iteration
                    for (int j = 0; j < VAMBE.get_ndof() * VAMBE.get_nnode(); j++)
                        U_new(j) = U(j) + dU(j);

                    //Error calculation
                    for (int j = 0; j < VAMBE.get_ndof() * VAMBE.get_nnode(); j++)
                        error(j) = abs(U_new(j) - U(j));

                    max = error(0);
                    for (int j = 0; j < VAMBE.get_ndof() * VAMBE.get_nnode(); j++)
                        if (max < error(j))
                            max = error(j);

                    //Assignment for next iteration
                    double alpha = 0;
                    for (int j = 0; j < VAMBE.get_ndof() * VAMBE.get_nnode(); j++)
                        U(j) = alpha * U(j) + (1 - alpha) * U_new(j);
                    iter++;

                    //After calculating the internal forces through 1D nonlinear beam analysis,
                    //the strains are updated using the constitutive relation.
                    /*Tau = VAMBE.Update_Strains(VAMBE, &U, file1);

                    file1 << "            Updates Strains After " << iter << " iteration" << std::endl;
                    for (int j = 0; j < VAMBE.get_nnode() * 6; j++)
                    {
                        file1 << "                " << Tau(j) << std::endl;
                    }*/
                    //std::cout << max << std::endl;
                } while (max > pow(10, -6) && iter < maxiter);
                if (iter < maxiter)
                {
                    fiter++;
                    file2 << VAMBE.get_load(3)[2] / pow(10, 6) << "     " << U(VAMBE.get_nnode() * 12 - 4) << "    " << U(VAMBE.get_nnode() * 12 - 6) << "   " << U(VAMBE.get_nnode() * 12 - 2) << std::endl;
                    //std::cout << VAMBE.get_load(3) << "     " << U(12 * (VAMBE.get_nnode() - 2) + 6 + 2) << "    " << iter << std::endl;
                    //file2 << "Converged Displacements for load " << VAMBE.get_load(3) / pow(10, 6) << std::endl;
                    for (int j = 0; j < VAMBE.get_nnode() * VAMBE.get_ndof(); j++)
                        std::cout << U(j) << std::endl;
                }
                else
                {
                    std::cout << "Code didn't converge for force   " << VAMBE.get_load(3) << std::endl;
                    break;
                }
            } while (fiter <= VAMBE.get_nls());
        }
        else
        {
            if (!file1.is_open())
                std::cout << "Couldn't open file Result_Log.txt" << std::endl;
            if (!file2.is_open())
                std::cout << "Couldn't open file Results.txt" << std::endl;
        }
        file1.close();
        file2.close();
    }
    //--------------------------------------------------------//
    //----------3D EulerBernouli Beam Element-----------------//
    //-----------------Large displacement---------------------//
    //-----------------Updated Lagrangian---------------------//
    //---------------------Quadratic--------------------------//
    //--------------------------------------------------------//
    else if (choice == 5)
    {
        //NLEBBE3D EBBE3D = EBBE3D.ReadINP(1)
        NLEBBE3D EBBE3D = NLEBBE3D("LINELAS");

        Eigen::VectorXd dGU = Eigen::VectorXd::Zero(EBBE3D.get_ndof() * EBBE3D.get_nnode());
        Eigen::VectorXd GU = Eigen::VectorXd::Zero(EBBE3D.get_ndof() * EBBE3D.get_nnode());
        Eigen::VectorXd error = Eigen::VectorXd::Zero(EBBE3D.get_ndof() * EBBE3D.get_nnode());
        Eigen::VectorXd GU_new = Eigen::VectorXd::Zero(EBBE3D.get_ndof() * EBBE3D.get_nnode());
        Eigen::VectorXd GR = Eigen::VectorXd::Zero(EBBE3D.get_ndof() * EBBE3D.get_nnode());
        Eigen::SparseMatrix<double, Eigen::ColMajor> GT(EBBE3D.get_ndof()* EBBE3D.get_nnode(), EBBE3D.get_ndof()* EBBE3D.get_nnode());

        double max;
        int iter = 1;
        int fiter = 1;
        int maxiter = 100;

        std::fstream file1, file2;
        file1.open("E:/Adhithya/MTech_Res_Thesis/Cpp/ThesisCode/TextileComposites/TextileComposites/Result_Log.txt", std::fstream::in | std::fstream::out);
        file2.open("E:/Adhithya/MTech_Res_Thesis/Cpp/ThesisCode/TextileComposites/TextileComposites/Results.txt", std::fstream::in | std::fstream::out);
    
        //declare and initialize material parameters and other input variables
        double* D = new double[7];
        D[0] = EBBE3D.get_matprop("E");
        D[1] = EBBE3D.get_matprop("nu");
        D[2] = EBBE3D.get_matprop("Bp");
        D[3] = EBBE3D.get_matprop("Hp");
        D[4] = EBBE3D.get_matprop("Zx");
        D[5] = EBBE3D.get_matprop("Zy");
        D[6] = EBBE3D.get_matprop("Zz");
   
        /*double* D = new double[14];
        D[0] = EBBE3D.get_matprop("E1");
        D[1] = EBBE3D.get_matprop("E2");
        D[2] = EBBE3D.get_matprop("E3");
        D[3] = EBBE3D.get_matprop("nu12");
        D[4] = EBBE3D.get_matprop("nu13");
        D[5] = EBBE3D.get_matprop("nu23");
        D[6] = EBBE3D.get_matprop("G12");
        D[7] = EBBE3D.get_matprop("G13");
        D[8] = EBBE3D.get_matprop("G23");
        D[9] = EBBE3D.get_matprop("B");
        D[10] = EBBE3D.get_matprop("H");
        D[11] = EBBE3D.get_matprop("Zx");
        D[12] = EBBE3D.get_matprop("Zy");
        D[13] = EBBE3D.get_matprop("Zz");*/

        /*double* D = new double[6];
        D[0] = EBBE3D.get_matprop("C10");
        D[1] = EBBE3D.get_matprop("Bp");
        D[2] = EBBE3D.get_matprop("Hp");
        D[3] = EBBE3D.get_matprop("Zx");
        D[4] = EBBE3D.get_matprop("Zy");
        D[5] = EBBE3D.get_matprop("Zz");*/

        /*double* D = new double[7];
        D[0] = EBBE3D.get_matprop("C10");
        D[1] = EBBE3D.get_matprop("C01");
        D[2] = EBBE3D.get_matprop("Bp");
        D[3] = EBBE3D.get_matprop("Hp");
        D[4] = EBBE3D.get_matprop("Zx");
        D[5] = EBBE3D.get_matprop("Zy");
        D[6] = EBBE3D.get_matprop("Zz");*/

        /*for (int i = 0; i < 7; i++)
            std::cout << D << std::endl;*/
        double load = 0;

        if (file1.is_open() && file2.is_open())
        {
            do
            {
                iter = 0;
                do
                {
                    //Initialize all matrices to zero
                    GU_new.setZero();
                    GR.setZero();
                    GT.setZero();
                    
                    //declare local tangent stiffness matrix and residual vector
                    double** T;
                    T = new double* [EBBE3D.get_ndof() * EBBE3D.get_nen()];
                    for (int i = 0; i < EBBE3D.get_ndof() * EBBE3D.get_nen(); i++)
                        T[i] = new double[EBBE3D.get_ndof() * EBBE3D.get_nen()];

                    double* R = new double[EBBE3D.get_ndof() * EBBE3D.get_nen()];

                    for (int i = 0; i < EBBE3D.get_nelem(); i++)
                    {
                        int a = (int)EBBE3D.get_connectivity(i, 1) - 1;
                        int b = (int)EBBE3D.get_connectivity(i, 2) - 1;
                        int c = (int)EBBE3D.get_connectivity(i, 3) - 1;

                        //declare and initialize position vectors
                        /*double* X[2];
                        for (int j = 0; j < 2; j++)
                            X[j] = new double[3];*/
                        double X[3][3];

                        for (int j = 0; j < 3; j++)
                        {
                            X[0][j] = EBBE3D.get_coordinates(a, j);
                            X[1][j] = EBBE3D.get_coordinates(b, j);
                            X[2][j] = EBBE3D.get_coordinates(c, j);
                        }
                        //double h = X[2][0] - X[0][0];

                        //declare and initialize dofs
                        /*double* U[2];
                        for (int j = 0; j < 2; j++)
                            U[j] = new double[6];*/
                        double U[3][6];

                        for (int j = 0; j < 6; j++)
                        {
                            U[0][j] = GU(6 * a + j);
                            U[1][j] = GU(6 * b + j);
                            U[2][j] = GU(6 * c + j);
                        }

                        for (int j = 0; j < EBBE3D.get_ndof() * EBBE3D.get_nen(); j++)
                        {
                            R[j] = 0;
                            for (int k = 0; k < EBBE3D.get_ndof() * EBBE3D.get_nen(); k++)
                                T[j][k] = 0;
                        }

                        if (EBBE3D.get_matmodel() == "LINELAS")
                            EBBE3D.RKt(D, X, U, T, R);
                        else if (EBBE3D.get_matmodel() == "COMP_NEOHOOKE")
                            EBBE3D.RKt_NHC(D, X, U, T, R);
                        else if (EBBE3D.get_matmodel() == "INCOMP_NEOHOOKE")
                            EBBE3D.RKt_NH(D, X, U, T, R);
                        else if (EBBE3D.get_matmodel() == "MOONEY RIVLIN")
                            EBBE3D.RKt_MooneyRivlin(D, X, U, T, R);
                        else if (EBBE3D.get_matmodel() == "ORTHO")
                            EBBE3D.RKt_Orthotropic(D, X, U, T, R);
                        
                        /*for (int i = 0; i < 18; i++)
                            std::cout << R[i] << std::endl;

                        for (int j = 0; j < 18; j++)
                        {
                            for (int k = 0; k < 18; k++)
                                std::cout << T[j][k] << "   ";
                            std::cout << std::endl;
                        }*/

                        //Assembly
                        for (int j = 0; j < EBBE3D.get_ndof() * EBBE3D.get_nen(); j++)
                        {
                            GR(EBBE3D.get_ndof() * a + j) += R[j];
                            for (int k = 0; k < EBBE3D.get_ndof() * EBBE3D.get_nen(); k++)
                                GT.coeffRef(EBBE3D.get_ndof() * a + j, EBBE3D.get_ndof() * a + k) += T[j][k];
                        }
                    }

                    Loading LOAD = Loading();

                    LOAD.SetLoad(&GR, &EBBE3D, fiter, EBBE3D.get_nls());

                    Boundary BOUND = Boundary();

                    BOUND.SetBoundary(&GT, &GR, &EBBE3D, fiter, EBBE3D.get_nls(), &GU);
                    
                    //std::cout << GR << std::endl;

                    //Solve the equation
                    GT.makeCompressed();

                    Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>> solver;
                    solver.analyzePattern(GT);
                    solver.factorize(GT); //LU decomposition

                    //std::cout<<J<<std::endl;

                    assert(solver.info() == Eigen::Success);

                    dGU = solver.solve(-GR);

                    //Calculate errors and update for next iteration
                    for (int j = 0; j < EBBE3D.get_nnode() * EBBE3D.get_ndof(); j++)
                        GU_new(j) = GU(j) + dGU(j);

                    //Error calculation
                    for (int j = 0; j < EBBE3D.get_nnode() * EBBE3D.get_ndof(); j++)
                        error(j) = abs(GU_new(j) - GU(j));

                    max = error(0);
                    for (int j = 0; j < EBBE3D.get_nnode() * EBBE3D.get_ndof(); j++)
                        if (max < error(j))
                            max = error(j);

                    //Assignment for next iteration
                    for (int j = 0; j < EBBE3D.get_nnode() * EBBE3D.get_ndof(); j++)
                        GU(j) = GU_new(j);
                    iter++;

                    //Free memory
                    for (int i = 0; i < EBBE3D.get_ndof() * EBBE3D.get_nen(); i++)
                        delete T[i];
                    delete[] R;

                    std::cout << max << std::endl;
                } while (max > pow(10, -6) && iter < maxiter);
                if (iter < maxiter)
                {
                    VTKGrid VTK1 = VTKGrid(EBBE3D, 16, fiter, "EBBE3D", GU, EBBE3D.get_ndim(), 9, EBBE3D.get_ndof(), 4, 0);
                    
                    /*for (int j = 0; j < EBBE3D.get_nnode() * EBBE3D.get_ndof(); j++)
                    {
                        //file2 << "Displacements of node " << j + 1 << std::endl;
                        //for (int k = 0; k < 3; k++)
                            //file2 << U(12 * (j - 1) + 6 + k) << std::endl;
                        file2 << GU(j) << std::endl;
                    }*/

                    /*std::cout << "Load" << fiter << std::endl;
                    std::cout << "X direction" << std::endl;
                    std::cout << GU(EBBE3D.get_nnode() * EBBE3D.get_ndof() - 6) << std::endl;
                    std::cout << "Y direction" << std::endl;
                    std::cout << GU(EBBE3D.get_nnode() * EBBE3D.get_ndof() - 5) << std::endl;
                    std::cout << "Z direction" << std::endl;
                    std::cout << GU(EBBE3D.get_nnode() * EBBE3D.get_ndof() - 4) << std::endl;
                    file2 << GU(EBBE3D.get_nnode() * EBBE3D.get_ndof() - 4) << std::endl;
                    std::cout << "Iteration" << std::endl;
                    std::cout << iter << std::endl;*/

                    //std::cout << "Load" << fiter << std::endl;
                    //std::cout << GU << std::endl;
                    /*std::cout << "X - direction" << std::endl;
                    std::cout << GU((EBBE3D.get_nnode() + 1) * EBBE3D.get_ndof() / 2) << std::endl;
                    std::cout << "Y - direction" << std::endl;
                    std::cout << GU((EBBE3D.get_nnode() + 1) * EBBE3D.get_ndof() / 2 + 1) << std::endl;
                    std::cout << "Z - direction" << std::endl;
                    std::cout << GU((EBBE3D.get_nnode() + 1) * EBBE3D.get_ndof() / 2 + 2) << std::endl;*/
                    fiter++;


                }
                else
                {
                    std::cout << "Code didn't converge" << std::endl;
                    break;
                }
            } while (fiter <= EBBE3D.get_nls() + 1);
        }
    }
    //--------------------------------------------------------//
    //----------3D EulerBernouli Beam Element-----------------//
    //-----------------Large displacement---------------------//
    //-----------------Updated Lagrangian---------------------//
    //--------------------Quadratic---------------------------//
    //---------------Node to Node Contact---------------------//
    //--------------------------------------------------------//
    else if (choice == 6)
    {

        const int nbeams = 8;

        NLEBBE3D* EBBE3D = new NLEBBE3D[nbeams];

        for (int i = 0; i < nbeams; i++)
        {
            //Textile simulation
            EBBE3D[i] = NLEBBE3D(nbeams, 4.5, 1, 4, "LINELAS", i + 1, 21, 3, 3);
            //Patch test
            //EBBE3D[i] = NLEBBE3D(40, 0.5, 10, 40, 3, 3, pow(10, 7), 0.001, i);
            //Benchmark problem
            //EBBE3D[i] = NLEBBE3D(20, 0.5, 0.3, 15, 3, 3, 2 * pow(10, 11), 0.298, i);
            //EBBE3D[i] = NLEBBE3D(i + 1, "LINELAS");

            std::cout << "New coordinates of beam " << i << std::endl;
            for (int j = 0; j < EBBE3D[i].get_nnode(); j++)
            {
                for (int k = 0; k < EBBE3D[i].get_ndim(); k++)
                    std::cout << EBBE3D[i].get_coordinates(j, k) << " ";
                std::cout << std::endl;
            }
            for (int j = 0; j < EBBE3D[i].get_nelem(); j++)
            {
                for (int k = 0; k < EBBE3D[i].get_nen() + 1; k++)
                    std::cout << EBBE3D[i].get_connectivity(j, k) << " ";
                std::cout << std::endl;
            }
        }

        int size = 0;
        for (int i = 0; i < nbeams; i++)
            size += EBBE3D[i].get_ndof() * EBBE3D[i].get_nnode();

        Eigen::VectorXd dGU = Eigen::VectorXd::Zero(size);
        Eigen::VectorXd GU = Eigen::VectorXd::Zero(size);
        Eigen::VectorXd error = Eigen::VectorXd::Zero(size);
        Eigen::VectorXd GU_new = Eigen::VectorXd::Zero(size);
        Eigen::VectorXd GR = Eigen::VectorXd::Zero(size);
        Eigen::SparseMatrix<double, Eigen::ColMajor> GT(size, size);

        BeamContact BCon{ nbeams, "NTN" };

        double** gap;
        gap = new double* [nbeams];
        for (int i = 0; i < nbeams; i++)
            gap[i] = new double[EBBE3D[i].get_nnode()];

        for (int i = 0; i < nbeams; i++)
            for (int j = 0; j < EBBE3D[i].get_nnode(); j++)
                gap[i][j] = 0;

        double max = 0;
        int iter = 1;
        int fiter = 0;
        int maxiter = 300;
        double conv = pow(10, -6);

        std::fstream file1, file2;
        file1.open("E:/Adhithya/MTech_Res_Thesis/Cpp/ThesisCode/TextileComposites/TextileComposites/Result_Log.txt", std::fstream::in | std::fstream::out);
        file2.open("E:/Adhithya/MTech_Res_Thesis/Cpp/ThesisCode/TextileComposites/TextileComposites/Results.txt", std::fstream::in | std::fstream::out);

        //declare and initialize material parameters and other input variables
        //master
        //std::unique_ptr<double[]> D1 = std::make_unique<double[]>(7);

        for (int i = 0; i < nbeams; i++)
            VTKGrid VTK1 = VTKGrid(EBBE3D[i], 16, fiter, "EBBE3D", GU, EBBE3D[i].get_ndim(), 9, EBBE3D[i].get_ndof(), 4, i, gap, BCon.get_penaltyparameter());

        double load = 0;

        if (file1.is_open() && file2.is_open())
        {
            do
            {
                iter = 0;
                //file2 << "Residual" << "    " << "Load step  " << fiter << std::endl;
                do
                {
                    //Initialize all matrices to zero
                    GU_new.setZero();
                    GR.setZero();
                    GT.setZero();

                    for (int m = 0; m < nbeams; m++)
                    {
                        int temp = 0;
                        for (int j = m - 1; j >= 0; j--)
                            temp += EBBE3D[j].get_nnode() * EBBE3D[j].get_ndof();

                        //declare element tangent stiffness matrix and residual vector
                        double** T;
                        T = new double* [EBBE3D[m].get_ndof() * EBBE3D[m].get_nen()];
                        for (int i = 0; i < EBBE3D[m].get_ndof() * EBBE3D[m].get_nen(); i++)
                            T[i] = new double[EBBE3D[m].get_ndof() * EBBE3D[m].get_nen()];
                        double* R = new double[EBBE3D[m].get_ndof() * EBBE3D[m].get_nen()];

                        double* D1 = new double[7];
                        D1[0] = EBBE3D[m].get_matprop("E");
                        D1[1] = EBBE3D[m].get_matprop("nu");
                        D1[2] = EBBE3D[m].get_matprop("Bp");
                        D1[3] = EBBE3D[m].get_matprop("Hp");
                        D1[4] = EBBE3D[m].get_matprop("Zx");
                        D1[5] = EBBE3D[m].get_matprop("Zy");
                        D1[6] = EBBE3D[m].get_matprop("Zz");

                        for (int i = 0; i < EBBE3D[m].get_nelem(); i++)
                        {
                            int a = (int)EBBE3D[m].get_connectivity(i, 1) - 1;
                            int b = (int)EBBE3D[m].get_connectivity(i, 2) - 1;
                            int c = (int)EBBE3D[m].get_connectivity(i, 3) - 1;

                            //std::cout << a << std::endl;
                            //std::cout << b << std::endl;
                            //std::cout << c << std::endl;

                            //declare and initialize position vector
                            double X[3][3];

                            for (int j = 0; j < 3; j++)
                            {
                                X[0][j] = EBBE3D[m].get_coordinates(a, j);
                                X[1][j] = EBBE3D[m].get_coordinates(b, j);
                                X[2][j] = EBBE3D[m].get_coordinates(c, j);
                            }
                            double h = X[2][0] - X[0][0];

                            /*for (int j = 0; j < 3; j++)
                            {
                                for (int k = 0; k < 3; k++)
                                    std::cout << X[j][k] << " ";
                                std::cout << std::endl;
                            }*/
                            //declare and initialize dofs
                            double U[3][6];

                            for (int j = 0; j < 6; j++)
                            {
                                U[0][j] = GU(EBBE3D[m].get_ndof() * a + j + temp);
                                U[1][j] = GU(EBBE3D[m].get_ndof() * b + j + temp);
                                U[2][j] = GU(EBBE3D[m].get_ndof() * c + j + temp);
                            }

                            /*for (int j = 0; j < 3; j++)
                            {
                                for (int k = 0; k < 6; k++)
                                    std::cout << U[j][k] << " ";
                                std::cout << std::endl;
                            }*/

                            for (int j = 0; j < EBBE3D[m].get_ndof() * EBBE3D[m].get_nen(); j++)
                            {
                                R[j] = 0;
                                for (int k = 0; k < EBBE3D[m].get_ndof() * EBBE3D[m].get_nen(); k++)
                                    T[j][k] = 0;
                            }

                            //Local Residual and tangent matrix
                            EBBE3D[m].RKt(D1, X, U, T, R);

                            //Assembly
                            for (int j = 0; j < EBBE3D[m].get_ndof() * EBBE3D[m].get_nen(); j++)
                            {
                                GR(EBBE3D[m].get_ndof() * a + j + temp) += R[j];
                                for (int k = 0; k < EBBE3D[m].get_ndof() * EBBE3D[m].get_nen(); k++)
                                    GT.coeffRef(EBBE3D[m].get_ndof() * a + j + temp, EBBE3D[m].get_ndof() * a + k + temp) += T[j][k];
                            }
                        }
                    }

                    Loading* LOAD = new Loading[nbeams];

                    //If start is zero that means there is no contact
                    if (BCon.get_startingpoint() == 0)
                    {
                        for (int i = 0; i < nbeams; i++)
                            LOAD[i] = Loading(i + 1, BCon.get_startingpoint());
                    }
                    else
                    {
                        for (int i = 0; i < nbeams; i++)
                            LOAD[i] = Loading(i + 1);
                    }

                    if (BCon.get_startingpoint() == 0)
                    {
                        for (int i = 0; i < nbeams; i++)
                            LOAD[i].SetLoad(&GR, EBBE3D, i, fiter + 1, BCon.get_precontiter());
                    }
                    else
                    {
                        for (int i = 0; i < nbeams; i++)
                            LOAD[i].SetLoad(&GR, EBBE3D, i, fiter - BCon.get_precontiter() + 1, EBBE3D[i].get_nls());
                    }

                    //std::cout << "Residual vector after applying load" << std::endl;
                    //std::cout << GR << std::endl;
                    //std::cout << "Before applying boundary condition" << std::endl;
                    //std::cout << GU << std::endl;

                    Boundary* BOUND = new Boundary[nbeams];

                    if (BCon.get_startingpoint() == 0)
                    {
                        for (int i = 0; i < nbeams; i++)
                            BOUND[i] = Boundary(i + 1, 0.505);
                    }
                    else
                    {
                        for (int i = 0; i < nbeams; i++)
                            BOUND[i] = Boundary(i + 1);
                    }
                    
                    if (BCon.get_startingpoint() == 0)
                    {
                        for (int i = 0; i < nbeams; i++)
                            BOUND[i].SetBoundary(&GT, &GR, EBBE3D, i, fiter + 1, BCon.get_precontiter(), &GU);
                    }
                    else
                    {
                        for (int i = 0; i < nbeams; i++)
                            BOUND[i].SetBoundary(&GT, &GR, EBBE3D, i, fiter - BCon.get_precontiter() + 1, EBBE3D[i].get_nls(), &GU);
                    }

                    if (BCon.get_startingpoint() == 1)
                    {
                        //-------------------------------------------------//
                        //-------------------------------------------------//
                        //--------------------Contact----------------------//
                        //-------------------------------------------------//
                        //-------------------------------------------------//
                        //BCon.GlobalContactSearch(EBBE3D);
                        //Contact Search 
                        std::vector<std::vector<int>> ConPair = BCon.LocalContactSearch_NTN(EBBE3D, &GU, nbeams);

                        if (ConPair.size() == 0)
                            conv = pow(10, -6);
                        else
                            conv = 0.0008;

                        for (int i = 0; i < ConPair.size(); i++)
                        {
                            for (int j = 0; j < ConPair[i].size(); j++)
                                std::cout << ConPair[i][j] << std::endl;
                            std::cout << std::endl;
                        }
                        //declare contact element tangent stiffness matrix and residual vector for contact
                        //Beam elements with different degrees of freedom cannot be used.
                        double** CT;
                        CT = new double* [EBBE3D[0].get_ndof() * BCon.get_nen()];
                        for (int i = 0; i < EBBE3D[0].get_ndof() * BCon.get_nen(); i++)
                            CT[i] = new double[EBBE3D[0].get_ndof() * BCon.get_nen()];
                        double* CR = new double[EBBE3D[0].get_ndof() * BCon.get_nen()];

                        //-------------------------------------------------//
                        //------------Apply contact constraint-------------//
                        //-------------------------------------------------// 
                        for (int i = 0; i < ConPair.size(); i++)
                        {
                            //Find the corresponding master element for the slave element.
                            //Find the corresponding master element for the slave element.
                            int mbeam = ConPair[i][0];
                            int sbeam = ConPair[i][1];
                            int mnode1 = ConPair[i][2];
                            int snode1 = ConPair[i][3];

                            int temp1 = 0;
                            for (int temp = mbeam; temp > 0; temp--)
                                temp1 += EBBE3D[temp].get_nnode() * EBBE3D[temp].get_ndof();

                            int temp2 = 0;
                            for (int temp = sbeam; temp > 0; temp--)
                                temp2 += EBBE3D[temp].get_nnode() * EBBE3D[temp].get_ndof();

                            double X1[3], u1[6], x1[3];

                            double Y1[3], v1[6], y1[3];

                            //Reference config
                            for (int j = 0; j < 3; j++)
                            {
                                X1[j] = EBBE3D[mbeam].get_coordinates(mnode1, j);
                                Y1[j] = EBBE3D[sbeam].get_coordinates(snode1, j);
                            }

                            //Displacements
                            for (int j = 0; j < 6; j++)
                            {
                                u1[j] = GU.coeffRef(EBBE3D[mbeam].get_ndof() * mnode1 + j + temp1);
                                v1[j] = GU.coeffRef(EBBE3D[sbeam].get_ndof() * snode1 + j + temp2);
                            }

                            for (int j = 0; j < 3; j++)
                            {
                                x1[j] = X1[j] + u1[j];
                                y1[j] = Y1[j] + v1[j];
                            }

                            double D[6];
                            D[0] = BCon.get_penaltyparameter() / (iter + 1);
                            //std::cout << D[0] << std::endl;
                            //Eigen::VectorXd n = BCon.set_normal(X1, Y1, BCon.get_normal());
                            Eigen::VectorXd n = BCon.set_normal(x1, y1);
                            //std::cout << n << std::endl;
                            //Eigen::VectorXd n(3);
                            //n(0) = 0;
                            //n(1) = 1;
                            //n(2) = 0;
                            for (int k = 0; k < 3; k++)
                                D[k + 1] = n(k);
                            //std::cout << EBBE3D[sbeam].get_diameter() << std::endl;
                            D[4] = EBBE3D[sbeam].get_diameter() / 2.0;
                            //std::cout << D[5] << std::endl;
                            D[5] = EBBE3D[mbeam].get_diameter() / 2.0;

                            //Initialize contact residual and contact tangent matrix
                            for (int l = 0; l < EBBE3D[sbeam].get_ndof() * BCon.get_nen(); l++)
                            {
                                CR[l] = 0;
                                for (int m = 0; m < EBBE3D[sbeam].get_ndof() * BCon.get_nen(); m++)
                                    CT[l][m] = 0;
                            }

                            double g;
                            BCon.Contact_NTN(D, X1, Y1, u1, v1, CR, CT, &g);
                            gap[mbeam][mnode1] = g;
                            gap[sbeam][snode1] = g;
                            //std::cout << g << std::endl;
                            /*std::cout << "Contact tangent matrix" << std::endl;
                            for (int l = 0; l < EBBE3D[mbeam].get_ndof() * BCon.get_nen(); l++)
                            {
                                for (int m = 0; m < EBBE3D[mbeam].get_ndof() * BCon.get_nen(); m++)
                                    std::cout << CT[l][m] << " ";
                                std::cout << std::endl;
                            }
                            std::cout << "Contact residual" << std::endl;
                            for (int l = 0; l < EBBE3D[mbeam].get_ndof() * BCon.get_nen(); l++)
                            {
                                std::cout << CR[l] << std::endl;
                            }*/
                            //Assembly
                            if (g < 0)
                            {
                                for (int l = 0; l < EBBE3D[sbeam].get_ndof(); l++)
                                {
                                    GR(EBBE3D[mbeam].get_ndof() * mnode1 + l + temp1) += CR[l];
                                    GR(EBBE3D[sbeam].get_ndof() * snode1 + l + temp2) += CR[l + 6];
                                    for (int m = 0; m < EBBE3D[mbeam].get_ndof(); m++)
                                    {
                                        GT.coeffRef(EBBE3D[mbeam].get_ndof() * mnode1 + l + temp1, EBBE3D[mbeam].get_ndof() * mnode1 + m + temp1) += CT[l][m];
                                        GT.coeffRef(EBBE3D[sbeam].get_ndof() * mnode1 + l + temp1, EBBE3D[mbeam].get_ndof() * snode1 + m + temp2) += CT[l][m + 6];
                                        GT.coeffRef(EBBE3D[mbeam].get_ndof() * snode1 + l + temp2, EBBE3D[sbeam].get_ndof() * mnode1 + m + temp1) += CT[l + 6][m];
                                        GT.coeffRef(EBBE3D[mbeam].get_ndof() * snode1 + l + temp2, EBBE3D[mbeam].get_ndof() * snode1 + m + temp2) += CT[l + 6][m + 6];
                                    }
                                }
                            }
                        }
                    }

                    //std::cout << "Residual Vector after applying boundary conditions" << std::endl;
                    //std::cout << GR << std::endl;
                    //std::cout << "After applying boundary condition" << std::endl;
                    //std::cout << GT << std::endl;

                    //Solve the equation
                    GT.makeCompressed();

                    Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>> solver;
                    solver.analyzePattern(GT);
                    solver.factorize(GT); //LU decomposition

                    assert(solver.info() == Eigen::Success);

                    dGU = solver.solve(-GR);

                    //std::cout << dGU << std::endl;
                    //Calculate errors and update for next iteration
                    for (int j = 0; j < size; j++)
                        GU_new(j) = GU(j) + dGU(j);

                    //Error calculation
                    for (int j = 0; j < size; j++)
                        error(j) = abs(GU_new(j) - GU(j));

                    max = error(0);
                    for (int j = 1; j < size; j++)
                        if (max < error(j))
                            max = error(j);

                    std::cout << max << std::endl;

                    double alpha = 1;
                    if (conv > pow(10, -5))
                        alpha = 1;
                    else
                        alpha = 1;
                    //Assignment for next iteration
                    for (int j = 0; j < size; j++)
                        GU(j) = alpha * GU_new(j) + (1 - alpha) * GU(j);
                    iter++;

                    //std::cout << GU << std::endl;
                    /*if (max < pow(10, -8))
                        std::cout << gap << std::endl;*/
                    //file2 << max << std::endl;
                        //Free memory
                        //delete[] TP, RP, CT, T, R, CR, ContactPairs;
                    
                } while ((max > conv && iter < maxiter) || iter == 1);
                if (iter < maxiter)
                {
                    fiter++;
                    //file2 << "Contact pressure  " << fiter << std::endl;
                    //for (int j = 0; j < EBBE3D[0].get_nnode(); j++)
                    //    file2 << gap[0][j] * BCon.get_penaltyparameter() << std::endl;
                    
                    //Change for benchmark problem of contact
                    //for(int j = 0; j < EBBE3D[1].get_nnode()* EBBE3D[1].get_ndof(); j++)
                    //    file2 << GU(size - j - 1) << std::endl;
                    int k;
                    if (BCon.get_startingpoint() == 1)
                    {
                        int size = 0;
                        file2 << "Iteration " << fiter + 1 << std::endl;
                        for (int i = 0; i < nbeams; i++)
                        {
                            double sum = 0;
                            for (int j = 0; j < 3; j++)
                                sum += pow(GU(size + EBBE3D[i].get_ndof() * (EBBE3D[i].get_nnode() - 1) + j), 2);
                            file2 << sqrt(sum) << std::endl;
                            size += EBBE3D[i].get_ndof() * EBBE3D[i].get_nnode();
                        }
                    }
                    //std::cout << fiter << std::endl;
                    /*file2 << "Converged values at load  " << load << std::endl;
                    file2 << "Values for beam 1 " << std::endl;
                    for (int j = 0; j < EBBE3D1.get_nnode() * EBBE3D1.get_ndof(); j++)
                        file2 << GU(j) << std::endl;*/
                    for (int i = 0; i < nbeams; i++)
                        VTKGrid VTK1 = VTKGrid(EBBE3D[i], 16, fiter, "EBBE3D", GU, EBBE3D[i].get_ndim(), 9, EBBE3D[i].get_ndof(), 4, i, gap, BCon.get_penaltyparameter());

                    if (fiter == BCon.get_precontiter())
                    {
                        BCon.set_startingpoint(1);
                        for (int i = 0; i < nbeams; i++)
                            EBBE3D[i].UpdateNodes(EBBE3D, &GU, i);
                        for (int i = 0; i < nbeams; i++)
                        {
                            std::cout << "New coordinates of beam " << i << std::endl;
                            for (int j = 0; j < EBBE3D[i].get_nnode(); j++)
                            {
                                for (int k = 0; k < EBBE3D[i].get_ndim(); k++)
                                    std::cout << EBBE3D[i].get_coordinates(j, k) << " ";
                                std::cout << std::endl;
                            }
                        }
                      
                        TextileMicrostructureGen(EBBE3D, nbeams);
                        conv = pow(10, -3);
                        GU.setZero();
                    }
                    //Print the current configuration in a file
                }
                else
                {
                    std::cout << "Code didn't converge" << std::endl;
                    break;
                }
            } while (fiter < EBBE3D[0].get_nls() + BCon.get_precontiter());
        }
        //delete[] D1, D2;
    }
    //--------------------------------------------------------//
    //----------3D EulerBernouli Beam Element-----------------//
    //-----------------Large displacement---------------------//
    //-----------------Updated Lagrangian---------------------//
    //--------------------Quadratic---------------------------//
    //------------Segment to Segment Contact------------------//
    else if (choice == 7)
    {
        const int nbeams = 2;

        NLEBBE3D* EBBE3D = new NLEBBE3D[nbeams];

        for (int i = 0; i < nbeams; i++)
            EBBE3D[i] = NLEBBE3D(i + 1, "LINELAS");
 
        int size = 0;
        for (int i = 0; i < nbeams; i++)
            size += EBBE3D[i].get_ndof() * EBBE3D[i].get_nnode();

        Eigen::VectorXd dGU = Eigen::VectorXd::Zero(size);
        Eigen::VectorXd GU = Eigen::VectorXd::Zero(size);
        Eigen::VectorXd error = Eigen::VectorXd::Zero(size);
        Eigen::VectorXd GU_new = Eigen::VectorXd::Zero(size);
        Eigen::VectorXd GR = Eigen::VectorXd::Zero(size);
        Eigen::SparseMatrix<double, Eigen::ColMajor> GT(size, size);

        BeamContact BCon{ nbeams, "STS" };

        double** gap;
        gap = new double* [nbeams];
        for (int i = 0; i < nbeams; i++)
            gap[i] = new double[EBBE3D[i].get_nnode()];

        for (int i = 0; i < nbeams; i++)
            for (int j = 0; j < EBBE3D[i].get_nnode(); j++)
                gap[i][j] = 0;

        double max = 0;
        int iter = 1;
        int fiter = 0;
        int maxiter = 100;
        double conv = pow(10, -6);

        std::fstream file1, file2;
        file1.open("E:/Adhithya/MTech_Res_Thesis/Cpp/ThesisCode/TextileComposites/TextileComposites/Result_Log.txt", std::fstream::in | std::fstream::out);
        file2.open("E:/Adhithya/MTech_Res_Thesis/Cpp/ThesisCode/TextileComposites/TextileComposites/Results.txt", std::fstream::in | std::fstream::out);

        //declare and initialize material parameters and other input variables
        //master
        //std::unique_ptr<double[]> D1 = std::make_unique<double[]>(7);

        double load = 0;

        if (file1.is_open() && file2.is_open())
        {
            do
            {
                iter = 0;
                do
                {
                    //Initialize all matrices to zero
                    GU_new.setZero();
                    GR.setZero();
                    GT.setZero();

                    for (int m = 0; m < nbeams; m++)
                    {
                        int temp = 0;
                        for (int j = m - 1; j >= 0; j--)
                            temp += EBBE3D[j].get_nnode() * EBBE3D[j].get_ndof();

                        //declare element tangent stiffness matrix and residual vector
                        double** T;
                        T = new double* [EBBE3D[m].get_ndof() * EBBE3D[m].get_nen()];
                        for (int i = 0; i < EBBE3D[m].get_ndof() * EBBE3D[m].get_nen(); i++)
                            T[i] = new double[EBBE3D[m].get_ndof() * EBBE3D[m].get_nen()];
                        double* R = new double[EBBE3D[m].get_ndof() * EBBE3D[m].get_nen()];

                        double* D1 = new double[7];
                        D1[0] = EBBE3D[m].get_matprop("E");
                        D1[1] = EBBE3D[m].get_matprop("nu");
                        D1[2] = EBBE3D[m].get_matprop("Bp");
                        D1[3] = EBBE3D[m].get_matprop("Hp");
                        D1[4] = EBBE3D[m].get_matprop("Zx");
                        D1[5] = EBBE3D[m].get_matprop("Zy");
                        D1[6] = EBBE3D[m].get_matprop("Zz");

                        for (int i = 0; i < EBBE3D[m].get_nelem(); i++)
                        {
                            int a = (int)EBBE3D[m].get_connectivity(i, 1) - 1;
                            int b = (int)EBBE3D[m].get_connectivity(i, 2) - 1;
                            int c = (int)EBBE3D[m].get_connectivity(i, 3) - 1;

                            //std::cout << a << std::endl;
                            //std::cout << b << std::endl;
                            //std::cout << c << std::endl;

                            //declare and initialize position vector
                            double X[3][3];


                            for (int j = 0; j < 3; j++)
                            {
                                X[0][j] = EBBE3D[m].get_coordinates(a, j);
                                X[1][j] = EBBE3D[m].get_coordinates(b, j);
                                X[2][j] = EBBE3D[m].get_coordinates(c, j);
                            }
                            double h = X[2][0] - X[0][0];

                            //declare and initialize dofs
                            double U[3][6];

                            for (int j = 0; j < 6; j++)
                            {
                                U[0][j] = GU(EBBE3D[m].get_ndof() * a + j + temp);
                                U[1][j] = GU(EBBE3D[m].get_ndof() * b + j + temp);
                                U[2][j] = GU(EBBE3D[m].get_ndof() * c + j + temp);
                            }

                            for (int j = 0; j < EBBE3D[m].get_ndof() * EBBE3D[m].get_nen(); j++)
                            {
                                R[j] = 0;
                                for (int k = 0; k < EBBE3D[m].get_ndof() * EBBE3D[m].get_nen(); k++)
                                    T[j][k] = 0;
                            }

                            //Local Residual and tangent matrix
                            EBBE3D[m].RKt(D1, X, U, T, R);

                            //Assembly
                            for (int j = 0; j < EBBE3D[m].get_ndof() * EBBE3D[m].get_nen(); j++)
                            {
                                GR(EBBE3D[m].get_ndof() * a + j + temp) += R[j];
                                for (int k = 0; k < EBBE3D[m].get_ndof() * EBBE3D[m].get_nen(); k++)
                                    GT.coeffRef(EBBE3D[m].get_ndof() * a + j + temp, EBBE3D[m].get_ndof() * a + k + temp) += T[j][k];
                            }
                        }
                    }

                    Loading* LOAD = new Loading[nbeams];

                    for (int i = 0; i < nbeams; i++)
                        LOAD[i] = Loading(i + 1, BCon.get_startingpoint());

                    if (fiter < BCon.get_precontiter())
                    {
                        for (int i = 0; i < nbeams; i++)
                            LOAD[i].SetLoad(&GR, EBBE3D, i, fiter + 1, BCon.get_precontiter());
                    }
                    else
                    {
                        for (int i = 0; i < nbeams; i++)
                            LOAD[i].SetLoad(&GR, EBBE3D, i, fiter - BCon.get_precontiter() + 1, EBBE3D[i].get_nls());
                    }

                    Boundary* BOUND = new Boundary[nbeams];

                    for (int i = 0; i < nbeams; i++)
                        BOUND[i] = Boundary(i + 1, BCon.get_startingpoint());

                    if (fiter < BCon.get_precontiter())
                    {
                        for (int i = 0; i < nbeams; i++)
                            BOUND[i].SetBoundary(&GT, &GR, EBBE3D, i, fiter + 1, BCon.get_precontiter(), &GU);
                    }
                    else
                    {
                        for (int i = 0; i < nbeams; i++)
                            BOUND[i].SetBoundary(&GT, &GR, EBBE3D, i, fiter - BCon.get_precontiter() + 1, EBBE3D[i].get_nls(), &GU);
                    }

                    if (BCon.get_startingpoint() == 1)
                    {
                        //-------------------------------------------------//
                        //-------------------------------------------------//
                        //--------------------Contact----------------------//
                        //-------------------------------------------------//
                        //-------------------------------------------------//
                        //BCon.GlobalContactSearch(EBBE3D);
                        //Contact Search 
                        std::vector<std::vector<int>> ConPair = BCon.LocalContactSearch_STS(EBBE3D, &GU, nbeams);

                        for (int i = 0; i < ConPair.size(); i++)
                        {
                            for (int j = 0; j < ConPair[i].size(); j++)
                                std::cout << ConPair[i][j] << std::endl;
                            std::cout << std::endl;
                        }
                        //declare contact element tangent stiffness matrix and residual vector for contact
                        //Beam elements with different degrees of freedom cannot be used.
                        double** CT;
                        CT = new double* [EBBE3D[0].get_ndof() * BCon.get_nen()];
                        for (int i = 0; i < EBBE3D[0].get_ndof() * BCon.get_nen(); i++)
                            CT[i] = new double[EBBE3D[0].get_ndof() * BCon.get_nen()];
                        double* CR = new double[EBBE3D[0].get_ndof() * BCon.get_nen()];

                        double** TP;
                        TP = new double* [2];
                        for (int i = 0; i < 2; i++)
                            TP[i] = new double[2];
                        double* RP = new double[2];

                        //-------------------------------------------------//
                        //------------Apply contact constraint-------------//
                        //-------------------------------------------------// 
                        for (int i = 0; i < ConPair.size(); i++)
                        {
                            //Find the corresponding master element for the slave element.
                            int mbeam = ConPair[i][0];
                            int sbeam = ConPair[i][1];
                            int mnode1 = ConPair[i][2];
                            int mnode2 = ConPair[i][3];
                            int snode1 = ConPair[i][4];
                            int snode2 = ConPair[i][5];

                            int temp1 = 0;
                            for (int temp = mbeam; temp > 0; temp--)
                                temp1 += EBBE3D[temp].get_nnode() * EBBE3D[temp].get_ndof();

                            int temp2 = 0;
                            for (int temp = sbeam; temp > 0; temp--)
                                temp2 += EBBE3D[temp].get_nnode() * EBBE3D[temp].get_ndof();

                            double X1[3], X2[3], u1[6], u2[6];

                            double Y1[3], Y2[3], v1[6], v2[6], v3[6];

                            for (int j = 0; j < 3; j++)
                            {
                                X1[j] = EBBE3D[mbeam].get_coordinates(mnode1, j);
                                X2[j] = EBBE3D[mbeam].get_coordinates(mnode2, j);
                                Y1[j] = EBBE3D[sbeam].get_coordinates(snode1, j);
                                Y2[j] = EBBE3D[sbeam].get_coordinates(snode2, j);
                            }

                            for (int j = 0; j < 6; j++)
                            {
                                u1[j] = GU.coeffRef(EBBE3D[mbeam].get_ndof() * mnode1 + j + temp1);
                                u2[j] = GU.coeffRef(EBBE3D[mbeam].get_ndof() * mnode2 + j + temp1);
                                v1[j] = GU.coeffRef(EBBE3D[sbeam].get_ndof() * snode1 + j + temp2);
                                v2[j] = GU.coeffRef(EBBE3D[sbeam].get_ndof() * snode2 + j + temp2);
                            }

                            double x1[3], x2[3], y1[3], y2[3];

                            //current config
                            for (int j = 0; j < 3; j++)
                            {
                                x1[j] = X1[j] + u1[j];
                                x2[j] = X2[j] + u2[j];
                                y1[j] = Y1[j] + v1[j];
                                y2[j] = Y2[j] + v2[j];
                            }

                            double delx[3], dely[3];
                            for (int j = 0; j < 3; j++)
                            {
                                delx[j] = x2[j] - x1[j];
                                dely[j] = y2[j] - y1[j];
                            }

                            //Eigen::VectorXd n = BCon.get_normal(delx, dely, delx, dely);
                            //Eigen::VectorXd n = Eigen::VectorXd::Zero(3);
                            //std::cout << n << std::endl;
                            Eigen::VectorXd n = BCon.set_normal(X1, Y1, BCon.get_normal());

                            double Data2[6];
                            Data2[0] = BCon.get_penaltyparameter();
                            Data2[1] = EBBE3D[mbeam].get_diameter() / 2.0;
                            Data2[2] = EBBE3D[sbeam].get_diameter() / 2.0;
                            Data2[3] = n(0);
                            Data2[4] = n(1);
                            Data2[5] = n(2);

                            for (int j = 0; j < EBBE3D[0].get_ndof() * BCon.get_nen(); j++)
                            {
                                CR[j] = 0;
                                for (int k = 0; k < EBBE3D[0].get_ndof() * BCon.get_nen(); k++)
                                    CT[j][k] = 0;
                            }
                            //segment n1 - n2 and m1 - m2
                            double g;
                            BCon.Contact_STS(Data2, X1, X2, Y1, Y2, u1, u2, v1, v2, CT, CR, RP, &g);
                            std::cout << g << std::endl;
                            std::cout << "Contact tangent matrix" << std::endl;
                            for (int l = 0; l < EBBE3D[mbeam].get_ndof() * BCon.get_nen(); l++)
                            {
                                for (int m = 0; m < EBBE3D[mbeam].get_ndof() * BCon.get_nen(); m++)
                                    std::cout << CT[l][m] << " ";
                                std::cout << std::endl;
                            }
                            std::cout << "Contact residual" << std::endl;
                            for (int l = 0; l < EBBE3D[mbeam].get_ndof() * BCon.get_nen(); l++)
                            {
                                std::cout << CR[l] << std::endl;
                            }
                            if (g < 0)
                            {
                                for (int l = 0; l < EBBE3D[mbeam].get_ndof(); l++)
                                {
                                    GR(EBBE3D[mbeam].get_ndof() * mnode1 + l + temp1) += CR[l];
                                    GR(EBBE3D[mbeam].get_ndof() * mnode2 + l + temp1) += CR[l + 6];
                                    GR(EBBE3D[sbeam].get_ndof() * snode1 + l + temp2) += CR[l + 12];
                                    GR(EBBE3D[sbeam].get_ndof() * snode2 + l + temp2) += CR[l + 18];
                                    for (int m = 0; m < EBBE3D[mbeam].get_ndof(); m++)
                                    {
                                        GT.coeffRef(EBBE3D[mbeam].get_ndof() * mnode1 + l + temp1, EBBE3D[mbeam].get_ndof() * mnode1 + m + temp1) += CT[l][m];
                                        GT.coeffRef(EBBE3D[mbeam].get_ndof() * mnode1 + l + temp1, EBBE3D[mbeam].get_ndof() * mnode2 + m + temp1) += CT[l][m + 6];
                                        GT.coeffRef(EBBE3D[mbeam].get_ndof() * mnode1 + l + temp1, EBBE3D[mbeam].get_ndof() * snode1 + m + temp2) += CT[l][m + 12];
                                        GT.coeffRef(EBBE3D[mbeam].get_ndof() * mnode1 + l + temp1, EBBE3D[mbeam].get_ndof() * snode2 + m + temp2) += CT[l][m + 18];

                                        GT.coeffRef(EBBE3D[mbeam].get_ndof() * mnode2 + l + temp1, EBBE3D[mbeam].get_ndof() * mnode1 + m + temp1) += CT[l + 6][m];
                                        GT.coeffRef(EBBE3D[mbeam].get_ndof() * mnode2 + l + temp1, EBBE3D[mbeam].get_ndof() * mnode2 + m + temp1) += CT[l + 6][m + 6];
                                        GT.coeffRef(EBBE3D[mbeam].get_ndof() * mnode2 + l + temp1, EBBE3D[mbeam].get_ndof() * snode1 + m + temp2) += CT[l + 6][m + 12];
                                        GT.coeffRef(EBBE3D[mbeam].get_ndof() * mnode2 + l + temp1, EBBE3D[mbeam].get_ndof() * snode2 + m + temp2) += CT[l + 6][m + 18];

                                        GT.coeffRef(EBBE3D[mbeam].get_ndof() * snode1 + l + temp2, EBBE3D[mbeam].get_ndof() * mnode1 + m + temp1) += CT[l + 12][m];
                                        GT.coeffRef(EBBE3D[mbeam].get_ndof() * snode1 + l + temp2, EBBE3D[mbeam].get_ndof() * mnode2 + m + temp1) += CT[l + 12][m + 6];
                                        GT.coeffRef(EBBE3D[mbeam].get_ndof() * snode1 + l + temp2, EBBE3D[mbeam].get_ndof() * snode1 + m + temp2) += CT[l + 12][m + 12];
                                        GT.coeffRef(EBBE3D[mbeam].get_ndof() * snode1 + l + temp2, EBBE3D[mbeam].get_ndof() * snode2 + m + temp2) += CT[l + 12][m + 18];

                                        GT.coeffRef(EBBE3D[mbeam].get_ndof() * snode2 + l + temp2, EBBE3D[mbeam].get_ndof() * mnode1 + m + temp1) += CT[l + 18][m];
                                        GT.coeffRef(EBBE3D[mbeam].get_ndof() * snode2 + l + temp2, EBBE3D[mbeam].get_ndof() * mnode2 + m + temp1) += CT[l + 18][m + 6];
                                        GT.coeffRef(EBBE3D[mbeam].get_ndof() * snode2 + l + temp2, EBBE3D[mbeam].get_ndof() * snode1 + m + temp2) += CT[l + 18][m + 12];
                                        GT.coeffRef(EBBE3D[mbeam].get_ndof() * snode2 + l + temp2, EBBE3D[mbeam].get_ndof() * snode2 + m + temp2) += CT[l + 18][m + 18];
                                    }
                                }
                            }
                        }
                        //Special case for contact between endpoints 
                        //Enforce constraint for closest enpoints
                        //closest points n1, m1
                        /*Contact_STS_Endpoints(Data3, X1, Y1, u1, v1, &GT, &GR, n1, m1);

                        //closest points n2, m2
                        Contact_STS_Endpoints(Data3, X2, Y2, u2, v2, &GT, &GR, n2, m2);

                        //closest points n3, m3
                        Contact_STS_Endpoints(Data3, X3, Y3, u3, v3, &GT, &GR, n3, m3);*/
                    }

                    //std::cout << "Residual vector before applying load" << std::endl;
                                        //std::cout << GR << std::endl;
                        
                    //std::cout << "Residual vector after applying load" << std::endl;
                    //std::cout << GR << std::endl;
                    //std::cout << "Before applying boundary condition" << std::endl;
                    //std::cout << GU << std::endl;

                            
                            
                            
                    //std::cout << "Residual Vector after applying boundary conditions" << std::endl;
                    //std::cout << GR << std::endl;
                    //std::cout << "After applying boundary condition" << std::endl;
                    //std::cout << GU << std::endl;

                    //Solve the equation
                    GT.makeCompressed();

                    Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>> solver;
                    solver.analyzePattern(GT);
                    solver.factorize(GT); //LU decomposition

                    assert(solver.info() == Eigen::Success);

                    dGU = solver.solve(-GR);

                    //std::cout << dGU << std::endl;

                    //Calculate errors and update for next iteration
                    for (int j = 0; j < size; j++)
                        GU_new(j) = GU(j) + dGU(j);

                    //Error calculation
                    for (int j = 0; j < size; j++)
                        error(j) = abs(GU_new(j) - GU(j));

                    max = error(0);
                    for (int j = 1; j < size; j++)
                        if (max < error(j))
                            max = error(j);

                    std::cout << max << std::endl;

                    //Assignment for next iteration
                    for (int j = 0; j < size; j++)
                        GU(j) = GU_new(j);
                    iter++;
                        
                    /*if (max < pow(10, -8))
                        std::cout << gap << std::endl;*/    

                    //Free memory
                    //delete[] TP, RP, CT, T, R, CR, ContactPairs;

                        
                } while ((max > conv && iter < maxiter )|| iter == 1);
                if (iter < maxiter)
                {
                    fiter++;
                    /*file2 << "Converged values at load  " << load << std::endl;
                    file2 << "Values for beam 1 " << std::endl;
                    for (int j = 0; j < EBBE3D1.get_nnode() * EBBE3D1.get_ndof(); j++)
                        file2 << GU(j) << std::endl;*/
                    for (int i = 0; i < nbeams; i++)
                        VTKGrid VTK1 = VTKGrid(EBBE3D[i], 16, fiter, "EBBE3D", GU, EBBE3D[i].get_ndim(), 9, EBBE3D[i].get_ndof(), 4, i, gap, BCon.get_penaltyparameter());

                    if (fiter == BCon.get_precontiter())
                    {
                        BCon.set_startingpoint(1);
                        for (int i = 0; i < nbeams; i++)
                            EBBE3D[i].UpdateNodes(EBBE3D, &GU, i);
                        for (int i = 0; i < nbeams; i++)
                        {
                            std::cout << "New coordinates of beam " << i << std::endl;
                            for (int j = 0; j < EBBE3D[i].get_nnode(); j++)
                            {
                                for (int k = 0; k < EBBE3D[i].get_ndim(); k++)
                                    std::cout << EBBE3D[i].get_coordinates(j, k) << " ";
                                std::cout << std::endl;
                            }
                            GU.setZero();
                        }
                    }
                }
                else
                {
                    std::cout << "Code didn't converge" << std::endl;
                    break;
                }
            } while (fiter < EBBE3D[0].get_nls() + BCon.get_precontiter());
        }
        //delete[] D1, D2;
    }
    //----------------------VAM Based-------------------------//
    //----------Geometrically Exact Beam Element--------------//
    //-----------------Large displacement---------------------//
    //-----------------Updated Lagrangian---------------------//
    //----------------------Linear----------------------------//
    //----------------Node to Node contact--------------------//
    else if (choice == 8)
    {
        const int nbeams = 4;

        VAMBeamElement* VAMBE = new VAMBeamElement[nbeams];

        ReadTextileInp(VAMBE, 10, 12, 3);

        for (int i = 0; i < nbeams; i++)
        {
            std::cout << VAMBE[i].get_nls() << std::endl;
            std::cout << VAMBE[i].get_nnode() << std::endl;
            std::cout << VAMBE[i].get_ndof() << std::endl;
            std::cout << VAMBE[i].get_ndim() << std::endl;
            std::cout << VAMBE[i].get_nelem() << std::endl;
            std::cout << VAMBE[i].get_diameter() << std::endl;
            std::cout << VAMBE[i].get_axis() << std::endl;
            std::cout << "New coordinates of beam " << i << std::endl;
            for (int j = 0; j < VAMBE[i].get_nnode(); j++)
            {
                for (int k = 0; k < VAMBE[i].get_ndim(); k++)
                    std::cout << VAMBE[i].get_coordinates(j, k) << " ";
                std::cout << std::endl;
            }
            for (int j = 0; j < VAMBE[i].get_nelem(); j++)
            {
                for (int k = 0; k < 3; k++)
                    std::cout << VAMBE[i].get_connectivity(j, k) << " ";
                std::cout << std::endl;
            }
        }

        //for (int i = 0; i < nbeams; i++)
        //    VAMBE[i] = VAMBeamElement(i + 1);

        int size = 0;
        for (int i = 0; i < nbeams; i++)
            size += VAMBE[i].get_ndof() * VAMBE[i].get_nnode();

        Eigen::VectorXd dU = Eigen::VectorXd::Zero(size);
        Eigen::VectorXd U = Eigen::VectorXd::Zero(size);
        Eigen::VectorXd error = Eigen::VectorXd::Zero(size);
        Eigen::VectorXd U_new = Eigen::VectorXd::Zero(size);
        Eigen::VectorXd R = Eigen::VectorXd::Zero(size);
        Eigen::SparseMatrix<double, Eigen::ColMajor> J(size, size);
        Eigen::VectorXd Tau = Eigen::VectorXd::Zero(6 * size / 12);

        BeamContact BCon{ nbeams, "NTN" };

        double** gap;
        gap = new double* [nbeams];
        for (int i = 0; i < nbeams; i++)
            gap[i] = new double[VAMBE[i].get_nnode()];

        for (int i = 0; i < nbeams; i++)
            for (int j = 0; j < VAMBE[i].get_nnode(); j++)
                gap[i][j] = 0;

        double max;
        int iter = 1;
        int fiter = 0;
        int maxiter = 100;
        std::fstream file1, file2;
        file1.open("E:/Adhithya/MTech_Res_Thesis/Cpp/ThesisCode/TextileComposites/TextileComposites/Result_Log.txt", std::fstream::in | std::fstream::out);
        file2.open("E:/Adhithya/MTech_Res_Thesis/Cpp/ThesisCode/TextileComposites/TextileComposites/Results.txt", std::fstream::in | std::fstream::out);

        //Don't know what Theta means but it is used in C_ab formula
       
        double conv = pow(10, -6);

        //Loop for force iterations.
        if (file1.is_open() && file2.is_open())
        {
            do
            {
                iter = 0;
                do
                {
                    //Initialize all matrices to zero
                    U_new.setZero();
                    R.setZero();
                    J.setZero();
                    Tau.setZero();

                    Eigen::VectorXd U1(6), F1(6), U2(6), F2(6);
                    Eigen::MatrixXd Seq_1 = Eigen::MatrixXd::Zero(6, 6);
                    Eigen::MatrixXd Seq_2 = Eigen::MatrixXd::Zero(6, 6);

                    for (int l = 0; l < nbeams; l++)
                    {
                        int temp = 0;
                        for (int i = 0; i < l; i++)
                            temp += VAMBE[i].get_nnode() * VAMBE[i].get_ndof();
                        Eigen::VectorXd axis = VAMBE[l].get_axis();
                        double Theta;
                        if (axis(0) == 1)
                            Theta = 0;
                        else if (axis(1) == 1)
                            Theta = 90 * M_PI / 180;

                        for (int i = 0; i < VAMBE[l].get_nelem(); i++)
                        {
                            int a = (int)VAMBE[l].get_connectivity(i, 1) - 1;
                            int b = (int)VAMBE[l].get_connectivity(i, 2) - 1;
                            Eigen::VectorXd xa(3), xb(3);
                            for (int j = 0; j < 3; j++)
                            {
                                xa(j) = VAMBE[l].get_coordinates(a, j);
                                xb(j) = VAMBE[l].get_coordinates(b, j);
                            }
                            //std::cout << xa << std::endl;
                            //std::cout << xb << std::endl;
                            double h = 0;
                            
                            //std::cout << axis << std::endl;
                            for (int j = 0; j < 3; j++)
                                h += abs(xb(j) * axis(j) - xa(j) * axis(j));
                            //std::cout << h << std::endl;
                            if (i < VAMBE[l].get_nelem() - 1)
                            {
                                if (i == 0)
                                {
                                    Eigen::VectorXd Element_R = Eigen::VectorXd::Zero(12);
                                    Eigen::MatrixXd Element_J = Eigen::MatrixXd::Zero(12, 18);

                                    Eigen::VectorXd U0(6);
                                    for (int j = 0; j < 6; j++)
                                    {
                                        U1(j) = U(a * 12 + 6 + j + temp);
                                        F1(j) = U(a * 12 + 12 + j + temp);
                                        //Fhat
                                        //Mhat
                                        U0(j) = U(j + temp);
                                    }

                                    Seq_1 = VAMBE[l].Equivalent_StiffnessMatrix_FirstOrder(file1);

                                    Element_R = VAMBE[l].Element_Residual(U1, F1, h, Seq_1, Theta, 0, U0);
                                    Element_J = VAMBE[l].Element_Jacobian(U1, F1, h, Seq_1, Theta, 0, file1);

                                    //Assembly of Node 0 into Global Stiffness matrix
                                    //fu1 - F1 = 0
                                    //fpsi1 - M1 = 0
                                    //fF1 - u1 = 0
                                    //fM1 - theta1 = 0
                                    for (int j = 0; j < 12; j++)
                                    {
                                        R.coeffRef(j + temp) += Element_R(j);
                                        for (int k = 0; k < 18; k++)
                                        {
                                            J.coeffRef(j + temp, k + temp) += Element_J(j, k);
                                        }
                                    }
                                }

                                for (int j = 0; j < 6; j++)
                                {
                                    U1(j) = U(a * 12 + 6 + j + temp);
                                    F1(j) = U(a * 12 + 12 + j + temp);
                                    U2(j) = U(b * 12 + 6 + j + temp);
                                    F2(j) = U(b * 12 + 12 + j + temp);
                                }

                                Seq_1 = VAMBE[l].Equivalent_StiffnessMatrix_FirstOrder(file1);
                                Seq_2 = VAMBE[l].Equivalent_StiffnessMatrix_FirstOrder(file1);

                                Eigen::VectorXd Element_R = Eigen::VectorXd::Zero(12);
                                Eigen::MatrixXd Element_J = Eigen::MatrixXd::Zero(12, 24);

                                //Element Residual
                                Element_R = VAMBE[l].Element_Residual(U1, U2, F1, F2, h, Seq_1, Seq_2, Theta);
                                //Element Jacobian
                                Element_J = VAMBE[l].Element_Jacobian(U1, U2, F1, F2, h, Seq_1, Seq_2, Theta, file1);

                                //Assembly
                                //Solve for intermediate nodes
                                //fui + fui+1
                                //fpsii + fpsii+1
                                //fFi + fFi+1
                                //fMi + fMi+1
                                for (int j = 0; j < 12; j++)
                                {
                                    R.coeffRef(12 * b + j + temp) += Element_R(j);
                                    for (int k = 0; k < 24; k++)
                                        J.coeffRef(12 * b + j + temp, 12 * a + 6 + k + temp) += Element_J(j, k);
                                }
                            }
                            else if (i == VAMBE[l].get_nelem() - 1)
                            {
                                Eigen::VectorXd Element_R = Eigen::VectorXd::Zero(12);
                                Eigen::MatrixXd Element_J = Eigen::MatrixXd::Zero(12, 18);

                                Eigen::VectorXd U0(6);
                                for (int j = 0; j < 6; j++)
                                {
                                    U1(j) = U(12 * a + 6 + j + temp);
                                    F1(j) = U(12 * a + 12 + j + temp);
                                    //FN+1
                                    //MN+1
                                    U0(j) = U(12 * b + 6 + j + temp);
                                }

                                Seq_1 = VAMBE[l].Equivalent_StiffnessMatrix_FirstOrder(file1);

                                Element_R = VAMBE[l].Element_Residual(U1, F1, h, Seq_1, Theta, VAMBE[l].get_nnode() - 1, U0);
                                Element_J = VAMBE[l].Element_Jacobian(U1, F1, h, Seq_1, Theta, VAMBE[l].get_nnode() - 1, file1);

                                //fuN - FN+1 = 0
                                //fpsiN - MN+1 = 0
                                //fFN + uN+1 = 0
                                //fMN + thetaN+1 = 0
                                for (int j = 0; j < 12; j++)
                                {
                                    R.coeffRef(12 * b + j + temp) += Element_R(j);
                                    for (int k = 0; k < 18; k++)
                                        J.coeffRef(12 * b + j + temp, 12 * a + 6 + k + temp) += Element_J(j, k);
                                }
                            }
                        }
                    }

                    //std::cout << R << std::setprecision(3) << std::endl;

                    //std::cout << J << std::setprecision(3) << std::endl;

                    Loading* LOAD = new Loading[nbeams];

                    for (int i = 0; i < nbeams; i++)
                        LOAD[i] = Loading(i + 1, BCon.get_startingpoint());

                    if (BCon.get_startingpoint() == 0)
                    {
                        for (int i = 0; i < nbeams; i++)
                            LOAD[i].SetLoad(&R, VAMBE, i, fiter + 1, BCon.get_precontiter());
                    }
                    else
                    {
                        for (int i = 0; i < nbeams; i++)
                            LOAD[i].SetLoad(&R, VAMBE, i, fiter - BCon.get_precontiter() + 1, VAMBE[i].get_nls());
                    }

                    if (BCon.get_startingpoint() == 1)
                    {
                        //-------------------------------------------------//
                        //-------------------------------------------------//
                        //--------------------Contact----------------------//
                        //-------------------------------------------------//
                        //-------------------------------------------------//
                        //BCon.GlobalContactSearch(EBBE3D);
                        //Contact Search 
                        std::vector<std::vector<int>> ConPair = BCon.LocalContactSearch_NTN(VAMBE, &U, nbeams);

                        if (ConPair.size() == 0)
                            conv = pow(10, -6);
                        else
                            conv = 15;

                        for (int i = 0; i < ConPair.size(); i++)
                        {
                            for (int j = 0; j < ConPair[i].size(); j++)
                                std::cout << ConPair[i][j] << std::endl;
                            std::cout << std::endl;
                        }
                        //declare contact element tangent stiffness matrix and residual vector for contact
                        //Beam elements with different degrees of freedom cannot be used.
                        double** CT;
                        CT = new double* [VAMBE[0].get_ndof() * BCon.get_nen()];
                        for (int i = 0; i < VAMBE[0].get_ndof() * BCon.get_nen(); i++)
                            CT[i] = new double[VAMBE[0].get_ndof() * BCon.get_nen()];
                        double* CR = new double[VAMBE[0].get_ndof() * BCon.get_nen()];

                        //-------------------------------------------------//
                        //------------Apply contact constraint-------------//
                        //-------------------------------------------------// 
                        for (int i = 0; i < ConPair.size(); i++)
                        {
                            //Find the corresponding master element for the slave element.
                            int mbeam = ConPair[i][0];
                            int sbeam = ConPair[i][1];
                            int mnode1 = ConPair[i][2];
                            int snode1 = ConPair[i][3];

                            int temp1 = 0;
                            for (int temp = mbeam; temp > 0; temp--)
                                temp1 += VAMBE[temp].get_nnode() * VAMBE[temp].get_ndof();

                            int temp2 = 0;
                            for (int temp = sbeam; temp > 0; temp--)
                                temp2 += VAMBE[temp].get_nnode() * VAMBE[temp].get_ndof();

                            double X1[3], u1[6];

                            double Y1[3], v1[6];

                            //Reference config
                            for (int j = 0; j < 3; j++)
                            {
                                X1[j] = VAMBE[mbeam].get_coordinates(mnode1, j);
                                Y1[j] = VAMBE[sbeam].get_coordinates(snode1, j);
                            }

                            //Displacements
                            for (int j = 0; j < 6; j++)
                            {
                                u1[j] = U.coeffRef(VAMBE[mbeam].get_ndof() * mnode1 + j + temp1 + 6);
                                v1[j] = U.coeffRef(VAMBE[sbeam].get_ndof() * snode1 + j + temp2 + 6);
                            }

                            double D[6];
                            D[0] = BCon.get_penaltyparameter();

                            Eigen::VectorXd n = BCon.set_normal(X1, Y1, BCon.get_normal());
                            //std::cout << n << std::endl;
                            //Eigen::VectorXd n(3);
                            //n(0) = 0;
                            //n(1) = 1;
                            //n(2) = 0;
                            for (int k = 0; k < 3; k++)
                                D[k + 1] = n(k);
                            //std::cout << EBBE3D[sbeam].get_diameter() << std::endl;
                            D[4] = VAMBE[sbeam].get_diameter() / 2.0;
                            //std::cout << D[5] << std::endl;
                            D[5] = VAMBE[mbeam].get_diameter() / 2.0;

                            //Initialize contact residual and contact tangent matrix
                            for (int l = 0; l < VAMBE[sbeam].get_ndof() * BCon.get_nen(); l++)
                            {
                                CR[l] = 0;
                                for (int m = 0; m < VAMBE[sbeam].get_ndof() * BCon.get_nen(); m++)
                                    CT[l][m] = 0;
                            }

                            double g;
                            BCon.Contact_NTN(D, X1, Y1, u1, v1, CR, CT, &g);
                            gap[mbeam][mnode1] = g;
                            gap[sbeam][snode1] = g;
                            std::cout << g << std::endl;
                            std::cout << "Contact tangent matrix" << std::endl;
                            for (int l = 0; l < 6 * BCon.get_nen(); l++)
                            {
                                for (int m = 0; m < 6 * BCon.get_nen(); m++)
                                    std::cout << CT[l][m] << " ";
                                std::cout << std::endl;
                            }
                            std::cout << "Contact residual" << std::endl;
                            for (int l = 0; l < 6 * BCon.get_nen(); l++)
                            {
                                std::cout << CR[l] << std::endl;
                            }
                            //Assembly
                            if (g < 0)
                            {
                                for (int l = 0; l < 6; l++)
                                {
                                    R(VAMBE[mbeam].get_ndof() * mnode1 + l + temp1 + 6) += CR[l];
                                    R(VAMBE[sbeam].get_ndof() * snode1 + l + temp2 + 6) += CR[l + 6];
                                    for (int m = 0; m < 6; m++)
                                    {
                                        J.coeffRef(VAMBE[mbeam].get_ndof() * mnode1 + l + temp1 + 6, VAMBE[mbeam].get_ndof() * mnode1 + m + temp1 + 6) += CT[l][m];
                                        J.coeffRef(VAMBE[sbeam].get_ndof() * mnode1 + l + temp1 + 6, VAMBE[mbeam].get_ndof() * snode1 + m + temp2 + 6) += CT[l][m + 6];
                                        J.coeffRef(VAMBE[mbeam].get_ndof() * snode1 + l + temp2 + 6, VAMBE[sbeam].get_ndof() * mnode1 + m + temp1 + 6) += CT[l + 6][m];
                                        J.coeffRef(VAMBE[mbeam].get_ndof() * snode1 + l + temp2 + 6, VAMBE[mbeam].get_ndof() * snode1 + m + temp2 + 6) += CT[l + 6][m + 6];
                                    }
                                }
                            }
                        }
                    }

                    //std::cout << "Residual vector before applying load" << std::endl;
                    //std::cout << GR << std::endl;

                    //Solve the equation
                    J.makeCompressed();

                    Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>> solver;
                    solver.analyzePattern(J);
                    solver.factorize(J); //LU decomposition

                    assert(solver.info() == Eigen::Success);

                    dU = solver.solve(-R);

                    //Calculate errors and update for next iteration
                    for (int j = 0; j < size; j++)
                        U_new(j) = U(j) + dU(j);

                    //Error calculation
                    for (int j = 0; j < size; j++)
                        error(j) = abs(U_new(j) - U(j));

                    max = error(0);
                    for (int j = 0; j < size; j++)
                        if (max < error(j))
                            max = error(j);

                    //Assignment for next iteration
                    double alpha = 0;
                    for (int j = 0; j < size; j++)
                        U(j) = alpha * U(j) + (1 - alpha) * U_new(j);
                    iter++;

                    //After calculating the internal forces through 1D nonlinear beam analysis,
                    //the strains are updated using the constitutive relation.
                    /*Tau = VAMBE[0].Update_Strains(VAMBE, nbeams, &U, file1);

                    file1 << "            Updates Strains After " << iter << " iteration" << std::endl;
                    for (int j = 0; j < VAMBE[0].get_nnode() * 6; j++)
                    {
                        file1 << "                " << Tau(j) << std::endl;
                    }*/
                    std::cout << max << std::endl;
                    /*for (int j = 0; j < size; j++)
                        std::cout << dU(j) << std::endl;*/
                } while (max > pow(10, -6) && iter < maxiter);
                if (iter < maxiter)
                {
                    fiter++;
                    /*file2 << VAMBE.get_load(3) << "     " << U((VAMBE.get_nnode() - 1) * 12 - 18) << "    " << iter << std::endl;
                    std::cout << VAMBE.get_load(3) << "     " << U(12 * (VAMBE.get_nnode() - 2) + 6 + 2) << "    " << iter << std::endl;*/
                    for (int j = 0; j < size; j++)
                        std::cout << U(j) << std::endl;

                    for (int i = 0; i < nbeams; i++)
                        VTKGrid VTK1 = VTKGrid(VAMBE[i], 16, fiter, "VAMBE", U, VAMBE[i].get_ndim(), 9, VAMBE[i].get_ndof(), 4, i, gap, BCon.get_penaltyparameter());

                }
                else
                {
                    file2 << "Code didn't converge for force   " << VAMBE[0].get_load(3) << std::endl;
                    break;
                }
            } while (fiter <= VAMBE[0].get_nls());
        }
        else
        {
            if (!file1.is_open())
                std::cout << "Couldn't open file Result_Log.txt" << std::endl;
            if (!file2.is_open())
                std::cout << "Couldn't open file Results.txt" << std::endl;
        }
        file1.close();
        file2.close();
    }
    else
        std::cout << "Wrong option" << std::endl;

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << duration.count() << " seconds" << std::endl;
}





#ifndef VARIABLES_H
#define VARIABLES_H
#include "Variables.h"
#include<fstream>
#include<limits>


int main()
{
    // 1 - Linear Bar Element
    // 2 - Nonlinear Bar Element
    // 3 - Nonlinear Euler Bernouli Beam Element (3 dofs per node)
    // 4 - GEBT(Geometrically Exact Beam theory) based beam element
    // 5 - Nonlinear Euler Bernouli Beam Element (6 dofs per node)(quadratic)
    // 6 - Nonlinear Euler Bernouli Beam Element with contact using Node to Node (2 bodies)
    // 7 - Nonlinear Euler Bernouli beam element with Segment to Segment contact (2 bodies)
    int choice = 5;
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
            }
            if (iter == maxiter)
            {
                EBBE2D.set_loadprop("vf", EBBE2D.get_loadprop("vf") - 0.5);
                std::cout << "Solution not converging for the given load" << std::endl;
                std::cout << "Trying with a load with less increment" << std::endl;
            }
            std::cout << "Load  " << EBBE2D.get_loadprop("vf") << std::endl;
            std::cout << "Displacements  " << std::endl << U << std::endl;

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
        VAMBeamElement VAMBE = ReadVAMBEFile();

        Eigen::VectorXd dU = Eigen::VectorXd::Zero(VAMBE.NDOF * VAMBE.NNODE);
        Eigen::VectorXd U = Eigen::VectorXd::Zero(VAMBE.NDOF * VAMBE.NNODE);
        Eigen::VectorXd error = Eigen::VectorXd::Zero(VAMBE.NDOF * VAMBE.NNODE);
        Eigen::VectorXd U_new = Eigen::VectorXd::Zero(VAMBE.NDOF * VAMBE.NNODE);
        Eigen::VectorXd R = Eigen::VectorXd::Zero(VAMBE.NDOF * VAMBE.NNODE);
        Eigen::SparseMatrix<double, Eigen::ColMajor> J(VAMBE.NNODE * VAMBE.NDOF, VAMBE.NNODE * VAMBE.NDOF);
        Eigen::VectorXd Tau = Eigen::VectorXd::Zero(6 * VAMBE.NNODE);
        double max;
        int iter = 1;
        int fiter = 0;
        int maxiter = 50;
        std::fstream file1, file2;
        file1.open("E:/Adhithya/MTech_Res_Thesis/Cpp/ThesisCode/TextileComposites/TextileComposites/Result_Log.txt", std::fstream::in | std::fstream::out);
        file2.open("E:/Adhithya/MTech_Res_Thesis/Cpp/ThesisCode/TextileComposites/TextileComposites/Results.txt", std::fstream::in | std::fstream::out);

        //Don't know what Theta means but it is used in C_ab formula
        double Theta = 90 * M_PI / 180;
        double load = VAMBE.B.aN(0);
        //Loop for force iterations.
        if (file1.is_open() && file2.is_open())
        {
            do
            {
                iter = 0;
                VAMBE.B.aN(0) = fiter * load / VAMBE.NLS;
                std::cout << VAMBE.B.aN(0) << std::endl;
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

                    for (int i = 0; i < VAMBE.NELEM; i++)
                    {
                        int a = (int)VAMBE.ELEM(i, 1) - 1;
                        int b = (int)VAMBE.ELEM(i, 2) - 1;
                        double xa = VAMBE.NODE(a, 0);
                        double xb = VAMBE.NODE(b, 0);
                        double h = xb - xa;
                        if (i < VAMBE.NELEM - 1)
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
                                Seq_1 = Equivalent_StiffnessMatrix_FirstOrder(Strain, VAMBE.inittwist, VAMBE.CS.Rect.width, file1);

                                Element_R = Element_Residual(U1, F1, h, Seq_1, Theta, 0, VAMBE.B, U0);
                                Element_J = Element_Jacobian(U1, F1, h, Seq_1, Theta, 0, file1);

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

                            Seq_1 = Equivalent_StiffnessMatrix_FirstOrder(Strain1, VAMBE.inittwist, VAMBE.CS.Rect.width, file1);
                            Seq_2 = Equivalent_StiffnessMatrix_FirstOrder(Strain2, VAMBE.inittwist, VAMBE.CS.Rect.width, file1);

                            Eigen::VectorXd Element_R = Eigen::VectorXd::Zero(12);
                            Eigen::MatrixXd Element_J = Eigen::MatrixXd::Zero(12, 24);

                            //Element Residual
                            Element_R = Element_Residual(U1, U2, F1, F2, h, Seq_1, Seq_2, Theta);
                            //Element Jacobian
                            Element_J = Element_Jacobian(U1, U2, F1, F2, h, Seq_1, Seq_2, Theta, file1);

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
                        else if (i == VAMBE.NELEM - 1)
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

                            Seq_1 = Equivalent_StiffnessMatrix_FirstOrder(Strain, VAMBE.inittwist, VAMBE.CS.Rect.width, file1);

                            Element_R = Element_Residual(U1, F1, h, Seq_1, Theta, VAMBE.NNODE - 1, VAMBE.B, U0);
                            Element_J = Element_Jacobian(U1, F1, h, Seq_1, Theta, VAMBE.NNODE - 1, file1);

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

                    //Calculate errors and update for next iteration
                    for (int j = 0; j < VAMBE.NDOF * VAMBE.NNODE; j++)
                        U_new(j) = U(j) + dU(j);

                    //Error calculation
                    for (int j = 0; j < VAMBE.NDOF * VAMBE.NNODE; j++)
                        error(j) = abs(U_new(j) - U(j));

                    max = error(0);
                    for (int j = 0; j < VAMBE.NDOF * VAMBE.NNODE; j++)
                        if (max < error(j))
                            max = error(j);

                    //Assignment for next iteration
                    for (int j = 0; j < VAMBE.NDOF * VAMBE.NNODE; j++)
                        U(j) = U_new(j);
                    iter++;

                    //After calculating the internal forces through 1D nonlinear beam analysis,
                    //the strains are updated using the constitutive relation.
                    Tau = Update_Strains(VAMBE, &U, file1);

                    file1 << "            Updates Strains After " << iter << " iteration" << std::endl;
                    for (int j = 0; j < VAMBE.NNODE * 6; j++)
                    {
                        file1 << "                " << Tau(j) << std::endl;
                    }
                    std::cout << max << std::endl;
                } while (max > pow(10, -6) && iter < maxiter);
                if (iter < maxiter)
                {
                    fiter++;
                    file2 << VAMBE.B.aN(1) << "     " << U((VAMBE.NNODE - 1) * 12 - 18) << "    " << iter << std::endl;
                    std::cout << VAMBE.B.aN(1) << "     " << U(12 * (VAMBE.NNODE - 2) + 6 + 2) << "    " << iter << std::endl;
                    file2 << "Converged Displacements for load " << VAMBE.B.aN(1) << std::endl;
                    for (int j = 0; j < VAMBE.NNODE * VAMBE.NDOF; j++)
                        file2 << U(j) << std::endl;
                }
                else
                {
                    file2 << "Code didn't converge for force   " << VAMBE.B.aN(0) << std::endl;
                    VAMBE.B.aN(0) -= 0.05;
                    file2 << "Trying with force" << VAMBE.B.aN(0) << std::endl;
                }
            } while (fiter <= VAMBE.NLS);
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
        const int nbeams = 1;

        NLEBBE3D EBBE3D = { NLEBBE3D::NLEBBE3D(1) };

        Eigen::VectorXd dGU = Eigen::VectorXd::Zero(EBBE3D.get_ndof() * EBBE3D.get_nnode());
        Eigen::VectorXd GU = Eigen::VectorXd::Zero(EBBE3D.get_ndof() * EBBE3D.get_nnode());
        Eigen::VectorXd error = Eigen::VectorXd::Zero(EBBE3D.get_ndof() * EBBE3D.get_nnode());
        Eigen::VectorXd GU_new = Eigen::VectorXd::Zero(EBBE3D.get_ndof() * EBBE3D.get_nnode());
        Eigen::VectorXd GR = Eigen::VectorXd::Zero(EBBE3D.get_ndof() * EBBE3D.get_nnode());
        Eigen::SparseMatrix<double, Eigen::ColMajor> GT(EBBE3D.get_ndof()* EBBE3D.get_nnode(), EBBE3D.get_ndof()* EBBE3D.get_nnode());

        double max;
        int iter = 1;
        int fiter = 1;
        int maxiter = 50;

        std::fstream file1, file2;
        file1.open("E:/Adhithya/MTech_Res_Thesis/Cpp/ThesisCode/TextileComposites/TextileComposites/Result_Log.txt", std::fstream::in | std::fstream::out);
        file2.open("E:/Adhithya/MTech_Res_Thesis/Cpp/ThesisCode/TextileComposites/TextileComposites/Results.txt", std::fstream::in | std::fstream::out);
    
        //declare and initialize material parameters and other input variables
        double *D = new double[7];
        D[0] = EBBE3D.get_matprop("E");
        D[1] = EBBE3D.get_matprop("nu");
        D[2] = EBBE3D.get_matprop("Bp");
        D[3] = EBBE3D.get_matprop("Hp");
        D[4] = EBBE3D.get_matprop("Zx");
        D[5] = EBBE3D.get_matprop("Zy");
        D[6] = EBBE3D.get_matprop("Zz");
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
                        double h = X[2][0] - X[0][0];

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

                        EBBE3D.RKt(D, X, U, T, R);

                        /*for (int i = 0; i < 12; i++)
                            std::cout << R[i] << std::endl;

                        for (int j = 0; j < 12; j++)
                        {
                            for (int k = 0; k < 12; k++)
                                std::cout << T[j][k] << "   ";
                            std::cout << std::endl;
                        }*/

                        //Assembly
                        for (int j = 0; j < EBBE3D.get_ndof() * EBBE3D.get_nen(); j++)
                        {
                            GR(6 * a + j) += R[j];
                            for (int k = 0; k < EBBE3D.get_ndof() * EBBE3D.get_nen(); k++)
                                GT.coeffRef(6 * a + j, 6 * a + k) += T[j][k];
                        }
                    }

                    //Assign Boundary conditions
                    for (int j = 0; j < 6; j++)
                    {
                        for (int k = 0; k < EBBE3D.get_nnode() * EBBE3D.get_ndof(); k++)
                            GT.coeffRef(j, k) = 0;
                        GR(j) = 0;
                        GT.coeffRef(j, j) = 1;
                    }
                    GR(EBBE3D.get_nnode()*EBBE3D.get_ndof() - 4) -= load;
                    
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
                    VTKGrid VTK1 = VTKGrid(EBBE3D, 8, load, "EBBE3D", GU, 3, 8, EBBE3D.get_ndof(), 4);
                    fiter++;
                    load += 4.167;
                    /*for (int j = 0; j < EBBE3D.get_nnode() * EBBE3D.get_ndof(); j++)
                    {
                        //file2 << "Displacements of node " << j + 1 << std::endl;
                        //for (int k = 0; k < 3; k++)
                            //file2 << U(12 * (j - 1) + 6 + k) << std::endl;
                        file2 << GU(j) << std::endl;
                    }*/
                    file2 << GU(EBBE3D.get_nnode() * EBBE3D.get_ndof() - 4) << std::endl;
                }
                else
                {
                    std::cout << "Code didn't converge" << std::endl;
                    break;
                }
            } while (fiter < EBBE3D.get_nls());
        }
    }
    //--------------------------------------------------------//
    //----------3D EulerBernouli Beam Element-----------------//
    //-----------------Large displacement---------------------//
    //-----------------Updated Lagrangian---------------------//
    //--------------------Quadratic---------------------------//
    //---------------Node to Node Contact---------------------//
   /* else if (choice == 6)
    {
        //Master 1
        //Slave 2
        double ms;
        //Read master data
        ms = 1;
        NLEBBE3D EBBE3D1 = ReadEBBE3DElement(ms);
        //Read slave data
        ms = 2;
        NLEBBE3D EBBE3D2 = ReadEBBE3DElement(ms);

        int** ContactPairs = new int* [EBBE3D2.get_nnode()];
        for (int j = 0; j < EBBE3D2.get_nnode(); j++)
            ContactPairs[j] = new int[EBBE3D1.get_nnode()];

        Eigen::VectorXd dGU = Eigen::VectorXd::Zero(EBBE3D1.get_ndof() * (EBBE3D1.get_nnode() + EBBE3D2.get_nnode()));
        Eigen::VectorXd GU = Eigen::VectorXd::Zero(EBBE3D1.get_ndof() * (EBBE3D1.get_nnode() + EBBE3D2.get_nnode()));
        Eigen::VectorXd error = Eigen::VectorXd::Zero(EBBE3D1.get_ndof() * (EBBE3D1.get_nnode() + EBBE3D2.get_nnode()));
        Eigen::VectorXd GU_new = Eigen::VectorXd::Zero(EBBE3D1.get_ndof() * (EBBE3D1.get_nnode() + EBBE3D2.get_nnode()));
        Eigen::VectorXd GR = Eigen::VectorXd::Zero(EBBE3D1.get_ndof() * (EBBE3D1.get_nnode() + EBBE3D2.get_nnode()));
        Eigen::SparseMatrix<double, Eigen::ColMajor> GT(EBBE3D1.get_ndof() * (EBBE3D1.get_nnode() + EBBE3D2.get_nnode()), EBBE3D1.get_ndof() * (EBBE3D1.get_nnode() + EBBE3D2.get_nnode()));

        double max;
        int iter = 1;
        int fiter = 1;
        int maxiter = 50;

        std::fstream file1, file2;
        file1.open("E:/Adhithya/MTech_Res_Thesis/Cpp/ThesisCode/TextileComposites/TextileComposites/Result_Log.txt", std::fstream::in | std::fstream::out);
        file2.open("E:/Adhithya/MTech_Res_Thesis/Cpp/ThesisCode/TextileComposites/TextileComposites/Results.txt", std::fstream::in | std::fstream::out);

        //declare and initialize material parameters and other input variables
        //master
        double* D1 = new double[7];
        D1[0] = EBBE3D1.get_matprop("E");
        D1[1] = EBBE3D1.get_matprop("nu");
        D1[2] = EBBE3D1.get_matprop("Bp");
        D1[3] = EBBE3D1.get_matprop("Hp");
        D1[4] = EBBE3D1.get_matprop("Zx");
        D1[5] = EBBE3D1.get_matprop("Zy");
        D1[6] = EBBE3D1.get_matprop("Zz");

        //slave
        double* D2 = new double[7];
        D2[0] = EBBE3D2.get_matprop("E");
        D2[1] = EBBE3D2.get_matprop("nu");
        D2[2] = EBBE3D2.get_matprop("Bp");
        D2[3] = EBBE3D2.get_matprop("Hp");
        D2[4] = EBBE3D2.get_matprop("Zx");
        D2[5] = EBBE3D2.get_matprop("Zy");
        D2[6] = EBBE3D2.get_matprop("Zz");

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
                    T = new double* [EBBE3D1.get_ndof() * EBBE3D1.get_nen()];
                    for (int i = 0; i < EBBE3D1.get_ndof() * EBBE3D1.get_nen(); i++)
                        T[i] = new double[EBBE3D1.get_ndof() * EBBE3D1.get_nen()];
                    double* R = new double[EBBE3D1.get_ndof() * EBBE3D1.get_nen()];

                    //declare local tangent stiffness matrix and residual vector for contact
                    double** CT;
                    CT = new double* [EBBE3D1.get_ndof() * 2];
                    for (int i = 0; i < EBBE3D1.get_ndof() * 2; i++)
                        CT[i] = new double[EBBE3D1.get_ndof() * 2];
                    double* CR = new double[EBBE3D1.get_ndof() * 2];


                    for (int i = 0; i < (EBBE3D1.get_nelem() + EBBE3D2.get_nelem()); i++)
                    {
                        //master
                        if (i < EBBE3D1.get_nelem())
                        {
                            int a = (int)EBBE3D1.get_connectivity(i, 1) - 1;
                            int b = (int)EBBE3D1.get_connectivity(i, 2) - 1;
                            int c = (int)EBBE3D1.get_connectivity(i, 3) - 1;

                            //declare and initialize position vectors
                            /*double* X[2];
                            for (int j = 0; j < 2; j++)
                                X[j] = new double[3];
                            double X[3][3];

                            for (int j = 0; j < 3; j++)
                            {
                                X[0][j] = EBBE3D1.get_coordinates(a, j);
                                X[1][j] = EBBE3D1.get_coordinates(b, j);
                                X[2][j] = EBBE3D1.get_coordinates(c, j);
                            }
                            double h = X[2][0] - X[0][0];

                            //declare and initialize dofs
                            double U[3][6];

                            for (int j = 0; j < 6; j++)
                            {
                                U[0][j] = GU(6 * a + j);
                                U[1][j] = GU(6 * b + j);
                                U[2][j] = GU(6 * c + j);
                            }


                            for (int j = 0; j < EBBE3D1.get_ndof() * EBBE3D1.get_nen(); j++)
                            {
                                R[j] = 0;
                                for (int k = 0; k < EBBE3D1.get_ndof() * EBBE3D1.get_nen(); k++)
                                    T[j][k] = 0;
                            }

                            //Local Residual and tangent matrix
                            RKt(D1, X, U, T, R);

                            /*for (int i = 0; i < 12; i++)
                                std::cout << R[i] << std::endl;

                            for (int j = 0; j < 12; j++)
                            {
                                for (int k = 0; k < 12; k++)
                                    std::cout << T[j][k] << "   ";
                                std::cout << std::endl;
                            }

                            //Assembly
                            for (int j = 0; j < EBBE3D1.get_ndof() * EBBE3D1.get_nen(); j++)
                            {
                                GR(6 * a + j) += R[j];
                                for (int k = 0; k < EBBE3D1.get_ndof() * EBBE3D1.get_nen(); k++)
                                    GT.coeffRef(6 * a + j, 6 * a + k) += T[j][k];
                            }
                        }
                        //slave
                        else
                        {
                            int a = (int)EBBE3D2.get_connectivity(i - EBBE3D1.get_nelem(), 1) - 1;
                            int b = (int)EBBE3D2.get_connectivity(i - EBBE3D1.get_nelem(), 2) - 1;
                            int c = (int)EBBE3D2.get_connectivity(i - EBBE3D1.get_nelem(), 3) - 1;

                            //std::cout << a << std::endl;
                            //std::cout << b << std::endl;
                            //std::cout << c << std::endl;

                            //declare and initialize position vectors
                            double X[3][3];

                            for (int j = 0; j < 3; j++)
                            {
                                X[0][j] = EBBE3D2.get_coordinates(a, j);
                                //std::cout << X[0][j] << std::endl;
                                X[1][j] = EBBE3D2.get_coordinates(b, j);
                                //std::cout << X[1][j] << std::endl;
                                X[2][j] = EBBE3D2.get_coordinates(c, j);
                                //std::cout << X[2][j] << std::endl;
                            }
                            double h = X[2][0] - X[0][0];

                            //declare and initialize dofs
                            double U[3][6];

                            for (int j = 0; j < 6; j++)
                            {
                                U[0][j] = GU(6 * (a + EBBE3D1.get_nnode()) + j);
                                U[1][j] = GU(6 * (b + EBBE3D1.get_nnode()) + j);
                                U[2][j] = GU(6 * (c + EBBE3D1.get_nnode()) + j);
                            }

                            for (int j = 0; j < EBBE3D2.get_ndof() * EBBE3D2.get_nen(); j++)
                            {
                                R[j] = 0;
                                for (int k = 0; k < EBBE3D2.get_ndof() * EBBE3D2.get_nen(); k++)
                                    T[j][k] = 0;
                            }

                            //Local Residual and tangent matrix
                            RKt(D2, X, U, T, R);

                            /*for (int i = 0; i < 12; i++)
                                std::cout << R[i] << std::endl;

                            for (int j = 0; j < 12; j++)
                            {
                                for (int k = 0; k < 12; k++)
                                    std::cout << T[j][k] << "   ";
                                std::cout << std::endl;
                            }

                            //Assembly
                            for (int j = 0; j < EBBE3D2.get_ndof() * EBBE3D2.get_nen(); j++)
                            {
                                GR(6 * a + j + EBBE3D1.get_nnode() * EBBE3D1.get_ndof()) += R[j];
                                for (int k = 0; k < EBBE3D2.get_ndof() * EBBE3D2.get_nen(); k++)
                                    GT.coeffRef(6 * a + j + EBBE3D1.get_nnode() * EBBE3D1.get_ndof(), 6 * a + k + EBBE3D1.get_nnode() * EBBE3D1.get_ndof()) += T[j][k];
                            }

                        }
                    }

                    //Contact Search 
                    ContactSearch(EBBE3D1, EBBE3D2, ContactPairs,"NTN");
                    //Apply contact constraint
                    //EBBE3D2 - slave
                    //EBBE3D1 - master
                    for (int i = 0; i < EBBE3D2.get_nnode(); i++)
                    {
                        //Find the corresponding master node for the slave node.
                        int mnode = ContactPairs[i][1] - 1;
                        int snode = ContactPairs[i][0] - 1;

                        //std::cout << mnode << std::endl;
                        //std::cout << snode << std::endl;

                        double X1[3], X2[3], u1[6], u2[6], D[4];

                        D[0] = 1e10;
                        D[1] = 0;
                        D[2] = 1;
                        D[3] = 0;

                        for (int j = 0; j < 3; j++)
                        {
                            X1[j] = EBBE3D2.get_coordinates(mnode, j);
                            X2[j] = EBBE3D1.get_coordinates(snode, j);
                        }

                        /*for (int j = 0; j < 3; j++)
                            std::cout << X1[j] << std::endl;

                        for (int j = 0; j < 3; j++)
                            std::cout << X2[j] << std::endl;

                        for (int j = 0; j < 6; j++)
                        {
                            u1[j] = GU.coeffRef(EBBE3D2.get_ndof() * (snode + EBBE3D1.get_nnode()) + j);
                            u2[j] = GU.coeffRef(EBBE3D1.get_ndof() * mnode + j);
                        }

                        //Initialize contact residual and contact tangent matrix
                        for (int j = 0; j < EBBE3D1.get_ndof() * 2; j++)
                        {
                            CR[j] = 0;
                            for (int k = 0; k < EBBE3D1.get_ndof() * 2; k++)
                                CT[j][k] = 0;
                        }
                        
                        double g;
                        Contact_NTN(D, X1, X2, u1, u2, CR, CT, &g);

                        /*std::cout << "Tangent matrix" << std::endl;
                        for (int j = 0; j < EBBE3D1.get_ndof() * 2; j++)
                        {
                            for (int k = 0; k < EBBE3D1.get_ndof() * 2; k++)
                                std::cout << CT[j][k] << "  ";
                            std::cout << std::endl;
                        }*/

                        /*std::cout << "Residual vector" << std::endl;
                        for (int j = 0; j < EBBE3D1.get_ndof() * 2; j++)
                            std::cout << CR[j] << std::endl;

                        //Assembly
                        std::cout << g << std::endl;
                        if (g < 0)
                        {
                            for (int j = 0; j < 6; j++)
                            {
                                GR(6 * mnode + j) += CR[j + 6];
                                GR(6 * (snode + EBBE3D1.get_nnode()) + j) += CR[j];
                                for (int k = 0; k < 6; k++)
                                {
                                    GT.coeffRef(6 * (snode + EBBE3D1.get_nnode()) + j, 6 * (snode + EBBE3D1.get_nnode()) + k) += CT[j][k];
                                    GT.coeffRef(6 * (snode + EBBE3D1.get_nnode()) + j, 6 * mnode + k) += CT[j][k + 6];
                                    GT.coeffRef(6 * mnode + j, 6 * (snode + EBBE3D1.get_nnode()) + k) += CT[j + 6][k];
                                    GT.coeffRef(6 * mnode + j, 6 * mnode + k) += CT[j + 6][k + 6];
                                }
                            }
                        }
                    }

                    //Assign Boundary conditions
                    //master
                    //fix zero end
                    for (int j = 0; j < 6; j++)
                    {
                        for (int k = 0; k < EBBE3D1.get_nnode() * EBBE3D1.get_ndof(); k++)
                            GT.coeffRef(j, k) = 0;
                        GR(j) = 0;
                        GT.coeffRef(j, j) = 1;
                    }
                    //fix last end
                    for (int j = (EBBE3D1.get_nnode() - 1) * EBBE3D1.get_ndof(); j < EBBE3D1.get_nnode() * EBBE3D1.get_ndof(); j++)
                    {
                        for (int k = 0; k < EBBE3D1.get_nnode() * EBBE3D1.get_ndof(); k++)
                            GT.coeffRef(j, k) = 0;
                        GR(j) = 0;
                        GT.coeffRef(j, j) = 1;
                    }

                    int p = EBBE3D1.get_loadnode() - 1;
                    GR(EBBE3D1.get_ndof() * p + 1) -= load;
                    //slave
                    for (int j = EBBE3D1.get_nnode() * EBBE3D1.get_ndof(); j < EBBE3D1.get_nnode() * EBBE3D1.get_ndof() + 6; j++)
                    {
                        for (int k = EBBE3D1.get_nnode() * EBBE3D1.get_ndof(); k < EBBE3D1.get_nnode() * EBBE3D1.get_ndof() + EBBE3D2.get_nnode() * EBBE3D2.get_ndof(); k++)
                            GT.coeffRef(j, k) = 0;
                        GR(j) = 0;
                        GT.coeffRef(j, j) = 1;
                    }
                    for (int j = EBBE3D1.get_nnode() * EBBE3D1.get_ndof() + (EBBE3D2.get_nnode() - 1) * EBBE3D2.get_ndof(); j < EBBE3D1.get_nnode() * EBBE3D1.get_ndof() + EBBE3D2.get_nnode() * EBBE3D2.get_ndof(); j++)
                    {
                        for (int k = EBBE3D1.get_nnode() * EBBE3D1.get_ndof(); k < EBBE3D1.get_nnode() * EBBE3D1.get_ndof() + EBBE3D2.get_nnode() * EBBE3D2.get_ndof(); k++)
                            GT.coeffRef(j, k) = 0;
                        GR(j) = 0;
                        GT.coeffRef(j, j) = 1;
                    }
                    int q = EBBE3D2.get_loadnode() - 1;
                    GR(EBBE3D1.get_nnode() * EBBE3D1.get_ndof() + q * EBBE3D2.get_ndof() + 1) -= -load;

                    //Solve the equation
                    GT.makeCompressed();

                    Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>> solver;
                    solver.analyzePattern(GT);
                    solver.factorize(GT); //LU decomposition

                    assert(solver.info() == Eigen::Success);

                    dGU = solver.solve(-GR);

                    //Calculate errors and update for next iteration
                    for (int j = 0; j < (EBBE3D1.get_nnode() + EBBE3D2.get_nnode()) * EBBE3D1.get_ndof(); j++)
                        GU_new(j) = GU(j) + dGU(j);

                    //Error calculation
                    for (int j = 0; j < (EBBE3D1.get_nnode() + EBBE3D2.get_nnode()) * EBBE3D1.get_ndof(); j++)
                        error(j) = abs(GU_new(j) - GU(j));

                    max = error(0);
                    for (int j = 0; j < (EBBE3D1.get_nnode() + EBBE3D2.get_nnode()) * EBBE3D1.get_ndof(); j++)
                        if (max < error(j))
                            max = error(j);

                    //Assignment for next iteration
                    for (int j = 0; j < (EBBE3D1.get_nnode() + EBBE3D2.get_nnode()) * EBBE3D1.get_ndof(); j++)
                        GU(j) = GU_new(j);
                    iter++;

                    //Free memory
                    for (int i = 0; i < EBBE3D1.get_ndof() * EBBE3D1.get_nen(); i++)
                        delete T[i];
                    delete[] R;

                    //std::cout << max << std::endl;
                } while (max > pow(10, -6) && iter < maxiter);
                if (iter < maxiter)
                {
                    fiter++;
                    /*file2 << "Converged values at load  " << load << std::endl;
                    file2 << "Values for beam 1 " << std::endl;
                    for (int j = 0; j < EBBE3D1.get_nnode() * EBBE3D1.get_ndof(); j++)
                        file2 << GU(j) << std::endl;
                    file2 << "Values for beam 2 " << std::endl;
                    for (int j = EBBE3D1.get_nnode() * EBBE3D1.get_ndof(); j < EBBE3D1.get_nnode() * EBBE3D1.get_ndof() + EBBE3D2.get_nnode() * EBBE3D2.get_ndof(); j++)
                        file2 << GU(j) << std::endl;
                    Eigen::VectorXd GU1 = Eigen::VectorXd::Zero(EBBE3D1.get_nnode() * EBBE3D1.get_ndof());
                    Eigen::VectorXd GU2 = Eigen::VectorXd::Zero(EBBE3D2.get_nnode() * EBBE3D2.get_ndof());
                    for (int j = 0; j < EBBE3D1.get_nnode() * EBBE3D1.get_ndof(); j++)
                        GU1(j) = GU(j);
                    for (int j = EBBE3D1.get_nnode() * EBBE3D1.get_ndof(); j < EBBE3D1.get_nnode() * EBBE3D1.get_ndof() + EBBE3D2.get_nnode() * EBBE3D2.get_ndof(); j++)
                        GU2(j - EBBE3D1.get_nnode() * EBBE3D1.get_ndof()) = GU(j);
                    PostProcessing(EBBE3D1.get_coordinates, GU1, load, "EBBE3D1", EBBE3D1.get_nnode(), 3, EBBE3D1.get_ndof());
                    PostProcessing(EBBE3D2.get_coordinates, GU2, load, "EBBE3D2", EBBE3D2.get_nnode(), 3, EBBE3D2.get_ndof());
                    load += 20;
                    std::cout << load << std::endl;
                    //Print the current configuration in a file
                }
                else
                {
                    std::cout << "Code didn't converge" << std::endl;
                    break;
                }
            } while (fiter < EBBE3D1.get_nls());
        }
    }
    //--------------------------------------------------------//
    //----------3D EulerBernouli Beam Element-----------------//
    //-----------------Large displacement---------------------//
    //-----------------Updated Lagrangian---------------------//
    //--------------------Quadratic---------------------------//
    //------------Segment to Segment Contact------------------//
    else if (choice == 7)
    {
        //Master 1
        //Slave 2
        const int nbeams = 2;

        NLEBBE3D EBBE3D[nbeams] = { NLEBBE3D::NLEBBE3D(0), NLEBBE3D::NLEBBE3D(1) };

        int size = 0;
        for (int i = 0; i < nbeams; i++)
            size += EBBE3D[i].get_ndof() * EBBE3D[i].get_nnode();

        Eigen::VectorXd dGU = Eigen::VectorXd::Zero(size);
        Eigen::VectorXd GU = Eigen::VectorXd::Zero(size);
        Eigen::VectorXd error = Eigen::VectorXd::Zero(size);
        Eigen::VectorXd GU_new = Eigen::VectorXd::Zero(size);
        Eigen::VectorXd GR = Eigen::VectorXd::Zero(size);
        Eigen::SparseMatrix<double, Eigen::ColMajor> GT(size, size);

        double max = 0;
        int iter = 1;
        int fiter = 1;
        int maxiter = 50;

        std::fstream file1, file2;
        file1.open("E:/Adhithya/MTech_Res_Thesis/Cpp/ThesisCode/TextileComposites/TextileComposites/Result_Log.txt", std::fstream::in | std::fstream::out);
        file2.open("E:/Adhithya/MTech_Res_Thesis/Cpp/ThesisCode/TextileComposites/TextileComposites/Results.txt", std::fstream::in | std::fstream::out);

        //declare and initialize material parameters and other input variables
        //master
        //std::unique_ptr<double[]> D1 = std::make_unique<double[]>(7);

        double load = 0;

        if (file1.is_open() && file2.is_open())
        {
            for (int i = 0; i < nbeams; i++)
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

                        //declare element tangent stiffness matrix and residual vector
                        double** T;
                        T = new double* [EBBE3D[i].get_ndof() * EBBE3D[i].get_nen()];
                        for (int i = 0; i < EBBE3D[i].get_ndof() * EBBE3D[i].get_nen(); i++)
                            T[i] = new double[EBBE3D[i].get_ndof() * EBBE3D[i].get_nen()];
                        double* R = new double[EBBE3D[i].get_ndof() * EBBE3D[i].get_nen()];

                        //initialize element tangent stiffness and element residual to zero. 
                        for (int i = 0; i < EBBE3D[i].get_ndof() * EBBE3D[i].get_nen(); i++)
                            R[i] = 0;
                        for (int i = 0; i < EBBE3D[i].get_ndof() * EBBE3D[i].get_nen(); i++)
                            for (int j = 0; j < EBBE3D[i].get_ndof() * EBBE3D[i].get_nen(); j++)
                                T[i][j] = 0;

                        for (int i = 0; i < EBBE3D[i].get_nelem(); i++)
                        {
                            double* D1 = new double[7];
                            D1[0] = EBBE3D[i].get_matprop("E");
                            D1[1] = EBBE3D[i].get_matprop("nu");
                            D1[2] = EBBE3D[i].get_matprop("Bp");
                            D1[3] = EBBE3D[i].get_matprop("Hp");
                            D1[4] = EBBE3D[i].get_matprop("Zx");
                            D1[5] = EBBE3D[i].get_matprop("Zy");
                            D1[6] = EBBE3D[i].get_matprop("Zz");

                            //master
                            if (i == 0)
                            {
                                int a = (int)EBBE3D[i].get_connectivity(i, 1) - 1;
                                int b = (int)EBBE3D[i].get_connectivity(i, 2) - 1;
                                int c = (int)EBBE3D[i].get_connectivity(i, 3) - 1;

                                double X[3][3];

                                for (int j = 0; j < 3; j++)
                                {
                                    X[0][j] = EBBE3D[i].get_coordinates(a, j);
                                    X[1][j] = EBBE3D[i].get_coordinates(b, j);
                                    X[2][j] = EBBE3D[i].get_coordinates(c, j);
                                }
                                double h = X[2][0] - X[0][0];

                                //declare and initialize dofs
                                double U[3][6];

                                for (int j = 0; j < 6; j++)
                                {
                                    U[0][j] = GU(6 * a + j);
                                    U[1][j] = GU(6 * b + j);
                                    U[2][j] = GU(6 * c + j);
                                }


                                for (int j = 0; j < EBBE3D[i].get_ndof() * EBBE3D[i].get_nen(); j++)
                                {
                                    R[j] = 0;
                                    for (int k = 0; k < EBBE3D[i].get_ndof() * EBBE3D[i].get_nen(); k++)
                                        T[j][k] = 0;
                                }

                                //Local Residual and tangent matrix
                                EBBE3D[i].RKt(D1, X, U, T, R);

                                //Assembly
                                for (int j = 0; j < EBBE3D[i].get_ndof() * EBBE3D[i].get_nen(); j++)
                                {
                                    GR(6 * a + j) += R[j];
                                    for (int k = 0; k < EBBE3D[i].get_ndof() * EBBE3D[i].get_nen(); k++)
                                        GT.coeffRef(6 * a + j, 6 * a + k) += T[j][k];
                                }
                            }
                            //slave
                            else
                            {
                                int a = (int)EBBE3D[i].get_connectivity(i - EBBE3D[i - 1].get_nelem(), 1) - 1;
                                int b = (int)EBBE3D[i].get_connectivity(i - EBBE3D[i - 1].get_nelem(), 2) - 1;
                                int c = (int)EBBE3D[i].get_connectivity(i - EBBE3D[i - 1].get_nelem(), 3) - 1;

                                //declare and initialize position vectors
                                double X[3][3];

                                for (int j = 0; j < 3; j++)
                                {
                                    X[0][j] = EBBE3D[i].get_coordinates(a, j);
                                    //std::cout << X[0][j] << std::endl;
                                    X[1][j] = EBBE3D[i].get_coordinates(b, j);
                                    //std::cout << X[1][j] << std::endl;
                                    X[2][j] = EBBE3D[i].get_coordinates(c, j);
                                    //std::cout << X[2][j] << std::endl;
                                }
                                double h = X[2][0] - X[0][0];

                                //declare and initialize dofs
                                double U[3][6];

                                for (int j = 0; j < 6; j++)
                                {
                                    U[0][j] = GU(6 * (a + EBBE3D[i - 1].get_nnode()) + j);
                                    U[1][j] = GU(6 * (b + EBBE3D[i - 1].get_nnode()) + j);
                                    U[2][j] = GU(6 * (c + EBBE3D[i - 1].get_nnode()) + j);
                                }

                                for (int j = 0; j < EBBE3D[i].get_ndof() * EBBE3D[i].get_nen(); j++)
                                {
                                    R[j] = 0;
                                    for (int k = 0; k < EBBE3D[i].get_ndof() * EBBE3D[i].get_nen(); k++)
                                        T[j][k] = 0;
                                }

                                //Local Residual and tangent matrix
                                EBBE3D[i].RKt(D1, X, U, T, R);

                                //Assembly
                                for (int j = 0; j < EBBE3D[i].get_ndof() * EBBE3D[i].get_nen(); j++)
                                {
                                    GR(6 * a + j + EBBE3D[i - 1].get_nnode() * EBBE3D[i - 1].get_ndof()) += R[j];
                                    for (int k = 0; k < EBBE3D[i].get_ndof() * EBBE3D[i].get_nen(); k++)
                                        GT.coeffRef(6 * a + j + EBBE3D[i - 1].get_nnode() * EBBE3D[i - 1].get_ndof(), 6 * a + k + EBBE3D[i - 1].get_nnode() * EBBE3D[i - 1].get_ndof()) += T[j][k];
                                }

                            }
                        }

                        //-------------------------------------------------//
                        //-------------------------------------------------//
                        //--------------------Contact----------------------//
                        //-------------------------------------------------//
                        //-------------------------------------------------//
                        BeamContact BCon{ nbeams, "STS"};

                        //Global Contact Search
                        BCon.GlobalContactSearch(EBBE3D);
                        //Local Contact Search 
                        BCon.LocalContactSearch(EBBE3D);

                        //declare contact element tangent stiffness matrix and residual vector for contact
                        //Beam elements with different degrees of freedom cannot be used.
                        double** CT;
                        CT = new double* [EBBE3D[0].get_ndof() * BCon.get_nen()];
                        for (int i = 0; i < EBBE3D[0].get_ndof() * BCon.get_nen(); i++)
                            CT[i] = new double[EBBE3D[0].get_ndof() * BCon.get_nen()];
                        double* CR = new double[EBBE3D[0].get_ndof() * BCon.get_nen()];

                        const Eigen::MatrixXd globalconn = BCon.get_Globalelem();
                        double** TP;
                        TP = new double* [2];
                        for (int i = 0; i < 2; i++)
                            TP[i] = new double[2];
                        double* RP = new double[2];

                        for (int i = 0; i < globalconn.rows(); i++)
                        {
                            int* b1 = new int [globalconn.cols()];
                            for (int j = 0; j < globalconn.cols(); j++)
                                b1[j] = globalconn(i, j);

                            for (int j = 0; j < EBBE3D[b1[0]].get_nelem(); j++)
                            {
                                //Find the corresponding master element for the slave element.
                                int melem = BCon.getContactPairs[i][1] - 1;
                                int selem = ContactPairs[i][0] - 1;
                            }

                            
                        }
                       
                        

                        

                        //-------------------------------------------------//
                        //------------Apply contact constraint-------------//
                        //-------------------------------------------------// 
                        //EBBE3D2 - slave
                        //EBBE3D1 - master
                        for (int i = 0; i < EBBE3D2.get_nelem(); i++)
                        {
                            //Find the corresponding master element for the slave element.
                            int melem = ContactPairs[i][1] - 1;
                            int selem = ContactPairs[i][0] - 1;

                            //Find the nodes for slave and master elements
                            int n1, n2, n3, m1, m2, m3;
                            n1 = EBBE3D1.get_connectivity(melem, 1) - 1;
                            n2 = EBBE3D1.get_connectivity(melem, 2) - 1;
                            n3 = EBBE3D1.get_connectivity(melem, 3) - 1;
                            m1 = EBBE3D2.get_connectivity(selem, 1) - 1;
                            m2 = EBBE3D2.get_connectivity(selem, 2) - 1;
                            m3 = EBBE3D2.get_connectivity(selem, 3) - 1;

                            double X1[3], X2[3], X3[3], u1[6], u2[6], u3[6];

                            double Y1[3], Y2[3], Y3[3], v1[6], v2[6], v3[6];

                            for (int j = 0; j < 3; j++)
                            {
                                X1[j] = EBBE3D1.get_coordinates(n1, j);
                                X2[j] = EBBE3D1.get_coordinates(n2, j);
                                X3[j] = EBBE3D1.get_coordinates(n3, j);
                                Y1[j] = EBBE3D2.get_coordinates(m1, j);
                                Y2[j] = EBBE3D2.get_coordinates(m2, j);
                                Y3[j] = EBBE3D2.get_coordinates(m3, j);
                            }

                            for (int j = 0; j < 6; j++)
                            {
                                u1[j] = GU.coeffRef(EBBE3D1.get_ndof() * n1 + j);
                                u2[j] = GU.coeffRef(EBBE3D1.get_ndof() * n2 + j);
                                u3[j] = GU.coeffRef(EBBE3D1.get_ndof() * n3 + j);

                                v1[j] = GU.coeffRef(EBBE3D2.get_ndof() * (m1 + EBBE3D1.get_nnode()) + j);
                                v2[j] = GU.coeffRef(EBBE3D2.get_ndof() * (m2 + EBBE3D1.get_nnode()) + j);
                                v3[j] = GU.coeffRef(EBBE3D2.get_ndof() * (m3 + EBBE3D1.get_nnode()) + j);
                            }

                            double Data2[4];
                            Data2[0] = EBBE3D1.epsilon;
                            Data2[1] = EBBE3D1.D;
                            Data2[2] = EBBE3D2.D;
                            Data2[3] = EBBE3D1.get_nnode();
                            //segment n1 - n2 and m1 - m2
                            Contact_STS(Data2, X1, X2, Y1, Y2, u1, u2, v1, v2, &GT, &GR, n1, n2, m1, m2);

                            //segment n2 - n3 and m2 - m3
                            Contact_STS(Data2, X2, X3, Y2, Y3, u2, u3, v2, v3, &GT, &GR, n2, n3, m2, m3);

                            double Data3[5];
                            Data3[0] = EBBE3D1.epsilon;
                            Data3[1] = EBBE3D1.D;
                            Data3[2] = EBBE3D2.D;
                            Data3[3] = EBBE3D1.get_nnode();
                            Data3[4] = EBBE3D2.get_nnode();

                            //Special case for contact between endpoints 
                            //Enforce constraint for closest enpoints
                            //closest points n1, m1
                            Contact_STS_Endpoints(Data3, X1, Y1, u1, v1, &GT, &GR, n1, m1);

                            //closest points n2, m2
                            Contact_STS_Endpoints(Data3, X2, Y2, u2, v2, &GT, &GR, n2, m2);

                            //closest points n2, m2
                            Contact_STS_Endpoints(Data3, X3, Y3, u3, v3, &GT, &GR, n3, m3);
                        }

                        //Assign Boundary conditions
                        //master
                        //fix zero end
                        for (int j = 0; j < 6; j++)
                        {
                            for (int k = 0; k < EBBE3D1.get_nnode() * EBBE3D1.get_ndof(); k++)
                                GT.coeffRef(j, k) = 0;
                            GR(j) = 0;
                            GT.coeffRef(j, j) = 1;
                        }
                        //fix last end
                        for (int j = (EBBE3D1.get_nnode() - 1) * EBBE3D1.get_ndof(); j < EBBE3D1.get_nnode() * EBBE3D1.get_ndof(); j++)
                        {
                            for (int k = 0; k < EBBE3D1.get_nnode() * EBBE3D1.get_ndof(); k++)
                                GT.coeffRef(j, k) = 0;
                            GR(j) = 0;
                            GT.coeffRef(j, j) = 1;
                        }

                        int p = EBBE3D1.load - 1;
                        GR(EBBE3D1.get_ndof() * p + 1) -= load;
                        //slave
                        for (int j = EBBE3D1.get_nnode() * EBBE3D1.get_ndof(); j < EBBE3D1.get_nnode() * EBBE3D1.get_ndof() + 6; j++)
                        {
                            for (int k = EBBE3D1.get_nnode() * EBBE3D1.get_ndof(); k < EBBE3D1.get_nnode() * EBBE3D1.get_ndof() + EBBE3D2.get_nnode() * EBBE3D2.get_ndof(); k++)
                                GT.coeffRef(j, k) = 0;
                            GR(j) = 0;
                            GT.coeffRef(j, j) = 1;
                        }
                        for (int j = EBBE3D1.get_nnode() * EBBE3D1.get_ndof() + (EBBE3D2.get_nnode() - 1) * EBBE3D2.get_ndof(); j < EBBE3D1.get_nnode() * EBBE3D1.get_ndof() + EBBE3D2.get_nnode() * EBBE3D2.get_ndof(); j++)
                        {
                            for (int k = EBBE3D1.get_nnode() * EBBE3D1.get_ndof(); k < EBBE3D1.get_nnode() * EBBE3D1.get_ndof() + EBBE3D2.get_nnode() * EBBE3D2.get_ndof(); k++)
                                GT.coeffRef(j, k) = 0;
                            GR(j) = 0;
                            GT.coeffRef(j, j) = 1;
                        }
                        int q = EBBE3D2.get_load() - 1;
                        GR(EBBE3D1.get_nnode() * EBBE3D1.get_ndof() + q * EBBE3D2.get_ndof() + 1) -= -load;

                        //Solve the equation
                        GT.makeCompressed();

                        Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>> solver;
                        solver.analyzePattern(GT);
                        solver.factorize(GT); //LU decomposition

                        assert(solver.info() == Eigen::Success);

                        dGU = solver.solve(-GR);

                        //control step size
                        max = dGU(0);
                        for (int j = 1; j < (EBBE3D1.get_nnode() + EBBE3D2.get_nnode()); j++)
                        {
                            if (max < dGU(j))
                                max = dGU(j);
                            else if (max < dGU(j + 1))
                                max = dGU(j);
                            else if (max < dGU(j + 2))
                                max = dGU(j);
                            else
                                continue;
                        }
                        //std::cout << max << std::endl;
                        while (max > EBBE3D1.D || max > EBBE3D2.D)
                        {
                            for (int j = 0; j < (EBBE3D1.get_nnode() + EBBE3D2.get_nnode()); j++)
                            {
                                dGU(j) = 0.5 * dGU(j);
                                dGU(j + 1) = 0.5 * dGU(j + 1);
                                dGU(j + 2) = 0.5 * dGU(j + 2);
                            }

                            max = dGU(0);
                            for (int j = 1; j < (EBBE3D1.get_nnode() + EBBE3D2.get_nnode()); j++)
                            {
                                if (max < dGU(j))
                                    max = dGU(j);
                                else if (max < dGU(j + 1))
                                    max = dGU(j);
                                else if (max < dGU(j + 2))
                                    max = dGU(j);
                                else
                                    continue;
                            }
                        }
                        /*if (load >= 80)
                        {
                            std::cout << "Change in displacements" << std::endl;
                            for (int j = 0; j < (EBBE3D1.get_nnode() + EBBE3D2.get_nnode()) * EBBE3D1.get_ndof(); j++)
                                std::cout << dGU(j) << std::endl;
                            std::cout << "Tangent stiffness matrix" << std::endl;
                            for (int j = 0; j < EBBE3D1.get_nnode() * EBBE3D1.get_ndof() + EBBE3D2.get_nnode() * EBBE3D2.get_ndof(); j++)
                            {
                                for (int k = 0; k < EBBE3D1.get_nnode() * EBBE3D1.get_ndof() + EBBE3D2.get_nnode() * EBBE3D2.get_ndof(); k++)
                                    std::cout << GT.coeffRef(j, k) << " ";
                                std::cout << std::endl;
                            }
                            std::cout << "Residual vector" << std::endl;
                            for (int j = 0; j < EBBE3D1.get_nnode() * EBBE3D1.get_ndof() + EBBE3D2.get_nnode() * EBBE3D2.get_ndof(); j++)
                                std::cout << GR.coeffRef(j) << std::endl;
                        }

                        //Calculate errors and update for next iteration
                        for (int j = 0; j < (EBBE3D1.get_nnode() + EBBE3D2.get_nnode()) * EBBE3D1.get_ndof(); j++)
                            GU_new(j) = GU(j) + dGU(j);

                        //Error calculation
                        for (int j = 0; j < (EBBE3D1.get_nnode() + EBBE3D2.get_nnode()) * EBBE3D1.get_ndof(); j++)
                            error(j) = abs(GU_new(j) - GU(j));

                        max = error(0);
                        for (int j = 0; j < (EBBE3D1.get_nnode() + EBBE3D2.get_nnode()) * EBBE3D1.get_ndof(); j++)
                            if (max < error(j))
                                max = error(j);

                        //Assignment for next iteration
                        for (int j = 0; j < (EBBE3D1.get_nnode() + EBBE3D2.get_nnode()) * EBBE3D1.get_ndof(); j++)
                            GU(j) = GU_new(j);
                        iter++;

                        //Free memory
                        delete[] TP, RP, CT, T, R, CR, ContactPairs;

                        std::cout << max << std::endl;
                    } while (max > pow(10, -8) && iter < maxiter);
                    if (iter < maxiter)
                    {
                        fiter++;
                        /*file2 << "Converged values at load  " << load << std::endl;
                        file2 << "Values for beam 1 " << std::endl;
                        for (int j = 0; j < EBBE3D1.get_nnode() * EBBE3D1.get_ndof(); j++)
                            file2 << GU(j) << std::endl;
                        file2 << "Values for beam 2 " << std::endl;
                        for (int j = EBBE3D1.get_nnode() * EBBE3D1.get_ndof(); j < EBBE3D1.get_nnode() * EBBE3D1.get_ndof() + EBBE3D2.get_nnode() * EBBE3D2.get_ndof(); j++)
                            file2 << GU(j) << std::endl;
                            //for (int j = 0; j < (EBBE3D1.get_nnode() + EBBE3D2.get_nnode()) * EBBE3D1.get_ndof(); j++)
                            //    std::cout << GU(j) << std::endl;
                        Eigen::VectorXd GU1 = Eigen::VectorXd::Zero(EBBE3D1.get_nnode() * EBBE3D1.get_ndof());
                        Eigen::VectorXd GU2 = Eigen::VectorXd::Zero(EBBE3D2.get_nnode() * EBBE3D2.get_ndof());
                        for (int j = 0; j < EBBE3D1.get_nnode() * EBBE3D1.get_ndof(); j++)
                            GU1(j) = GU(j);
                        for (int j = EBBE3D1.get_nnode() * EBBE3D1.get_ndof(); j < EBBE3D1.get_nnode() * EBBE3D1.get_ndof() + EBBE3D2.get_nnode() * EBBE3D2.get_ndof(); j++)
                            GU2(j - EBBE3D1.get_nnode() * EBBE3D1.get_ndof()) = GU(j);
                        //PostProcessing(EBBE3D1.get_coordinates, GU1, load, "EBBE3D1", EBBE3D1.get_nnode(), 3, EBBE3D1.get_ndof());
                        //PostProcessing(EBBE3D2.get_coordinates, GU2, load, "EBBE3D2", EBBE3D2.get_nnode(), 3, EBBE3D2.get_ndof());
                        VTKGrid VTK1 = VTKGrid(EBBE3D1, 8, load, "EBBE3D1", GU1);
                        VTKGrid VTK2 = VTKGrid(EBBE3D2, 8, load, "EBBE3D2", GU2);
                        load += 10;
                        std::cout << load << std::endl;
                        //Print the current configuration in a file
                    }
                    else
                    {
                        Eigen::VectorXd GU1 = Eigen::VectorXd::Zero(EBBE3D1.get_nnode() * EBBE3D1.get_ndof());
                        Eigen::VectorXd GU2 = Eigen::VectorXd::Zero(EBBE3D2.get_nnode() * EBBE3D2.get_ndof());
                        for (int j = 0; j < EBBE3D1.get_nnode() * EBBE3D1.get_ndof(); j++)
                            GU1(j) = GU(j);
                        for (int j = EBBE3D1.get_nnode() * EBBE3D1.get_ndof(); j < EBBE3D1.get_nnode() * EBBE3D1.get_ndof() + EBBE3D2.get_nnode() * EBBE3D2.get_ndof(); j++)
                            GU2(j - EBBE3D1.get_nnode() * EBBE3D1.get_ndof()) = GU(j);
                        PostProcessing(EBBE3D1.get_coordinates, GU1, load, "EBBE3D1", EBBE3D1.get_nnode(), 3, EBBE3D1.get_ndof());
                        PostProcessing(EBBE3D2.get_coordinates, GU2, load, "EBBE3D2", EBBE3D2.get_nnode(), 3, EBBE3D2.get_ndof());
                        std::cout << "Code didn't converge" << std::endl;
                        break;
                    }
                } while (fiter < EBBE3D1.get_nls());
            }
        }
        delete[] D1, D2;
    }
    else
        std::cout << "Wrong option" << std::endl;*/
}

#endif




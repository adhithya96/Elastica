
#ifndef VARIABLES_H
#define VARIABLES_H
#include "Variables.h"
#include<fstream>
#include<limits>


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
        NonLinearEulerBernouliBeamElement NLEBBE = ReadNLEBBEFile();
        Eigen::VectorXd dU = Eigen::VectorXd::Zero(NLEBBE.NDOF * NLEBBE.NNODE);
        Eigen::VectorXd U = Eigen::VectorXd::Zero(NLEBBE.NDOF * NLEBBE.NNODE);
        Eigen::VectorXd error = Eigen::VectorXd::Zero(NLEBBE.NDOF * NLEBBE.NNODE);
        Eigen::MatrixXd k = Eigen::MatrixXd::Zero(6, 6);
        Eigen::MatrixXd t = Eigen::MatrixXd::Zero(6, 6);
        Eigen::VectorXd f = Eigen::VectorXd::Zero(6);
        Eigen::VectorXd F = Eigen::VectorXd::Zero(NLEBBE.NDOF * NLEBBE.NNODE);
        Eigen::VectorXd U_new = Eigen::VectorXd::Zero(NLEBBE.NDOF * NLEBBE.NNODE);
        Eigen::VectorXd R = Eigen::VectorXd::Zero(NLEBBE.NDOF * NLEBBE.NNODE);
        Eigen::SparseMatrix<double, Eigen::ColMajor> K(NLEBBE.NNODE * NLEBBE.NDOF, NLEBBE.NNODE * NLEBBE.NDOF),
            T(NLEBBE.NNODE * NLEBBE.NDOF, NLEBBE.NNODE * NLEBBE.NDOF);
        //Eigen::VectorXd X_ref = Eigen::VectorXd::Zero(NLEBBE.NDM * NLEBBE.NNODE);
        //Eigen::VectorXd X_curr = Eigen::VectorXd::Zero(NLEBBE.NDM * NLEBBE.NNODE);
     
        //Eigen::MatrixXd LOAD;
        double max = 1;
        int iter = 1;
        int fiter = 1;
        int maxiter = 20;
        //Initialize Reference and current configuration
        //for (int i = 0; i < (int)NLEBBE.NNODE; i++)
        //{
        //    X_ref(NLEBBE.NDM * i) = NLEBBE.NODE(i, 0);
        //    X_curr(NLEBBE.NDM* i) = NLEBBE.NODE(i, 0);
        //    X_ref(NLEBBE.NDM * i + 1) = NLEBBE.NODE(i, 1);
        //    X_curr(NLEBBE.NDM* i + 1) = NLEBBE.NODE(i, 1);
        //}
        //Loop for force iterations.
        do
        {
            //Newton Raphson loop
            //Update Reference configuration
            //for (int i = 0; i < (int)NLEBBE.NNODE; i++)
            //{
            //    X_ref(NLEBBE.NDM * i) += U(NLEBBE.NDOF * i);
            //    X_ref(NLEBBE.NDM * i + 1) += U(NLEBBE.NDOF * i + 1);
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
                

                for (int i = 0; i < (int)NLEBBE.NELEM; i++)
                {
                    int a = (int)NLEBBE.ELEM(i, 1) - 1;
                    int b = (int)NLEBBE.ELEM(i, 2) - 1;
                    double xa = NLEBBE.NODE((int)a, 0);
                    double xb = NLEBBE.NODE((int)b, 0);
                    double E = NLEBBE.MAT((int)NLEBBE.ELEM(i, 0) - 1, 0);
                    double nu = NLEBBE.MAT((int)NLEBBE.ELEM(i, 0) - 1, 1);

                    double A = Area(&NLEBBE.CS);

                    double I = MomentOfInertia(&NLEBBE.CS);

                    //Local stiffness matrix
                    k = StiffnessMatrix_NLEBBE(xa, xb, E, nu, A, I, U, a, b);

                    //Tangent Stiffness matrix
                    t = TangentStiffnessMatrix_NLEBBE(k, xa, xb, E, nu, A, I, U, a, b);

                    //Local Force Vector
                    f = LocalFoceVec_NLEBBE(xa, xb, a, b, NLEBBE.vf, NLEBBE.af,U, E, nu);

                    RearrangeElementStiffness_NLEBBE(k, t, f);

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
                ApplyConstraints_NLEBBE(T, U, NLEBBE.CNODE, NLEBBE.NNODE, R);

                //---------------------------------------------------------------//
                //---------------------------------------------------------------//
                //-----------------------CONTACT--------------------------------//
                //--------------------------------------------------------------//
                //--------------------------------------------------------------//     
                // CONTACT SEARCH
                // GLOBAL SEARCH
                // 1. Create 2 seperate element arrays for 2 beams
                /*int k = 0, l = 0;
                Eigen::MatrixXd Beam1 = Eigen::MatrixXd::Zero(NLEBBE.NELEM, 2);
                Eigen::MatrixXd Beam2 = Eigen::MatrixXd::Zero(NLEBBE.NELEM, 2);
                for (int i = 0; i < NLEBBE.NELEM; i++)
                {
                    if (NLEBBE.ELEM(i, 0) == 1)
                    {
                        Beam1(k, 0) = NLEBBE.ELEM(i, 1);
                        Beam1(k, 1) = NLEBBE.ELEM(i, 2);
                        k++;
                    }
                    else
                    {
                        Beam2(l, 0) = NLEBBE.ELEM(i, 1);
                        Beam2(l, 1) = NLEBBE.ELEM(i, 2);
                        l++;
                    }
                }
                Beam1.conservativeResize(k, 2);
                Beam2.conservativeResize(l, 2);
                std::cout << "Beam 1 connectivity" << std::endl;
                std::cout << Beam1 << std::endl;
                std::cout << "Beam 2 connectivity" << std::endl;
                std::cout << Beam2 << std::endl;
                // 2. Calculate midpoint of each element
                Eigen::MatrixXd Mid1 = Eigen::MatrixXd::Zero(k, 2);
                Eigen::MatrixXd Mid2 = Eigen::MatrixXd::Zero(l, 2);
                std::cout << Beam1.rows();
                for (int i = 0; i < Beam1.rows(); i++)
                {
                    int n1 = Beam1(i, 0) - 1;
                    int n2 = Beam1(i, 1) - 1;
                    double x1 = X_curr(2 * n1);
                    double y1 = X_curr(2 * n1 + 1);
                    double x2 = X_curr(2 * n2);
                    double y2 = X_curr(2 * n2 + 1);
                    Mid1(i, 0) = (x1 + x2) / 2.0;
                    Mid1(i, 1) = (y1 + y2) / 2.0;
                }

                for (int i = 0; i < Beam2.rows(); i++)
                {
                    int n1 = Beam2(i, 0) - 1;
                    int n2 = Beam2(i, 1) - 1;
                    double x1 = X_curr(2 * n1);
                    double y1 = X_curr(2 * n1 + 1);
                    double x2 = X_curr(2 * n2);
                    double y2 = X_curr(2 * n2 + 1);
                    Mid2(i, 0) = (x1 + x2) / 2.0;
                    Mid2(i, 1) = (y1 + y2) / 2.0;
                }
                std::cout << "Beam 1 mid point node data" << std::endl;
                std::cout << Mid1 << std::endl;
                std::cout << "Beam 2 mid point node data" << std::endl;
                std::cout << Mid2 << std::endl;
                // 3. Calculate the minimum distance among the midpoints
                double min = std::numeric_limits<double>::infinity();
                int ele1 = -1, ele2 = -1;
                for (int i = 0; i < Beam1.rows(); i++)
                {
                    double x1 = Mid1(i, 0);
                    double y1 = Mid1(i, 1);

                    for (int j = 0; j < Beam2.rows(); j++)
                    {
                        double x2 = Mid2(j, 0);
                        double y2 = Mid2(j, 1);
                        double dist = sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2));
                        std::cout << dist << std::endl;
                        if (dist < min)
                        {
                            min = dist;
                            ele1 = i;
                            ele2 = j;
                        }
                    }
                }
                std::cout << ele1 << std::endl;
                std::cout << ele2 << std::endl;

                //LOCAL SEARCH
                Eigen::VectorXd x1(2);
                int n1 = Beam1(ele1, 0) - 1;
                int n2 = Beam1(ele1, 1) - 1;
                x1(0) = NLEBBE.NODE(n1, 0);
                x1(1) = NLEBBE.NODE(n1, 1);
                Eigen::VectorXd x2(2);
                x2(0) = NLEBBE.NODE(n2, 0);
                x2(1) = NLEBBE.NODE(n2, 1);

                n1 = Beam2(ele2, 0) - 1;
                n2 = Beam2(ele2, 1) - 1;
                Eigen::VectorXd y1(2);
                y1(0) = NLEBBE.NODE(n1, 0);
                y1(1) = NLEBBE.NODE(n1, 1);
                Eigen::VectorXd y2(2);
                y2(0) = NLEBBE.NODE(n2, 0);
                y2(1) = NLEBBE.NODE(n2, 1);
                std::cout << "Point 1 - Beam 1" << std::endl;
                std::cout << x1 << std::endl;
                std::cout << "Point 2 - Beam 1" << std::endl;
                std::cout << x2 << std::endl;
                std::cout << "Point 1 - Beam 2" << std::endl;
                std::cout << y1 << std::endl;
                std::cout << "Point 2 - Beam 2" << std::endl;
                std::cout << y2 << std::endl;

                Eigen::VectorXd exi_vec(2);
                exi_vec = ContactPoints(x1, x2, y1, y2, NLEBBE.NBODIES);
                double exi = exi_vec(0);
                double exib = exi_vec(1);
                if ((exi > -1) && (exi < 1) && (exib > -1) && (exib < -1))
                {
                    std::cout << exi << std::endl;
                    std::cout << exib << std::endl;
                }
                //Check if contact constraint is satisfied
                double gN, rm, rs;
                rm = NLEBBE.CS.Cir.radius;
                rs = NLEBBE.CS.Cir.radius;
                //Evaluate minimum distance d
                double d = MinimumDistance(x1, x2, y1, y2, exi_vec);
                gN = d - rm - rs;

                //Evaluate Tangent Stiffness matrix
                //Evaluate HTilde Matrix
                Eigen::MatrixXd HTilde(2, 12);
                HTilde = EvaluateHTildeMatrix(x1, x2, y1, y2, exi_vec);
                //Evaluate A matrix
                Eigen::MatrixXd A(2, 2);
                A = EvaluateAMatrix(x1, x2, y1, y2, exi_vec);
                //Evaluate B matrix
                Eigen::MatrixXd B(2, 4);
                B = EvaluateBMatrix(x1, x2, y1, y2, exi_vec);
                //Evaluate C matrix
                Eigen::MatrixXd C(2, 4);
                C = EvaluateCMatrix(x1, x2, y1, y2, exi_vec);
                //Evaluate D matrix
                Eigen::MatrixXd D(2, 12);
                D = EvaluateDMatrix(A, B, C, x1, x2, y1, y2, exi_vec);
                Eigen::VectorXd dvec(12), dbvec(12);
                for (int i = 0; i < 12; i++)
                {
                    dvec(i) = D(0, i);
                    dbvec(i) = D(1, i);
                }
                //Evaluate normal vector
                Eigen::VectorXd n(2);
                n = EvaluateNormalVector(x1, x2, y1, y2, exi_vec);
                //Evaluate E matrix
                Eigen::MatrixXd E(12, 12);
                E = EvaluateEMatrix(x1, x2, y1, y2, exi_vec, n, dvec, dbvec);
                //Evaluate F matrix
                Eigen::MatrixXd F(12, 12);
                F = EvaluateFMatrix(x1, x2, y1, y2, exi_vec, n, dvec, dbvec);
                //Evaluate G matrix
                Eigen::MatrixXd G(12, 12);
                G = EvaluateGMatrix(x1, x2, y1, y2, exi_vec, n, dvec, dbvec, HTilde);

                //Evaluate Contact Local Tangent stiffness matrix
                Eigen::MatrixXd Kcon(12, 12);
                Kcon = HTilde.transpose() * n * n.transpose() * HTilde + gN * (E + E.transpose() + F + G);

                //Evaluate Contact Element Residual Vector 
                Eigen::VectorXd Rcon(12);
                Rcon = gN * HTilde.transpose() * n;*/

                //Assembly


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
                for(int i=0;i<NLEBBE.NNODE*NLEBBE.NDOF;i++)
                    U_new(i) = U(i) + dU(i);

                //Error calculation
                for (int i = 0; i < NLEBBE.NNODE * NLEBBE.NDOF; i++)
                    error(i) = abs(U_new(i)-U(i));
                max = error(0);
                for (int i = 0; i < NLEBBE.NNODE * NLEBBE.NDOF; i++)
                    if (max < error(i))
                        max = error(i);

                //Assignment for next iteration
                for (int i = 0; i < NLEBBE.NNODE * NLEBBE.NDOF; i++)
                    U(i) = U_new(i);
                iter++;

                //Update current configuration
                //for (int i = 0; i < (int)NLEBBE.NNODE; i++)
                //{
                //   X_curr(NLEBBE.NDM * i) = X_curr(NLEBBE.NDM * i) + U(NLEBBE.NDOF * i);
                //    X_curr(NLEBBE.NDM * i + 1) = X_ref(NLEBBE.NDM * i + 1) + U(NLEBBE.NDOF * i + 1);
                //}
                
            } while (max > pow(10, -3) && iter < maxiter);
            if (fiter < NLEBBE.NLS)
            {
                PostProcessing(U, NLEBBE, fiter);
                fiter++;
            }
            if (iter == maxiter)
            {
                NLEBBE.vf = NLEBBE.vf - 0.5;
                std::cout << "Solution not converging for the given load" << std::endl;
                std::cout << "Trying with a load with less increment" << std::endl;
            }
            std::cout << "Load  " << NLEBBE.vf << std::endl;
            std::cout << "Displacements  " << std::endl << U << std::endl;

            NLEBBE.vf += 1;
            //Load step increment
        } while (fiter < NLEBBE.NLS);

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
                                Seq_1 = Equivalent_StiffnessMatrix_FirstOrder(Strain, VAMBE.CMP.inittwist, VAMBE.CS.Rect.width, file1);

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

                            Seq_1 = Equivalent_StiffnessMatrix_FirstOrder(Strain1, VAMBE.CMP.inittwist, VAMBE.CS.Rect.width, file1);
                            Seq_2 = Equivalent_StiffnessMatrix_FirstOrder(Strain2, VAMBE.CMP.inittwist, VAMBE.CS.Rect.width, file1);

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

                            Seq_1 = Equivalent_StiffnessMatrix_FirstOrder(Strain, VAMBE.CMP.inittwist, VAMBE.CS.Rect.width, file1);

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
    //--------------------------------------------------------//
    //----------3D EulerBernouli Beam Element-----------------//
    //-----------------Large displacement---------------------//
    //-----------------Updated Lagrangian---------------------//
    //---------------------Quadratic--------------------------//
    else if (choice == 5)
    {
        NonLinearEulerBernouliBeamElement3D EBBE3D = ReadEBBE3DElement();

        Eigen::VectorXd dGU = Eigen::VectorXd::Zero(EBBE3D.NDOF * EBBE3D.NNODE);
        Eigen::VectorXd GU = Eigen::VectorXd::Zero(EBBE3D.NDOF * EBBE3D.NNODE);
        Eigen::VectorXd error = Eigen::VectorXd::Zero(EBBE3D.NDOF * EBBE3D.NNODE);
        Eigen::VectorXd GU_new = Eigen::VectorXd::Zero(EBBE3D.NDOF * EBBE3D.NNODE);
        Eigen::VectorXd GR = Eigen::VectorXd::Zero(EBBE3D.NDOF * EBBE3D.NNODE);
        Eigen::SparseMatrix<double, Eigen::ColMajor> GT(EBBE3D.NDOF* EBBE3D.NNODE, EBBE3D.NDOF* EBBE3D.NNODE);

        double max;
        int iter = 1;
        int fiter = 1;
        int maxiter = 50;

        std::fstream file1, file2;
        file1.open("E:/Adhithya/MTech_Res_Thesis/Cpp/ThesisCode/TextileComposites/TextileComposites/Result_Log.txt", std::fstream::in | std::fstream::out);
        file2.open("E:/Adhithya/MTech_Res_Thesis/Cpp/ThesisCode/TextileComposites/TextileComposites/Results.txt", std::fstream::in | std::fstream::out);
    
        //declare and initialize material parameters and other input variables
        double *D = new double[7];
        D[0] = EBBE3D.E;
        D[1] = EBBE3D.nu;
        D[2] = EBBE3D.Bp;
        D[3] = EBBE3D.Hp;
        D[4] = EBBE3D.Zx;
        D[5] = EBBE3D.Zy;
        D[6] = EBBE3D.Zz;
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
                    T = new double* [EBBE3D.NDOF*EBBE3D.NEN];
                    for (int i = 0; i < EBBE3D.NDOF * EBBE3D.NEN; i++)
                        T[i] = new double[EBBE3D.NDOF * EBBE3D.NEN];

                    double* R = new double[EBBE3D.NDOF * EBBE3D.NEN];

                    for (int i = 0; i < EBBE3D.NELEM; i++)
                    {
                        int a = (int)EBBE3D.ELEM(i, 1) - 1;
                        int b = (int)EBBE3D.ELEM(i, 2) - 1;
                        int c = (int)EBBE3D.ELEM(i, 3) - 1;

                        //declare and initialize position vectors
                        /*double* X[2];
                        for (int j = 0; j < 2; j++)
                            X[j] = new double[3];*/
                        double X[3][3];

                        for (int j = 0; j < 3; j++)
                        {
                            X[0][j] = EBBE3D.NODE(a, j);
                            X[1][j] = EBBE3D.NODE(b, j);
                            X[2][j] = EBBE3D.NODE(c, j);
                        }
                        double h = X[1][0] - X[0][0];

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

                        for (int j = 0; j < EBBE3D.NDOF * EBBE3D.NEN; j++)
                        {
                            R[j] = 0;
                            for (int k = 0; k < EBBE3D.NDOF * EBBE3D.NEN; k++)
                                T[j][k] = 0;
                        }
                        RKt(D, X, U, T, R);

                        /*for (int i = 0; i < 12; i++)
                            std::cout << R[i] << std::endl;

                        for (int j = 0; j < 12; j++)
                        {
                            for (int k = 0; k < 12; k++)
                                std::cout << T[j][k] << "   ";
                            std::cout << std::endl;
                        }*/

                        //Assembly
                        for (int j = 0; j < EBBE3D.NDOF * EBBE3D.NEN; j++)
                        {
                            GR(6 * a + j) += R[j];
                            for (int k = 0; k < EBBE3D.NDOF * EBBE3D.NEN; k++)
                                GT.coeffRef(6 * a + j, 6 * a + k) += T[j][k];
                        }
                    }

                    //Assign Boundary conditions
                    for (int j = 0; j < 6; j++)
                    {
                        for (int k = 0; k < EBBE3D.NNODE * EBBE3D.NDOF; k++)
                            GT.coeffRef(j, k) = 0;
                        GR(j) = 0;
                        GT.coeffRef(j, j) = 1;
                    }
                    GR(EBBE3D.NNODE*EBBE3D.NDOF - 4) -= load;
                    
                    //Solve the equation
                    GT.makeCompressed();

                    Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>> solver;
                    solver.analyzePattern(GT);
                    solver.factorize(GT); //LU decomposition

                    //std::cout<<J<<std::endl;

                    assert(solver.info() == Eigen::Success);

                    dGU = solver.solve(-GR);

                    //Calculate errors and update for next iteration
                    for (int j = 0; j < EBBE3D.NNODE * EBBE3D.NDOF; j++)
                        GU_new(j) = GU(j) + dGU(j);

                    //Error calculation
                    for (int j = 0; j < EBBE3D.NNODE * EBBE3D.NDOF; j++)
                        error(j) = abs(GU_new(j) - GU(j));

                    max = error(0);
                    for (int j = 0; j < EBBE3D.NNODE * EBBE3D.NDOF; j++)
                        if (max < error(j))
                            max = error(j);

                    //Assignment for next iteration
                    for (int j = 0; j < EBBE3D.NNODE * EBBE3D.NDOF; j++)
                        GU(j) = 0.5 * GU_new(j) + 0.5 * GU(j);
                    iter++;

                    //Free memory
                    for (int i = 0; i < EBBE3D.NDOF * EBBE3D.NEN; i++)
                        delete T[i];
                    delete[] R;

                    std::cout << max << std::endl;
                } while (max > pow(10, -6) && iter < maxiter);
                if (iter < maxiter)
                {
                    fiter++;
                    load += 4.167;
                    /*for (int j = 0; j < EBBE3D.NNODE * EBBE3D.NDOF; j++)
                    {
                        //file2 << "Displacements of node " << j + 1 << std::endl;
                        //for (int k = 0; k < 3; k++)
                            //file2 << U(12 * (j - 1) + 6 + k) << std::endl;
                        file2 << GU(j) << std::endl;
                    }*/
                    file2 << GU(EBBE3D.NNODE * EBBE3D.NDOF - 4) << std::endl;
                }
                else
                {
                    std::cout << "Code didn't converge" << std::endl;
                    break;
                }
            } while (fiter < EBBE3D.NLS);
        }
    }
    //--------------------------------------------------------//
    //----------3D EulerBernouli Beam Element-----------------//
    //-----------------Large displacement---------------------//
    //-----------------Updated Lagrangian---------------------//
    //--------------------Quadratic---------------------------//
    //---------------Node to Node Contact---------------------//
    else if (choice == 6)
    {
        //Master 1
        //Slave 2
        double ms;
        //Read master data
        ms = 1;
        NonLinearEulerBernouliBeamElement3D EBBE3D1 = ReadEBBE3DElement(ms);
        //Read slave data
        ms = 2;
        NonLinearEulerBernouliBeamElement3D EBBE3D2 = ReadEBBE3DElement(ms);

        Eigen::VectorXd dGU = Eigen::VectorXd::Zero(EBBE3D1.NDOF * (EBBE3D1.NNODE + EBBE3D2.NNODE));
        Eigen::VectorXd GU = Eigen::VectorXd::Zero(EBBE3D1.NDOF * (EBBE3D1.NNODE + EBBE3D2.NNODE));
        Eigen::VectorXd error = Eigen::VectorXd::Zero(EBBE3D1.NDOF * (EBBE3D1.NNODE + EBBE3D2.NNODE));
        Eigen::VectorXd GU_new = Eigen::VectorXd::Zero(EBBE3D1.NDOF * (EBBE3D1.NNODE + EBBE3D2.NNODE));
        Eigen::VectorXd GR = Eigen::VectorXd::Zero(EBBE3D1.NDOF * (EBBE3D1.NNODE + EBBE3D2.NNODE));
        Eigen::SparseMatrix<double, Eigen::ColMajor> GT(EBBE3D1.NDOF * (EBBE3D1.NNODE + EBBE3D2.NNODE), EBBE3D1.NDOF * (EBBE3D1.NNODE + EBBE3D2.NNODE));

        double max;
        int iter = 1;
        int fiter = 1;
        int maxiter = 20;

        std::fstream file1, file2;
        file1.open("E:/Adhithya/MTech_Res_Thesis/Cpp/ThesisCode/TextileComposites/TextileComposites/Result_Log.txt", std::fstream::in | std::fstream::out);
        file2.open("E:/Adhithya/MTech_Res_Thesis/Cpp/ThesisCode/TextileComposites/TextileComposites/Results.txt", std::fstream::in | std::fstream::out);

        //declare and initialize material parameters and other input variables
        //master
        double* D1 = new double[7];
        D1[0] = EBBE3D1.E;
        D1[1] = EBBE3D1.nu;
        D1[2] = EBBE3D1.Bp;
        D1[3] = EBBE3D1.Hp;
        D1[4] = EBBE3D1.Zx;
        D1[5] = EBBE3D1.Zy;
        D1[6] = EBBE3D1.Zz;

        //slave
        double* D2 = new double[7];
        D2[0] = EBBE3D2.E;
        D2[1] = EBBE3D2.nu;
        D2[2] = EBBE3D2.Bp;
        D2[3] = EBBE3D2.Hp;
        D2[4] = EBBE3D2.Zx;
        D2[5] = EBBE3D2.Zy;
        D2[6] = EBBE3D2.Zz;

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
                    T = new double* [12];
                    for (int i = 0; i < 12; i++)
                        T[i] = new double[12];
                    double* R = new double[12];

                    for (int i = 0; i < (EBBE3D1.NELEM + EBBE3D2.NELEM); i++)
                    {
                        //master
                        if (i < EBBE3D1.NELEM)
                        {
                            int a = (int)EBBE3D1.ELEM(i, 1) - 1;
                            int b = (int)EBBE3D1.ELEM(i, 2) - 1;

                            //declare and initialize position vectors
                            /*double* X[2];
                            for (int j = 0; j < 2; j++)
                                X[j] = new double[3];*/
                            double X[2][3];

                            for (int j = 0; j < 3; j++)
                            {
                                X[0][j] = EBBE3D1.NODE(a, j);
                                X[1][j] = EBBE3D1.NODE(b, j);
                            }
                            double h = X[1][0] - X[0][0];

                            //declare and initialize dofs
                            /*double* U[2];
                            for (int j = 0; j < 2; j++)
                                U[j] = new double[6];*/
                            double U[2][6];

                            for (int j = 0; j < 6; j++)
                            {
                                U[0][j] = GU(6 * a + j);
                                U[1][j] = GU(6 * b + j);
                            }

                            for (int j = 0; j < 12; j++)
                            {
                                R[j] = 0;
                                for (int k = 0; k < 12; k++)
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
                            }*/

                            //Apply contact constraint

                            //Assembly
                            for (int j = 0; j < 12; j++)
                            {
                                GR(6 * a + j) += R[j];
                                for (int k = 0; k < 12; k++)
                                    GT.coeffRef(6 * a + j, 6 * a + k) += T[j][k];
                            }
                        }
                        //slave
                        else
                        {
                            int a = (int)EBBE3D2.ELEM(i, 1) - 1;
                            int b = (int)EBBE3D2.ELEM(i, 2) - 1;

                            //declare and initialize position vectors
                            /*double* X[2];
                            for (int j = 0; j < 2; j++)
                                X[j] = new double[3];*/
                            double X2[2][3], X1[2][3];

                            for (int j = 0; j < 3; j++)
                            {
                                X1[0][j] = EBBE3D1.NODE(a, j);
                                X1[1][j] = EBBE3D1.NODE(b, j);
                                X2[0][j] = EBBE3D2.NODE(a, j);
                                X2[1][j] = EBBE3D2.NODE(b, j);
                            }
                            double h = X2[1][0] - X2[0][0];

                            //declare and initialize dofs
                            /*double* U[2];
                            for (int j = 0; j < 2; j++)
                                U[j] = new double[6];*/
                            double U2[2][6], U1[2][6];

                            for (int j = 0; j < 6; j++)
                            {
                                U2[0][j] = GU(6 * a + j);
                                U2[1][j] = GU(6 * b + j);
                                U1[0][j] = GU(6 * a + j);
                                U1[1][j] = GU(6 * b + j);
                            }

                            for (int j = 0; j < 12; j++)
                            {
                                R[j] = 0;
                                for (int k = 0; k < 12; k++)
                                    T[j][k] = 0;
                            }

                            //Local Residual and tangent matrix
                            //RKt(D1, X, U, T, R);

                            /*for (int i = 0; i < 12; i++)
                                std::cout << R[i] << std::endl;

                            for (int j = 0; j < 12; j++)
                            {
                                for (int k = 0; k < 12; k++)
                                    std::cout << T[j][k] << "   ";
                                std::cout << std::endl;
                            }*/

                            //Apply contact constraint


                            //Assembly
                            for (int j = 0; j < 12; j++)
                            {
                                GR(6 * a + j) += R[j];
                                for (int k = 0; k < 12; k++)
                                    GT.coeffRef(6 * a + j, 6 * a + k) += T[j][k];
                            }

                        }
                    }

                    //Assign Boundary conditions
                    for (int j = 0; j < 6; j++)
                    {
                        for (int k = 0; k < EBBE3D1.NNODE * EBBE3D1.NDOF; k++)
                            GT.coeffRef(j, k) = 0;
                        GR(j) = 0;
                        GT.coeffRef(j, j) = 1;
                    }
                    GR(EBBE3D1.NNODE * EBBE3D1.NDOF - 4) -= load;

                    //Solve the equation
                    GT.makeCompressed();

                    Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>> solver;
                    solver.analyzePattern(GT);
                    solver.factorize(GT); //LU decomposition

                    //std::cout<<J<<std::endl;

                    assert(solver.info() == Eigen::Success);

                    dGU = solver.solve(-GR);

                    //Calculate errors and update for next iteration
                    for (int j = 0; j < EBBE3D1.NNODE * EBBE3D1.NDOF; j++)
                        GU_new(j) = GU(j) + dGU(j);

                    //Error calculation
                    for (int j = 0; j < EBBE3D1.NNODE * EBBE3D1.NDOF; j++)
                        error(j) = abs(GU_new(j) - GU(j));

                    max = error(0);
                    for (int j = 0; j < EBBE3D1.NNODE * EBBE3D1.NDOF; j++)
                        if (max < error(j))
                            max = error(j);

                    //Assignment for next iteration
                    for (int j = 0; j < EBBE3D1.NNODE * EBBE3D1.NDOF; j++)
                        GU(j) = 0.5 * GU_new(j) + 0.5 * GU(j);
                    iter++;

                    //Free memory
                    for (int i = 0; i < 12; i++)
                        delete T[i];
                    delete[] R;

                    std::cout << max << std::endl;
                } while (max > pow(10, -6) && iter < maxiter);
                if (iter < maxiter)
                {
                    fiter++;
                    load += 4.167;
                    /*for (int j = 0; j < EBBE3D.NNODE * EBBE3D.NDOF; j++)
                    {
                        //file2 << "Displacements of node " << j + 1 << std::endl;
                        //for (int k = 0; k < 3; k++)
                            //file2 << U(12 * (j - 1) + 6 + k) << std::endl;
                        file2 << GU(j) << std::endl;
                    }*/
                    file2 << GU(EBBE3D1.NNODE * EBBE3D1.NDOF - 4) << std::endl;
                }
                else
                {
                    std::cout << "Code didn't converge" << std::endl;
                    break;
                }
            } while (fiter < EBBE3D1.NLS);
        }
        }
    else
        std::cout << "Wrong option" << std::endl;
}

#endif




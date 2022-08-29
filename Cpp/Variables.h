#pragma once
//Seperate structures for each of the elements in the FEM code solves the problem without much hassle.
//Even in the readinp file we can read accordingly
#include<Eigen/Sparse>
#include<Eigen/Dense>
#include<Eigen/SparseCholesky>
#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<math.h>
#include<cassert>


struct Stiffness
{
    Eigen::MatrixXd Sl;
    Eigen::MatrixXd Sln;
    Eigen::MatrixXd Sn;
};

struct Bar
{
    double A11;

    double B11;
    double B12;
    double B16;

    double D11;
    double D12;
    double D22;
    double D16;
    double D26;
    double D66;
};

struct DoubleBar
{
    double A11;

    double B11;
    double B16;

    double D11;
    double D16;
    double D66;
};

struct Variables
{
    Bar B;
    DoubleBar BB;
};
/*
struct ELSET
{
    char* elset;
    Eigen::VectorXd elements;
};*/

struct Orthotropic
{
    double E11;
    double E22;
    double E33;
    double nu12;
    double nu13;
    double nu23;
    double G12;
    double G13;
    double G23;
};

struct Material
{
    double E;
    double nu;
};

struct CrossSection
{
    double height;
    double width;
};

struct CompositeMaterialProperties
{
    Eigen::VectorXd Orient;
    Eigen::VectorXd inittwist;
    double np;
    double E11;
    double E22;
    double E33;
    double G12;
    double G13;
    double nu12;
    double nu13;
    double nu23;
    double G23;

};

struct Loading
{
    double* cload;
};

struct Boundary
{
    Eigen::VectorXd F1, FN, M1, MN;
    Eigen::VectorXd u1, uN, theta1, thetaN;
};

struct LinearBarElement
{
    int NMAT;

    int NBODIES;

    int NDOF;

    Eigen::MatrixXd NODE;

    Eigen::MatrixXd ELEM;

    int NNODE;

    int NELEM;

    Eigen::MatrixXd CNODE;

    Eigen::MatrixXd LOAD;

    Eigen::MatrixXd MAT;

    //Material* MAT;

    CrossSection CS;
};

struct NonLinearBarElement
{
    int NMAT;

    int NBODIES;

    int NDOF;

    int NLS;

    Eigen::MatrixXd NODE;

    Eigen::MatrixXd ELEM;

    int NNODE;

    int NELEM;

    Eigen::MatrixXd CNODE;

    Eigen::MatrixXd MAT;

    CrossSection CS;
};

struct NonLinearEulerBernouliBeamElement
{
    int NMAT;

    int NBODIES;

    int NDOF;

    int NLS;

    Eigen::MatrixXd NODE;

    Eigen::MatrixXd ELEM;

    int NNODE;

    int NELEM;

    Eigen::MatrixXd CNODE;

    Eigen::MatrixXd MAT;

    CrossSection CS;

    double vf;
    double af;
};

struct VAMBeamElement
{
    //Number of Materials
    int NMAT;
    //No. of load steps
    int NLS;
    //Number of bodies
    int NBODIES;
    //Number of degrees of freedom
    int NDOF;
    //Number of cross sections
    int NCS;
    //Nodal Coordinates
    Eigen::MatrixXd NODE;
    //Number of Nodes
    int NNODE;
    //Number of Elements
    int NELEM;
    //Material constants
    CompositeMaterialProperties CMP;
    //Geometry
    CrossSection CS;
    Eigen::MatrixXd ELEM;
    //Loading
    Boundary B;
    //struct LOAD load;
    //struct ELSET elset;
};

//Functions in ReadInpFile.cpp
LinearBarElement ReadLBEFile();
NonLinearBarElement ReadNLBEFile();
NonLinearEulerBernouliBeamElement ReadNLEBBEFile();
VAMBeamElement ReadVAMBEFile();

//Functions in BarElement_Linear_1D.cpp
void StiffnessMatrix_LBE(double A, double E, double xa, double xb, Eigen::MatrixXd& k);
void LocalForceVec_LBE(double f1, double f2, Eigen::VectorXd& f);
void ApplyConstraints_LBE(Eigen::MatrixXd CNODE, int NNODE, Eigen::VectorXd F,
    Eigen::SparseMatrix<double, Eigen::ColMajor>& K);

//Functions in BarElement_NonLinear_1D.cpp
void TangentStiffnessMatrix_NLBE(Eigen::MatrixXd& t, Eigen::VectorXd U, int StartNode, int EndNode, double h);
void StiffnessMatrix_NLBE(Eigen::MatrixXd& k, Eigen::VectorXd U, int StartNode, int EndNode, double h);
void LocalForceVec_NLBE(Eigen::VectorXd& f, int StartNode, int EndNode, double h);
void ApplyConstraints_NLBE(Eigen::SparseMatrix<double, Eigen::RowMajor>& T, Eigen::VectorXd& U, Eigen::VectorXd& R, int NNODE, Eigen::MatrixXd CNODE);

//Functions in EulerBernoulli2D.cpp
void StiffnessMatrix_NLEBBE(Eigen::MatrixXd& k, double xa, double xb, double E, double nu, double base, double height,
    Eigen::VectorXd& U, int a, int b);
void LocalFoceVec_NLEBBE(Eigen::VectorXd& f, double xa, double xb, int a, int b, double vf, double af);
void TangentStiffnessMatrix_NLEBBE(Eigen::MatrixXd& t, Eigen::MatrixXd& k, double xa, double xb, double E, double nu,
    double base, double height, Eigen::VectorXd& U, int a, int b);
void ApplyConstraints_NLEBBE(Eigen::SparseMatrix<double, Eigen::ColMajor>& T, Eigen::VectorXd& U, Eigen::MatrixXd& CNODE,
    int NNODE, Eigen::VectorXd& R);
void RearrangeElementStiffness_NLEBBE(Eigen::MatrixXd& k, Eigen::MatrixXd& t, Eigen::VectorXd& f);
//GaussPoints GaussQuadraturePoints();

//Functions for VAM Beam Element
Stiffness StiffnessMatrix_VAM(double E1, double E2, double nu12, double  G12, double nu21, double width, double height,
    double k1, Eigen::VectorXd Orient, int np, std::fstream& file1);
Eigen::MatrixXd Equivalent_ClassicalStiffnessModel_VAM(Stiffness S, Eigen::VectorXd Strain, std::fstream& file1);

//Functions from GEBT.cpp
Eigen::VectorXd Element_Residual(Eigen::VectorXd U1, Eigen::VectorXd U2, Eigen::VectorXd F1, Eigen::VectorXd F2, double h,
    Eigen::MatrixXd S1, Eigen::MatrixXd S2, double Theta);
Eigen::VectorXd Element_Residual(Eigen::VectorXd U1, Eigen::VectorXd F1, double h,
    Eigen::MatrixXd S1, double Theta, int nodenum, Boundary B);

Eigen::MatrixXd Element_Jacobian(Eigen::VectorXd U1, Eigen::VectorXd U2, Eigen::VectorXd F1, Eigen::VectorXd F2, double h,
    Eigen::MatrixXd S1, Eigen::MatrixXd S2, double Theta, std::fstream& file1);
Eigen::MatrixXd Element_Jacobian(Eigen::VectorXd U1, Eigen::VectorXd F1, double h, Eigen::MatrixXd S1, double Theta, int nodenum,
    std::fstream& file1);
void Update_Strains(Eigen::VectorXd& Tau, int nnode, Stiffness S, Eigen::VectorXd U, std::fstream &file1);
// VARIABLES_H

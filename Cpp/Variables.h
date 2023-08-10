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
#include<memory>
#include "sms.h"

#define M_PI 3.14159265358979323846

struct Boundary
{
    Eigen::VectorXd a1, b1;
    int x1, y1;
    Eigen::VectorXd aN, bN;
    int xN, yN;
};

struct Rectangle
{
    double height;
    double width;
};

struct Circle
{
    double radius;
};

struct CrossSection
{
    struct Rectangle Rect;
    struct Circle Cir;
    std::string choice;
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

class NLEBBE2D
{
private:
    int NMAT;
    int NDOF;
    int NLS;
    int NNODE;
    int NELEM;
    int NEN;
    double E, nu, b, h, vf, af;
    Eigen::MatrixXd NODE;
    Eigen::MatrixXd ELEM;
    Eigen::MatrixXd CNODE;
public:
    NLEBBE2D();
    int get_nnode();
    int get_nelem();
    int get_ndof();
    int get_nls();
    int get_nen();
    int get_connectivity(int i, int j);
    double get_coordinates(int i, int j);
    double get_modelprop(std::string str);
    double get_loadprop(std::string str);
    void set_loadprop(std::string str, double val);
    int get_cnode(int i, int j);
    Eigen::MatrixXd StiffnessMatrix_NLEBBE(double xa, double xb, double E, double nu, double base, double height, Eigen::VectorXd& U, int a, int b);
    Eigen::VectorXd LocalFoceVec_NLEBBE(double xa, double xb, int a, int b, double vf, double af, Eigen::VectorXd& U, double E, double nu);
    Eigen::MatrixXd TangentStiffnessMatrix_NLEBBE(Eigen::MatrixXd& k, double xa, double xb, double E, double nu, double base, double height, Eigen::VectorXd& U, int a, int b);
    void ApplyConstraints_NLEBBE(Eigen::SparseMatrix<double, Eigen::ColMajor>& T, Eigen::VectorXd& U, int NNODE, Eigen::VectorXd& R);
    void RearrangeElementStiffness_NLEBBE(Eigen::MatrixXd& k, Eigen::MatrixXd& t, Eigen::VectorXd& f);
};

class NLEBBE3D
{
private:
    int NEN;
    int NDOF;
    int NNODE;
    int NELEM;
    int NLS;
    int NBEAMS;
    Eigen::MatrixXd NODE;
    Eigen::MatrixXd ELEM;
    double E, nu, Bp, Hp, Zx, Zy, Zz;
    int loadnode;
    double DIA;

public:
    NLEBBE3D(int nbeams);
    //getters
    int get_nen();
    int get_ndof();
    int get_nnode();
    int get_nelem();
    int get_nls();
    int get_coordinates(int i, int j);
    int get_connectivity(int i, int j);
    int get_nbeams();
    double get_matprop(std::string str);
    int get_loadnode();
    void RKt(double D[7], double X[3][3], double U[3][6], double** T, double* R);
    void RKtLin(double D[7], double X[2][3], double U[2][6], double** T, double* R);
    double get_diameter();
};

class BeamContact
{
private:
    double epsilon;
    double NEN;
    Eigen::MatrixXd ELEM;
    int NBEAMS;
    Eigen::MatrixXd GlobalELEM;
public:
    double get_penaltyparameter();
    int get_nen();
    Eigen::MatrixXd get_Globalelem();
    void LocalContactSearch(NLEBBE3D* EBBE3D);
    int get_slavesegment();
    BeamContact(const int nbeams, std::string str);
    void Contact_NTN(double D[4], double X1[3], double X2[3], double u1[6], double u2[6], double* CR, double** CT, double* g);
    void GlobalContactSearch(NLEBBE3D* EBBE3D);
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
    //Geometry
    CrossSection CS;
    Eigen::MatrixXd ELEM;
    //Material constants
    Eigen::VectorXd inittwist;
    //Boundary
    Boundary B;
};

class VTKGrid
{
private:
    Eigen::MatrixXd Points;
    Eigen::MatrixXd Cells;
    int nelem, nnode, ntype, eletype;
    std::string dim;

public:
    VTKGrid(NLEBBE3D EBBE3D, int refine, double load, std::string BeamName, Eigen::MatrixXd GU, int ndim, int eletype, int ndof, int ncells);
    VTKGrid(NLEBBE2D EBBE2D, double load, std::string BeamName, Eigen::MatrixXd GU, int ndim, int eletype, int ndof, int ncells);
    void WriteVTK(int load, std::string BeamName, std::string filepath);
    Eigen::MatrixXd FindSurfacePoints(Eigen::VectorXd n1, Eigen::VectorXd node, double radius, int refine);
    Eigen::MatrixXd FindCells(int p1, int p2, int refine);
};

//Functions in ReadInpFile.cpp
LinearBarElement ReadLBEFile();
NonLinearBarElement ReadNLBEFile();
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




//Functions for VAM Beam Element
Eigen::MatrixXd Equivalent_StiffnessMatrix_FirstOrder(Eigen::VectorXd Strain, Eigen::VectorXd inittwist, double b, std::fstream& file1);

//Functions from GEBT.cpp
Eigen::VectorXd Element_Residual(Eigen::VectorXd U1, Eigen::VectorXd U2, Eigen::VectorXd F1, Eigen::VectorXd F2, double h,
    Eigen::MatrixXd S1, Eigen::MatrixXd S2, double Theta);
Eigen::VectorXd Element_Residual(Eigen::VectorXd U1, Eigen::VectorXd F1, double h,
    Eigen::MatrixXd S1, double Theta, int nodenum, Boundary B, Eigen::VectorXd U0);

Eigen::MatrixXd Element_Jacobian(Eigen::VectorXd U1, Eigen::VectorXd U2, Eigen::VectorXd F1, Eigen::VectorXd F2, double h,
    Eigen::MatrixXd S1, Eigen::MatrixXd S2, double Theta, std::fstream& file1);
Eigen::MatrixXd Element_Jacobian(Eigen::VectorXd U1, Eigen::VectorXd F1, double h, Eigen::MatrixXd S1, double Theta, int nodenum,
    std::fstream& file1);
Eigen::VectorXd Update_Strains(VAMBeamElement VAMBE, Eigen::VectorXd* U, std::fstream& file1);
// VARIABLES_H

//Functions from 3DBeamElement_Contact_NodeToNode.cpp
void ContactSearch(NLEBBE3D EBBE3D1, NLEBBE3D EBBE3D2, int** ContactPairs, std::string choice);
void Contact_NTN(double D[4], double X1[3], double X2[3], double u1[6], double u2[6], double* CR, double** CT, double* g);
void PostProcessing(Eigen::MatrixXd X, Eigen::VectorXd U, double load, std::string BeamElement, int nnode, int ndm, int ndof);

//Functions from 3DBeamElement_EulerBernouli_Contact_STS
void Contact_STS(double D[6], double X1[3], double X2[3]
    , double X3[3], double X4[3], double u1[6], double u2[6], double u3[6], double u4[6],
    Eigen::SparseMatrix<double, Eigen::ColMajor>* GT, Eigen::VectorXd* GR, int mnode1, int mnode2,
    int snode1, int snode2);
void Contact_STS_Endpoints(double D[5], double X1[3]
    , double X2[3], double u1[6], double u2[6],
    Eigen::SparseMatrix<double, Eigen::ColMajor>* GT, Eigen::VectorXd* GR, int mnode, int snode);


//Functions from 3DBeamElement_NonLinear_EulerBernoulli_QuadraticInterpolation
double GaussIntegrationPoints(int i, int j);
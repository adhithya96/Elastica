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
#include<iomanip>

#define M_PI 3.14159265358979323846

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

class VAMBeamElement
{
private:
    int NLS;
    int NDOF;
    int NNODE;
    int NELEM;
    int NDM;
    double D;
    Eigen::MatrixXd NODE;
    Eigen::MatrixXd ELEM;
    Eigen::VectorXd inittwist;
    double b;
    Eigen::VectorXd a1, b1, aN, bN, axis;
public:
    VAMBeamElement();
    VAMBeamElement(int beamnum);
    int get_nls();
    int get_ndof();
    int get_nnode();
    int get_nelem();
    int get_ndim();
    Eigen::VectorXd get_axis();
    double get_diameter(); 
    Eigen::VectorXd get_inittwist();
    double get_width();
    double get_coordinates(int i, int j);
    int get_connectivity(int i, int j);
    Eigen::VectorXd get_load(int choice);
    void set_load(int choice, int dir, double value);
    //Functions from GEBT.cpp
    Eigen::VectorXd Element_Residual(Eigen::VectorXd U1, Eigen::VectorXd U2, Eigen::VectorXd F1, Eigen::VectorXd F2, double h,
        Eigen::MatrixXd S1, Eigen::MatrixXd S2, double Theta);
    Eigen::VectorXd Element_Residual(Eigen::VectorXd U1, Eigen::VectorXd F1, double h,
        Eigen::MatrixXd S1, double Theta, int nodenum, Eigen::VectorXd a1, Eigen::VectorXd b1, Eigen::VectorXd aN,
        Eigen::VectorXd bN, Eigen::VectorXd U0);
    Eigen::VectorXd Element_Residual(Eigen::VectorXd U1, Eigen::VectorXd F1, double h,
        Eigen::MatrixXd S1, double Theta, int nodenum, Eigen::VectorXd U0);

    Eigen::MatrixXd Element_Jacobian(Eigen::VectorXd U1, Eigen::VectorXd U2, Eigen::VectorXd F1, Eigen::VectorXd F2, double h,
        Eigen::MatrixXd S1, Eigen::MatrixXd S2, double Theta, std::fstream& file1);
    Eigen::MatrixXd Element_Jacobian(Eigen::VectorXd U1, Eigen::VectorXd F1, double h, Eigen::MatrixXd S1, double Theta, int nodenum,
        std::fstream& file1);
    Eigen::MatrixXd Equivalent_StiffnessMatrix_FirstOrder(std::fstream& file1);
    VAMBeamElement(int nls, int ndof, int ndim, int nnode, int nelem, std::vector<std::vector<double>>& node,
        std::vector<std::vector<int>>& conn, double dia, Eigen::VectorXd axis);

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
    int NDM;
    Eigen::MatrixXd NODE;
    Eigen::MatrixXd ELEM;
    double E, nu, Bp, Hp, Zx, Zt, Zz, C10, C01;
    double E1, E2, E3, nu12, nu13, nu23, G12, G13, G23, B, H, Zy;
    int loadnode;
    double DIA;
    std::string MAT;

public:
    NLEBBE3D() {};
    NLEBBE3D(int choice);
    NLEBBE3D(int choice, std::string str);
    NLEBBE3D(std::string str);
    NLEBBE3D(int yn, double ys, double yd, int refine, std::string str, int beamnum, int nls, int nen, int ndm);
    NLEBBE3D(int nele, double bd, double bl, int nls, int nen, int ndm, double E, double nu, int beamnum);
    //getters
    int get_nen();
    int get_ndof();
    int get_nnode();
    int get_nelem();
    int get_ndim();
    int get_nls();
    double get_coordinates(int i, int j);
    int get_connectivity(int i, int j);
    int get_nbeams();
    double get_matprop(std::string str);
    int get_loadnode();
    std::string get_matmodel();
    void UpdateNodes(NLEBBE3D* EBBE3D, Eigen::VectorXd* GU, int m);
    void RKt(double D[7], double X[3][3], double U[3][6], double** T, double* R);
    void RKtLin(double D[7], double X[2][3], double U[2][6], double** T, double* R);
    double get_diameter();
    void RKt_NHC(double D[7], double X[3][3], double U[3][6]
        , double** T, double* R);
    void RKt_NH(double D[6], double X[3][3], double U[3][6], double** T, double* R);
    void RKt_MooneyRivlin(double D[7], double X[3][3], double U[3][6], double **T, double *R);
    void RKt_Orthotropic(double D[14], double X[3][3], double U[3][6], double** T, double *R);

};

class BeamContact
{
private:
    double epsilon;
    double NEN;
    //Eigen::MatrixXd ELEM;
    int NBEAMS;
    Eigen::MatrixXd GlobalELEM;
    Eigen::VectorXd n;
    double start;
    double precontactiter;
public:
    Eigen::VectorXd set_normal(double* a, double* b);
    double get_precontiter() { return precontactiter; }
    double get_startingpoint() { return start; }
    void set_startingpoint(double value) { start = value; }
    BeamContact(const int nbeams, std::string str);
    double get_penaltyparameter();
    int get_nen();
    int get_Globalelem(int i, int j);
    Eigen::VectorXd set_normal(double* a, double* b, Eigen::VectorXd n);
    Eigen::VectorXd get_normal() { return n; }
    void Contact_NTN(double D[6], double X1[3], double X2[3], double u1[6], double u2[6], double* CR, double** CT, double* g);
    std::vector<std::vector<int>> LocalContactSearch_STS(NLEBBE3D* EBBE3D, Eigen::VectorXd* GU, int nbeams);
    void Contact_STS(double D[4], double X1[3], double X2[3]
        , double X3[3], double X4[3], double u1[6], double u2[6], double u3[6], double u4[6],
        double** CT, double* CR, double* RP, double* g);
    std::vector<std::vector<int>> LocalContactSearch_NTN(NLEBBE3D* EBBE3D, Eigen::VectorXd* GU, int nbeams);
    std::vector<std::vector<int>> LocalContactSearch_STS(VAMBeamElement* VAMBE, Eigen::VectorXd* GU, int nbeams);
    std::vector<std::vector<int>> LocalContactSearch_NTN(VAMBeamElement* VAMBE, Eigen::VectorXd* GU, int nbeams);
};

class Boundary
{
private:
    Eigen::MatrixXd ELEM;
public:
    Boundary();
    Boundary(int beamnum, int contactswitch);
    Boundary(int beamnum, double val);
    Boundary(int beamnum);
    void SetBoundary(Eigen::SparseMatrix<double, Eigen::ColMajor>* GT, Eigen::VectorXd* GR, NLEBBE3D* EBBE3D, int beamnum, int fiter, int continter, Eigen::VectorXd* GU);
    void SetBoundary(Eigen::SparseMatrix<double, Eigen::ColMajor>* GT, Eigen::VectorXd* GR, NLEBBE3D* EBBE3D, int fiter, int coniter, Eigen::VectorXd* GU);
};

class Loading
{
private:
    Eigen::MatrixXd ELEM;
public:
    Loading();
    Loading(int beamnum, int contactswitch);
    Loading(int beamnum);
    void SetLoad(Eigen::VectorXd* GR, NLEBBE3D* EBBE3D, int beamnum, int fiter, int preconiter);
    void SetLoad(Eigen::VectorXd* GR, NLEBBE3D* EBBE3D, int fiter, int coniter);
    void SetLoad(Eigen::VectorXd* GR, VAMBeamElement* VAMBE, int beamnum, int fiter, int coniter);
};

class VTKGrid
{
private:
    Eigen::MatrixXd Points;
    Eigen::MatrixXd Cells;
    int nelem, nnode, ntype, eletype;
    std::string dim;
    Eigen::MatrixXd U, theta, F, M;
    Eigen::MatrixXd CP;
    Eigen::MatrixXd gap;
public:
    VTKGrid(NLEBBE3D EBBE3D, int refine, double load, std::string BeamName, Eigen::MatrixXd GU, int ndim, int eletype, int ndof, int ncells, int beamnum, double ** g, double epsilon);
    VTKGrid(NLEBBE2D EBBE2D, double load, std::string BeamName, Eigen::MatrixXd GU, int ndim, int eletype, int ndof, int ncells);
    VTKGrid(NLEBBE3D EBBE3D, int refine, double load, std::string BeamName, Eigen::MatrixXd GU, int ndim, int eletype, int ndof, int ncells, int beamnum);
    VTKGrid(VAMBeamElement VAMBE, int refine, double load, std::string BeamName, Eigen::VectorXd U, int ndim, int eletype, int dof, int ncells, int beamnum, double** g, double epsilon);
    void WriteVTK(int load, std::string BeamName, std::string filepath, int beamnum, double epsilon);
    void WriteVTK_VAM(int load, std::string BeamName, std::string filepath, int beamnum, double epsilon);
    Eigen::MatrixXd FindSurfacePoints(Eigen::VectorXd n1, Eigen::VectorXd node, double radius, int refine);
    Eigen::MatrixXd FindCells(int p1, int p2, int refine);
};

//Functions in BarElement_Linear_1D.cpp
LinearBarElement ReadLBEFile();
void StiffnessMatrix_LBE(double A, double E, double xa, double xb, Eigen::MatrixXd& k);
void LocalForceVec_LBE(double f1, double f2, Eigen::VectorXd& f);
void ApplyConstraints_LBE(Eigen::MatrixXd CNODE, int NNODE, Eigen::VectorXd F,
    Eigen::SparseMatrix<double, Eigen::ColMajor>& K);

//Functions in BarElement_NonLinear_1D.cpp
NonLinearBarElement ReadNLBEFile();
void TangentStiffnessMatrix_NLBE(Eigen::MatrixXd& t, Eigen::VectorXd U, int StartNode, int EndNode, double h);
void StiffnessMatrix_NLBE(Eigen::MatrixXd& k, Eigen::VectorXd U, int StartNode, int EndNode, double h);
void LocalForceVec_NLBE(Eigen::VectorXd& f, int StartNode, int EndNode, double h);
void ApplyConstraints_NLBE(Eigen::SparseMatrix<double, Eigen::RowMajor>& T, Eigen::VectorXd& U, Eigen::VectorXd& R, int NNODE, Eigen::MatrixXd CNODE);

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

//Textile microstructure 
void TextileMicrostructureGen(NLEBBE3D* EBBE3D, int nbeams);
void ReadTextileInp(VAMBeamElement* VAMBE, int nls, int ndof, int ndim);
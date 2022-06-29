#include<Eigen/Dense>
#include<Eigen/Sparse>

struct LOAD
{
    char* name[50];
    double loadval;
};

struct ELSET
{
    char* elset;
    Eigen::VectorXd elements;
};

struct Mesh
{
//Number of Materials
    int NMAT;
//No. of load steps
    int NLS;
//Number of bodies
    int NBODIES;
//Number of degrees of freedom
    int NDOF;
//Nodal Coordinates
    Eigen::MatrixXd NODE;
//Number of Nodes
    int NNODE;
//Number of Elements
    int NELEM;
//Dirichlet Boundary conditions
    Eigen::MatrixXd CNODE;
//Neumann Boundary conditions
    Eigen::MatrixXd FORCE;
//Material constants
    Eigen::MatrixXd MAT;
//Geometry
    double b1,h1;
//Geometry
    double b2,h2;
    Eigen::MatrixXd ELEM;
    //Element used
    //1. Euler Benoulli Beam
    int choice;
//Loading
    //struct LOAD load;
    //struct ELSET elset;
    double vf;
    double af;
};


/*struct GaussPoints
{
    double x1;
    double x2;
    double x3;
    double w1;
    double w2;
    double w3;
};*/

//Functions in ReadInpFile.cpp
Mesh ReadInpFile();


//Functions in EulerBernoulli2D.cpp
double Jacobian(double xa, double xb);
double HC(double exi, int num, double xa, double xb);
double DDHC(double exi, int num, double xa, double xb);
double SDHC(double exi, int num, double xa, double xb);
double HF(double exi, double num);
double SDHF(double num);
void RearrangeElementMatrices(Eigen::MatrixXd &k, Eigen::MatrixXd &t, Eigen::VectorXd &f);
double dw0dx(Eigen::VectorXd &U, double x, int a, int b,double xa,double xb);
double dU0dx(Eigen::VectorXd &U, int a, int b);
void TangentStiffnessEulerBernoulli(Eigen::MatrixXd &t, Eigen::MatrixXd &k, double xa, double xb,double E, double nu, double base, double height, Eigen::VectorXd &U, int a, int b);
void StiffnessEulerBernoulli(Eigen::MatrixXd &k, double xa, double xb,double E, double nu, double base, double height, Eigen::VectorXd &U, int a, int b);
void ForceVec(Eigen::VectorXd &f, double xa, double xb, int a, int b, double vf, double af);
void ApplyConstraints2DEB(Eigen::SparseMatrix<double, Eigen::ColMajor>& T, Eigen::VectorXd& U, Eigen::MatrixXd &CNODE, int NNODE, Eigen::VectorXd& R);
//GaussPoints GaussQuadraturePoints();

//Functions in BarElement.cpp
void NonLinearStiffnessMatrix(Eigen::MatrixXd& k,Eigen::VectorXd U, int StartNode, int EndNode, double h);
void TangentStiffnessMatrix(Eigen::MatrixXd& t, Eigen::VectorXd U, int StartNode, int EndNode, double h);
void ForceVector(Eigen::VectorXd& f, int StartNode, int EndNode, double h);
void ApplyConstraints(Eigen::SparseMatrix<double, Eigen::RowMajor>& T, Eigen::VectorXd &U, Eigen::VectorXd &R, int NNODE, int cNODE, double value);

//Functions for VAM Beam Element


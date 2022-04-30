#include<Eigen/Dense>
#include<Eigen/Sparse>

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
//Max no. of iterations
    int maxiter;
    int iter, fiter;
};

Mesh ReadInpFile();
void StiffnessEulerBernoulli(Eigen::MatrixXd &k, double xa, double xb,double E, double nu, double b1, double h1, Eigen::VectorXd U, int a, int b);
void TangentStiffnessEulerBernoulli(Eigen::MatrixXd &t, Eigen::MatrixXd &k, double xa, double xb,double E, double nu, double base, double height, Eigen::VectorXd U, int a, int b);
void RearrangeElementStiffness(Eigen::MatrixXd &k);

//void BarElement(double A, double E, double l, Eigen::SparseMatrix<double, Eigen::RowMajor>& K,std::vector<Eigen::Triplet<double>>matTriplets, int StartNode, int EndNode);
//f Vector
void ForceVector(Eigen::VectorXd& F, int StartNode, int EndNode, double h);
//Local Stiffness matrix (non linear case)
void NonLinearStiffnessMatrix(Eigen::MatrixXd& k,Eigen::VectorXd U, int StartNode, int EndNode, double h);
//Apply boundary conditions
void ApplyConstraints(Eigen::SparseMatrix<double, Eigen::RowMajor>& T, Eigen::VectorXd &U, Eigen::VectorXd &R, int NNODE, int cNODE, double value);
//Tangent Stiffness Matrix
void TangentStiffnessMatrix(Eigen::MatrixXd& t, Eigen::VectorXd U, int StartNode, int EndNode, double h);

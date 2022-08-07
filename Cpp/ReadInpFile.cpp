#ifndef VARIABLES_H
#define VARIABLES_H
#include "Variables.h"

//Mesh ReadInpFile()
//{
//    struct Mesh M;
//    M.NBODIES = 2;
//    M.NNODE = 12;
//    M.NELEM = 10;
//    M.NMAT = 2;
//    M.NDOF = 3;
//    M.NLS = 20;

//    M.NODE = Eigen::MatrixXd::Zero(M.NNODE,3);
//    M.ELEM = Eigen::MatrixXd::Zero(M.NELEM, 3);
//    M.MAT = Eigen::MatrixXd::Zero(M.NMAT, M.NBODIES);
//    M.CNODE = Eigen::MatrixXd::Zero(6,2);
//    M.LOAD = Eigen::MatrixXd::Zero(2,3);
////    Nodal Information
//    M.NODE(0,0) = 0;
//    M.NODE(0,1) = 0;
//    M.NODE(0,2) = 0;
//    M.NODE(1,0) = 0.508/2;
//    M.NODE(1,1) = 0;
//    M.NODE(1,2) = 0;
//    M.NODE(2,0) = 1.016/2;
//    M.NODE(2,1) = 0;
//    M.NODE(2,2) = 0;
//    M.NODE(3,0) = 1.524/2;
//    M.NODE(3,1) = 0;
//    M.NODE(3,2) = 0;
//    M.NODE(4,0) = 2.032/2;
//    M.NODE(4,1) = 0;
//    M.NODE(4,2) = 0;
//    M.NODE(5,0) = 2.54/2;
//    M.NODE(5,1) = 0;
//    M.NODE(5,2) = 0;
//    M.NODE(6,0) = 0;
//    M.NODE(6,1) = 1;
//    M.NODE(6,2) = 0;
//    M.NODE(7,0) = 0.508/2;
//    M.NODE(7,1) = 1;
//    M.NODE(7,2) = 0;
//    M.NODE(8,0) = 1.016/2;
//    M.NODE(8,1) = 1;
//    M.NODE(8,2) = 0;
//    M.NODE(9,0) = 1.524/2;
//    M.NODE(9,1) = 1;
//    M.NODE(9,2) = 0;
//    M.NODE(10,0) = 2.032/2;
//    M.NODE(10,1) = 1;
//    M.NODE(10,2) = 0;
//    M.NODE(11,0) = 2.54/2;
//    M.NODE(11,1) = 1;
//    M.NODE(11,2) = 0;
////    Connectivity information
//    M.ELEM(0,0) = 1;
//    M.ELEM(0,1) = 1;
//    M.ELEM(0,2) = 2;
//    M.ELEM(1,0) = 1;
//    M.ELEM(1,1) = 2;
//    M.ELEM(1,2) = 3;
//    M.ELEM(2,0) = 1;
//    M.ELEM(2,1) = 3;
//    M.ELEM(2,2) = 4;
//    M.ELEM(3,0) = 1;
//    M.ELEM(3,1) = 4;
//    M.ELEM(3,2) = 5;
//    M.ELEM(4,0) = 1;
//    M.ELEM(4,1) = 5;
//    M.ELEM(4,2) = 6;
//    M.ELEM(5,0) = 2;
//    M.ELEM(5,1) = 7;
//    M.ELEM(5,2) = 8;
//    M.ELEM(6,0) = 2;
//    M.ELEM(6,1) = 8;
//    M.ELEM(6,2) = 9;
//    M.ELEM(7,0) = 2;
//    M.ELEM(7,1) = 9;
//    M.ELEM(7,2) = 10;
//    M.ELEM(8,0) = 2;
//    M.ELEM(8,1) = 10;
//    M.ELEM(8,2) = 11;
//    M.ELEM(9,0) = 2;
//    M.ELEM(9,1) = 11;
//    M.ELEM(9,2) = 12;
////    Beam material information
//    M.MAT(0,0) = 2.068423e+11;
//    M.MAT(0,1) = 0.33;
//    M.MAT(1,0) = 2.068423e+11;
//    M.MAT(1,1) = 0.33;
////    Beam cross section information
//    M.b1 = 0.0254;
//    M.h1 = 0.0254;
//    M.b2 = 0.0254;
//    M.h2 = 0.0254;
////    Boundary conditions and constrainst
//    //There are 3 dof for each node. Numbering followed is as follows:
//    //Dispalcement
//        //0. Horizontal deflection
//        //1. Vertical Deflection
//        //2. Angle of rotation
//    //Force/Neumann
//        //1. Axial Force
//        //2. Vertical Force
//        //3. Moment
//    //Displacement Boundary
//    M.CNODE(0,0) = 1;
//    M.CNODE(0,1) = 1;
//    M.CNODE(1,0) = 6;
//    M.CNODE(1,1) = 0;
//    M.CNODE(2,0) = 6;
//    M.CNODE(2,1) = 2;
//    M.CNODE(3,0) = 7;
//    M.CNODE(3,1) = 1;
//    M.CNODE(4,0) = 12;
//    M.CNODE(4,1) = 0;
//    M.CNODE(5,0) = 12;
//    M.CNODE(5,1) = 2;
////Loading conditions
//    //Element sets have to be written in inp file for distributed load
////    Load can act in 3 ways depending on the no. of degrees of freedom
////    0. Axial Force
////    1. Vertical Force
////    2. Moment
////    Element sets are described for each load indicating the nodes on which load act
//  /*  M.LOAD(0,0) = 0;
//    M.LOAD(0,1) = -175.126835;
//    M.LOAD(0,2) = 0;
//    M.LOAD(1,0) = 0;
//    M.LOAD(1,1) = -175.126835;
//    M.LOAD(1,2) = 0;
//    M.choice = 3;
////Element Sets
//    M.ELSET(0,0) = 1;
//    M.ELSET(0,1) = 2;
//    M.ELSET(0,2) = 3;
//    M.ELSET(0,3) = 4;
//    M.ELSET(0,4) = 5;
//    M.ELSET(0,5) = 6;
//    M.ELSET(1,0) = 7;
//    M.ELSET(1,1) = 8;
//    M.ELSET(1,2) = 9;
//    M.ELSET(1,3) = 10;
//    M.ELSET(1,4) = 11;
//    M.ELSET(1,5) = 12;*/
//    M.choice = 3;
//    return M;
//}

/*void LoadVector(Eigen::MatrixXd force, Eigen::MatrixXd& ELSET, Eigen::MatrixXd &LOAD, int NNODE)
{
    //Initialie load vector to zero
    for(int i=0;i<NNODE;i++)
    {
        LOAD(i,0) = 0;
        LOAD(i,1) = 0;
        LOAD(i,2) = 0;
    }
//    initialize load vector with the values of load given by the user
    for(int j=0;j<ELSET.rows();j++)
        for(int i=0;i<ELSET.cols();i++)
        {
            LOAD((int)ELSET(j,i)-1,0) = force(j,0);
            LOAD((int)ELSET(j,i)-1,1) = force(j,1);
            LOAD((int)ELSET(j,i)-1,2) = force(j,2);
        }
}
*/
/*
GaussPoints GaussQuadraturePoints()
{
    struct GaussPoints GP;
    //Three Point
    GP.w1 = 5/9;
    GP.w2 = 8/9;
    GP.w3 = 5/9;
    GP.x1 = -sqrt(3/5);
    GP.x2 = 0;
    GP.x3 = sqrt(3/5);

    return GP;
}*/

/*
LinearBarElement ReadLBEFile()
{
    struct LinearBarElement LBE;

    LBE.NBODIES = 1;
    LBE.NNODE = 5;
    LBE.NELEM = 4;
    LBE.NMAT = 2;
    LBE.NDOF = 3;
    \
    LBE.NODE = Eigen::MatrixXd::Zero(LBE.NNODE,3);
    LBE.ELEM = Eigen::MatrixXd::Zero(LBE.NELEM, 3);

    LBE.E = 30e6;
    //LBE.nu = 0.33;

    LBE.CS.height = 1;
    LBE.CS.height = 1;

    LBE.NODE(0,0) = 0;
    LBE.NODE(0,1) = 0;
    LBE.NODE(0,2) = 0;
    LBE.NODE(1,0) = 63.5e-3;
    LBE.NODE(1,1) = 0;
    LBE.NODE(1,2) = 0;
    LBE.NODE(2,0) = 127e-3;
    LBE.NODE(2,1) = 0;
    LBE.NODE(2,2) = 0;
    LBE.NODE(3,0) = 190.5e-3;
    LBE.NODE(3,1) = 0;
    LBE.NODE(3,2) = 0;
    LBE.NODE(4,0) = 254e-3;;
    LBE.NODE(4,1) = 0;
    LBE.NODE(4,2) = 0;

    LBE.ELEM(0,0) = 1;
    LBE.ELEM(0,1) = 1;
    LBE.ELEM(0,2) = 2;
    LBE.ELEM(1,0) = 1;
    LBE.ELEM(1,1) = 2;
    LBE.ELEM(1,2) = 3;
    LBE.ELEM(2,0) = 1;
    LBE.ELEM(2,1) = 3;
    LBE.ELEM(2,2) = 4;
    LBE.ELEM(3,0) = 2;
    LBE.ELEM(3,1) = 4;
    LBE.ELEM(3,2) = 5;

    LBE.CNODE(0) = 1;

    LBE.LOAD.cload[0] = 0;
    LBE.LOAD.cload[1] = 0;
    LBE.LOAD.cload[2] = 0;
    LBE.LOAD.cload[3] = 0;
    LBE.LOAD.cload[4] = 1;

    return LBE;
}

NonLinearBarElement ReadNLBEFile()
{
    struct NonLinearBarElement NLBE;

    NLBE.NBODIES = 1;
    NLBE.NNODE = 5;
    NLBE.NELEM = 4;
    NLBE.NMAT = 2;
    NLBE.NDOF = 3;
    \
    NLBE.NODE = Eigen::MatrixXd::Zero(NLBE.NNODE,3);
    NLBE.ELEM = Eigen::MatrixXd::Zero(NLBE.NELEM, 3);

    NLBE.E = 30e6;
    NLBE.MAT[0].nu = 0.33;
    NLBE.MAT[1].E = 200e9;
    NLBE.MAT[1].nu = 0.3;

    NLBE.CS.height = 1;
    NLBE.CS.height = 1;

    NLBE.NODE(0,0) = 0;
    NLBE.NODE(0,1) = 0;
    NLBE.NODE(0,2) = 0;
    NLBE.NODE(1,0) = 63.5e-3;
    NLBE.NODE(1,1) = 0;
    NLBE.NODE(1,2) = 0;
    NLBE.NODE(2,0) = 127e-3;
    NLBE.NODE(2,1) = 0;
    NLBE.NODE(2,2) = 0;
    NLBE.NODE(3,0) = 190.5e-3;
    NLBE.NODE(3,1) = 0;
    NLBE.NODE(3,2) = 0;
    NLBE.NODE(4,0) = 254e-3;;
    NLBE.NODE(4,1) = 0;
    NLBE.NODE(4,2) = 0;

    NLBE.ELEM(0,0) = 1;
    NLBE.ELEM(0,1) = 1;
    NLBE.ELEM(0,2) = 2;
    NLBE.ELEM(1,0) = 1;
    NLBE.ELEM(1,1) = 2;
    NLBE.ELEM(1,2) = 3;
    NLBE.ELEM(2,0) = 1;
    NLBE.ELEM(2,1) = 3;
    NLBE.ELEM(2,2) = 4;
    NLBE.ELEM(3,0) = 2;
    NLBE.ELEM(3,1) = 4;
    NLBE.ELEM(3,2) = 5;

    NLBE.CNODE(0) = 1;

    NLBE.LOAD.cload[0] = 0;
    NLBE.LOAD.cload[1] = 0;
    NLBE.LOAD.cload[2] = 0;
    NLBE.LOAD.cload[3] = 0;
    NLBE.LOAD.cload[4] = 1;

    return NLBE;
}

NonLinearEulerBernouliBeamElement ReadNLEBBEFile()
{
    struct NonLinearEulerBernouliBeamElement NLEBBE;

    NLEBBE.NBODIES = 1;
    NLEBBE.NNODE = 5;
    NLEBBE.NELEM = 4;
    NLEBBE.NMAT = 2;
    NLEBBE.NDOF = 3;
    \
    NLEBBE.NODE = Eigen::MatrixXd::Zero(NLEBBE.NNODE,3);
    NLEBBE.ELEM = Eigen::MatrixXd::Zero(NLEBBE.NELEM, 3);

    NLEBBE.MAT[0].E = 30e6;
    NLEBBE.MAT[0].nu = 0.33;
    NLEBBE.MAT[1].E = 200e9;
    NLEBBE.MAT[1].nu = 0.3;

    NLEBBE.CS.height = 1;
    NLEBBE.CS.width = 1;

    NLEBBE.NODE(0,0) = 0;
    NLEBBE.NODE(0,1) = 0;
    NLEBBE.NODE(0,2) = 0;
    NLEBBE.NODE(1,0) = 63.5e-3;
    NLEBBE.NODE(1,1) = 0;
    NLEBBE.NODE(1,2) = 0;
    NLEBBE.NODE(2,0) = 127e-3;
    NLEBBE.NODE(2,1) = 0;
    NLEBBE.NODE(2,2) = 0;
    NLEBBE.NODE(3,0) = 190.5e-3;
    NLEBBE.NODE(3,1) = 0;
    NLEBBE.NODE(3,2) = 0;
    NLEBBE.NODE(4,0) = 254e-3;;
    NLEBBE.NODE(4,1) = 0;
    NLEBBE.NODE(4,2) = 0;

    NLEBBE.ELEM(0,0) = 1;
    NLEBBE.ELEM(0,1) = 1;
    NLEBBE.ELEM(0,2) = 2;
    NLEBBE.ELEM(1,0) = 1;
    NLEBBE.ELEM(1,1) = 2;
    NLEBBE.ELEM(1,2) = 3;
    NLEBBE.ELEM(2,0) = 1;
    NLEBBE.ELEM(2,1) = 3;
    NLEBBE.ELEM(2,2) = 4;
    NLEBBE.ELEM(3,0) = 2;
    NLEBBE.ELEM(3,1) = 4;
    NLEBBE.ELEM(3,2) = 5;

    NLEBBE.CNODE(0) = 1;

    NLEBBE.LOAD.cload[0] = 0;
    NLEBBE.LOAD.cload[1] = 0;
    NLEBBE.LOAD.cload[2] = 0;
    NLEBBE.LOAD.cload[3] = 0;
    NLEBBE.LOAD.cload[4] = 1;

    return NLEBBE;
}*/

VAMBeamElement ReadVAMBEFile()
{
    struct VAMBeamElement M;
    M.NBODIES = 1;
    M.NNODE = 5;
    M.NELEM = 4;
    M.NMAT = 1;
    M.NDOF = 4;
    M.NLS = 20;
    M.CMP.np = 6;
    M.NCS = 1;

    M.NODE = Eigen::MatrixXd::Zero(M.NNODE,3);
    M.ELEM = Eigen::MatrixXd::Zero(M.NELEM, 3);
    M.CMP.Orient = Eigen::VectorXd::Zero(M.CMP.np);
    M.CMP.inittwist = Eigen::VectorXd::Zero(3);
    //M.LOAD = Eigen::MatrixXd::Zero(1,3);
//    Nodal Information
    M.NODE(0,0) = 0;
    M.NODE(0,1) = 0;
    M.NODE(0,2) = 0;
    M.NODE(1,0) = 63.5e-3;
    M.NODE(1,1) = 0;
    M.NODE(1,2) = 0;
    M.NODE(2,0) = 127e-3;
    M.NODE(2,1) = 0;
    M.NODE(2,2) = 0;
    M.NODE(3,0) = 190.5e-3;
    M.NODE(3,1) = 0;
    M.NODE(3,2) = 0;
    M.NODE(4,0) = 254e-3;;
    M.NODE(4,1) = 0;
    M.NODE(4,2) = 0;
//    Connectivity information
    M.ELEM(0,0) = 1;
    M.ELEM(0,1) = 1;
    M.ELEM(0,2) = 2;
    M.ELEM(1,0) = 1;
    M.ELEM(1,1) = 2;
    M.ELEM(1,2) = 3;
    M.ELEM(2,0) = 1;
    M.ELEM(2,1) = 3;
    M.ELEM(2,2) = 4;
    M.ELEM(3,0) = 1;
    M.ELEM(3,1) = 4;
    M.ELEM(3,2) = 5;
//    Beam material information
    //M.MAT(0,0) = 30e6;
    //M.MAT(0,1) = 0.33;
    M.CMP.E11 = 135.6e9;
    M.CMP.E22 = 9.9e9;
    M.CMP.E33 = 9.9e9;
    M.CMP.G12 = 4.2e9;
    M.CMP.G13 = 4.2e9;
    M.CMP.G23 = 3.3e9;
    M.CMP.nu12 = 0.3;
    M.CMP.nu13 = 0.3;
    M.CMP.nu23 = 0.5;

    M.CMP.Orient(0) = 20;
    M.CMP.Orient(1) = -70;
    M.CMP.Orient(2) = 20;
    M.CMP.Orient(3) = -20;
    M.CMP.Orient(4) = 70;
    M.CMP.Orient(5) = -20;

    M.CMP.inittwist(0) = -0.20;
    M.CMP.inittwist(1) = 0;
    M.CMP.inittwist(2) = 0;

//    Beam cross section information
    M.CS.width = 25.4e-3;
    M.CS.height = 1.168e-3;

    //Loading
    M.B.F1 = Eigen::VectorXd::Zero(3);
    M.B.FN = Eigen::VectorXd::Zero(3);
    M.B.M1 = Eigen::VectorXd::Zero(3);
    M.B.MN = Eigen::VectorXd::Zero(3);

    M.B.u1 = Eigen::VectorXd::Zero(3);
    M.B.uN = Eigen::VectorXd::Zero(3);
    M.B.theta1 = Eigen::VectorXd::Zero(3);
    M.B.thetaN = Eigen::VectorXd::Zero(3);

    for(int i=0;i<3;i++)
    {
        M.B.F1(i) = 0;
        M.B.M1(i) = 0;
        M.B.MN(i) = 0;
        M.B.u1(i) = 0;


    }

    return M;
}

/*
Mesh ReadInpFile(std::string filename)
{
    struct Mesh M;
    std::fstream file;
    file.open(filename.c_str());
    if(file.is_open())
    {
        //Read the file
        while(!file.eof())
        {
            std::string str;
            getline(file,str);
            std::stringstream s(str);
            s>>str;
            //Type of analysis
            if(std::strcmp(str.c_str(),"*STATIC"))
            {
                getline(file,str);
                M.choice = stoi(str);
            }
            else
            {
                std::cout<<"You seem to have forgotten to enter the kind of analysis that you want me to do."<<std::endl;
                std::cout<<"Use the numbers listed below to tell me about the same"<<std::endl;
                std::cout<<"1. Linear Static Analysis using Bar Elements"<<std::endl;
                std::cout<<"2. Non Linear Static Anlaysis using Bar Elements"<<std::endl;
                std::cout<<"3. Non Linear Static Analysis using Euler Bernoulli Beam Elements"<<std::endl;
                std::cout<<"4. Non Linear Static Analysis using VAM Beam Elements"<<std::endl;
                std::cout<<"Check the sample .inp file to know more."<<std::endl;

            }
            //Number of Bodies
            if(std::strcmp(str.c_str(),"*NBODIES"))
            {
                getline(file, str);
                M.NBODIES = stoi(str);
            }
            else
            {
                std::cout<<"You seem to have forgotten to enter the number of bodies to taken for me to solve."<<std::endl;
                std::cout<<"I'm assuming there's only one body"<<std::endl;
                M.NBODIES = 1;
            }
            //Get nodal data.
            if(std::strcmp(str.c_str(),"*NODE"))
            {
                std::cout<<"Accessing Nodal Data"<<std::endl;
                int i=0;
                do
                {
                    getline(file, str);
                    std::stringstream s(str);
                    s>>M.NODE(i,0)>>M.NODE(i,1)>>M.NODE(i,3);
                    i++;
                } while (!std::strcmp(str.c_str()," "));
            }
            else
            {
                std::cout<<"I dont' see any nodal data"<<std::endl;
                std::cout<<"Check the sample .inp file to know more"<<std::endl;
            }
            //Get Elemental Data
            if(std::strcmp(str.c_str(), "*ELEMENT"))
            {
                std::cout<<"Accessing Connectivity Data"<<std::endl;
                int i=0;
                do
                {
                    getline(file, str);
                    std::stringstream s(str);
                    s>>M.ELEM(i,0)>>M.ELEM(i,1)>>M.ELEM(i,2);
                    i++;
                } while (!std::strcmp(str.c_str()," "));

            }
            else
            {
                std::cout<<"I dont' see any connectivity data"<<std::endl;
                std::cout<<"Check the sample .inp file to know more"<<std::endl;
            }
            //Get Material Parameters
            if(std::strcmp(str.c_str(),"*MATERIAL"))
            {
                do
                {
                    if(std::strcmp(str.c_str(),"NMAT"))
                    {
                        getline(file, str);
                        std::stringstream s(str);
                        s>>M.NMAT;
                    }
                    std::string temp;
                    getline(file, str);
                    std::stringstream s(str);
                    s>>temp;
                    if(strcmp(temp.c_str(),"*ELASTIC"))
                    {
                        s>>temp>>temp>>temp;
                        if(strcmp(temp.c_str(),"ISOTROPIC"))
                        {
                            int i=0;
                            do
                            {
                                getline(file, str);
                                std::stringstream s(str);
                                s>>M.Mat[i].Iso.E>>M.Mat[i].Iso.nu;
                                i++;
                            }
                            while (!std::strcmp(str.c_str()," "));

                        }
                        else if(strcmp(temp.c_str(),"ORTHOTROPIC"))
                        {
                            getline(file, str);
                            std::stringstream s(str);
                        }

                    }


                }while (!std::strcmp(str.c_str()," "));
            }
            else
            {

            }
            //Get Loading and Boundary conditions
            //Element Sets

        }
    }
    else
    {
        std::cout<<"File not found"<<std::endl;
        std::cout<<"Maybe there's a spelling error"<<std::endl;
        std::cout<<"Maybe the slash is reversed"<<std::endl;
    }

}*/

#endif

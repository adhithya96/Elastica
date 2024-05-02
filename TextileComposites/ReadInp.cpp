#ifndef VARIABLES_H
#define VARIABLES_H
#include "Variables.h"

LinearBarElement ReadLBEFile()
{
    struct LinearBarElement LBE;

    LBE.NBODIES = 1;
    LBE.NNODE = 5;
    LBE.NELEM = 4;
    LBE.NMAT = 2;
    LBE.NDOF = 1;

    LBE.NODE = Eigen::MatrixXd::Zero(LBE.NNODE, LBE.NDOF);
    LBE.ELEM = Eigen::MatrixXd::Zero(LBE.NELEM, 3);
    LBE.MAT = Eigen::MatrixXd::Zero(LBE.NMAT, 2);

    LBE.MAT(0, 0) = 30e6;
    LBE.MAT(0, 1) = 0.33;
    LBE.MAT(1, 0) = 30e6;
    LBE.MAT(1, 1) = 0.33;

    LBE.CS.choice = "RECT";
    LBE.CS.Rect.height = 1;
    LBE.CS.Rect.height = 1;

    LBE.NODE(0, 0) = 0;
    LBE.NODE(1, 0) = 0.25;
    LBE.NODE(2, 0) = 0.5;
    LBE.NODE(3, 0) = 0.75;
    LBE.NODE(4, 0) = 1;

    LBE.ELEM(0, 0) = 1;
    LBE.ELEM(0, 1) = 1;
    LBE.ELEM(0, 2) = 2;
    LBE.ELEM(1, 0) = 1;
    LBE.ELEM(1, 1) = 2;
    LBE.ELEM(1, 2) = 3;
    LBE.ELEM(2, 0) = 1;
    LBE.ELEM(2, 1) = 3;
    LBE.ELEM(2, 2) = 4;
    LBE.ELEM(3, 0) = 2;
    LBE.ELEM(3, 1) = 4;
    LBE.ELEM(3, 2) = 5;

    LBE.CNODE = Eigen::MatrixXd::Zero(1, 2);
    LBE.CNODE(0, 0) = 1;
    LBE.CNODE(0, 1) = 0;


    LBE.LOAD = Eigen::MatrixXd::Zero(1, 2);
    LBE.LOAD(0, 0) = 5;
    LBE.LOAD(0, 1) = 1;

    return LBE;
}

NonLinearBarElement ReadNLBEFile()
{
    struct NonLinearBarElement NLBE;

    NLBE.NBODIES = 1;
    NLBE.NNODE = 3;
    NLBE.NELEM = 2;
    NLBE.NMAT = 1;
    NLBE.NDOF = 1;

    NLBE.NODE = Eigen::MatrixXd::Zero(NLBE.NNODE, 2);
    NLBE.ELEM = Eigen::MatrixXd::Zero(NLBE.NELEM, 2);

    NLBE.CS.choice = "RECT";
    NLBE.CS.Rect.height = 1;
    NLBE.CS.Rect.width = 1;

    NLBE.NODE(0, 0) = 0;
    NLBE.NODE(0, 1) = 0;
    NLBE.NODE(1, 0) = 0.5;
    NLBE.NODE(1, 1) = 0;
    NLBE.NODE(2, 0) = 1;
    NLBE.NODE(2, 1) = 0;

    NLBE.ELEM(0, 0) = 1;
    NLBE.ELEM(0, 1) = 2;
    NLBE.ELEM(1, 0) = 2;
    NLBE.ELEM(1, 1) = 3;

    NLBE.CNODE = Eigen::MatrixXd::Zero(1, 2);
    NLBE.CNODE(0, 0) = 3;
    NLBE.CNODE(0, 1) = sqrt(2);

    return NLBE;
}

//Validation case for GEBT (12 dofs per node) axial load 0/0 ply.
/*VAMBeamElement ReadVAMBEFile()
{
    struct VAMBeamElement M;
    M.NBODIES = 1;
    M.NNODE = 3;
    M.NELEM = 2;
    M.NMAT = 1;
    M.NDOF = 12;
    M.NLS = 10;

    M.NODE = Eigen::MatrixXd::Zero(M.NNODE, 3);
    M.ELEM = Eigen::MatrixXd::Zero(M.NELEM, 3);

    //    Nodal Information
    M.NODE(0, 0) = 0;
    M.NODE(0, 1) = 0;
    M.NODE(0, 2) = 0;
    M.NODE(1, 0) = 0.127;
    M.NODE(1, 1) = 0;
    M.NODE(1, 2) = 0;
    M.NODE(2, 0) = 0.254;
    M.NODE(2, 1) = 0;
    M.NODE(2, 2) = 0;

    //Connectivity information
    M.ELEM(0, 0) = 1;
    M.ELEM(0, 1) = 1;
    M.ELEM(0, 2) = 2;
    M.ELEM(1, 0) = 1;
    M.ELEM(1, 1) = 2;
    M.ELEM(1, 2) = 3;

    M.inittwist = Eigen::VectorXd::Zero(3);
    M.inittwist(0) = 0;
    M.inittwist(1) = 0;
    M.inittwist(2) = 0;

    M.B.a1 = Eigen::VectorXd::Zero(3);
    M.B.b1 = Eigen::VectorXd::Zero(3);
    M.B.aN = Eigen::VectorXd::Zero(3);
    M.B.bN = Eigen::VectorXd::Zero(3);

    M.B.x1 = 1;
    M.B.y1 = 2;
    M.B.xN = 3;
    M.B.yN = 4;

    M.B.aN(0) = 600;

    return M;
}*/







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

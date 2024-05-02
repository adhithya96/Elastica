#pragma once
#include "Variables.h"

NLEBBE3D::NLEBBE3D(int yn, double ys, double yd, int refine, std::string str, int beamnum, int nls, int nen, int ndm)
{
    double yl = (yn + 1) * ys;
    double delh = ys / refine;
    int nnode = (yn / 2 + 1) * refine + 1;

    this->NNODE = nnode;
    this->NELEM = (nnode - 1) / 2;
    this->NDOF = 6;
    this->NLS = nls;
    this->NEN = nen;
    this->NDM = ndm;

    this->NODE = Eigen::MatrixXd::Zero(nnode, ndm);
    this->ELEM = Eigen::MatrixXd::Zero(this->NELEM, nen + 1);

    this->E = 48383.10;
    this->nu = 0.2;
    this->Bp = 1;
    this->Hp = 1;
    this->Zx = 0;
    this->Zy = 0;
    this->Zz = 1;

    this->DIA = 1;
    if (beamnum <= yn / 2)
    {
        double xcoord;
        if (yn % 4 == 0)
        {
            int sign;
            if ((beamnum - yn / 4) > 0)
                sign = 1;
            else
                sign = -1;
            if (sign == -1)
                xcoord = ((beamnum - yn / 4) * ys + sign * ys / 2);
            else
                xcoord = ((beamnum - yn / 4) * ys - sign * ys / 2);
        }
        else
        {
            xcoord = (beamnum - (yn) / 2 + 1) * ys;
        }

        for (int i = 0; i < nnode; i++)
        {
            this->NODE(i, 0) = xcoord;
            this->NODE(i, 1) = (i - (nnode - 1) / 2) * delh;
            this->NODE(i, 2) = 0;
        }
        for (int i = 0; i < this->NELEM; i++)
        {
            this->ELEM(i, 0) = 1;
            this->ELEM(i, 1) = 2 * i + 1;
            this->ELEM(i, 2) = 2 * i + 2;
            this->ELEM(i, 3) = 2 * i + 3;
        }
    }
    else
    {
        double ycoord;
        if (yn % 4 == 0)
        {
            int sign;
            if ((beamnum - yn / 4) > 0)
                sign = 1;
            else
                sign = -1;
            if (sign == -1)
                ycoord = ((beamnum - yn / 2 - yn / 4) * ys + sign * ys / 2);
            else
                ycoord = ((beamnum - yn / 2 - yn / 4) * ys - sign * ys / 2);
        }
        else
        {
            ycoord = (beamnum - 3 - (yn) / 2 + 1) * ys;
        }

        for (int i = 0; i < nnode; i++)
        {
            this->NODE(i, 0) = (i - (nnode - 1) / 2) * delh;;
            this->NODE(i, 1) = ycoord;
            this->NODE(i, 2) = 0;
        }
        for (int i = 0; i < this->NELEM; i++)
        {
            this->ELEM(i, 0) = 1;
            this->ELEM(i, 1) = 2 * i + 1;
            this->ELEM(i, 2) = 2 * i + 2;
            this->ELEM(i, 3) = 2 * i + 3;
        }
    }
}

void TextileMicrostructureGen(NLEBBE3D* EBBE3D, int nbeams, Eigen::VectorXd GU)
{
    std::fstream file1, file2;
    file1.open("E:/Adhithya/MTech_Res_Thesis/Cpp/ThesisCode/TextileComposites/TextileComposites/TextileMicrostructure.txt", std::fstream::in | std::fstream::out);
    file2.open("E:/Adhithya/MTech_Res_Thesis/Cpp/ThesisCode/TextileComposites/TextileComposites/Displacements.txt", std::fstream::in | std::fstream::out);
    if (file1.is_open())
    {
        file1 << "NBEAMS" << std::endl;
        file1 << nbeams << std::endl;
        int size = 0;
        for (int i = 0; i < nbeams; i++)
        {
            file1 << "BEAM " << i + 1 << std::endl;
            file1 << "DIA " << EBBE3D[i].get_diameter() << std::endl;
            file1 << "NODE " << EBBE3D[i].get_nnode() << std::endl;
            for (int j = 0; j < EBBE3D[i].get_nnode(); j++)
                file1 << EBBE3D[i].get_coordinates(j, 0) << " " << EBBE3D[i].get_coordinates(j, 1) << " " << EBBE3D[i].get_coordinates(j, 2) << std::endl;
            file1 << "ELEM " << EBBE3D[i].get_nnode() - 1 << std::endl;
            for (int j = 0; j < EBBE3D[i].get_nnode() - 1; j++)
                file1 << 1 << " " << j + 1 << " " << j + 2 << std::endl;
        }
    }
    else
    {
        std::cout << "Textile microstructure file is not open" << std::endl;
    }
    if (file2.is_open())
    {
        int size = 0;
        for(int i = 0; i < nbeams; i++)
            size += EBBE3D[i].get_nnode() * EBBE3D[i].get_ndof();
        file2 << "U " << size << std::endl;
        for (int j = 0; j < size; j++)
            file2 << GU(j) << std::endl;
    }
    else
    {
        std::cout << "Displacements are not printed" << std::endl;
    }
}

void ReadTextileInp(VAMBeamElement* VAMBE, int nls, int ndof, int ndim)
{
    std::fstream file1;
    file1.open("E:/Adhithya/MTech_Res_Thesis/Cpp/ThesisCode/TextileComposites/TextileComposites/TextileMicrostructure_2cross2.txt", std::fstream::in | std::fstream::out);
    int beamnum = 0;
    double dia = 0;
    int nbeams, nelem = 0, nnode = 0;
    std::string str;
    getline(file1, str);
    std::stringstream s(str);
    s >> str;
    //Type of analysis
    if (!std::strcmp(str.c_str(), "NBEAMS"))
    {
        getline(file1, str);
        nbeams = stoi(str);
    }
    std::vector<std::vector<int>> conn;
    std::vector<std::vector<double>> node;
    Eigen::VectorXd axis(3);
    while (!(file1.eof()))
    {
        getline(file1, str);
        std::stringstream s(str);
        s >> str;
        //set diameter
        if (!std::strcmp(str.c_str(), "DIA"))
        {
            s >> str;
            dia = stod(str);
            continue;
        }
        if (!std::strcmp(str.c_str(), "AXIS"))
        {
            getline(file1, str);
            std::stringstream s(str);
            for (int i = 0; i < 3; i++)
            {
                s >> str;
                axis(i) = stod(str);
            }
            continue;
        }
        //set nodal data
        if (!std::strcmp(str.c_str(), "NODE"))
        {
            s >> str;
            nnode = stoi(str);
            for (int i = 0; i < nnode; i++)
            {
                getline(file1, str);
                std::stringstream s(str);
                node.push_back(std::vector<double>());
                for (int j = 0; j < 3; j++)
                {
                    s >> str;
                    node[i].push_back(stod(str));
                }
            }
            continue;
        }
        //set connectivity
        if (!std::strcmp(str.c_str(), "ELEM"))
        {
            s >> str;
            nelem = stoi(str);
            for (int i = 0; i < nelem; i++)
            {
                getline(file1, str);
                std::stringstream s(str);
                conn.push_back(std::vector<int>());
                for (int j = 0; j < 3; j++)
                {
                    s >> str;
                    conn[i].push_back(stoi(str));
                }
            }
            continue;
        }
        if (!std::strcmp(str.c_str(), "U"))
        {
            getline(file1, str);
        }
        if (!std::strcmp(str.c_str(), "END"))
        {
            VAMBE[beamnum] = VAMBeamElement(nls, ndof, ndim, nnode, nelem, node, conn, dia, axis);
            beamnum++;
            node.resize(0);
            conn.resize(0);
        }
    }
}

void ReadDisplacements(Eigen::VectorXd* U, int nnode)
{
    std::fstream file1;
    file1.open("E:/Adhithya/MTech_Res_Thesis/Cpp/ThesisCode/TextileComposites/TextileComposites/Displacements_2cross2.txt", std::fstream::in | std::fstream::out);
    std::string str;
    getline(file1, str);
    std::stringstream s(str);
    s >> str;
    int i = 0;
    if (!std::strcmp(str.c_str(), "U"))
    {
        while (!(file1.eof()))
        {
            for (int j = 0; j < 6; j++)
            {
                getline(file1, str);
                (*U)(i * 12 + 6 + j) = stod(str);
            }
            i++;
        }
    }
    std::cout << *U << std::endl;
}

VAMBeamElement::VAMBeamElement(int nls, int ndof, int ndim, int nnode, int nelem, std::vector<std::vector<double>>& node,
    std::vector<std::vector<int>>& conn, double dia, Eigen::VectorXd axis)
{
    this->NNODE = nnode;
    this->NLS = nls;
    this->NDOF = ndof;
    this->NDM = ndim;
    this->NELEM = nelem;
    this->D = dia;

    this->NODE = Eigen::MatrixXd::Zero(nnode, 3);
    this->ELEM = Eigen::MatrixXd::Zero(nelem, 3);
    this->axis = Eigen::VectorXd::Zero(3);

    for (int i = 0; i < node.size(); i++)
        for (int j = 0; j < node[i].size(); j++)
            this->NODE(i, j) = node[i][j];
    for (int i = 0; i < conn.size(); i++)
        for (int j = 0; j < conn[i].size(); j++)
            this->ELEM(i, j) = conn[i][j];

    this->axis = axis;

}
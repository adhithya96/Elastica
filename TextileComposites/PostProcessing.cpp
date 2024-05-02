#pragma once
#include "Variables.h"
#include<iomanip>

void PostProcessing(Eigen::MatrixXd X, Eigen::VectorXd U, double load, std::string BeamElement, int nnode, int ndm, int ndof)
{
    std::fstream file1;
    std::string str;
    std::ostringstream str1;
    str1 << load;
    str = "E:/Adhithya/MTech_Res_Thesis/Cpp/ThesisCode/TextileComposites/TextileComposites/3DBeam/" + BeamElement + "_" + str1.str() + ".txt";
    file1.open(str.c_str(), std::ios::out);

    //std::cout << X.size() << std::endl;
    //std::cout << str << std::endl;

    Eigen::MatrixXd x = Eigen::MatrixXd::Zero(nnode, ndm);

    if (file1.is_open())
    {
        for (int i = 0; i < nnode; i++)
        {
            for (int j = 0; j < ndm; j++)
            {
                x(i, j) = X(i, j) + U(ndof * i + j);
                if (x(i, j) < 1e-10)
                    x(i, j) = 0;
                file1 << std::setprecision(5) << x(i, j) << " ";
            }
            file1 << std::endl;
        }
    }
    else
    {
        std::cout << "Unable to open file for saving current configuration " << std::endl;
    }
}

Eigen::VectorXd CrossProduct(Eigen::VectorXd A, Eigen::VectorXd B)
{
    Eigen::VectorXd C(3);
    C(0) = A(1) * B(2) - A(2) * B(1);
    C(1) = -(A(0) * B(2) - A(2) * B(0));
    C(2) = A(0) * B(1) - A(1) * B(0);

    return C;
}

VTKGrid::VTKGrid(VAMBeamElement VAMBE, int refine, double load, std::string BeamName, Eigen::VectorXd U, int ndim, int eletype, int ndof, int ncells, int beamnum, double** g, double epsilon)
{
    this->nnode = VAMBE.get_nnode();
    this->nelem = VAMBE.get_nelem();
    this->eletype = eletype;
    this->Points = Eigen::MatrixXd::Zero(this->nnode * refine, ndim);
    this->Cells = Eigen::MatrixXd::Zero(this->nelem * refine, ncells);
    this->U = Eigen::MatrixXd::Zero(this->nnode * refine, ndim);
    this->F = Eigen::MatrixXd::Zero(this->nnode * refine, ndim);
    this->M = Eigen::MatrixXd::Zero(this->nnode * refine, ndim);
    this->gap = Eigen::VectorXd::Zero(this->nnode * refine);
    this->theta = Eigen::MatrixXd::Zero(this->nnode * refine, ndim);
    int temp1 = 0;
    for (int i = beamnum; i > 0; i--)
        temp1 += this->nnode * ndof;
    //Find nodes on the surface
    for (int i = 0; i < VTKGrid::nelem; i++)
    {
        int p1 = VAMBE.get_connectivity(i, 1) - 1;
        int p2 = VAMBE.get_connectivity(i, 2) - 1;
        //std::cout << p1 << std::endl;
        //std::cout << p2 << std::endl;
        //std::cout << p3 << std::endl;
        //normal vectors
        Eigen::VectorXd a(3), b(3);
        for (int k = 0; k < 3; k++)
        {
            a(k) = VAMBE.get_coordinates(p1, k) + U(ndof * p1 + k + temp1 + 6);
            b(k) = VAMBE.get_coordinates(p2, k) + U(ndof * p2 + k + temp1 + 6);
        }

        for (int k = 0; k < refine; k++)
            for (int j = 0; j < 3; j++)
            {
                VTKGrid::F(p1 * refine + k, j) = U(ndof * p1 + j + temp1);
                VTKGrid::M(p1 * refine + k, j) = U(ndof * p1 + j + temp1 + 3);
                VTKGrid::U(p1 * refine + k, j) = U(ndof * p1 + j + temp1 + 6);
                VTKGrid::theta(p1 * refine + k, j) = U(ndof * p1 + j + temp1 + 9);
                VTKGrid::F(p2 * refine + k, j) = U(ndof * p2 + j + temp1);
                VTKGrid::M(p2 * refine + k, j) = U(ndof * p2 + j + temp1 + 3);
                VTKGrid::U(p2 * refine + k, j) = U(ndof * p2 + j + temp1 + 6);
                VTKGrid::theta(p2 * refine + k, j) = U(ndof * p2 + j + temp1 + 9);
            }

        for (int k = 0; k < refine; k++)
        {
            VTKGrid::gap(p1 * refine + k) = g[beamnum][p1];
            VTKGrid::gap(p2 * refine + k) = g[beamnum][p2];
        }

        //std::cout << a << std::endl;
        //std::cout << b << std::endl;
        //std::cout << c << std::endl;
        Eigen::VectorXd n1(3), n2(3);
        n1 = b - a;
        //std::cout << n1 << std::endl;
        //std::cout << n2 << std::endl;
        double r = VAMBE.get_diameter() / 2.0;
        if (i != VTKGrid::nelem - 1)
        {
            Eigen::MatrixXd temp = VTKGrid::FindSurfacePoints(n1, a, r, refine);
            for (int k = 0; k < refine; k++)
                for (int j = 0; j < 3; j++)
                    VTKGrid::Points(p1 * refine + k, j) = temp(k, j);
                    //std::cout << VTKGrid::Points(p1 * refine, j) << std::endl;
        }
        else
        {
            Eigen::MatrixXd temp = VTKGrid::FindSurfacePoints(n1, a, r, refine);
            for (int k = 0; k < refine; k++)
                for (int j = 0; j < 3; j++)
                    VTKGrid::Points(p1 * refine + k, j) = temp(k, j);
            temp = VTKGrid::FindSurfacePoints(n1, b, r, refine);
            for (int k = 0; k < refine; k++)
                for (int j = 0; j < 3; j++)
                    VTKGrid::Points(p2 * refine + k, j) = temp(k, j);
        }
        Eigen::MatrixXd temp = VTKGrid::FindCells(p1, p2, refine);
        for (int k = 0; k < refine; k++)
        {
            for (int j = 0; j < 4; j++)
                VTKGrid::Cells(p1 * refine + k, j) = temp(k, j);
        }
    }

    VTKGrid::WriteVTK_VAM(load, BeamName, ".\\3DBeam\\", beamnum, epsilon);
}

void VTKGrid::WriteVTK_VAM(int load, std::string BeamName, std::string filepath, int beamnum, double epsilon)
{
    std::fstream fname;
    std::ostringstream str1;
    str1 << load;
    std::string str = filepath + BeamName + std::to_string(beamnum) + "_" + str1.str() + ".vtu";
    fname.open(str.c_str(), std::ios::out);

    if (fname.is_open())
    {
        fname << "# vtk DataFile Version 2.0" << std::endl;
        fname << "Beam Element " << BeamName << " load " << str1.str() << std::endl;
        fname << "ASCII " << std::endl;
        fname << "DATASET " << "UNSTRUCTURED_GRID " << std::endl;
        fname << "POINTS " << VTKGrid::Points.rows() << " double" << std::endl;
        for (int i = 0; i < VTKGrid::Points.rows(); i++)
        {
            for (int j = 0; j < VTKGrid::Points.cols(); j++)
            {
                fname << VTKGrid::Points(i, j) << " ";
                //std::cout << VTKGrid::Points(i, j) << std::endl;
            }
            fname << std::endl;
        }
        fname << "CELLS " << VTKGrid::Cells.rows() << " " << VTKGrid::Cells.rows() * (VTKGrid::Cells.cols() + 1) << std::endl;
        for (int i = 0; i < VTKGrid::Cells.rows(); i++)
        {
            fname << 4 << " ";
            for (int j = 0; j < VTKGrid::Cells.cols(); j++)
                fname << VTKGrid::Cells(i, j) << " ";
            fname << std::endl;
        }
        fname << "CELL_TYPES " << VTKGrid::Cells.rows() << std::endl;
        for (int i = 0; i < VTKGrid::Cells.rows(); i++)
            fname << VTKGrid::eletype << std::endl;

        fname << "POINT_DATA " << VTKGrid::Points.rows() << std::endl;
        fname << "VECTORS " << "U " << "double" << std::endl;
        for (int i = 0; i < VTKGrid::U.rows(); i++)
        {
            for (int j = 0; j < VTKGrid::U.cols(); j++)
                fname << VTKGrid::U(i, j) << " ";
            fname << std::endl;
        }
        fname << "VECTORS " << "Theta " << "double" << std::endl;
        for (int i = 0; i < VTKGrid::theta.rows(); i++)
        {
            for (int j = 0; j < VTKGrid::theta.cols(); j++)
                fname << VTKGrid::theta(i, j) << " ";
            fname << std::endl;
        }
        fname << "VECTORS " << "Force " << "double" << std::endl;
        for (int i = 0; i < VTKGrid::F.rows(); i++)
        {
            for (int j = 0; j < VTKGrid::F.cols(); j++)
                fname << VTKGrid::F(i, j) << " ";
            fname << std::endl;
        }
        fname << "VECTORS " << "Moment " << "double" << std::endl;
        for (int i = 0; i < VTKGrid::M.rows(); i++)
        {
            for (int j = 0; j < VTKGrid::M.cols(); j++)
                fname << VTKGrid::M(i, j) << " ";
            fname << std::endl;
        }
        fname << "SCALARS " << "CP " << "double " << 1 << std::endl;
        fname << "LOOKUP_TABLE " << "default " << std::endl;
        for (int i = 0; i < VTKGrid::gap.size(); i++)
            fname << VTKGrid::gap(i) * epsilon << std::endl;
    }
    else
    {
        std::cout << "Unable to open file for saving current configuration " << std::endl;
    }

}




void VTKGrid::WriteVTK(int load, std::string BeamName, std::string filepath, int beamnum, double epsilon)
{
    std::fstream fname;
    std::ostringstream str1;
    str1 << load;
    std::string str = filepath + BeamName + std::to_string(beamnum) + "_" + str1.str() + ".vtu";
    fname.open(str.c_str(), std::ios::out);
    
    if (fname.is_open())
    {
        fname << "# vtk DataFile Version 2.0" << std::endl;
        fname << "Beam Element " << BeamName << " load " << str1.str() << std::endl;
        fname << "ASCII " << std::endl;
        fname << "DATASET " << "UNSTRUCTURED_GRID " << std::endl;
        fname << "POINTS " << VTKGrid::Points.rows() << " double" << std::endl;
        for (int i = 0; i < VTKGrid::Points.rows(); i++)
        {
            for (int j = 0; j < VTKGrid::Points.cols(); j++)
            {
                fname << VTKGrid::Points(i, j) << " ";
                //std::cout << VTKGrid::Points(i, j) << std::endl;
            }
            fname << std::endl;
        }
        fname << "CELLS " << VTKGrid::Cells.rows() << " " << VTKGrid::Cells.rows() * (VTKGrid::Cells.cols() + 1) << std::endl;
        for (int i = 0; i < VTKGrid::Cells.rows(); i++)
        {
            fname << 4 << " ";
            for (int j = 0; j < VTKGrid::Cells.cols(); j++)
                fname << VTKGrid::Cells(i, j) << " ";
            fname << std::endl;
        }
        fname << "CELL_TYPES " << VTKGrid::Cells.rows() << std::endl;
        for (int i = 0; i < VTKGrid::Cells.rows(); i++)
            fname << VTKGrid::eletype << std::endl;

        fname << "POINT_DATA " << VTKGrid::Points.rows() << std::endl;
        fname << "VECTORS " << "U " << "double" << std::endl;
        for (int i = 0; i < VTKGrid::U.rows(); i++)
        {
            for (int j = 0; j < VTKGrid::U.cols(); j++)
                fname << VTKGrid::U(i, j) << " ";
            fname << std::endl;
        }
        fname << "VECTORS " << "theta " << "double" << std::endl;;
        for (int i = 0; i < VTKGrid::theta.rows(); i++)
        {
            for (int j = 0; j < VTKGrid::theta.cols(); j++)
                fname << VTKGrid::theta(i, j) << " ";
            fname << std::endl;
        }
        fname << "SCALARS " << "CP " << "double " << 1 << std::endl;
        fname << "LOOKUP_TABLE " << "default " << std::endl;
        for (int i = 0; i < VTKGrid::gap.size(); i++)
            fname << std::setprecision(7) << VTKGrid::gap(i) * epsilon << std::endl;
    }
    else
    {
        std::cout << "Unable to open file for saving current configuration " << std::endl;
    }

}

VTKGrid::VTKGrid(NLEBBE3D EBBE3D, int refine, double load, std::string BeamName, Eigen::MatrixXd GU, int ndim, int eletype, int ndof, int ncells, int beamnum)
{
    this->nnode = EBBE3D.get_nnode();
    this->nelem = EBBE3D.get_nelem();
    this->eletype = eletype;
    this->Points = Eigen::MatrixXd::Zero(this->nnode * refine, ndim);
    this->Cells = Eigen::MatrixXd::Zero(this->nelem * refine * 2, ncells);
    this->U = Eigen::MatrixXd::Zero(this->nnode * refine, ndim);
    this->gap = Eigen::VectorXd::Zero(this->nnode * refine);
    this->theta = Eigen::MatrixXd::Zero(this->nnode * refine, ndim);
    int temp1 = 0;
    for (int i = beamnum - 1; i >= 0; i--)
        temp1 += this->nnode * ndof;
    //Find nodes on the surface
    for (int i = 0; i < VTKGrid::nelem; i++)
    {
        int p1 = EBBE3D.get_connectivity(i, 1) - 1;
        int p2 = EBBE3D.get_connectivity(i, 2) - 1;
        int p3 = EBBE3D.get_connectivity(i, 3) - 1;
        //std::cout << p1 << std::endl;
        //std::cout << p2 << std::endl;
        //std::cout << p3 << std::endl;
        //normal vectors
        Eigen::VectorXd a(3), b(3), c(3);
        for (int k = 0; k < 3; k++)
        {
            a(k) = EBBE3D.get_coordinates(p1, k) + GU(ndof * p1 + k + temp1);
            b(k) = EBBE3D.get_coordinates(p2, k) + GU(ndof * p2 + k + temp1);
            c(k) = EBBE3D.get_coordinates(p3, k) + GU(ndof * p3 + k + temp1);
        }

        for (int k = 0; k < refine; k++)
            for (int j = 0; j < 3; j++)
            {
                VTKGrid::U(p1 * refine + k, j) = GU(ndof * p1 + j + temp1);
                VTKGrid::theta(p1 * refine + k, j) = GU(ndof * p1 + j + temp1 + 3);
                VTKGrid::U(p2 * refine + k, j) = GU(ndof * p2 + j + temp1);
                VTKGrid::theta(p2 * refine + k, j) = GU(ndof * p2 + j + temp1 + 3);
                VTKGrid::U(p3 * refine + k, j) = GU(ndof * p3 + j + temp1);
                VTKGrid::theta(p3 * refine + k, j) = GU(ndof * p3 + j + temp1 + 3);
            }

        //std::cout << a << std::endl;
        //std::cout << b << std::endl;
        //std::cout << c << std::endl;
        Eigen::VectorXd n1(3), n2(3);
        n1 = b - a;
        n2 = c - b;
        double r = EBBE3D.get_diameter() / 2;
        if (i != VTKGrid::nelem - 1)
        {
            Eigen::MatrixXd temp = VTKGrid::FindSurfacePoints(n1, a, r, refine);
            for (int k = 0; k < refine; k++)
                for (int j = 0; j < 3; j++)
                {
                    VTKGrid::Points(p1 * refine + k, j) = temp(k, j);
                    //std::cout << VTKGrid::Points(p1 * refine, j) << std::endl;
                }
            temp = VTKGrid::FindSurfacePoints(n2, b, r, refine);
            for (int k = 0; k < refine; k++)
                for (int j = 0; j < 3; j++)
                {
                    VTKGrid::Points(p2 * refine + k, j) = temp(k, j);
                    //std::cout << VTKGrid::Points(p1 * refine, j) << std::endl;
                }
        }
        else
        {
            Eigen::MatrixXd temp = VTKGrid::FindSurfacePoints(n1, a, r, refine);
            for (int k = 0; k < refine; k++)
                for (int j = 0; j < 3; j++)
                    VTKGrid::Points(p1 * refine + k, j) = temp(k, j);
            temp = VTKGrid::FindSurfacePoints(n2, b, r, refine);
            for (int k = 0; k < refine; k++)
                for (int j = 0; j < 3; j++)
                    VTKGrid::Points(p2 * refine + k, j) = temp(k, j);
            temp = VTKGrid::FindSurfacePoints(n2, c, r, refine);
            for (int k = 0; k < refine; k++)
                for (int j = 0; j < 3; j++)
                    VTKGrid::Points(p3 * refine + k, j) = temp(k, j);
        }
        Eigen::MatrixXd temp = VTKGrid::FindCells(p1, p2, refine);
        for (int k = 0; k < refine; k++)
        {
            for (int j = 0; j < 4; j++)
                VTKGrid::Cells(p1 * refine + k, j) = temp(k, j);
        }
        temp = VTKGrid::FindCells(p2, p3, refine);
        for (int k = 0; k < refine; k++)
        {
            for (int j = 0; j < 4; j++)
                VTKGrid::Cells(p2 * refine + k, j) = temp(k, j);
        }


    }

    VTKGrid::WriteVTK(load, BeamName, ".\\3DBeam\\", beamnum, 0);
}

VTKGrid::VTKGrid(NLEBBE3D EBBE3D, int refine, double load, std::string BeamName, Eigen::MatrixXd GU, int ndim, int eletype, int ndof, int ncells, int beamnum, double ** g, double epsilon)
{
    this->nnode = EBBE3D.get_nnode();
    this->nelem = EBBE3D.get_nelem();
    this->eletype = eletype;
    this->Points = Eigen::MatrixXd::Zero(this->nnode * refine, ndim);
    this->Cells = Eigen::MatrixXd::Zero(this->nelem * refine * 2, ncells);
    this->U = Eigen::MatrixXd::Zero(this->nnode * refine, ndim);
    this->gap = Eigen::VectorXd::Zero(this->nnode * refine);
    this->theta = Eigen::MatrixXd::Zero(this->nnode * refine, ndim);
    int temp1 = 0;
    for (int i = beamnum; i > 0; i--)
        temp1 += this->nnode * ndof;
    //Find nodes on the surface
    for (int i = 0; i < VTKGrid::nelem; i++)
    {
        int p1 = EBBE3D.get_connectivity(i, 1) - 1;
        int p2 = EBBE3D.get_connectivity(i, 2) - 1;
        int p3 = EBBE3D.get_connectivity(i, 3) - 1;
        //std::cout << p1 << std::endl;
        //std::cout << p2 << std::endl;
        //std::cout << p3 << std::endl;
        //normal vectors
        Eigen::VectorXd a(3), b(3), c(3);
        for (int k = 0; k < 3; k++)
        {
            a(k) = EBBE3D.get_coordinates(p1, k) + GU(ndof * p1 + k + temp1);
            b(k) = EBBE3D.get_coordinates(p2, k) + GU(ndof * p2 + k + temp1);
            c(k) = EBBE3D.get_coordinates(p3, k) + GU(ndof * p3 + k + temp1);
        }

        for(int k = 0; k < refine; k++)
            for (int j = 0; j < 3; j++)
            {
                VTKGrid::U(p1 * refine + k, j) = GU(ndof * p1 + j + temp1);
                VTKGrid::theta(p1 * refine + k, j) = GU(ndof * p1 + j + temp1 + 3);
                VTKGrid::U(p2 * refine + k, j) = GU(ndof * p2 + j + temp1);
                VTKGrid::theta(p2 * refine + k, j) = GU(ndof * p2 + j + temp1 + 3);
                VTKGrid::U(p3 * refine + k, j) = GU(ndof * p3 + j + temp1);
                VTKGrid::theta(p3 * refine + k, j) = GU(ndof * p3 + j + temp1 + 3);
            }

        for (int k = 0; k < refine; k++)
        {
            VTKGrid::gap(p1 * refine + k) = g[beamnum][p1];
            VTKGrid::gap(p2 * refine + k) = g[beamnum][p2];
            VTKGrid::gap(p3 * refine + k) = g[beamnum][p3];
        }

        //std::cout << a << std::endl;
        //std::cout << b << std::endl;
        //std::cout << c << std::endl;
        Eigen::VectorXd n1(3), n2(3);
        n1 = b - a;
        //std::cout << n1 << std::endl;
        n2 = c - b;
        //std::cout << n2 << std::endl;
        double r = EBBE3D.get_diameter() / 2.0;
        if (i != VTKGrid::nelem - 1)
        {
            Eigen::MatrixXd temp = VTKGrid::FindSurfacePoints(n1, a, r, refine);
            for (int k = 0; k < refine; k++)
                for (int j = 0; j < 3; j++)
                {
                    VTKGrid::Points(p1 * refine + k, j) = temp(k, j);
                    //std::cout << VTKGrid::Points(p1 * refine, j) << std::endl;
                }
            temp = VTKGrid::FindSurfacePoints(n2, b, r, refine);
            for (int k = 0; k < refine; k++)
                for (int j = 0; j < 3; j++)
                {
                    VTKGrid::Points(p2 * refine + k, j) = temp(k, j);
                    //std::cout << VTKGrid::Points(p1 * refine, j) << std::endl;
                }
        }
        else
        {
            Eigen::MatrixXd temp = VTKGrid::FindSurfacePoints(n1, a, r, refine);
            for (int k = 0; k < refine; k++)
                for (int j = 0; j < 3; j++)
                    VTKGrid::Points(p1 * refine + k, j) = temp(k, j);
            temp = VTKGrid::FindSurfacePoints(n2, b, r, refine);
            for (int k = 0; k < refine; k++)
                for (int j = 0; j < 3; j++)
                    VTKGrid::Points(p2 * refine + k, j) = temp(k, j);
            temp = VTKGrid::FindSurfacePoints(n2, c, r, refine);
            for (int k = 0; k < refine; k++)
                for (int j = 0; j < 3; j++)
                    VTKGrid::Points(p3 * refine + k, j) = temp(k, j);
        }
        Eigen::MatrixXd temp = VTKGrid::FindCells(p1, p2, refine);
        for (int k = 0; k < refine; k++)
        {
            for (int j = 0; j < 4; j++)
                VTKGrid::Cells(p1 * refine + k, j) = temp(k, j);
        }
        temp = VTKGrid::FindCells(p2, p3, refine);
        for (int k = 0; k < refine; k++)
        {
            for (int j = 0; j < 4; j++)
                VTKGrid::Cells(p2 * refine + k, j) = temp(k, j);
        }


    }

    VTKGrid::WriteVTK(load, BeamName, ".\\3DBeam\\", beamnum, epsilon);
}

VTKGrid::VTKGrid(NLEBBE2D EBBE2D, double load, std::string BeamName, Eigen::MatrixXd GU, int ndim, int eletype, int ndof, int ncells)
{
    VTKGrid::nnode = EBBE2D.get_nnode();
    VTKGrid::nelem = EBBE2D.get_nelem();
    VTKGrid::eletype = eletype;
    VTKGrid::Points = Eigen::MatrixXd::Zero(VTKGrid::nnode * 3, 3);
    VTKGrid::Cells = Eigen::MatrixXd::Zero(VTKGrid::nelem * 2, ncells);
    int count = 0;

    for (int i = 0; i < VTKGrid::nelem; i++)
    {
        int p1 = EBBE2D.get_connectivity(i, 1) - 1;
        int p2 = EBBE2D.get_connectivity(i, 2) - 1;

        Eigen::VectorXd a = Eigen::VectorXd::Zero(3);
        Eigen::VectorXd b = Eigen::VectorXd::Zero(3);
        for (int k = 0; k < 2; k++)
        {
            a(k) = EBBE2D.get_coordinates(p1, k) + GU(ndof * p1 + k);
            b(k) = EBBE2D.get_coordinates(p2, k) + GU(ndof * p2 + k);
        }
        double h = EBBE2D.get_modelprop("h");

        Eigen::VectorXd n1(3);
        n1 = b - a;
        Eigen::VectorXd n2(3);
        n2(0) = 0;
        n2(1) = 0;
        n2(2) = 1;
        Eigen::VectorXd n3 = CrossProduct(n1, n2);
        Eigen::VectorXd n3unit = n3 / n3.norm();
        //Points
        if (i != VTKGrid::nelem - 1)
        {
            for (int j = 0; j < ndim; j++)
            {
                VTKGrid::Points(p1 * 3 + 0, j) = a(j);
                VTKGrid::Points(p1 * 3 + 1, j) = a(j) + n3unit(j) * h;
                VTKGrid::Points(p1 * 3 + 2, j) = a(j) - n3unit(j) * h;
            }
        }
        else
        {
            for (int j = 0; j < ndim; j++)
            {
                VTKGrid::Points(p1 * 3 + 0, j) = a(j);
                VTKGrid::Points(p1 * 3 + 1, j) = a(j) + n3unit(j) * h;
                VTKGrid::Points(p1 * 3 + 2, j) = a(j) - n3unit(j) * h;
            }
            for (int j = 0; j < ndim; j++)
            {
                VTKGrid::Points(p2 * 3 + 0, j) = b(j);
                VTKGrid::Points(p2 * 3 + 1, j) = b(j) + n3unit(j) * h;
                VTKGrid::Points(p2 * 3 + 2, j) = b(j) - n3unit(j) * h;
            }
        }
        //Cell 1
        VTKGrid::Cells(count, 0) = p1 * 3;
        VTKGrid::Cells(count, 1) = p2 * 3;
        VTKGrid::Cells(count, 2) = p2 * 3 + 1;
        VTKGrid::Cells(count, 3) = p1 * 3 + 1;
        count++;
        //Cell 2
        VTKGrid::Cells(count, 0) = p1 * 3 + 2;
        VTKGrid::Cells(count, 1) = p2 * 3 + 2;
        VTKGrid::Cells(count, 2) = p2 * 3;
        VTKGrid::Cells(count, 3) = p1 * 3;
        count++;
    }

    VTKGrid::WriteVTK(load, BeamName, ".\\2DBeam\\", 1, 0);
}



Eigen::MatrixXd VTKGrid::FindSurfacePoints(Eigen::VectorXd n1, Eigen::VectorXd centre, double radius, int refine)
{
    //std::cout << radius << std::endl;
    //Second normal
    Eigen::VectorXd n2(3);
    n2(0) = -n1(1);
    n2(1) = n1(0);
    n2(2) = 0;
    if (n2(0) == 0 && n2(1) == 0 && n2(2) == 0)
    {
        n2(0) = 0;
        n2(1) = -n1(2);
        n2(2) = n1(1);
    }
    //std::cout << n2 << std::endl;
    //std::cout << n1 << std::endl;
    //normalize
    //std::cout << n2.norm() << std::endl;
    Eigen::VectorXd n2unit = n2 / n2.norm();
    //std::cout << n2unit << std::endl;
    //Evaluate cross product to find the third vector
    Eigen::VectorXd n3 = CrossProduct(n1, n2unit);
    //normalize
    Eigen::VectorXd n3unit = n3 / n3.norm();
    //std::cout << n3unit << std::endl;
    //locus of circle
    Eigen::MatrixXd Points = Eigen::MatrixXd::Zero(refine, 3);
    //std::cout << M_PI << std::endl;
    for (int i = 0; i < refine; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Points(i, j) = centre(j) + radius * cos((i * 2.0 * M_PI) / refine) * n2unit(j) + radius * sin((i * 2.0 * M_PI) / refine) * n3unit(j);
            //std::cout << Points(i, j) << std::endl;

        }
    }
    //std::cout << Points << std::endl;
    return Points;
}

Eigen::MatrixXd VTKGrid::FindCells(int p1, int p2, int refine)
{
    Eigen::MatrixXd Cells(refine, 4);
    int i;
    for (i = 0; i < refine - 1; i++)
    {
        Cells(i, 0) = p1 * refine + i;    
        Cells(i, 1) = p2 * refine + i;
        Cells(i, 2) = p2 * refine + i + 1;
        Cells(i, 3) = p1 * refine + i + 1;
    }
    Cells(refine - 1, 0) = p1 * refine + i;
    Cells(refine - 1, 1) = p2 * refine + i;
    Cells(refine - 1, 2) = p2 * refine;
    Cells(refine - 1, 3) = p1 * refine;



    return Cells;
}
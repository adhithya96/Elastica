#pragma once

#include "Variables.h"

//Hyerelastic validation
Loading::Loading()
{
	this->ELEM = Eigen::MatrixXd::Zero(1, 4);

	//node
	this->ELEM(0, 0) = 21;
	//type
	this->ELEM(0, 1) = -2;
	//direction
	this->ELEM(0, 2) = 1;
	//value
	this->ELEM(0, 3) = 100;
}

// -1 fixed boundary condition
// -2 hinged boundary condition
// 0 free boundary condition
Boundary::Boundary()
{
	this->ELEM = Eigen::MatrixXd::Zero(2, 2);
	this->ELEM(0, 0) = 1;
	this->ELEM(0, 1) = -1;
}


// -1 fixed boundary condition
// -2 hinged boundary condition
// 0 free boundary condition
/*Boundary::Boundary()
{
	ELEM = Eigen::MatrixXd::Zero(2, 2);
	ELEM(0, 0) = 1;
	ELEM(0, 1) = -1;
}


//Linear elastic validation
Loading::Loading()
{
	this->ELEM = Eigen::MatrixXd::Zero(1, 4);
	//node
	this->ELEM(0, 0) = 17;
	//type
	this->ELEM(0, 1) = -2;
	//direction
	this->ELEM(0, 2) = 3;
	//value
	this->ELEM(0, 3) = 833;
}*/

/*Boundary::Boundary(int beamnum, int contactswitch)
{
	if (contactswitch == 1)
	{
		if (beamnum == 1)
		{
			this->ELEM = Eigen::MatrixXd::Zero(2, 2);
			this->ELEM(0, 0) = 1;
			this->ELEM(0, 1) = -1;
			this->ELEM(1, 0) = 5;
			this->ELEM(1, 1) = -1;
		}
		else if (beamnum == 2)
		{
			//this->ELEM = Eigen::MatrixXd::Zero(0, 0);
			this->ELEM = Eigen::MatrixXd(2, 2);
			this->ELEM(0, 0) = 1;
			this->ELEM(0, 1) = -1;
			this->ELEM(1, 0) = 5;
			this->ELEM(1, 1) = -1;
		}
	}
	else
	{
		if (beamnum == 1)
		{
			this->ELEM = Eigen::MatrixXd::Zero(2, 2);
			this->ELEM(0, 0) = 1;
			this->ELEM(0, 1) = -1;
			this->ELEM(1, 0) = 5;
			this->ELEM(1, 1) = -1;
		}
		else if (beamnum == 2)
		{
			this->ELEM = Eigen::MatrixXd(2, 2);
			this->ELEM(0, 0) = 1;
			this->ELEM(0, 1) = -1;
			this->ELEM(1, 0) = 5;
			this->ELEM(1, 1) = -1;
		}
	}
}*/
//Loading GEBT Multibeam
void Loading::SetLoad(Eigen::VectorXd* GR, VAMBeamElement* VAMBE, int beamnum, int fiter, int coniter)
{
	for (int i = 0; i < this->ELEM.rows(); i++)
	{
		int lnode = this->ELEM(i, 0) - 1;
		int type = this->ELEM(i, 1);
		int dir = abs(this->ELEM(i, 2)) - 1;
		double value = this->ELEM(i, 3);
		int sign = this->ELEM(i, 2) / abs(this->ELEM(i, 2));

		//distributed load
		if (type == -1)
		{
			int iter = 0;
			for (int j = 0; j < beamnum; j++)
				iter += VAMBE[j].get_nnode() * VAMBE[j].get_ndof();
			for (int j = 0; j < VAMBE[beamnum].get_nnode(); j++)
				(*GR)(iter + VAMBE[beamnum].get_ndof() * j + dir) -= value * (fiter - 1) * sign / (coniter - 1);
		}
		//concentrated load
		else if (type == -2)
		{
			int iter = 0;
			for (int j = 0; j < beamnum; j++)
				iter += VAMBE[j].get_nnode() * VAMBE[j].get_ndof();

			iter += lnode * VAMBE[beamnum].get_ndof() + dir;
			//std::cout << iter << std::endl;

			(*GR)(iter) -= value * (fiter - 1) * sign / (coniter - 1);
			std::cout << value * (fiter - 1) * sign / (coniter - 1) << std::endl;
		}
	}
}






//Boundary condition
//Multibeam case
void Boundary::SetBoundary(Eigen::SparseMatrix<double, Eigen::ColMajor>* GT, Eigen::VectorXd* GR, NLEBBE3D* EBBE3D, int beamnum, int fiter, int coniter, Eigen::VectorXd* GU)
{
	for (int i = 0; i < this->ELEM.rows(); i++)
	{
		int bnode = this->ELEM(i, 0) - 1;
		int type = this->ELEM(i, 1);
		//fixed boundary condition
		if (type == -1)
		{
			int iter1 = 0, iter2 = 0;
			for (int j = 0; j < beamnum; j++)
				iter1 += EBBE3D[j].get_nnode() * EBBE3D[j].get_ndof();
			iter2 = iter1 + bnode * EBBE3D[beamnum].get_ndof();
			for (int j = iter2; j < iter2 + EBBE3D[beamnum].get_ndof(); j++)
			{
				for (int k = iter1; k < iter1 + EBBE3D[beamnum].get_nnode() * EBBE3D[beamnum].get_ndof(); k++)
					GT->coeffRef(j, k) = 0;
				GT->coeffRef(j, j) = 1;
				(*GR)(j) = 0;
			}
		}
		else if (type == -2)
		{
			int iter1 = 0, iter2 = 0;
			for (int j = 0; j < beamnum; j++)
				iter1 += EBBE3D[j].get_nnode() * EBBE3D[j].get_ndof();
			iter2 = iter1 + bnode * EBBE3D[beamnum].get_ndof();
			for (int j = iter2; j < iter2 + 3; j++)
			{
				for (int k = iter1; k < iter1 + EBBE3D[beamnum].get_nnode() * EBBE3D[beamnum].get_ndof(); k++)
					GT->coeffRef(j, k) = 0;
				GT->coeffRef(j, j) = 1;
				(*GR)(j) = 0;
			}
		}
		//Displacement boundary condition
		else if (type == -3)
		{
			int sign = ELEM(i, 2) / abs(ELEM(i, 2));
			int dir = abs(ELEM(i, 2)) - 1 ;
			double value = ELEM(i, 3);
			int iter1 = 0, iter2 = 0;
			for (int j = 0; j < beamnum; j++)
				iter1 += EBBE3D[j].get_nnode() * EBBE3D[j].get_ndof();
			iter2 = iter1 + bnode * EBBE3D[beamnum].get_ndof() + dir;
				
			for (int k = iter1; k < iter1 + EBBE3D[beamnum].get_nnode() * EBBE3D[beamnum].get_ndof(); k++)
				GT->coeffRef(iter2, k) = 0;
			GT->coeffRef(iter2, iter2) = 1;
			(*GR)(iter2) = 0;
			(*GU)(iter2) = value * sign * (fiter  - 1) * 1.0 / (coniter - 1);
			std::cout << value * sign * (fiter - 1) * 1.0 / (coniter - 1) << std::endl;
		}
		else if (type = -4)
		{
			int iter1 = 0;
			for (int j = 0; j < beamnum; j++)
				iter1 += EBBE3D[j].get_nnode() * EBBE3D[j].get_ndof();
			for (int j = iter1; j < iter1 + EBBE3D[beamnum].get_nnode() * EBBE3D[beamnum].get_ndof(); j++)
			{
				for (int k = iter1; k < iter1 + EBBE3D[beamnum].get_nnode() * EBBE3D[beamnum].get_ndof(); k++)
					GT->coeffRef(j, k) = 0;
				GT->coeffRef(j, j) = 1;
				(*GR)(j) = 0;
			}
		}
	}
	//std::cout << ELEM.rows() << std::endl;	
}

//Multi beam case
//Loading condition
void Loading::SetLoad(Eigen::VectorXd* GR, NLEBBE3D* EBBE3D, int beamnum, int fiter, int coniter)
{
	for (int i = 0; i < this->ELEM.rows(); i++)
	{
		int lnode = this->ELEM(i, 0) - 1;
		int type = this->ELEM(i, 1);
		int dir = abs(this->ELEM(i, 2)) - 1;
		double value = this->ELEM(i, 3);
		int sign = this->ELEM(i, 2) / abs(this->ELEM(i, 2));

		//distributed load
		if (type == -1)
		{
			int iter = 0;
			for (int j = 0; j < beamnum; j++)
				iter += EBBE3D[j].get_nnode() * EBBE3D[j].get_ndof();
			for (int j = 0; j < EBBE3D[beamnum].get_nnode(); j++)
				(*GR)(iter + EBBE3D[beamnum].get_ndof() * j + dir) -= value * (fiter - 1) * sign / (coniter - 1);

		}
		//concentrated load
		else if (type == -2)
		{
			int iter = 0;
			for (int j = 0; j < beamnum; j++)
				iter += EBBE3D[j].get_nnode() * EBBE3D[j].get_ndof();
			iter += lnode * EBBE3D[beamnum].get_ndof() + dir;
			(*GR)(iter) -= value * (fiter - 1) * sign / (coniter - 1);
			std::cout << value * (fiter - 1) * sign / (coniter - 1) << std::endl;
		}
	}
}

//Single beam case
//Loading condition
void Loading::SetLoad(Eigen::VectorXd* GR, NLEBBE3D* EBBE3D, int fiter, int coniter)
{
	for (int i = 0; i < this->ELEM.rows(); i++)
	{
		int lnode = this->ELEM(i, 0) - 1;
		int type = this->ELEM(i, 1);
		int dir = abs(this->ELEM(i, 2)) - 1;
		double value = this->ELEM(i, 3);
		int sign = this->ELEM(i, 2) / abs(this->ELEM(i, 2));
		//std::cout << type << std::endl;
		//std::cout << dir << std::endl;
		//std::cout << value << std::endl;
		//distributed load
		if (type == -1)
		{
			for (int j = 0; j < (*EBBE3D).get_nnode(); j++)
				(*GR)((*EBBE3D).get_ndof() * j + abs(dir)) -= (double)value * (fiter - 1) * sign / (coniter - 2);

		}
		//concentrated load
		else if (type == -2)
		{
			(*GR)(lnode * (*EBBE3D).get_ndof() + abs(dir)) -= (double)value * (fiter - 1) * sign / (coniter - 2);
		}
	}
}

//Boundary condition
//Single beam case
void Boundary::SetBoundary(Eigen::SparseMatrix<double, Eigen::ColMajor>* GT, Eigen::VectorXd* GR, NLEBBE3D* EBBE3D, int fiter, int coniter, Eigen::VectorXd* GU)
{
	for (int i = 0; i < this->ELEM.rows(); i++)
	{
		int bnode = this->ELEM(i, 0) - 1;
		int type = this->ELEM(i, 1);
		//fixed boundary condition
		if (type == -1)
		{
			for (int j = bnode * (*EBBE3D).get_ndof(); j < (bnode + 1) * (*EBBE3D).get_ndof(); j++)
			{
				for (int k = 0; k < (*EBBE3D).get_nnode() * (*EBBE3D).get_ndof(); k++)
					GT->coeffRef(j, k) = 0;
				GT->coeffRef(j, j) = 1;
				(*GR)(j) = 0;
			}
		}
		else if (type == -2)
		{
			for (int j = bnode * (*EBBE3D).get_ndof(); j < bnode * (*EBBE3D).get_ndof() + 3; j++)
			{
				for (int k = 0; k < (*EBBE3D).get_nnode() * (*EBBE3D).get_ndof(); k++)
					GT->coeffRef(j, k) = 0;
				GT->coeffRef(j, j) = 1;
				(*GR)(j) = 0;
			}
		}
		//displacement boundary condition
		else if (type == -3)
		{
			int dir = ELEM(i, 2);
			int value = ELEM(i, 3);
			int j = bnode * (*EBBE3D).get_ndof() + dir;
			for (int k = 0; k < (*EBBE3D).get_nnode() * (*EBBE3D).get_ndof(); k++)
				GT->coeffRef(j, k) = 0;
			GT->coeffRef(j, j) = 1;

			(*GR)(j) = 0;
			(*GU)(j) = (double)value * (fiter - 1) / (coniter - 2);
			//std::cout << (double)value * (fiter - 1) / (coniter - 2) << std::endl;
		}
	}
}
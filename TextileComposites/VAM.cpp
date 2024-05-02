
#ifndef VARIABLES_H
#define VARIABLES_H


#include "Variables.h"


int VAMBeamElement::get_nls()
{
    return this->NLS;
}

int VAMBeamElement::get_ndof()
{
    return this->NDOF;
}

int VAMBeamElement::get_nnode()
{
    return this->NNODE;
}

int VAMBeamElement::get_nelem()
{
    return this->NELEM;
}

double VAMBeamElement::get_coordinates(int i, int j)
{
    return this->NODE(i, j);
}

int VAMBeamElement::get_connectivity(int i, int j)
{
    return this->ELEM(i, j);
}

Eigen::VectorXd VAMBeamElement::get_inittwist()
{
    return this->inittwist;
}

double VAMBeamElement::get_width()
{
    return b;
}

Eigen::VectorXd VAMBeamElement::get_load(int choice)
{
    if (choice == 1)
        return this->a1;
    else if (choice == 2)
        return this->b1;
    else if (choice == 3)
        return this->aN;
    else if (choice == 4)
        return this->bN;
}

void VAMBeamElement::set_load(int choice, int dir, double value)
{
    int temp = dir - 1;
    if (choice == 1)
        this->a1(temp) = value;
    else if (choice == 2)
        this->b1(temp) = value;
    else if (choice == 3)
    {
        this->aN(temp) = value;
        std::cout << value << std::endl;
    }
    else if (choice == 4)
        this->bN(temp) = value;

    
}

#endif
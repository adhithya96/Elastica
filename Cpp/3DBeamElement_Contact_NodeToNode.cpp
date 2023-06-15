#ifndef VARIABLES_H
#define VARIABLES_H
#include "Variables.h"
#include<iomanip>

/*************************************************************
* AceGen    7.505 Windows (16 Aug 22)                        *
*           Co. J. Korelc  2020           19 May 23 14:24:13 *
**************************************************************
User     : Full professional version
Notebook : ContactElement
Evaluation time                 : 1 s     Mode  : Debug
Number of formulae              : 52      Method: Automatic
Subroutine                      : Beam Contact size: 554
Total size of Mathematica  code : 554 subexpressions
Total size of C code            : 3181 bytes */
#include "sms.h"

/******************* S U B R O U T I N E *********************/
//D[1] - penalty parameter
//{D[2],D[3],D[4]} - normal vector
//X1 - position vector of master node
//X2 - position vector of slave node
//u1 - displacement vector of master node
//u2 - displacement vector of slave node
//CR - residual vector for contact
//CT - tangent matrix for contact
void Beam_Contact(double D[4], double X1[3], double X2[3], double u1[6], double u2[6], double* CR, double** CT)
{
    double v[251];

    int i38, i46, b34, b35;
    /* 1 = X01_1 */
    v[1] = X1[0];
    /* 2 = X01_2 */
    v[2] = X1[1];
    /* 3 = X01_3 */
    v[3] = X1[2];
    /* 4 = X02_1 */
    v[4] = X2[0];
    /* 5 = X02_2 */
    v[5] = X2[1];
    /* 6 = X02_3 */
    v[6] = X2[2];
    /* 7 = [u01_1][u_1] */
    v[7] = u1[0];
    /* 8 = [u01_2][u_2] */
    v[8] = u1[1];
    /* 9 = [u01_3][u_3] */
    v[9] = u1[2];
    /* 10 = [u01_4][u_4] */
    v[10] = u1[3];
    /* 11 = [u01_5][u_5] */
    v[11] = u1[4];
    /* 12 = [u01_6][u_6] */
    v[12] = u1[5];
    /* 13 = [u02_1][u_7] */
    v[13] = u2[0];
    /* 14 = [u02_2][u_8] */
    v[14] = u2[1];
    /* 15 = [u02_3][u_9] */
    v[15] = u2[2];
    /* 16 = [u02_4][u_10] */
    v[16] = u2[3];
    /* 17 = [u02_5][u_11] */
    v[17] = u2[4];
    /* 18 = [u02_6][u_12] */
    v[18] = u2[5];
    v[73] = v[7];
    v[74] = v[8];
    v[75] = v[9];
    v[76] = v[10];
    v[77] = v[11];
    v[78] = v[12];
    v[79] = v[13];
    v[80] = v[14];
    v[81] = v[15];
    v[82] = v[16];
    v[83] = v[17];
    v[84] = v[18];
    v[47] = 0e0;/*debug*/
    /* 19 = Data_1 */
    v[19] = D[0];
    /* 20 = Data_2 */
    v[20] = D[1];
    /* 21 = Data_3 */
    v[21] = D[2];
    /* 22 = Data_4 */
    v[22] = D[3];
    /* 23 = \[Epsilon] */
    v[23] = v[19];
    /* 24 = n_1 */
    v[24] = v[20];
    /* 25 = n_2 */
    v[25] = v[21];
    /* 26 = n_3 */
    v[26] = v[22];
    /* 27 = x1_1 */
    v[27] = v[1] + v[7];
    /* 28 = x1_2 */
    v[28] = v[2] + v[8];
    /* 29 = x1_3 */
    v[29] = v[3] + v[9];
    /* 30 = x2_1 */
    v[30] = v[13] + v[4];
    /* 31 = x2_2 */
    v[31] = v[14] + v[5];
    /* 32 = x2_3 */
    v[32] = v[15] + v[6];
    /* 33 = gN */
    v[33] = v[24] * (v[27] - v[30]) + v[25] * (v[28] - v[31]) + v[26] * (v[29] - v[32]);
    v[42] = 1e0 * v[23] * v[26] * v[33];
    v[41] = 1e0 * v[23] * v[25] * v[33];
    v[40] = 1e0 * v[23] * v[24] * v[33];
    v[85] = v[40];
    v[86] = v[41];
    v[87] = v[42];
    v[88] = 0e0;
    v[89] = 0e0;
    v[90] = 0e0;
    v[91] = -v[40];
    v[92] = -v[41];
    v[93] = -v[42];
    v[94] = 0e0;
    v[95] = 0e0;
    v[96] = 0e0;
    v[43] = 0e0;/*debug*/
    b34 = v[33] > 0e0;
    v[34] = b34;
    b35 = b34;
    v[35] = b35;
    if (b35) {
        return;
        v[36] = 0e0;/*debug*/
    }
    else {
    };
    /* 37 = \[CapitalPi] */
    v[37] = 0.5e0 * v[23] * (v[33] * v[33]);
    for (i38 = 1; i38 <= 12; i38++) {
        v[38] = i38;
        /* 39 = \[DoubleStruckCapitalG]_i */
        v[39] = v[72 + i38];
        /* 44 = Rgi */
        v[44] = v[84 + i38];
        /* 49 = \[OverBracket]_gN_(\[Yen]|Rgi)_42 */
        v[104] = 1e0;
        v[105] = 0e0;
        v[106] = 0e0;
        v[107] = 0e0;
        v[108] = 0e0;
        v[109] = 0e0;
        v[110] = -1e0;
        v[111] = 0e0;
        v[112] = 0e0;
        v[113] = 0e0;
        v[114] = 0e0;
        v[115] = 0e0;
        v[116] = 0e0;
        v[117] = 1e0;
        v[118] = 0e0;
        v[119] = 0e0;
        v[120] = 0e0;
        v[121] = 0e0;
        v[122] = 0e0;
        v[123] = -1e0;
        v[124] = 0e0;
        v[125] = 0e0;
        v[126] = 0e0;
        v[127] = 0e0;
        v[128] = 0e0;
        v[129] = 0e0;
        v[130] = 1e0;
        v[131] = 0e0;
        v[132] = 0e0;
        v[133] = 0e0;
        v[134] = 0e0;
        v[135] = 0e0;
        v[136] = -1e0;
        v[137] = 0e0;
        v[138] = 0e0;
        v[139] = 0e0;
        v[49] = 1e0 * v[103 + i38] * v[23] * v[24] + 1e0 * v[115 + i38] * v[23] * v[25] + 1e0 * v[127 + i38] * v[23] * v[26];
        v[52] = v[26] * v[49];
        v[51] = v[25] * v[49];
        v[50] = v[24] * v[49];
        v[140] = v[50];
        v[141] = v[51];
        v[142] = v[52];
        v[143] = 0e0;
        v[144] = 0e0;
        v[145] = 0e0;
        v[146] = -v[50];
        v[147] = -v[51];
        v[148] = -v[52];
        v[149] = 0e0;
        v[150] = 0e0;
        v[151] = 0e0;
        v[53] = 0e0;/*debug*/
        CR[i38 - 1] += v[44];
        v[45] = 0e0;/*debug*/
        for (i46 = i38; i46 <= 12; i46++) {
            v[46] = i46;
            /* 48 = \[DoubleStruckCapitalG]_j */
            v[48] = v[72 + i46];
            /* 54 = Kgij */
            v[54] = v[139 + i46];
            CT[i38 - 1][i46 - 1] += v[54];
            v[55] = 0e0;/*debug*/
        };/* end for */
    };/* end for */

    for (int i = 0; i < 12; i++)
    {
        for (int j = 0; j < 12; j++)
            if (j != i)
                CT[j][i] = CT[i][j];
    }
};

void ContactSearch(NonLinearEulerBernouliBeamElement3D EBBE3D1, NonLinearEulerBernouliBeamElement3D EBBE3D2, int** ContactPairs)
{
    for (int i = 0; i < EBBE3D2.NNODE; i++)
    {
        ContactPairs[i][0] = i + 1;
        ContactPairs[i][1] = i + 1;
    }
}

void PostProcessing(Eigen::MatrixXd X, Eigen::VectorXd U, double load, std::string BeamElement, int nnode, int ndm, int ndof)
{
    std::fstream file1;
    std::string str;
    std::ostringstream str1;
    str1 << load;
    str = "E:/Adhithya/MTech_Res_Thesis/Cpp/ThesisCode/TextileComposites/TextileComposites/3DBeam/" + BeamElement + "_" + str1.str() + ".txt";
    file1.open(str.c_str(), std::ios::out);

    std::cout << X.size() << std::endl;
    std::cout << str << std::endl;

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

#endif

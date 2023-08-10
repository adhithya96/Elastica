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


double BeamContact::get_penaltyparameter()
{
    return BeamContact::epsilon;
}

int BeamContact::get_nen()
{
    return BeamContact::NEN;
}

Eigen::MatrixXd BeamContact::get_Globalelem()
{
    return BeamContact::GlobalELEM;
}

/*void BeamContact::LocalContactSearch(NLEBBE3D* EBBE3D)
{
    if (choice == "NTN")
    {
        for (int i = 0; i < EBBE3D2.NNODE; i++)
        {
            ContactPairs[i][0] = i + 1;
            ContactPairs[i][1] = i + 1;
        }
    }
    else if (choice == "STS")
    {
        for (int i = 0; i < EBBE3D2.NELEM; i++)
        {
            ContactPairs[i][0] = i + 1;
            ContactPairs[i][1] = i + 1;
        }
    }
    else
    {
        std::cout << "Wrong choice" << std::endl;
        std::cout << "Choices for contact FE discretisation" << std::endl << "1. NTN  -  Node to Node" << std::endl << "2. STS  -  Segment to Segment" << std::endl;
    }
}*/


/******************* S U B R O U T I N E *********************/
//D[1] - penalty parameter
//{D[2],D[3],D[4]} - normal vector
//X1 - position vector of master node
//X2 - position vector of slave node
//u1 - displacement vector of master node
//u2 - displacement vector of slave node
//CR - residual vector for contact
//CT - tangent matrix for contact
//Tangent and residual for node to node contact
void BeamContact::Contact_NTN(double D[4], double X1[3], double X2[3], double u1[6], double u2[6], double* CR, double** CT, double* g)
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
    *g = v[33];
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
        for (i46 = i38; i46 <= 12; i46++)
        {
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

    for (int i = 0; i < 12; i++)
        std::cout << CR[i] << std::endl;
    for (int i = 0; i < 12; i++)
    {
        for (int j = 0; j < 12; j++)
            std::cout << CT[i][j] << "  ";
        std::cout << std::endl;
    }
}


//Tangent and residual for segment to segment contact
/******************* S U B R O U T I N E *********************/
void Contact_STS(double D[4], double X1[3], double X2[3]
    , double X3[3], double X4[3], double u1[6], double u2[6], double u3[6], double u4[6],
    Eigen::SparseMatrix<double, Eigen::ColMajor>* GT, Eigen::VectorXd* GR, int mnode1, int mnode2, 
    int snode1, int snode2)
{
    double RP[2];

    double exi = D[1];
    double eta = D[2];
    //master current config 
    double x1[3], x2[3];
    //slave current config
    double y1[3], y2[3];

    for (int j = 0; j < 3; j++)
    {
        x1[j] = X1[j] + u1[j];
        x2[j] = X2[j] + u2[j];

        y1[j] = X3[j] + u3[j];
        y2[j] = X4[j] + u4[j];
    }
    Eigen::VectorXd bm(3), tm(3), ts(3), bs(3);
    for (int j = 0; j < 3; j++)
    {
        bm(j) = x2[j] + x1[j];
        bs(j) = y2[j] + y1[j];
        tm(j) = x2[j] - x1[j];
        ts(j) = y2[j] - y1[j];
    }
    Eigen::VectorXd exp(3), exp2(3);
    exp2 = (tm * (tm.dot(ts)) - ts * (tm.dot(tm))) / (tm.dot(tm) * ts.dot(ts) - tm.dot(ts) * tm.dot(ts));
    exp = (ts * (tm.dot(ts)) - tm * (ts.dot(ts))) / (tm.dot(tm) * ts.dot(ts) - tm.dot(ts) * tm.dot(ts));

    RP[0] = (bm - bs).dot(exp);
    RP[1] = -(bm - bs).dot(exp2);

    if ((RP[0]<1 && RP[1]>-1) && (RP[1]<1 && RP[1]>-1))
    {
        double CR[24], CT[24][24];
        //Initialize contact residual and contact tangent matrix
        for (int j = 0; j < 24; j++)
        {
            CR[j] = 0;
            for (int k = 0; k < 24; k++)
                CT[j][k] = 0;
        }

        double v[406];
        int i64, i73, b61, b62;
        /* 1 = X1_1 */
        v[1] = X1[0];
        /* 2 = X1_2 */
        v[2] = X1[1];
        /* 3 = X1_3 */
        v[3] = X1[2];
        /* 4 = X2_1 */
        v[4] = X2[0];
        /* 5 = X2_2 */
        v[5] = X2[1];
        /* 6 = X2_3 */
        v[6] = X2[2];
        /* 7 = X3_1 */
        v[7] = X3[0];
        /* 8 = X3_2 */
        v[8] = X3[1];
        /* 9 = X3_3 */
        v[9] = X3[2];
        /* 10 = X4_1 */
        v[10] = X4[0];
        /* 11 = X4_2 */
        v[11] = X4[1];
        /* 12 = X4_3 */
        v[12] = X4[2];
        /* 13 = [u1_1][u_1] */
        v[13] = u1[0];
        /* 14 = [u1_2][u_2] */
        v[14] = u1[1];
        /* 15 = [u1_3][u_3] */
        v[15] = u1[2];
        /* 16 = [u1_4][u_4] */
        v[16] = u1[3];
        /* 17 = [u1_5][u_5] */
        v[17] = u1[4];
        /* 18 = [u1_6][u_6] */
        v[18] = u1[5];
        /* 19 = [u2_1][u_7] */
        v[19] = u2[0];
        /* 20 = [u2_2][u_8] */
        v[20] = u2[1];
        /* 21 = [u2_3][u_9] */
        v[21] = u2[2];
        /* 22 = [u2_4][u_10] */
        v[22] = u2[3];
        /* 23 = [u2_5][u_11] */
        v[23] = u2[4];
        /* 24 = [u2_6][u_12] */
        v[24] = u2[5];
        /* 25 = [u3_1][u_13] */
        v[25] = u3[0];
        /* 26 = [u3_2][u_14] */
        v[26] = u3[1];
        /* 27 = [u3_3][u_15] */
        v[27] = u3[2];
        /* 28 = [u3_4][u_16] */
        v[28] = u3[3];
        /* 29 = [u3_5][u_17] */
        v[29] = u3[4];
        /* 30 = [u3_6][u_18] */
        v[30] = u3[5];
        /* 31 = [u4_1][u_19] */
        v[31] = u4[0];
        /* 32 = [u4_2][u_20] */
        v[32] = u4[1];
        /* 33 = [u4_3][u_21] */
        v[33] = u4[2];
        /* 34 = [u4_4][u_22] */
        v[34] = u4[3];
        /* 35 = [u4_5][u_23] */
        v[35] = u4[4];
        /* 36 = [u4_6][u_24] */
        v[36] = u4[5];
        v[108] = v[13];
        v[109] = v[14];
        v[110] = v[15];
        v[111] = v[16];
        v[112] = v[17];
        v[113] = v[18];
        v[114] = v[19];
        v[115] = v[20];
        v[116] = v[21];
        v[117] = v[22];
        v[118] = v[23];
        v[119] = v[24];
        v[120] = v[25];
        v[121] = v[26];
        v[122] = v[27];
        v[123] = v[28];
        v[124] = v[29];
        v[125] = v[30];
        v[126] = v[31];
        v[127] = v[32];
        v[128] = v[33];
        v[129] = v[34];
        v[130] = v[35];
        v[131] = v[36];
        v[74] = 0e0;/*debug*/
        /* 37 = Data_1 */
        v[37] = D[0];
        /* 38 = [\[Xi]][Data_2] */
        v[38] = RP[0];
        /* 39 = [\[Eta]][Data_3] */
        v[39] = RP[1];
        /* 40 = [ds][Data_4] */
        v[40] = D[1];
        /* 41 = [dm][Data_5] */
        v[41] = D[2];
        /* 42 = \[Epsilon] */
        v[42] = v[37];
        /* 43 = Nm_1 */
        v[43] = (1e0 - v[38]) / 2e0;
        /* 44 = Nm_2 */
        v[44] = (1e0 + v[38]) / 2e0;
        /* 45 = Ns_1 */
        v[45] = (1e0 - v[39]) / 2e0;
        /* 46 = Ns_2 */
        v[46] = (1e0 + v[39]) / 2e0;
        /* 47 = Xm_1 */
        v[47] = v[1] * v[43] + v[4] * v[44];
        /* 48 = Xm_2 */
        v[48] = v[2] * v[43] + v[44] * v[5];
        /* 49 = Xm_3 */
        v[49] = v[3] * v[43] + v[44] * v[6];
        /* 50 = Xs_1 */
        v[50] = v[10] * v[46] + v[45] * v[7];
        /* 51 = Xs_2 */
        v[51] = v[11] * v[46] + v[45] * v[8];
        /* 52 = Xs_3 */
        v[52] = v[12] * v[46] + v[45] * v[9];
        /* 53 = xm_1 */
        v[53] = v[13] * v[43] + v[19] * v[44] + v[47];
        //std::cout << v[53] << std::endl;
        /* 54 = xm_2 */
        v[54] = v[14] * v[43] + v[20] * v[44] + v[48];
        //std::cout << v[54] << std::endl;
        /* 55 = xm_3 */
        v[55] = v[15] * v[43] + v[21] * v[44] + v[49];
        //std::cout << v[55] << std::endl;
        /* 56 = xs_1 */
        v[56] = v[25] * v[45] + v[31] * v[46] + v[50];
        v[67] = -v[53] + v[56];
        //std::cout << v[56] << std::endl;
        /* 57 = xs_2 */
        v[57] = v[26] * v[45] + v[32] * v[46] + v[51];
        v[68] = -v[54] + v[57];
        //std::cout << v[57] << std::endl;
        /* 58 = xs_3 */
        v[58] = v[27] * v[45] + v[33] * v[46] + v[52];
        v[69] = -v[55] + v[58];
        //std::cout << v[58] << std::endl;
        v[66] = sqrt((v[67] * v[67]) + (v[68] * v[68]) + (v[69] * v[69]));
        v[86] = 1e0 / v[66];
        v[84] = 1e0 / v[66];
        v[78] = 1e0 / v[66];
        v[76] = 1e0 / (v[66] * v[66]);
        /* 59 = gN */
        v[59] = -v[40] - v[41] + v[66];
        double g = v[59];
        //std::cout << v[59] << std::endl;
        v[83] = (1e0 * v[42] * v[46] * v[59]) / v[66];
        v[82] = (1e0 * v[42] * v[45] * v[59]) / v[66];
        v[81] = (-1e0 * v[42] * v[44] * v[59]) / v[66];
        v[80] = (-1e0 * v[42] * v[43] * v[59]) / v[66];
        v[132] = (-1e0 * v[42] * v[43] * v[59] * v[67]) / v[66];
        v[133] = (-1e0 * v[42] * v[43] * v[59] * v[68]) / v[66];
        v[134] = (-1e0 * v[42] * v[43] * v[59] * v[69]) / v[66];
        v[135] = 0e0;
        v[136] = 0e0;
        v[137] = 0e0;
        v[138] = (-1e0 * v[42] * v[44] * v[59] * v[67]) / v[66];
        v[139] = (-1e0 * v[42] * v[44] * v[59] * v[68]) / v[66];
        v[140] = (-1e0 * v[42] * v[44] * v[59] * v[69]) / v[66];
        v[141] = 0e0;
        v[142] = 0e0;
        v[143] = 0e0;
        v[144] = (1e0 * v[42] * v[45] * v[59] * v[67]) / v[66];
        v[145] = (1e0 * v[42] * v[45] * v[59] * v[68]) / v[66];
        v[146] = (1e0 * v[42] * v[45] * v[59] * v[69]) / v[66];
        v[147] = 0e0;
        v[148] = 0e0;
        v[149] = 0e0;
        v[150] = (1e0 * v[42] * v[46] * v[59] * v[67]) / v[66];
        v[151] = (1e0 * v[42] * v[46] * v[59] * v[68]) / v[66];
        v[152] = (1e0 * v[42] * v[46] * v[59] * v[69]) / v[66];
        v[153] = 0e0;
        v[154] = 0e0;
        v[155] = 0e0;
        v[70] = 0e0;/*debug*/
        /* 60 = \[CapitalPi] */
        v[60] = 0.5e0 * v[42] * (v[59] * v[59]);
        b61 = v[59] > 0e0;
        //std::cout << b61 << std::endl;
        v[61] = b61;
        b62 = b61;
        v[62] = b62;
        if (b62) {
            return;
            v[63] = 0e0;/*debug*/
        }
        else {
        };
        for (i64 = 1; i64 <= 24; i64++) {
            v[64] = i64;
            /* 65 = \[DoubleStruckCapitalG]_i */
            v[65] = v[107 + i64];
            /* 71 = Rgi */
            v[71] = v[131 + i64];
            /* 77 = \[OverBracket]_\[Yen]_66(gN|Rgi) */
            v[163] = (-1e0 * v[42] * v[43] * v[67]) / v[66];
            v[164] = (-1e0 * v[42] * v[43] * v[68]) / v[66];
            v[165] = (-1e0 * v[42] * v[43] * v[69]) / v[66];
            v[166] = 0e0;
            v[167] = 0e0;
            v[168] = 0e0;
            v[169] = (-1e0 * v[42] * v[44] * v[67]) / v[66];
            v[170] = (-1e0 * v[42] * v[44] * v[68]) / v[66];
            v[171] = (-1e0 * v[42] * v[44] * v[69]) / v[66];
            v[172] = 0e0;
            v[173] = 0e0;
            v[174] = 0e0;
            v[175] = (1e0 * v[42] * v[45] * v[67]) / v[66];
            v[176] = (1e0 * v[42] * v[45] * v[68]) / v[66];
            v[177] = (1e0 * v[42] * v[45] * v[69]) / v[66];
            v[178] = 0e0;
            v[179] = 0e0;
            v[180] = 0e0;
            v[181] = 1e0 * v[42] * v[46] * v[67] * v[86];
            v[182] = 1e0 * v[42] * v[46] * v[68] * v[84];
            v[183] = 1e0 * v[42] * v[46] * v[69] * v[78];
            v[184] = 0e0;
            v[185] = 0e0;
            v[186] = 0e0;
            v[187] = 1e0 * v[42] * v[43] * v[59] * v[67] * v[76];
            v[188] = 1e0 * v[42] * v[43] * v[59] * v[68] * v[76];
            v[189] = 1e0 * v[42] * v[43] * v[59] * v[69] * v[76];
            v[190] = 0e0;
            v[191] = 0e0;
            v[192] = 0e0;
            v[193] = 1e0 * v[42] * v[44] * v[59] * v[67] * v[76];
            v[194] = 1e0 * v[42] * v[44] * v[59] * v[68] * v[76];
            v[195] = 1e0 * v[42] * v[44] * v[59] * v[69] * v[76];
            v[196] = 0e0;
            v[197] = 0e0;
            v[198] = 0e0;
            v[199] = -1e0 * v[42] * v[45] * v[59] * v[67] * v[76];
            v[200] = -1e0 * v[42] * v[45] * v[59] * v[68] * v[76];
            v[201] = -1e0 * v[42] * v[45] * v[59] * v[69] * v[76];
            v[202] = 0e0;
            v[203] = 0e0;
            v[204] = 0e0;
            v[205] = -1e0 * v[42] * v[46] * v[59] * v[67] * v[76];
            v[206] = -1e0 * v[42] * v[46] * v[59] * v[68] * v[76];
            v[207] = -1e0 * v[42] * v[46] * v[59] * v[69] * v[76];
            v[208] = 0e0;
            v[209] = 0e0;
            v[210] = 0e0;
            v[77] = v[162 + i64] + v[186 + i64];
            /* 79 = \[OverBracket]_\[Yen]_69(\[Yen]|Rgi)_66 */
            v[211] = 0e0;
            v[212] = 0e0;
            v[213] = v[80];
            v[214] = 0e0;
            v[215] = 0e0;
            v[216] = 0e0;
            v[217] = 0e0;
            v[218] = 0e0;
            v[219] = v[81];
            v[220] = 0e0;
            v[221] = 0e0;
            v[222] = 0e0;
            v[223] = 0e0;
            v[224] = 0e0;
            v[225] = v[82];
            v[226] = 0e0;
            v[227] = 0e0;
            v[228] = 0e0;
            v[229] = 0e0;
            v[230] = 0e0;
            v[231] = v[83];
            v[232] = 0e0;
            v[233] = 0e0;
            v[234] = 0e0;
            v[79] = v[210 + i64] + v[69] * v[77] * v[78];
            /* 85 = \[OverBracket]_\[Yen]_68(\[Yen]|Rgi)_66 */
            v[235] = 0e0;
            v[236] = v[80];
            v[237] = 0e0;
            v[238] = 0e0;
            v[239] = 0e0;
            v[240] = 0e0;
            v[241] = 0e0;
            v[242] = v[81];
            v[243] = 0e0;
            v[244] = 0e0;
            v[245] = 0e0;
            v[246] = 0e0;
            v[247] = 0e0;
            v[248] = v[82];
            v[249] = 0e0;
            v[250] = 0e0;
            v[251] = 0e0;
            v[252] = 0e0;
            v[253] = 0e0;
            v[254] = v[83];
            v[255] = 0e0;
            v[256] = 0e0;
            v[257] = 0e0;
            v[258] = 0e0;
            v[85] = v[234 + i64] + v[68] * v[77] * v[84];
            /* 87 = \[OverBracket]_\[Yen]_67(\[Yen]|Rgi)_66 */
            v[259] = v[80];
            v[260] = 0e0;
            v[261] = 0e0;
            v[262] = 0e0;
            v[263] = 0e0;
            v[264] = 0e0;
            v[265] = v[81];
            v[266] = 0e0;
            v[267] = 0e0;
            v[268] = 0e0;
            v[269] = 0e0;
            v[270] = 0e0;
            v[271] = v[82];
            v[272] = 0e0;
            v[273] = 0e0;
            v[274] = 0e0;
            v[275] = 0e0;
            v[276] = 0e0;
            v[277] = v[83];
            v[278] = 0e0;
            v[279] = 0e0;
            v[280] = 0e0;
            v[281] = 0e0;
            v[282] = 0e0;
            v[87] = v[258 + i64] + v[67] * v[77] * v[86];
            v[283] = -(v[43] * v[87]);
            v[284] = -(v[43] * v[85]);
            v[285] = -(v[43] * v[79]);
            v[286] = 0e0;
            v[287] = 0e0;
            v[288] = 0e0;
            v[289] = -(v[44] * v[87]);
            v[290] = -(v[44] * v[85]);
            v[291] = -(v[44] * v[79]);
            v[292] = 0e0;
            v[293] = 0e0;
            v[294] = 0e0;
            v[295] = v[45] * v[87];
            v[296] = v[45] * v[85];
            v[297] = v[45] * v[79];
            v[298] = 0e0;
            v[299] = 0e0;
            v[300] = 0e0;
            v[301] = v[46] * v[87];
            v[302] = v[46] * v[85];
            v[303] = v[46] * v[79];
            v[304] = 0e0;
            v[305] = 0e0;
            v[306] = 0e0;
            v[88] = 0e0;/*debug*/
            CR[i64 - 1] += v[71];
            v[72] = 0e0;/*debug*/
            for (i73 = 1; i73 <= 24; i73++) {
                v[73] = i73;
                /* 75 = \[DoubleStruckCapitalG]_j */
                v[75] = v[107 + i73];
                /* 89 = Kgij */
                v[89] = v[282 + i73];
                CT[i64 - 1][i73 - 1] += v[89];
                v[90] = 0e0;/*debug*/
            };/* end for */
        };/* end for */

        //Assembly
        int MNNODE = D[3];
        std::cout << g << std::endl;
        if (g < 0)
        {
            for (int j = 0; j < 6; j++)
            {
                GR->coeffRef(6 * mnode1 + j) += CR[j];
                GR->coeffRef(6 * mnode2 + j) += CR[j];
                GR->coeffRef(6 * (snode1 + MNNODE) + j) += CR[j + 12];
                GR->coeffRef(6 * (snode2 + MNNODE) + j) += CR[j + 12];
                for (int k = 0; k < 6; k++)
                {
                    GT->coeffRef(6 * mnode1 + j, 6 * mnode1 + k) += CT[j][k];
                    GT->coeffRef(6 * mnode1 + j, 6 * mnode2 + k) += CT[j][k + 6];
                    GT->coeffRef(6 * mnode1 + j, 6 * snode1 + k) += CT[j][k + 12];
                    GT->coeffRef(6 * mnode1 + j, 6 * snode2 + k) += CT[j][k + 18];

                    GT->coeffRef(6 * mnode2 + j, 6 * mnode1 + k) += CT[j + 6][k];
                    GT->coeffRef(6 * mnode2 + j, 6 * mnode2 + k) += CT[j + 6][k + 6];
                    GT->coeffRef(6 * mnode2 + j, 6 * snode1 + k) += CT[j + 6][k + 12];
                    GT->coeffRef(6 * mnode2 + j, 6 * snode2 + k) += CT[j + 6][k + 18];

                    GT->coeffRef(6 * snode1 + j, 6 * mnode1 + k) += CT[j + 12][k];
                    GT->coeffRef(6 * snode1 + j, 6 * mnode2 + k) += CT[j + 12][k + 6];
                    GT->coeffRef(6 * snode1 + j, 6 * snode1 + k) += CT[j + 12][k + 12];
                    GT->coeffRef(6 * snode1 + j, 6 * snode2 + k) += CT[j + 12][k + 18];

                    GT->coeffRef(6 * mnode1 + j, 6 * mnode1 + k) += CT[j + 18][k];
                    GT->coeffRef(6 * mnode1 + j, 6 * mnode2 + k) += CT[j + 18][k + 6];
                    GT->coeffRef(6 * mnode1 + j, 6 * snode1 + k) += CT[j + 18][k + 12];
                    GT->coeffRef(6 * mnode1 + j, 6 * snode2 + k) += CT[j + 18][k + 18];
                }
            }
        }

    }
    else 
        return;
};

//Interaction between endpoints for segment to segment contact
//Function to generate tangent stiffness matrix and residual vector
void Contact_STS_Endpoints(double D[5], double X1[3]
    , double X2[3], double u1[6], double u2[6],
    Eigen::SparseMatrix<double, Eigen::ColMajor> *GT, Eigen::VectorXd *GR, int mnode, int snode)
{
    double CR[12];
    double CT[12][12];

    //Initialize contact residual and contact tangent matrix
    for (int j = 0; j < 12; j++)
    {
        CR[j] = 0;
        for (int k = 0; k < 12; k++)
            CT[j][k] = 0;
    }

    /*for (int j = 0; j < 3; j++)
    std::cout << X1[j] << " ";
    std::cout << std::endl;
    for (int j = 0; j < 3; j++)
        std::cout << X2[j] << " ";
    std::cout << std::endl;

    for (int j = 0; j < 6; j++)
        std::cout << u1[j] << " ";
    std::cout << std::endl;
    for (int j = 0; j < 6; j++)
        std::cout << u2[j] << " ";
    std::cout << std::endl;

    for (int j = 0; j < 5; j++)
        std::cout << D[j] << " ";
    std::cout << std::endl;*/

    double v[249];
    int i30, i43;
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
    /* 7 = u01_1 */
    v[7] = u1[0];
    /* 8 = u01_2 */
    v[8] = u1[1];
    /* 9 = u01_3 */
    v[9] = u1[2];
    /* 10 = u01_4 */
    v[10] = u1[3];
    /* 11 = u01_5 */
    v[11] = u1[4];
    /* 12 = u01_6 */
    v[12] = u1[5];
    /* 13 = u02_1 */
    v[13] = u2[0];
    /* 14 = u02_2 */
    v[14] = u2[1];
    /* 15 = u02_3 */
    v[15] = u2[2];
    /* 16 = u02_4 */
    v[16] = u2[3];
    /* 17 = u02_5 */
    v[17] = u2[4];
    /* 18 = u02_6 */
    v[18] = u2[5];
    v[77] = v[7];
    v[78] = v[8];
    v[79] = v[9];
    v[80] = v[10];
    v[81] = v[11];
    v[82] = v[12];
    v[83] = v[13];
    v[84] = v[14];
    v[85] = v[15];
    v[86] = v[16];
    v[87] = v[17];
    v[88] = v[18];
    v[44] = 0e0;/*debug*/
    /* 19 = [\[Epsilon]][Data_1] */
    v[19] = D[0];
    /* 20 = [dm][Data_2] */
    v[20] = D[1];
    /* 21 = [ds][Data_3] */
    v[21] = D[2];
    /* 22 = xm_1 */
    v[22] = v[1] + v[7];
    /* 23 = xm_2 */
    v[23] = v[2] + v[8];
    /* 24 = xm_3 */
    v[24] = v[3] + v[9];
    /* 25 = xs_1 */
    v[25] = v[13] + v[4];
    v[33] = v[22] - v[25];
    /* 26 = xs_2 */
    v[26] = v[14] + v[5];
    v[36] = v[23] - v[26];
    /* 27 = xs_3 */
    v[27] = v[15] + v[6];
    v[38] = v[24] - v[27];
    v[32] = sqrt((v[33] * v[33]) + (v[36] * v[36]) + (v[38] * v[38]));
    v[55] = 1e0 / v[32];
    v[53] = 1e0 / v[32];
    v[51] = 1e0 / v[32];
    v[49] = 1e0 / (v[32] * v[32]);
    v[35] = 1e0 / v[32];
    /* 28 = gN */
    v[28] = -v[20] - v[21] + v[32];
    double g = v[28];
    if (g > 0)
        return;
    v[39] = 1e0 * v[19] * v[28] * v[35] * v[38];
    v[37] = 1e0 * v[19] * v[28] * v[35] * v[36];
    v[34] = 1e0 * v[19] * v[28] * v[33] * v[35];
    v[89] = v[34];
    v[90] = v[37];
    v[91] = v[39];
    v[92] = 0e0;
    v[93] = 0e0;
    v[94] = 0e0;
    v[95] = -v[34];
    v[96] = -v[37];
    v[97] = -v[39];
    v[98] = 0e0;
    v[99] = 0e0;
    v[100] = 0e0;
    v[40] = 0e0;/*debug*/
    /* 29 = \[CapitalPi] */
    v[29] = 0.5e0 * v[19] * (v[28] * v[28]);
    for (i30 = 1; i30 <= 12; i30++) {
        v[30] = i30;
        /* 31 = \[DoubleStruckCapitalG]_i */
        v[31] = v[76 + i30];
        /* 41 = Rgi */
        v[41] = v[88 + i30];
        /* 46 = \[OverBracket]_\[Yen]_34(Rgi|Rgi) */
        v[108] = 1e0;
        v[109] = 0e0;
        v[110] = 0e0;
        v[111] = 0e0;
        v[112] = 0e0;
        v[113] = 0e0;
        v[114] = -1e0;
        v[115] = 0e0;
        v[116] = 0e0;
        v[117] = 0e0;
        v[118] = 0e0;
        v[119] = 0e0;
        v[46] = v[107 + i30];
        /* 47 = \[OverBracket]_\[Yen]_37(Rgi|Rgi) */
        v[120] = 0e0;
        v[121] = 1e0;
        v[122] = 0e0;
        v[123] = 0e0;
        v[124] = 0e0;
        v[125] = 0e0;
        v[126] = 0e0;
        v[127] = -1e0;
        v[128] = 0e0;
        v[129] = 0e0;
        v[130] = 0e0;
        v[131] = 0e0;
        v[47] = v[119 + i30];
        /* 48 = \[OverBracket]_\[Yen]_39(Rgi|Rgi) */
        v[132] = 0e0;
        v[133] = 0e0;
        v[134] = 1e0;
        v[135] = 0e0;
        v[136] = 0e0;
        v[137] = 0e0;
        v[138] = 0e0;
        v[139] = 0e0;
        v[140] = -1e0;
        v[141] = 0e0;
        v[142] = 0e0;
        v[143] = 0e0;
        v[48] = v[131 + i30];
        /* 50 = \[OverBracket]_\[Yen]_32(gN|Rgi) */
        v[50] = (1e0 * v[19] * v[33] * v[46]) / v[32] + (1e0 * v[19] * v[36] * v[47]) / v[32]
            - 1e0 * v[19] * v[28] * v[33] * v[46] * v[49] - 1e0 * v[19] * v[28] * v[36] * v[47] * v[49]
            - 1e0 * v[19] * v[28] * v[38] * v[48] * v[49] + 1e0 * v[19] * v[38] * v[48] * v[51];
        /* 52 = \[OverBracket]_\[Yen]_38(\[Yen]|Rgi)_32 */
        v[52] = v[38] * v[50] * v[51] + 1e0 * v[19] * v[28] * v[48] * v[53];
        /* 54 = \[OverBracket]_\[Yen]_36(\[Yen]|Rgi)_32 */
        v[54] = v[36] * v[50] * v[53] + 1e0 * v[19] * v[28] * v[47] * v[55];
        /* 56 = \[OverBracket]_\[Yen]_33(\[Yen]|Rgi)_32 */
        v[56] = (1e0 * v[19] * v[28] * v[46]) / v[32] + v[33] * v[50] * v[55];
        v[144] = v[56];
        v[145] = v[54];
        v[146] = v[52];
        v[147] = 0e0;
        v[148] = 0e0;
        v[149] = 0e0;
        v[150] = -v[56];
        v[151] = -v[54];
        v[152] = -v[52];
        v[153] = 0e0;
        v[154] = 0e0;
        v[155] = 0e0;
        v[57] = 0e0;/*debug*/
        CR[i30 - 1] += v[41];
        v[42] = 0e0;/*debug*/
        for (i43 = 1; i43 <= 12; i43++) {
            v[43] = i43;
            /* 45 = \[DoubleStruckCapitalG]_j */
            v[45] = v[76 + i43];
            /* 58 = Kgij */
            v[58] = v[143 + i43];
            CT[i30 - 1][i43 - 1] += v[58];
            v[59] = 0e0;/*debug*/
        };/* end for */
    };/* end for */

    int MNNODE = D[3];
    int SNNODE = D[4];
    std::cout << g << std::endl;
    for (int j = 0; j < 6; j++)
    {
        GR->coeffRef(6 * mnode + j) += CR[j];
        GR->coeffRef(6 * (snode + MNNODE) + j) += CR[j + 6];
        for (int k = 0; k < 6; k++)
        {
            GT->coeffRef(6 * mnode + j, 6 * mnode + k) += CT[j][k];
            GT->coeffRef(6 * mnode + j, 6 * (snode + MNNODE) + k) += CT[j][k + 6];
            GT->coeffRef(6 * (snode + MNNODE) + j, 6 * mnode + k) += CT[j + 6][k];
            GT->coeffRef(6 * (snode + MNNODE) + j, 6 * (snode + MNNODE) + k) += CT[j + 6][k + 6];
        }
    }
};







//Closest point projection for segment to segment contact 
/******************* S U B R O U T I N E *********************/
void CPP_STS(double D[3], double X1[3], double X2[3]
    , double X3[3], double X4[3], double u1[6], double u2[6], double u3[6]
    , double u4[6], double* RP, double** TP)
{

    for (int j = 0; j < 2; j++)
    {
        RP[j] = 0;
        for (int k = 0; k < 2; k++)
            TP[j][k] = 0;
    }

    double exi = D[1];
    double eta = D[2];
    //master current config 
    double x1[3], x2[3];
    //slave current config
    double y1[3], y2[3];
    for (int j = 0; j < 3; j++)
    {
        x1[j] = X1[j] + u1[j];
        x2[j] = X2[j] + u2[j];

        y1[j] = X3[j] + u3[j];
        y2[j] = X4[j] + u4[j];
    }
    Eigen::VectorXd bm(3), tm(3), ts(3), bs(3);
    for (int j = 0; j < 3; j++)
    {
        bm(j) = x2[j] + x1[j];
        bs(j) = y2[j] + y1[j];
        tm(j) = x2[j] - x1[j];
        ts(j) = y2[j] - y1[j];
    }
    Eigen::VectorXd exp(3), exp2(3);
    exp2 = (tm * (tm.dot(ts)) - ts * (tm.dot(tm))) / (tm.dot(tm) * ts.dot(ts) - tm.dot(ts) * tm.dot(ts));
    exp = (ts * (tm.dot(ts)) - tm * (ts.dot(ts))) / (tm.dot(tm) * ts.dot(ts) - tm.dot(ts) * tm.dot(ts));

    RP[0] = (bm - bs).dot(exp);
    RP[1] = -(bm - bs).dot(exp2);
};

#endif
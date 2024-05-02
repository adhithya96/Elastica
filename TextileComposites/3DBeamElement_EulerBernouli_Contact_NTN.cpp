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

/******************* S U B R O U T I N E *********************/
//D[1] - penalty parameter
//{D[2],D[3],D[4]} - normal vector
// D[5] - slave beam radius
// D[6] - master beam radius
//X1 - position vector of slave node
//X2 - position vector of master node
//u1 - displacement vector of slave node
//u2 - displacement vector of master node
//CR - residual vector for contact
//CT - tangent matrix for contact
//Tangent and residual for node to node contact
/******************* S U B R O U T I N E *********************/
//void BeamContact::Contact_NTN( double D[6], double X1[3]
//    , double X2[3], double u1[6], double u2[6], double* CR, double** CT, double *g)
//{
//    double v[255];
//    int i42, i50, b38, b39;
//    /* 1 = X01_1 */
//    v[1] = X1[0];
//    /* 2 = X01_2 */
//    v[2] = X1[1];
//    /* 3 = X01_3 */
//    v[3] = X1[2];
//    /* 4 = X02_1 */
//    v[4] = X2[0];
//    /* 5 = X02_2 */
//    v[5] = X2[1];
//    /* 6 = X02_3 */
//    v[6] = X2[2];
//    /* 7 = u01_1 */
//    v[7] = u1[0];
//    /* 8 = u01_2 */
//    v[8] = u1[1];
//    /* 9 = u01_3 */
//    v[9] = u1[2];
//    /* 10 = u01_4 */
//    v[10] = u1[3];
//    /* 11 = u01_5 */
//    v[11] = u1[4];
//    /* 12 = u01_6 */
//    v[12] = u1[5];
//    /* 13 = u02_1 */
//    v[13] = u2[0];
//    /* 14 = u02_2 */
//    v[14] = u2[1];
//    /* 15 = u02_3 */
//    v[15] = u2[2];
//    /* 16 = u02_4 */
//    v[16] = u2[3];
//    /* 17 = u02_5 */
//    v[17] = u2[4];
//    /* 18 = u02_6 */
//    v[18] = u2[5];
//    v[77] = v[7];
//    v[78] = v[8];
//    v[79] = v[9];
//    v[80] = v[10];
//    v[81] = v[11];
//    v[82] = v[12];
//    v[83] = v[13];
//    v[84] = v[14];
//    v[85] = v[15];
//    v[86] = v[16];
//    v[87] = v[17];
//    v[88] = v[18];
//    v[51] = 0e0;/*debug*/
//    /* 19 = Data_1 */
//    v[19] = D[0];
//    /* 20 = Data_2 */
//    v[20] = D[1];
//    /* 21 = Data_3 */
//    v[21] = D[2];
//    /* 22 = Data_4 */
//    v[22] = D[3];
//    /* 23 = Data_5 */
//    v[23] = D[4];
//    /* 24 = Data_6 */
//    v[24] = D[5];
//    /* 25 = \[Epsilon] */
//    v[25] = v[19];
//    /* 26 = n_1 */
//    v[26] = v[20];
//    /* 27 = n_2 */
//    v[27] = v[21];
//    /* 28 = n_3 */
//    v[28] = v[22];
//    /* 29 = rs */
//    v[29] = v[23];
//    /* 30 = rm */
//    v[30] = v[24];
//    /* 31 = x1_1 */
//    v[31] = v[1] + v[7];
//    /* 32 = x1_2 */
//    v[32] = v[2] + v[8];
//    /* 33 = x1_3 */
//    v[33] = v[3] + v[9];
//    /* 34 = x2_1 */
//    v[34] = v[13] + v[4];
//    /* 35 = x2_2 */
//    v[35] = v[14] + v[5];
//    /* 36 = x2_3 */
//    v[36] = v[15] + v[6];
//    /* 37 = gN */
//    v[37] = -v[29] - v[30] + v[26] * (-v[31] + v[34]) + v[27] * (-v[32] + v[35]) + v[28] * (-v[33] + v[36]);
//    *g = v[37];
//    v[46] = 1e0 * v[25] * v[28] * v[37];
//    v[45] = 1e0 * v[25] * v[27] * v[37];
//    v[44] = 1e0 * v[25] * v[26] * v[37];
//    v[89] = -v[44];
//    v[90] = -v[45];
//    v[91] = -v[46];
//    v[92] = 0e0;
//    v[93] = 0e0;
//    v[94] = 0e0;
//    v[95] = v[44];
//    v[96] = v[45];
//    v[97] = v[46];
//    v[98] = 0e0;
//    v[99] = 0e0;
//    v[100] = 0e0;
//    v[47] = 0e0;/*debug*/
//    b38 = v[37] > 0e0;
//    v[38] = b38;
//    b39 = b38;
//    v[39] = b39;
//    if (b39) {
//        return;
//        v[40] = 0e0;/*debug*/
//    }
//    else {
//    };
//    /* 41 = \[CapitalPi] */
//    v[41] = 0.5e0 * v[25] * (v[37] * v[37]);
//    for (i42 = 1; i42 <= 12; i42++) {
//        v[42] = i42;
//        /* 43 = \[DoubleStruckCapitalG]_i */
//        v[43] = v[76 + i42];
//        /* 48 = Rgi */
//        v[48] = v[88 + i42];
//        /* 53 = \[OverBracket]_gN_(\[Yen]|Rgi)_46 */
//        v[108] = -1e0;
//        v[109] = 0e0;
//        v[110] = 0e0;
//        v[111] = 0e0;
//        v[112] = 0e0;
//        v[113] = 0e0;
//        v[114] = 1e0;
//        v[115] = 0e0;
//        v[116] = 0e0;
//        v[117] = 0e0;
//        v[118] = 0e0;
//        v[119] = 0e0;
//        v[120] = 0e0;
//        v[121] = -1e0;
//        v[122] = 0e0;
//        v[123] = 0e0;
//        v[124] = 0e0;
//        v[125] = 0e0;
//        v[126] = 0e0;
//        v[127] = 1e0;
//        v[128] = 0e0;
//        v[129] = 0e0;
//        v[130] = 0e0;
//        v[131] = 0e0;
//        v[132] = 0e0;
//        v[133] = 0e0;
//        v[134] = -1e0;
//        v[135] = 0e0;
//        v[136] = 0e0;
//        v[137] = 0e0;
//        v[138] = 0e0;
//        v[139] = 0e0;
//        v[140] = 1e0;
//        v[141] = 0e0;
//        v[142] = 0e0;
//        v[143] = 0e0;
//        v[53] = 1e0 * v[107 + i42] * v[25] * v[26] + 1e0 * v[119 + i42] * v[25] * v[27] + 1e0 * v[131 + i42] * v[25] * v[28];
//        v[56] = v[28] * v[53];
//        v[55] = v[27] * v[53];
//        v[54] = v[26] * v[53];
//        v[144] = -v[54];
//        v[145] = -v[55];
//        v[146] = -v[56];
//        v[147] = 0e0;
//        v[148] = 0e0;
//        v[149] = 0e0;
//        v[150] = v[54];
//        v[151] = v[55];
//        v[152] = v[56];
//        v[153] = 0e0;
//        v[154] = 0e0;
//        v[155] = 0e0;
//        v[57] = 0e0;/*debug*/
//        CR[i42 - 1] += v[48];
//        v[49] = 0e0;/*debug*/
//        for (i50 = 1; i50 <= 12; i50++) {
//            v[50] = i50;
//            /* 52 = \[DoubleStruckCapitalG]_j */
//            v[52] = v[76 + i50];
//            /* 58 = Kgij */
//            v[58] = v[143 + i50];
//            CT[i42 - 1][i50 - 1] += v[58];
//            v[59] = 0e0;/*debug*/
//        };/* end for */
//    };/* end for */
//};


//void BeamContact::Contact_NTN(double D[6], double X1[3]
//        , double X2[3], double u1[6], double u2[6], double* CR, double** CT, double* g)
//{
//    double v[255];
//    int i43, i56, b39, b40;
//    /* 1 = X01_1 */
//    v[1] = X1[0];
//    /* 2 = X01_2 */
//    v[2] = X1[1];
//    /* 3 = X01_3 */
//    v[3] = X1[2];
//    /* 4 = X02_1 */
//    v[4] = X2[0];
//    /* 5 = X02_2 */
//    v[5] = X2[1];
//    /* 6 = X02_3 */
//    v[6] = X2[2];
//    /* 7 = u01_1 */
//    v[7] = u1[0];
//    /* 8 = u01_2 */
//    v[8] = u1[1];
//    /* 9 = u01_3 */
//    v[9] = u1[2];
//    /* 10 = u01_4 */
//    v[10] = u1[3];
//    /* 11 = u01_5 */
//    v[11] = u1[4];
//    /* 12 = u01_6 */
//    v[12] = u1[5];
//    /* 13 = u02_1 */
//    v[13] = u2[0];
//    /* 14 = u02_2 */
//    v[14] = u2[1];
//    /* 15 = u02_3 */
//    v[15] = u2[2];
//    /* 16 = u02_4 */
//    v[16] = u2[3];
//    /* 17 = u02_5 */
//    v[17] = u2[4];
//    /* 18 = u02_6 */
//    v[18] = u2[5];
//    v[87] = v[7];
//    v[88] = v[8];
//    v[89] = v[9];
//    v[90] = v[10];
//    v[91] = v[11];
//    v[92] = v[12];
//    v[93] = v[13];
//    v[94] = v[14];
//    v[95] = v[15];
//    v[96] = v[16];
//    v[97] = v[17];
//    v[98] = v[18];
//    v[57] = 0e0;/*debug*/
//    /* 19 = Data_1 */
//    v[19] = D[0];
//    /* 20 = Data_2 */
//    v[20] = D[1];
//    /* 21 = Data_3 */
//    v[21] = D[2];
//    /* 22 = Data_4 */
//    v[22] = D[3];
//    /* 23 = Data_5 */
//    v[23] = D[4];
//    /* 24 = Data_6 */
//    v[24] = D[5];
//    /* 25 = \[Epsilon] */
//    v[25] = v[19];
//    /* 26 = n_1 */
//    v[26] = v[20];
//    /* 27 = n_2 */
//    v[27] = v[21];
//    /* 28 = n_3 */
//    v[28] = v[22];
//    /* 29 = rs */
//    v[29] = v[23];
//    /* 30 = rm */
//    v[30] = v[24];
//    /* 31 = x1_1 */
//    v[31] = v[1] + v[7];
//    /* 32 = x1_2 */
//    v[32] = v[2] + v[8];
//    /* 33 = x1_3 */
//    v[33] = v[3] + v[9];
//    /* 34 = x2_1 */
//    v[34] = v[13] + v[4];
//    v[47] = -v[31] + v[34];
//    /* 35 = x2_2 */
//    v[35] = v[14] + v[5];
//    v[49] = -v[32] + v[35];
//    /* 36 = x2_3 */
//    v[36] = v[15] + v[6];
//    v[51] = -v[33] + v[36];
//    v[37] = (v[47] * v[47]) + (v[49] * v[49]) + (v[51] * v[51]);
//    v[45] = 1e0 / sqrt(v[37]);
//    /* 38 = gN */
//    v[38] = -v[29] - v[30] + v[37] * v[45];
//    *g = v[38];
//    b39 = v[38] > 0e0;
//    v[39] = b39;
//    b40 = b39;
//    v[40] = b40;
//    if (b40) {
//        return;
//        v[41] = 0e0;/*debug*/
//    }
//    else {
//    };
//    /* 42 = \[CapitalPi] */
//    v[42] = 0.5e0 * v[25] * (v[38] * v[38]);
//    /* 46 = \[OverBracket]_\[Yen]_37(\[CapitalPi]|\[CapitalPi]) */
//    v[46] = 0.5e0 * v[25] * v[38] * v[45];
//    v[52] = 2e0 * v[46] * v[51];
//    v[50] = 2e0 * v[46] * v[49];
//    v[48] = 2e0 * v[46] * v[47];
//    v[99] = -v[48];
//    v[100] = -v[50];
//    v[101] = -v[52];
//    v[102] = 0e0;
//    v[103] = 0e0;
//    v[104] = 0e0;
//    v[105] = v[48];
//    v[106] = v[50];
//    v[107] = v[52];
//    v[108] = 0e0;
//    v[109] = 0e0;
//    v[110] = 0e0;
//    v[53] = 0e0;/*debug*/
//    for (i43 = 1; i43 <= 12; i43++) {
//        v[43] = i43;
//        /* 44 = \[DoubleStruckCapitalG]_i */
//        v[44] = v[86 + i43];
//        /* 54 = Rgi */
//        v[54] = v[98 + i43];
//        /* 59 = \[OverBracket]_\[Yen]_48(Rgi|Rgi) */
//        v[118] = -1e0;
//        v[119] = 0e0;
//        v[120] = 0e0;
//        v[121] = 0e0;
//        v[122] = 0e0;
//        v[123] = 0e0;
//        v[124] = 1e0;
//        v[125] = 0e0;
//        v[126] = 0e0;
//        v[127] = 0e0;
//        v[128] = 0e0;
//        v[129] = 0e0;
//        v[59] = v[117 + i43];
//        /* 60 = \[OverBracket]_\[Yen]_50(Rgi|Rgi) */
//        v[130] = 0e0;
//        v[131] = -1e0;
//        v[132] = 0e0;
//        v[133] = 0e0;
//        v[134] = 0e0;
//        v[135] = 0e0;
//        v[136] = 0e0;
//        v[137] = 1e0;
//        v[138] = 0e0;
//        v[139] = 0e0;
//        v[140] = 0e0;
//        v[141] = 0e0;
//        v[60] = v[129 + i43];
//        /* 61 = \[OverBracket]_\[Yen]_52(Rgi|Rgi) */
//        v[142] = 0e0;
//        v[143] = 0e0;
//        v[144] = -1e0;
//        v[145] = 0e0;
//        v[146] = 0e0;
//        v[147] = 0e0;
//        v[148] = 0e0;
//        v[149] = 0e0;
//        v[150] = 1e0;
//        v[151] = 0e0;
//        v[152] = 0e0;
//        v[153] = 0e0;
//        v[61] = v[141 + i43];
//        /* 62 = \[OverBracket]_\[OverBracket]_\[Yen]_37(\[CapitalPi]|\[CapitalPi])(\[Yen]|Rgi)_52 */
//        v[62] = 2e0 * v[47] * v[59] + 2e0 * v[49] * v[60] + 2e0 * v[51] * v[61];
//        /* 63 = \[OverBracket]_\[Yen]_37(\[Yen]|Rgi)_45 */
//        v[63] = -0.5e0 * (v[45] * (0.5e0 * v[25] * v[38] + 0.5e0 * v[25] * v[37] * v[45]) * v[62]) / v[37] + 0.5e0 * v[25] *
//            (v[45] * v[45]) * v[62];
//        /* 64 = \[OverBracket]_\[Yen]_51(\[Yen]|Rgi)_37 */
//        v[64] = 2e0 * v[46] * v[61] + 2e0 * v[51] * v[63];
//        /* 65 = \[OverBracket]_\[Yen]_49(\[Yen]|Rgi)_37 */
//        v[65] = 2e0 * v[46] * v[60] + 2e0 * v[49] * v[63];
//        /* 66 = \[OverBracket]_\[Yen]_47(\[Yen]|Rgi)_37 */
//        v[66] = 2e0 * v[46] * v[59] + 2e0 * v[47] * v[63];
//        v[154] = -v[66];
//        v[155] = -v[65];
//        v[156] = -v[64];
//        v[157] = 0e0;
//        v[158] = 0e0;
//        v[159] = 0e0;
//        v[160] = v[66];
//        v[161] = v[65];
//        v[162] = v[64];
//        v[163] = 0e0;
//        v[164] = 0e0;
//        v[165] = 0e0;
//        v[67] = 0e0;/*debug*/
//        CR[i43 - 1] += v[54];
//        v[55] = 0e0;/*debug*/
//        for (i56 = 1; i56 <= 12; i56++) {
//            v[56] = i56;
//            /* 58 = \[DoubleStruckCapitalG]_j */
//            v[58] = v[86 + i56];
//            /* 68 = Kgij */
//            v[68] = v[153 + i56];
//            CT[i43 - 1][i56 - 1] += v[68];
//            v[69] = 0e0;/*debug*/
//        };/* end for */
//    };/* end for */
//};


void BeamContact::Contact_NTN(double D[6], double X1[3]
    , double X2[3], double u1[6], double u2[6], double *CR, double **CT, double* g)
{
    double v[255];
    int i42, i50, b38, b39;
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
    v[51] = 0e0;/*debug*/
    /* 19 = Data_1 */
    v[19] = D[0];
    /* 20 = Data_2 */
    v[20] = D[1];
    /* 21 = Data_3 */
    v[21] = D[2];
    /* 22 = Data_4 */
    v[22] = D[3];
    /* 23 = Data_5 */
    v[23] = D[4];
    /* 24 = Data_6 */
    v[24] = D[5];
    /* 25 = \[Epsilon] */
    v[25] = v[19];
    /* 26 = n_1 */
    v[26] = v[20];
    /* 27 = n_2 */
    v[27] = v[21];
    /* 28 = n_3 */
    v[28] = v[22];
    /* 29 = rs */
    v[29] = v[23];
    /* 30 = rm */
    v[30] = v[24];
    /* 31 = x1_1 */
    v[31] = v[1] + v[7];
    /* 32 = x1_2 */
    v[32] = v[2] + v[8];
    /* 33 = x1_3 */
    v[33] = v[3] + v[9];
    /* 34 = x2_1 */
    v[34] = v[13] + v[4];
    /* 35 = x2_2 */
    v[35] = v[14] + v[5];
    /* 36 = x2_3 */
    v[36] = v[15] + v[6];
    /* 37 = gN */
    v[37] = -v[29] - v[30] + v[26] * (-v[31] + v[34]) + v[27] * (-v[32] + v[35]) + v[28] * (-v[33] + v[36]);
    *g = v[37];
    v[46] = 1e0 * v[25] * v[28] * v[37];
    v[45] = 1e0 * v[25] * v[27] * v[37];
    v[44] = 1e0 * v[25] * v[26] * v[37];
    v[89] = -v[44];
    v[90] = -v[45];
    v[91] = -v[46];
    v[92] = 0e0;
    v[93] = 0e0;
    v[94] = 0e0;
    v[95] = v[44];
    v[96] = v[45];
    v[97] = v[46];
    v[98] = 0e0;
    v[99] = 0e0;
    v[100] = 0e0;
    v[47] = 0e0;/*debug*/
    b38 = v[37] > 0e0;
    v[38] = b38;
    b39 = b38;
    v[39] = b39;
    if (b39) {
        return;
        v[40] = 0e0;/*debug*/
    }
    else {
    };
    /* 41 = \[CapitalPi] */
    v[41] = 0.5e0 * v[25] * (v[37] * v[37]);
    for (i42 = 1; i42 <= 12; i42++) {
        v[42] = i42;
        /* 43 = \[DoubleStruckCapitalG]_i */
        v[43] = v[76 + i42];
        /* 48 = Rgi */
        v[48] = v[88 + i42];
        /* 53 = \[OverBracket]_gN_(\[Yen]|Rgi)_46 */
        v[108] = -1e0;
        v[109] = 0e0;
        v[110] = 0e0;
        v[111] = 0e0;
        v[112] = 0e0;
        v[113] = 0e0;
        v[114] = 1e0;
        v[115] = 0e0;
        v[116] = 0e0;
        v[117] = 0e0;
        v[118] = 0e0;
        v[119] = 0e0;
        v[120] = 0e0;
        v[121] = -1e0;
        v[122] = 0e0;
        v[123] = 0e0;
        v[124] = 0e0;
        v[125] = 0e0;
        v[126] = 0e0;
        v[127] = 1e0;
        v[128] = 0e0;
        v[129] = 0e0;
        v[130] = 0e0;
        v[131] = 0e0;
        v[132] = 0e0;
        v[133] = 0e0;
        v[134] = -1e0;
        v[135] = 0e0;
        v[136] = 0e0;
        v[137] = 0e0;
        v[138] = 0e0;
        v[139] = 0e0;
        v[140] = 1e0;
        v[141] = 0e0;
        v[142] = 0e0;
        v[143] = 0e0;
        v[53] = 1e0 * v[107 + i42] * v[25] * v[26] + 1e0 * v[119 + i42] * v[25] * v[27] + 1e0 * v[131 + i42] * v[25] * v[28];
        v[56] = v[28] * v[53];
        v[55] = v[27] * v[53];
        v[54] = v[26] * v[53];
        v[144] = -v[54];
        v[145] = -v[55];
        v[146] = -v[56];
        v[147] = 0e0;
        v[148] = 0e0;
        v[149] = 0e0;
        v[150] = v[54];
        v[151] = v[55];
        v[152] = v[56];
        v[153] = 0e0;
        v[154] = 0e0;
        v[155] = 0e0;
        v[57] = 0e0;/*debug*/
        CR[i42 - 1] += v[48];
        v[49] = 0e0;/*debug*/
        for (i50 = 1; i50 <= 12; i50++) {
            v[50] = i50;
            /* 52 = \[DoubleStruckCapitalG]_j */
            v[52] = v[76 + i50];
            /* 58 = Kgij */
            v[58] = v[143 + i50];
            CT[i42 - 1][i50 - 1] += v[58];
            v[59] = 0e0;/*debug*/
        };/* end for */
    };/* end for */
};


/*NLEBBE3D::NLEBBE3D(int choice)
{
    if (choice == 1)
    {
        this->NNODE = 5;
        this->NELEM = 2;
        this->NDOF = 6;
        this->NLS = 30;
        this->NEN = 3;
        this->NDM = 3;

        this->NODE = Eigen::MatrixXd::Zero(this->NNODE, this->NDM);
        this->ELEM = Eigen::MatrixXd::Zero(this->NELEM, this->NEN + 1);

        this->NODE(0, 0) = 0;
        this->NODE(0, 1) = 0;
        this->NODE(0, 2) = 0;
        this->NODE(1, 0) = 2.5;
        this->NODE(1, 1) = 0;
        this->NODE(1, 2) = 0;
        this->NODE(2, 0) = 5;
        this->NODE(2, 1) = 0;
        this->NODE(2, 2) = 0;
        this->NODE(3, 0) = 7.5;
        this->NODE(3, 1) = 0;
        this->NODE(3, 2) = 0;
        this->NODE(4, 0) = 10;
        this->NODE(4, 1) = 0;
        this->NODE(4, 2) = 0;

        this->ELEM(0, 0) = 1;
        this->ELEM(0, 1) = 1;
        this->ELEM(0, 2) = 2;
        this->ELEM(0, 3) = 3;
        this->ELEM(1, 0) = 1;
        this->ELEM(1, 1) = 3;
        this->ELEM(1, 2) = 4;
        this->ELEM(1, 3) = 5;

        this->E = 1e7;
        this->nu = 0.0001;
        this->Bp = 1;
        this->Hp = 1;
        this->Zx = 0;
        this->Zy = 0;
        this->Zz = 1;

        this->DIA = 0.5;
    }
    else if (choice == 2)
    {
        this->NNODE = 5;
        this->NELEM = 2;
        this->NDOF = 6;
        this->NLS = 30;
        this->NEN = 3;
        this->NDM = 3;

        this->NODE = Eigen::MatrixXd::Zero(this->NNODE, this->NDM);
        this->ELEM = Eigen::MatrixXd::Zero(this->NELEM, this->NEN + 1);

        this->NODE(0, 0) = 0;
        this->NODE(0, 1) = 2;
        this->NODE(0, 2) = 0;
        this->NODE(1, 0) = 2.5;
        this->NODE(1, 1) = 2;
        this->NODE(1, 2) = 0;
        this->NODE(2, 0) = 5;
        this->NODE(2, 1) = 2;
        this->NODE(2, 2) = 0;
        this->NODE(3, 0) = 7.5;
        this->NODE(3, 1) = 2;
        this->NODE(3, 2) = 0;
        this->NODE(4, 0) = 10;
        this->NODE(4, 1) = 2;
        this->NODE(4, 2) = 0;

        this->ELEM(0, 0) = 1;
        this->ELEM(0, 1) = 1;
        this->ELEM(0, 2) = 2;
        this->ELEM(0, 3) = 3;
        this->ELEM(1, 0) = 1;
        this->ELEM(1, 1) = 3;
        this->ELEM(1, 2) = 4;
        this->ELEM(1, 3) = 5;

        this->E = 1e7;
        this->nu = 0.0001;
        this->Bp = 1;
        this->Hp = 1;
        this->Zx = 0;
        this->Zy = 0;
        this->Zz = 1;

        this->DIA = 0.5;
    }
    else if (choice == 3)
    {
        this->NNODE = 5;
        this->NELEM = 2;
        this->NDOF = 6;
        this->NLS = 30;
        this->NEN = 3;
        this->NDM = 3;

        this->NODE = Eigen::MatrixXd::Zero(this->NNODE, this->NDM);
        this->ELEM = Eigen::MatrixXd::Zero(this->NELEM, this->NEN + 1);

        this->NODE(0, 0) = 0;
        this->NODE(0, 1) = 0;
        this->NODE(0, 2) = 2;
        this->NODE(1, 0) = 2.5;
        this->NODE(1, 1) = 0;
        this->NODE(1, 2) = 2;
        this->NODE(2, 0) = 5;
        this->NODE(2, 1) = 0;
        this->NODE(2, 2) = 2;
        this->NODE(3, 0) = 7.5;
        this->NODE(3, 1) = 0;
        this->NODE(3, 2) = 2;
        this->NODE(4, 0) = 10;
        this->NODE(4, 1) = 0;
        this->NODE(4, 2) = 2;

        this->ELEM(0, 0) = 1;
        this->ELEM(0, 1) = 1;
        this->ELEM(0, 2) = 2;
        this->ELEM(0, 3) = 3;
        this->ELEM(1, 0) = 1;
        this->ELEM(1, 1) = 3;
        this->ELEM(1, 2) = 4;
        this->ELEM(1, 3) = 5;

        this->E = 1e7;
        this->nu = 0.0001;
        this->Bp = 1;
        this->Hp = 1;
        this->Zx = 0;
        this->Zy = 0;
        this->Zz = 1;

        this->DIA = 0.5;
    }
    else if (choice == 4)
    {
        this->NNODE = 5;
        this->NELEM = 2;
        this->NDOF = 6;
        this->NLS = 30;
        this->NEN = 3;
        this->NDM = 3;

        this->NODE = Eigen::MatrixXd::Zero(this->NNODE, this->NDM);
        this->ELEM = Eigen::MatrixXd::Zero(this->NELEM, this->NEN + 1);

        this->NODE(0, 0) = 0;
        this->NODE(0, 1) = 2;
        this->NODE(0, 2) = 2;
        this->NODE(1, 0) = 2.5;
        this->NODE(1, 1) = 2;
        this->NODE(1, 2) = 2;
        this->NODE(2, 0) = 5;
        this->NODE(2, 1) = 2;
        this->NODE(2, 2) = 2;
        this->NODE(3, 0) = 7.5;
        this->NODE(3, 1) = 2;
        this->NODE(3, 2) = 2;
        this->NODE(4, 0) = 10;
        this->NODE(4, 1) = 2;
        this->NODE(4, 2) = 2;

        this->ELEM(0, 0) = 1;
        this->ELEM(0, 1) = 1;
        this->ELEM(0, 2) = 2;
        this->ELEM(0, 3) = 3;
        this->ELEM(1, 0) = 1;
        this->ELEM(1, 1) = 3;
        this->ELEM(1, 2) = 4;
        this->ELEM(1, 3) = 5;

        this->E = 1e7;
        this->nu = 0.0001;
        this->Bp = 1;
        this->Hp = 1;
        this->Zx = 0;
        this->Zy = 0;
        this->Zz = 1;

        this->DIA = 0.5;
    }
}*/

/*NLEBBE3D::NLEBBE3D(int choice, std::string str)
{
    if (choice == 1)
    {
        NLEBBE3D::MAT = str;

        this->NNODE = 5;
        this->NELEM = 2;
        this->NDOF = 6;
        this->NLS = 11;
        this->NEN = 3;
        this->NDM = 3;

        this->NODE = Eigen::MatrixXd::Zero(this->NNODE, this->NDM);
        this->ELEM = Eigen::MatrixXd::Zero(this->NELEM, this->NEN + 1);

        this->NODE(0, 0) = 0;
        this->NODE(0, 1) = 0;
        this->NODE(0, 2) = 0;
        this->NODE(1, 0) = 2.5;
        this->NODE(1, 1) = 0;
        this->NODE(1, 2) = 0;
        this->NODE(2, 0) = 5;
        this->NODE(2, 1) = 0;
        this->NODE(2, 2) = 0;
        this->NODE(3, 0) = 7.5;
        this->NODE(3, 1) = 0;
        this->NODE(3, 2) = 0;
        this->NODE(4, 0) = 10;
        this->NODE(4, 1) = 0;
        this->NODE(4, 2) = 0;

        this->ELEM(0, 0) = 1;
        this->ELEM(0, 1) = 1;
        this->ELEM(0, 2) = 2;
        this->ELEM(0, 3) = 3;
        this->ELEM(1, 0) = 1;
        this->ELEM(1, 1) = 3;
        this->ELEM(1, 2) = 4;
        this->ELEM(1, 3) = 5;

        this->E = 1e7;
        this->nu = 0.0001;
        this->Bp = 1;
        this->Hp = 1;
        this->Zx = 0;
        this->Zy = 0;
        this->Zz = 1;

        this->DIA = 0.5;
    }
    else if (choice == 2)
    {
        NLEBBE3D::MAT = str;

        this->NNODE = 5;
        this->NELEM = 2;
        this->NDOF = 6;
        this->NLS = 11;
        this->NEN = 3;
        this->NDM = 3;

        this->NODE = Eigen::MatrixXd::Zero(this->NNODE, this->NDM);
        this->ELEM = Eigen::MatrixXd::Zero(this->NELEM, this->NEN + 1);

        this->NODE(0, 0) = 0;
        this->NODE(0, 1) = 0.50;
        this->NODE(0, 2) = 0;
        this->NODE(1, 0) = 2.5;
        this->NODE(1, 1) = 0.50;
        this->NODE(1, 2) = 0;
        this->NODE(2, 0) = 5;
        this->NODE(2, 1) = 0.50;
        this->NODE(2, 2) = 0;
        this->NODE(3, 0) = 7.5;
        this->NODE(3, 1) = 0.50;
        this->NODE(3, 2) = 0;
        this->NODE(4, 0) = 10;
        this->NODE(4, 1) = 0.50;
        this->NODE(4, 2) = 0;

        this->ELEM(0, 0) = 1;
        this->ELEM(0, 1) = 1;
        this->ELEM(0, 2) = 2;
        this->ELEM(0, 3) = 3;
        this->ELEM(1, 0) = 1;
        this->ELEM(1, 1) = 3;
        this->ELEM(1, 2) = 4;
        this->ELEM(1, 3) = 5;

        this->E = 1e7;
        this->nu = 0.0001;
        this->Bp = 1;
        this->Hp = 1;
        this->Zx = 0;
        this->Zy = 0;
        this->Zz = 1;

        this->DIA = 0.5;
    }
}*/




#endif

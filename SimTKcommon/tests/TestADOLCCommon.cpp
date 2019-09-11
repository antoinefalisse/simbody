/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-17 Stanford University and the Authors.        *
 * Authors: Chris Dembia                                                      *
 * Contributors: Antoine Falisse, Gil Serrancoli                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

/*
 * These tests ensure that using ADOL-C's adouble scalar type works as intended.
 */

#include "SimTKcommon.h"
#include "SimTKcommon/Testing.h"

#include <iostream>
using std::cout;
using std::endl;
using std::cin;


using namespace SimTK;

void testVectorOperations() {

    // Vectors

    Vector_<adouble> v(2);
    v[0] = 4;
    v[1] = 6;
    SimTK_TEST(v[0] == 4);
    SimTK_TEST(v[1] == 6);
    double b = 3.5;
    int b_int = 3;
    adouble b_ad = 3.5;
    float b_fl = 3.5;

    /// multiplication
    Vector_<adouble> v2 = b * v;
	Vector_<adouble> v3 = v * b;
    Vector_<adouble> v2_int = b_int * v;
    Vector_<adouble> v3_int = v * b_int;
    Vector_<adouble> v2_ad = b_ad * v;
    Vector_<adouble> v3_ad = v * b_ad;
    Vector_<adouble> v2_fl = b_fl * v;
    Vector_<adouble> v3_fl = v * b_fl;
    /// Ensure that v is not modified simply by using it in an operation.
    SimTK_TEST(v[0] == 4);
    SimTK_TEST(v[1] == 6);
    /// Ensure that the correct operator* is used.
    SimTK_TEST(v2[0] == 14);
    SimTK_TEST(v2[1] == 21);
	SimTK_TEST(v3[0] == 14);
	SimTK_TEST(v3[1] == 21);
    SimTK_TEST(v2_int[0] == 12);
    SimTK_TEST(v2_int[1] == 18);
    SimTK_TEST(v3_int[0] == 12);
    SimTK_TEST(v3_int[1] == 18);
    SimTK_TEST(v2_ad[0] == 14);
    SimTK_TEST(v2_ad[1] == 21);
    SimTK_TEST(v3_ad[0] == 14);
    SimTK_TEST(v3_ad[1] == 21);
    SimTK_TEST(v2_fl[0] == 14);
    SimTK_TEST(v2_fl[1] == 21);
    SimTK_TEST(v3_fl[0] == 14);
    SimTK_TEST(v3_fl[1] == 21);

	/// addition
	Vector_<adouble> v4 = b + v;
	Vector_<adouble> v5 = v + b;
    Vector_<adouble> v4_int = b_int + v;
    Vector_<adouble> v5_int = v + b_int;
    Vector_<adouble> v4_ad = b_ad + v;
    Vector_<adouble> v5_ad = v + b_ad;
    Vector_<adouble> v4_fl = b_fl + v;
    Vector_<adouble> v5_fl = v + b_fl;
	/// Ensure that v is not modified simply by using it in an operation.
	SimTK_TEST(v[0] == 4);
	SimTK_TEST(v[1] == 6);
	/// Ensure that the correct operator+ is used.
	SimTK_TEST(v4[0] == 7.5);
	SimTK_TEST(v4[1] == 9.5);
	SimTK_TEST(v5[0] == 7.5);
	SimTK_TEST(v5[1] == 9.5);
    SimTK_TEST(v4_int[0] == 7);
    SimTK_TEST(v4_int[1] == 9);
    SimTK_TEST(v5_int[0] == 7);
    SimTK_TEST(v5_int[1] == 9);
    SimTK_TEST(v4_ad[0] == 7.5);
    SimTK_TEST(v4_ad[1] == 9.5);
    SimTK_TEST(v5_ad[0] == 7.5);
    SimTK_TEST(v5_ad[1] == 9.5);
    SimTK_TEST(v4_fl[0] == 7.5);
    SimTK_TEST(v4_fl[1] == 9.5);
    SimTK_TEST(v5_fl[0] == 7.5);
    SimTK_TEST(v5_fl[1] == 9.5);

	/// substraction
	Vector_<adouble> v6 = b - v;
	Vector_<adouble> v7 = v - b;
    Vector_<adouble> v6_int = b_int - v;
    Vector_<adouble> v7_int = v - b_int;
    Vector_<adouble> v6_ad = b_ad - v;
    Vector_<adouble> v7_ad = v - b_ad;
    Vector_<adouble> v6_fl = b_fl - v;
    Vector_<adouble> v7_fl = v - b_fl;
	/// Ensure that v is not modified simply by using it in an operation.
	SimTK_TEST(v[0] == 4);
	SimTK_TEST(v[1] == 6);
    /// Ensure that the correct operator- is used.
	SimTK_TEST(v6[0] == -0.5);
	SimTK_TEST(v6[1] == -2.5);
	SimTK_TEST(v7[0] == 0.5);
	SimTK_TEST(v7[1] == 2.5);
    SimTK_TEST(v6_int[0] == -1);
    SimTK_TEST(v6_int[1] == -3);
    SimTK_TEST(v7_int[0] == 1);
    SimTK_TEST(v7_int[1] == 3);
    SimTK_TEST(v6_ad[0] == -0.5);
    SimTK_TEST(v6_ad[1] == -2.5);
    SimTK_TEST(v7_ad[0] == 0.5);
    SimTK_TEST(v7_ad[1] == 2.5);
    SimTK_TEST(v6_fl[0] == -0.5);
    SimTK_TEST(v6_fl[1] == -2.5);
    SimTK_TEST(v7_fl[0] == 0.5);
    SimTK_TEST(v7_fl[1] == 2.5);

	/// division
    Vector_<adouble> v1(2);
    Vector_<adouble> v1_int(2);
    v1[0] = 7.5;
    v1[1] = 11.25;
    v1_int[0] = 8;
    v1_int[1] = 12;
    SimTK_TEST(v1[0] == 7.5);
    SimTK_TEST(v1[1] == 11.25);
    SimTK_TEST(v1_int[0] == 8);
    SimTK_TEST(v1_int[1] == 12);
    double b2 = 2.5;
    int b2_int = 2;
    adouble b2_ad = 2.5;
    float b2_fl = 2.5;

	Vector_<adouble> v8 = v1 / b2;
    Vector_<adouble> v8_int = v1_int / b2_int;
    Vector_<adouble> v8_ad = v1 / b2_ad;
    Vector_<adouble> v8_fl = v1 / b2_fl;
	/// Ensure that v1 is not modified simply by using it in an operation.
	SimTK_TEST(v1[0] == 7.5);
	SimTK_TEST(v1[1] == 11.25);
    SimTK_TEST(v1_int[0] == 8);
    SimTK_TEST(v1_int[1] == 12);
	/// Ensure that the correct operator/ is used.
    SimTK_TEST(v8[0] == 3);
    SimTK_TEST(v8[1] == 4.5);
    SimTK_TEST(v8_int[0] == 4);
    SimTK_TEST(v8_int[1] == 6);
    SimTK_TEST(v8_ad[0] == 3);
    SimTK_TEST(v8_ad[1] == 4.5);
    SimTK_TEST(v8_fl[0] == 3);
    SimTK_TEST(v8_fl[1] == 4.5);


    // Matrices
    Matrix_<adouble> m(2, 3);
    m[0][0] = 6;
    m[1][0] = 2;
    m[0][1] = 3.5;
    m[1][1] = 1.5;
    m[0][2] = 2.5;
    m[1][2] = 0.5;
    Matrix_<double> md(2, 3);
    md[0][0] = 6;
    md[1][0] = 2;
    md[0][1] = 3.5;
    md[1][1] = 1.5;
    md[0][2] = 2.5;
    md[1][2] = 0.5;
    SimTK_TEST(m[0][0] == 6);
    SimTK_TEST(m[1][0] == 2);
    SimTK_TEST(m[0][1] == 3.5);
    SimTK_TEST(m[1][1] == 1.5);
    SimTK_TEST(m[0][2] == 2.5);
    SimTK_TEST(m[1][2] == 0.5);

    /// multiplication
    Matrix_<adouble> m2 = b * m;
    Matrix_<adouble> m3 = m * b;
    Matrix_<adouble> m2_int = b_int * m;
    Matrix_<adouble> m3_int = m * b_int;
    Matrix_<adouble> m2_ad = b_ad * m;
    Matrix_<adouble> m3_ad = m * b_ad;
    Matrix_<adouble> m2_fl = b_fl * m;
    Matrix_<adouble> m3_fl = m * b_fl;
    SimTK_TEST(m2[0][0] == 21);
    SimTK_TEST(m2[1][0] == 7);
    SimTK_TEST(m2[0][1] == 12.25);
    SimTK_TEST(m2[1][1] == 5.25);
    SimTK_TEST(m2[0][2] == 8.75);
    SimTK_TEST(m2[1][2] == 1.75);
    SimTK_TEST(m2_int[0][0] == 18);
    SimTK_TEST(m2_int[1][0] == 6);
    SimTK_TEST(m2_int[0][1] == 10.5);
    SimTK_TEST(m2_int[1][1] == 4.5);
    SimTK_TEST(m2_int[0][2] == 7.5);
    SimTK_TEST(m2_int[1][2] == 1.5);
    SimTK_TEST(m2_ad[0][0] == 21);
    SimTK_TEST(m2_ad[1][0] == 7);
    SimTK_TEST(m2_ad[0][1] == 12.25);
    SimTK_TEST(m2_ad[1][1] == 5.25);
    SimTK_TEST(m2_ad[0][2] == 8.75);
    SimTK_TEST(m2_ad[1][2] == 1.75);
    SimTK_TEST(m2_fl[0][0] == 21);
    SimTK_TEST(m2_fl[1][0] == 7);
    SimTK_TEST(m2_fl[0][1] == 12.25);
    SimTK_TEST(m2_fl[1][1] == 5.25);
    SimTK_TEST(m2_fl[0][2] == 8.75);
    SimTK_TEST(m2_fl[1][2] == 1.75);
    SimTK_TEST(m3[0][0] == 21);
    SimTK_TEST(m3[1][0] == 7);
    SimTK_TEST(m3[0][1] == 12.25);
    SimTK_TEST(m3[1][1] == 5.25);
    SimTK_TEST(m3[0][2] == 8.75);
    SimTK_TEST(m3[1][2] == 1.75);
    SimTK_TEST(m3_int[0][0] == 18);
    SimTK_TEST(m3_int[1][0] == 6);
    SimTK_TEST(m3_int[0][1] == 10.5);
    SimTK_TEST(m3_int[1][1] == 4.5);
    SimTK_TEST(m3_int[0][2] == 7.5);
    SimTK_TEST(m3_int[1][2] == 1.5);
    SimTK_TEST(m3_ad[0][0] == 21);
    SimTK_TEST(m3_ad[1][0] == 7);
    SimTK_TEST(m3_ad[0][1] == 12.25);
    SimTK_TEST(m3_ad[1][1] == 5.25);
    SimTK_TEST(m3_ad[0][2] == 8.75);
    SimTK_TEST(m3_ad[1][2] == 1.75);
    SimTK_TEST(m3_fl[0][0] == 21);
    SimTK_TEST(m3_fl[1][0] == 7);
    SimTK_TEST(m3_fl[0][1] == 12.25);
    SimTK_TEST(m3_fl[1][1] == 5.25);
    SimTK_TEST(m3_fl[0][2] == 8.75);
    SimTK_TEST(m3_fl[1][2] == 1.75);

    // division
    double b3 = 0.5;
    int b3_int = 2;
    float b3_fl = 0.5;
    adouble b3_ad = 0.5;

    Matrix_<adouble> m4 = m / b3;
    Matrix_<adouble> m4_int = m / b3_int;
    Matrix_<adouble> m4_ad = m / b3_ad;
    Matrix_<adouble> m4_fl = m / b3_fl;
    SimTK_TEST(m4[0][0] == 12);
    SimTK_TEST(m4[1][0] == 4);
    SimTK_TEST(m4[0][1] == 7);
    SimTK_TEST(m4[1][1] == 3);
    SimTK_TEST(m4[0][2] == 5);
    SimTK_TEST(m4[1][2] == 1);
    SimTK_TEST(m4_int[0][0] == 3);
    SimTK_TEST(m4_int[1][0] == 1);
    SimTK_TEST(m4_int[0][1] == 1.75);
    SimTK_TEST(m4_int[1][1] == 0.75);
    SimTK_TEST(m4_int[0][2] == 1.25);
    SimTK_TEST(m4_int[1][2] == 0.25);
    SimTK_TEST(m4_ad[0][0] == 12);
    SimTK_TEST(m4_ad[1][0] == 4);
    SimTK_TEST(m4_ad[0][1] == 7);
    SimTK_TEST(m4_ad[1][1] == 3);
    SimTK_TEST(m4_ad[0][2] == 5);
    SimTK_TEST(m4_ad[1][2] == 1);
    SimTK_TEST(m4_fl[0][0] == 12);
    SimTK_TEST(m4_fl[1][0] == 4);
    SimTK_TEST(m4_fl[0][1] == 7);
    SimTK_TEST(m4_fl[1][1] == 3);
    SimTK_TEST(m4_fl[0][2] == 5);
    SimTK_TEST(m4_fl[1][2] == 1);

    /// addition
    //Matrix_<adouble> m5 = b + m;
    //Matrix_<adouble> m6 = m + b;
    //Matrix_<adouble> m5_int = b_int +m;
    //Matrix_<adouble> m6_int = m + b_int;
    //Matrix_<adouble> m5_ad = b_ad +m;
    //Matrix_<adouble> m6_ad = m + b_ad;
    //Matrix_<adouble> m5_fl = b_fl +m;
    //Matrix_<adouble> m6_fl = m + b_fl;
    //std::cout << m5 << std::endl;
    //Matrix_<double> m5d = b + md;
    //std::cout << m5d << std::endl;


    //SimTK_TEST(m5[0][0] == 9.5);
    //SimTK_TEST(m5[1][0] == 5.5);
    //SimTK_TEST(m5[0][1] == 7);
    //SimTK_TEST(m5[1][1] == 5);
    //SimTK_TEST(m5[0][2] == 6);
    //SimTK_TEST(m5[1][2] == 5);
    //SimTK_TEST(m5_int[0][0] == 9);
    //SimTK_TEST(m5_int[1][0] == 5);
    //SimTK_TEST(m5_int[0][1] == 6.5);
    //SimTK_TEST(m5_int[1][1] == 4.5);
    //SimTK_TEST(m5_int[0][2] == 5.5);
    //SimTK_TEST(m5_int[1][2] == 3.5);
    //SimTK_TEST(m5_ad[0][0] == 9.5);
    //SimTK_TEST(m5_ad[1][0] == 5.5);
    //SimTK_TEST(m5_ad[0][1] == 7);
    //SimTK_TEST(m5_ad[1][1] == 5);
    //SimTK_TEST(m5_ad[0][2] == 6);
    //SimTK_TEST(m5_ad[1][2] == 5);
    //SimTK_TEST(m5_fl[0][0] == 9.5);
    //SimTK_TEST(m5_fl[1][0] == 5.5);
    //SimTK_TEST(m5_fl[0][1] == 7);
    //SimTK_TEST(m5_fl[1][1] == 5);
    //SimTK_TEST(m5_fl[0][2] == 6);
    //SimTK_TEST(m5_fl[1][2] == 5);
    //SimTK_TEST(m6[0][0] == 9.5);
    //SimTK_TEST(m6[1][0] == 7);
    //SimTK_TEST(m6[0][1] == 12.25);
    //SimTK_TEST(m6[1][1] == 5.25);
    //SimTK_TEST(m6[0][2] == 8.75);
    //SimTK_TEST(m6[1][2] == 1.75);
    //SimTK_TEST(m6_int[0][0] == 9);
    //SimTK_TEST(m6_int[1][0] == 5);
    //SimTK_TEST(m6_int[0][1] == 6.5);
    //SimTK_TEST(m6_int[1][1] == 4.5);
    //SimTK_TEST(m6_int[0][2] == 5.5);
    //SimTK_TEST(m6_int[1][2] == 3.5);
    //SimTK_TEST(m6_ad[0][0] == 9.5);
    //SimTK_TEST(m6_ad[1][0] == 5.5);
    //SimTK_TEST(m6_ad[0][1] == 7);
    //SimTK_TEST(m6_ad[1][1] == 5);
    //SimTK_TEST(m6_ad[0][2] == 6);
    //SimTK_TEST(m6_ad[1][2] == 5);
    //SimTK_TEST(m6_fl[0][0] == 9.5);
    //SimTK_TEST(m6_fl[1][0] == 5.5);
    //SimTK_TEST(m6_fl[0][1] == 7);
    //SimTK_TEST(m6_fl[1][1] == 5);
    //SimTK_TEST(m6_fl[0][2] == 6);
    //SimTK_TEST(m6_fl[1][2] == 5);

    // RowVectors
    RowVector_<adouble> vr(2);
    vr[0] = 4;
    vr[1] = 6;
    SimTK_TEST(vr[0] == 4);
    SimTK_TEST(vr[1] == 6);
    /// multiplication
    RowVector_<adouble> vr2 = b * vr;
    RowVector_<adouble> vr3 = vr * b;
    RowVector_<adouble> vr2_int = b_int * vr;
    RowVector_<adouble> vr3_int = vr * b_int;
    RowVector_<adouble> vr2_ad = b_ad * vr;
    RowVector_<adouble> vr3_ad = vr * b_ad;
    RowVector_<adouble> vr2_fl = b_fl * vr;
    RowVector_<adouble> vr3_fl = vr * b_fl;
    /// Ensure that v is not modified simply by using it in an operation.
    SimTK_TEST(vr[0] == 4);
    SimTK_TEST(vr[1] == 6);
    /// Ensure that the correct operator* is used.
    SimTK_TEST(vr2[0] == 14);
    SimTK_TEST(vr2[1] == 21);
    SimTK_TEST(vr3[0] == 14);
    SimTK_TEST(vr3[1] == 21);
    SimTK_TEST(vr2_int[0] == 12);
    SimTK_TEST(vr2_int[1] == 18);
    SimTK_TEST(vr3_int[0] == 12);
    SimTK_TEST(vr3_int[1] == 18);
    SimTK_TEST(vr2_ad[0] == 14);
    SimTK_TEST(vr2_ad[1] == 21);
    SimTK_TEST(vr3_ad[0] == 14);
    SimTK_TEST(vr3_ad[1] == 21);
    SimTK_TEST(vr2_fl[0] == 14);
    SimTK_TEST(vr2_fl[1] == 21);
    SimTK_TEST(vr3_fl[0] == 14);
    SimTK_TEST(vr3_fl[1] == 21);
    /// division
    RowVector_<adouble> vr4 = vr / b3;
    RowVector_<adouble> vr4_int = vr / b3_int;
    RowVector_<adouble> vr4_ad = vr / b3_ad;
    RowVector_<adouble> vr4_fl = vr / b3_fl;
    /// Ensure that vr is not modified simply by using it in an operation.
    SimTK_TEST(vr[0] == 4);
    SimTK_TEST(vr[1] == 6);
    /// Ensure that the correct operator/ is used.
    SimTK_TEST(vr4[0] == 8);
    SimTK_TEST(vr4[1] == 12);
    SimTK_TEST(vr4_int[0] == 2);
    SimTK_TEST(vr4_int[1] == 3);
    SimTK_TEST(vr4_ad[0] == 8);
    SimTK_TEST(vr4_ad[1] == 12);
    SimTK_TEST(vr4_fl[0] == 8);
    SimTK_TEST(vr4_fl[1] == 12);
    /// addition
    RowVector_<adouble> vr6 = vr + b3;
    RowVector_<adouble> vr5 = b3 + vr;
    RowVector_<adouble> vr6_int = vr + b3_int;
    RowVector_<adouble> vr5_int = b3_int + vr;
    RowVector_<adouble> vr6_ad = vr + b3_ad;
    RowVector_<adouble> vr5_ad = b3_ad + vr;
    RowVector_<adouble> vr6_fl = vr + b3_fl;
    RowVector_<adouble> vr5_fl = b3_fl + vr;
    /// Ensure that vr is not modified simply by using it in an operation.
    SimTK_TEST(vr[0] == 4);
    SimTK_TEST(vr[1] == 6);
    /// Ensure that the correct operator/ is used.
    SimTK_TEST(vr6[0] == 4.5);
    SimTK_TEST(vr6[1] == 6.5);
    SimTK_TEST(vr6_int[0] == 6);
    SimTK_TEST(vr6_int[1] == 8);
    SimTK_TEST(vr6_ad[0] == 4.5);
    SimTK_TEST(vr6_ad[1] == 6.5);
    SimTK_TEST(vr6_fl[0] == 4.5);
    SimTK_TEST(vr6_fl[1] == 6.5);
    SimTK_TEST(vr5[0] == 4.5);
    SimTK_TEST(vr5[1] == 6.5);
    SimTK_TEST(vr5_int[0] == 6);
    SimTK_TEST(vr5_int[1] == 8);
    SimTK_TEST(vr5_ad[0] == 4.5);
    SimTK_TEST(vr5_ad[1] == 6.5);
    SimTK_TEST(vr5_fl[0] == 4.5);
    SimTK_TEST(vr5_fl[1] == 6.5);
    /// substraction
    RowVector_<adouble> vr7 = vr - b3;
    RowVector_<adouble> vr8 = b3 - vr;
    RowVector_<adouble> vr7_int = vr - b3_int;
    RowVector_<adouble> vr8_int = b3_int - vr;
    RowVector_<adouble> vr7_ad = vr - b3_ad;
    RowVector_<adouble> vr8_ad = b3_ad - vr;
    RowVector_<adouble> vr7_fl = vr - b3_fl;
    RowVector_<adouble> vr8_fl = b3_fl - vr;
    /// Ensure that vr is not modified simply by using it in an operation.
    SimTK_TEST(vr[0] == 4);
    SimTK_TEST(vr[1] == 6);
    /// Ensure that the correct operator/ is used.
    SimTK_TEST(vr7[0] == 3.5);
    SimTK_TEST(vr7[1] == 5.5);
    SimTK_TEST(vr7_int[0] == 2);
    SimTK_TEST(vr7_int[1] == 4);
    SimTK_TEST(vr7_ad[0] == 3.5);
    SimTK_TEST(vr7_ad[1] == 5.5);
    SimTK_TEST(vr7_fl[0] == 3.5);
    SimTK_TEST(vr7_fl[1] == 5.5);
    SimTK_TEST(vr8[0] == -3.5);
    SimTK_TEST(vr8[1] == -5.5);
    SimTK_TEST(vr8_int[0] == -2);
    SimTK_TEST(vr8_int[1] == -4);
    SimTK_TEST(vr8_ad[0] == -3.5);
    SimTK_TEST(vr8_ad[1] == -5.5);
    SimTK_TEST(vr8_fl[0] == -3.5);
    SimTK_TEST(vr8_fl[1] == -5.5);
}

void testAdouble() {
    double* xp = new double[1];
    xp[0] = -2.3;
    adouble* x = new adouble[1];
    trace_on(1);
    x[0] <<= xp[0];
    adouble y[1];
    y[0] = 3 * pow(x[0], 3) + cos(x[0]) + 1;
    double y0[1];
    y[0] >>= y0[0];
    trace_off();

    double** J;
    J = myalloc(1, 1);
    jacobian(1, 1, 1, xp, J);
    SimTK_TEST(J[0][0] == 48.355705212176716);
}

int main() {
    SimTK_START_TEST("TestAutodiffADOLCCommon");
        SimTK_SUBTEST(testVectorOperations);
        SimTK_SUBTEST(testAdouble);
    SimTK_END_TEST();
}

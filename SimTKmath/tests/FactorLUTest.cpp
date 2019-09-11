/* -------------------------------------------------------------------------- *
 *                          Simbody(tm): SimTKmath                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-12 Stanford University and the Authors.        *
 * Authors: Jack Middleton                                                    *
 * Contributors:                                                              *
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

/**@file
 * This is a test program which uses the FactorLU  class to do an LU 
 * factorization on a system of linear equations and then use the 
 * factored LU matrix to solve for a particular right hand side 
 */

/*
The data for this test is from an example FORTRAN  program from the
Numerical Algorithms Group (NAG)
URL:http://www.nag.com/lapack-ex/lapack-ex.html


Solves:

Ax = B,

where A is the general matrix


      1.80   2.88   2.05   -0.89           9.52
A =   5.25  -2.95  -0.95   -3.80   and B = 24.35
      1.58  -2.69  -2.90   -1.04           0.77
     -1.11  -0.66  -0.59    0.80          -6.22


   Solution
 x =     1.0000    -1.0000     3.0000    -5.0000


 LU factorization:
             1          2          3          4
 1      5.2500    -2.9500    -0.9500    -3.8000
 2      0.3429     3.8914     2.3757     0.4129
 3      0.3010    -0.4631    -1.5139     0.2948
 4     -0.2114    -0.3299     0.0047     0.1314

 Pivot indices
             2          2          3          4

*/

#include "SimTKmath.h"

#include <iostream>
using std::cout; using std::endl;

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

using namespace SimTK;

#ifndef SimTK_REAL_IS_ADOUBLE
	Real A[16] = { 1.80,   2.88,   2.05,   -0.89,
				   5.25,  -2.95,  -0.95,   -3.80,
				   1.58,  -2.69,  -2.90,   -1.04,
				  -1.11,  -0.66,  -0.59,    0.80  };

	Real B[4] =  { 9.52, 24.35,  0.77, -6.22 };

	Real X[4] =  { 1.,   -1.,    3.,   -5.   };
#else
	double A[16] = { 1.80,   2.88,   2.05,   -0.89,
	5.25,  -2.95,  -0.95,   -3.80,
	1.58,  -2.69,  -2.90,   -1.04,
	-1.11,  -0.66,  -0.59,    0.80 };

	double B[4] = { 9.52, 24.35,  0.77, -6.22 };

	double X[4] = { 1.,   -1.,    3.,   -5. };
#endif

int main () {
    try { 
            // Default precision (Real, normally double) test.
		#ifndef SimTK_REAL_IS_ADOUBLE
			Matrix a(4,4, A);
			Vector b(4, B);
			Vector x_right(4, X);
			Vector x; // should get sized automatically to 4 by solve()
		#else
			Matrix_<double> a(4, 4, A);
			Vector_<double> b(4, B);
			Vector_<double> x_right(4, X);
			Vector_<double> x; // should get sized automatically to 4 by solve()
		#endif

        FactorLU lu(a);  // perform LU factorization 

        lu.solve( b, x );  // solve for x given a right hand side 

        cout << " Real SOLUTION: " << x << "  errnorm=" << (x-x_right).norm() << endl;
        ASSERT((x-x_right).norm() < 10*SignificantReal);

            // float test
		
        Matrix_<float> af(4,4); for (int i=0; i<4; ++i) for (int j=0; j<4; ++j) af(i,j)=(float)a(i,j);
        Vector_<float> bf(4); for (int i=0; i<4; ++i) bf[i] = (float)b[i];
        Vector_<float> xf_right(4); for (int i=0; i<4; ++i) xf_right[i] = (float)x_right[i];
        Vector_<float> xf; // should get sized automatically to 4 by solve()

        FactorLU luf;
        luf.factor(af);
        luf.solve(bf, xf);

        cout << " float SOLUTION: " << xf << "  errnorm=" << (xf-xf_right).norm() << endl;
        const float SignificantFloat = NTraits<float>::getSignificant();
        ASSERT((xf-xf_right).norm() < 10*SignificantFloat);

        luf.factor(a);
        lu.solve( b, x );  // solve for x given a right hand side 
        cout << " Real SOLUTION: " << x << "  errnorm=" << (x-x_right).norm() << endl;
        ASSERT((x-x_right).norm() < 10*SignificantReal);
        
		#ifndef SimTK_REAL_IS_ADOUBLE
			Real C[4] = { 1.0,   2.0,
						  1.0,   3.0  };
			Matrix c(2,2, C);
			FactorLU clu(c);
			Matrix invC;
			clu.inverse(invC);
		#else
				double C[4] = { 1.0,   2.0,
					1.0,   3.0 };
				Matrix_<double> c(2, 2, C);
				FactorLU clu(c);
				Matrix_<double> invC;
				clu.inverse(invC);
		#endif
        
		cout << "Inverse c: " << endl;
        cout << invC[0] << endl;
        cout << invC[1] << endl;
        
		#ifndef SimTK_REAL_IS_ADOUBLE
			Real Z[4] = { 0.0,   0.0,
						 0.0,   0.0  };
			Matrix z(2,2, Z);
			FactorLU zlu(z);
			Vector_<SimTK::Real> xz;
			Vector_<SimTK::Real> bz(2);
		#else
			double Z[4] = { 0.0,   0.0,
				0.0,   0.0 };
			Matrix_<double> z(2, 2, Z);
			FactorLU zlu(z);
			Vector_<double> xz;
			Vector_<double> bz(2);
		#endif
        bz(1) = bz(0) = 0.0;
        zlu.solve( bz, xz );
        cout << " solve with mat all zeros : " << endl;
        for(int i=0;i<xz.size();i++) printf("%f ", xz(i) );  printf("\n");
   
        try {
			#ifndef SimTK_REAL_IS_ADOUBLE
				Matrix_<SimTK::Real> z0;
				FactorLU z0lu(z0);
				Vector_<SimTK::Real> bz0(0);
			#else
				Matrix_<double> z0;
				FactorLU z0lu(z0);
				Vector_<double> bz0(0);
			#endif
            z0lu.solve( bz0, xz );
            cout << " solve with mat(0,0) : " << endl;
            for(int i=0;i<xz.size();i++) printf("%f ", xz(i) );  printf("\n");
        } catch (const std::exception& e) {
             cout << "(EXPECTED EXCEPTION) NULL matrix test: " 
                 << e.what() << endl;
        }
    } 
    catch (const std::exception& e) {
        std::printf("FAILED: %s\n", e.what());
        return 1;
    }

    return 0;
}



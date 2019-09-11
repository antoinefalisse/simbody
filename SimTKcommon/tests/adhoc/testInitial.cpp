/* -------------------------------------------------------------------------- *
 *                       SimTK Simbody: SimTKcommon                           *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
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

//#include <simbody.h>
#include <SimTKcommon.h>
#include <adolc.h>
#include <adolc_sparse.h>
//
#include <cstdlib> // for rand()
#include <ctime>
#include <string>
#include <complex>
#include <iostream>
#include <iterator>
#include <random>
#include <cassert>
#include <algorithm>
////
////
////
//////
using std::cout;
using std::endl;
using std::complex;

using namespace SimTK;

void print_and_clear(int n, int m, double* const * J) {
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			std::cout << J[i][j] << " ";
			// Clear for the next method.
			J[i][j] = 0;
		}
		std::cout << std::endl;
	}
}
//
///// This is a dense function R^n -> R^m. It is templatized so that we can
///// evaluate it on both double and adouble.
void constraint_function_dense(int n, int m, SimTK::Vector_<SimTK::Real> x, SimTK::Vector_<SimTK::Real> y) {
	for (int j = 0; j < m; ++j) {
		y[j] = 0; // Clean up the given memory.
		for (int i = 0; i < n; ++i) {
			y[j] += log(j + 2) * cos(x[i]) * exp(-x[(n - 1) - i]);
		}
	}
}

int main() {
	int n = 10; int m = 10000; 
	double** J=myalloc(m,n);
	std::vector<double> px(n,1);
	

	trace_on(0);
	SimTK::Real a;
	typedef SimTK::Vector_<SimTK::Real> vector;
	vector x(n), y(m);

	SimTK::Vector_<double> py(m);
	const double* px_i = px.data();

	// Indicate independent variables.
	for (int i = 0; i < n; ++i) {
		x[i] <<= px[i];
	}
	//// Evaluate function.
	constraint_function_dense(n, m, x, y);

	//// Indicate dependent variables.
	for (int j = 0; j < m; ++j) y[j] >>= py[j];

	//// Stop recording.
	trace_off();

	// Use the recorded tape to compute the jacobian.
	using namespace std::chrono;
	auto start = steady_clock::now();

	
	int success = jacobian(0, m, n, px_i, J);

		
	auto end = steady_clock::now();

	/*std::cout << "J(x):" << std::endl;
	print_and_clear(x.size(), y.size(), J);*/
	
	auto dur = duration_cast<milliseconds>(end - start).count();
	std::cout << "Duration: " << dur << " milliseconds" << std::endl;
	return 0;

}

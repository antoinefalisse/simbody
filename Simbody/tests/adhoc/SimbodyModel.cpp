#include "Simbody.h"
#include "SimTKcommon.h"
#include <adolc.h>
#include <adolc_sparse.h>
#include "SimTKmath.h"

#include <iostream>
#include <iterator>
#include <random>
#include <cassert>
#include <algorithm>
#include <fstream>

using namespace SimTK;

int main() {

	// Define the system.
	MultibodySystem* system = new MultibodySystem();
	SimbodyMatterSubsystem* matter = new SimbodyMatterSubsystem(*system);
	GeneralForceSubsystem* forces = new GeneralForceSubsystem(*system);
	Force::UniformGravity* gravity = new Force::UniformGravity(*forces, *matter, Vec3(0, Real(-9.81), 0));

	// Describe mass and visualization properties for a generic body.
	Body::Rigid* bodyInfo = new Body::Rigid(MassProperties(1.0, Vec3(0, 0, 0), UnitInertia(0)));

	// Create the moving (mobilized) bodies of the pendulum.
	MobilizedBody::Pin* pendulum1 = new MobilizedBody::Pin(matter->Ground(), Transform(Vec3(0)),
		*bodyInfo, Transform(Vec3(0, 1, 0)));
	MobilizedBody::Pin* pendulum2 = new MobilizedBody::Pin(*pendulum1, Transform(Vec3(0)),
		*bodyInfo, Transform(Vec3(0, 1, 0)));

	MobilizedBodyIndex* indxpd1 = new MobilizedBodyIndex(pendulum1->getMobilizedBodyIndex());
	MobilizedBodyIndex* indxpd2 = new MobilizedBodyIndex(pendulum2->getMobilizedBodyIndex());

	// Initialize the system and state.
	State* Modelstate = new State(system->realizeTopology());
	adouble* p = new adouble[1];

	/// Option1	
	//double xp[5];
	//xp[0] = 2;
	//xp[1] = 4;
	//xp[2] = 3;
	//xp[3] = 5;
	//xp[4] = 0;

	//// Tape start
	//trace_on(1, 1, 8000000, 8000000, 8000000, 8000000, 1);

	//Vector& Q = Modelstate->updQ();
	//Vector& U = Modelstate->updU();

	//// Mark independent variables
	//Q[0] <<= xp[0];
	//U[0] <<= xp[1];
	//Q[1] <<= xp[2];
	//U[1] <<= xp[3];
	//p[0] <<= xp[4];
	/// end Option1

	/// Option2, using new and delete
	double* xp = new double[5];
	xp[0] = 2;
	xp[1] = 4;
	xp[2] = 3;
	xp[3] = 5;
	xp[4] = 0;

	trace_on(1, 1, 8000000, 8000000, 8000000, 8000000, 1);

	adouble* Q = new adouble[2];
	adouble* U = new adouble[2];

	Q[0] <<= xp[0];
	U[0] <<= xp[1];
	Q[1] <<= xp[2];
	U[1] <<= xp[3];
	p[0] <<= xp[4];

	pendulum1->setOneQ(*Modelstate, 0, Q[0]);
	pendulum1->setOneU(*Modelstate, 0, U[0]);
	pendulum2->setOneQ(*Modelstate, 0, Q[1]);
	pendulum2->setOneU(*Modelstate, 0, U[1]);
	/// end option2

	system->realize(*Modelstate, Stage::Velocity);

	// Computation residual forces
	Vector appliedMobilityForces(2);
	appliedMobilityForces.setToZero();
	double grav = 9.81;
	double m1 = pendulum1->getBodyMass(*Modelstate).getValue();
	double m2 = pendulum2->getBodyMass(*Modelstate).getValue();

	/// Length bars = length to center of mass
	Transform Trl1 = pendulum1->getDefaultOutboardFrame();
	Vec3 Vecl1 = Trl1.updP();
	double l1 = Vecl1.get(1).getValue();
	Transform Trl2 = pendulum2->getDefaultOutboardFrame();
	Vec3 Vecl2 = Trl2.updP();
	double l2 = Vecl2.get(1).getValue();

	appliedMobilityForces[0] = p[0] * grav * (l2*m2*cos(Q[0] + Q[1]) + l1*m2*cos(Pi - Q[0]) + l1*m1*cos(Pi - Q[0]));
	appliedMobilityForces[1] = p[0] * grav * (l2*m2*cos(Q[0] + Q[1]));

	Vector_<SpatialVec> appliedBodyForces;
	appliedBodyForces.resize(3);
	appliedBodyForces.setToZero();
	Vector knownUdot(2);
	knownUdot.setToZero();
	Vector residualMobilityForces(2);
	residualMobilityForces.setToZero();
	Vec3 g = gravity->getGravity();

	matter->addInStationForce(*Modelstate, *indxpd1, pendulum1->getBodyMassCenterStation(*Modelstate), pendulum1->getBodyMass(*Modelstate)*g, appliedBodyForces);
	matter->addInStationForce(*Modelstate, *indxpd2, pendulum2->getBodyMassCenterStation(*Modelstate), pendulum2->getBodyMass(*Modelstate)*g, appliedBodyForces);

	matter->calcResidualForceIgnoringConstraints(*Modelstate, appliedMobilityForces, appliedBodyForces, knownUdot, residualMobilityForces);

	printf("Nominal output \n %f %f \n \n",
		residualMobilityForces[0].getValue(), residualMobilityForces[1].getValue());

	adouble y[2];
	y[0] = residualMobilityForces[0];
	y[1] = residualMobilityForces[1];

	double residualMobilityForces_f0;
	double residualMobilityForces_f1;

	y[0] >>= residualMobilityForces_f0;
	y[1] >>= residualMobilityForces_f1;

	/// Option2
	delete[] Q;
	delete[] U;
	delete[] p;

	trace_off();

	size_t tape_stats[STAT_SIZE];
	tapestats(1, tape_stats);

	std::ofstream myfile1;
	myfile1.open("test1.txt");
	myfile1 << tape_stats[0]; myfile1 << "  ";
	myfile1 << tape_stats[1]; myfile1 << "  ";
	myfile1 << tape_stats[2]; myfile1 << "  ";
	myfile1 << tape_stats[3]; myfile1 << "  ";
	myfile1 << tape_stats[4]; myfile1 << "  ";
	myfile1 << tape_stats[5]; myfile1 << "  ";
	myfile1 << tape_stats[6]; myfile1 << "  ";
	myfile1 << tape_stats[7]; myfile1 << "  ";
	myfile1 << tape_stats[8]; myfile1 << "  ";
	myfile1 << tape_stats[9]; myfile1 << "  ";
	myfile1 << tape_stats[10]; myfile1 << "  ";
	myfile1 << tape_stats[11]; myfile1 << "  ";
	myfile1 << tape_stats[12]; myfile1 << "  ";
	myfile1 << tape_stats[13]; myfile1 << "  ";
	myfile1 << tape_stats[14]; myfile1 << "  ";
	myfile1 << tape_stats[15]; myfile1 << "  ";
	myfile1 << tape_stats[16]; myfile1 << "  ";
	myfile1 << tape_stats[17]; myfile1 << "  ";
	myfile1.close();

	/*size_t tape_stats[STAT_SIZE];
	tapestats(1, tape_stats);
	cout << "number of independent variables " << tape_stats[0] << "\n";
	cout << "number of dependent variables " << tape_stats[1] << "\n";
	cout << "max number of live active variables " << tape_stats[2] << "\n";
	cout << "size taylor stack " << tape_stats[3] << "\n";
	cout << "buffer size " << tape_stats[4] << "\n";
	cout << "total number of operations recorded " << tape_stats[5] << "\n";*/

	double** J;
	J = myalloc(2, 5);
	jacobian(1, 2, 5, xp, J);
	printf("Jacobian \n %f %f %f %f %f \n  %f %f %f %f %f \n \n",
		J[0][0], J[0][1], J[0][2], J[0][3], J[0][4],
		J[1][0], J[1][1], J[1][2], J[1][3], J[1][4]);

	double xp_dot[5];
	xp_dot[0] = 3;
	xp_dot[1] = 5;
	xp_dot[2] = 8;
	xp_dot[3] = 9;
	xp_dot[4] = 0;

	double yy[2];
	yy[0] = y[0].getValue();
	yy[1] = y[1].getValue();

	double y_dot[2];
	fos_forward(1, 2, 5, 0, xp, xp_dot, yy, y_dot);
	printf("Forward sensitivities \n %f %f \n \n",
		y_dot[0], y_dot[1]);

	double yy_bar[2];
	yy_bar[0] = 2;
	yy_bar[1] = 3;

	double xp_bar[5];
	fos_reverse(1, 2, 5, yy_bar, xp_bar);
	printf("Reverse sensitivities \n %f %f %f %f %f \n \n",
		xp_bar[0], xp_bar[1], xp_bar[2], xp_bar[3], xp_bar[4]);

	Visualizer viz(*system);

	viz.report(*Modelstate);
}
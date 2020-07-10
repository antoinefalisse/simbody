#include "Simbody.h"

#include <iostream>
#include <iterator>
#include <random>
#include <cassert>
#include <algorithm>
#include <vector>
#include <fstream>

using namespace SimTK;

constexpr int n_in = 3;
constexpr int n_out = 1;

constexpr int NX = 4; // states
constexpr int NU = 2; // controls
constexpr int NP = 1; // perturbation
constexpr int NR = 2; // residual forces

// Function F
int F_generic(const double** arg, double** res) {

    MultibodySystem* systemmodel;
    SimbodyMatterSubsystem* matter;
    GeneralForceSubsystem* forces;
    Force::UniformGravity* gravity;
    Body::Rigid* bodyInfo1;
    Body::Rigid* bodyInfo2;
    MobilizedBody::Pin* pendulum1;
    MobilizedBody::Pin* pendulum2;
    MobilizedBodyIndex* indxpd1;
    MobilizedBodyIndex* indxpd2;
    State* Modelstate;

    double grav = 9.81;
    double m1;
    double m2;
    Transform Trl1;
    Vec3 Vecl1;
    double l1;
    Transform Trl2;
    Vec3 Vecl2;
    double l2;
    Vec3 g;

    // Define the system.
	systemmodel = new MultibodySystem();
	matter = new SimbodyMatterSubsystem(*systemmodel);
	forces = new GeneralForceSubsystem(*systemmodel);
	gravity = new Force::UniformGravity(*forces, *matter, Vec3(0, Real(-9.81), 0));

	// Describe mass properties for a generic body.
	bodyInfo1 = new Body::Rigid(MassProperties(26.0, Vec3(0, 0.55, 0), UnitInertia(1.4 / 26 + 0.55*0.55))); // leg 1.4
	bodyInfo2 = new Body::Rigid(MassProperties(46.0, Vec3(0, 0.365, 0), UnitInertia(2.9 / 46 + 0.365*0.365))); // torso 2.9;

	// Create the moving (mobilized) bodies of the pendulum.
	pendulum1 = new MobilizedBody::Pin(matter->Ground(), Transform(Vec3(0)),
		*bodyInfo1, Transform(Vec3(0)));
	pendulum2 = new MobilizedBody::Pin(*pendulum1, Transform(Vec3(0, 0.853, 0)),
		*bodyInfo2, Transform(Vec3(0)));

	indxpd1 = new MobilizedBodyIndex(pendulum1->getMobilizedBodyIndex());
	indxpd2 = new MobilizedBodyIndex(pendulum2->getMobilizedBodyIndex());

	// Initialize the system and state.
	Modelstate = new State(systemmodel->realizeTopology());

    #ifndef SimTK_REAL_IS_ADOUBLE
	    m1 = pendulum1->getBodyMass(*Modelstate);
	    m2 = pendulum2->getBodyMass(*Modelstate);
	    Trl1 = pendulum2->getDefaultInboardFrame();
	    Vecl1 = Trl1.updP();
	    l1 = Vecl1.get(1);
    #else
	    m1 = pendulum1->getBodyMass(*Modelstate).getValue();
	    m2 = pendulum2->getBodyMass(*Modelstate).getValue();
	    Trl1 = pendulum2->getDefaultInboardFrame();
	    Vecl1 = Trl1.updP();
	    l1 = Vecl1.get(1).getValue();
    #endif

	g = gravity->getGravity();

    // Read inputs
	std::vector<double> x(arg[0], arg[0] + NX);
	std::vector<double> u(arg[1], arg[1] + NU);

    static double pert[NP];
	static double ua[NU];

	Vector& Q = Modelstate->updQ();
	Vector& U = Modelstate->updU();

	Q[0] = x[0];
	U[0] = x[1];
	Q[1] = x[2];
	U[1] = x[3];
	ua[0] = u[0];
	ua[1] = u[1];

    systemmodel->realize(*Modelstate, Stage::Acceleration);

	// Computation residual forces
	/// Applied mobility forces
	static Vector appliedMobilityForces(2);
	appliedMobilityForces.setToZero();

	Vec3 r1_v = pendulum1->getBodyMassCenterStation(*Modelstate);
	Vec3 r2_v = pendulum2->getBodyMassCenterStation(*Modelstate);
    #ifndef SimTK_REAL_IS_ADOUBLE
        double r1 = r1_v[1];
	    double r2 = r2_v[1];
    #else
	    double r1 = r1_v[1].getValue();
	    double r2 = r2_v[1].getValue();
    #endif
	appliedMobilityForces[0] = grav * (r2*m2*cos(Q[0] + Q[1]) + r1*m1*cos(Q[0]) + m2*l1*cos(Q[0]));
	appliedMobilityForces[1] = grav * (r2*m2*cos(Q[0] + Q[1]));
	/// Applied mobility forces
	static Vector_<SpatialVec> appliedBodyForces;
	appliedBodyForces.resize(3);
	appliedBodyForces.setToZero();
	static Vector knownUdot(2);
	knownUdot[0] = ua[0];
	knownUdot[1] = ua[1];
	static Vector residualMobilityForces(2);
	residualMobilityForces.setToZero();
	matter->addInStationForce(*Modelstate, *indxpd1, pendulum1->getBodyMassCenterStation(*Modelstate), pendulum1->getBodyMass(*Modelstate)*g, appliedBodyForces);
	matter->addInStationForce(*Modelstate, *indxpd2, pendulum2->getBodyMassCenterStation(*Modelstate), pendulum2->getBodyMass(*Modelstate)*g, appliedBodyForces);
	/// Residual forces
	matter->calcResidualForceIgnoringConstraints(*Modelstate, appliedMobilityForces, appliedBodyForces, knownUdot, residualMobilityForces);

    Vec3 test = pendulum1->getBodyAngularVelocity(*Modelstate);

    SpatialVec temp = pendulum1->getBodyAcceleration(*Modelstate);




    #ifndef SimTK_REAL_IS_ADOUBLE
        if (res[0]) for (int i = 0; i < NR; ++i) res[0][i] = (residualMobilityForces[i]);
    #else
        if (res[0]) for (int i = 0; i < NR; ++i) res[0][i] = (residualMobilityForces[i].getValue());
    #endif

	return 0;


}

int main() {

    double x[NX];
    double u[NU];
    double tau[NR];

    for (int i = 0; i < NX; ++i) x[i] = 1;
    for (int i = 0; i < NU; ++i) u[i] = 1;

    const double* Recorder_arg[n_in] = { x,u };
    double* Recorder_res[n_out] = { tau };

    F_generic(Recorder_arg, Recorder_res);

    double res[NR];
    for (int i = 0; i < NR; ++i) {
        Recorder_res[0][i];
        std::cout << Recorder_res[0][i] << std::endl;
    }

    return 0;

}
/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
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

#include "SimTKsimbody.h"

using namespace SimTK;
using namespace std;

static int is_initialised = false;
static int is_initialised_Fadolc = false;

/// Declare static Simbody objects
static MultibodySystem* systemmodel;
static SimbodyMatterSubsystem* matter;
static GeneralForceSubsystem* forces;
static Force::UniformGravity* gravity;
static HuntCrossleyForce_smooth* hc_heel;
static HuntCrossleyForce_smooth* hc_front;
static Body::Rigid* bodyInfo0;
static Body::Rigid* bodyInfo1;
static Body::Rigid* bodyInfo2;
static MobilizedBody::Planar* footground;
static MobilizedBody::Pin* pendulum1;
static MobilizedBody::Pin* pendulum2;
static MobilizedBodyIndex* indxpd0;
static MobilizedBodyIndex* indxpd1;
static MobilizedBodyIndex* indxpd2;
static State* Modelstate;

static double grav = 9.81;
static double m0;
static double m1;
static double m2;
static Transform Trl1;
static Vec3 Vecl1;
static double l1;
static Transform Trl2;
static Vec3 Vecl2;
static double l2;
static Vec3 g;

void model_setup() {

	// Define the system.
	systemmodel = new MultibodySystem();
	matter = new SimbodyMatterSubsystem(*systemmodel);
	forces = new GeneralForceSubsystem(*systemmodel);
	gravity = new Force::UniformGravity(*forces, *matter, Vec3(0, Real(-9.81), 0));
	hc_heel = new HuntCrossleyForce_smooth(*forces);
	hc_front = new HuntCrossleyForce_smooth(*forces);

	// Describe mass properties for a generic body.
	bodyInfo0 = new Body::Rigid(MassProperties(2.0, Vec3(0.05, 0, 0), UnitInertia(0.01 / 2.0 + 0.05*0.05))); // foot
	bodyInfo1 = new Body::Rigid(MassProperties(24.0, Vec3(0, 0.54, 0), UnitInertia(1.4 / 24 + 0.54*0.54))); // leg 1.4
	bodyInfo2 = new Body::Rigid(MassProperties(46.0, Vec3(0, 0.365, 0), UnitInertia(2.9 / 46 + 0.365*0.365))); // torso 2.9;
																								   //pendulum1 = new MobilizedBody::Pin(matter->Ground(), Transform(Vec3(0)),
																											   //	*bodyInfo1, Transform(Vec3(0)));
	footground = new MobilizedBody::Planar(matter->Ground(), Transform(Vec3(0)),
		*bodyInfo0, Transform(Vec3(0)));
	

	pendulum1 = new MobilizedBody::Pin(*footground, Transform(Vec3(0)),
		*bodyInfo1, Transform(Vec3(0)));
	pendulum2 = new MobilizedBody::Pin(*pendulum1, Transform(Vec3(0, 0.853, 0)),
		*bodyInfo2, Transform(Vec3(0)));

	indxpd0 = new MobilizedBodyIndex(footground->getMobilizedBodyIndex());
	indxpd1 = new MobilizedBodyIndex(pendulum1->getMobilizedBodyIndex());
	indxpd2 = new MobilizedBodyIndex(pendulum2->getMobilizedBodyIndex());

	// Initialize the system and state.
	Modelstate = new State(systemmodel->realizeTopology());

	m0 = footground->getBodyMass(*Modelstate).getValue();
	m1 = pendulum1->getBodyMass(*Modelstate).getValue();
	m2 = pendulum2->getBodyMass(*Modelstate).getValue();
	Trl1 = pendulum2->getDefaultInboardFrame();
	Vecl1 = Trl1.updP();
	l1 = Vecl1.get(1).getValue();

	g = gravity->getGravity();

	Real static_friction = 0.8;
	Real dynamic_friction = 0.8;
	Real viscous_friction = 0.5;
	Real transition_velocity = 0.1;
	Real stiffness = 3777162;
	Real dissipation = 2;
	Real loc_heel_x = -0.05;
	Real radius_heel = 0.02;
	Real loc_front_x = 0.15;
	Real radius_front = 0.02;

	hc_heel = new HuntCrossleyForce_smooth(*forces);
	hc_heel->setParameters(stiffness, dissipation, static_friction, dynamic_friction, viscous_friction, transition_velocity);
	hc_heel->setLocSphere(Vec3(loc_heel_x, 0, 0));
	hc_heel->setBodySphere(matter->getMobilizedBody(*indxpd0));
	hc_heel->setRadiusSphere(radius_heel);

	hc_front = new HuntCrossleyForce_smooth(*forces);
	hc_front->setParameters(stiffness, dissipation, static_friction, dynamic_friction, viscous_friction, transition_velocity);
	hc_front->setLocSphere(Vec3(loc_front_x, 0, 0));
	hc_front->setBodySphere(matter->getMobilizedBody(*indxpd0));
	hc_front->setRadiusSphere(radius_front);
}

void main() {
	if (!is_initialised) {
		model_setup();
		is_initialised = true;
	}

	static Real pert[1];
	static Real ua[5];

	Vector& Q = Modelstate->updQ();
	Vector& U = Modelstate->updU();

	Q[0] = 1;
	U[0] = 0;
	Q[1] = 1;
	U[1] = 0;
	Q[2] = 1;
	U[2] = 0;
	Q[3] = 1;
	U[3] = 0;
	Q[4] = 1;
	U[4] = 0;
	ua[0] = 0;
	ua[1] = 0;
	ua[2] = 0;
	ua[3] = 0;
	ua[4] = 0;
	pert[0] = 0;
	std::cout << "Q=" << Q << " U=" << U << std::endl;
	//Transform bodyT=matter->getMobilizedBody(*indxpd0).getBodyTransform(*Modelstate);

	hc_heel->enable(*Modelstate);
	hc_front->enable(*Modelstate);
	systemmodel->realize(*Modelstate, Stage::Velocity);

	// Computation residual forces
	/// Applied mobility forces
	static Vector appliedMobilityForces(5);
	appliedMobilityForces.setToZero();

	Vec3 r0_v = footground->getBodyMassCenterStation(*Modelstate);
	Vec3 r1_v = pendulum1->getBodyMassCenterStation(*Modelstate);
	Vec3 r2_v = pendulum2->getBodyMassCenterStation(*Modelstate);
	double r0 = r0_v[1].getValue();
	double r1 = r1_v[1].getValue();
	double r2 = r2_v[1].getValue();
	appliedMobilityForces[0] = pert[0] * grav*(r0*m0*sin(Q[0]) + r1*m1*cos(Q[0] + Q[1]) + m2*(l1*cos(Q[0] + Q[1]) + r2*cos(Q[0] + Q[1] + Q[2])));
	appliedMobilityForces[1] = -pert[0] * grav*(m0 + m1 + m2);
	appliedMobilityForces[2] = 0;
	appliedMobilityForces[3] = pert[0] * grav * (r1*m1*cos(Q[0] + Q[1]) + m2*(l1*cos(Q[0] + Q[1]) + r2*cos(Q[0] + Q[1] + Q[2])));
	appliedMobilityForces[4] = pert[0] * grav * (r2*m2*cos(Q[0] + Q[1] + Q[2]));

	/// Applied body forces
	static Vector_<SpatialVec> appliedBodyForces;
	appliedBodyForces.resize(4);
	appliedBodyForces.setToZero();
	static Vector knownUdot(5);
	knownUdot[0] = ua[0];
	knownUdot[1] = ua[1];
	knownUdot[2] = ua[2];
	knownUdot[3] = ua[3];
	knownUdot[4] = ua[4];
	static Vector residualMobilityForces(5);
	residualMobilityForces.setToZero();

	matter->addInStationForce(*Modelstate, *indxpd0, footground->getBodyMassCenterStation(*Modelstate), footground->getBodyMass(*Modelstate)*g, appliedBodyForces);
	matter->addInStationForce(*Modelstate, *indxpd1, pendulum1->getBodyMassCenterStation(*Modelstate), pendulum1->getBodyMass(*Modelstate)*g, appliedBodyForces);
	matter->addInStationForce(*Modelstate, *indxpd2, pendulum2->getBodyMassCenterStation(*Modelstate), pendulum2->getBodyMass(*Modelstate)*g, appliedBodyForces);

	//// Computation of GRF forces
	//double loc_heel_x = -0.05;
	//double radius_heel = 0.02;
	//double loc_front_x = 0.15;
	//double radius_front = 0.02;
	//std::vector<adouble> GRF_spheres = Compute_GRF_spheres(Q, U,loc_heel_x,radius_heel,loc_front_x,radius_front);

	Vec3 pos_heelC_in_body(0);
	Vec3 pos_frontC_in_body(0);
	/*pos_heelC_in_body[0] = loc_heel_x;
	pos_heelC_in_body[1] = -radius_heel;
	pos_frontC_in_body[0] = loc_front_x;
	pos_frontC_in_body[1] = -radius_front;*/
	Vector_<Vec3> particleForces;
	particleForces.setToZero();
	Vector_<SpatialVec> GRF_heel;
	GRF_heel.resize(4);
	GRF_heel.setToZero();
	Vector_<SpatialVec> GRF_front;
	GRF_front.resize(4);
	GRF_front.setToZero();

	particleForces.resize(matter->getNumParticles());
	particleForces.setToZero();
	std::cout << "start calc forces" << std::endl;
	//forces->setForceIsDisabled(*Modelstate, hc_heel->getForceIndex(), 0);
	
	hc_heel->calcForceContribution(*Modelstate, appliedBodyForces, particleForces, appliedMobilityForces);
	std::cout << "isdisabled in test " << hc_heel->isDisabled(*Modelstate) << std::endl;
	pos_heelC_in_body = hc_heel->getContactPointInBody(*Modelstate);
	std::cout << "appliedBody forces after calling hc_heel=" << appliedBodyForces << std::endl;
	
	//forces->setForceIsDisabled(*Modelstate, hc_front->getForceIndex(), 0);
	hc_front->calcForceContribution(*Modelstate, appliedBodyForces, particleForces, appliedMobilityForces);
	pos_frontC_in_body = hc_front->getContactPointInBody(*Modelstate);
	std::cout << "appliedBody forces after calling hc_front=" << appliedBodyForces << std::endl;

	//matter->addInStationForce(*Modelstate, *indxpd0, pos_heelC_in_body, GRF_heel[, appliedBodyForces);
	//matter->addInStationForce(*Modelstate, *indxpd0, pos_frontC_in_body, GRF_front, appliedBodyForces);

	/// Residual forces
	matter->calcResidualForceIgnoringConstraints(*Modelstate, appliedMobilityForces, appliedBodyForces, knownUdot, residualMobilityForces);
	std::cout << "residualMobilityForces" << residualMobilityForces << std::endl;
}
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

const Real TOL = 2e-6;

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

#ifdef SimTK_REAL_IS_ADOUBLE
	template <class adouble>
	void assertEqual(adouble val1, adouble val2) {
		ASSERT(fabs(val1.getValue() - val2.getValue()) < TOL);
	}
	template <int N>
	void assertEqual(Vec<N> val1, Vec<N> val2) {
		for (int i = 0; i < N; ++i) {
			ASSERT(fabs(val1[i] - val2[i]) < TOL);
		}
	}
#else
	<class T>
	void assertEqual(T val1, T val2) {
		ASSERT(fabs(val1-val2) < TOL);
	}
	template <int N>
	void assertEqual(Vec<N> val1, Vec<N> val2) {
		for (int i = 0; i < N; ++i) {
			ASSERT(fabs(val1[i] - val2[i]) < TOL);
			
		}
	}
#endif


Real step5(Real x) {
	Real out = 10 * pow(x, 3) - 15 * pow(x, 4) + 6 * pow(x, 5);
	return out;
}

Real stribeck(Real us, Real ud,
	Real uv, Real vrel) {
	Real mu_wet = uv*vrel;
	Real mu_dry;
	Real mu_dry3 = ud;
	Real mu_dry2 = us - (us - ud)*step5((vrel - 1) / 2);
	Real mu_dry1 = us*step5(vrel);
	Real mu_dry_aux = mu_dry2 + ((1 + NTraits<Real>::tanh(1 * (vrel - 3))) / 2)*(mu_dry3 - mu_dry2);
	mu_dry = mu_dry1 + ((1 + NTraits<Real>::tanh(1 * (vrel - 1))) / 2)*(mu_dry_aux - mu_dry1);
	/*if (vrel > 3) {
		mu_dry = ud;
	}
	else if ((vrel >= 1) && (vrel<3)) {
		mu_dry = us - (us - ud)*step5((vrel - 1) / 2);
	}
	else if (vrel<1) {
		mu_dry = us*step5(vrel);
	}*/
	Real mufriction = mu_wet + mu_dry;
	return mufriction;
}

void testForces() {
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    const Vec3 gravity = Vec3(0, -9.8, 0);
    Force::UniformGravity(forces, matter, gravity, 0);
    const Real radius = 0.8;
	Vec3 locSphereInBody(-0.1, 0.1, 0.2);
	Real stiffness = 3e5;
	Real E = stiffness*pow(0.5, 1.5);
	Real k = (4.0 / 3.0)*sqrt(radius)*E;
    const Real c = 0.5;
    const Real us = 0.8;
    const Real ud = 0.8;
    const Real uv = 0.1;
	const Real vt = 0.001;

	double eps = 1e-8;
	double bd = 1500;
	double bv = 3.5;
	double bn = 1;

    Random::Uniform random(0.0, 1.0);
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    MobilizedBody::Translation Mbody(matter.updGround(), Transform(), body, Transform());
    HuntCrossleyForce_smooth hc(forces);
    
	hc.setParameters(stiffness, c, us, ud, uv, vt);
	hc.setLocSphere(locSphereInBody);
	Real offset = 0.0;
	Vec3 normal = Vec3(0, 1, 0);
	hc.setGroundPlane(normal, offset);
	hc.setBodySphere(Mbody);
	hc.setRadiusSphere(radius);

    /*assertEqual(vt, hc.getTransitionVelocity());*/
    State state = system.realizeTopology();
    
    // Position the sphere at a variety of positions and check the normal force.
	Vector VecQi(3, 0.0);
	Vector VecUi(3, 0.0);
	Real depthVel = VecUi[1];
    for (Real height = radius+0.2; height >= 0; height -= 0.1) {
		VecQi[1] = height;
        Mbody.setQFromVector(state, VecQi);
		Mbody.setUFromVector(state, VecUi);
		system.realize(state, Stage::Dynamics);
		const Real depth = radius-(height+locSphereInBody[1]) +offset;

		Real fp = pow((pow(depth, 2) + eps), 3.0 / 4.0);
		Real fv = (1.0 + (3.0 / 2.0)*c*depthVel);
		/*Real fn_nosmooth = k*fp*fv;*/

		//Remove negative parts
		Real fp_nonneg = fp*(1.0 / 2.0 + (1.0 / 2.0)*NTraits<Real>::tanh(1e7*depth));
		Real fv_nonneg = fv*(1.0 / 2.0 + (1.0 / 2.0)*NTraits<Real>::tanh(1e7*(depthVel + (2.0 / (3.0*c)))));

		//Decreasing exponential when there is no contact
		Real auxfun_p = (NTraits<Real>::exp((depth - 0.01) / 0.1))*bn*(1.0 / 2.0 - (1.0 / 2.0)*NTraits<Real>::tanh(bd*(depth - 0.01)));

		//Smooth curves
		Real fn = k*fp*fv*(1.0 / 2.0 + (1.0 / 2.0)*NTraits<Real>::tanh(bd*depth))*(1.0 / 2.0 + (1.0 / 2.0)*NTraits<Real>::tanh(bv*(depthVel + (2.0 / (3.0*c)))));
		std::cout << "k=" << k << " \n fp_nonneg=" << fp_nonneg << "\n fv_nonneg=" << fv_nonneg << std::endl;

		//Add auxiliar function
		fn = auxfun_p + ((1 + NTraits<Real>::tanh(2000 * (depth - 0.002))) / 2)*(fn - auxfun_p);
		std::cout << "height" << height << std::endl;
		std::cout << "Indentation=" << depth << std::endl;
		std::cout << "calculated= " << gravity + Vec3(0, fn, 0) << std::endl;
		std::cout << "from rigidbody=" << system.getRigidBodyForces(state, Stage::Dynamics)[Mbody.getMobilizedBodyIndex()][1] << std::endl;
		assertEqual(system.getRigidBodyForces(state, Stage::Dynamics)[Mbody.getMobilizedBodyIndex()][1], gravity + Vec3(0, fn, 0));

	}
	printf("\n");

    // Now do it with a vertical velocity and see if the dissipation force is correct.
	
	VecQi(3,0.0);
	for (Real height = radius+0.2; height > 0; height -= 0.1) {
		//VecQi[1] = height - locSphereInBody[1];
		VecQi[1] = height ;

		Vector VecUi(3, 0.0);
		depthVel = VecUi[1] ;
		Mbody.setQFromVector(state, VecQi);
        const Real depth = radius-(height+locSphereInBody[1]) +offset;

		Real fh = k *pow((pow(depth, 2) + eps), 3.0 / 4.0);
        for (Real v = -1.0; v <= 1.0; v += 0.1) {
			VecUi[1] = v;
			depthVel = -v;
	
            Mbody.setUFromVector(state, VecUi);
            system.realize(state, Stage::Dynamics);
			Real fp = pow((pow(depth, 2) + eps), 3.0 / 4.0);
			Real fv = (1.0 + (3.0 / 2.0)*c*depthVel);
			/*Real fn_nosmooth = k*fp*fv;*/

			//Remove negative parts
			Real fp_nonneg = fp*(1.0 / 2.0 + (1.0 / 2.0)*NTraits<Real>::tanh(1e7*depth));
			Real fv_nonneg = fv*(1.0 / 2.0 + (1.0 / 2.0)*NTraits<Real>::tanh(1e7*(depthVel + (2.0 / (3.0*c)))));

			//Decreasing exponential when there is no contact
			Real auxfun_p = (NTraits<Real>::exp((depth - 0.01) / 0.1))*bn*(1.0 / 2.0 - (1.0 / 2.0)*NTraits<Real>::tanh(bd*(depth - 0.01)));

			//Smooth curves
			Real fn = k*fp*fv*(1.0 / 2.0 + (1.0 / 2.0)*NTraits<Real>::tanh(bd*depth))*(1.0 / 2.0 + (1.0 / 2.0)*NTraits<Real>::tanh(bv*(depthVel + (2.0 / (3.0*c)))));

			//Add auxiliar function
			fn = auxfun_p + ((1 + NTraits<Real>::tanh(2000 * (depth - 0.002))) / 2)*(fn - auxfun_p);
			assertEqual(system.getRigidBodyForces(state, Stage::Dynamics)[Mbody.getMobilizedBodyIndex()][1], gravity+Vec3(0, fn, 0));
        }
    }
    
    // Do it with a horizontal velocity and see if the friction force is correct.
	VecQi(3, 0.0);
    Vector_<SpatialVec> expectedForce(matter.getNumBodies());
    for (Real height = radius+0.2; height > 0; height -= 0.1) {
		std::cout << "height=" << height << std::endl;
		for (Real v = -1.0; v <= 1.0; v += 0.1) {
			
			VecQi[1] = height ;
			Vector VecUi(3, 0.0);
			VecUi[0] = v;
			Mbody.setQFromVector(state, VecQi);
			
			const Real depth = radius-(height+ locSphereInBody[1])+offset;
			depthVel = 2;
			VecUi[1] = -depthVel;
			
			Mbody.setUFromVector(state, VecUi);
			Vec3 vec3v(0);
			for (int i = 0; i < 3; i++) {
				vec3v[i] = VecUi[i];
			}
			const Real vnormal = dot(vec3v, normal);
			const Vec3 vtangent = vec3v - vnormal*normal;

			Real fp = pow((pow(depth, 2) + eps), 3.0 / 4.0);
			Real fv = (1.0 + (3.0 / 2.0)*c*depthVel);
			//Real fn_nosmooth = k*fp*fv;

			//Remove negative parts
			Real fp_nonneg = fp*(1.0 / 2.0 + (1.0 / 2.0)*NTraits<Real>::tanh(1e7*depth));
			Real fv_nonneg = fv*(1.0 / 2.0 + (1.0 / 2.0)*NTraits<Real>::tanh(1e7*(depthVel + (2.0 / (3.0*c)))));

			Real auxfun_p = (NTraits<Real>::exp((depth - 0.01) / 0.1))*bn*(1.0 / 2.0 - (1.0 / 2.0)*NTraits<Real>::tanh(bd*(depth - 0.01)));


			//Smooth curves
			Real fn = k*fp*fv*(1.0 / 2.0 + (1.0 / 2.0)*NTraits<Real>::tanh(bd*depth))*(1.0 / 2.0 + (1.0 / 2.0)*NTraits<Real>::tanh(bv*(depthVel + (2.0 / (3.0*c)))));

			//Add auxiliar function
			fn = auxfun_p + ((1 + NTraits<Real>::tanh(2000 * (depth - 0.002))) / 2)*(fn - auxfun_p);

		/*for (Real v = -1.0; v <= 1.0; v += 0.1) {
			VecUi[0]=v;
			Mbody.setUFromVector(state, VecUi);*/
			system.realize(state, Stage::Dynamics);
			//Real aux = pow(vtangent[0], 2) + pow(vtangent[1], 2) + pow(vtangent[2], 2) + eps;
			const Real vslip = sqrt((pow(v, 2) + eps));

            const Real vrel = vslip / vt;
            
			Real mu_friction = stribeck(us, ud, uv*vt, vrel);
			std::cout << "vrel=" << vrel << std::endl;
			std::cout << "mu_friction=" << mu_friction << std::endl;

			Real ff = fn*mu_friction;
			std::cout << "ff=" << ff << std::endl;
			std::cout << "v=" << v << std::endl;
			std::cout << "MBody.getU= " << Mbody.getU(state) << std::endl;
			/*Vec3 force = fn*normal;*/
			
			Vec3 vectortangforce = ff*(-vtangent) / vslip;
			Vec3 vectornormforce = fn*normal;

			const Vec3 totalForce = gravity+ vectortangforce+ vectornormforce;
			Vec3 contactForce = vectortangforce + vectornormforce;

            expectedForce = SpatialVec(Vec3(0), Vec3(0));
			Vec3 LocSphereInGround = Mbody.findStationLocationInGround(state, locSphereInBody);
			Vec3 locContactInGround(0);
			locContactInGround = LocSphereInGround - Vec3(0.0, radius, 0.0);
			Vec3 contactPointInSphere = Mbody.findStationAtGroundPoint(state, locContactInGround);
			std::cout << "depth=" << depth << std::endl;
			std::cout << "contact force=" << contactForce << std::endl;
            //Mbody.applyForceToBodyPoint(state, contactPointInSphere, Vec3(ff, fn, 0), expectedForce);
			Mbody.applyForceToBodyPoint(state, contactPointInSphere, contactForce, expectedForce);
			std::cout << "expectedForce0 " << expectedForce << std::endl;
			Mbody.applyForceToBodyPoint(state, Mbody.getBodyMassCenterStation(state), 1.0*gravity, expectedForce);
			std::cout << "expectedForce1 " << expectedForce << std::endl;

            SpatialVec actualForce = system.getRigidBodyForces(state, Stage::Dynamics)[Mbody.getMobilizedBodyIndex()];

			std::cout << "actualForce[0]= " << actualForce[0] << " \n expectedForce[0]=  " << expectedForce[Mbody.getMobilizedBodyIndex()][0] << "\n" << std::endl;
			std::cout << "actualForce[1]= " << actualForce[1] << " \n expectedForce[1]=  " << expectedForce[Mbody.getMobilizedBodyIndex()][1] << "\n\n" << std::endl;
            assertEqual(actualForce[0], expectedForce[Mbody.getMobilizedBodyIndex()][0]);
			std::cout << "actualForce[1]= " << actualForce[1] << " \n expectedForce[1]=  " << expectedForce[Mbody.getMobilizedBodyIndex()][1] << "\n\n" << std::endl;
			assertEqual(actualForce[1], expectedForce[Mbody.getMobilizedBodyIndex()][1]);
        }
    }

	// Test vertical force with an offset of 2 cm
	std::cout << "Test vertical force with an offset of 2 cm" << std::endl;
	offset = 0.02;
	hc.setGroundPlane(Vec3(0, 1, 0), offset);
	VecQi=Vector(3,0.0);
	VecUi=Vector(3,0.0);
	depthVel = VecUi[1];
	for (Real height = radius + 0.2; height > 0; height -= 0.1) {
		VecQi[1] = height;
		Mbody.setQFromVector(state, VecQi);
		Mbody.setUFromVector(state, VecUi);
		system.realize(state, Stage::Dynamics);
		const Real depth = radius - (height+ locSphereInBody[1]) + offset;
		Real fp = pow((pow(depth, 2) + eps), 3.0 / 4.0);
		Real fv = (1.0 + (3.0 / 2.0)*c*depthVel);

		//Remove negative parts
		Real fp_nonneg = fp*(1.0 / 2.0 + (1.0 / 2.0)*NTraits<Real>::tanh(1e7*depth));
		Real fv_nonneg = fv*(1.0 / 2.0 + (1.0 / 2.0)*NTraits<Real>::tanh(1e7*(depthVel + (2.0 / (3.0*c)))));

		//Decreasing exponential when there is no contact
		Real auxfun_p = (NTraits<Real>::exp((depth - 0.01) / 0.1))*bn*(1.0 / 2.0 - (1.0 / 2.0)*NTraits<Real>::tanh(bd*(depth - 0.01)));

		//Smooth curves
		Real fn = k*fp*fv*(1.0 / 2.0 + (1.0 / 2.0)*NTraits<Real>::tanh(bd*depth))*(1.0 / 2.0 + (1.0 / 2.0)*NTraits<Real>::tanh(bv*(depthVel + (2.0 / (3.0*c)))));

		//Add auxiliar function
		fn = auxfun_p + ((1 + NTraits<Real>::tanh(2000 * (depth - 0.002))) / 2)*(fn - auxfun_p);

		std::cout << "actual force=" << gravity + Vec3(0, fn, 0) << std::endl;
		std::cout << "expected force= " << system.getRigidBodyForces(state, Stage::Dynamics)[Mbody.getMobilizedBodyIndex()][1] << std::endl;

		assertEqual(system.getRigidBodyForces(state, Stage::Dynamics)[Mbody.getMobilizedBodyIndex()][1], gravity + Vec3(0, fn, 0));

	}


}

int main() {
    /*try {*/
        testForces();
    /*}
    catch(const std::exception& e) {*/
    /*    cout << "exception: " << e.what() << endl;
        return 1;
    }*/
    cout << "Done" << endl;
    return 0;
}

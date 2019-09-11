#ifndef SimTK_SIMBODY_HUNT_CROSSLEY_FORCE_IMPL_H_
#define SimTK_SIMBODY_HUNT_CROSSLEY_FORCE_IMPL_H_

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

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/HuntCrossleyForce_smooth.h"
#include "ForceImpl.h"

namespace SimTK {

class HuntCrossleyForceImpl_smooth : public ForceImpl {
public:
    class Parameters {
    public:
        Parameters() : stiffness(0), dissipation(0), staticFriction(0), 
            dynamicFriction(0), viscousFriction(0), transitionVelocity(0) {
        }
        Parameters(Real stiffness, Real dissipation,
			Real staticFriction, Real dynamicFriction, Real viscousFriction,
			Real trans_velocity)
            : stiffness(stiffness), dissipation(dissipation), 
			staticFriction(staticFriction), dynamicFriction(dynamicFriction),
			viscousFriction(viscousFriction), 
			transitionVelocity(trans_velocity) {}
		Real stiffness, dissipation, staticFriction,
            dynamicFriction, viscousFriction, transitionVelocity;
    };
	struct InstanceVars {
		InstanceVars(const MobilizedBody defbody, const Real defradius, const Vec3 deflocSphere,
			const Parameters defparameters)
			: BodySphere(defbody), RadiusSphere(defradius), LocSphere(deflocSphere), 
			parameters(defparameters) {}
		MobilizedBody BodySphere;
		Real RadiusSphere;
		Vec3 LocSphere;
		Parameters parameters;
	};

	Real            RadiusSphere;
	Vec3            LocSphere;
	MobilizedBody   BodySphere;
	Parameters		parameters;
	Plane			GroundPlane;

	HuntCrossleyForceImpl_smooth(GeneralForceSubsystem& subsystem);
	
	HuntCrossleyForceImpl_smooth* clone() const {
		return new HuntCrossleyForceImpl_smooth(*this);
	}
	
    /** Set body parameters as individual entries. The order is stiffness, 
    * dissipation, staticFriction, dynamicFriction, viscous friction and 
    * trans_velocity */
    void setParameters(Real stiffness, Real dissipation,
        Real staticFriction, Real dynamicFriction, Real viscousFriction, 
        Real transitionVelocity);

    /** get parameters */
    const Parameters& getParameters() const;
    /** obtain a writtable version of Parameters */
    Parameters& updParameters();
    /** Set stiffness */
    void setStiffness(Real stiffness);
    /** Set dissipation */
    void setDissipation(Real dissipation); 
    /** Set static friction */
    void setStaticFriction(Real staticFriction);
    /** Set dynamic friction */
    void setDynamicFriction(Real dynamicFriction);
    /** Set viscous friction */
    void setViscousFriction(Real viscousFriction);
    /** Set transition velocity */
    void setTransitionVelocity(Real transitionVelocity);
	/** Set Ground Plane */
	void setGroundPlane(Vec3 normal, Real offset);

    /** Set Mobilized Body where the sphere is attached */
	void setBodySphere(MobilizedBody bodyInput);
	/** Set location of the sphere using the frame of the body where it is
	attached*/
	void setLocSphere(Vec3 locSphere);
	/** Set radius of the spehere */
	void setRadiusSphere(Real radius);
	/** Get Mobilized Body where the sphere is attached */
    MobilizedBody getBodySphere();
	/** Get location of the sphere using the frame of the body where it is
	attached*/
	Vec3 getLocSphere();
	/** Get radius of the spehere */
	Real getRadiusSphere();

	/** Calculate contact force given the state, outputs the result as bodyForces 
   (forces in Ground frame) or as mobilityforces (generalized forces) */
	void calcForce(const State& state, Vector_<SpatialVec>& bodyForces,
		Vector_<Vec3>& particleForces, Vector& mobilityForces) const override;
	Real calcPotentialEnergy(const State& state) const;
	/** get the location of the contact in body's frame*/
	Vec3 getContactPointInBody(const State& state);

    void realizeTopology(State& state) const;

	/*void realizeTopology(State& s) const {
		HuntCrossleyForceImpl_smooth* mThis = const_cast<HuntCrossleyForceImpl_smooth*>(this);

		const InstanceVars iv(defX_B1F, defX_B2M, defK, defC);
		mThis->instanceVarsIx = getForceSubsystem()
			.allocateDiscreteVariable(s, Stage::Instance,
				new Value<InstanceVars>(iv));

		Vector einit(1, Real(0));
		mThis->dissipatedEnergyIx = getForceSubsystem().allocateZ(s, einit);

		mThis->positionCacheIx = getForceSubsystem().allocateCacheEntry(s,
			Stage::Position, Stage::Infinity, new Value<PositionCache>());
		mThis->potEnergyCacheIx = getForceSubsystem().allocateCacheEntry(s,
			Stage::Position, Stage::Infinity, new Value<Real>(NaN));
		mThis->velocityCacheIx = getForceSubsystem().allocateCacheEntry(s,
			Stage::Velocity, Stage::Infinity, new Value<VelocityCache>());
		mThis->forceCacheIx = getForceSubsystem().allocateCacheEntry(s,
			Stage::Velocity, Stage::Infinity, new Value<ForceCache>());
	}*/

private:
	const GeneralForceSubsystem&          subsystem;
	mutable CacheEntryIndex               energyCacheIndex;

	//const GeneralForceSubsystem&          subsystem;
    void getContactPointSphere(const State& state, 
		Vec3& contactPointPos) const;
	Real step5(Real x) const;
	Real stribeck(Real us, Real ud, Real uv, Real vrel) const;
	
	
	
	//bool isPotentialEnergyValid(const State& s) const
	//{
	//	return getForceSubsystem().isCacheValueRealized(s, potEnergyCacheIx);
	//}
	//void ensurePotentialEnergyValid(const State&) const;
	//const InstanceVars& getInstanceVars(const State& s) const
	//{
	//	return Value<InstanceVars>::downcast
	//	(getForceSubsystem().getDiscreteVariable(s, instanceVarsIx));
	//}
	//InstanceVars& updInstanceVars(State& s) const
	//{
	//	return Value<InstanceVars>::updDowncast
	//	(getForceSubsystem().updDiscreteVariable(s, instanceVarsIx));
	//}
	//// TOPOLOGY CACHE
	//DiscreteVariableIndex           instanceVarsIx;
	//ZIndex                          dissipatedEnergyIx;
	//CacheEntryIndex                 positionCacheIx;
	//CacheEntryIndex                 potEnergyCacheIx;
	//CacheEntryIndex                 velocityCacheIx;
	//CacheEntryIndex                 forceCacheIx;

    /** Set plane giving orientation and location respect to the body where it
    is attached (in its local frame) */
    //void setPlane(Body body, Vec3 location_plane, Rotation rotation_plane); 

    /** get the Mobilized body where the sphere is attached and its position
    respect to the frame of this body*/
    /** get the Mobilized body where the plane is attached and its position and
    orientation respect to the frame of this body */
    //void getPlane(MobilizedBody body, Vec3& location_plane, 
    //    Rotation& rotation_plane);
    
    //Vec3									LocSphereBody2;
    //Real									RadiusSphere;
    //Transform								TransPlane;
    //MobilizedBody						    Body1; // body containing the plane (normally ground);
    //MobilizedBody							Body2; // body containing the sphere

};



} // namespace SimTK

#endif // SimTK_SIMBODY_HUNT_CROSSLEY_FORCE_IMPL_H_

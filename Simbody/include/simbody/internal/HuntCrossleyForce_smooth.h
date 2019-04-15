#ifndef SimTK_SIMBODY_HUNT_CROSSLEY_FORCE_SMOOTH_H_
#define SimTK_SIMBODY_HUNT_CROSSLEY_FORCE_SMOOTH_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
 * Authors: Antoine Falisse, Gil Serrancoli                                                     *
 * Contributors:                                                             *
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
#include "simbody/internal/Force.h"

namespace SimTK {

class HuntCrossleyForceImpl_smooth;

/**
 * This class models the forces generated by simple point contacts between a
 * sphere and a half space. The contact model is a smooth (i.e., twice
 * continuously differentiable) approximation of the Hunt-Crossley contact
 * model available in Simbody. The proposed implementation differs from the
 * original implementation as it does not rely on a GeneralContactSubsystem.
 * Instead, it assumes contact between a sphere and a fixed half space that
 * can have an offset and an inclinationi with respect to the ground.
 *
 * This smooth contact model was primarily developed for use in gradient-based
 * optimization (e.g., in predictive simulations of walking to simulate the
 * contacts between the feet and the ground).
 *
 * The smoothing of the if statements in the original contact model is
 * performed by using hyperbolic tangent functions. For example, the following
 * if statement:
 *      y = 0, if x < d
 *      y = a, if x >= d
 * can be approximated by:
 *      f = 0.5 + 0.5*tanh(b(x-d))
 *      y = a*f
 * where b is a parameter determining the smoothness of the transition.
 *
 * Similar to the original implementation, this contact model includes
 * components for the normal restoring force, dissipation in the material,
 * and surface friction. This force is only applied to point contacts.
 *
 * <h1>Normal Force Components</h1>
 *
 * The force in the normal direction is based on a model due to Hunt & Crossley
 * : K. H. Hunt and F. R. E. Crossley, "Coefficient of Restitution Interpreted
 * as Damping in Vibroimpact,"ASME Journal of Applied Mechanics, pp. 440-445,
 * June 1975. This is a continuous model based on Hertz elastic contact theory,
 * which correctly reproduces the empirically observed dependence on velocity
 * of coefficient of restitution, where e=(1-cv) for (small) impact velocity v
 * and a material property c with units 1/v. Note that c can be measured right
 * off the coefficient of restitution-vs.-velocity curves: it is the absolute
 * value of the slope at low velocities.
 *
 * Given a collision between a sphere and a plane, we can generate a contact
 * force from this equation
 *      f = kx^n(1 + 3/2 cv)
 * where k is a stiffness constant incorporating material properties
 * and geometry (to be defined below), x is penetration depth and
 * v = dx/dt is penetration rate (positive during penetration and
 * negative during rebound). Exponent n depends on the surface
 * geometry. For Hertz contact where the geometry can be approximated
 * by sphere (or sphere-plane) interactions, which is all we are
 * currently handling here, n=3/2.
 *
 * Stiffness k is defined in terms of the relative radius of curvature R and
 * effective plane-strain modulus E, each of which is a combination of
 * the description of the two individual contacting elements:
 *
 * <pre>
 *          R1*R2                                         E2^(2/3)
 *     R = -------,  E = (s1 * E1^(2/3))^(3/2),  s1= -------------------
 *         R1 + R2                                   E1^(2/3) + E2^(2/3)
 * </pre>
 *
 *     c = c1*s1 + c2*(1-s1)
 *     k = (4/3) sqrt(R) E
 *     f = k x^(3/2) (1 + 3/2 c xdot)
 *     pe = 2/5 k x^(5/2)
 * Also, we can calculate the contact patch radius a as
 *     a = sqrt(R*x)
 *
 * In the above, E1 and E2 are the *plane strain* moduli. If you have instead
 * Young's modulus Y1 and Poisson's ratio p1, then E1=Y1/(1-p1^2). The interface
 * to this subsystem asks for E1 (pressure/%strain) and c1 (1/velocity), and
 * E2,c2 only.
 *
 * <h1>Friction Force</h1>
 *
 * The friction force is based on a model by Michael Hollars:
 *
 * f = fn*[min(vs/vt,1)*(ud+2(us-ud)/(1+(vs/vt)^2))+uv*vs]
 *
 * where fn is the normal force at the contact point, vs is the slip (tangential)
 * velocity of the two bodies at the contact point, vt is a transition velocity
 * (see below), and us, ud, and uv are the coefficients of static, dynamic, and
 * viscous friction respectively.  Each of the three friction coefficients is calculated
 * based on the friction coefficients of the two bodies in contact:
 *
 * u = 2*u1*u2/(u1+u2)
 *
 * Because the friction force is a continuous function of the slip velocity, this
 * model cannot represent stiction; as long as a tangential force is applied, the
 * two bodies will move relative to each other.  There will always be a nonzero
 * drift, no matter how small the force is.  The transition velocity vt acts as an
 * upper limit on the drift velocity.  By setting vt to a sufficiently small value,
 * the drift velocity can be made arbitrarily small, at the cost of making the
 * equations of motion very stiff.  The default value of vt is 0.01.
 */
class SimTK_SIMBODY_EXPORT HuntCrossleyForce_smooth : public Force {
public:
    /**
     * Create a Hunt-Crossley contact model.
     *
     * @param forces the subsystem that will own this HuntCrossleyForce element
     */
    HuntCrossleyForce_smooth(GeneralForceSubsystem& forces);
    /**
     * Set the material parameters for a surface in the contact set.
     *
     * @param stiffness       the stiffness constant (k) for the body
     * @param dissipation     the dissipation coefficient (c) for the body
     * @param staticFriction  the coefficient of static friction (us) for the body
     * @param dynamicFriction the coefficient of dynamic friction (ud) for the body
     * @param viscousFriction the coefficient of viscous friction (uv) for the body
     */
    void setParameters(Real stiffness, Real dissipation, Real staticFriction,
       Real dynamicFriction, Real viscousFriction, Real transitionVelocity);
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
	/** Set ground plane */
	void setGroundPlane(Vec3 normal, Real offset);
	/** Get location of the contact point in the body frame */
	Vec3 getContactPointInBody(const State& state);
	/** Set Mobilized Body the sphere is attached to */
	void setBodySphere(MobilizedBody bodyInput);
	/** Set location of the sphere in the body frame */
	void setLocSphere(Vec3 locSphere);
	/** Set radius of the sphere */
	void setRadiusSphere(Real radius);
	/** Get Mobilized Body the sphere is attached to */
	MobilizedBody getBodySphere();
	/** Get location of the sphere in the body frame */
	Vec3 getLocSphere();
	/** Set radius of the sphere */
	Real setRadiusSphere();

    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(HuntCrossleyForce_smooth,
        HuntCrossleyForceImpl_smooth, Force);
};

} // namespace SimTK

#endif // SimTK_SIMBODY_HUNT_CROSSLEY_FORCE_SMOOTH_H_

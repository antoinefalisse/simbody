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
 * Authors: Peter Eastman                                                     *
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
 * can have an offset and a slope with respect to the ground.
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
 * To use it, do the following:
 *
 * <ol>
 * <li>Add a GeneralForceSubsystem to a MultibodySystem.</li>
 * <li>Add a HuntCrossleyForce_smooth to the GeneralForceSubsystem, call
 * setParameters(), setContactSphere(), setRadiusContactSphere(),
 * setLocContactSphere(), and setContactPlane(). </li>
 * </ol>
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
 * The original Hertz force between a sphere and a plane is given by:
 *      fh = 4/3*k*x*(R*k*x)^(1/2)
 * where k = 0.5*stiffness^(2/3) where stiffness is the effective Young's
 * modulus, which is assumed identic for both contacting materials (i.e.,
 * sphere and plane), x is penetration depth, and R is sphere radius.
 * In the smooth approximation, we use the expression:
 *      fh_pos = 4/3*k*(R*k*)^(1/2)*((x^2+eps)^(1/2))^(3/2)
 *      fh_smooth = fh_pos*(1./2.+(1./2.)*tanh(bd*x));
 * where eps=1e-5 enforces a small force even when there is no contact between
 * the sphere and the plane, and bd=300 determines the smoothness of the tanh
 * transition.
 *
 * The original Hunt-Crossley force is given by:
 *      f = fh*(1+3/2*c*v)
 * where c is dissipation and v is penetration rate.
 * In the smooth approximation, we use the expression:
 *      f_pos = fh_smooth*(1.+(3./2.)*c*v);
 *      f_smooth = f_pos*(1./2.+(1./2.)*tanh(bv*(v+(2./(3.*c)))));
 * where bv=50 determines the smoothness of the tanh transition.
 *
 * <h1>Friction Force</h1>
 *
 * The friction force is based on a model by Michael Hollars:
 *
 * f = f_smooth*[min(vs/vt,1)*(ud+2(us-ud)/(1+(vs/vt)^2))+uv*vs]
 *
 * where f_smooth is the smooth normal force at the contact point, vs is the
 * slip (tangential) velocity of the two bodies at the contact point, vt is a
 * transition velocity (see below), and us, ud, and uv are the coefficients of
 * static, dynamic, and viscous friction respectively. Each of the three
 * friction coefficients is calculated based on the friction coefficients of
 * the two bodies in contact:
 *      u = 2*u1*u2/(u1+u2)
 * In the smooth approximation, we assume the same coefficients for both
 * contacting materials.
 *
 * Because the friction force is a continuous function of the slip velocity,
 * this model cannot represent stiction; as long as a tangential force is
 * applied, the two bodies will move relative to each other. There will always
 * be a nonzero drift, no matter how small the force is. The transition
 * velocity vt acts as an upper limit on the drift velocity. By setting vt to a
 * sufficiently small value, the drift velocity can be made arbitrarily small,
 * at the cost of making the equations of motion very stiff. The default value
 * of vt is 0.01.
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
     * Set the contact material parameters
     *
     * @param stiffness       the stiffness constant
     * @param dissipation     the dissipation coefficient (c)
     * @param staticFriction  the coefficient of static friction (us)
     * @param dynamicFriction the coefficient of dynamic friction (ud)
     * @param viscousFriction the coefficient of viscous friction (uv)
     */
    void setParameters(Real stiffness, Real dissipation, Real staticFriction,
       Real dynamicFriction, Real viscousFriction, Real transitionVelocity);
    /** Set the stiffness constant */
    void setStiffness(Real stiffness);
    /** Set the dissipation coefficient */
    void setDissipation(Real dissipation);
    /** Set the coefficient of static friction */
    void setStaticFriction(Real staticFriction);
    /** Set the coefficient of dynamic friction */
    void setDynamicFriction(Real dynamicFriction);
    /** Set the coefficient of viscous friction */
    void setViscousFriction(Real viscousFriction);
    /** Set the transition velocity */
    void setTransitionVelocity(Real transitionVelocity);
    /**
    * Set contact plane
    *
    * @param normal     direction of the normal to the plane of contact
    * @param offset     distance to the ground origin along the normal
    */
    void setContactPlane(Vec3 normal, Real offset);
    /** Set the Mobilized Body to which the contact sphere is attached */
    void setContactSphere(MobilizedBody bodyInput);
    /** Set the location of the contact sphere in the body frame */
    void setLocContactSphere(Vec3 locSphere);
    /** Set the radius of the contact sphere */
    void setRadiusContactSphere(Real radius);
    /** Get the Mobilized Body to which the contact sphere is attached */
    MobilizedBody getBodySphere();
    /** Get the location of the contact sphere in the body frame */
    Vec3 getLocContactSphere();
    /** Set the radius of the sphere */
    Real setRadiusContactSphere();

    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(HuntCrossleyForce_smooth,
        HuntCrossleyForceImpl_smooth, Force);
};

} // namespace SimTK

#endif // SimTK_SIMBODY_HUNT_CROSSLEY_FORCE_SMOOTH_H_

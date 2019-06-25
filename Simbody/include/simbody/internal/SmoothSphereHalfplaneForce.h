#ifndef SimTK_SIMBODY_SMOOTH_SPHERE_HALFPLANE_FORCE_H_
#define SimTK_SIMBODY_SMOOTH_SPHERE_HALFPLANE_FORCE_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-19 Stanford University and the Authors.        *
 * Authors: Antoine Falisse, Gil Serrancoli                                   *
 * Contributors: Peter Eastman                                                *
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

class SmoothSphereHalfplaneForceImpl;

/**
 * This class models the forces generated by simple point contacts between a
 * sphere and a half space. The contact model is a smooth (i.e., twice
 * continuously differentiable) approximation of the HuntCrossleyForce already
 * available in Simbody. The proposed implementation was designed for use
 * with gradient-based optimization algorithms and algorithmic/automatic
 * differentiation. To this aim, conditional if statements were approximated by
 * using hyperbolic tangent functions. For example, the following if statement:
 * <pre>     y = 0, if x < d </pre>
 * <pre>     y = a, if x >= d </pre>
 * can be approximated by:
 * <pre>     f = 0.5 + 0.5 tanh(b(x-d)) </pre>
 * <pre>     y = a f </pre>
 * where b is a parameter determining the smoothness of the transition.
 *
 * The proposed implementation does not rely on a GeneralContactSubsystem.
 * Instead, it assumes contact between a sphere and a fixed half space that
 * can have an offset and a slope with respect to the ground. The contact model
 * includes components for the normal restoring force, dissipation in the
 * material, and surface friction. The force is only applied to point contacts.
 *
 * To use it, do the following:
 *
 * <ol>
 * <li>Add a GeneralForceSubsystem to a MultibodySystem.</li>
 * <li>Add a SmoothSphereHalfplaneForce to the GeneralForceSubsystem, call
 * setParameters(), setContactSphereInBody(), setContactSphereRadius(),
 * setContactSphereLocationInBody(), and setContactPlane(). </li>
 * </ol>
 *
 */
class SimTK_SIMBODY_EXPORT SmoothSphereHalfplaneForce : public Force {
public:
    /**
     * @param forces the subsystem that will own this
        SmoothSphereHalfplaneForce element
     */
    SmoothSphereHalfplaneForce(GeneralForceSubsystem& forces);
    /**
     * Set the contact material parameters.
     *
     * @param stiffness the stiffness constant, default is 1
     * @param dissipation the dissipation coefficient, default is 0
     * @param staticFriction the coefficient of static friction, default is 0
     * @param dynamicFriction the coefficient of dynamic friction, default is 0
     * @param viscousFriction the coefficient of viscous friction, default is 0
     * @param transitionVelocity the transition velocity, default is 0.01
     * @param cf the constant that enforces a small contact force even when
           there is no contact between the sphere and the plane to ensure
           differentiability of the model, default is 1e-5
     * @param bd the parameter that determines the smoothness of the transition
           of the tanh used to smooth the Hertz force, default is 300
     * @param bv the parameter that determines the smoothness of the transition
           of the tanh used to smooth the Hunt-Crossley force, default is 50
     */
    void setParameters(Real stiffness, Real dissipation, Real staticFriction,
       Real dynamicFriction, Real viscousFriction, Real transitionVelocity,
       Real eps, Real bd, Real bv);
    /** Set the stiffness constant. */
    void setStiffness(Real stiffness);
    /** Set the dissipation coefficient. */
    void setDissipation(Real dissipation);
    /** Set the coefficient of static friction. */
    void setStaticFriction(Real staticFriction);
    /** Set the coefficient of dynamic friction. */
    void setDynamicFriction(Real dynamicFriction);
    /** Set the coefficient of viscous friction. */
    void setViscousFriction(Real viscousFriction);
    /** Set the transition velocity. */
    void setTransitionVelocity(Real transitionVelocity);
    /** Set the constant that enforces a small contact force even when there
        is no contact between the sphere and the plane. */
    void setConstantContactForce(Real cf);
    /** Set the parameter that determines the smoothness of the transition
        of the tanh used to smooth the Hertz force. */
    void setParameterTanhHertzForce(Real bd);
    /** Set the parameter that determines the smoothness of the transition
        of the tanh used to smooth the Hunt-Crossley force. */
    void setParameterTanhHuntCrossleyForce(Real bv);
    /**
    * Set the contact plane.
    *
    * @param normal     direction of the normal to the plane of contact
    * @param offset     distance to the ground origin along the normal
    */
    void setContactPlane(Vec3 normal, Real offset);
    /** Set the MobilizedBody to which the contact sphere is attached. */
    void setContactSphereInBody(MobilizedBody bodyInput);
    /** Set the location of the contact sphere in the body frame. */
    void setContactSphereLocationInBody(Vec3 locationSphere);
    /** Set the radius of the contact sphere. */
    void setContactSphereRadius(Real radius);
    /** Get the MobilizedBody to which the contact sphere is attached. */
    MobilizedBody getBodySphere();
    /** Get the location of the contact sphere in the body frame. */
    Vec3 getContactSphereLocationInBody();
    /** Get the radius of the contact sphere. */
    Real getContactSphereRadius();

    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(SmoothSphereHalfplaneForce,
        SmoothSphereHalfplaneForceImpl, Force);
};

} // namespace SimTK

#endif // SimTK_SIMBODY_SMOOTH_SPHERE_HALFPLANE_FORCE_H_

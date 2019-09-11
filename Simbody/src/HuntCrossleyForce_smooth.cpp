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

#include <iostream>

#include "SimTKmath.h"

#include "simbody/internal/common.h"
#include "simbody/internal/GeneralContactSubsystem.h"
#include "simbody/internal/MobilizedBody.h"

#include "HuntCrossleyForceImpl_smooth.h"
#include <fstream>



namespace SimTK {

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(HuntCrossleyForce_smooth,
    HuntCrossleyForceImpl_smooth, Force);

HuntCrossleyForce_smooth::HuntCrossleyForce_smooth(GeneralForceSubsystem& forces) :
    Force(new HuntCrossleyForceImpl_smooth(forces)) {
    updImpl().setForceSubsystem(forces, forces.adoptForce(*this));
}

HuntCrossleyForceImpl_smooth::HuntCrossleyForceImpl_smooth(GeneralForceSubsystem& subsystem)  :
    subsystem(subsystem){
}

void HuntCrossleyForceImpl_smooth::realizeTopology(State& state) const {
    energyCacheIndex = state.allocateCacheEntry(subsystem.getMySubsystemIndex(), Stage::Dynamics, new Value<Real>());
    Stage st = state.getSubsystemStage(subsystem.getMySubsystemIndex());
}

void HuntCrossleyForce_smooth::setParameters
(Real stiffness, Real dissipation, Real staticFriction, Real dynamicFriction,
    Real viscousFriction, Real transitionVelocity) {
    updImpl().setParameters(stiffness, dissipation, staticFriction,
        dynamicFriction, viscousFriction, transitionVelocity);
}

void HuntCrossleyForceImpl_smooth::setParameters
(Real stiffness, Real dissipation, 	Real staticFriction, Real dynamicFriction,
    Real viscousFriction, Real transitionVelocity) {
    updParameters() = Parameters(stiffness, dissipation, staticFriction,
        dynamicFriction, viscousFriction, transitionVelocity);
}

void HuntCrossleyForce_smooth::setStiffness(Real stiffness) {
    updImpl().parameters.stiffness= stiffness;
}

void HuntCrossleyForceImpl_smooth::setStiffness(Real stiffness) {
    parameters.stiffness = stiffness;
}

void HuntCrossleyForce_smooth::setDissipation(Real dissipation) {
    updImpl().parameters.dissipation = dissipation;
}

void HuntCrossleyForceImpl_smooth::setDissipation(Real dissipation) {
    parameters.dissipation = dissipation;
}

void HuntCrossleyForce_smooth::setStaticFriction(Real staticFriction) {
    updImpl().parameters.staticFriction = staticFriction;
}

void HuntCrossleyForceImpl_smooth::setStaticFriction(Real staticFriction) {
    parameters.staticFriction = staticFriction;
}

void HuntCrossleyForce_smooth::setDynamicFriction(Real dynamicFriction) {
    updImpl().parameters.dynamicFriction = dynamicFriction;
}

void HuntCrossleyForceImpl_smooth::setDynamicFriction(Real dynamicFriction) {
    parameters.dynamicFriction = dynamicFriction;
}

void HuntCrossleyForce_smooth::setViscousFriction(Real viscousFriction) {
    updImpl().parameters.viscousFriction = viscousFriction;
}

void HuntCrossleyForceImpl_smooth::setViscousFriction(Real viscousFriction) {
    parameters.viscousFriction = viscousFriction;
}

void HuntCrossleyForce_smooth::setTransitionVelocity(Real transitionVelocity) {
    updImpl().parameters.transitionVelocity = transitionVelocity;
}

void HuntCrossleyForceImpl_smooth::setTransitionVelocity(Real transitionVelocity) {
    parameters.transitionVelocity = transitionVelocity;
}

void HuntCrossleyForce_smooth::setGroundPlane(Vec3 normal, Real offset) {
    updImpl().setGroundPlane(normal, offset);
}

void HuntCrossleyForceImpl_smooth::setGroundPlane(Vec3 normal, Real offset) {
    GroundPlane=Plane(normal, offset);
}

void HuntCrossleyForce_smooth::setBodySphere(MobilizedBody bodyInput) {
    updImpl().BodySphere = bodyInput;
}

void HuntCrossleyForceImpl_smooth::setBodySphere(MobilizedBody bodyInput) {
    BodySphere = bodyInput;
}

void HuntCrossleyForce_smooth::setLocSphere(Vec3 locSphere) {
    updImpl().LocSphere = locSphere;
}

void HuntCrossleyForceImpl_smooth::setLocSphere(Vec3 locSphere) {
    LocSphere = locSphere;
}

void HuntCrossleyForce_smooth::setRadiusSphere(Real radius) {
    updImpl().RadiusSphere = radius;
}

void HuntCrossleyForceImpl_smooth::setRadiusSphere(Real radius){
    RadiusSphere = radius;
}

MobilizedBody HuntCrossleyForce_smooth::getBodySphere() {
    return updImpl().BodySphere;
}

MobilizedBody HuntCrossleyForceImpl_smooth::getBodySphere() {
    return BodySphere;
}

Vec3 HuntCrossleyForce_smooth::getLocSphere() {
    return updImpl().LocSphere;
}

Vec3 HuntCrossleyForceImpl_smooth::getLocSphere() {
    return LocSphere;
}

Real HuntCrossleyForce_smooth::setRadiusSphere() {
    return updImpl().RadiusSphere;
}

Real HuntCrossleyForceImpl_smooth::getRadiusSphere() {
    return RadiusSphere;
}

const HuntCrossleyForceImpl_smooth::Parameters& HuntCrossleyForceImpl_smooth::
getParameters() const {
    return parameters;
}

HuntCrossleyForceImpl_smooth::Parameters& HuntCrossleyForceImpl_smooth::
updParameters() {
   return parameters;
}

void HuntCrossleyForceImpl_smooth::getContactPointSphere(const State& state,
    Vec3& contactPointPos) const {

    Vec3 posSphereInGround=
        BodySphere.findStationLocationInGround(state, LocSphere);

    contactPointPos = posSphereInGround - RadiusSphere*GroundPlane.getNormal();
}

Vec3 HuntCrossleyForce_smooth::getContactPointInBody(const State& state)
    {
    return updImpl().getContactPointInBody(state);
}

Vec3 HuntCrossleyForceImpl_smooth::getContactPointInBody(const State& state)
{
    Vec3 posSphereInGround =
        BodySphere.findStationLocationInGround(state, LocSphere);
    Vec3 contactPointPos = posSphereInGround - Vec3(0, RadiusSphere, 0);

    return BodySphere.findStationAtGroundPoint(state, contactPointPos);
}

void HuntCrossleyForceImpl_smooth::calcForce(const State& state,
    Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces,
    Vector& mobilityForces) const {

    Vec3 contactPointPos;

    getContactPointSphere(state, contactPointPos);
    Real Indentation = - GroundPlane.getDistance(contactPointPos);

    // Adjust the contact location based on the relative stiffness of the two
    // materials. Here we assume, as in the original Simbody Hunt-Crossley
    // contact model, that both materials have the same relative stiffnesses.
    // As described in Sherman(2011), the point of contact will then be
    // located midway between the two surfaces. We therefore need to add half
    // of the penetration to contactPointPos that was determined as the
    // location of the contact sphere center minus its radius.
    Vec3 normal = GroundPlane.getNormal();
    Vec3 contactPointPosAdj = contactPointPos+Real(1./2.)*Indentation*normal;

    Vec3 contactPointPosAdjInB =
        BodySphere.findStationAtGroundPoint(state, contactPointPosAdj);
    Vec3 contactPointVel =
        BodySphere.findStationVelocityInGround(state, contactPointPosAdjInB);

    const Vec3 v = contactPointVel;
    const Real vnormal = dot(v, normal);
    const Vec3 vtangent = v - vnormal*normal;
    Real IndentationVel = -vnormal;

    Parameters parameters = getParameters();
    Real stiffness = parameters.stiffness;
    Real dissipation = parameters.dissipation;
    const Real vt = parameters.transitionVelocity;
    const Real us = parameters.staticFriction;
    const Real ud = parameters.dynamicFriction;
    const Real uv = parameters.viscousFriction;

    double eps = 1e-5;
    double bv = 50;
    double bd = 300;
    // relationship between stiffness and elastic modulus (s1 = 0.5)
    Real k = (1./2.)*NTraits<Real>::pow(stiffness, (2./3.));
    Real fH = (4./3.)*k*NTraits<Real>::sqrt(RadiusSphere*k)*NTraits<Real>::pow(NTraits<Real>::sqrt(Indentation*Indentation+eps),(3./2.));
    Real c = dissipation;
    Real fHd = fH*(1.+(3./2.)*c*IndentationVel);
    Real fn = fHd*(1./2.+(1./2.)*NTraits<Real>::tanh(bd*Indentation))*(1./2.+(1./2.)*NTraits<Real>::tanh(bv*(IndentationVel+(2./(3.*c)))));
    Vec3 force = fn*normal;
    Real aux = pow(vtangent[0],2)+pow(vtangent[1],2)+pow(vtangent[2],2)+eps;
    const Real vslip = pow(aux,1./2.);
    Real vrel = vslip / vt;
    const Real ffriction = fn*(NTraits<Real>::min(vrel,Real(1))*(ud+2*(us-ud)/(1+vrel*vrel))+uv*vslip);
    force += ffriction*(-vtangent) / vslip;
    // Apply the force to the bodies.
    BodySphere.applyForceToBodyPoint(state, contactPointPosAdjInB, force, bodyForces);
}

Real HuntCrossleyForceImpl_smooth::calcPotentialEnergy(const State& state)
    const {
    Real PotentialEnergy = 0;
    return PotentialEnergy;
}

} // namespace SimTK


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

#include "SimTKmath.h"

#include "simbody/internal/common.h"
#include "simbody/internal/MobilizedBody.h"

#include "SmoothSphereHalfplaneForceImpl.h"

namespace SimTK {

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(SmoothSphereHalfplaneForce,
    SmoothSphereHalfplaneForceImpl, Force);

SmoothSphereHalfplaneForce::SmoothSphereHalfplaneForce
    (GeneralForceSubsystem& forces) :
    Force(new SmoothSphereHalfplaneForceImpl(forces)) {
    updImpl().setForceSubsystem(forces, forces.adoptForce(*this));
}

SmoothSphereHalfplaneForceImpl::SmoothSphereHalfplaneForceImpl
    (GeneralForceSubsystem& subsystem) :
    subsystem(subsystem){
}

void SmoothSphereHalfplaneForceImpl::realizeTopology(State& state) const {
    energyCacheIndex=state.allocateCacheEntry(subsystem.getMySubsystemIndex(),
        Stage::Dynamics, new Value<Real>());
}

void SmoothSphereHalfplaneForce::setParameters
    (Real stiffness, Real dissipation, Real staticFriction,
    Real dynamicFriction, Real viscousFriction, Real transitionVelocity, Real
    cf, Real bd, Real bv) {
    updImpl().setParameters(stiffness, dissipation, staticFriction,
        dynamicFriction, viscousFriction, transitionVelocity, cf, bd, bv);
}

void SmoothSphereHalfplaneForceImpl::setParameters
    (Real stiffness, Real dissipation,  Real staticFriction,
    Real dynamicFriction, Real viscousFriction, Real transitionVelocity, Real
    cf, Real bd, Real bv) {
    updParameters() = Parameters(stiffness, dissipation, staticFriction,
        dynamicFriction, viscousFriction, transitionVelocity, cf, bd, bv);
}

void SmoothSphereHalfplaneForce::setStiffness(Real stiffness) {
    updImpl().parameters.stiffness = stiffness;
}

void SmoothSphereHalfplaneForceImpl::setStiffness(Real stiffness) {
    parameters.stiffness = stiffness;
}

void SmoothSphereHalfplaneForce::setDissipation(Real dissipation) {
    updImpl().parameters.dissipation = dissipation;
}

void SmoothSphereHalfplaneForceImpl::setDissipation(Real dissipation) {
    parameters.dissipation = dissipation;
}

void SmoothSphereHalfplaneForce::setStaticFriction(Real staticFriction) {
    updImpl().parameters.staticFriction = staticFriction;
}

void SmoothSphereHalfplaneForceImpl::setStaticFriction(Real staticFriction) {
    parameters.staticFriction = staticFriction;
}

void SmoothSphereHalfplaneForce::setDynamicFriction(Real dynamicFriction) {
    updImpl().parameters.dynamicFriction = dynamicFriction;
}

void SmoothSphereHalfplaneForceImpl::setDynamicFriction(Real dynamicFriction) {
    parameters.dynamicFriction = dynamicFriction;
}

void SmoothSphereHalfplaneForce::setViscousFriction(Real viscousFriction) {
    updImpl().parameters.viscousFriction = viscousFriction;
}

void SmoothSphereHalfplaneForceImpl::setViscousFriction(Real viscousFriction) {
    parameters.viscousFriction = viscousFriction;
}

void SmoothSphereHalfplaneForce::setTransitionVelocity(
    Real transitionVelocity) {
    updImpl().parameters.transitionVelocity = transitionVelocity;
}

void SmoothSphereHalfplaneForceImpl::setTransitionVelocity
    (Real transitionVelocity) {
    parameters.transitionVelocity = transitionVelocity;
}

void SmoothSphereHalfplaneForce::setConstantContactForce(Real cf) {
    updImpl().parameters.cf = cf;
}

void SmoothSphereHalfplaneForceImpl::setConstantContactForce(Real cf) {
    parameters.cf = cf;
}

void SmoothSphereHalfplaneForce::setParameterTanhHertzForce(Real bd) {
    updImpl().parameters.bd = bd;
}

void SmoothSphereHalfplaneForceImpl::setParameterTanhHertzForce(Real bd) {
    parameters.bd = bd;
}

void SmoothSphereHalfplaneForce::setParameterTanhHuntCrossleyForce(Real bv) {
    updImpl().parameters.bv = bv;
}

void SmoothSphereHalfplaneForceImpl::setParameterTanhHuntCrossleyForce(
    Real bv) {
    parameters.bv = bv;
}

void SmoothSphereHalfplaneForce::setContactPlane(Vec3 normal, Real offset) {
    updImpl().setContactPlane(normal, offset);
}

void SmoothSphereHalfplaneForceImpl::setContactPlane(Vec3 normal,Real offset) {
    contactPlane = Plane(normal, offset);
}

void SmoothSphereHalfplaneForce::setContactSphereInBody(
    MobilizedBody bodyInput) {
    updImpl().bodySphere = bodyInput;
}

void SmoothSphereHalfplaneForceImpl::setContactSphereInBody(
    MobilizedBody bodyInput){
    bodySphere = bodyInput;
}

void SmoothSphereHalfplaneForce::setContactSphereLocationInBody(
    Vec3 LocContactSphere) {
    updImpl().locationContactSphere = LocContactSphere;
}

void SmoothSphereHalfplaneForceImpl::setContactSphereLocationInBody
    (Vec3 locationContactSphere) {
    locationContactSphere = locationContactSphere;
}

void SmoothSphereHalfplaneForce::setContactSphereRadius(Real radius) {
    updImpl().radiusContactSphere = radius;
}

void SmoothSphereHalfplaneForceImpl::setContactSphereRadius(Real radius) {
    radiusContactSphere = radius;
}

MobilizedBody SmoothSphereHalfplaneForce::getBodySphere() {
    return updImpl().bodySphere;
}

MobilizedBody SmoothSphereHalfplaneForceImpl::getBodySphere() {
    return bodySphere;
}

Vec3 SmoothSphereHalfplaneForce::getContactSphereLocationInBody() {
    return updImpl().locationContactSphere;
}

Vec3 SmoothSphereHalfplaneForceImpl::getContactSphereLocationInBody() {
    return locationContactSphere;
}

Real SmoothSphereHalfplaneForce::getContactSphereRadius() {
    return updImpl().radiusContactSphere;
}

Real SmoothSphereHalfplaneForceImpl::getContactSphereRadius() {
    return radiusContactSphere;
}

const SmoothSphereHalfplaneForceImpl::Parameters&
    SmoothSphereHalfplaneForceImpl::getParameters() const {
    return parameters;
}

SmoothSphereHalfplaneForceImpl::Parameters& SmoothSphereHalfplaneForceImpl::
    updParameters() {
    return parameters;
}

void SmoothSphereHalfplaneForceImpl::getContactPointSphere(const State& state,
    Vec3& contactPointPos) const {
    Vec3 posSphereInGround =
        bodySphere.findStationLocationInGround(state, locationContactSphere);
    contactPointPos = posSphereInGround -
        radiusContactSphere*contactPlane.getNormal();
}

void SmoothSphereHalfplaneForceImpl::calcForce(const State& state,
    Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces,
    Vector& mobilityForces) const {
    // Calculate the indentation based on the contact point location.
    Vec3 contactPointPos;
    getContactPointSphere(state, contactPointPos);
    const Real indentation = - contactPlane.getDistance(contactPointPos);
    // Initialize the potential energy.
    Real& pe = Value<Real>::updDowncast(state.updCacheEntry(
        subsystem.getMySubsystemIndex(), energyCacheIndex)).upd();
    pe = 0.0;
    // Adjust the contact location based on the relative stiffness of the two
    // materials. Here we assume, as in the original Simbody Hunt-Crossley
    // contact model, that both materials have the same relative stiffness.
    // As described in Sherman(2011), the point of contact will then be
    // located midway between the two surfaces. We therefore need to add half
    // the indentation to the contact location that was determined as the
    // location of the contact sphere center minus its radius.
    const Vec3 normal = contactPlane.getNormal();
    const Vec3 contactPointPosAdj =
        contactPointPos+Real(1./2.)*indentation*normal;
    const Vec3 contactPointPosAdjInB =
        bodySphere.findStationAtGroundPoint(state, contactPointPosAdj);
    // Calculate the contact point velocity.
    const Vec3 contactPointVel =
        bodySphere.findStationVelocityInGround(state, contactPointPosAdjInB);
    // Calculate the tangential and indentation velocities.
    const Vec3 v = contactPointVel;
    const Real vnormal = dot(v, normal);
    const Vec3 vtangent = v - vnormal*normal;
    const Real indentationVel = -vnormal;
    // Get the contact model parameters.
    const Parameters parameters = getParameters();
    const Real stiffness = parameters.stiffness;
    const Real dissipation = parameters.dissipation;
    const Real vt = parameters.transitionVelocity;
    const Real us = parameters.staticFriction;
    const Real ud = parameters.dynamicFriction;
    const Real uv = parameters.viscousFriction;
    const Real cf = parameters.cf;
    const Real bd = parameters.bd;
    const Real bv = parameters.bv;
    // Calculate the Hertz force.
    const Real k = (1./2.)*std::pow(stiffness, (2./3.));
    const Real fH = (4./3.)*k*std::sqrt(radiusContactSphere*k)*
        std::pow(std::sqrt(indentation*indentation+cf),(3./2.));
    pe += Real(2./5.)*fH*indentation;
    // Calculate the Hunt-Crossley force.
    const Real c = dissipation;
    const Real fHd = fH*(1.+(3./2.)*c*indentationVel);
    const Real fn = fHd*(1./2.+(1./2.)*std::tanh(bd*indentation))*
        (1./2.+(1./2.)*std::tanh(bv*(indentationVel+(2./(3.*c)))));
    Vec3 force = fn*normal;
    // Calculate the friction force.
    const Real aux = vtangent.normSqr()+cf;
    const Real vslip = pow(aux,1./2.);
    const Real vrel = vslip / vt;
    const Real ffriction = fn*(std::min(vrel,Real(1))*
        (ud+2*(us-ud)/(1+vrel*vrel))+uv*vslip);
    force += ffriction*(-vtangent) / vslip;
    // Apply the force to the bodies.
    bodySphere.applyForceToBodyPoint(state, contactPointPosAdjInB,
        force, bodyForces);
}

Real SmoothSphereHalfplaneForceImpl::calcPotentialEnergy(const State& state)
    const { return Value<Real>::downcast(state.getCacheEntry(
        subsystem.getMySubsystemIndex(), energyCacheIndex)).get();
}

} // namespace SimTK


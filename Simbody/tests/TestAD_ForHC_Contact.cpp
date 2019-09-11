#include "Simbody.h"
#include "SimTKcommon.h"
#include <adolc.h>
#include <adolc_sparse.h>
#include "SimTKmath.h"

using namespace SimTK;

void main() {

	MultibodySystem* systemmodel=new MultibodySystem;
	SimbodyMatterSubsystem* matter=new SimbodyMatterSubsystem(*systemmodel);
	GeneralForceSubsystem* forces=new GeneralForceSubsystem(*systemmodel);
	Force::UniformGravity* gravity=new Force::UniformGravity(*forces, *matter, Vec3(0, Real(-9.81), 0), 0);
	GeneralContactSubsystem* contacts=new GeneralContactSubsystem(*systemmodel);

	// add contact geometries
	ContactGeometry::Sphere ball1(0.8);
	ContactGeometry::HalfSpace plane;
	const Real radius = 0.8;
	const Real k1 = 1.0;
	const Real k2 = 2.0;
	const Real stiffness1 = pow(k1, 2.0 / 3.0);
	const Real stiffness2 = pow(k2, 2.0 / 3.0);
	const Real dissipation1 = 0.5;
	const Real dissipation2 = 1.0;
	const Real us1 = 1.0;
	const Real us2 = 0.7;
	const Real ud1 = 0.5;
	const Real ud2 = 0.2;
	const Real uv1 = 0.1;
	const Real uv2 = 0.05;
	//Describe mass properties for a generic body
	Body::Rigid* ball1Body= new Body::Rigid(MassProperties(1.0, Vec3(0), Inertia(1)));
	ContactSetIndex CSetI = contacts->createContactSet();
	//ContactSetIndex CSetI = systemmodel->updContactSubsystem().createContactSet();

	MobilizedBody::Bushing* ball1MBody= new MobilizedBody::Bushing(matter->updGround(), Transform(Vec3(0)), *ball1Body, Transform(Vec3(0)));
	//MobilizedBody::Translation ball1MBody(matter.updGround(), Transform(), ball1Body, Transform());
	ContactGeometry::HalfSpace HS;
	std::cout << HS.getNormal() << std::endl;
	Transform Tid, Tid2;
	Tid.setToZero();
	contacts->addBody(CSetI, *ball1MBody, ContactGeometry::Sphere(0.8), Tid);
	contacts->addBody(CSetI, matter->updGround(), HS, Transform(Rotation(-0.5*Pi, ZAxis), Vec3(0))); // y < 0
	HuntCrossleyForce HCF(*forces, *contacts, CSetI);
	//HuntCrossleyForce HCF(*forces, *contacts, CSetI);

	HCF.setBodyParameters(ContactSurfaceIndex(0), k1, dissipation1, us1, ud1, uv1);
	HCF.setBodyParameters(ContactSurfaceIndex(1), k2, dissipation2, us2, ud2, uv2);
	//	systemmodel->setContactSubsystem(*contacts);
	HCF.setTransitionVelocity(0.001);
	std::cout << "Tran vel =" << HCF.getTransitionVelocity() << std::endl << std::endl;
	State* state = new State(systemmodel->realizeTopology());
	Vector Qvec(6);

	//Qvec[0] = 0.002;
	//Qvec[1] = 0.003;
	//Qvec[2] = 0.004;
	//Qvec[3] = 0.001;
	//Qvec[4] = 0.001;
	//Qvec[5] = 0.006;
	
	double* xp = new double[6];
	xp[0] = 0.002;
	xp[1] = 0.003;
	xp[2] = 0.004;
	xp[3] = 0.001;
	xp[4] = 0.001;
	xp[5] = 0.006;
	
	/*adouble* Q = new adouble[6];*/
		
	// Start trace
	trace_on(1, 1);
		
	for (int i = 0; i < 6; i++) {
		Qvec[i] <<= xp[i];
	}
	
	std::cout << "vec Qvec " << Qvec << std::endl;
	std::cout << "Mobilized body index for ball1= " << ball1MBody->getMobilizedBodyIndex() << std::endl;
	printf("There are %d forces \n", forces->getNumForces());
	std::cout << "num contact sets " << contacts->getNumContactSets() << std::endl;
	std::cout << "num bodies = " << matter->getNumBodies() << std::endl;
	std::cout << "Num Mobilities =" << matter->getNumMobilities() << std::endl;
	std::cout << "Nq= " << state->getNQ() << " Nu= " << state->getNU() << " Ny=" << state->getNY() << std::endl;
	matter->updMobilizedBody(ball1MBody->getMobilizedBodyIndex()).setQFromVector(*state,Qvec);
	matter->updMobilizedBody(ball1MBody->getMobilizedBodyIndex()).setUFromVector(*state,Vector(6, 0.0));

	/*ball1MBody.setQ(state, Qvec);*/
	//ball1MBody->setUFromVector(*state, Vector(6, 0.0));
	std::cout << "start realizing" << std::endl;
	systemmodel->realize(*state, Stage::Dynamics);
	std::cout << "model realized" << std::endl;

	std::cout << state->getSystemStage() << std::endl;
	/*ContactSnapshot PredCont=tracker.getPredictedContacts(state);*/
	Vector Q=state->getQ();
	Vector U=state->getU();
	std::cout << "Q Vector =" << Q << std::endl;
	
	printf("There are %d forces \n", forces->getNumForces());
		
	//if (nc != 0) {

	Vector_<SpatialVec> RBF = systemmodel->getRigidBodyForces(*state, Stage::Dynamics);
	for (int i = 0; i < systemmodel->getMatterSubsystem().getNumBodies(); ++i) {
		std::cout << "Rigid body forces for body " << i << " =" << RBF[i] << std::endl;

	}
	printf("CalcForceCont done");
	//}
	
	adouble y[6];
	y[0] = RBF[1][0][0];
	y[1] = RBF[1][0][1];
	y[2] = RBF[1][0][2];
	y[3] = RBF[1][1][0];
	y[4] = RBF[1][1][1];
	y[5] = RBF[1][1][2];
	
	double vals_0[6];
	y[0] >>= vals_0[0];
	y[1] >>= vals_0[1];
	y[2] >>= vals_0[2];
	y[3] >>= vals_0[3];
	y[4] >>= vals_0[4];
	y[5] >>= vals_0[5];
	
	//delete Q;
	trace_off();



	system("pause");
	

}
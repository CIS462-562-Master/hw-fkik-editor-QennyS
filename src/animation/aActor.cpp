#include "aActor.h"
#pragma warning(disable : 4018)



/****************************************************************
*
*    	    Actor functions
*
****************************************************************/

AActor::AActor() 
{
	m_pInternalSkeleton = new ASkeleton();
	m_pSkeleton = m_pInternalSkeleton;

	m_BVHController = new BVHController();
	m_BVHController->setActor(this);

	m_IKController = new IKController();
	m_IKController->setActor(this);

	// code to update additional Actor data goes here
	resetGuide();

}

AActor::AActor(const AActor* actor)
{
	*this = *actor;
}

AActor& AActor::operator = (const AActor& actor)
{
	// Performs a deep copy
	if (&actor == this)
	{
		return *this;
	}
	m_pSkeleton = actor.m_pSkeleton;

	// code to update additional Actor data goes here


	return *this;
}

AActor::~AActor()
{
	 delete m_IKController;
	 delete m_BVHController;
	 delete m_pInternalSkeleton;

}

void AActor::clear()
{
	// looks like it is clearing more times than the number of actors.  as a result, m_pSkeleton is not defined for last case.
	m_pSkeleton->clear();  

	// code to update additional Actor data goes here
}

void AActor::update()
{
	if (!m_pSkeleton->getRootNode() )
		 return; // Nothing loaded
	else m_pSkeleton->update();

	// code to update additional Actor data goes here

}

ASkeleton* AActor::getSkeleton()
{
	return m_pSkeleton;
}

void AActor::setSkeleton(ASkeleton* pExternalSkeleton)
{
	m_pSkeleton = pExternalSkeleton;
}

void AActor::resetSkeleton()
{
	m_pSkeleton = m_pInternalSkeleton;
}

BVHController* AActor::getBVHController()
{
	return m_BVHController;
}

IKController* AActor::getIKController()
{
	return m_IKController;
}

void AActor::updateGuideJoint(vec3 guideTargetPos)
{
	if (!m_pSkeleton->getRootNode()) { return; }

	// TODO: 
	// 1.	Set the global position of the guide joint to the global position of the root joint
	// 2.	Set the y component of the guide position to 0
	// 3.	Set the global rotation of the guide joint towards the guideTarget
	vec3 rootPos = m_Guide.getGlobalTranslation() + m_Guide.getGlobalRotation() * m_pSkeleton->getRootNode()->getGlobalTranslation();
	rootPos[1] = 0;
	guideTargetPos[1] = 0;
	m_Guide.setGlobalTranslation(rootPos);

	vec3 z = vec3(0, 0, 1);
	vec3 forward = (guideTargetPos - m_Guide.getGlobalTranslation()).Normalize();
	vec3 axis = z ^ forward;
	double angle = acos(z * forward);
	mat3 rotation;
	rotation.FromAxisAngle(axis, angle);
	m_Guide.setGlobalRotation(rotation);
}

void AActor::solveFootIK(float leftHeight, float rightHeight, bool rotateLeft, bool rotateRight, vec3 leftNormal, vec3 rightNormal)
{
	if (!m_pSkeleton->getRootNode()) { return; }
	AJoint* leftFoot = m_pSkeleton->getJointByID(m_IKController->mLfootID);
	AJoint* rightFoot = m_pSkeleton->getJointByID(m_IKController->mRfootID);

	// TODO: 
	// The normal and the height given are in the world space

	// 1.	Update the local translation of the root based on the left height and the right height
	vec3 pos = m_pSkeleton->getRootNode()->getGlobalTranslation();
	pos[1] = pos[1] + std::max(leftHeight, rightHeight);
	m_pSkeleton->getRootNode()->setLocalTranslation(pos);
	m_pSkeleton->update();

	// 2.	Update the character with Limb-based IK 
	
	// Rotate Foot
	if (rotateLeft)
	{
		//AIKchain IKChain = m_IKController->createIKchain(leftFoot->getID(), -1, m_pSkeleton);
		
		// Update the local orientation of the left foot based on the left normal
		vec3 leftNormalLocal = (leftFoot->getLocal2Global().m_rotation.Inverse() * leftNormal).Normalize();
		float x = leftNormalLocal[0];
		float y = leftNormalLocal[1];
		float z = leftNormalLocal[2];
		mat3 rotation(vec3(y / sqrt(x * x + y * y), x * z / sqrt(x * x + y * y), x), vec3(-x / sqrt(x * x + y * y), y * z / sqrt(x * x + y * y), y),
			vec3(0, -sqrt(x * x + y * y), z));
		ATarget target;
		vec3 leftFootPos = leftFoot->getGlobalTranslation();
		leftFootPos[1] = leftHeight;
		target.setGlobalTranslation(leftFootPos);
		target.setGlobalRotation(rotation);
		m_IKController->IKSolver_Limb(leftFoot->getID(), target);
		//m_IKController->computeLimbIK(target, IKChain, vec3(0, 1, 0), m_pSkeleton);
	}
	if (rotateRight)
	{	
		//AIKchain IKChain = m_IKController->createIKchain(rightFoot->getID(), -1, m_pSkeleton);
		// Update the local orientation of the right foot based on the right normal
		vec3 rightNormalLocal = (rightFoot->getLocal2Global().m_rotation.Inverse() * rightNormal).Normalize();
		float x = rightNormalLocal[0];
		float y = rightNormalLocal[1];
		float z = rightNormalLocal[2];
		mat3 rotation(vec3(y / sqrt(x * x + y * y), (x * z) / sqrt(x * x + y * y), x), vec3(-x / sqrt(x * x + y * y), (y * z) / sqrt(x * x + y * y), y),
			vec3(0, -sqrt(x * x + y * y), z));
		ATarget target;
		vec3 rightFootPos = rightFoot->getGlobalTranslation();
		rightFootPos[1] = rightHeight;
		target.setGlobalTranslation(rightFootPos);
		target.setGlobalRotation(rotation);
		m_IKController->IKSolver_Limb(rightFoot->getID(), target);
		//m_IKController->computeLimbIK(target, IKChain, vec3(0, 1, 0), m_pSkeleton);
	}
	m_pSkeleton->update();
}

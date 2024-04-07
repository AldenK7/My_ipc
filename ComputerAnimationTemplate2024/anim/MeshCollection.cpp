#include "MeshCollection.h"

MeshCollection::MeshCollection(int max) :
	maxSize(max),
	numMeshes(0)
{
	meshes = new Mesh*[max];
}


void MeshCollection::addMesh(Mesh* newMesh) {
	animTcl::OutputMessage("adding mesh");
	meshes[numMeshes] = newMesh;
	numMeshes++;
}

int MeshCollection::getSize() {
	return numMeshes;
}

void MeshCollection::updatePos(Eigen::MatrixXd qt) {
	for (int i = 0; i < numMeshes; i++) {
		Vector pos;
		setVector(
			pos,
			qt(i, 0),
			qt(i, 1),
			qt(i, 2)
		);

		meshes[i]->setTranslation(pos);
	}
}

void MeshCollection::updateRot(Eigen::MatrixXd Qt) {
	for (int i = 0; i < numMeshes; i++) {
		Vector rot;
		setVector(
			rot,
			Qt(i, 0),
			Qt(i, 1),
			Qt(i, 2)
		);

		meshes[i]->setRotation(rot);
	}
}

double MeshCollection::getBarrier(Eigen::Vector3d trans, Eigen::Vector3d rot) {
	// Mesh[1] is stationary
	return meshes[0]->barrierWithMesh(meshes[1], trans, rot);
}

Eigen::MatrixXd MeshCollection::getTranslations() {
	Eigen::MatrixXd translations(numMeshes, 3);

	for (int i = 0; i < numMeshes; i++) {
		translations(i, 0) = meshes[i]->m_pos[0];
		translations(i, 1) = meshes[i]->m_pos[1];
		translations(i, 2) = meshes[i]->m_pos[2];
	}

	return translations;
}

Eigen::MatrixXd MeshCollection::getRotations() {
	Eigen::MatrixXd rotations(numMeshes, 3);

	for (int i = 0; i < numMeshes; i++) {
		rotations(i, 0) = meshes[i]->m_rot[0];
		rotations(i, 1) = meshes[i]->m_rot[1];
		rotations(i, 2) = meshes[i]->m_rot[2];
	}

	return rotations;
}

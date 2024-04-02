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
		pos[0] = qt(i, 0);
		pos[1] = qt(i, 1);
		pos[2] = qt(i, 2);

		meshes[i]->setTranslation(pos);
	}
}

double MeshCollection::getBarrier() {
	return meshes[0]->barrierWithMesh(meshes[1]);
}

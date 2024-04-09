#include "MeshCollection.h"

MeshCollection::MeshCollection(int max) :
	maxSize(max),
	numMeshes(0)
{
	meshes = new Mesh*[max];
	barrier = -1;
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
	Vector pos;
	setVector(
		pos,
		pos[0] = qt(0, 0),
		pos[1] = qt(0, 1),
		pos[2] = qt(0, 2)
	);

	//animTcl::OutputMessage("%d [%f %f %f]", 16, meshes[0]->init_vertices[16](0), meshes[0]->init_vertices[16](1), meshes[0]->init_vertices[16](2));
	//animTcl::OutputMessage("%d [%f %f %f]", 16, meshes[1]->init_vertices[16](0), meshes[1]->init_vertices[16](1), meshes[1]->init_vertices[16](2));

	//animTcl::OutputMessage("%d [%f %f %f]", 16, meshes[0]->vertices[16](0), meshes[0]->vertices[16](1), meshes[0]->vertices[16](2));
	//animTcl::OutputMessage("%d [%f %f %f]", 16, meshes[1]->vertices[16](0), meshes[1]->vertices[16](1), meshes[1]->vertices[16](2));

	//animTcl::OutputMessage("%d [%f %f %f]", 16, meshes[0]->triangles[0].p1(0), meshes[0]->triangles[0].p1(1), meshes[0]->triangles[0].p1(2));
	//animTcl::OutputMessage("%d [%f %f %f]", 16, meshes[1]->triangles[0].p1(0), meshes[1]->triangles[0].p1(1), meshes[1]->triangles[0].p1(2));

	for (int v = 0; v < meshes[0]->num_vertices; v++) {
		//animTcl::OutputMessage("%d [%f %f %f]", v, meshes[0]->vertices[v](0), meshes[0]->vertices[v](1), meshes[0]->vertices[v](2));
	}

	meshes[0]->setTranslation(pos);
}

double boundary(double d, double dhat) {

	if (d >= dhat) {
		return 0;
	}
	else {
		return -pow(d - dhat, 2) * log(d / dhat);
	}
}

double MeshCollection::getBarrier(Eigen::Vector3d trans) {
	return meshes[0]->barrierWithMesh(meshes[1], trans);
}


double MeshCollection::getDistance(Eigen::Vector3d trans) {
	return meshes[0]->minDistance(meshes[1], trans);
}


double MeshCollection::CCD() {
	return 0;
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
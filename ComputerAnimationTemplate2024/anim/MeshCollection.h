#ifndef MESH_COLLECTION_H
#define MESH_COLLECTION_H

#include <shared/defs.h>
#include <util/util.h>
#include <Eigen/Core>
#include <GLmodel/GLmodel.h>
#include "Mesh.h"
#include "AABB.h"
#include "Eigen/Dense"

class MeshCollection
{

public:
	MeshCollection(int max);
		
	int getSize();
	void updatePos(Eigen::MatrixXd qt);
	void updateRot(Eigen::MatrixXd Qt);

	void addMesh(Mesh* newMesh);

	double getBarrier(Eigen::Vector3d trans, Eigen::Vector3d rot);

	Eigen::MatrixXd getTranslations();
	Eigen::MatrixXd getRotations();

protected:
	int maxSize;
	int numMeshes;

	Mesh** meshes;
};
#endif

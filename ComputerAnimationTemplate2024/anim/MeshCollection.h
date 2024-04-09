#ifndef MESH_COLLECTION_H
#define MESH_COLLECTION_H

#include <shared/defs.h>
#include <util/util.h>
#include <Eigen/Core>
#include <GLmodel/GLmodel.h>
#include "Mesh.h"
#include "AABB.h"

class MeshCollection
{

public:
	MeshCollection(int max);
		
	int getSize();
	void updatePos(Eigen::MatrixXd qt);

	void addMesh(Mesh* newMesh);

	double getBarrier(Eigen::Vector3d trans);

	double getDistance(Eigen::Vector3d trans);

	double CCD();

	Eigen::MatrixXd getTranslations();

protected:
	int maxSize;
	int numMeshes;
	
	Mesh** meshes;

	double barrier;

	double dhat = 0.5;
};
#endif

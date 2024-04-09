#pragma once
#ifndef MY_MESH_H
#define MY_MESH_H

#include "BaseSystem.h"
#include <shared/defs.h>
#include <util/util.h>
#include "animTcl.h"
#include <GLmodel/GLmodel.h>
#include "AABB.h"
#include "Eigen/Dense"

#include "shared/opengl.h"

struct Triangle {
	Eigen::Vector3d p1;
	Eigen::Vector3d p2;
	Eigen::Vector3d p3;
};

// a sample system
class Mesh : public BaseSystem
{

public:
	Mesh(const std::string& name, float size, char* mesh, int id);
	virtual void getState(double* p);
	virtual void setState(double* p);
	void reset(double time);

	double minXYZ(double arr[3]);
	double maxXYZ(double arr[3]);

	void setupPrimitives();

	void createBVH();
	void updateBVH();

	double barrierWithMesh(Mesh* mesh, Eigen::Vector3d trans);
	double minDistance(Mesh* mesh, Eigen::Vector3d trans);

	void setTranslation(double* p);
	void setRotation(double* p);

	void getTranslation(double* p);
	void getRotation(double* p);

	void display(GLenum mode = GL_RENDER);

	void readModel(char* fname) { m_model.ReadOBJ(fname); }
	void flipNormals(void) { glmReverseWinding(&m_model); }
	int command(int argc, myCONST_SPEC char** argv);

	GLMmodel m_model;
	aabb::Tree bvhTree;

	Eigen::Vector3d* vertices;
	Eigen::Vector3d* init_vertices;
	int num_vertices;

	Triangle* triangles;
	Triangle* init_triangles;
	int num_triangles;

	Eigen::Vector3d m_pos;
	Eigen::Vector3d m_rot;

	int id;

protected:

	float m_sx;
	float m_sy;
	float m_sz;

	float size;

	double dhat = 0.05;
};
#endif
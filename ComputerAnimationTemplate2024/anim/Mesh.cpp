#include "Mesh.h"

using namespace Eigen;

Mesh::Mesh(const std::string& name, float size, char* mesh, int id) :
	BaseSystem(name),
	m_sx(1.0f),
	m_sy(1.0f),
	m_sz(1.0f),
	size(size),
	id(id)
{
	m_model.ReadOBJ(mesh);
	m_pos << 0.0, 0.0, 0.0;
	m_rot << 0.0, 0.0, 0.0;

	setupPrimitives();

}

void Mesh::getState(double* p)
{

	//VecCopy(p, m_pos);

}

void Mesh::setState(double* p)
{

	//VecCopy(m_pos, p);

}

void Mesh::reset(double time)
{

	//setVector(m_pos, 0, 0, 0);

}


int Mesh::command(int argc, myCONST_SPEC char** argv)
{
	if (argc < 1)
	{
		animTcl::OutputMessage("system %s: wrong number of params.", m_name.c_str());
		return TCL_ERROR;
	}
	else if (strcmp(argv[0], "read") == 0)
	{
		if (argc == 2)
		{
			m_model.ReadOBJ(argv[1]);
			glmFacetNormals(&m_model);
			glmVertexNormals(&m_model, 90);
			return TCL_OK;
		}
		else
		{
			animTcl::OutputMessage("Usage: read <file_name>");
			return TCL_ERROR;
		}
	}
	else if (strcmp(argv[0], "scale") == 0)
	{
		if (argc == 4)
		{
			m_sx = (float)atof(argv[1]);
			m_sy = (float)atof(argv[2]);
			m_sz = (float)atof(argv[3]);
		}
		else
		{
			animTcl::OutputMessage("Usage: scale <sx> <sy> <sz> ");
			return TCL_ERROR;

		}
	}
	else if (strcmp(argv[0], "pos") == 0)
	{
		if (argc == 4)
		{
			m_pos[0] = atof(argv[1]);
			m_pos[1] = atof(argv[2]);
			m_pos[2] = atof(argv[3]);
		}
		else
		{
			animTcl::OutputMessage("Usage: pos <x> <y> <z> ");
			return TCL_ERROR;

		}
	}
	else if (strcmp(argv[0], "flipNormals") == 0)
	{
		flipNormals();
		return TCL_OK;

	}
	else if (strcmp(argv[0], "reset") == 0)
	{
		double p[3] = { 0,0,0 };
		setState(p);
	}

	glutPostRedisplay();
	return TCL_OK;

}	// SampleParticle::command


double Mesh::minXYZ(double arr[3]) {
	return min(arr[0], min(arr[1], arr[2]));
}

double Mesh::maxXYZ(double arr[3]) {
	return max(arr[0], max(arr[1], arr[2]));
}


void Mesh::setupPrimitives() {
	num_vertices = m_model.numvertices;
	vertices = new Vector3d[num_vertices];
	init_vertices = new Vector3d[num_vertices];

	std::vector<double> vecs;
	for (int i = 0; i < num_vertices * 3; i++) {
		vecs.push_back(m_model.vertices[i]);
	}

	for (int v = 0; v < num_vertices; v++) {
		/* Base coords */
		init_vertices[v] <<
			m_model.vertices[(v * 3) + 3],
			m_model.vertices[(v * 3) + 4],
			m_model.vertices[(v * 3) + 5];

		/* Translated coords */
		vertices[v][0] = init_vertices[v][0];
		vertices[v][1] = init_vertices[v][1];
		vertices[v][2] = init_vertices[v][2];
	}

	num_triangles = m_model.numtriangles;
	triangles = new Triangle[num_triangles];
	init_triangles = new Triangle[num_triangles];
	int v1;
	int v2;
	int v3;

	for (int t = 0; t < num_triangles; t++) {
		v1 = m_model.triangles[t].vindices[0];
		v2 = m_model.triangles[t].vindices[1];
		v3 = m_model.triangles[t].vindices[2];

		/* Base coords */
		init_triangles[t].p1 <<
			vertices[v1 - 1][0],
			vertices[v1 - 1][1],
			vertices[v1 - 1][2];

		init_triangles[t].p2 <<
			vertices[v2 - 1][0],
			vertices[v2 - 1][1],
			vertices[v2 - 1][2];

		init_triangles[t].p3 <<
			vertices[v3 - 1][0],
			vertices[v3 - 1][1],
			vertices[v3 - 1][2];

		/* Translated coords */
		triangles[t].p1[0] = init_triangles[t].p1[0];
		triangles[t].p1[1] = init_triangles[t].p1[1];
		triangles[t].p1[2] = init_triangles[t].p1[2];

		triangles[t].p2[0] = init_triangles[t].p2[0];
		triangles[t].p2[1] = init_triangles[t].p2[1];
		triangles[t].p2[2] = init_triangles[t].p2[2];

		triangles[t].p3[0] = init_triangles[t].p3[0];
		triangles[t].p3[1] = init_triangles[t].p3[1];
		triangles[t].p3[2] = init_triangles[t].p3[2];

	}
}

void Mesh::createBVH() {
	bvhTree = aabb::Tree(3, dhat/2.0, num_vertices + num_triangles, true);

	double radius = dhat;

	for (int v = 0; v < num_vertices; v++) {
		/* Create small bounding boxes around each vertex */
		std::vector<double> lowerBound({
				vertices[v][0] - radius,
				vertices[v][1] - radius,
				vertices[v][2] - radius
			});

		std::vector<double> upperBound({
				vertices[v][0] + radius,
				vertices[v][1] + radius,
				vertices[v][2] + radius
			});

		bvhTree.insertParticle(v, lowerBound, upperBound);
	}

	for (int t = 0; t < num_triangles; t++) {
		/*
		* Create bounding box from triangle.
		* Lower bound = min x,y,z from all points
		* Upper bound = max x,y,z from all points
		*/
		double x[3] = { triangles[t].p1[0], triangles[t].p2[0], triangles[t].p3[0] };
		double y[3] = { triangles[t].p1[1], triangles[t].p2[1], triangles[t].p3[1] };
		double z[3] = { triangles[t].p1[2], triangles[t].p2[2], triangles[t].p3[2] };

		std::vector<double> lowerBound({
				minXYZ(x),
				minXYZ(y),
				minXYZ(z)
			});

		std::vector<double> upperBound({
				maxXYZ(x),
				maxXYZ(y),
				maxXYZ(z)
			});

		bvhTree.insertParticle(num_vertices + t, lowerBound, upperBound);
	}
}

void Mesh::updateBVH() {
	double radius = dhat;

	for (int v = 0; v < num_vertices; v++) {
		/* Create small bounding boxes around each vertex */
		std::vector<double> lowerBound({
				vertices[v][0] - radius,
				vertices[v][1] - radius,
				vertices[v][2] - radius
			});

		std::vector<double> upperBound({
				vertices[v][0] + radius,
				vertices[v][1] + radius,
				vertices[v][2] + radius
			});

		bvhTree.updateParticle(v, lowerBound, upperBound);
	}

	for (int t = 0; t < num_triangles; t++) {
		/*
		* Create bounding box from triangle.
		* Lower bound = min x,y,z from all points
		* Upper bound = max x,y,z from all points
		*/
		double x[3] = { triangles[t].p1[0], triangles[t].p2[0], triangles[t].p3[0] };
		double y[3] = { triangles[t].p1[1], triangles[t].p2[1], triangles[t].p3[1] };
		double z[3] = { triangles[t].p1[2], triangles[t].p2[2], triangles[t].p3[2] };

		std::vector<double> lowerBound({
				minXYZ(x),
				minXYZ(y),
				minXYZ(z)
			});

		std::vector<double> upperBound({
				maxXYZ(x),
				maxXYZ(y),
				maxXYZ(z)
			});

		bvhTree.updateParticle(num_vertices + t, lowerBound, upperBound);
	}
}

double b(double d, double dhat) {

	if (d >= dhat) {
		return 0;
	}
	else {
		return -pow(d - dhat, 2) * log(d / dhat);
	}
}

double Mesh::barrierWithMesh(Mesh* mesh, Vector3d trans) {
	double barrier = 0.0;

	double t[3] = { trans[0], trans[1], trans[2] };
	setTranslation(t);

	/*double r[3] = { rot[0], rot[1], rot[2] };
	setRotation(r);*/

	/* Points ------------------------------------------------------- */
	for (int v = 0; v < num_vertices; v++) {
		aabb::AABB box = bvhTree.getAABB(v);
		std::vector<unsigned int> contacts = mesh->bvhTree.query(box);

		for (int col = 0; col < contacts.size(); col++) {
			int ind = contacts[col];

			// Point Point
			if (ind < num_triangles) {
				Vector3d sub = vertices[v] - mesh->vertices[ind];

				// animTcl::OutputMessage("POINT POINT");

				double dist = sub.norm();
				if (dist == 0.0)
					dist = 1.0e-20;

				barrier += b(sub.norm(), dhat);
			}

			// Point Triangle
			else {
				// ind is the index of the bvh tree, not the triangles array
				ind -= num_vertices;

				Vector3d sub = vertices[v] - mesh->triangles[ind].p1;
				Vector3d cross =
					(mesh->triangles[ind].p2 - mesh->triangles[ind].p1).cross(mesh->triangles[ind].p3 - mesh->triangles[ind].p1);

				// animTcl::OutputMessage("POINT TRIANGLE");

				double dist = abs(sub.dot(cross.normalized()));
				if (dist == 0.0)
					dist = 1.0e-20;

				barrier += b(dist, dhat);
			}
		}
	}

	/* Triangles ----------------------------------------------------- */
	for (int t = 0; t < num_triangles; t++) {
		aabb::AABB box = bvhTree.getAABB(t + num_vertices);
		std::vector<unsigned int> contacts = mesh->bvhTree.query(box);

		for (int col = 0; col < contacts.size(); col++) {
			int ind = contacts[col];

			// Triangle Point
			if (ind < mesh->num_vertices) {
				Vector3d sub = mesh->vertices[ind] - triangles[t].p1;
				Vector3d cross =
					(triangles[t].p2 - mesh->vertices[ind]).cross(triangles[t].p3 - mesh->vertices[ind]);

				//animTcl::OutputMessage("TRIANGLE POINT");

				double dist = abs(sub.dot(cross.normalized()));
				if (dist == 0.0)
					dist = 1.0e-20;

				barrier += b(dist, dhat);
			}
		}
	}

	return barrier;
}

double Mesh::minDistance(Mesh* mesh, Vector3d trans) {

	double min_dist = DBL_MAX;

	double t[3] = { trans[0], trans[1], trans[2] };
	setTranslation(t);

	/* Points ------------------------------------------------------- */
	for (int v = 0; v < num_vertices; v++) {
		aabb::AABB box = bvhTree.getAABB(v);
		std::vector<unsigned int> contacts = mesh->bvhTree.query(box);

		for (int col = 0; col < contacts.size(); col++) {
			int ind = contacts[col];

			// Point Point
			if (ind < num_triangles) {
				Vector3d sub = vertices[v] - mesh->vertices[ind];

				// animTcl::OutputMessage("POINT POINT");

				double dist = sub.norm();
				if (dist < min_dist) {
					min_dist = dist;
				}
			}

			// Point Triangle
			else {
				// ind is the index of the bvh tree, not the triangles array
				ind -= num_vertices;

				Vector3d sub = vertices[v] - mesh->triangles[ind].p1;
				Vector3d cross =
					(mesh->triangles[ind].p2 - mesh->triangles[ind].p1).cross(mesh->triangles[ind].p3 - mesh->triangles[ind].p1);

				// animTcl::OutputMessage("POINT TRIANGLE");

				double dist = abs(sub.dot(cross.normalized()));
				if (dist < min_dist) {
					min_dist = dist;
				}
			}
		}
	}

	/* Triangles ----------------------------------------------------- */
	for (int t = 0; t < num_triangles; t++) {
		aabb::AABB box = bvhTree.getAABB(t + num_vertices);
		std::vector<unsigned int> contacts = mesh->bvhTree.query(box);

		for (int col = 0; col < contacts.size(); col++) {
			int ind = contacts[col];

			// Triangle Point
			if (ind < mesh->num_vertices) {
				Vector3d sub = mesh->vertices[ind] - triangles[t].p1;
				Vector3d cross =
					(triangles[t].p2 - mesh->vertices[ind]).cross(triangles[t].p3 - mesh->vertices[ind]);

				//animTcl::OutputMessage("TRIANGLE POINT");

				double dist = abs(sub.dot(cross.normalized()));
				if (dist < min_dist) {
					min_dist = dist;
				}
			}
		}
	}

	return min_dist;
}

void Mesh::setTranslation(double* p) {
	m_pos(0) = p[0];
	m_pos(1) = p[1];
	m_pos(2) = p[2];

	for (int v = 0; v < num_vertices; v++) {
		vertices[v] = init_vertices[v] + m_pos;
	}

	for (int t = 0; t < num_triangles; t++) {
		triangles[t].p1 = init_triangles[t].p1 + m_pos;
		triangles[t].p2 = init_triangles[t].p2 + m_pos;
		triangles[t].p3 = init_triangles[t].p3 + m_pos;
	}

	updateBVH();
}

void Mesh::setRotation(double* p) {
	m_rot(0) = p[0];
	m_rot(1) = p[1];
	m_rot(2) = p[2];

	double theta = m_rot.norm() * PI / 180;

	if (theta == 0.0)
		return;

	Vector3d e = m_rot.normalized();

	for (int v = 0; v < num_vertices; v++) {
		vertices[v] =
			(vertices[v] * cos(theta)) +
			(e.cross(vertices[v])) +
			((1 - cos(theta)) * e.dot(vertices[v]) * e);
	}

	for (int t = 0; t < num_triangles; t++) {
		triangles[t].p1 =
			(triangles[t].p1 * cos(theta)) +
			(e.cross(triangles[t].p1)) +
			((1 - cos(theta)) * e.dot(triangles[t].p1) * e);

		triangles[t].p2 =
			(triangles[t].p2 * cos(theta)) +
			(e.cross(triangles[t].p2)) +
			((1 - cos(theta)) * e.dot(triangles[t].p2) * e);

		triangles[t].p3 =
			(triangles[t].p3 * cos(theta)) +
			(e.cross(triangles[t].p3)) +
			((1 - cos(theta)) * e.dot(triangles[t].p3) * e);
	}

	updateBVH();
}

void Mesh::getTranslation(double* p) {
	p[0] = m_pos[0];
	p[1] = m_pos[1];
	p[2] = m_pos[2];
}

void Mesh::getRotation(double* p) {
	p[0] = m_rot[0];
	p[1] = m_rot[1];
	p[2] = m_rot[2];
}


void Mesh::display(GLenum mode)
{
	glEnable(GL_LIGHTING);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glEnable(GL_COLOR_MATERIAL);

	if (id == 0) {
		glColor3f(1.0, 0.0, 0.0);
	}
	else {
		glColor3f(1.0, 1.0, 1.0);
	}

	glTranslated(m_pos[0], m_pos[1], m_pos[2]);
	glRotated(m_rot.norm(), m_rot[0], m_rot[1], m_rot[2]);
	glScalef(m_sx, m_sy, m_sz);

	glmDraw(&m_model, GLM_SMOOTH);

	glPopMatrix();
	glPopAttrib();

}
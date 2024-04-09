#include "Ipc.h"

using namespace Eigen;
using namespace std;
using namespace fd;

Ipc::Ipc(const std::string& name, MeshCollection* collection) :
	BaseSimulator(name),
	collection(collection),
	prevTime(0.0),
	h(0.0),
	m(1.0),
	use_file(false)
{
	GlobalResourceManager::use()->setSimulationTime(0.01);

	// size = collection->getSize();
	size = 1;

	qt.resize(size, 3);
	qprev.resize(size, 3);
	qdot.resize(size, 3);
	fe.resize(size, 3);

	qt = collection->getTranslations();

	qprev = collection->getTranslations();

	qdot << 0.0, 0.0, 0.0;

	fe << 0.0, 0.0, 0.0;

	g << 0.0, -9.8, 0.0;

	/* File stuff */
	string filename = "../sim.txt";
	remove(filename.c_str());
}

Ipc::Ipc(const std::string& name, MeshCollection* collection, string file) :
	BaseSimulator(name),
	collection(collection),
	file(file),
	use_file(true),
	line(0)
{
	ifstream in(file);

	vector<Vector3d> coords_vec;

	double x, y, z;
	while (in >> x >> y >> z) {
		Vector3d temp;
		temp << x, y, z;

		coords_vec.push_back(temp);
	}

	in.close();

	line_nums = coords_vec.size();
	coords = new Vector3d[line_nums];

	for (int i = 0; i < line_nums; i++) {
		coords[i] = coords_vec[i];
	}

}

VectorXd Ipc::vecCopy(VectorXd v1) {
	VectorXd v2(size * 3);

	for (int i = 0; i < size * 3; i++) {
		v2(i) = v1(i);
	}

	return v2;
}

MatrixXd Ipc::matCopy(MatrixXd m1) {
	MatrixXd m2(size, 3);

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < 3; j++) {
			m2(i, j) = m1(i, j);
		}
	}

	return m2;
}

double Ipc::norm(VectorXd v) {
	return abs(v.maxCoeff());
}

VectorXd Ipc::flatten(MatrixXd m) {
	VectorXd v(size * 3);

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < 3; j++) {
			v((i * 3) + j) = m(i, j);
		}
	}

	return v;
}

Eigen::MatrixXd Ipc::unflatten(VectorXd v) {
	MatrixXd m(size, 3);

	for (int i = 0; i < size * 3; i++) {
		m((i / 3) + (i % 3)) = v(i);
	}

	return m;
}


double Ipc::E_q(VectorXd x) {
	double E = 0.0;

	for (int i = 0; i < size; i++) {
		Vector3d qt_i;
		Vector3d qdot_i;
		Vector3d fe_i;
		Vector3d q_i;
		qt_i << qt(i, 0), qt(i, 1), qt(i, 2);
		qdot_i << qdot(i, 0), qdot(i, 1), qdot(i, 2);
		fe_i << fe(i, 0), fe(i, 1), fe(i, 2);
		q_i << x((i * 3) + 0), x((i * 3) + 1), x((i * 3) + 2);

		Vector3d q_tilde = qt_i + (h * qdot_i) + (pow(h, 2) * (g + (fe_i/m)));

		E += (0.5 * m * q_i.dot(q_i)) - (m * q_i.dot(q_tilde));
	}

	Vector3d trans;
	trans << x(0), x(1), x(2);

	E += 1.0e10 * collection->getBarrier(trans);

	return E;
}


int Ipc::step(double time)
{
	if (use_file) {

		if (line < line_nums) {
			MatrixXd pos(1, 3);
			pos(0, 0) = coords[line][0];
			pos(0, 1) = coords[line][1];
			pos(0, 2) = coords[line][2];

			collection->updatePos(pos);

			line++;
		}

		return 0;
	}

	h = (time - prevTime);

	for (int i = 0; i < qt.rows(); i++) {
		qdot = (1 / h) * (qt - qprev);
	}

	qprev = matCopy(qt);

	/* Newton */
	VectorXd grad(size*3);
	MatrixXd hess(size*3, size*3);
	MatrixXd hess_inv(size * 3, size * 3);
	VectorXd p(size*3);
	VectorXd x(size*3);
	VectorXd xprev(size*3);

	x = flatten(qt);
	xprev = vecCopy(x);

	double offset = 1.0;
	int iterations = 0;

	do {
		auto E = [&](VectorXd x) {
			return this->E_q(x);
		};

		finite_gradient(x, E, grad);
		finite_hessian(x, E, hess);

		if (hess.determinant() == 0.0 || iterations > 50) {
			hess += offset * MatrixXd::Identity(size*3, size*3);
			offset *= 10;
			iterations = 0;
		}

		hess_inv = hess.inverse();

		p = -1 * hess.inverse() * grad;

		x = xprev + p;

		xprev = vecCopy(x);

		iterations++;

	} while ((1 / h) * norm(p) > 1.0e-10);

	animTcl::OutputMessage("[%f %f %f]", x(0), x(1), x(2));
	animTcl::OutputMessage("e: %f", E_q(x));

	ofstream outfile;
	outfile.open("../sim.txt", ios_base::app); // append instead of overwrite
	outfile << x(0) << " " << x(1) << " " << x(2) << "\n";
	outfile.close();

	qt = unflatten(x);

	collection->updatePos(qt);

	Vector3d trans;
	trans << qt(0, 0), qt(0, 1), qt(0, 2);
	animTcl::OutputMessage("dist: %f", collection->getDistance(trans));

	prevTime = time;

	return 0;

}

#include "Ipc.h"

using namespace Eigen;
using namespace std;
using namespace fd;

Ipc::Ipc(const std::string& name, MeshCollection* collection) :
	BaseSimulator(name),
	collection(collection),
	prevTime(0.0),
	h(0.0),
	m(1.0)
{
	GlobalResourceManager::use()->setSimulationTime(0.01);

	size = collection->getSize();

	qt.resize(size, 3);
	qprev.resize(size, 3);
	qdot.resize(size, 3);
	fe.resize(size, 3);

	qt << 1.0, 2.5, 0.0,
		  0.0, 0.0, 0.0;

	qprev << 1.0, 2.5, 0.0,
		     0.0, 0.0, 0.0;

	qdot << 0.0, 0.0, 0.0,
		    0.0, 0.0, 0.0;

	fe << 0.0, 0.0, 0.0,
		  0.0, 9.8, 0.0;

	g << 0.0, -9.8, 0.0;
}

double Ipc::E_q(VectorXd q) {
	double E = 0.0;

	for (int i = 0; i < qt.rows(); i++) {
		Vector3d qt_i;
		Vector3d qdot_i;
		Vector3d fe_i;
		Vector3d q_i;
		qt_i << qt(i, 0), qt(i, 1), qt(i, 2);
		qdot_i << qdot(i, 0), qdot(i, 1), qdot(i, 2);
		fe_i << fe(i, 0), fe(i, 1), fe(i, 2);
		q_i << q((i * 3) + 0), q((i * 3) + 1), q((i * 3) + 2);

		Vector3d q_tilde = qt_i + (h * qdot_i) + (pow(h, 2) * (g + (fe_i/m)));

		E += (0.5 * m * q_i.dot(q_i)) - (m * q_i.dot(q_tilde));

		E += 1000.0 * collection->getBarrier();
	}

	return E;
}


void vecCopy(VectorXd v1, VectorXd &v2) {
	for (int i = 0; i < v1.size(); i++) {
		v2(i) = v1(i);
	}
}


double norm(VectorXd v) {
	return abs(v.maxCoeff());
}


int Ipc::step(double time)
{
	h = (time - prevTime);

	for (int i = 0; i < qt.rows(); i++) {
		qdot(i, 0) = (1 / h) * (qt(i, 0) - qprev(i, 0));
		qdot(i, 1) = (1 / h) * (qt(i, 1) - qprev(i, 1));
		qdot(i, 2) = (1 / h) * (qt(i, 2) - qprev(i, 2));
	}

	/* Newton */
	VectorXd grad(size*3);
	MatrixXd hess(size*3, size*3);
	MatrixXd hess_inv(size * 3, size * 3);
	VectorXd p(size*3);
	VectorXd x(size*3);
	VectorXd xprev(size*3);

	x = flatten(qt);
	vecCopy(x, xprev);

	do {
		auto E = [&](VectorXd q) {
			return this->E_q(q);
		};

		finite_gradient(x, E, grad);
		finite_hessian(x, E, hess);

		hess_inv = hess.inverse();

		p = -1 * hess.inverse() * grad;

		x = xprev + p;


		vecCopy(x, xprev);

	} while ((1 / h) * norm(p) > 0.01);

	animTcl::OutputMessage("[%f %f %f]", x(0), x(1), x(2));
	animTcl::OutputMessage("b: %f", collection->getBarrier());

	qt = unflatten(x, qt.cols());
	qprev = unflatten(xprev, qprev.cols());

	collection->updatePos(qt);

	prevTime = time;

	return 0;

}

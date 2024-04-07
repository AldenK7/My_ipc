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

	Qt.resize(size, 3);
	Qprev.resize(size, 3);
	Qdot.resize(size, 3);

	fe.resize(size, 3);

	/* Translations */
	qt = collection->getTranslations();

	qprev = collection->getTranslations();

	qdot << 0.0, 0.0, 0.0,
		    0.0, 0.0, 0.0;

	/* Rotation */
	Qt = collection->getRotations();

	Qprev = collection->getRotations();

	Qdot << 0.0, 0.0, 0.0,
		    0.0, 0.0, 0.0;

	/* Forces */
	fe << 15.0, 15.0, 0.0,
		  0.0, 9.8*m, 0.0;

	g << 0.0, -9.8, 0.0;


	int h = 2;
	int d = 2;
	int w = 2;

	I << m * (pow(h, 2) + pow(d, 2)) / 12,
		 m * (pow(w, 2) + pow(d, 2)) / 12,
		 m * (pow(w, 2) + pow(h, 2)) / 12;

	J << -I[0] + I[1] + I[2], 0, 0,
		 0, I[0] - I[1] + I[2], 0,
		 0, 0, I[0] + I[1] - I[2];

	J /= 2.0;
}

double Ipc::E_q(VectorXd x) {
	double E = 0.0;

	for (int i = 0; i < size; i++) {
		/* Translation */
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

		/* Rotation */
		Vector3d vQt_i;
		Vector3d vQdot_i;
		Vector3d vQ_i;
		vQt_i << Qt(i, 0), Qt(i, 1), Qt(i, 2);
		vQdot_i << Qdot(i, 0), Qdot(i, 1), Qdot(i, 2);
		vQ_i << x((i * 3) + 3 + 0), x((i * 3) + 3 + 1), x((i * 3) + 3 + 2);

		Matrix3d Qt_i = rodrigues(vQt_i);
		Matrix3d Qdot_i = rodrigues(vQdot_i);
		Matrix3d Q_i = rodrigues(vQ_i);

		// Add torque at end here
		//Matrix3d Q_tilde = Qt_i + (h * Qdot_i) + (pow(h, 2) * J.inverse());
		Matrix3d Q_tilde = Qt_i + (h * Qdot_i);

		E += (0.5 * (Q_i * J * Q_i.transpose()).trace()) - (Q_i * J * Q_tilde.transpose()).trace();
	}

	Vector3d trans;
	trans << x(0), x(1), x(2);

	Vector3d rot;
	rot << x(3), x(4), x(5);

	/* Barrier */
	// E += 10.0 * collection->getBarrier(trans, rot);

	return E;
}



double sinc(double x) {
	if (x == 0.0) {
		return 1.0;
	}
	else {
		return (sin(x) / x);
	}
}


Matrix3d skew(Vector3d v) {
	Matrix3d skewed;
	skewed <<
		0.0, -v[2], v[1],
		v[2], 0.0, -v[0],
		-v[1], v[0], 0.0;

	return skewed;
}


Matrix3d Ipc::rodrigues(Vector3d theta) {
	Matrix3d rot_mat;
	Matrix3d skewed = skew(theta);

	rot_mat = Matrix3d::Identity() + 
			  (sinc(theta.norm()*PI/180) * skewed) +
			  (2 * pow(sinc(theta.norm()*PI/180/2), 2) * skewed * skewed);
	
	return rot_mat;
}


void vecCopy(VectorXd v1, VectorXd &v2) {
	for (int i = 0; i < v1.size(); i++) {
		v2(i) = v1(i);
	}
}


void matCopy(MatrixXd m1, MatrixXd &m2) {
	for (int i = 0; i < m1.rows(); i++) {
		for (int j = 0; j < m2.cols(); j++) {
			m2(i, j) = m1(i, j);
		}
	}
}


double norm(VectorXd v) {
	return abs(v.maxCoeff());
}


/*
* < t1_x, t1_y, t1_z, r1_x, r1_y, r1_z, ... , tn_x, ... >
*/
VectorXd flatten(MatrixXd trans, MatrixXd rot) {
	assert(trans.rows() == rot.rows());

	VectorXd flat((trans.rows() + rot.rows()) * 3);

	for (int row = 0; row < trans.rows(); row++) {
		for (int col = 0; col < 3; col++) {
			flat((row * 6) + col) = trans(row, col);
			flat((row * 6) + 3 + col) = rot(row, col);
		}
	}

	return flat;
}


void unflatten(VectorXd x, MatrixXd& t, MatrixXd& r) {

	/*for (int i = 0; i < x.size(); i++) {
		animTcl::OutputMessage("x [%f]", x[i]);
	}*/

	for (int i = 0; i < x.size()/6; i++) {
		t(i, 0) = x[(i * 6)];
		t(i, 1) = x[(i * 6) + 1];
		t(i, 2) = x[(i * 6) + 2];

		r(i, 0) = x[(i * 6) + 3];
		r(i, 1) = x[(i * 6) + 4];
		r(i, 2) = x[(i * 6) + 5];
	}

	/*for (int i = 0; i < 2; i++) {
		animTcl::OutputMessage("Q [%.2f %.2f %.2f]",
			r(i, 0), r(i, 1), r(i, 2));
	}*/
}


int Ipc::step(double time)
{
	h = (time - prevTime);
	h = 0.1;

	for (int i = 0; i < qt.rows(); i++) {
		qdot = (1 / h) * (qt - qprev);
		Qdot = (1 / h) * (Qt - Qprev);
	}

	matCopy(qt, qprev);
	matCopy(Qt, Qprev);

	/* Newton */
	VectorXd grad(size*6);
	MatrixXd hess(size*6, size*6);
	MatrixXd hess_inv(size*6, size*6);
	VectorXd p(size*6);
	VectorXd x(size*6);
	VectorXd xprev(size*6);

	x = flatten(qt, Qt);
	vecCopy(x, xprev);

	double offset = 10.0;
	int iterations = 0;

	do {
		auto E = [&](VectorXd x) {
			return this->E_q(x);
		};

		finite_gradient(x, E, grad);
		finite_hessian(x, E, hess);
		
		if (hess.determinant() == 0.0 || iterations > 50) {
			hess += offset * MatrixXd::Identity(12, 12);
			offset *= 10;
		}

		// animTcl::OutputMessage("[%f]", hess.determinant());

		
		/*for (int i = 0; i < 12; i++) {	
			animTcl::OutputMessage("[%.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f]",
				hess(i, 0), hess(i, 1), hess(i, 2), hess(i, 3), hess(i, 4), hess(i, 5), hess(i, 6), 
				hess(i, 7), hess(i, 8), hess(i, 9), hess(i, 10), hess(i, 11));
		}*/

		hess_inv = hess.inverse();

		p = -1 * hess.inverse() * grad;

		x = xprev + p;

		vecCopy(x, xprev);

		iterations++;

	} while ((1 / h) * norm(p) > 1.0e-10);

	unflatten(x, qt, Qt);
	// unflatten(xprev, qprev, Qprev);

	for (int i = 0; i < 2; i++) {
		animTcl::OutputMessage("q [%.4f %.4f %.4f]",
			qt(i, 0), qt(i, 1), qt(i, 2));
	}

	for (int i = 0; i < 2; i++) {
		animTcl::OutputMessage("Q [%.4f %.4f %.4f]",
			Qt(i, 0), Qt(i, 1), Qt(i, 2));
	}

	animTcl::OutputMessage("e: %f", E_q(x));

	collection->updatePos(qt);
	collection->updateRot(Qt);

	prevTime = time;

	return 0;

}

#ifndef MY_IPC_H
#define MY_IPC_H

#include <GLModel/GLModel.h>
#include <shared/defs.h>
#include <util/util.h>
#include "animTcl.h"
#include "BaseSimulator.h"
#include "Mesh.h"
#include "MeshCollection.h"
#include "finitediff.h"
#include "GlobalResourceManager.h"
#include "Eigen/Dense"

#include <string>


class Ipc : public BaseSimulator
{
public:

	Ipc(const std::string& name, MeshCollection* collection);

	double E_q(Eigen::VectorXd x);

	Eigen::Matrix3d rodrigues(Eigen::Vector3d);

	int step(double time);
	int init(double time)
	{
		return 0;
	};

	int command(int argc, myCONST_SPEC char** argv) { return TCL_OK; }

protected:

	Eigen::Vector3d g;

	double prevTime;
	int size;

	double h;
	double dhat;
	double m;

	Eigen::MatrixXd qt;
	Eigen::MatrixXd qprev;
	Eigen::MatrixXd qdot;

	Eigen::MatrixXd Qt;
	Eigen::MatrixXd Qprev;
	Eigen::MatrixXd Qdot;

	Eigen::MatrixXd fe;

	Eigen::Vector3d I;
	Eigen::Matrix3d J;

	MeshCollection* collection;
};


#endif
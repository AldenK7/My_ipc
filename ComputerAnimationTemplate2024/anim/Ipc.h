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
#include "Eigen/Core"

#include <string>


class Ipc : public BaseSimulator
{
public:

	Ipc(const std::string& name, MeshCollection* collection);

	double E_q(Eigen::VectorXd q);

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
	Eigen::MatrixXd fe;

	MeshCollection* collection;
};


#endif
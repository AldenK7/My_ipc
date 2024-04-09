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

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdio>

class Ipc : public BaseSimulator
{
public:

	Ipc(const std::string& name, MeshCollection* collection);
	Ipc(const std::string& name, MeshCollection* collection, std::string file);

	double E_q(Eigen::VectorXd x);

	Eigen::VectorXd vecCopy(Eigen::VectorXd v1);
	Eigen::MatrixXd matCopy(Eigen::MatrixXd m1);
	double norm(Eigen::VectorXd v);
	Eigen::VectorXd flatten(Eigen::MatrixXd m);
	Eigen::MatrixXd unflatten(Eigen::VectorXd m);

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

	bool use_file;
	std::string file;
	Eigen::Vector3d* coords;
	int line;
	int line_nums;
};


#endif
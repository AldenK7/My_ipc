////////////////////////////////////////////////////
// // Template code for  CSC 473
////////////////////////////////////////////////////

#ifdef WIN32
#include <windows.h>
#endif


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>

#include <shared/defs.h>

#include "shared/opengl.h"

#include <string.h>
#include <util/util.h>
#include <GLModel/GLModel.h>
#include "anim.h"
#include "animTcl.h"
#include "myScene.h"
#include "SampleParticle.h"
#include "SampleGravitySimulator.h"
#include "Mesh.h"
#include "Ipc.h"
#include "MeshCollection.h"

//#include <util/jama/tnt_stopwatch.h>
//#include <util/jama/jama_lu.h>

// register a sample variable with the shell.
// Available types are:
// - TCL_LINK_INT 
// - TCL_LINK_FLOAT

int g_testVariable = 10;

SETVAR myScriptVariables[] = {
	"testVariable", TCL_LINK_INT, (char *) &g_testVariable,
	"",0,(char *) NULL
};


//---------------------------------------------------------------------------------
//			Hooks that are called at appropriate places within anim.cpp
//---------------------------------------------------------------------------------

// start or end interaction
void myMouse(int button, int state, int x, int y)
{

	// let the global resource manager know about the new state of the mouse 
	// button
	GlobalResourceManager::use()->setMouseButtonInfo( button, state );

	if( button == GLUT_LEFT_BUTTON && state == GLUT_DOWN )
	{
		animTcl::OutputMessage(
			"My mouse received a mouse button press event\n");

	}
	if( button == GLUT_LEFT_BUTTON && state == GLUT_UP )
	{
		animTcl::OutputMessage(
			"My mouse received a mouse button release event\n") ;
	}
}	// myMouse

// interaction (mouse motion)
void myMotion(int x, int y)
{

	GLMouseButtonInfo updatedMouseButtonInfo = 
		GlobalResourceManager::use()->getMouseButtonInfo();

	if( updatedMouseButtonInfo.button == GLUT_LEFT_BUTTON )
	{
		animTcl::OutputMessage(
			"My mouse motion callback received a mousemotion event\n") ;
	}

}	// myMotion


void MakeScene(void)
{

	bool success;

	int scene = 0;

	MeshCollection* collection;
	Mesh* obj1;
	Mesh* obj2;

	Vector pos1;
	Vector pos2;

	switch (scene) {
		case 0:
			obj1 = new Mesh("1", 1.0, "data/cube.obj", 0);
			success = GlobalResourceManager::use()->addSystem(obj1, true);
			assert(success);

			obj1->setupPrimitives();

			pos1;
			setVector(pos1, 0.0, 2.0, 0.0);

			obj1->createBVH();
			obj1->setTranslation(pos1);

			obj2 = new Mesh("2", 1.0, "data/floor.obj", 1);
			success = GlobalResourceManager::use()->addSystem(obj2, true);
			assert(success);

			pos2;
			setVector(pos2, 0.0, -0.9, 0.0);

			obj2->createBVH();
			obj2->setTranslation(pos2);

			collection = new MeshCollection(2);
			collection->addMesh(obj1);
			collection->addMesh(obj2);

			break;

		case 1:
			obj1 = new Mesh("1", 1.0, "data/cone.obj", 0);
			success = GlobalResourceManager::use()->addSystem(obj1, true);
			assert(success);

			obj1->setupPrimitives();

			pos1;
			setVector(pos1, 0.0, 2.0, 0.0);

			obj1->createBVH();
			obj1->setTranslation(pos1);

			obj2 = new Mesh("2", 1.0, "data/torus_small.obj", 1);
			success = GlobalResourceManager::use()->addSystem(obj2, true);
			assert(success);

			pos2;
			setVector(pos2, 0.0, -0.9, 0.0);

			obj2->createBVH();
			obj2->setTranslation(pos2);

			collection = new MeshCollection(2);
			collection->addMesh(obj1);
			collection->addMesh(obj2);

			break;

		case 2:
			obj1 = new Mesh("1", 1.0, "data/cone.obj", 0);
			success = GlobalResourceManager::use()->addSystem(obj1, true);
			assert(success);

			obj1->setupPrimitives();

			pos1;
			setVector(pos1, 0.0, 1.0, 0.0);

			obj1->createBVH();
			obj1->setTranslation(pos1);

			obj2 = new Mesh("2", 1.0, "data/floor.obj", 1);
			success = GlobalResourceManager::use()->addSystem(obj2, true);
			assert(success);

			pos2;
			setVector(pos2, 0.0, -0.9, 0.0);

			obj2->createBVH();
			obj2->setTranslation(pos2);

			collection = new MeshCollection(2);
			collection->addMesh(obj1);
			collection->addMesh(obj2);

			break;

		case 3:
			obj1 = new Mesh("1", 1.0, "data/cube.obj", 0);
			success = GlobalResourceManager::use()->addSystem(obj1, true);
			assert(success);

			obj1->setupPrimitives();

			pos1;
			setVector(pos1, 0.0, 1.0, 0.0);

			obj1->createBVH();
			obj1->setTranslation(pos1);

			obj2 = new Mesh("2", 1.0, "data/torus.obj", 1);
			success = GlobalResourceManager::use()->addSystem(obj2, true);
			assert(success);

			pos2;
			setVector(pos2, 0.0, -0.9, 0.0);

			obj2->createBVH();
			obj2->setTranslation(pos2);

			collection = new MeshCollection(2);
			collection->addMesh(obj1);
			collection->addMesh(obj2);

			break;
	}

	// Ipc* ipc = new Ipc("ipc", collection);
	Ipc* ipc = new Ipc("ipc", collection, "../sim.txt");
	success = GlobalResourceManager::use()->addSimulator(ipc);
	assert( success );

}	// MakeScene

// OpenGL initialization
void myOpenGLInit(void)
{
	animTcl::OutputMessage("Initialization routine was called.");

}	// myOpenGLInit

void myIdleCB(void)
{
	
	return;

}	// myIdleCB

void myKey(unsigned char key, int x, int y)
{
	 animTcl::OutputMessage("My key callback received a key press event\n");
	return;

}	// myKey

static int testGlobalCommand(ClientData clientData, Tcl_Interp *interp, int argc, myCONST_SPEC char **argv)
{
	 animTcl::OutputMessage("This is a test command!");
    animTcl::OutputResult("100") ;
	return TCL_OK;

}	// testGlobalCommand

void mySetScriptCommands(Tcl_Interp *interp)
{

	// here you can register additional generic (they do not belong to any object) 
	// commands with the shell

	Tcl_CreateCommand(interp, "test", testGlobalCommand, (ClientData) NULL,
					  (Tcl_CmdDeleteProc *)	NULL);

}	// mySetScriptCommands

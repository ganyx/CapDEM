//#include <sdtlib.h>
#include <iostream> //library for C++
#include <unistd.h>
//#include <for>
//#include <QList>	/**<Qt containers*/
//#include <QThread>  /**<Qt multicore*/

#include <omp.h>	// OpenMP
#define NTHREADS 1	// Number of CPUs for one job
#include <time.h>

#include <cstdlib>
#include <algorithm>

#include <sys/resource.h>

using namespace std; // for a easier use of cin and cout

#include <math.h>     /**<used for pow(),sqrt()*/
#include <sys/stat.h> /**<used for mkdir();*/

#include <dirent.h>   /**<used for opendir();*/
#include <string> /**<Used to handle string.*/
#include <fstream> /**<Used to open files.*/ 
#include <sstream> /**<Used to convert integer into string.*/

#define DIM 3  /**< Space dimension. */ 
#define PI 3.1415926535 
#define MESH_SIZE 2.0  /** ratio of mesh size and the maximum particle diameter */
#define INIT_BOND 0
#define BOND_STRENGTH 20.0
#define SHEAR_LIMIT 0.5
#define DILAT_LIMIT 0.0
#define WALL_KD	5e-4
#define WALL_DAMP 0.001
#define GLOBAL_DAMPING	1.0

//#define COMP_FRACTION 0.10 /**< Composite fraction, number fraction */
bool LIQUID_TRANSFER;
#define LIQUID_DENSITY 0.1	//exp
#define MAX_CAP_LENGTH 0.5	//exp
//#define CONTACT_ANGLE 3.1416/4.0
#define WATER_K 1000.0
#define AIR_K 0.01
#define MIN_SATURATION 0.01
#define MAX_SATURATION 0.995
#define MAX_SATURATION_AIR 0.95 // S_a
#define MIN_S 0.0
#define MAX_S 1.0
#define MAX_DWATER 0.05 // limit the max rate of water transport per contact
#define WATER_CONDUCTION_RATIO	0.05

#define CONTACT_ANGLE_MAX 	60.0/180.0*PI // PI/3.0
#define CONTACT_ANGLE_MIN 	5.0/180.0*PI // PI/18.0
#define INLET_OUTLET_PRESSURE	0.01	// Water inlet/outlet rate
#define WETTING_RATE	0.0005

bool Voronoi_Update;

#define PSEUDO_2D 0 /**< 1= yes 0 = no (then 3D)*/	     
	     
class Cvector;
class Cmatrix;
class Cwall;
class Cstat;
class Cpair_mobility_line;
class Cparticle;
class Ccontact;
class Cconfig;
class Cin_out;

class Cparameter;
class Cevent;
class Cpost_process;
class Crun;

class Ccell;
class Cgrid;
class Cnode;
class Cbox;
class Cmesh;
class Cprofile;

// Voronoi tesselation from Voro++
#define NMAX 100	// The maximum number of Voronoi neighbors
#include "voro++/voro++.cc"
//using namespace voro;

#include "matrix.h"
#include "inout.h"
#include "stat.h"
#include "contact.h"
#include "particle.h"
#include "parameter.h"
#include "cell.h"   
#include "mesh.h" 
#include "config.h"
#include "profile.h"

//#include "post_process/grid.h"
//#include "post_process/node.h"
//#include "post_process/distribution.h"
//#include "post_process/post_process.h"	     
#include "run.h"
#include "post_process.h"

#include "matrix.cpp"
#include "stat.cpp"
#include "contact.cpp"
#include "particle.cpp"
#include "inout.cpp"
#include "parameter.cpp"  
#include "cell.cpp"
#include "mesh.cpp"
#include "config.cpp"
#include "profile.cpp"
#include "run.cpp"
#include "post_process.cpp"
//#include "post_process/spatial.cpp"
//#include "post_process/grid.cpp"
//#include "post_process/node.cpp"
//#include "post_process/distribution.cpp"
//#include "post_process/post_process.cpp"	    
	
	
/*! \mainpage 

\author Yixiang Gan
\version CapDEM: DEM with Capillary interactions, based on Soft Dyanamics
\date    FEB-2013

 \section intro_sec Introduction
 This version of the CapDEM can simulate the motion, the temperature of grains, and water movement between grains, including:

  \li  elastic-frictional contact, normal viscous dissipation, rolling and twist torque and friction;
  \li  heat conduction through contacts;
  \li  shear flow;
  \li  flow down a slope under gravity.
  \li  capillary interaction between grains
  \li  water movement

\section licence_sec Licence

The CapDEM is a free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation.
   
    It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License in the file 'COPYING' along with the Soft Dynamics package. If not, see <http://www.gnu.org/licenses/>.

\section install_sec Installation-compilation

The CapDEM is written in C++ and include other GPL library (see file soft_dynamics.h). It uses the Qt library, and can be easily complile using the qmake facility, for instance through Qdevelop editor.
*/     

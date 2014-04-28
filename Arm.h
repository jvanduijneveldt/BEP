#ifndef ARM_H
#define ARM_H
#include <cstdlib>
#include <iostream>
#include "dlib\matrix.h"

using namespace std;
using namespace dlib;

int const NUMBER_OF_JOINTS = 4;
int const NUMBER_OF_FRAMES = NUMBER_OF_JOINTS + 2; 
int const LAST_FRAME = NUMBER_OF_FRAMES - 1;
enum joint_type {ROTATIONAL, PRISMATIC};

//delete in ros
typedef int C3mxlROS;
//

class Arm{
    private:
        // arm properies
        matrix<joint_type,1,NUMBER_OF_JOINTS>           E;              // prisatic/rotational configuration 
        double                                          v_Max;          // maximum linear speed of the end effector
        double                                          a_Max;          // maximum linear acceleration of the end effector
        matrix<double,NUMBER_OF_JOINTS,1>               omega_Max;      // maximum angular speed of each joint
        matrix<double,NUMBER_OF_JOINTS,1>               alpha_Max;      // maximum angular acceleration of each joint
        matrix<double,4,4>                              T01;            // Transformation from zeroth to first frame
        matrix<double,4,4>                              Tend;           // Transformation from second last to last frame
        matrix<double,NUMBER_OF_JOINTS,NUMBER_OF_JOINTS>jointDecoupling;// Decoupling of motor angles and joint angles
        matrix<double,NUMBER_OF_JOINTS,NUMBER_OF_JOINTS>motorDecoupling;// Decoupling of motor angles and joint angles
        matrix<C3mxlROS*,NUMBER_OF_JOINTS,1>            motorAdresses;  // adresses of each motor used to operate the joints
        
        //dynamic model
        matrix<double,3,NUMBER_OF_JOINTS>               C;  // Positions of centers of masses expressed in the frame of the respective joints
        matrix<double,1,NUMBER_OF_JOINTS>               mc; // Masses of each joint
        matrix<matrix<double,3,3>,1,NUMBER_OF_JOINTS>   Ic; // Rotational inertia of each joint expressed in the frame of the joint
        
        //state variables
        matrix<double,NUMBER_OF_JOINTS,4>                               DH;             // DH parameters describing the system of the system
        matrix<double,NUMBER_OF_JOINTS,1>                               q;              // Current state of the system
        matrix<matrix<double,4,4>,NUMBER_OF_FRAMES,NUMBER_OF_FRAMES>    Transform;      // Saved frame transformations
        matrix<bool,NUMBER_OF_FRAMES,NUMBER_OF_FRAMES>                  TransformExists;// matrix used to check if transforms need to be recalculated.
        
        //functions
        matrix<double,4,4>                                  getT(int a, int b); //returns the transformation matrix from frame b to frame a P_a  = Tab*P_b
        matrix<double,6,NUMBER_OF_JOINTS>                   getJ();             //returns the jacobian of the system
        matrix<double,3,NUMBER_OF_JOINTS>                   getJv(int i);       //returns the linear part of the jacobian 
        matrix<double,3,NUMBER_OF_JOINTS>                   getJw(int i);       //returns the rotational part of the jacobian
        matrix<double,3,NUMBER_OF_FRAMES>                   getO();             //returns the position of the origins of each frame expressed in frame 0
        matrix<double,3,1>                                  getO(int frame);    //returns the position of the origin of the selected frame expressed in frame 0
        matrix<double,3,NUMBER_OF_FRAMES>                   getZ();             //returns the coordinates of the Z vector of each frame expressed in frame 0
        matrix<double,3,1>                                  getZ(int frame);    //returns the coordinates of the Z vector of the selected frame expressed in frame 0 
        matrix<double,NUMBER_OF_JOINTS,NUMBER_OF_JOINTS>    getM();             //returns the mass matrix of the system
        matrix<double,3,3>                                  getIc(int joint);   //returns the rotational inertia of the selected joint expressed in frame 0
        matrix<double,3,1>                                  getC(int joint);    //returns the position of the center of mass of selected joint expressed in frame 0
        matrix<double,3,1>                                  cross(matrix<double,3,1> A , matrix<double,3,1> B); //returns the crossproduct between A and B
        void                                                updateDH();         //updates the DH variable the match the current configurationn
        void                                                moveTo(matrix<double,2,1>);
     
    public:
        Arm();                                              //constructor
        void init(matrix<C3mxlROS*,NUMBER_OF_JOINTS,1>);    //initialize variables for first usage
        
};
#endif

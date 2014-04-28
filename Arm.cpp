#include "Arm.h"

Arm::Arm(){

}

void Arm::init(matrix<C3mxlROS*,NUMBER_OF_JOINTS,1> newAdresses){
        
    // Define joint types (4 rotational joints)        
    E = ROTATIONAL,ROTATIONAL,ROTATIONAL,ROTATIONAL;
       
    // initialize DH parameters
    //      alpha   L       theta   d  
    DH=     0,      0,      0,      0,
            0,      0.6,    0,      0,
            0,      0.5,    0,      0,
            pi/2,   0,      0.2,    0; 
                    
    motorAdresses = newAdresses;
    Arm::updateDH();
    
     
    // no begin and end transforms
    T01=    1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1;
    
    Tend=   1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1;
    
    // motorangle= jointdecoupling * jointangles
    jointDecoupling =   1,0,0,0,
                        1,1,0,0,
                        0,0,1,0,
                        0,0,0,1;
    motorDecoupling = inv(jointDecoupling);                    
                        
    
    //set centre of mass positions 
    //these are expressed in the frame of their respective joints
    C= 0;
    // first joint
    C(0,0) = 0.25;
    C(1,0) =    0;
    C(2,0) =    0; 
           
    // second joint            
    C(0,1) =  0.2;
    C(1,1) =    0;
    C(2,1) =    0; 
                         
    // third joint            
    C(0,2) = 0.08;
    C(1,2) =    0;
    C(2,2) =    0; 
    // fourth joint            
    C(0,3) =  0.1;
    C(1,3) =    0;
    C(2,3) =    0; 
    
    // joint masses
    mc=1;
    
    //joint rotational inertia
    matrix<double,3,3> Itest;
    
    Itest  =    1,0,0,
                0,2,0,
                0,0,3;
    
    Ic(0)=Itest;
    Ic(1)=Itest;
    Ic(2)=Itest;
    Ic(3)=Itest;
    
    

    
    // No transforms known     
    TransformExists = false;
        
       


    
    
}


void Arm::moveTo(matrix<double,2,1> Target){
    //begin definitions
    matrix<double, 2, 2>
        J,
        subDecouple;
    matrix<double, 2, 1>
        O_effector,
        direction,
        v_vect,
        omega,
        omega_motor,
        k;
    double
        distance,
        v_linear,scale;
    //end definitions
    
    subDecouple=subm(jointDecoupling, range(0,1), range(0,1));
    
    while(1){

    
        J = subm(getJ(), range(0,1), range(0,1));
        O_effector = rowm(getO(LAST_FRAME), range(0,1));

        // Direction determination
        direction = Target - O_effector;
        distance = length(direction);
        v_linear=min(v_Max,sqrt(2*a_Max*distance));
        v_vect = v_linear*normalize(direction);

        // Angular ration calculation
        omega = inv(J)*v_vect;
        omega_motor= subDecouple*omega;
        
        
        scale = min(1.0,min(pointwise_multiply(rowm(omega_Max,range(0,1)),reciprocal(omega_motor))));
    
        omega_motor = omega_motor*scale;
       

        if(distance < 0.01){
            break;
        }

    }
    
}

void Arm::updateDH(){
    matrix<double,NUMBER_OF_JOINTS,1> q_motor,q_joint;
    for(int joint = 1;joint<=NUMBER_OF_JOINTS;joint++){
        //motorAdresses(joint-1)->getPos();  
        //q_motor(joint-1) = motorAdresses(joint-1)->presentPos();
    }
    q_joint=motorDecoupling*q_motor;
    for(int joint = 1;joint<=NUMBER_OF_JOINTS;joint++){
        if (E(joint-1)==ROTATIONAL) { 
            //rotational
            DH(joint-1,3)=q_joint(joint-1);            
        }else /* == PRISMATIC */{
            //prismatic 
            DH(joint-1,2)=q_joint(joint-1);
        }
    }
    TransformExists = false;  
}


matrix<double,4,4> Arm::getT(int a, int b){
    matrix<double,4,4> Rt, Pd, Ra, Pa, T;
    
    if (TransformExists(a,b)) return Transform(a,b);
    
    if (a==b)                                   T = identity_matrix<double>(4);
    else if (a==0 && b==1)                      T = T01;
    else if (a==LAST_FRAME-1 && b==LAST_FRAME)  T = Tend;
    else if (a==1 && b==0)                      T = inv(T01);
    else if (a==LAST_FRAME && b==LAST_FRAME-1)  T = inv(Tend);
    else if (a<b){
        if (a+1==b){
            // translation along X    
            Pa= 1, 0, 0, DH(a,1),
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1;
             // rotation of Z about X
            Ra= 1, 0, 0, 0, 
                0, cos(DH(a,0)), -sin(DH(a,0)), 0,
                0, sin(DH(a,0)),  cos(DH(a,0)), 0,
                0, 0, 0, 1;
            // translation along Z
            Pd= 1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, DH(a,2),
                0, 0, 0, 1;
            //rotation of X about Z
            Rt= cos(DH(a,3)), -sin(DH(a,3)), 0, 0,
                sin(DH(a,3)),  cos(DH(a,3)), 0, 0,
                0 ,0 ,1 ,0,
                0, 0, 0, 1;
            T=Pa*Ra*Pd*Rt;
        } else {
            T=getT(a,b-1)*getT(b-1,b);
        }      
    } else {
        if (a==b+1) {
            //rotation of X about Z
            Rt= cos(-DH(b,3)), -sin(-DH(b,3)), 0, 0,
                sin(-DH(b,3)),  cos(-DH(b,3)), 0, 0,
                0 ,0 ,1 ,0,
                0, 0, 0, 1;
            // translation along Z
            Pd= 1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, -DH(b,2),
                0, 0, 0, 1;
            // rotation of Z about X
            Ra= 1, 0, 0, 0, 
                0, cos(-DH(b,0)), -sin(-DH(b,0)), 0,
                0, sin(-DH(b,0)),  cos(-DH(b,0)), 0,
                0, 0, 0, 1;
            // translation along X    
            Pa= 1, 0, 0, -DH(b,1),
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1;
            T=Rt*Pd*Ra*Pa;
        } else {
            T=getT(a,b+1)*getT(b+1,b);
        }   
    }
    Transform(a,b)=T;
    TransformExists(a,b)=true;
    return T;      
}

matrix<double,3,1> Arm::getO(int frame){
    // returns x,y,z coordinates of the origin of the selected frame
    return subm(getT(0,frame),range(0,2),range(3,3));    
}


matrix<double,3,NUMBER_OF_FRAMES> Arm::getO(){
    // returns x,y,z coordinates of the origines 
    matrix<double,3,NUMBER_OF_FRAMES> O;
    for (int frame = 0 ; frame<NUMBER_OF_FRAMES ;frame++){  
        set_colm(O,frame)=getO(frame);
    }
    return O;    
}

matrix<double,3,1> Arm::getZ(int frame){
    return subm(getT(0,frame),range(0,2),range(2,2));    
}


matrix<double,3,NUMBER_OF_FRAMES> Arm::getZ(){
    // returns coordinates of the z axes of each frame expressed in frame 0
    matrix<double,3,NUMBER_OF_FRAMES> Z;
    for (int frame = 0 ; frame<NUMBER_OF_FRAMES ;frame++){  
        set_colm(Z,frame)=getZ(frame);
    }
    return Z;        
}

matrix<double,6,NUMBER_OF_JOINTS> Arm::getJ(){
    // returns the jacobian of the arm
    matrix<double,6,NUMBER_OF_JOINTS> J;
    
    // J is calculated    
    for (int joint = 1;joint<=NUMBER_OF_JOINTS;joint++){
        if (E(joint-1)==ROTATIONAL) {
            // rotational joint
            set_subm(J,range(0,2),range(joint-1,joint-1))=cross(getZ(joint),getO(LAST_FRAME)-getO(joint));
            set_subm(J,range(3,5),range(joint-1,joint-1))=getZ(joint);
            
        }
        else{
            //prismatic joint
            set_subm(J,range(0,2),range(joint-1,joint-1))=getZ(joint);
            set_subm(J,range(3,5),range(joint-1,joint-1))=0;
        }    
    }  
    return J; 
}
 
 
matrix<double,3,1> Arm::cross(matrix<double,3,1> A , matrix<double,3,1> B){
    // cross product implementation for matrixes
    matrix<double,3,1> C;
    C(0)=A(1)*B(2)-A(2)*B(1);
    C(1)=A(2)*B(0)-A(0)*B(2);
    C(2)=A(0)*B(1)-A(1)*B(0);
    return C;  
}


matrix<double,3,NUMBER_OF_JOINTS> Arm::getJv(int i){ 
    // 1  <= i <= NUMBER_OF_JOINTS 
    matrix<double,3,NUMBER_OF_JOINTS> Jv;
    Jv=0;
    
    for(int joint=1;joint<=i;joint++){
        if (E(joint-1)==ROTATIONAL){
            //rotational
            set_colm(Jv,joint-1)=cross(getZ(joint),getC(i)-getO(joint));
        } else {
            //prismatic
            set_colm(Jv,joint-1)=getZ(joint);
            
        }
    }
    return Jv;
}
        


matrix<double,3,NUMBER_OF_JOINTS> Arm::getJw(int i){
    // 1  <= i <= NUMBER_OF_JOINTS 
    matrix<double,3,NUMBER_OF_JOINTS> Jw;
    Jw =0;
    for(int joint=1;joint<=i;joint++){
        if (E(joint-1)==ROTATIONAL){
            //rotational
            set_colm(Jw,joint-1)=getZ(joint);
            
        } else {
            //prismatic
            //do nothing column is allready 0
        }
    }
    return Jw;
}      

matrix<double,3,1> Arm::getC(int joint){
    //gets Cjoint expressed in frame 0
    matrix<double,4,1> Ci;
    Ci=1;
    set_rowm(Ci,range(0,2))=colm(C,joint-1);
    return rowm(getT(0,joint)*Ci,range(0,2));
}


matrix<double,3,3> Arm::getIc(int joint){
    //gets Ic of the joint expressed in frame 0
    
    
    matrix<double,4,4> Ii;
    Ii=identity_matrix<double>(4); 
    set_subm(Ii,range(0,2),range(0,2))=Ic(joint-1);
    return subm(getT(0,joint)*Ii*getT(joint,0),range(0,2),range(0,2));
    
    
}


matrix<double,NUMBER_OF_JOINTS,NUMBER_OF_JOINTS> Arm::getM(){
    
    matrix<double,NUMBER_OF_JOINTS,NUMBER_OF_JOINTS> M;
    matrix<double,3,NUMBER_OF_JOINTS> Jv,Jw;
    M = 0;
    
    for (int joint = 1; joint <= NUMBER_OF_JOINTS;joint++){
        Jv = getJv(joint);
        Jw = getJw(joint);
        M=M+mc(joint-1)*trans(Jv)*Jv + trans(Jw)*getIc(joint)*Jw;
    }
    return M;    
}
    




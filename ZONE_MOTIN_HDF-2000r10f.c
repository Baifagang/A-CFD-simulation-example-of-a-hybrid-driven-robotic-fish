#include"udf.h"

#define rot_speed 2000
#define dis_x 0.1253876
#define dis_y 0.2305749
#define dis_z 0.014386
#define fr 1
#define Pai 3.141592654

#define head_ID 270
#define joint1_ID 286
#define joint2_ID 290
#define tail_ID 294
#define left_prop_ID 205
#define right_prop_ID 201

#define t_act 0.2
#define mass 48.099

#define l1 0.351
#define l2 0.221
#define l3 0.244

static real dx = 0.0;
static real x = 0.0;
static real ve = 0.0;
static real a = 0.0;
static real y = 0.0;


DEFINE_EXECUTE_AT_END(execute_at_end)
{
	
	real dtime = CURRENT_TIMESTEP;
	real time = CURRENT_TIME;
	
	real f_glob[3] = { 0,0,0 };
	real m_glob[3] = { 0,0,0 };

	real f_head[3] = { 0,0,0 };
	real m_head[3] = { 0,0,0 };
	
    real f_joint1[3] = { 0,0,0 };
	real m_joint1[3] = { 0,0,0 };
	
	real f_joint2[3] = { 0,0,0 };
	real m_joint2[3] = { 0,0,0 };

	real f_tail[3] = { 0,0,0 };
	real m_tail[3] = { 0,0,0 };
	
	real f_left_prop[3] = { 0,0,0 };
	real m_left_prop[3] = { 0,0,0 };

	real f_right_prop[3] = { 0,0,0 };
	real m_right_prop[3] = { 0,0,0 };
	
	real head_cg[3] = { 0,0,0 };
	real joint1_cg[3] = { 0,0,0 };
	real joint2_cg[3] = { 0,0,0 };
	real tail_cg[3] = { 0,0,0 };
	real left_prop_cg[3] = { 0,0,0 };
	real right_prop_cg[3] = { 0,0,0 };
	
	head_cg[0] = x;
	head_cg[1] = y;
	
	joint1_cg[0] = x-l1;
	joint1_cg[1] = y+l1*5.757*sin(5.76*fr*time+1.635)*(1-exp(-time/t_act))*(Pai/180);
	
	joint2_cg[0] = x-(l1+l2);
	joint2_cg[1] = y+(l1*5.757*sin(5.76*fr*time+1.635)+l2*4.353*sin(5.76*fr*time+2.031))*(1-exp(-time/t_act))*(Pai/180);
	
	tail_cg[0] = x-(l1+l2+l3);
	tail_cg[1] = y+(l1*5.757*sin(5.76*fr*time+1.635)+l2*4.353*sin(5.76*fr*time+2.031)+l3*7.679*sin(5.76*fr*time+2.356))*(1-exp(-time/t_act))*(Pai/180);
	
	left_prop_cg[0] = x-dis_x;
	left_prop_cg[1] = y+dis_y;
	
	right_prop_cg[0] = x-dis_x;
	right_prop_cg[1] = y-dis_y;
	
	


#if RP_NODE

	if (!Data_Valid_P())
		return;
	Domain *domain1 = Get_Domain(1);  //return fluid domain
	Thread *tf1 = Lookup_Thread(domain1, head_ID);
	Compute_Force_And_Moment(domain1, tf1, head_cg, f_head, m_head, FALSE);
	
	Domain *domain2 = Get_Domain(1);  //return fluid domain
	Thread *tf2 = Lookup_Thread(domain2, joint1_ID);
	Compute_Force_And_Moment(domain2, tf2, joint1_cg, f_joint1, m_joint1, FALSE);
	
	Domain *domain3 = Get_Domain(1);  //return fluid domain
	Thread *tf3 = Lookup_Thread(domain3, joint2_ID);
	Compute_Force_And_Moment(domain3, tf3, joint2_cg, f_joint2, m_joint2, FALSE);
	
	Domain *domain4 = Get_Domain(1);  //return fluid domain
	Thread *tf4 = Lookup_Thread(domain4, tail_ID);
	Compute_Force_And_Moment(domain4, tf4, tail_cg, f_tail, m_tail, FALSE);	
	
	Domain *domain5 = Get_Domain(1);  //return fluid domain
	Thread *tf5 = Lookup_Thread(domain5, left_prop_ID);
	Compute_Force_And_Moment(domain5, tf5, left_prop_cg, f_left_prop, m_left_prop, FALSE);
	
	Domain *domain6 = Get_Domain(1);  //return fluid domain
	Thread *tf6 = Lookup_Thread(domain6, right_prop_ID);
	Compute_Force_And_Moment(domain6, tf6, right_prop_cg, f_right_prop, m_right_prop, FALSE);	
					

	real ve_before = ve, a_before = a;
	
	f_glob[0] = f_head[0]+f_joint1[0]+f_joint2[0]+f_tail[0]+f_left_prop[0]+f_right_prop[0];
	
	a = f_glob[0]/mass;
	ve = ve + (a + a_before)*dtime / 2;
	dx = (ve + ve_before)*dtime / 2;
	x = x + dx;


#endif
	node_to_host_real_4(a, ve, x,dx);

#if RP_HOST

	FILE *fpx = NULL;
	fpx = fopen("positionx.txt", "a");
	fprintf(fpx, "%.32f %.32f %.32f %.32f %.32f  \n", time, x, ve, a,dx);
	fclose(fpx);


#endif
}

DEFINE_ZONE_MOTION(moving, omega, axis, origin, velocity, time, dtime)
{	
	*omega = 0.0;
	origin[0] = x;
	origin[1] = 0.0;
	N3V_D(axis,=,0.0,0.0,1.0);
	velocity[0] = ve;
    velocity[1] = 0.0;
return;
}

DEFINE_ZONE_MOTION(leftprop_rot,omega,axis,origin,vel,time,dtime)
{
	real ox1,oy1,oz1;
	vel[0]= ve;
	vel[1]= 0;
	vel[2]= 0;
	ox1= -dis_x+vel[0]*(time-dtime);
	oy1= dis_y+vel[1]*(time-dtime);
	oz1= -dis_z+vel[2]*(time-dtime);
	origin[0]=ox1;
	origin[1]=oy1;
	origin[2]=oz1;
	axis[0]=1;
	axis[1]=0;
	axis[2]=0;
	*omega=rot_speed*Pai/30;
	return;
}


DEFINE_ZONE_MOTION(rightprop_rot,omega,axis,origin,vel,time,dtime)
{
	real ox1,oy1,oz1;
	vel[0]= ve;
	vel[1]= 0;
	vel[2]= 0;
	ox1= -dis_x+vel[0]*(time-dtime);
	oy1= -dis_y+vel[1]*(time-dtime);
	oz1= -dis_z+vel[2]*(time-dtime);
	origin[0]=ox1;
	origin[1]=oy1;
	origin[2]=oz1;
	axis[0]=1;
	axis[1]=0;
	axis[2]=0;
	*omega=-rot_speed*Pai/30;
	return;
}


DEFINE_CG_MOTION(Joint1,dt,vel,omega,time,dtime)
{
	NV_S(vel, =, 0.0);
	NV_S(omega, =, 0.0);
	vel[0]= ve;
	omega[2] =(1-exp(-time/t_act))*5.76*fr*5.757*cos(5.76*fr*time+1.635)*(Pai/180);
}
DEFINE_CG_MOTION(Joint2,dt,vel,omega,time,dtime)
{
	NV_S(vel, =, 0.0);
	NV_S(omega, =, 0.0);
	vel[0]= 0.0;
	omega[2] =(1-exp(-time/t_act))*5.76*fr*4.353*cos(5.76*fr*time+2.031)*(Pai/180);
}
DEFINE_CG_MOTION(Tail,dt,vel,omega,time,dtime)
{
	NV_S(vel, =, 0.0);
	NV_S(omega, =, 0.0);
	vel[0]= 0.0;
	omega[2] =(1-exp(-time/t_act))*5.76*fr*7.679*cos(5.76*fr*time+2.356)*(Pai/180);
}


DEFINE_ON_DEMAND(assignment)
{
	real x_date[4] = { 6.63582218871586704267429013270885, 1.29842501077505279027946016867645, -0.09709791289928049184965175300022, 0.01298959443089301932228174507600 };
	x = x_date[0];
	ve = x_date[1];
	a = x_date[2];
	dx = x_date[3];
	Message("x=%f,vx=%f,a=%f,dx=%f \n",x,ve,a,dx);

} 


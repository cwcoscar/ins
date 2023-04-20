#include <string>
#include <iostream>
#include <math.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#define PI 3.141592652
typedef struct
{
   double Lat;
   double Lon;
   double High;
}llh_InitTypeDef;
typedef struct
{
   double xx;
   double yy;
   double zz;
}xyz_InitTypeDef;
typedef struct
{
   double e;
   double n;
   double u;
}enu_InitTypeDef;
typedef struct
{
   double n;
   double e;
   double d;
}ned_InitTypeDef;


xyz_InitTypeDef Radii_of_curvature(double Lat){

   xyz_InitTypeDef Rmn;
   double R_0 = 6378137;
   double e = 0.0818191908425;
   double temp = 1 - pow((e * sin(Lat*PI/180)),2); 
   Rmn.xx = R_0 * (1 - e*e) / pow(temp,1.5);//RM
   Rmn.yy = R_0 / sqrt(temp);//Rn
   return Rmn;
}

llh_InitTypeDef lever_arm(double L1, double L2, double H1,double roll,double pitch,double yaw,double lever_arm_x,double lever_arm_y,double lever_arm_z)
{

   llh_InitTypeDef LLA;
   static xyz_InitTypeDef radii = Radii_of_curvature(L1);

   roll = -roll;
   pitch = - pitch;
   yaw = -yaw;

   Eigen::Vector3f b2a_in_base(lever_arm_x, lever_arm_y, lever_arm_z); //base_link(IMU) to antenna vector in vehicle-frame(base_link)
   Eigen::Matrix3f tf_baseton;
   Eigen::Matrix3f D_inv;
   D_inv << 1/(radii.xx + H1),                               0.0,   0.0,
                          0.0,  1/(radii.yy + H1)/cos(L1*PI/180),   0.0,
                          0.0,                               0.0,  -1.0;
   Eigen::Vector3f b2a_in_n; //base_link(IMU) to antenna vector in navigation-frame

   Eigen::AngleAxisf rot_x((roll + 180)*PI/180, Eigen::Vector3f::UnitX());    
   Eigen::AngleAxisf rot_y(pitch*PI/180, Eigen::Vector3f::UnitY());
   Eigen::AngleAxisf rot_z(yaw*PI/180, Eigen::Vector3f::UnitZ());


   tf_baseton = (rot_x*rot_y*rot_z).matrix(); //the angle between vehicle-frame(base_link) and map-frame along z-axis
   b2a_in_n = D_inv * tf_baseton * b2a_in_base;

   LLA.Lat = L1 - b2a_in_n(0,0) * 180/PI;
   LLA.Lon = L2 - b2a_in_n(1,0) * 180/PI;
   LLA.High = H1 - b2a_in_n(2,0); 



   return (LLA);
}


xyz_InitTypeDef lla2xyz(double L1, double L2, double H1)
{
   double a, b, e;
   double sinphi, cosphi, coslam, sinlam, tan2phi;
   double tmp, tmpden, tmp2;
   xyz_InitTypeDef xyzC;

   L1=(L1*PI)/180;
   L2=(L2*PI)/180;

   a = 6378137.0000;
   b = 6356752.3142;
   e = sqrt(1-(b/a)*(b/a));  

   sinphi = sin(L1);
   cosphi = cos(L1);
   coslam = cos(L2);
   sinlam = sin(L2);
   tan2phi = (tan(L1))*(tan(L1));
   tmp = 1 - e*e;
   tmpden = sqrt( 1 + tmp*tan2phi );

   xyzC.xx = (a*coslam)/tmpden + H1*coslam*cosphi;

   xyzC.yy = (a*sinlam)/tmpden + H1*sinlam*cosphi;

   tmp2 = sqrt(1 - e*e*sinphi*sinphi);
   xyzC.zz = (a*tmp*sinphi)/tmp2 + H1*sinphi;

   return (xyzC);
}


enu_InitTypeDef lla2enu(double lat0,double lon0,double alt0,double lat1,double lon1,double alt1)
{
   xyz_InitTypeDef xe2 = lla2xyz(lat0,lon0,alt0);
   xyz_InitTypeDef xe1 = lla2xyz(lat1,lon1,alt1);

   enu_InitTypeDef enuC;
   double a, b, c;
   double phi, lam, sinphi, cosphi, sinlam, coslam;

   a=xe1.xx-xe2.xx;
   b=xe1.yy-xe2.yy;
   c=xe1.zz-xe2.zz;

   phi=(lat0*PI)/180;
   lam=(lon0*PI)/180;
   sinphi=sin(phi);
   cosphi=cos(phi);
   sinlam=sin(lam);
   coslam=cos(lam);

   enuC.e=(-sinlam)*a+(coslam)*b+(0)*c;
   enuC.n=(-sinphi*coslam)*a+(-sinphi*sinlam)*b+(cosphi)*c;
   enuC.u=(cosphi*coslam)*a+(cosphi*sinlam)*b+(sinphi)*c;


   return (enuC);
}


ned_InitTypeDef lla2ned(double lat0,double lon0,double alt0,double lat1,double lon1,double alt1)
{

   xyz_InitTypeDef xe2 = lla2xyz(lat0,lon0,alt0);
   xyz_InitTypeDef xe1 = lla2xyz(lat1,lon1,alt1);  
   ned_InitTypeDef nedC;
   double a, b, c;
   double phi, lam, sinphi, cosphi, sinlam, coslam;

   a=xe1.xx-xe2.xx;
   b=xe1.yy-xe2.yy;
   c=xe1.zz-xe2.zz;

    phi=(lat0*PI)/180;
   lam=(lon0*PI)/180;
   sinphi=sin(phi);
   cosphi=cos(phi);
   sinlam=sin(lam);
   coslam=cos(lam);

   nedC.e=(-sinlam)*a+(coslam)*b+(0)*c;
   nedC.n=(-sinphi*coslam)*a+(-sinphi*sinlam)*b+(cosphi)*c;
   nedC.d=-((cosphi*coslam)*a+(cosphi*sinlam)*b+(sinphi)*c);

   return (nedC);
}
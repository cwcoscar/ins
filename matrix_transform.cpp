#include "mec_transformation.h"

#include <Eigen/Dense>
#include <iostream>
#include <cmath>

// #include <sensors_msgs/Imu.h>

#define EARTH_SEMIMAJOR 6378137 // meter
#define EARTH_ROTATION_RATE 0.00007292115 // rad/s
#define FLATTENING 1/298.25722563
#define ECCENTRICITY sqrt(FLATTENING*(2-FLATTENING))
#define NCKUEE_LATITUDE 0.401367412 // 22.99665875 degree
#define NCKUEE_LONGITUDE 2.09827994 // 120.222584889 degree
#define NCKUEE_HEIGHT 98.211 // meter

typedef struct ins_mechanization_parameter{
    Eigen::VectorXd f_b; // measurement: accelerometer
    Eigen::MatrixXd omega_ib_b; // measurement: gyroscope
    Eigen::MatrixXd omega_il_b;
    Eigen::MatrixXd omega_ie_b;
    Eigen::MatrixXd omega_el_b;
    Eigen::MatrixXd D_inverse
    Eigen::VectorXd r_dot_l;
    Eigen::VectorXd v_dot_l;
    Eigen::MatrixXd R_dot_b_l;
    double g_l;
}mechanization;

typedef struct state_vector{
    Eigen::VectorXd r_l; // lat long h
    Eigen::VectorXd v_l; // enu
    Eigen::MatrixXd R_b_l; 

}state;

Eigen::MatrixXd Skew_symmetric_transform(const Eigen::VectorXd& vector, int dim){
    int length = vector.size();
    Eigen::MatrixXd result(length,length);
    if (dim == 3){
        result << 0,          -vector(2), vector(1),
                  vector(2),  0,          -vector(0),
                  -vector(1), vector(0),  0;
    }
    else if (dim == 4){
        result << 0,          vector(2),  -vector(1), vector(0),
                  -vector(2), 0,          vector(0),  vector(1),
                  vector(1),  -vector(0), 0,          vector(2),
                  -vector(0), -vector(1), -vector(2), 0;
    }
    else{
        std::cout << "The dimension of transformation is not supported." << std::endl;
    }
    return result;
}

double gravity(double lat_rad){
    return 9.7803267715*(1+0.0052790414*pow(sin(lat_rad),2)+0.0000232718*pow(sin(lat_rad),4));
}

double earth_radius_along_meridian(double lat_rad){
    double m = (EARTH_SEMIMAJOR*(1-ECCENTRICITY*ECCENTRICITY))/pow(1-pow(ECCENTRICITY*sin(lat_rad),2),1.5);
    return m;
}

double earth_radius_along_prime_vertical(double lat_rad){
    double n = EARTH_SEMIMAJOR/sqrt(1-pow(ECCENTRICITY*sin(lat_rad),2));
    return n;
}

Eigen::VectorXd w_el_l_update(Eigen::VectorXd vel_enu, double lat_rad, double h){
    Eigen::VectorXd result(3);
    result << -vel_enu(1)/(earth_radius_along_meridian(lat_rad)+h), 
              vel_enu(0)/(earth_radius_along_prime_vertical(lat_rad)+h), 
              (vel_enu(0)*tan(lat_rad))/(earth_radius_along_prime_vertical(lat_rad)+h);
    return result;
}

Eigen::VectorXd w_ie_l_update(double lat_rad){
    Eigen::VectorXd result(3);
    result << 0, 
              EARTH_ROTATION_RATE*cos(lat_rad),
              EARTH_ROTATION_RATE*sin(lat_rad);
    return result;
}

Eigen::MatrixXd omega_el_l_update(Eigen::VectorXd vel_enu, double lat_rad, double h){
    return Skew_symmetric_transform(w_el_l_update(vel_enu,lat_rad,h), 3);
}

Eigen::MatrixXd omega_ie_l_update(double lat_rad){
    return Skew_symmetric_transform(w_ie_l_update(lat_rad), 3);
}

Eigen::MatrixXd omega_il_b_update(Eigen::VectorXd vel_enu, double lat_rad, double h, Eigen::MatrixXd R_b_l){
    Eigen::VectorXd tmp(3);
    tmp = w_l_el_update(vel_enu,lat_rad,h) + w_l_ie_update(lat_rad);
    return Skew_symmetric_transform(R_b_l.transpose()*tmp, 3);
}

Eigen::MatrixXd D_inverse_update(double lat_rad, double h){
    Eigen::MatrixXd D_inverse(3, 3);
    D_inverse << 0,                                                                 1/(earth_radius_along_meridian(vec.r_l(0))+vec.r_l(2)), 0,
                 1/((earth_radius_along_prime_vertical(vec.r_l(0))+h)*cos(vec(0))), 0,                                                      0,
                 0,                                                                 0,                                                      1;
    return D_inverse;
}


// Think about how to update variables with vec as input !!!!!!



void R_dot_l_b_update(mechanization& variables, state& vec){
    variables.omega_il_b = omega_il_b_update(vec.v_l, vec.r_l(0), vec.r_l(2), vec.R_b_l);
    variables.R_dot_b_l = vec.R_b_l * (omega_ib_b - omega_il_b);
}

void v_dot_l(mechanization& variables, state& vec){
    variables.omega_ie_l = omega_ie_l_update(vec.r_l(0));
    variables.omega_el_l = omega_el_l_update(vec.v_l, vec.r_l(0), vec.r_l(2));
    variables.g_l = gravity(vec.r_l(0));
    variables.v_dot_l = vec.R_b_l * variables.f_b - (2*variables.omega_ie_l + variables.omega_el_l)*vec.v_l + variables.g_l;
}

void r_dot_l(mechanization& variables, state& vec){
    variables.D_inverse = D_inverse_update(vec.r_l(0), vec.r_l(2))
    variables.r_dot_l = variables.D_inverse * vec.v_l;
}

int main(int argc, char **argv) {
    Eigen::VectorXd v(3);
    v << 10, 10, 0;
    Eigen::MatrixXd result = Skew_symmetric_transform(v,3);

    Eigen::Matrix<double, 1, 4> a = {1, 2, 3, 4};

    std::cout << omega_il_b_update(v, NCKUEE_LATITUDE, NCKUEE_HEIGHT) << std::endl;
    return 0;
}
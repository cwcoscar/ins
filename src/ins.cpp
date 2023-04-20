#include "ins.h"
#include "coordinate_mat_transformation.h"

using namespace INS;

static bool fix_flag = false;
static bool vel_flag = false;
static bool att_flag = false;

Eigen::Vector3d a_b(-0.348513345960806, -0.26021251227717, 0.132337782341719);
Eigen::Vector3d a_s(-0.00426378745053201, 0.000725755116990245, 0.00263712381843959);
Eigen::Vector3d g_b(0.00710659149934062, 0.00211909908717263, -0.0000592951099686292);
Eigen::Vector3d g_s(-2.36707629594559, -0.490347919324706, -0.686283178454847);
Eigen::Vector3d nckuee(NCKUEE_LATITUDE, NCKUEE_LONGITUDE, NCKUEE_HEIGHT);

Ins_mechanization::Ins_mechanization(ros::Publisher pub_fix, ros::Publisher pub_vel, ros::Publisher pub_att) 
: pub_ins_fix(pub_fix), pub_ins_vel(pub_vel), pub_ins_att(pub_att){
    imu_correction.acc_bias = a_b;
    imu_correction.acc_scale = a_s;
    imu_correction.gyro_bias = g_b;
    imu_correction.gyro_scale = g_s;
}

inline double Ins_mechanization::gravity(double lat_rad){
    return 9.7803267715*(1+0.0052790414*pow(sin(lat_rad),2)+0.0000232718*pow(sin(lat_rad),4));
}

inline double Ins_mechanization::earth_radius_along_meridian(double lat_rad){
    double m = (EARTH_SEMIMAJOR*(1-ECCENTRICITY*ECCENTRICITY))/pow(1-pow(ECCENTRICITY*sin(lat_rad),2),1.5);
    return m;
}

inline double Ins_mechanization::earth_radius_along_prime_vertical(double lat_rad){
    double n = EARTH_SEMIMAJOR/sqrt(1-pow(ECCENTRICITY*sin(lat_rad),2));
    return n;
}

Eigen::VectorXd Ins_mechanization::w_el_l_update(Eigen::VectorXd vel_enu, double lat_rad, double h){
    Eigen::VectorXd result(3);
    result << -vel_enu(1)/(earth_radius_along_meridian(lat_rad)+h), 
            vel_enu(0)/(earth_radius_along_prime_vertical(lat_rad)+h), 
            (vel_enu(0)*tan(lat_rad))/(earth_radius_along_prime_vertical(lat_rad)+h);
    return result;
}

Eigen::VectorXd Ins_mechanization::w_ie_l_update(double lat_rad){
    Eigen::VectorXd result(3);
    result << 0, 
            EARTH_ROTATION_RATE*cos(lat_rad),
            EARTH_ROTATION_RATE*sin(lat_rad);
    return result;
}

Eigen::MatrixXd Ins_mechanization::omega_el_l_update(state& vector){
    return coordinate_mat_transformation::Skew_symmetric_transform(w_el_l_update(vector.v_l, vector.r_l(0), vector.r_l(2)), 3);
}

Eigen::MatrixXd Ins_mechanization::omega_ie_l_update(state& vector){
    return coordinate_mat_transformation::Skew_symmetric_transform(w_ie_l_update(vector.r_l(0)), 3);
}

Eigen::MatrixXd Ins_mechanization::omega_il_b_update(state& vector, int dim){ //dim = 3 or 4
    Eigen::VectorXd w_il_b(3);
    w_il_b = w_el_l_update(vector.v_l,vector.r_l(0),vector.r_l(2)) + w_ie_l_update(vector.r_l(0));
    return coordinate_mat_transformation::Skew_symmetric_transform(vector.R_b_l.transpose()*w_il_b, dim);
}

Eigen::MatrixXd Ins_mechanization::D_inverse_update(state& vector){
    Eigen::MatrixXd D_inverse(3, 3);
    D_inverse << 0,                                   1/(earth_radius_along_meridian(vector.r_l(0))+vector.r_l(2)),       0,
                1/((earth_radius_along_prime_vertical(vector.r_l(0))+vector.r_l(2))*cos(vector.r_l(0)*DEG_TO_RAD)), 0,    0,
                0,                                    0,                                                                  1;
    return D_inverse;
}

void Ins_mechanization::R_dot_l_b_update(mechanization& variables, state& vec){
    Eigen::MatrixXd omega_il_b;
    omega_il_b = omega_il_b_update(vec,3);

    Eigen::MatrixXd omega_ib_b;
    omega_ib_b = coordinate_mat_transformation::Skew_symmetric_transform(variables.w_ib_b, 3);
    variables.R_dot_b_l_past = variables.R_dot_b_l_now;
    variables.R_dot_b_l_now = vec.R_b_l * (omega_ib_b - omega_il_b);
}

void Ins_mechanization::Q_dot_update(mechanization& variables, state& vec){
    Eigen::MatrixXd omega_il_b;
    omega_il_b = omega_il_b_update(vec,4);
    // std::cout << "w_il_b" << std::endl <<  w_el_l_update(vec.v_l,vec.r_l(0),vec.r_l(2)) + w_ie_l_update(vec.r_l(0)) << std::endl;
    Eigen::MatrixXd omega_ib_b;
    omega_ib_b = coordinate_mat_transformation::Skew_symmetric_transform(variables.w_ib_b, 4);
    // std::cout << "w_ib_b" << std::endl << variables.w_ib_b << std::endl;
    Eigen::VectorXd Q;
    Q = coordinate_mat_transformation::att_Q_transform(vec.att_l);
    variables.Q_dot_past = variables.Q_dot_now;
    variables.Q_dot_now = 0.5 * (omega_ib_b - omega_il_b) * Q;
}

void Ins_mechanization::v_dot_l_update(mechanization& variables, state& vec){
    Eigen::MatrixXd omega_ie_l;
    omega_ie_l = omega_ie_l_update(vec);
    Eigen::MatrixXd omega_el_l;
    omega_el_l = omega_el_l_update(vec);
    
    Eigen::Vector3d g_l(0, 0, -gravity(vec.r_l(0)));
    variables.v_dot_l_past = variables.v_dot_l_now;
    variables.v_dot_l_now = vec.R_b_l * variables.f_b - (2*omega_ie_l + omega_el_l)*vec.v_l + vec.R_b_l * g_l;

    // std::cout << "att" << std::endl << vec.att_l << std::endl;
    // std::cout << "vec.R_b_l" << std::endl << vec.R_b_l << std::endl;
    // std::cout << "variables.f_b" << std::endl << variables.f_b << std::endl;
    // std::cout << "vec.R_b_l * variables.f_b" << std::endl << vec.R_b_l * variables.f_b << std::endl;
}

void Ins_mechanization::r_dot_l_update(mechanization& variables, state& vec){
    Eigen::MatrixXd D_inverse;
    D_inverse = D_inverse_update(vec); 
    variables.r_dot_l_past = variables.r_dot_l_now;
    variables.r_dot_l_now = D_inverse * vec.v_l;
    Eigen::Vector3d r_dot_l_enu = coordinate_mat_transformation::lla2enu(variables.r_dot_l_now, nckuee);
}

void Ins_mechanization::GNSSfixcallback(const sensor_msgs::NavSatFix& msg){
    Eigen::Vector3d pos(msg.latitude, msg.longitude, msg.altitude);
    state_vector.r_l = pos;
    fix_flag = true;
}

void Ins_mechanization::GNSSvelcallback(const geometry_msgs::TwistWithCovarianceStamped& msg){
    if (att_flag == true){
        Eigen::Vector3d vel_l(msg.twist.twist.linear.x, msg.twist.twist.linear.y, msg.twist.twist.linear.z);
        state_vector.v_l = vel_l;
        vel_flag = true;
    }
    else{
        std::cout << "\033[33m" << "Waiting for initial attitude !" << "\033[0m" << std::endl;
    }
}
            
void Ins_mechanization::GNSSattcallback(const ublox_msgs::NavATT& msg){ 
    Eigen::Vector3d att(msg.roll/100000, msg.pitch/100000, msg.heading/100000);
    // NED -> ENU
    att(2) = 360 - att(2);
    att(1) = att(0);
    att(0) = att(1);

    // 0 ~ 2PI -> -PI ~ PI
    for (int i = 0; i < 3; i++){
        if (att(i) > 180 && att(i) < 360){
            att(i) = att(i) - 360;
        }
    }
    att = att*DEG_TO_RAD; 
    // std::cout << "att" << std::endl << att << std::endl;
    state_vector.R_b_l = coordinate_mat_transformation::Rotation_matrix(att);
    state_vector.att_l = att;
    att_flag = true;
}

void Ins_mechanization::Imucallback(const sensor_msgs::Imu& msg){
    // if body frame is enu
    double current_time_tag = msg.header.stamp.toSec();
    Eigen::Vector3d acc_raw(msg.linear_acceleration.x, msg.linear_acceleration.y, msg.linear_acceleration.z);
    Eigen::Vector3d gyro_raw(msg.angular_velocity.x, msg.angular_velocity.y, msg.angular_velocity.z);
    Imu_data_calibration(acc_raw, gyro_raw);
    // acc first then gyro
    if ((current_time_tag - mech_variables.compute_time_tag) > 0.005){
        compute();
        mech_variables.compute_time_tag = current_time_tag;
    }
}

void Ins_mechanization::Imu_data_calibration(Eigen::Vector3d acc_raw, Eigen::Vector3d gyro_raw){
    // save last raw data
    Eigen::Vector3d acc_calibrated;
    Eigen::Matrix3d acc_scale;
    acc_scale << imu_correction.acc_scale(0),                          0,                           0,
                                           0,imu_correction.acc_scale(1),                           0,
                                           0,                          0, imu_correction.acc_scale(2);   
    acc_calibrated = (Eigen::MatrixXd::Identity(3,3) + acc_scale).inverse() * (acc_raw - imu_correction.acc_bias);
    acc_calibrated = (acc_raw - imu_correction.acc_bias);
    Eigen::Vector3d gyro_calibrated;
    Eigen::Matrix3d gyro_scale;
    gyro_scale << imu_correction.gyro_scale(0),                           0,                            0,
                                             0,imu_correction.gyro_scale(1),                            0,
                                             0,                           0, imu_correction.gyro_scale(2);   
    gyro_calibrated = (Eigen::MatrixXd::Identity(3,3) + gyro_scale).inverse() * (gyro_raw - imu_correction.gyro_bias);
    // if body frame is enu
    // mech_variables.f_b = acc_calibrated;
    // mech_variables.w_ib_b = gyro_calibrated;

    // if body frame is nwu
    mech_variables.f_b(0) = -acc_calibrated(1);
    mech_variables.f_b(1) = acc_calibrated(0);
    mech_variables.f_b(2) = acc_calibrated(2);
    mech_variables.w_ib_b(0) = -gyro_calibrated(1);
    mech_variables.w_ib_b(1) = gyro_calibrated(0);
    mech_variables.w_ib_b(2) = gyro_calibrated(2);
    // std::cout << "f_b" << std::endl << mech_variables.f_b << std::endl;
    // std::cout << "w_ib_b" << std::endl << mech_variables.w_ib_b << std::endl;
}

void Ins_mechanization::Initialize_state(){
    Eigen::Vector3d unit(0,0,0);
    Eigen::Vector4d unit_q(0,0,0,1);

    state_vector.r_l = unit;
    state_vector.v_l = unit;
    state_vector.R_b_l = Eigen::MatrixXd::Identity(3,3);
    state_vector.att_l = unit;
    mech_variables.f_b = unit;
    mech_variables.w_ib_b = unit;
    mech_variables.r_dot_l_now = unit;
    mech_variables.r_dot_l_past = unit;
    mech_variables.v_dot_l_now = unit;
    mech_variables.v_dot_l_past = unit;
    mech_variables.R_dot_b_l_now = Eigen::MatrixXd::Identity(3,3);
    mech_variables.R_dot_b_l_past = Eigen::MatrixXd::Identity(3,3);
    mech_variables.Q_dot_now = unit_q;
    mech_variables.Q_dot_past = unit_q;
}

void Ins_mechanization::Attitude_update(){
    // inference the current attitude
    Q_dot_update(mech_variables, state_vector);
    // std::cout << "att_l" << std::endl << state_vector.att_l << std::endl;
    // std::cout << "q" << std::endl << coordinate_mat_transformation::att_Q_transform(state_vector.att_l) << std::endl;
    // std::cout << "q_dot" << std::endl << mech_variables.Q_dot_past * SAMPLING_TIME << std::endl;
    Eigen::VectorXd Q_now = coordinate_mat_transformation::att_Q_transform(state_vector.att_l) + mech_variables.Q_dot_past * SAMPLING_TIME;
    // std::cout << "Q_now" << std::endl << Q_now << std::endl;

    // std::cout << "Q_dot_past" << std::endl << mech_variables.Q_dot_past * SAMPLING_TIME << std::endl;
    // std::cout << "Q_now" << std::endl << Q_now << std::endl;

    // update current attitude and rotation matrix
    Eigen::VectorXd att_now = coordinate_mat_transformation::Q_att_transform(Q_now);
    state_vector.att_l = att_now;
    state_vector.R_b_l = coordinate_mat_transformation::Rotation_matrix(att_now);
    // std::cout << "att_now" << std::endl << att_now << std::endl;
}

void Ins_mechanization::Velocity_update(){
    // inference the current velocity
    v_dot_l_update(mech_variables, state_vector);
    Eigen::VectorXd vel_now = state_vector.v_l + 0.5 * (mech_variables.v_dot_l_now + mech_variables.v_dot_l_past) * SAMPLING_TIME;

    // update current velocity
    state_vector.v_l = vel_now;
    // std::cout << "vel_now" << std::endl << vel_now << std::endl;
}

void Ins_mechanization::Position_update(){
    // inference the current position
    Eigen::VectorXd r_now;
    r_dot_l_update(mech_variables, state_vector);
    r_now = state_vector.r_l + 0.5 * ((mech_variables.r_dot_l_now + mech_variables.r_dot_l_past)*RAD_TO_DEG) * SAMPLING_TIME;

    // update current position
    state_vector.r_l = r_now;
    Eigen::Vector3d pos_enu = coordinate_mat_transformation::lla2enu(r_now, nckuee);
    std::cout << std::fixed << std::setprecision(10);
    // std::cout << "ENU now: " << std::endl << pos_enu << std::endl;
}

void Ins_mechanization::Publish_ins(){
    sensor_msgs::NavSatFix llh;
    geometry_msgs::TwistWithCovarianceStamped vel_l;
    ins::InsATT att_l;
    ros::Time now = ros::Time::now();

    llh.header.stamp = now;
    llh.header.frame_id = "local";
    llh.latitude = state_vector.r_l(0);
    llh.longitude = state_vector.r_l(1);
    llh.altitude = state_vector.r_l(2);

    vel_l.header.stamp = now;
    vel_l.header.frame_id = "local";
    vel_l.twist.twist.linear.x = state_vector.v_l(0);
    vel_l.twist.twist.linear.y = state_vector.v_l(1);
    vel_l.twist.twist.linear.z = state_vector.v_l(2);

    att_l.header.stamp = now;
    att_l.header.frame_id = "local";
    att_l.att_e = state_vector.att_l(0)*RAD_TO_DEG;
    att_l.att_n = state_vector.att_l(1)*RAD_TO_DEG;
    att_l.att_u = state_vector.att_l(2)*RAD_TO_DEG;

    pub_ins_fix.publish(llh);
    pub_ins_vel.publish(vel_l);
    pub_ins_att.publish(att_l);
}

void Ins_mechanization::compute(){
    if(fix_flag == false || vel_flag == false || att_flag == false){
        std::cout << "\033[33m" << "Waiting for initialization !" << "\033[0m" << std::endl;
    }
    else{
        Attitude_update();
        Velocity_update();
        Position_update();
        Publish_ins();
    }
}
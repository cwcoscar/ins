#include "ins.h"
#include "coordinate_mat_transformation.h"

using namespace INS;

static bool fix_flag = false;
static bool vel_flag = false;
static bool att_flag = false;

/* IMU Calibration parameters */
Eigen::Vector3d a_b(-0.348513345960806, -0.26021251227717, 0.132337782341719);
Eigen::Vector3d a_s(-0.00426378745053201, 0.000725755116990245, 0.00263712381843959);
Eigen::Vector3d g_b(0.00710659149934062, 0.00211909908717263, -0.0000592951099686292);
Eigen::Vector3d nckuee(NCKUEE_LATITUDE, NCKUEE_LONGITUDE, NCKUEE_HEIGHT);

Ins_mechanization::Ins_mechanization(ros::Publisher pub_fix, config ins_config) 
: pub_ins_fix_(pub_fix), ins_config_(ins_config){
    imu_correction_.acc_bias = a_b;
    imu_correction_.acc_scale = a_s;
    imu_correction_.gyro_bias = g_b;
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
    return coordinate_mat_transformation::Skew_symmetric_transform(w_el_l_update(vector.v_l, vector.r_l(0)*DEG_TO_RAD, vector.r_l(2)), 3);
}

Eigen::MatrixXd Ins_mechanization::omega_ie_l_update(state& vector){
    return coordinate_mat_transformation::Skew_symmetric_transform(w_ie_l_update(vector.r_l(0)*DEG_TO_RAD), 3);
}

/* dim = 3 or 4 */
Eigen::MatrixXd Ins_mechanization::omega_il_b_update(state& vector, int dim){ 
    Eigen::VectorXd w_il_b(3);
    w_il_b = w_el_l_update(vector.v_l,vector.r_l(0)*DEG_TO_RAD,vector.r_l(2)) + w_ie_l_update(vector.r_l(0)*DEG_TO_RAD);
    return coordinate_mat_transformation::Skew_symmetric_transform(vector.R_b_l.transpose()*w_il_b, dim);
}

Eigen::MatrixXd Ins_mechanization::D_inverse_update(state& vector){
    Eigen::MatrixXd D_inverse(3, 3);
    D_inverse << 0,                                   1/(earth_radius_along_meridian(vector.r_l(0)*DEG_TO_RAD)+vector.r_l(2)),       0,
                1/((earth_radius_along_prime_vertical(vector.r_l(0)*DEG_TO_RAD)+vector.r_l(2))*cos(vector.r_l(0)*DEG_TO_RAD)), 0,    0,
                0,                                    0,                                                                             1;
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
    /* omega_il_b and omega_il_b are supposed to be small compared with measurement */
    Eigen::MatrixXd omega_il_b;
    omega_il_b = omega_il_b_update(vec,4);
    Eigen::MatrixXd omega_ib_b;
    omega_ib_b = coordinate_mat_transformation::Skew_symmetric_transform(variables.w_ib_b, 4);
    Eigen::VectorXd Q;
    Q = coordinate_mat_transformation::att_Q_transform(vec.att_l);
    variables.Q_dot_past = variables.Q_dot_now;
    variables.Q_dot_now = 0.5 * (omega_ib_b - omega_il_b) * Q;
}

void Ins_mechanization::v_dot_l_update(mechanization& variables, state& vec){
    static bool flag_static = false;
    double amplify = 1000;
    static std::vector<Eigen::Vector3d> f_b_window;

    Eigen::MatrixXd omega_ie_l;
    omega_ie_l = omega_ie_l_update(vec);
    Eigen::MatrixXd omega_el_l;
    omega_el_l = omega_el_l_update(vec);
    
    Eigen::Vector3d g_l(0, 0, -gravity(vec.r_l(0)*DEG_TO_RAD));

    if(flag_static == true){
        /* If the vehicle is static, velocity is constraint by a scaler */
        state_vector_.v_l = 0.1*state_vector_.v_l;
        /* Resume last v_dot */
        variables.v_dot_l_past = amplify*variables.v_dot_l_now;
    }
    else{
        variables.v_dot_l_past = variables.v_dot_l_now;
    }

    /* Use Sliding window of f_b to check if vehicle is static */
    /* 
    f_b_mean is the mean acceleration in sliding window which has length @ 100
    g_b_mean is the mean of closest 5 measurements
    */
    Eigen::Vector3d f_b_mean(0,0,0);
    Eigen::Vector3d g_b_mean(0,0,0);
    for(int i = 0; i < f_b_window.size(); i++){
        f_b_mean = f_b_mean + f_b_window[i];
        if(f_b_window.size() >= 5 && i < f_b_window.size() && i >= f_b_window.size()-5){
            g_b_mean = g_b_mean + f_b_window[i];
        }
    }
    f_b_mean = f_b_mean/f_b_window.size();
    g_b_mean = g_b_mean/5;

    /* gravity compensation */
    /* 
    - original method: 
        approximate gravity vector by attitude
    - Revised method:
        Due to unsteady of attitude, gravity can't be derived from attitude accurately
        Calculate mean acceleration of a short period (5 measurements) to approximate the gravity.
        (first 5 measurement use original method)
    */
    if(f_b_window.size() >= 5){
        variables.v_dot_l_now = vec.R_b_l*(variables.f_b - g_b_mean) - (2*omega_ie_l + omega_el_l)*vec.v_l;
    }
    else{
        variables.v_dot_l_now = (vec.R_b_l * variables.f_b - (2*omega_ie_l + omega_el_l)*vec.v_l + vec.R_b_l * g_l);
    }

    /* When vehicle is static, jitter of IMU measurements is small. */
    if((variables.f_b-f_b_mean).norm() < 0.1){
        variables.v_dot_l_now = (1/amplify)*variables.v_dot_l_now;
        flag_static = true;
    }
    else flag_static = false;

    /* Sliding window of f_b is constraint @ length 100 */
    if(f_b_window.size() < 100){
        f_b_window.push_back(variables.f_b);
    }
    else if (f_b_window.size() == 100){
        f_b_window.push_back(variables.f_b);
        f_b_window.erase(f_b_window.begin());
    }

    // std::cout << "------" << std::endl;
    // std::cout << "norm of diff" << std::endl << (variables.f_b-f_b_mean).norm() << std::endl;
    // std::cout << "flag_static" << std::endl << flag_static << std::endl;
    // std::cout << "variables.f_b - g_b_mean: " << std::endl << variables.f_b - g_b_mean << std::endl;
    // std::cout << "v_dot_l_now: " << std::endl << variables.v_dot_l_now << std::endl;
    // std::cout << "------" << std::endl;
}

void Ins_mechanization::r_dot_l_update(mechanization& variables, state& vec){
    Eigen::MatrixXd D_inverse;
    D_inverse = D_inverse_update(vec); 
    variables.r_dot_l_past = variables.r_dot_l_now;
    variables.r_dot_l_now = D_inverse * vec.v_l;
    // std::cout << "r_dot_l_past: " << std::endl << variables.r_dot_l_past << std::endl;
    Eigen::Vector3d r_dot_l_enu = coordinate_mat_transformation::lla2enu(variables.r_dot_l_now, nckuee);
}

void Ins_mechanization::GNSSfixcallback(const sensor_msgs::NavSatFix& msg){
    if(ins_config_.fix_type == 3){
        /* Update once  */
        if(fix_flag == false){
            Eigen::Vector3d pos(msg.latitude, msg.longitude, msg.altitude);
            state_vector_.r_l = pos;
        }
    }
    else if(ins_config_.mode == 0 || fix_flag == false){
        Eigen::Vector3d pos(msg.latitude, msg.longitude, msg.altitude);
        state_vector_.r_l = pos;
    }
    else if(ins_config_.mode == 1){
        /* Don't update  */
    }
    if(ins_config_.mode == 0 && att_flag == true && ins_config_.fix_type != 3) send_tf(state_vector_.r_l, state_vector_.att_l, "ublox");
    fix_flag = true;
    ins_config_.transform2baselink = false;
}

// /ublox_f9k/fix_velocity:ENU
void Ins_mechanization::GNSSvelcallback(const geometry_msgs::TwistWithCovarianceStamped& msg){
    if(ins_config_.fix_type == 3){
        //*Update once  */
        if(ins_config_.odometer == 0 && vel_flag == false){
            Eigen::Vector3d vel_l(msg.twist.twist.linear.x, msg.twist.twist.linear.y, msg.twist.twist.linear.z);
            state_vector_.v_l = vel_l;
            vel_flag = true;
            // std::cout << "\033[33m" << "vel_l" << std::endl << vel_l << "\033[0m" << std::endl;
        }
    }
    else if(ins_config_.odometer == 0 && (ins_config_.mode == 0 || vel_flag == false)){
        Eigen::Vector3d vel_l(msg.twist.twist.linear.x, msg.twist.twist.linear.y, msg.twist.twist.linear.z);
        state_vector_.v_l = vel_l;
        vel_flag = true;
    }
    else if(ins_config_.mode == 1){
        /* Don't update */
    }
}

// /ublox_f9k/navpvt
void Ins_mechanization::GNSSspeedcallback(const ublox_msgs::NavPVT& msg){ 
    if(ins_config_.mode == 0 || vel_flag == false){
        state_vector_.v_forward = (double)msg.gSpeed/1000;
        // std::cout << "\033[33m" << "v_forward" << std::endl << state_vector_.v_forward << "\033[0m" << std::endl;
        vel_flag = true;
    }
    else if(ins_config_.mode == 1){
        state_vector_.v_forward = (double)msg.gSpeed/1000;
        vel_flag = true;
    }
}
            
void Ins_mechanization::GNSSattcallback(const ublox_msgs::NavATT& msg){ 
    if(ins_config_.fix_type == 3){
        /* Update once to initilized*/
        if(att_flag == false){
            Eigen::Vector3d att(msg.roll/100000, msg.pitch/100000, msg.heading/100000);
            /* NED -> ENU */
            att(2) = 360 - att(2);
            double tmp = att(1);
            att(1) = att(0);
            att(0) = tmp;

            /* 0 ~ 360 -> -180 ~ 180 */
            for (int i = 0; i < 3; i++){
                if (att(i) > 180 && att(i) < 360){
                    att(i) = att(i) - 360;
                }
            }
            att = att*DEG_TO_RAD; 
            // std::cout << "att" << std::endl << att << std::endl;
            state_vector_.R_b_l = coordinate_mat_transformation::Rotation_matrix(att);
            state_vector_.att_l = att;
        }
    }
    else if(ins_config_.mode == 0 || att_flag == false){
        Eigen::Vector3d att(msg.roll/100000, msg.pitch/100000, msg.heading/100000);
        /* NED -> ENU */
        att(2) = 360 - att(2);
        double tmp = att(1);
        att(1) = att(0);
        att(0) = tmp;

        /* 0 ~ 360 -> -180 ~ 180 */
        for (int i = 0; i < 3; i++){
            if (att(i) > 180 && att(i) < 360){
                att(i) = att(i) - 360;
            }
        }
        att = att*DEG_TO_RAD; 
        // std::cout << "att" << std::endl << att << std::endl;
        state_vector_.R_b_l = coordinate_mat_transformation::Rotation_matrix(att);
        state_vector_.att_l = att;
    }
    else if(ins_config_.mode == 1){/*Don't update */}
    if(ins_config_.mode == 0 && fix_flag == true && ins_config_.fix_type != 3) send_tf(state_vector_.r_l, state_vector_.att_l, "ublox");
    att_flag = true;
}

void Ins_mechanization::Novatelfixcallback(const novatel_gps_msgs::Inspva& msg){
    if(ins_config_.mode == 0 || fix_flag == false){
        if (ins_config_.novatel_count >= 100){
            ins_config_.novatel_count = 0;
            Eigen::Vector3d pos(msg.latitude, msg.longitude, msg.height);
            Eigen::Vector3d vel_l(msg.east_velocity, msg.north_velocity, msg.up_velocity);
            Eigen::Vector3d att(msg.roll, msg.pitch, msg.azimuth);
            /* NED -> ENU */
            att(2) = 360 - att(2);
            double tmp = att(1);
            att(1) = att(0);
            att(0) = tmp;
            /* 0 ~ 360 -> -180 ~ 180 */
            for (int i = 0; i < 3; i++){
                if (att(i) > 180 && att(i) < 360){
                    att(i) = att(i) - 360;
                }
            }
            att = att*DEG_TO_RAD; 

            /* If odometer is on, velocity isn't updated with Novatel */
            if(ins_config_.odometer == 0){
                state_vector_.v_l = vel_l;
                vel_flag = true;
            }
            state_vector_.r_l = pos;
            state_vector_.R_b_l = coordinate_mat_transformation::Rotation_matrix(att);
            state_vector_.att_l = att;
            fix_flag = att_flag = true;
            ins_config_.transform2baselink = false;

            // std::cout << "\033[33m" << "att" << std::endl << att << "\033[0m" << std::endl;
            if(ins_config_.mode == 0) send_tf(state_vector_.r_l, state_vector_.att_l, "novatel");
        }
        else ins_config_.novatel_count++;
    }
    else if(ins_config_.mode == 1){/* Don't update */}
}
void Ins_mechanization::uwbfixcallback(const uwb_ins_eskf_msgs::uwbFIX& msg){
    static int control_rate = 3;
    if(control_rate == 3){
        if(ins_config_.mode == 0 || fix_flag == false){
            Eigen::Vector3d pos(msg.latitude, msg.longitude, msg.altitude);
            Eigen::Vector3d vel_l(msg.velocity_e, msg.velocity_n, msg.velocity_u);
            Eigen::Vector3d att(msg.att_e, msg.att_n, msg.att_u);

            /* 0 ~ 360 -> -180 ~ 180 */
            for (int i = 0; i < 3; i++){
                if (att(i) > 180 && att(i) < 360){
                    att(i) = att(i) - 360;
                }
            }
            att = att*DEG_TO_RAD; 

            /* If odometer is on, velocity isn't updated with UWB */
            if(ins_config_.odometer == 0){
                state_vector_.v_l = vel_l;
                vel_flag = true;
            }
            state_vector_.r_l = pos;
            state_vector_.R_b_l = coordinate_mat_transformation::Rotation_matrix(att);
            state_vector_.att_l = att;
            fix_flag = att_flag = true;
            ins_config_.transform2baselink = false;
        }
        else if(ins_config_.mode == 1){
            /* Don't update */
        }
    }
    if(control_rate-- == 3) control_rate = 3;
}

void Ins_mechanization::fusionfixcallback(const uwb_ins_eskf_msgs::fusionFIX& msg){
    fix_flag = vel_flag = att_flag = true;

    Eigen::Vector3d pos(msg.latitude, msg.longitude, msg.altitude);
    Eigen::Vector3d vel_l(msg.velocity_e, msg.velocity_n, msg.velocity_u);
    Eigen::Vector3d att(msg.att_e, msg.att_n, msg.att_u);

    /* 0 ~ 360 -> -180 ~ 180 */
    for (int i = 0; i < 3; i++){
        if (att(i) > 180 && att(i) < 360){
            att(i) = att(i) - 360;
        }
    }
    att = att*DEG_TO_RAD; 

    state_vector_.r_l = pos;
    state_vector_.v_l = vel_l;
    state_vector_.R_b_l = coordinate_mat_transformation::Rotation_matrix(att);
    state_vector_.att_l = att;
    // std::cout << "\033[33m" << "att" << std::endl << att << "\033[0m" << std::endl;
}

void Ins_mechanization::Imucallback(const sensor_msgs::Imu& msg){
    double current_time_tag = msg.header.stamp.toSec();
    Eigen::Vector3d acc_raw(msg.linear_acceleration.x, msg.linear_acceleration.y, msg.linear_acceleration.z);
    Eigen::Vector3d gyro_raw(msg.angular_velocity.x, msg.angular_velocity.y, msg.angular_velocity.z);
    Imu_data_calibration(acc_raw, gyro_raw);

    /* /imu_meas topic is 250 hz or 100 hz */
    static int count = 1;
    if ((current_time_tag - mech_variables_.compute_time_tag) > 0.005 && count >= 1){
        count = 1;
        compute();
        mech_variables_.compute_time_tag = current_time_tag;
    }
    else count++;
}

void Ins_mechanization::Odometercallback(const geometry_msgs::TwistStamped& msg){
    if(ins_config_.mode == 0 || vel_flag == false){
        state_vector_.v_forward = msg.twist.linear.x;
        vel_flag = true;
    }
    else if(ins_config_.mode == 1){
        state_vector_.v_forward = msg.twist.linear.x;
        vel_flag = true;
        /* Don't update */
    }
}


void Ins_mechanization::transform2baselink(){
    Eigen::Vector3d ref_lla(NCKUEE_LATITUDE, NCKUEE_LONGITUDE, NCKUEE_HEIGHT);
    Eigen::Vector3d r_b;
    if(ins_config_.fix_type == 0){
        r_b = ins_config_.uwb_b;
    }
    else if(ins_config_.fix_type == 1 || ins_config_.fix_type == 2){
        r_b = ins_config_.gnss_b;
    }

    Eigen::Vector3d r_enu = coordinate_mat_transformation::lla2enu(state_vector_.r_l ,ref_lla);
    Eigen::Vector3d baselink_enu = r_enu - state_vector_.R_b_l*r_b;
    Eigen::VectorXd baselink_lla = coordinate_mat_transformation::enu2Geodetic(baselink_enu ,ref_lla);
    state_vector_.r_l = baselink_lla;
    ins_config_.transform2baselink = true;
    // std::cout << "r_enu: " << std::endl << r_enu << std::endl;
    // std::cout << "baselink_enu: " << std::endl << baselink_enu << std::endl;
}

void Ins_mechanization::Imu_data_calibration(Eigen::Vector3d acc_raw, Eigen::Vector3d gyro_raw){
    Eigen::Vector3d acc_calibrated;
    Eigen::Matrix3d acc_scale;
    acc_scale << imu_correction_.acc_scale(0),                          0,                           0,
                                           0,imu_correction_.acc_scale(1),                           0,
                                           0,                          0, imu_correction_.acc_scale(2);   
    acc_calibrated = (Eigen::MatrixXd::Identity(3,3) + acc_scale).inverse() * (acc_raw - imu_correction_.acc_bias);
    // acc_calibrated = (acc_raw - imu_correction_.acc_bias);

    Eigen::Vector3d gyro_calibrated;
    Eigen::Matrix3d gyro_scale;
    gyro_calibrated = (gyro_raw - imu_correction_.gyro_bias);

    if(ins_config_.b_frame == 0){
        /* If body frame is ENU */
        mech_variables_.f_b = acc_calibrated;
        mech_variables_.w_ib_b = gyro_calibrated;
    }
    else if(ins_config_.b_frame == 1){
        /* If body frame is NWU */
        mech_variables_.f_b(0) = -acc_calibrated(1);
        mech_variables_.f_b(1) = acc_calibrated(0);
        mech_variables_.f_b(2) = acc_calibrated(2);
        mech_variables_.w_ib_b(0) = -gyro_calibrated(1);
        mech_variables_.w_ib_b(1) = gyro_calibrated(0);
        mech_variables_.w_ib_b(2) = gyro_calibrated(2);
    }
    
}

void Ins_mechanization::Initialize_state(){
    Eigen::Vector3d unit(0,0,0);
    Eigen::Vector4d unit_q(0,0,0,1);

    state_vector_.r_l = unit;
    state_vector_.v_l = unit;
    state_vector_.R_b_l = Eigen::MatrixXd::Identity(3,3);
    state_vector_.att_l = unit;
    mech_variables_.f_b = unit;
    mech_variables_.w_ib_b = unit;
    mech_variables_.r_dot_l_now = unit;
    mech_variables_.r_dot_l_past = unit;
    mech_variables_.v_dot_l_now = unit;
    mech_variables_.v_dot_l_past = unit;
    mech_variables_.R_dot_b_l_now = Eigen::MatrixXd::Identity(3,3);
    mech_variables_.R_dot_b_l_past = Eigen::MatrixXd::Identity(3,3);
    mech_variables_.Q_dot_now = unit_q;
    mech_variables_.Q_dot_past = unit_q;
}

void Ins_mechanization::Attitude_update(){
    /* inference the current attitude */
    Q_dot_update(mech_variables_, state_vector_);
    // std::cout << "att_l" << std::endl << state_vector_.att_l << std::endl;
    // std::cout << "q" << std::endl << coordinate_mat_transformation::att_Q_transform(state_vector_.att_l) << std::endl;
    // std::cout << "q_dot" << std::endl << mech_variables_.Q_dot_past * SAMPLING_TIME << std::endl;
    Eigen::VectorXd Q_now = coordinate_mat_transformation::att_Q_transform(state_vector_.att_l) + mech_variables_.Q_dot_past * SAMPLING_TIME;
    // std::cout << "Q_now" << std::endl << Q_now << std::endl;
    // std::cout << "Q_dot_past" << std::endl << mech_variables_.Q_dot_past * SAMPLING_TIME << std::endl;

    /* update current attitude and rotation matrix */
    Eigen::VectorXd att_now = coordinate_mat_transformation::Q_att_transform(Q_now);
    state_vector_.att_l = att_now;
    state_vector_.R_b_l = coordinate_mat_transformation::Rotation_matrix(att_now);
    // std::cout << "att_now" << std::endl << att_now << std::endl;
}

void Ins_mechanization::Velocity_update(){
    /* inference the current velocity */
    v_dot_l_update(mech_variables_, state_vector_);
    Eigen::VectorXd vel_now = state_vector_.v_l + 0.5 * (mech_variables_.v_dot_l_now + mech_variables_.v_dot_l_past) * SAMPLING_TIME;

    /* update current velocity */
    state_vector_.v_l = vel_now;
    // std::cout << "vel_now" << std::endl << vel_now << std::endl;
}

void Ins_mechanization::Position_update(){
    /* inference the current position*/
    Eigen::Vector3d r_now;
    r_dot_l_update(mech_variables_, state_vector_);
    r_now(0) = state_vector_.r_l(0) + 0.5 * ((mech_variables_.r_dot_l_now(0) + mech_variables_.r_dot_l_past(0))*RAD_TO_DEG) * SAMPLING_TIME; // latitude
    r_now(1) = state_vector_.r_l(1) + 0.5 * ((mech_variables_.r_dot_l_now(1) + mech_variables_.r_dot_l_past(1))*RAD_TO_DEG) * SAMPLING_TIME; // longitude
    r_now(2) = state_vector_.r_l(2) + 0.5 * (mech_variables_.r_dot_l_now(2) + mech_variables_.r_dot_l_past(02)) * SAMPLING_TIME; // height

    /* update current position */
    state_vector_.r_l = r_now;
    Eigen::Vector3d pos_enu = coordinate_mat_transformation::lla2enu(r_now, nckuee);
    std::cout << std::fixed << std::setprecision(10);
    // std::cout << "ENU now: " << std::endl << pos_enu << std::endl;
}

void Ins_mechanization::send_tf(){
    tf::Transform transform;
    tf::Quaternion current_q;
    ros::Time now = ros::Time::now();
    
    Eigen::Vector3d now_enu = coordinate_mat_transformation::lla2enu(state_vector_.r_l, nckuee);
    // std::cout << "\033[33m" << "ENU:" << std::endl << now_enu << "\033[0m" << std::endl;
    // std::cout << "\033[33m" << "yaw: " << -state_vector_.att_l(2) << "\033[0m" << std::endl;

    current_q.setRPY(state_vector_.att_l(0), 
                    state_vector_.att_l(1), 
                    state_vector_.att_l(2)); // refer to "map" frame (ENU)

    transform.setOrigin(tf::Vector3(now_enu(0), now_enu(1), now_enu(2)));
    transform.setRotation(current_q);
    br_.sendTransform(tf::StampedTransform(transform, now, "/map", "ins_baselink"));
}

void Ins_mechanization::send_tf(Eigen::Vector3d now_lla, Eigen::Vector3d now_att, std::string frame){
    tf::Transform transform;
    tf::Quaternion current_q;
    ros::Time now = ros::Time::now();
    
    Eigen::Vector3d now_enu = coordinate_mat_transformation::lla2enu(now_lla, nckuee);
    // std::cout << "\033[33m" << "ENU:" << std::endl << now_enu << "\033[0m" << std::endl;

    current_q.setRPY(now_att(0), now_att(1), now_att(2)); // refer to "map" frame (ENU)

    transform.setOrigin(tf::Vector3(now_enu(0), now_enu(1), now_enu(2)));
    transform.setRotation(current_q);
    br_.sendTransform(tf::StampedTransform(transform, now, "/map", frame));
}

void Ins_mechanization::Publish_ins(){
    uwb_ins_eskf_msgs::InsFIX msg;
    ros::Time now = ros::Time::now();

    msg.header.stamp = now;
    msg.header.frame_id = "local";
    msg.latitude = state_vector_.r_l(0);
    msg.longitude = state_vector_.r_l(1);
    msg.altitude = state_vector_.r_l(2);
    msg.velocity_e = state_vector_.v_l(0);
    msg.velocity_n = state_vector_.v_l(1);
    msg.velocity_u = state_vector_.v_l(2);
    // 0 ~ 360 -> -180 ~ 180
    Eigen::Vector3d att = state_vector_.att_l*RAD_TO_DEG;
    for (int i = 0; i < 3; i++){
        if (att(i) > 180 && att(i) < 360){
            att(i) = att(i) - 360;
        }
    }
    msg.att_e = att(0);
    msg.att_n = att(1);
    msg.att_u = att(2);
    Eigen::Vector3d f_l = state_vector_.R_b_l * mech_variables_.f_b;
    msg.f_e = f_l(0);
    msg.f_n = f_l(1);
    msg.f_u = f_l(2);

    pub_ins_fix_.publish(msg);
}

void Ins_mechanization::compute(){
    if(fix_flag == false || vel_flag == false || att_flag == false){
        std::cout << "\033[33m" << "Waiting for initialization !" << "\033[0m" << std::endl;
    }
    else{
        Attitude_update();
        /* 
        If Odometer is on
        Update velocity with GNSS speed / Odometer measuenment
        */
        if(ins_config_.odometer != 0){
            Eigen::Vector3d vel_b(0, state_vector_.v_forward, 0); //enu
            state_vector_.v_l = state_vector_.R_b_l*vel_b;
            // std::cout << "\033[33m" << "state_vector_.v_l" << std::endl << state_vector_.v_l << "\033[0m" << std::endl;
            // std::cout << "\033[33m" << "v_forward" << std::endl << state_vector_.v_forward << "\033[0m" << std::endl;
        }
        Velocity_update();
        /* Check if position is tranform to baselink */
        if(!ins_config_.transform2baselink){
            transform2baselink();
            std::cout << "\033[33m" << "Transform to baselink !" << "\033[0m" << std::endl;
        }
        Position_update();
        Publish_ins();
        send_tf();
    }
}
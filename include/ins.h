#ifndef INS_H_
#define INS_H_
#include "ros/ros.h"
#include <sensor_msgs/Imu.h>
#include <geometry_msgs/TwistWithCovarianceStamped.h>
#include <geometry_msgs/TwistStamped.h>
#include <ublox_msgs/NavATT.h>
#include <ublox_msgs/NavPVT.h>
#include <sensor_msgs/NavSatFix.h>
#include <novatel_gps_msgs/Inspva.h>
#include <uwb_ins_eskf_msgs/fusionFIX.h>
#include <uwb_ins_eskf_msgs/InsFIX.h>
#include <uwb_ins_eskf_msgs/uwbFIX.h>

#include <tf/tf.h>
#include <tf/transform_broadcaster.h>
#include <tf/transform_datatypes.h>

#define _USE_MATH_DEFINES
#include <cmath>
#include <Eigen/Dense>

#define EARTH_SEMIMAJOR 6378137 // meter
#define EARTH_ROTATION_RATE 0.00007292115 // rad/s
#define FLATTENING 1/298.25722563
#define ECCENTRICITY sqrt(FLATTENING*(2-FLATTENING))
#define NCKUEE_LATITUDE 22.99665875 //  degree
#define NCKUEE_LONGITUDE 120.222584889 //  degree
#define NCKUEE_HEIGHT 98.211 // meter
#define ACC_BIAS_X -0.348513345960806
#define ACC_BIAS_Y -0.26021251227717
#define ACC_BIAS_Z 0.132337782341719
#define ACC_SCALE_X -0.00426378745053201
#define ACC_SCALE_Y 0.000725755116990245
#define ACC_SCALE_Z 0.00263712381843959
#define GYRO_BIAS_X 0.00710659149934062
#define GYRO_BIAS_Y 0.00211909908717263
#define GYRO_BIAS_Z -0.0000592951099686292
#define GYRO_SCALE_X -2.36707629594559
#define GYRO_SCALE_Y -0.490347919324706
#define GYRO_SCALE_Z -0.686283178454847
#define SAMPLING_TIME 0.02 //sec
#define DEG_TO_RAD M_PI/180
#define RAD_TO_DEG 180/M_PI

// l: local-level frame (navigation frame) ENU
// b: body frame
// e: earth-center-earth-fixed frame
// i: inertia frame

namespace INS
{
    typedef struct INS_config{
        int mode;
        int fix_type;
        int odometer;
        int novatel_count = 100;
        int b_frame;
        Eigen::Vector3d gnss_b;
        Eigen::Vector3d uwb_b;
        bool transform2baselink = false;
    }config;

    typedef struct imu_calibration_parameter{
        Eigen::VectorXd acc_bias;
        Eigen::VectorXd acc_scale;
        Eigen::VectorXd gyro_bias;
        Eigen::VectorXd gyro_scale;
    }imu_calibration;

    typedef struct ins_mechanization_variable{
        double compute_time_tag = 0; // Time tag of computing ins result
        Eigen::VectorXd f_b; // measurement: accelerometer
        Eigen::VectorXd w_ib_b; // measurement: gyroscope
        Eigen::VectorXd r_dot_l_now; // at present
        Eigen::VectorXd r_dot_l_past; // last moment
        Eigen::VectorXd v_dot_l_now;
        Eigen::VectorXd v_dot_l_past;
        Eigen::MatrixXd R_dot_b_l_now;
        Eigen::MatrixXd R_dot_b_l_past;
        Eigen::VectorXd Q_dot_now;
        Eigen::VectorXd Q_dot_past;
    }mechanization;

    typedef struct ins_state_vector{
        Eigen::VectorXd r_l; // latitude longitude height (deg)
        Eigen::VectorXd v_l; // east north up (m/s)
        double v_forward; // (m/s)
        Eigen::MatrixXd R_b_l; //rotation matrix
        Eigen::VectorXd att_l; // roll pitch yaw (rad)
    }state;

    class Ins_mechanization
    {
        private:
            imu_calibration imu_correction_;
            state state_vector_;
            mechanization mech_variables_;
            ros::Publisher pub_ins_fix_;
            tf::TransformBroadcaster br_;
            config ins_config_;
        public:
            Ins_mechanization(ros::Publisher pub_fix, config ins_config);
            double gravity(double lat_rad);
            double earth_radius_along_meridian(double lat_rad);
            double earth_radius_along_prime_vertical(double lat_rad);
            Eigen::VectorXd w_el_l_update(Eigen::VectorXd vel_enu, double lat_rad, double h);
            Eigen::VectorXd w_ie_l_update(double lat_rad);
            Eigen::MatrixXd omega_el_l_update(state& vector);
            Eigen::MatrixXd omega_ie_l_update(state& vector);
            Eigen::MatrixXd omega_il_b_update(state& vector, int dim);
            Eigen::MatrixXd D_inverse_update(state& vector);
            void R_dot_l_b_update(mechanization& variables, state& vec);
            void Q_dot_update(mechanization& variables, state& vec);
            void v_dot_l_update(mechanization& variables, state& vec);
            void r_dot_l_update(mechanization& variables, state& vec);

            void GNSSfixcallback(const sensor_msgs::NavSatFix& msg);
            void GNSSvelcallback(const geometry_msgs::TwistWithCovarianceStamped& msg);
            void GNSSspeedcallback(const ublox_msgs::NavPVT& msg);
            void GNSSattcallback(const ublox_msgs::NavATT& msg);
            void Novatelfixcallback(const novatel_gps_msgs::Inspva& msg);
            void uwbfixcallback(const uwb_ins_eskf_msgs::uwbFIX& msg);
            void fusionfixcallback(const uwb_ins_eskf_msgs::fusionFIX& msg);
            void Imucallback(const sensor_msgs::Imu& msg);
            void Odometercallback(const geometry_msgs::TwistStamped& msg);
            
            void transform2baselink();
            void Imu_data_calibration(Eigen::Vector3d acc_raw, Eigen::Vector3d gyro_raw);
            void Initialize_state();
            void Attitude_update();
            void Velocity_update();
            void Position_update();
            void send_tf();
            void send_tf(Eigen::Vector3d now_lla, Eigen::Vector3d now_att, std::string frame);
            void Publish_ins();
            void compute();
    };
}
#endif
#include "ins.h"

int main(int argc, char **argv) {
    ros::init(argc, argv, "ins");

    ros::NodeHandle n;
    ros::NodeHandle nh("~");

    // Publishers of ins solution
    ros::Publisher pub_ins_fix = n.advertise<ins::InsFIX>("/ins/fix", 1);

    INS::Ins_mechanization ublox_ins(pub_ins_fix);

    // Subscribers of GNSS data
    ros::Subscriber sub_att = n.subscribe("/ublox_f9k/navatt", 1, &INS::Ins_mechanization::GNSSattcallback, &ublox_ins);
    ros::Subscriber sub_pos = n.subscribe("/ublox_f9k/fix", 1, &INS::Ins_mechanization::GNSSfixcallback, &ublox_ins);
    ros::Subscriber sub_vel = n.subscribe("/ublox_f9k/fix_velocity", 1, &INS::Ins_mechanization::GNSSvelcallback, &ublox_ins);
    ros::Subscriber sub_imu = n.subscribe("/ublox_f9k/imu_meas", 1, &INS::Ins_mechanization::Imucallback, &ublox_ins);

    ublox_ins.Initialize_state();
    ros::spin();

    return 0;
}
#include "ins.h"

int main(int argc, char **argv) {
    ros::init(argc, argv, "ins");
    INS::config ins_config;

    ros::NodeHandle n;
    ros::NodeHandle nh("~");

    nh.param("gnss_fix_type", ins_config.fix_type, 0);
    nh.param("fix_mode", ins_config.mode, 0);

    // Publishers of ins solution
    ros::Publisher pub_ins_fix = n.advertise<uwb_ins_eskf_msgs::InsFIX>("/ins/fix", 1);

    INS::Ins_mechanization ublox_ins(pub_ins_fix, ins_config);

    // Subscribers of GNSS data
    ros::Subscriber sub_att;
    ros::Subscriber sub_pos;
    ros::Subscriber sub_vel;
    ros::Subscriber sub_novatel;
    ros::Subscriber sub_eskf;
    if(ins_config.fix_type == 0){
        sub_att = n.subscribe("/ublox_f9k/navatt", 1, &INS::Ins_mechanization::GNSSattcallback, &ublox_ins);
        sub_pos = n.subscribe("/ublox_f9k/fix", 1, &INS::Ins_mechanization::GNSSfixcallback, &ublox_ins);
        sub_vel = n.subscribe("/ublox_f9k/fix_velocity", 1, &INS::Ins_mechanization::GNSSvelcallback, &ublox_ins);
    }
    else if(ins_config.fix_type == 1){
        sub_novatel = n.subscribe("/novatel/inspva", 1, &INS::Ins_mechanization::Novatelfixcallback, &ublox_ins);
    }
    
    if(ins_config.mode == 1){
        sub_eskf = n.subscribe("/fusion/fix", 1, &INS::Ins_mechanization::fusionfixcallback, &ublox_ins);
    }
    
    ros::Subscriber sub_imu = n.subscribe("/ublox_f9k/imu_meas", 1, &INS::Ins_mechanization::Imucallback, &ublox_ins);

    ublox_ins.Initialize_state();

    ros::spin();

    return 0;
}
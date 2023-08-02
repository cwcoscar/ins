#include "ins.h"

int main(int argc, char **argv) {
    ros::init(argc, argv, "ins");
    INS::config ins_config;

    ros::NodeHandle n;
    ros::NodeHandle nh("~");

    nh.param("gnss_fix_type", ins_config.fix_type, 0);
    nh.param("odometer", ins_config.odometer, 0);
    nh.param("fix_mode", ins_config.mode, 0);
    nh.param("body_frame", ins_config.b_frame, 1);
    nh.param("antenna_x", ins_config.gnss_b(0), 0.00);
    nh.param("antenna_y", ins_config.gnss_b(1), 0.25);
    nh.param("antenna_z", ins_config.gnss_b(2), 0.43);
    nh.param("A0_x", ins_config.uwb_b(0), 0.00);
    nh.param("A0_y", ins_config.uwb_b(1), 0.05);
    nh.param("A0_z", ins_config.uwb_b(2), 0.43);

    /* Publishers of ins solution */
    ros::Publisher pub_ins_fix = n.advertise<uwb_ins_eskf_msgs::InsFIX>("/ins/fix", 1);

    INS::Ins_mechanization ublox_ins(pub_ins_fix, ins_config);

    /* Subscribers of GNSS data */
    ros::Subscriber sub_ublox_att;
    ros::Subscriber sub_ublox_pos;
    ros::Subscriber sub_ublox_vel;
    ros::Subscriber sub_ublox_speed;
    ros::Subscriber sub_odo_speed;
    ros::Subscriber sub_novatel;
    ros::Subscriber sub_uwb;
    ros::Subscriber sub_eskf;
    if(ins_config.fix_type == 0){
        sub_uwb = n.subscribe("/uwb_position/A0", 1, &INS::Ins_mechanization::uwbfixcallback, &ublox_ins);
    }
    else if(ins_config.fix_type == 1){
        sub_ublox_att = n.subscribe("/ublox_f9k/navatt", 1, &INS::Ins_mechanization::GNSSattcallback, &ublox_ins);
        sub_ublox_pos = n.subscribe("/ublox_f9k/fix", 1, &INS::Ins_mechanization::GNSSfixcallback, &ublox_ins);
        sub_ublox_vel = n.subscribe("/ublox_f9k/fix_velocity", 1, &INS::Ins_mechanization::GNSSvelcallback, &ublox_ins);
    }
    else if(ins_config.fix_type == 2){
        sub_novatel = n.subscribe("/novatel/inspva", 1, &INS::Ins_mechanization::Novatelfixcallback, &ublox_ins);
    }
    else if(ins_config.fix_type == 3){
        sub_ublox_att = n.subscribe("/ublox_f9k/navatt", 1, &INS::Ins_mechanization::GNSSattcallback, &ublox_ins);
        sub_ublox_pos = n.subscribe("/ublox_f9k/fix", 1, &INS::Ins_mechanization::GNSSfixcallback, &ublox_ins);
        sub_ublox_vel = n.subscribe("/ublox_f9k/fix_velocity", 1, &INS::Ins_mechanization::GNSSvelcallback, &ublox_ins);
    }

    if(ins_config.odometer == 0){

    }
    else if (ins_config.odometer == 1){
        sub_odo_speed = n.subscribe("/can_velocity", 1, &INS::Ins_mechanization::Odometercallback, &ublox_ins);
    }
    else if (ins_config.odometer == 2){
        sub_ublox_speed = n.subscribe("/ublox_f9k/navpvt", 1, &INS::Ins_mechanization::GNSSspeedcallback, &ublox_ins);
    }
    
    if(ins_config.mode == 1){
        sub_eskf = n.subscribe("/fusion/baselink/fix", 1, &INS::Ins_mechanization::fusionfixcallback, &ublox_ins);
    }
    
    ros::Subscriber sub_imu = n.subscribe("/ublox_f9k/imu_meas", 1, &INS::Ins_mechanization::Imucallback, &ublox_ins);

    ublox_ins.Initialize_state();

    ros::spin();

    return 0;
}
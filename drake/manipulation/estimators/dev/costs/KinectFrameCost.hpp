#ifndef KINECT_FRAME_COST_H
#define KINECT_FRAME_COST_H

#include <stdexcept>
#include <iostream>
#include "ManipulationTrackerCost.hpp"
#include "drake/multibody/rigid_body_tree.h"
#include <lcm/lcm-cpp.hpp>
#include <memory>
#include <mutex>
#include "yaml-cpp/yaml.h"

#include "lcmtypes/bot_core/rigid_transform_t.hpp"
#include "lcmtypes/bot_core/raw_t.hpp"
#include "lcmtypes/kinect/frame_msg_t.hpp"
#include "lcmtypes/bot_core/rigid_transform_t.hpp"
#include "lcmtypes/bot_core/image_t.hpp"
#include <kinect/kinect-utils.h>
#include <mutex>
#include <bot_lcmgl_client/lcmgl.h>
#include <bot_frames/bot_frames.h>
#include <bot_param/param_client.h>
#include <lcmtypes/bot_core/pointcloud_t.hpp>
#include <lcmtypes/bot_core/image_t.hpp>

#include "drake/systems/sensors/camera_info.h"

#include "drake/systems/sensors/image.h"

using drake::systems::sensors::CameraInfo;

/**
 * Handling Kinect cost (point cloud and depth image)
 * HACK: This is going to be disassociated from a Kinect to simplify processing.
 */
class KinectFrameCost : public ManipulationTrackerCost {
public:
  typedef Eigen::Matrix3Xd PointCloud;
  typedef drake::systems::sensors::ImageDepth32F DepthImage;

  KinectFrameCost(std::shared_ptr<RigidBodyTreed> robot_,
                  std::shared_ptr<lcm::LCM> lcm_,
                  YAML::Node config,
                  const CameraInfo* camera_info);
  ~KinectFrameCost() {};
  bool constructCost(ManipulationTracker * tracker, const Eigen::VectorXd x_old, Eigen::MatrixXd& Q, Eigen::VectorXd& f, double& K);

  void initBotConfig(const char* filename);
  int get_trans_with_utime(std::string from_frame, std::string to_frame,
                               long long utime, Eigen::Isometry3d & mat);
  void handleSavePointcloudMsg(const lcm::ReceiveBuffer* rbuf,
                           const std::string& chan,
                           const bot_core::raw_t* msg);
  void readDepthImageAndPointCloud(const DepthImage& depth_image,
                                   const PointCloud& point_cloud);

  void handleCameraOffsetMsg(const lcm::ReceiveBuffer* rbuf,
                           const std::string& chan,
                           const bot_core::rigid_transform_t* msg);

   // bounds to cut down point cloud, in world coords
  struct BoundingBox
  {
      double xMin = -100.;
      double xMax = 100.;
      double yMin = -100.;
      double yMax = 100.;
      double zMin = -100.;
      double zMax = 100.;
  };
  void setBounds(BoundingBox bounds) { pointcloud_bounds = bounds; }

private:
  double downsample_amount = 10.0;
  int input_num_pixel_cols = 640;
  int input_num_pixel_rows = 480;
  int num_pixel_cols, num_pixel_rows;

  double icp_var = INFINITY;
  double free_space_var = INFINITY;
  double max_considered_icp_distance = 0.075;
  double min_considered_joint_distance = 0.03;
  double timeout_time = 0.5;
  double max_scan_dist = 10.0;
  bool verbose = false;
  bool verbose_lcmgl = false;

  // We operate in one of four modes:
  // 1) Kinect has a body whose pose is learned. Have_camera_body is true. have_camera_body=true,world_frame=true
  // 2) Kinect listens to pose on LCM, have_camera_body=false but world_frame=true
  // 3) Kinect has hard coded world frame if world_frame=true. have_hardcoded_kinect2world = true
  // 3) Kinect frame is world frame. have_camera_body and world_frame = false
  bool have_camera_body_ = false;
  std::string camera_body_name_;
  int camera_body_ind_;
  bool world_frame = true;

  std::shared_ptr<lcm::LCM> lcm;
  std::shared_ptr<RigidBodyTreed> robot;
  KinematicsCache<double> robot_kinematics_cache;
  int nq;

  bot_lcmgl_t* lcmgl_lidar_ = NULL;
  bot_lcmgl_t* lcmgl_icp_ = NULL;
  bot_lcmgl_t* lcmgl_measurement_model_ = NULL;
  BotParam* botparam_ = NULL;
  BotFrames* botframes_ = NULL;

  BoundingBox pointcloud_bounds;

  std::mutex latest_cloud_mutex;
  std::mutex camera_offset_mutex;
  Eigen::Isometry3d kinect2world_;
  bool have_hardcoded_kinect2world_ = false;
  Eigen::Isometry3d hardcoded_kinect2world_;

  Eigen::Matrix3Xd latest_cloud;
//  KinectCalibration* kcal;
  Eigen::MatrixXd latest_depth_image;
  Eigen::Matrix3Xd latest_color_image;
  Eigen::Matrix3Xd raycast_endpoints;
  const CameraInfo* camera_info_;

  double lastReceivedTime;
  double last_got_kinect_frame;
};

#endif

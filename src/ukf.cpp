#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  n_x_ = 5;


  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

   H_ = MatrixXd(2, 5);
      H_ << 1, 0, 0, 0, 0,
            0, 1, 0, 0, 0;


  R_ = MatrixXd(2, 2);
      R_ << std_laspx_*std_laspx_, 0,
            0, std_laspy_*std_laspy_ ;

  
    x_ << 1,1,1,1,1;
    P_ = Eigen::MatrixXd(5,5);
  
    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0, 
          0, 0 , 1, 0, 0,
          0, 0,  0,  0.0225, 0,
          0, 0,  0,  0, .0225;

  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
    // Radar measurement noise standard deviation radius change in m/s
  
   // Augmented state dimension
   n_aug_ = 7;


  // Weights of sigma points
   weights_ = VectorXd(2*n_aug_+1);

  // State dimension
 
 
  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

  time_us_ = 0.0;

  std::cout << " C'tor initialize" << std::endl;

  is_initialized_ = false;

}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  if (meas_package.sensor_type_ == MeasurementPackage::LASER || 
      meas_package.sensor_type_ == MeasurementPackage::RADAR
  )
  if (!is_initialized_) {

    if (meas_package.sensor_type_ == meas_package.SensorType::LASER)
    {
        x_ << meas_package.raw_measurements_(0), 
                meas_package.raw_measurements_(1), 
                0, 
                0,
                0;

    } else {
        double rho = meas_package.raw_measurements_(0);
        double si =  meas_package.raw_measurements_(1);
        double rhodot = meas_package.raw_measurements_(2);
        double vx = rhodot*cos(si);
        double vy = rhodot*sin(si);
        x_ << rho*cos(si), 
        rho*sin(si),
        0,
        0,
        0;
    }
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    std::cout <<  " Initialized " <<std::endl;
    return;
 
  }
  float dt = (meas_package.timestamp_ - time_us_)/1000000.0;
  time_us_ = meas_package.timestamp_;
  Prediction(dt);
  if (meas_package.sensor_type_ == meas_package.SensorType::LASER)
  {
    UpdateLidar(meas_package);
  } else {
    UpdateRadar(meas_package);
  }
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
  // Initialize

  Eigen::VectorXd x_aug = Eigen::VectorXd(7);

  // create augmented state covariance
  Eigen::MatrixXd P_aug = Eigen::MatrixXd(n_aug_, n_aug_);

  // create sigma point matrix
  Eigen::MatrixXd Xsig_aug = Eigen::MatrixXd(n_aug_, 2 * n_aug_ + 1);

  /**
   * Student part begin
   */
 
  // create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  // create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  // create square root matrix
  Eigen::MatrixXd L = P_aug.llt().matrixL();

  // create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; ++i) {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  //  Predict sigma points for the augmented points
  // augmented includes the process noise.
  for (int i = 0; i< 2*n_aug_+1; ++i) {
    // extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    // predicted state values
    double px_p, py_p;

    // avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    } else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    // add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    // write predicted sigma point into right column
    Xsig_pred(0,i) = px_p;
    Xsig_pred(1,i) = py_p;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
  }


  // Predict mean and covariance
  // create vector for weights
  Eigen::VectorXd weights = Eigen::VectorXd(2*n_aug_+1);
   // create vector for predicted state
  Eigen::VectorXd x = Eigen::VectorXd(n_x_);

  // create covariance matrix for prediction
  Eigen::MatrixXd P = Eigen::MatrixXd(n_x_, n_x_);

  double weight_0 = lambda_/(lambda_ + n_aug_);
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; ++i) {  // 2n+1 weights
    double weight = 0.5/(n_aug_ + lambda_);
    weights(i) = weight;
  }

  // predicted state mean
  x.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
    x = x + weights(i) * Xsig_pred.col(i);
  }

  // predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
    // state difference
    Eigen::VectorXd x_diff = Xsig_pred.col(i) - x;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P = P + weights(i) * x_diff * x_diff.transpose() ;
  }
  
  x_ = x;
  P_ = P;

  std::cout << " Prediction completed   " << std::endl;
  

}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  int n_z = 2;
  Eigen::VectorXd z = Eigen::VectorXd(n_z);
  z << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1);
  std::cout << " Update Lidar " << std::endl;

  Eigen::VectorXd z_pred = H_ * x_;
  Eigen::VectorXd y = z - z_pred;
  Eigen::MatrixXd Ht = H_.transpose();
  Eigen::MatrixXd S = H_ * P_ * Ht + R_;
  Eigen::MatrixXd K = P_*Ht *S.inverse();
 
  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  while (x_(3)> M_PI) x_(3)-=2.*M_PI;
  while (x_(3)<-M_PI) x_(3)+=2.*M_PI;
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;

std::cout << " Complete update lidar " << std::endl;
  
   
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

  std::cout << " In Update Radar" << std::endl;
  int n_z = 3;

  // define spreading parameter
  double lambda = 3 - n_aug_;

  // set vector for weights
  VectorXd weights = VectorXd(2*n_aug_+1);
  double weight_0 = lambda/(lambda+n_aug_);
  double weight = 0.5/(lambda+n_aug_);
  weights(0) = weight_0;

  for (int i=1; i<2*n_aug_+1; ++i) {  
    weights(i) = weight;
  }

  


// create matrix for sigma points in measurement space
  Eigen::MatrixXd Zsig = Eigen::MatrixXd(n_z, 2 * n_aug_ + 1);

  // mean predicted measurement
  Eigen::VectorXd z_pred = VectorXd(n_z);
  
  // measurement covariance matrix S
  Eigen::MatrixXd S = MatrixXd(n_z,n_z);
  Eigen::VectorXd z = meas_package.raw_measurements_;
 
  // transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
    // extract values for better readability
    // Xsig_pred is comming from generate sigma points and then
    // predict sigma points in the prediction step.
   
    double p_x = Xsig_pred(0,i);
    double p_y = Xsig_pred(1,i);
    double v  = Xsig_pred(2,i);
    double yaw = Xsig_pred(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                       // r
    Zsig(1,i) = atan2(p_y,p_x);                                // phi
    Zsig(2,i) = (p_x*v1 + p_y*v2) / std::max(0.0001,sqrt(p_x*p_x + p_y*p_y));   // r_dot
  }

  // mean predicted measurement
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; ++i) {
    z_pred = z_pred + weights(i) * Zsig.col(i);
  }

  // innovation covariance matrix S
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
    // residual
    Eigen::VectorXd z_diff = Zsig.col(i) - z_pred;

    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  
    S = S + weights(i) * z_diff * z_diff.transpose();
  }

  // add measurement noise covariance matrix
  Eigen::MatrixXd R = Eigen::MatrixXd(n_z,n_z);
  R <<  std_radr_*std_radr_, 0, 0,
        0, std_radphi_*std_radphi_, 0,
        0, 0,std_radrd_*std_radrd_;
  S = S + R;

  
 Eigen::MatrixXd Tc = Eigen::MatrixXd(n_x_, n_z);
 for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
    // residual
    Eigen::VectorXd z_diff = Zsig.col(i) - z_pred;
   
    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    Eigen::VectorXd x_diff = Xsig_pred.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights(i) * x_diff * z_diff.transpose();
    
  }

  // Kalman gain K;
  Eigen::MatrixXd K = Tc * S.inverse();

  // residual
  Eigen::VectorXd z_diff = z - z_pred;

  // angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  // update state mean and covariance matrix
  x_ = x_ + K * z_diff;



  P_ = P_ - K*S*K.transpose();

  // print result
 

}
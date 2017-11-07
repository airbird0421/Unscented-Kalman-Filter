#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 15, 0, 0, 
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;
  // Process noise standard deviation longitudinal acceleration in m/s^2

  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2

  std_yawdd_ = 0.4;

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_ = false;

  n_x_ = 5;
 
  n_aug_ = 7;

  lambda_ = 3 - n_aug_;

  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < 2 * n_aug_ + 1; i++) {
    weights_(i) = 0.5 / (n_aug_ + lambda_);
  } 

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  gt0_352 = gt7_815 = gt5_991 = gt0_103 = 0;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  if (!is_initialized_) {
    if(meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_ << meas_package.raw_measurements_(0),
            meas_package.raw_measurements_(1), 0, 0, 0;
    } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);
      x_ << rho * cos(phi), rho * sin(phi), 0, 0, 0;
    }

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;  
    return;
  }

  // call prediction function
  Prediction((meas_package.timestamp_ - time_us_ ) / 1000000.0);
  time_us_ = meas_package.timestamp_;

  // call update function
  if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } 
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  // generate sigma points
  VectorXd x(n_aug_);
  x << x_, 0, 0;
  
  MatrixXd Q(2, 2);
  Q << std_a_ * std_a_, 0, 0, std_yawdd_ * std_yawdd_;
  
  MatrixXd P_aug(n_aug_, n_aug_);
  P_aug << P_, MatrixXd::Zero(n_x_, n_aug_ - n_x_), 
           MatrixXd::Zero(n_aug_ - n_x_, n_x_), Q;
  
  MatrixXd A = P_aug.llt().matrixL();
  MatrixXd Xsig(n_aug_, 2 * n_aug_ + 1);

  Xsig << x, x.replicate(1, n_aug_) + sqrt(lambda_ + n_aug_) * A,
          x.replicate(1, n_aug_) - sqrt(lambda_ + n_aug_) * A;

  // predict
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    double px = Xsig(0, i), py = Xsig(1, i), v = Xsig(2, i);
    double yaw = Xsig(3, i), yawd = Xsig(4, i), a = Xsig(5, i);
    double yawdd = Xsig(6, i);

    double px_p, py_p;

    if (fabs(yawd) > 0.001) {
      px_p = px + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
      py_p = py + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    } else {
      px_p = px + v * delta_t * cos(yaw);
      py_p = py + v * delta_t * sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    // add noise
    px_p += 0.5 * a * delta_t * delta_t;
    py_p += 0.5 * a * delta_t * delta_t;
    v_p += a * delta_t;

    yaw_p += 0.5 * yawdd * delta_t * delta_t;
    yawd_p += yawdd * delta_t;

    Xsig_pred_.col(i) << px_p, py_p, v_p, yaw_p, yawd_p;
  }
  
  // update state and covariance
  //mean
  x_.setZero();
  for(int i = 0; i < 2 * n_aug_ + 1; i++) {
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  //Covariance
  P_.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd t = Xsig_pred_.col(i) - x_;

    while(t(3) > M_PI) t(3) -= 2 * M_PI;
    while(t(3) < -M_PI) t(3) += 2 * M_PI;
    
    P_ += weights_(i) * t * t.transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  MatrixXd H(2, n_x_);
  MatrixXd R(2, 2);
  R << std_laspx_ * std_laspx_, 0, 0, std_laspy_ * std_laspy_;
  H << 1, 0, 0, 0, 0,
       0, 1, 0, 0, 0;

  VectorXd y = meas_package.raw_measurements_ - H * x_;
  MatrixXd S = H * P_ * H.transpose() + R;
  MatrixXd K = P_ * H.transpose() * S.inverse();

  x_ = x_ + K * y;
  int x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H) * P_;

  //NIS
  double nis = y.transpose() * S.inverse() * y;
  if (nis > 0.103) gt0_103++;
  if (nis > 5.991) gt5_991++;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  // transform to measurement space
  MatrixXd Zsig(3, 2 * n_aug_ + 1);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    double px = Xsig_pred_(0, i), py = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i), yaw = Xsig_pred_(3, i);

    double rho = sqrt(px * px + py * py);
    if (rho < 0.001) {
      cout << "devide by zero when calculating rho\n";
      return;
    }
    Zsig(0, i) = rho;
    Zsig(1, i) = atan2(py, px);
    Zsig(2, i) = (px * v * cos(yaw) + py * v * sin(yaw)) / rho;
  }

  //mean
  VectorXd z_pred(3);
  z_pred.setZero();
  for(int i = 0; i < 2 * n_aug_ + 1; i++) {
    z_pred += weights_(i) * Zsig.col(i);
  }

  //Covariance
  MatrixXd S(3, 3), R(3, 3);
  S.setZero();
  R << std_radr_ * std_radr_, 0, 0,
       0, std_radphi_ * std_radphi_, 0,
       0, 0, std_radrd_ * std_radrd_;

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd t = Zsig.col(i) - z_pred;

    while(t(1) > M_PI) t(1) -= 2 * M_PI;
    while(t(1) < -M_PI) t(1) += 2 * M_PI;

    S += weights_(i) * t * t.transpose();
  }

  S += R;

  //cross correlation matrix Tc
  MatrixXd Tc(n_x_, 3);

  Tc.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd tz = Zsig.col(i) - z_pred;

    while(tz(1) > M_PI) tz(1) -= 2 * M_PI;
    while(tz(1) < -M_PI) tz(1) += 2 * M_PI;

    VectorXd tx = Xsig_pred_.col(i) - x_;

    while(tx(3) > M_PI) tx(3) -= 2 * M_PI;
    while(tx(3) < -M_PI) tx(3) += 2 * M_PI;

    Tc += weights_(i) * tx * tz.transpose();
  }

  VectorXd tm = meas_package.raw_measurements_ - z_pred;
  while(tm(1) > M_PI) tm(1) -= 2 * M_PI;
  while(tm(1) < -M_PI) tm(1) += 2 * M_PI;

  MatrixXd K = Tc * S.inverse();

  //update x and P
  x_ += K * tm;
  P_ -= K * S * K.transpose();

  //NIS
  double nis = tm.transpose() * S.inverse() * tm;
  if (nis > 0.352) gt0_352++;
  if (nis > 7.815) gt7_815++;
}

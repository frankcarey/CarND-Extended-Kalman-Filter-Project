#include "kalman_filter.h"
#include <iostream>
//#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {

  // Get predicted x_ values.
  double px = x_(0);
  double py = x_(1);
  double vx = x_(2);
  double vy = x_(3);

  double px2_py2 = px*px + py*py;

  // Check division by zero
  if (fabs(px2_py2) < 0.0001) {
    std::cout << "Err: Division by Zero";
    return;
  }

  // Convert old x_ values to polar coordinates to match radar measurement.
  double rho = sqrt(px2_py2);
  double phi = atan2(py,px);
  double rho_dot = (px*vx + py*vy)/rho;

  // This is called h(x) in the lessons.
  // We use this instead of H_ in EKF.
  VectorXd hx_ = VectorXd(3);
  hx_ << rho, phi, rho_dot;
  VectorXd y = z - hx_;

  y(1) = atan2(sin(y(1)), cos(y(1)));
  // TODO: Is it better to use this approach of addition and subtraction instead of normalizing like above?
  //  if (y(1) < -M_PI) {
  //    std::cout << "Update() -- PHI < -3.14!!!!  PHI= " << y(1);
  //    y(1) = y(1) + 2*M_PI;
  //    std::cout << " .....Now it's: " << y(1) << std::endl;
  //  }
  //  else if (y(1) > M_PI) {
  //    std::cout << "PHI > 3.14!!!!  PHI= " << y(1);
  //    y(1) = y(1) - 2*M_PI;
  //    std::cout << " .....Now it's: " << y(1) << std::endl;
  //  }

  // Use Hj (Jacobian) to calculate S, K, P
  // H_ is set to Hj_ in FusionEKF::ProcessMeasurement() before calling this function
  // R_ is set to R_radar_ in FusionEKF::ProcessMeasurement() before calling this function
  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  MatrixXd K = P_ * H_.transpose() * S.inverse();

  // Update state x and new P with a new estimate
  x_ = x_ + (K * y);

  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

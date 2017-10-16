#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices

  // R_ is the measurement covariance matrix.
  // R_laser is the R_ used for laser measurements and doesn't change.
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << 0.0225, 0,
      0, 0.0225;

  // R_radar is the R_ used for radar measurements and doesn't change.
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << 0.09, 0, 0,
      0, 0.0009, 0,
      0, 0, 0.09;


  // H_ is the measurement matrix.
  // H_laser is the H_ used for laser measurements and doesn't change.
  H_laser_ = MatrixXd(2, 4);
  H_laser_<< 1,0,0,0,
      0,1,0,0;

  // Hj_ is calculated via Tools::CalculateJacobian() on every
  // RADAR UPDATE in ProcessMeasurement below so no further
  // initialization is needed.
  Hj_ = MatrixXd(3, 4);

  // Kalman Filter is where the Predict(), Update(), and UpdateEKF() steps happen.
  // TODO: Is this step needed?
  ekf_ = KalmanFilter();

  //Create a 4D state vector, we don't know yet the values of the x state
  // TODO: Is this step needed?
  ekf_.x_ = VectorXd(4);

  // P_ is the State Covariance matrix.
  // This is updated on every Predict() and Update() step.
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1000, 0,
      0, 0, 0, 1000;

  // F_ is the Transition Matrix.
  // (0,2) and (1,3) get updated during the prediction step and are set to delta_t.
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
      0, 1, 0, 1,
      0, 0, 1, 0,
      0, 0, 0, 1;

  // Q_ is the Process Uncertainty matrix.
  // It is updated on every prediction step based on the acceleration noise and time delta.
  // We initialize to zero to start with.
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << 0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0;

  //NOTE: ekf_.H_ and ekf_.R_ get replaced with either laser or radar versions on each measurement,
  // so we don't need to initialize them here.

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates and initialize state.
      double rho = measurement_pack.raw_measurements_[0];
      double theta = measurement_pack.raw_measurements_[1];
      // Can we use ro-hat to set the predicted velocity? Per "tips" we don't get enough data for that,
      // so set x and y, but leave vx and vy at zero.
      ekf_.x_ << rho * cos(theta), rho * sin(theta), 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // Initialize state but laser only gives position data, so leave vx and vy at zero.
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */


  // Compute the time elapsed between the current and previous measurements.
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  double dt_2 = dt * dt;
  double dt_3 = dt_2 * dt;
  double dt_4 = dt_3 * dt;

  // Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  // Set the acceleration noise components.
  // TODO: Where do these come from?
  double noise_ax = 9;
  double noise_ay = 9;

  // Update the process noise covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
          0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
          dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
          0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  // RADAR UPDATE
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    if (Hj_.isZero()) {
      // If there was some error calculating the Jacobian, it will return all zeros.
      // If that happens, then just skip the update step.
      return;
    }
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    // Radar updates need to use Extended KF Updates.
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  // LASER UPDATE
  } else {
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    // Laser update can use regular KF Updates.
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}

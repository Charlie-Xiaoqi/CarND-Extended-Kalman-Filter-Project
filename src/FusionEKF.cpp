#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>
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
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
	//create a 4D state vector, we don't know yet the values of the x state
  H_laser_ << 1, 0, 0, 0,
          0, 1, 0, 0;


  ekf_.x_ = VectorXd(4);

    //state covariance matrix P
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
              0, 1, 0, 0,
              0, 0, 10, 0,
              0, 0, 0, 10;

    //the initial transition matrix F_
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
              0, 1, 0, 1,
              0, 0, 1, 0,
              0, 0, 0, 1;

  ekf_.Q_ = MatrixXd(4, 4);
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
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */

     float meas_rho = measurement_pack.raw_measurements_[0];
     float meas_phi = measurement_pack.raw_measurements_[1];
     if (meas_phi> 3.14159265358979323846 || meas_phi < -3.14159265358979323846 ) {
        while (meas_phi> 3.14159265358979323846){
            meas_phi = meas_phi - 2*3.14159265358979323846;
        }
        while (meas_phi< -3.14159265358979323846){
            meas_phi = meas_phi + 2*3.14159265358979323846;
        }
     }

     float meas_rho_dot = measurement_pack.raw_measurements_[2];

     float  px = meas_rho * cos(meas_phi);
     float  py = - meas_rho * sin(meas_phi);
     ekf_.x_ << px, py, 0, 0;
     /*
     Hj_ = tools.CalculateJacobian( ekf_.x_ );
     ekf_.Init(ekf_.x_, ekf_.P_, ekf_.F_,
                        H_laser_, R_radar_, ekf_.Q_)
     */
	 previous_timestamp_ = measurement_pack.timestamp_;


    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
     ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
     /*
     ekf_.Init(ekf_.x_, ekf_.P_, ekf_.F_,
                        H_laser_, R_laser_, ekf_.Q_)
	 */
	 previous_timestamp_ = measurement_pack.timestamp_;

    }

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
    //compute the time elapsed between the current and previous measurements
	float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
	previous_timestamp_ = measurement_pack.timestamp_;

    float dt_2 = dt * dt;
    float dt_3 = dt_2 * dt;
    float dt_4 = dt_3 * dt;
    ekf_.F_(0,2) = dt;
    ekf_.F_(1,3) = dt;
    float noise_ax = 9;
    float noise_ay = 9;
    ekf_.Q_ << dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
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

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    Hj_ = tools.CalculateJacobian( ekf_.x_ );
    ekf_.Init(ekf_.x_, ekf_.P_, ekf_.F_,
                        Hj_, R_radar_, ekf_.Q_);
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);


  } else {
    ekf_.Init(ekf_.x_, ekf_.P_, ekf_.F_,
                H_laser_, R_laser_, ekf_.Q_);
    // Laser updates
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}

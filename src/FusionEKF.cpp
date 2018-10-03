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

  H_laser_ << 1, 0, 0, 0,
			  0, 1, 0, 0;

  Hj_ << 1, 0, 0, 0,
			  0, 1, 0, 0,
        0, 0, 0, 0;

  //ekf_.Init(x_in, P_in, F_in, H_in, R_in, Q_in) 


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
      * 
    */

    // first measurement
    cout << "EKF: " << endl;
    VectorXd x = VectorXd(4);
    //x << 1, 1, 1, 1;

    MatrixXd H; MatrixXd R;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho = measurement_pack.raw_measurements_[0];
      float theta = measurement_pack.raw_measurements_[1];
      float x1 = rho* sin(theta);
      float x2 = rho* cos(theta);
      x << x1,x2, 0, 0;
      MatrixXd H = tools.CalculateJacobian(x);
      MatrixXd R = R_radar_;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
      MatrixXd H = H_laser_;
      MatrixXd R = R_laser_;
    } 

       	//state covariance matrix P
	  MatrixXd P = MatrixXd(4, 4);
	  P << 1, 0, 0, 0,
			  0, 1, 0, 0,
			  0, 0, 1000, 0,
			  0, 0, 0, 1000;

    MatrixXd Q = MatrixXd(4, 4);
	  Q << 1e-3, 0, 0, 0,
			  0, 1e-3, 0, 0,
			  0, 0, 1e-3, 0,
			  0, 0, 0, 1e-3;

    MatrixXd F = MatrixXd(4, 4);
    F << 1, 0, 1, 0,
			  0, 1, 0, 1,
			  0, 0, 1, 0,
			  0, 0, 0, 1;


    ekf_.Init(x, P, F, H, R, Q);

    // done initializing, no need to predict or update
    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    // cout << "finished initialisation of FusionEKF " << endl;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  //cout << "starting Prediction " << endl;
  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds. */
    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
	  previous_timestamp_ = measurement_pack.timestamp_;

    /* Update the process noise covariance matrix. */
    float dt_2 = dt * dt;
	  float dt_3 = dt_2 * dt;
	  float dt_4 = dt_3 * dt;

	  //Modify the F matrix so that the time is integrated
	  ekf_.F_(0, 2) = dt;
	  ekf_.F_(1, 3) = dt;

    //cout << ekf_.F_ << endl;

    /* Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
    */
    float noise_ax = 9;
	  float noise_ay = 9;

	  //set the process covariance matrix Q
	  ekf_.Q_ = MatrixXd(4, 4);
	  ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
			   0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
			   dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
			   0, dt_3/2*noise_ay, 0, dt_2*noise_ay;
     

  ekf_.Predict();
  //cout << "finished Prediction, starting Update " << endl;

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */
  //cout << "meas type" << measurement_pack.sensor_type_ << endl;
  
    VectorXd z_3d(3);
    z_3d << 1.,1.,1.;

    VectorXd z_2d(2);
    z_2d << 1.,1.;

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    //cout << " its a radar type measurement" << endl;
    // Radar updates
    //cout << measurement_pack.raw_measurements_[0] << " (rho) " << measurement_pack.raw_measurements_[1] << " (phi) " << measurement_pack.raw_measurements_[2] << " (rhodot)" << endl;
 
    //cout << "still alive, size of z_3d: " << z_3d.size() << " size of x: " << ekf_.x_.size() << endl;
    z_3d << measurement_pack.raw_measurements_[0],measurement_pack.raw_measurements_[1],measurement_pack.raw_measurements_[2];
    //cout << "just before Jacobian in Radar Update " << endl;
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    //cout << "finished Jacobian in Radar Update " << endl;
    ekf_.R_ = R_radar_;
    //cout << "just before UpdateEKF in Radar Update " << endl;
    ekf_.UpdateEKF(z_3d);
    //cout << "just after UpdateEKF in Radar Update " << endl;
  } else {
    //cout << " its a laser type measurement" << endl;
    // Laser updates
      
      z_2d << measurement_pack.raw_measurements_[0],measurement_pack.raw_measurements_[1];
      ekf_.H_ = H_laser_;
      ekf_.R_ = R_laser_;
      //cout << "just before Laser Update " << endl;
    	ekf_.Update(z_2d);
  }

  // print the output
  //cout << "x_ = " << ekf_.x_ << endl;
  //cout << "P_ = " << ekf_.P_ << endl;
}

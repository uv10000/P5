#include "kalman_filter.h"
#include <iostream>

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
  /**
  TODO:
    * predict the state
  */
    x_ = F_ * x_ ; // + u;
		MatrixXd Ft = F_.transpose();
		P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
    VectorXd y = z - H_ * x_;
		MatrixXd Ht = H_.transpose();
		MatrixXd S = H_ * P_ * Ht + R_;
		MatrixXd Si = S.inverse();
		MatrixXd K =  P_ * Ht * Si;

    //create identity matrix 
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);

    //new state
		x_ = x_ + (K * y);
		P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
    

    float px = x_(0);
	  float py = x_(1);
	  float vx = x_(2);
	  float vy = x_(3);

	  float rho = sqrt(px*px+py*py);
    if (rho < 1e-4 || px*px < 1e-8 || py*py < 1e-8) {cout << rho << px << py << " (rho,px,py), at least one nearly zero" << endl;}
    float theta = atan2(py,px); // possible angle jumps will be dealt with below
    float rhodot = (px*vx + py*vy)/rho; 
    VectorXd zp(3); // nonlinear measurement
    zp << rho, theta, rhodot;

    VectorXd y = z - zp; //H_ * x_;
    if (y(1) < - M_PI  ) {cout << "angle lt PI" << endl; y(1) = y(1) + 2 * M_PI;}      //handle jumps originating from atan2
    else if(y(1) > M_PI  ) {cout << "angle gt PI" << endl; y(1) = y(1)  - 2 * M_PI;}   //handle jumps originating from atan2
		MatrixXd Ht = H_.transpose();
		MatrixXd S = H_ * P_ * Ht + R_;
		MatrixXd Si = S.inverse();
		MatrixXd K =  P_ * Ht * Si;

		//new state
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
		x_ = x_ + (K * y);
		P_ = (I - K * H_) * P_;
}

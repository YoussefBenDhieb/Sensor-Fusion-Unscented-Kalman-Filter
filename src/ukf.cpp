#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd::Zero(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_<< 1,0,0,0,0,
  	   0,1,0,0,0,
  	   0,0,10,0,0,
  	   0,0,0,10,0,
  	   0,0,0,0,10;
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  ///* initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;

  //set state dimension
  n_x_ = 5;

  //set augmented state dimension
  n_aug_ = 7;

  //define spreading parameter
  lambda_ = 3 - n_aug_;

  //create vector for weights
  weights_ = VectorXd::Zero(2*n_aug_+1);

  //set weights
  weights_(0)= lambda_/ (lambda_ + n_aug_);
  for (int i=1; i<2*n_aug_+1; i++){
      weights_(i) = 0.5/(lambda_+n_aug_);
  }

  ///* predicted sigma points matrix
  Xsig_pred_ = MatrixXd::Zero(n_x_, 2 * n_aug_ + 1);

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_pack) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  if (!is_initialized_) {
	if (meas_pack.sensor_type_ == MeasurementPackage::RADAR) {
	  /**
	  Convert radar from polar to cartesian coordinates and initialize state.
	  */
		//set the state with the initial location and zero velocity
		x_ << meas_pack.raw_measurements_[0]*sin(meas_pack.raw_measurements_[1]), meas_pack.raw_measurements_[0]*cos(meas_pack.raw_measurements_[1]),
					 meas_pack.raw_measurements_[2]*sin(meas_pack.raw_measurements_[1]), meas_pack.raw_measurements_[2]*cos(meas_pack.raw_measurements_[1]), 1;

		time_us_ = meas_pack.timestamp_;

	}
	else if (meas_pack.sensor_type_ == MeasurementPackage::LASER) {
	  /**
	  Initialize state.
	  */
		//set the state with the initial location and zero velocity
		x_ << meas_pack.raw_measurements_[0], meas_pack.raw_measurements_[1], 1, 1, 1;

		time_us_ = meas_pack.timestamp_;
	}
	// done initializing, no need to predict or update
	is_initialized_ = true;
  }


  //compute the time elapsed between the current and previous measurements
  double dt = (meas_pack.timestamp_ - time_us_) / 1000000.0;		//dt - expressed in seconds
  time_us_= meas_pack.timestamp_;
  //Prediction step
  Prediction(dt);

  // Update step
  if (meas_pack.sensor_type_ == MeasurementPackage::LASER and use_laser_) {
  	UpdateLidar(meas_pack);

	}
  else if (meas_pack.sensor_type_ == MeasurementPackage::RADAR and use_radar_) {
	UpdateRadar(meas_pack);

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
  //create augmented mean vector
  VectorXd x_aug = VectorXd::Zero(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd::Zero(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd::Zero(n_aug_, 2 * n_aug_ + 1);
  
  //create augmented mean state
  x_aug.head(n_x_) = x_;
  
  //create augmented covariance matrix
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_ *std_yawdd_;
  
  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();
  
  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for (int i=0; i<n_aug_; i++){
      Xsig_aug.col(i+1) = x_aug + sqrt(lambda_+n_aug_)*A.col(i);
      Xsig_aug.col(i+n_aug_+1) = x_aug - sqrt(lambda_+n_aug_)*A.col(i);
  }

  //predict sigma points
  for (int i=0; i<2*n_aug_+1; i++){
      //avoid division by zero
      if (Xsig_aug(4,i) != 0){
          Xsig_pred_(0,i) = Xsig_aug(0,i)+(Xsig_aug(2,i)/Xsig_aug(4,i))*(sin(Xsig_aug(3,i)+Xsig_aug(4,i)*delta_t)-sin(Xsig_aug(3,i)))+0.5*delta_t*delta_t*cos(Xsig_aug(3,i))*Xsig_aug(5,i);
          Xsig_pred_(1,i) = Xsig_aug(1,i)+(Xsig_aug(2,i)/Xsig_aug(4,i))*(-cos(Xsig_aug(3,i)+Xsig_aug(4,i)*delta_t)+cos(Xsig_aug(3,i)))+0.5*delta_t*delta_t*sin(Xsig_aug(3,i))*Xsig_aug(5,i);
      }else{
          Xsig_pred_(0,i) = Xsig_aug(0,i)+Xsig_aug(2,i)*cos(Xsig_aug(3,i))*delta_t+0.5*delta_t*delta_t*cos(Xsig_aug(3,i))*Xsig_aug(5,i);
          Xsig_pred_(1,i) = Xsig_aug(1,i)+Xsig_aug(2,i)*sin(Xsig_aug(3,i))*delta_t+0.5*delta_t*delta_t*sin(Xsig_aug(3,i))*Xsig_aug(5,i);
      }
      //write predicted sigma points into right column
      Xsig_pred_(2,i) = Xsig_aug(2,i) + delta_t * Xsig_aug(5,i);
      Xsig_pred_(3,i) = Xsig_aug(3,i) + Xsig_aug(4,i) * delta_t + 0.5 * delta_t * delta_t * Xsig_aug(6,i);
      Xsig_pred_(4,i) = Xsig_aug(4,i) + delta_t * Xsig_aug(6,i);
  }
  
  //predict state mean
  x_.fill(0.0);
  for (int i = 0; i<2*n_aug_+1; i++){
      x_ += Xsig_pred_.col(i)*weights_(i);
  }
  
  //predict state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i<2*n_aug_+1; i++){
      VectorXd N = Xsig_pred_.col(i)-x_;
      while(N(3)> M_PI) N(3) -= 2.*M_PI;
      while(N(3)< -M_PI) N(3) += 2.*M_PI;
      P_ += (N*N.transpose())*weights_(i);
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
  //set measurement dimension, lidar can measure px, and py
  int n_z = 2;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd::Zero(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd::Zero(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd::Zero(n_z,n_z);

  //measurement vector
  VectorXd z = VectorXd::Zero(n_z);
  z = meas_package.raw_measurements_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);

  //transform sigma points into measurement space
  Zsig = Xsig_pred_.topRows(n_z); 
  
  //calculate mean predicted measurement
  for (int i=0; i<2*n_aug_+1; i++){
      z_pred += Zsig.col(i)*weights_(i); 
  }
  //calculate innovation covariance matrix S
  S.fill(0.0);
  for (int i=0; i<2*n_aug_+1; i++){
      MatrixXd N = Zsig.col(i)-z_pred;
      S += N*N.transpose()*weights_(i);
  }
  MatrixXd R = MatrixXd::Zero(n_z,n_z);
  R << std_laspx_*std_laspx_, 0,
       0, std_laspy_*std_laspy_;
  S += R;
  
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i=0; i< 2*n_aug_+1; i++){
      MatrixXd N = (Zsig.col(i)-z_pred);
      Tc += weights_(i)*(Xsig_pred_.col(i)-x_)*N.transpose();
  }
  //calculate Kalman gain K;
  MatrixXd K = Tc*S.inverse();
  
  //update state mean and covariance matrix
  x_ += K*(z-z_pred);
  P_ += -K*S*K.transpose();

  //Calculate NIS
  NIS = (z-z_pred)*S.inverse()*(z-z_pred);
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
  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd::Zero(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd::Zero(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd::Zero(n_z,n_z);

  //measurement vector
  VectorXd z = VectorXd::Zero(n_z);
  z = meas_package.raw_measurements_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);

  //transform sigma points into measurement space
  for(int i=0; i<2*n_aug_+1; i++){
      //Calculate rho
      Zsig(0,i) = sqrt(Xsig_pred_(0,i)*Xsig_pred_(0,i)+Xsig_pred_(1,i)*Xsig_pred_(1,i)); 
      //Calculate phi
      Zsig(1,i) = atan2(Xsig_pred_(1,i),Xsig_pred_(0,i)); 

      //Calculate rho_dot
      Zsig(2,i) = (Xsig_pred_(0,i)*Xsig_pred_(2,i)*cos(Xsig_pred_(3,i))+Xsig_pred_(1,i)*Xsig_pred_(2,i)*sin(Xsig_pred_(3,i)))/Zsig(0,i);
  }
  //calculate mean predicted measurement
  for (int i=0; i<2*n_aug_+1; i++){
      z_pred += Zsig.col(i)*weights_(i); 
  }

  //calculate innovation covariance matrix S
  for (int i=0; i<2*n_aug_+1; i++){
      MatrixXd N = Zsig.col(i)-z_pred;
      //Normalize phi
      while (N(1)>M_PI) N(1) -= 2*M_PI;
      while (N(1)<-M_PI) N(1) += 2*M_PI;
      S += N*N.transpose()*weights_(i);
  }
  cout<<S<<"\n";
  MatrixXd R = MatrixXd::Zero(n_z,n_z);
  R << std_radr_*std_radr_, 0, 0,
       0, std_radphi_*std_radphi_, 0,
       0, 0, std_radrd_*std_radrd_;
  S += R;
 
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i=0; i< 2*n_aug_+1; i++){
    MatrixXd z_diff = (Zsig.col(i)-z_pred);
    //Normalize angle
    while (z_diff(1)>M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
    Tc += weights_(i)*x_diff*z_diff.transpose();
  }
  //calculate Kalman gain K;
  MatrixXd K = Tc*S.inverse();
  
  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  //Calculate NIS
  NIS = (z-z_pred)*S.inverse()*(z-z_pred);

}

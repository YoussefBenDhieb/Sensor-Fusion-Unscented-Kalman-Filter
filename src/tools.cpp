#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
 	/**
  TODO:
    * Calculate the RMSE here.
  */
	VectorXd rmse(4);
	rmse << 0,0,0,0;
    VectorXd temp_diff(4);
    
    // check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
		if(estimations.size()==0)
	    {
	        cout<<"the estimation vector's size is equal to 0!";
	        return rmse;
	    }
	    else if (estimations.size()!=ground_truth.size())
	    {
	        cout<<"the estimation vector's size is not equal to the ground truth vector's size!";
	        return rmse;
	    }

	//accumulate squared residuals
	
	for(int i=0; i < estimations.size(); ++i){
       
        temp_diff = estimations[i]-ground_truth[i];
		rmse = rmse.array() + temp_diff.array()*temp_diff.array();
	}

	//calculate the mean
    rmse /= estimations.size();
	//calculate the squared root
    rmse = rmse.array().sqrt();
	//return the result
	return rmse;
}
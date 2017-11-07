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
  VectorXd rmse(4), diff(4);
  int size = estimations.size();
  rmse << 0.0, 0.0, 0.0, 0.0;

  if(size != ground_truth.size() || size == 0) {
    cout << "size error to calculate RMSE!";
    return rmse;
  }

  for(int i = 0; i < size; i++) {
    diff = estimations[i] - ground_truth[i];
    rmse += (VectorXd)(diff.array() * diff.array());
  }

  rmse /= size;
  rmse = rmse.array().sqrt();
  cout << "rmse" << size << "\n" << rmse << "\n";
  
  return rmse;

}

void Tools::NisStats(int size, int gt0_103, int gt5_991, int gt0_352, int gt7_815)
{
  // NIS statistics
  cout << "lidar > 0.103:" << (float)gt0_103 / size << "\n";
  cout << "lidar > 5.991:" << (float)gt5_991 / size << "\n";
  cout << "radar > 0.352:" << (float)gt0_352 / size << "\n";
  cout << "radar > 7.815:" << (float)gt7_815 / size << "\n";
}

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
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

    MatrixXd Hj(3,4);
    //recover state parameters
    double px = x_state(0);
    double py = x_state(1);
    double vx = x_state(2);
    double vy = x_state(3);

    //compute the Jacobian matrix
    double px2_py2 = pow(px, 2) + pow(py, 2);
    double px2_py2_32 = sqrt(pow(px2_py2, 3));
    double vx_py = vx * py;
    double vy_px = vy * px;

    //check division by zero-ish
    if (fabs(px2_py2) < 0.0001) {
        cout << "TODO: Division by Zero";
    }

    // Implement all the derivatives into the matrix.
    Hj << px / sqrt(px2_py2),                py / sqrt(px2_py2),                0,                  0,
          - py / px2_py2,                    px / px2_py2,                      0,                  0,
          py * (vx_py - vy_px) / px2_py2_32, py * (vy_px - vx_py) / px2_py2_32, px / sqrt(px2_py2), py / sqrt(px2_py2);

    return Hj;
}

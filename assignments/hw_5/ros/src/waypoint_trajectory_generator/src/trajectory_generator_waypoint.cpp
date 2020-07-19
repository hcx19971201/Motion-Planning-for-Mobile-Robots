#include "trajectory_generator_waypoint.h"
#include <stdio.h>
#include <ros/ros.h>
#include <ros/console.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;    
using namespace Eigen;

TrajectoryGeneratorWaypoint::TrajectoryGeneratorWaypoint(){}
TrajectoryGeneratorWaypoint::~TrajectoryGeneratorWaypoint(){}

//define factorial function, input i, output i!
int TrajectoryGeneratorWaypoint::Factorial(int x)
{
    int fac = 1;
    for(int i = x; i > 0; i--)
        fac = fac * i;
    return fac;
}
/*

    STEP 2: Learn the "Closed-form solution to minimum snap" in L5, then finish this PolyQPGeneration function

    variable declaration: input       const int d_order,                    // the order of derivative
                                      const Eigen::MatrixXd &Path,          // waypoints coordinates (3d)
                                      const Eigen::MatrixXd &Vel,           // boundary velocity
                                      const Eigen::MatrixXd &Acc,           // boundary acceleration
                                      const Eigen::VectorXd &Time)          // time allocation in each segment
                          output      MatrixXd PolyCoeff(m, 3 * p_num1d);   // position(x,y,z), so we need (3 * p_num1d) coefficients

*/

Eigen::MatrixXd TrajectoryGeneratorWaypoint::PolyQPGeneration(
            const int d_order,                    // the order of derivative
            const Eigen::MatrixXd &Path,          // waypoints coordinates (3d)
            const Eigen::MatrixXd &Vel,           // boundary velocity
            const Eigen::MatrixXd &Acc,           // boundary acceleration
            const Eigen::VectorXd &Time)          // time allocation in each segment
{
    // enforce initial and final velocity and accleration, for higher order derivatives, just assume them be 0;
    int p_order   = 2 * d_order - 1;              // the order of polynomial
    int p_num1d   = p_order + 1;                  // the number of variables in each segment

    int m = Time.size();                          // the number of segments
    MatrixXd PolyCoeff = MatrixXd::Zero(m, 3 * p_num1d);           // position(x,y,z), so we need (3 * p_num1d) coefficients
    VectorXd Px(p_num1d * m), Py(p_num1d * m), Pz(p_num1d * m);

    /*   Produce Mapping Matrix A to the entire trajectory, A is a mapping matrix that maps polynomial coefficients to derivatives.   */    
    // Q matrix
    MatrixXd Q = MatrixXd::Zero(m * p_num1d, m * p_num1d);

    for (int k = 0; k < m; k++)
    {
        MatrixXd Qk = MatrixXd::Zero(p_num1d, p_num1d);
        
        double T = Time(k);
        for (int i = d_order; i < p_num1d; i++)
        {
            for (int j = d_order; j < p_num1d; j++)
            {
                Qk(i, j) = Factorial(i) / Factorial(i - d_order) * Factorial(j) / Factorial(j - d_order) / (i + j - p_order) * pow(T, i + j - p_order);
            }
        }

        Q.block(k*p_num1d, k*p_num1d, p_num1d, p_num1d) = Qk;
    }

    // M matrix
    MatrixXd M = MatrixXd::Zero(m * p_num1d, m * p_num1d);
    
    for (int k = 0; k < m; k++)
    {
        MatrixXd Mk = MatrixXd::Zero(p_num1d, p_num1d);

        for (int i = 0; i < d_order; i++)
        {
            Mk(i, i) = Factorial(i);
        }
        
        double T = Time(k);
        for (int i = d_order; i < p_num1d; i++)
        {
            for (int j = i - d_order; j < p_num1d; j++)
            {
                Mk(i, j) = Factorial(j) / Factorial(j-(i-d_order)) * pow(T, j-(i-d_order));
            }
        }

        M.block(k*p_num1d, k*p_num1d, p_num1d, p_num1d) = Mk;
    }

    // C Matrix
    MatrixXd Ct = MatrixXd::Zero(m * p_num1d, (m+1) * d_order);

    /// constrained variables
    MatrixXd C0 = MatrixXd::Identity(d_order, d_order);
    MatrixXd CK = MatrixXd::Identity(d_order, d_order);

    Ct.block(0, 0, d_order, d_order) = C0;

    VectorXd fblock = VectorXd::Zero(2*d_order);
    fblock(0) = 1;
    fblock(d_order) = 1;

    for (int k = 0; k < m-1; k++)
    {
        Ct.block(d_order+2*d_order*k, d_order+k, 2*d_order, 1) = fblock;
    }

    Ct.block(d_order*(2*m-1), d_order+m-1, d_order, d_order) = CK;

    /// free variables
    MatrixXd pblock = MatrixXd::Zero(p_num1d, d_order-1);
    pblock.block(1, 0, d_order-1, d_order-1) = MatrixXd::Identity(d_order-1, d_order-1);
    pblock.block(d_order+1, 0, d_order-1, d_order-1) = MatrixXd::Identity(d_order-1, d_order-1);

    for (int k = 0; k < m-1; k++)
    {
        Ct.block(d_order+(p_num1d)*k, (d_order-1)*k+2*d_order+(m-1)*1, p_num1d, d_order-1) = pblock;
    }

    /*   Produce the dereivatives in X, Y and Z axis directly.  */
    MatrixXd M_inv = M.inverse();
    MatrixXd R = Ct.transpose() * M_inv.transpose() * Q * M_inv * Ct;

    int num_dF = 2*d_order+(m-1)*1;
    int num_dP = (m-1)*(d_order-1);

    MatrixXd R_fp = R.topRightCorner(num_dF, num_dP);
    MatrixXd R_pp = R.bottomRightCorner(num_dP, num_dP);
    
    for (int i = 0; i < 3; i++)
    {
        VectorXd path = Path.col(i);

        VectorXd dF = VectorXd::Zero(num_dF);
        VectorXd dP = VectorXd::Zero(num_dP);

        dF(0) = path(0);
        dF(1) = Vel(0, i);
        dF(2) = Acc(0, i);

        for (int k = 0; k < m-1; k++)
        {
            dF(d_order+k) = path(k+1);
        }
        
        dF(d_order+m-1) = path(m);
        dF(d_order+m  ) = Vel(1, i);
        dF(d_order+m+1) = Acc(1, i);

        dP = -1.0 * R_pp.inverse() * R_fp.transpose() * dF;

        VectorXd d = VectorXd::Zero(num_dF+num_dP);
        d << dF, dP;

        VectorXd poly_coef = M_inv * Ct * d;
        MatrixXd poly_coef_row = poly_coef.transpose();

        for (int k = 0; k < m; k++)
        {
            PolyCoeff.block(k, i*p_num1d, 1, p_num1d) = poly_coef_row.block(0, k*p_num1d, 1, p_num1d);
        }
    }

    
    return PolyCoeff;
}

#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// TODO: Set the timestep length and duration
size_t N = 10;
double dt = 0.1;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

double ref_v = 60;
size_t x_start = 0;
size_t y_start = x_start + N;
size_t psi_start = y_start + N;
size_t v_start = psi_start + N;
size_t cte_start = v_start + N;
size_t epsi_start = cte_start + N;
size_t delta_start = epsi_start + N;
size_t a_start = delta_start + N -1;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.
    
    // The cost value is stored in the first element of 'fg'.
    // Any additions to the cost should be added to 'fg[0]'.
    fg[0] = 0;

    // Reference State Cost
    for (int i = 0;i < N; i++){
      fg[0] += 4000*CppAD::pow(vars[cte_start + i],2);
      fg[0] += 3500*CppAD::pow(vars[epsi_start + i],2);
      fg[0] += CppAD::pow(vars[v_start + i] - ref_v,2);
    }

    for (int i = 0;i < N-1;i++){
      //fg[0] += CppAD::pow(vars[delta_start + i],2);
      //fg[0] += CppAD::pow(vars[a_start + i],2);
      // add penalty for speed and steer
      //fg[0] += 700*CppAD::pow(vars[delta_start + i] * vars[v_start+i],2);
      fg[0] += 60000 * CppAD::pow(vars[delta_start + i],2);
      fg[0] += CppAD::pow(vars[a_start],2);
    }

    for (int i = 0;i < N-2;i++){
      fg[0] += 200*CppAD::pow(vars[delta_start + i + 1] - vars[delta_start + i],2);
      fg[0] += 10*CppAD::pow(vars[a_start + i + 1] - vars[a_start + i],2);
    }

    // Setup Constraints
    fg[1 + x_start] = vars[x_start];
    fg[1 + y_start] = vars[y_start];
    fg[1 + psi_start] = vars[psi_start];
    fg[1 + v_start] = vars[v_start];
    fg[1 + cte_start] = vars[cte_start];
    fg[1 + epsi_start] = vars[epsi_start];

    // For the rest of the constraints
    for (int i = 1;i < N; i++){
      AD<double> x_1 = vars[x_start + i];
      AD<double> x_0 = vars[x_start + i -1];
      AD<double> y_1 = vars[y_start + i];
      AD<double> y_0 = vars[y_start + i -1];
      AD<double> psi_1 = vars[psi_start + i];
      AD<double> psi_0 = vars[psi_start + i -1];
      AD<double> v_1 = vars[v_start + i];
      AD<double> v_0 = vars[v_start + i -1];
      AD<double> cte_1 = vars[cte_start + i];
      AD<double> cte_0 = vars[cte_start + i - 1];
      AD<double> epsi_1 = vars[epsi_start + i];
      AD<double> epsi_0 = vars[epsi_start + i -1];
      AD<double> a = vars[a_start + i -1];
      AD<double> delta = vars[delta_start + i -1];
      if(i > 1){ //use previous actuations
        a = vars[a_start + i -2];
        delta = vars[delta_start + i -2];
      }
      AD<double> f_0 = coeffs[0] + coeffs[1]*x_0 + coeffs[2] * CppAD::pow(x_0,2) + coeffs[3] * CppAD::pow(x_0,3);
      AD<double> psides_0 = CppAD::atan(coeffs[1] + 2 * coeffs[2] * x_0 + 3 * coeffs[3] * CppAD::pow(x_0,2));

      // Setup the rest of the model constraints
      fg[1 + x_start + i] = x_1 - (x_0 + v_0 * CppAD::cos(psi_0) * dt);
      fg[1 + y_start + i] = y_1 - (y_0 + v_0 * CppAD::sin(psi_0) * dt);
      fg[1 + psi_start + i] = psi_1 - (psi_0 - v_0/Lf * delta * dt);
      fg[1 + v_start + i] = v_1 - (v_0 + a * dt);
      fg[1 + cte_start + i] = cte_1 - ((f_0 - y_0) + (v_0 * CppAD::sin(epsi_0) * dt));
      fg[1 + epsi_start + i] = epsi_1 - ((psi_0 - psides_0) - v_0/Lf * delta * dt);
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double v = state[3];
  double cte = state[4];
  double epsi = state[5];

  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  size_t n_vars = 0;
  n_vars = N * 6 + (N -1) * 2;
  // TODO: Set the number of constraints
  size_t n_constraints = 0;
  n_constraints = N * 6;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // TODO: Set lower and upper limits for variables.
  
  // Set the initial variable values
  vars[x_start] = x;
  vars[y_start] = y;
  vars[psi_start] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[epsi_start] = epsi;

  // Set all non-actuators upper and lowerlimits
  // to the max negative and positive values.
  for (int i = 0; i < delta_start; i++){
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  for (int i = delta_start; i < a_start; i++){
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }

  // Acceleration/decelaration upper and lower limits.
  for (int i = a_start;i < n_vars;i++){
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[psi_start] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[epsi_start] = epsi;

  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[psi_start] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = epsi;

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  
  vector<double> result;

  result.push_back(solution.x[delta_start]);
  result.push_back(solution.x[a_start]);

  for (int i = 0; i < N; i++) {
    result.push_back(solution.x[x_start + i]);
    result.push_back(solution.x[y_start + i]);
  }

  return result;
}

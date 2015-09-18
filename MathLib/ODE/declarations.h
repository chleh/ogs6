#pragma once

namespace MathLib
{

// maybe use Eigen::Map here
// and use std::function
typedef void (* const Function)(const double t, double const*const y, double *const ydot);
}

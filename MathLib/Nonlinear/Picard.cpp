/*
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <cassert>

#include "Picard.h"

namespace
{

std::unique_ptr<MathLib::Nonlinear::Picard>
createPicard(const BaseLib::ConfigTreeNew &config)
{
    config.checkConfParam<std::string>("type", "Picard");

    std::unique_ptr<MathLib::Nonlinear::Picard> pic(new MathLib::Nonlinear::Picard);

    auto const abs_tol = config.getConfParamOptional<double>("absolute_tolerance");
    if (abs_tol) pic->setAbsTolerance(*abs_tol);

    auto const rel_tol = config.getConfParamOptional<double>("relative_tolerance");
    if (rel_tol) pic->setRelTolerance(*rel_tol);

    auto const max_iter = config.getConfParamOptional<std::size_t>("max_iterations");
    if (max_iter) pic->setMaxIterations(*max_iter);

    auto const print = config.getConfParamOptional<bool>("print_errors");
    if (print) pic->printErrors(*print);

    auto const norm = config.getConfParamOptional<std::string>("norm_type");
    if (norm) {
        if      (*norm == "Norm1")   pic->setNormType(MathLib::VecNormType::NORM1);
        else if (*norm == "Norm2")   pic->setNormType(MathLib::VecNormType::NORM2);
        else if (*norm == "NormMax") pic->setNormType(MathLib::VecNormType::INFINITY_N);
        else {
            ERR("invalid vector norm type: `%s'", norm->c_str());
            std::abort();
            return nullptr;
        }
    }

    return pic;
}

}

namespace MathLib
{
namespace Nonlinear
{
Picard::Picard()
    : _normType(VecNormType::INFINITY_N),
      _abs_tol(std::numeric_limits<double>::max()),
      _rel_tol(1e-6),
      _max_itr(25),
      _printErrors(false),
      _n_iterations(0),
      _abs_error(.0),
      _rel_error(.0)
{
}


std::unique_ptr<Picard>
createNonlinearSolver(BaseLib::ConfigTreeNew const& config)
{
    auto const type = config.getConfParam<std::string>("type");

    if (type == "Picard") {
        return createPicard(config);
    } else {
        ERR("nonlinear solver type `%s' unknown", type.c_str());
        std::abort();
    }

    return nullptr;
}

}  // namespace Nonlinear

}  // namespace MathLib

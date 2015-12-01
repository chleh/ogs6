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
createPicard(const BaseLib::ConfigTree &config)
{
    assert(config.get<std::string>("type", "") == "Picard");

    std::unique_ptr<MathLib::Nonlinear::Picard> pic(new MathLib::Nonlinear::Picard);

    auto const abs_tol = config.get_optional<double>("absolute_tolerance");
    if (abs_tol) pic->setAbsTolerance(*abs_tol);

    auto const rel_tol = config.get_optional<double>("relative_tolerance");
    if (rel_tol) pic->setRelTolerance(*rel_tol);

    auto const max_iter = config.get_optional<std::size_t>("max_iterations");
    if (max_iter) pic->setMaxIterations(*max_iter);

    auto const print = config.get_optional<bool>("print_errors");
    if (print) pic->printErrors(*print);

    auto const norm = config.get_optional<std::string>("norm_type");
    if (norm) {
        if      (*norm == "Norm1")   pic->setNormType(MathLib::VecNormType::NORM1);
        else if (*norm == "Norm1")   pic->setNormType(MathLib::VecNormType::NORM2);
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
createNonlinearSolver(const BaseLib::ConfigTree &config)
{
    auto const type = config.get_optional<std::string>("type");

    if (type) {
        if (*type == "Picard") {
            return createPicard(config);



        } else {
            ERR("nonlinear solver type `%s' unknown", type->c_str());
            std::abort();
        }
    } else {
        ERR("nonlinear solver type not given.");
        std::abort();
    }

    return nullptr;
}

}  // namespace Nonlinear

}  // namespace MathLib

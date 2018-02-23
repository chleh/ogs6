/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <functional>
#include <vector>

#include "ProcessLib/BoundaryCondition/ActiveBoundaryConditionCriterion.h"

namespace ProcessLib
{
class ActiveBCCriterionThNonnen final : public ActiveBoundaryConditionCriterion
{
public:
    struct Criterion
    {
        Criterion(unsigned cid, std::function<bool(double, double)>&& rel,
                  double v, std::function<double(double, double)>&& agg,
                  double agg_start_)
            : component_id(cid),
              relation(rel),
              value(v),
              aggregation(agg),
              agg_start(agg_start_)
        {
        }

        unsigned const component_id;
        std::function<bool(double, double)> const relation;
        double const value;
        std::function<double(double, double)> const aggregation;
        double const agg_start;
    };

    ActiveBCCriterionThNonnen(std::vector<Criterion>&& criteria)
        : _criteria(std::move(criteria))
    {
    }

    std::size_t getActiveBC(const std::size_t current_active_bc,
                            const GlobalVector& x, const double t,
                            const NumLib::LocalToGlobalIndexMap& dof_table,
                            const std::size_t mesh_id,
                            const std::vector<std::size_t>& node_ids) override;

private:
    std::vector<Criterion> const _criteria;
};

std::unique_ptr<ActiveBCCriterionThNonnen> createActiveBCCriterionThNonnen(
    BaseLib::ConfigTree const& config);

}  // ProcessLib

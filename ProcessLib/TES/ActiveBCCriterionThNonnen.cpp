/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ActiveBCCriterionThNonnen.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"

namespace ProcessLib
{
std::size_t ActiveBCCriterionThNonnen::getActiveBC(
    const std::size_t current_active_bc, const GlobalVector& x,
    const double /*t*/, const NumLib::LocalToGlobalIndexMap& dof_table,
    const std::size_t mesh_id, const std::vector<std::size_t>& node_ids)
{
    if (current_active_bc >= _criteria.size())
        OGS_FATAL("Active BC index (%u) exceeds maximum (%d)",
                  current_active_bc, _criteria.size() - 1);

    auto const& crit = _criteria[current_active_bc];

    if (crit.component_id >= dof_table.getNumberOfComponents())
        OGS_FATAL("Component id too large (%u, maximum is %u).",
                  crit.component_id, dof_table.getNumberOfComponents() - 1);

    if (node_ids.empty())
        OGS_FATAL("No node for measuring provided.");

    MeshLib::Location loc(mesh_id, MeshLib::MeshItemType::Node, 0);
    // bool first_iter = true;
    double agg_value = crit.agg_start;

    DBUG("Measuring component %u.", crit.component_id);
    for (auto const node_id : node_ids)
    {
        loc.item_id = node_id;
        auto const dof_idx = dof_table.getLocalIndex(
            loc, crit.component_id, x.getRangeBegin(), x.getRangeEnd());
        (void)dof_idx;
        DBUG("node: %4lu\tdof idx: %4lu\tvalue: %g",
             node_id,
             dof_idx,
             x.get(dof_idx));
        agg_value = crit.aggregation(agg_value, x.get(dof_idx));
    }

    INFO("Aggregated value for component %u: %g", crit.component_id, agg_value);

    if (crit.relation(agg_value, crit.value))
    {
        INFO("Switching criterion fulfilled.");
        return current_active_bc + 1;
    }

    return current_active_bc;
}

std::unique_ptr<ActiveBCCriterionThNonnen> createActiveBCCriterionThNonnen(
    BaseLib::ConfigTree const& config)
{
    config.checkConfigParameter("type", "ThNonnen");

    std::vector<ActiveBCCriterionThNonnen::Criterion> criteria;

    auto const criteria_config = config.getConfigSubtree("criteria");
    for (auto criterion : criteria_config.getConfigParameterList("criterion"))
    {
        auto const property =
            criterion.getConfigAttribute<std::string>("property");
        auto const relation =
            criterion.getConfigAttribute<std::string>("relation");
        auto const value = criterion.getConfigAttribute<double>("value");
        auto const aggregation =
            criterion.getConfigAttribute<std::string>("aggregation");

        unsigned component_id;
        if (property == "pressure")
            component_id = 0;
        else if (property == "temperature")
            component_id = 1;
        else if (property == "vapour_mass_fraction")
            component_id = 2;
        else
            OGS_FATAL("Unknown property: `%s'.", property.c_str());

        std::function<bool(double, double)> rel;
        if (relation == "greater_than")
            rel = std::greater<double>{};
        else if (relation == "less_than")
            rel = std::less<double>{};
        else
            OGS_FATAL("Unknown relation: `%s'.", relation.c_str());

        std::function<double(double, double)> agg;
        double agg_start;
        if (aggregation == "min")
        {
            agg = [](double a, double b) { return std::min(a, b); };
            agg_start = std::numeric_limits<double>::max();
        }
        else if (aggregation == "max")
        {
            agg = [](double a, double b) { return std::max(a, b); };
            agg_start = std::numeric_limits<double>::lowest();
        }
        else
        {
            OGS_FATAL("Unknown aggregation: `%s'.", aggregation.c_str());
        }

        criteria.emplace_back(component_id, std::move(rel), value,
                              std::move(agg), agg_start);
    }

    return std::unique_ptr<ActiveBCCriterionThNonnen>(
        new ActiveBCCriterionThNonnen(std::move(criteria)));
}

}  // ProcessLib

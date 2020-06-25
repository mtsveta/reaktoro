//
// Created by root on 24/09/19.
//

#include "SmartKineticResult.hpp"

namespace Reaktoro{

auto SmartKineticTiming::operator+=(const SmartKineticTiming &other) -> SmartKineticTiming&
{
    initialize += other.initialize;

    solve += other.solve;

    estimate += other.estimate;
    estimate_search += other.estimate_search;
    estimate_error_control += other.estimate_error_control;
    estimate_taylor += other.estimate_taylor;
    estimate_database_priority_update += other.estimate_database_priority_update;

    learn += other.learn;
    learn_storage += other.learn_storage;
    learn_integration += other.learn_integration;
    learn_equilibration += other.learn_equilibration;
    learn_reaction_rates += other.learn_reaction_rates;
    learn_chemical_properties += other.learn_chemical_properties;
    learn_error_control_matrices += other.learn_error_control_matrices;

    equilibrate += other.equilibrate;

    return *this;
}
auto SmartKineticResult::operator+=(const Reaktoro::SmartKineticResult &other) -> SmartKineticResult&
{
    timing += other.timing;
    equilibrium += other.equilibrium;
    smart_equilibrium += other.smart_equilibrium;

    return *this;
}

} // namespace Reaktoro
//
// Created by root on 24/09/19.
//

#include "SmartKineticResult.hpp"

namespace Reaktoro{

auto SmartKineticTiming::operator+=(const SmartKineticTiming &other) -> SmartKineticTiming&
{
    solve += other.solve;
    estimate += other.estimate;
    learn += other.learn;
    return *this;
}
auto SmartKineticResult::operator+=(const Reaktoro::SmartKineticResult &other) -> SmartKineticResult&
{
    timing += other.timing;
    return *this;
}

}
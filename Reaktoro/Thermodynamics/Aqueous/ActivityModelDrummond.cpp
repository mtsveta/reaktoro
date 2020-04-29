// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#include "ActivityModelDrummond.hpp"

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Aqueous/AqueousMixture.hpp>

namespace Reaktoro {

using std::log;

ActivityModelDrummond::ActivityModelDrummond()
{}

ActivityModelDrummond::ActivityModelDrummond(String gas)
: gas(gas)
{}

ActivityModelDrummond::ActivityModelDrummond(String gas, Params params)
: gas(gas), params(params)
{}

auto ActivityModelDrummond::build(const SpeciesList& species) const -> ActivityPropsFn
{
    // The index of the dissolved gas in the aqueous phase.
    const auto igas = species.indexWithFormula(gas);

    ActivityPropsFn fn = [=](ActivityPropsRef props, ActivityArgs args)
    {
        // The aqueous mixture and its state exported by a base aqueous activity model.
        const auto& mixture = std::any_cast<AqueousMixture>(args.extra[0]);
        const auto& state = std::any_cast<AqueousMixtureState>(args.extra[1]);

        const auto& [a1, a2, a3, a4, a5] = params;
        const auto& T = state.T;
        const auto& I = state.Is;
        const auto c1 = a1 + a2*T + a3/T;
        const auto c2 = a4 + a5*T;
        props.ln_g[igas] = c1 * I - c2 * I/(I + 1);
        props.ln_a[igas] = props.ln_g[igas] + log(state.m[igas]);
    };

    return fn;
}

} // namespace Reaktoro

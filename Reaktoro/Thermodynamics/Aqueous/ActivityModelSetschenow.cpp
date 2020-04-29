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

#include "ActivityModelSetschenow.hpp"

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Aqueous/AqueousMixture.hpp>

namespace Reaktoro {

auto ActivityModelSetschenow(const AqueousMixture& mixture, double b) -> AqueousActivityModel
{
    // The value of ln(10)
    const double ln10 = 2.30258509299;

    AqueousActivityModel f = [=](const AqueousMixtureState& state)
    {
        // The effective ionic strength of the aqueous mixture
        const auto& I = state.Ie;

        // The activity coefficient of the given species (in molality scale)
        real ln_gi = ln10 * b * I;

        return ln_gi;
    };

    return f;
}

} // namespace Reaktoro
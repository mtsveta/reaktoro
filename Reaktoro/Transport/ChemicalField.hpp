// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

#pragma once

// C++ includes
#include <memory>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/StringList.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Math/Matrix.hpp>
#include <Reaktoro/Core/ChemicalOutput.hpp>

namespace Reaktoro {

/// Use this class for storing the collection of states and their properties
class ChemicalField
{
public:
    /// Construct a ChemicalField instance when the number of the states(cells) and system is provided.
    ChemicalField(Index size, const ChemicalSystem& system);

    /// Construct a ChemicalField instance when the number of the states(cells) and an instance of one state.
    ChemicalField(Index size, const ChemicalState& state);

    /// Construct a copy of a ChemicalField instance.
    ChemicalField(const ChemicalField& other);

    /// Destroy this ChemicalField instance.
    virtual ~ChemicalField();

    /// Assign a copy of an ChemicalField instance.
    auto operator=(ChemicalField other) -> ChemicalField&;

    /// Return the number of states in chamical field
    auto size() const -> Index;

    /// Set all the states of the field by the one provided
    auto set(const ChemicalState& state) -> void;

    /// Fetch the vector with temperatures in the states
    auto temperature(VectorRef values) -> void;

    /// Fetch the vector with pressures in the states
    auto pressure(VectorRef values) -> void;

    /// Fetch the vector with element amounts in the states
    auto elementAmounts(VectorRef values) -> void;

    auto output(const std::string& filename, const StringList& quantities) -> void;

    // Return the collection with properties
    auto properties() -> std::vector<ChemicalProperties>&;

    auto states() -> std::vector<ChemicalState>&;

private:

    struct Impl;

    std::unique_ptr<Impl> pimpl;

};

} // namespace Reaktoro

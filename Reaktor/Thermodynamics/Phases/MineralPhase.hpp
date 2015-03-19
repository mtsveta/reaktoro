// Reaktor is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#pragma once

// C++ includes
#include <tuple>
#include <vector>

// Reaktor includes
#include <Reaktor/Thermodynamics/Activity/MineralActivity.hpp>
#include <Reaktor/Thermodynamics/Mixtures/MineralMixture.hpp>

namespace Reaktor {

// Reaktor forward declarations
class Phase;

/// Class that defines an mineral phase
class MineralPhase : public MineralMixture
{
public:
    /// Construct a default MineralPhase instance.
    MineralPhase();

    /// Construct an MineralPhase instance with given species.
    explicit MineralPhase(const std::vector<MineralSpecies>& species);

    /// Set the activity model of a species.
    /// @param species The name of the species
    /// @param activity The activity function
    /// @see MineralActivity
    auto setActivityModel(const std::string& species, const MineralActivity& activity) -> void;

    /// Set the activity model of the species to be the ideal one.
    /// @param species The name of species to have its activity model set
    auto setActivityModelIdeal(const std::string& species) -> void;

    /// Calculate the concentrations of the mineral species.
    /// @param n The molar abundance of the species
    /// @return The concentrations of the mineral species
    auto concentrations(const Vector& n) const -> Vector;

    /// Calculate the activities of the mineral species and its molar derivatives.
    /// @param T The temperature used for the calculation (in units of K)
    /// @param P The pressure used for the calculation (in units of bar)
    /// @param n The molar composition of the mineral phase
    /// @return The activities of the mineral species and their molar derivatives
    auto activities(double T, double P, const Vector& n) const -> ChemicalVector;

private:
    /// The mineral activity functions
    std::vector<MineralActivity> activities$;
};

/// Create a Phase instance from a MineralPhase instance
/// @param phase The MineralPhase instance
/// @return A Phase instance created from the given mineral phase
auto createPhase(const MineralPhase& phase) -> Phase;

} // namespace Reaktor

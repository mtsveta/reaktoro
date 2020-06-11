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
#include <string>

// Reaktoro includes
#include <Reaktoro/Math/Matrix.hpp>
#include <Reaktoro/Common/ScalarTypes.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalProperties;
class ChemicalSystem;
class EquilibriumProperties;

/// Provides a computational representation of the state of a multiphase chemical system.
/// The chemical state of a multiphase system is defined by its temperature @f$(T)@f$,
/// pressure @f$(P)@f$, and molar composition @f$(\mathbf{n})@f$.
///
/// **Usage**
///
/// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// // Create a ChemicalState instance, where system is a ChemicalSystem instance
/// ChemicalState state(system);
///
/// // Set the temperature and pressure states
/// state.setTemperature(60.0, "celsius");
/// state.setPressure(  180.0, "bar");
///
/// // Set the amount of some species
/// state.set( "H2O(l)",  1.0, "kg");
/// state.set(    "Na+",  1.0, "mol");
/// state.set(    "Cl-",  1.0, "mol");
/// state.set("CO2(aq)",  0.5, "mol");
/// state.set("Calcite", 10.0, "g");
///
/// // Output the chemical state instance
/// std::cout << state << std::endl;
/// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// @see ChemicalSystem
/// @ingroup Core
class ChemicalState
{
public:

    /// Construct a ChemicalState instance with standard conditions
    /// This constructor creates an instance of ChemicalState with
    /// temperature 25 °C, pressure 1 bar, and zero species amounts.
    /// @param system The chemical system instance
    explicit ChemicalState(const ChemicalSystem& system);

    /// Construct a copy of a ChemicalState instance
    ChemicalState(const ChemicalState& other);

    /// Destroy this ChemicalState instance
    virtual ~ChemicalState();

    /// Assign a ChemicalState instance to this instance
    auto operator=(ChemicalState other) -> ChemicalState&;

    /// Set the temperature of the chemical state (in units of K)
    auto setTemperature(double val) -> void;

    /// Set the temperature of the chemical state with given units
    auto setTemperature(double val, std::string units) -> void;

    /// Set the pressure of the chemical state (in units of Pa)
    auto setPressure(double val) -> void;

    /// Set the pressure of the chemical state with given units
    auto setPressure(double val, std::string units) -> void;

    /// Set the molar amounts of the species with a single value (in units of mol)
    /// @param val The single molar amounts of the species
    auto setSpeciesAmounts(double val) -> void;

    /// Set the molar amounts of the species (in units of mol)
    /// @param n The vector of molar amounts of the species
    auto setSpeciesAmounts(VectorConstRef n) -> void;

    /// Set the molar amounts of species given by their indices (in units of mol)
    /// @param n The vector of molar amounts of the species
    /// @param indices The indices of the species to be set
    auto setSpeciesAmounts(VectorConstRef n, const Indices& indices) -> void;

    /// Set the molar amount of a species (in units of mol)
    /// @param index The index of the species
    /// @param amount The molar amount of the species
    auto setSpeciesAmount(Index index, double amount) -> void;

    /// Set the molar amount of a species (in units of mol)
    /// @param name The name of the species
    /// @param amount The amount of the species
    auto setSpeciesAmount(std::string name, double amount) -> void;

    /// Set the amount of a species with given units
    /// @param index The index of the species
    /// @param amount The amount of the species
    /// @param units The units of the amount (must be convertible to either mol or gram)
    auto setSpeciesAmount(Index index, double amount, std::string units) -> void;

    /// Set the amount of a species with given units
    /// @param name The name of the species
    /// @param amount The amount of the species
    /// @param units The units of the amount (must be convertible to either mol or gram)
    auto setSpeciesAmount(std::string name, double amount, std::string units) -> void;

    /// Set the mass of a species (in units of kg)
    /// @param index The index of the species
    /// @param mass The mass of the species
    auto setSpeciesMass(Index index, double mass) -> void;

    /// Set the mass of a species (in units of kg)
    /// @param name The name of the species
    /// @param mass The mass of the species
    auto setSpeciesMass(std::string name, double mass) -> void;

    /// Set the mass of a species with given units
    /// @param index The index of the species
    /// @param mass The mass of the species
    /// @param units The units of the mass
    auto setSpeciesMass(Index index, double mass, std::string units) -> void;

    /// Set the mass of a species with given units
    /// @param name The name of the species
    /// @param mass The mass of the species
    /// @param units The units of the mass
    auto setSpeciesMass(std::string name, double mass, std::string units) -> void;

    /// Set the dual potentials of the species (in units of J/mol).
    ///
    /// The dual potentials of the species are the Lagrange multipliers with
    /// respect to the positive bound constraints on the molar amounts of the
    /// species in a chemical equilibrium calculation. They can be seen as
    /// measures of stability of a species at equilibrium, with values closer
    /// to zero meaning more stability.
    /// @param z The Lagrange multipliers with respect to the positive constraints.
    [[deprecated("Use state.equilibrium().setSpeciesStabilities() instead, where state is a ChemicalState object.")]]
    auto setSpeciesDualPotentials(VectorConstRef z) -> void;

    /// Set the dual potentials of the elements (in units of J/mol).
    ///
    /// The dual potentials of the elements are the Lagrange multipliers with
    /// respect to the balance constraints on the molar amounts of the elements.
    /// They can be seen as a dual chemical potential of elements.
    /// @param values The Lagrange multipliers with respect to the balance constraints.
    [[deprecated("Use state.equilibrium().setElementChemicalPotentials() instead, where state is a ChemicalState object.")]]
    auto setElementDualPotentials(VectorConstRef y) -> void;

    /// Scale the molar amounts of the species by a given scalar.
    /// @param scalar The scale factor of the molar amounts
    auto scaleSpeciesAmounts(double scalar) -> void;

    /// Scale the molar amounts of the species in a phase by a given scalar.
    /// @param index The index of the phase
    /// @param scalar The scale factor of the molar amounts
    auto scaleSpeciesAmountsInPhase(Index index, double scalar) -> void;

    /// Scale the volume of a phase by adjusting the molar amounts of its species.
    /// @param index The index of the phase
    /// @param volume The volume of the phase (in units of m3)
    auto scalePhaseVolume(Index index, double volume) -> void;

    /// Scale the volume of a phase by adjusting the molar amounts of its species.
    /// @param index The index of the phase
    /// @param volume The volume of the phase
    /// @param units The units of the volume of the phase
    auto scalePhaseVolume(Index index, double volume, std::string units) -> void;

    /// Scale the volume of a phase by adjusting the molar amounts of its species.
    /// @param name The name of the phase
    /// @param volume The volume of the phase (in units of m3)
    auto scalePhaseVolume(std::string name, double volume) -> void;

    /// Scale the volume of a phase by adjusting the molar amounts of its species.
    /// @param name The name of the phase
    /// @param volume The volume of the phase
    /// @param units The units of the volume of the phase
    auto scalePhaseVolume(std::string name, double volume, std::string units) -> void;

    /// Scale the fluid volume of the chemical system.
    /// This method scales each fluid phase by the same factor.
    /// This factor is defined as the prescribed volume divided by
    /// the current volume of the fluid phases.
    /// @param volume The scaled fluid volume (in units of m3)
    auto scaleFluidVolume(double volume) -> void;

    /// Scale the fluid volume of the chemical system with given units.
    /// This method scales each fluid phase by the same factor.
    /// This factor is defined as the prescribed volume divided by
    /// the current volume of the fluid phases.
    /// @param volume The scaled fluid volume
    /// @param units The volume units
    auto scaleFluidVolume(double volume, std::string units) -> void;

    /// Scale the solid volume of the chemical system.
    /// This method scales each solid phase by the same factor.
    /// This factor is defined as the prescribed volume divided by
    /// the current volume of the solid phases.
    /// @param volume The scaled solid volume (in units of m3)
    auto scaleSolidVolume(double volume) -> void;

    /// Scale the solid volume of the chemical system with given units.
    /// This method scales each solid phase by the same factor.
    /// This factor is defined as the prescribed volume divided by
    /// the current volume of the solid phases.
    /// @param volume The scaled solid volume
    /// @param units The volume units
    auto scaleSolidVolume(double volume, std::string units) -> void;

    /// Scale the volume of the chemical system by adjusting the molar amounts of all species equally.
    /// @param volume The volume of the chemical system (in units of m3)
    auto scaleVolume(double volume) -> void;

    /// Scale the volume of the chemical system by adjusting the molar amounts of all species equally.
    /// @param volume The volume of the chemical system
    /// @param units The volume units
    auto scaleVolume(double volume, std::string units) -> void;

    /// Return the chemical system instance
    auto system() const -> const ChemicalSystem&;

    /// Return the temperature of the chemical state (in units of K)
    auto temperature() const -> double;

    /// Return the pressure of the chemical state (in units of Pa)
    auto pressure() const -> double;

    /// Return the molar amounts of the species (in units of mol)
    auto speciesAmounts() const -> VectorConstRef;

    /// Return the molar amounts of given species (in units of mol)
    /// @param indices The indices of the species
    auto speciesAmounts(const Indices& indices) const -> Vector;

    /// Return the molar mass of a chemical species (in units of mol)
    /// @param index The index of the species
    auto speciesAmount(Index index) const -> double;

    /// Return the molar amount of a chemical species (in units of mol)
    /// @param name The name of the species
    auto speciesAmount(std::string name) const -> double;

    /// Return the amount of a chemical species with given molar units
    /// @param index The index of the species
    /// @param units The units of the species amount
    auto speciesAmount(Index index, std::string units) const -> double;

    /// Return the amount of a chemical species with given molar units
    /// @param name The name of the species
    /// @param units The units of the species amount
    auto speciesAmount(std::string name, std::string units) const -> double;

    /// Return the dual potentials of the species (in units of J/mol)
    [[deprecated("Use state.equilibrium().speciesStabilities() instead, where state is a ChemicalState object.")]]
    auto speciesDualPotentials() const -> VectorConstRef;

    /// Return the molar amounts of the elements (in units of mol)
    auto elementAmounts() const -> Vector;

    /// Return the molar amounts of the elements in a phase (in units of mol)
    /// @param index The index of the phase
    auto elementAmountsInPhase(Index index) const -> Vector;

    /// Return the molar amounts of the elements in a set of species (in units of mol)
    /// @param indices The indices of the species
    auto elementAmountsInSpecies(const Indices& indices) const -> Vector;

    /// Return the molar amount of an element (in units of mol)
    /// @param index The index of the element
    auto elementAmount(Index index) const -> double;

    /// Return the molar amount of an element (in units of mol)
    /// @param name The name of the element
    auto elementAmount(std::string name) const -> double;

    /// Return the amount of an element with given units.
    /// @param index The index of the element
    /// @param units The units of the element amount
    auto elementAmount(Index index, std::string units) const -> double;

    /// Return the amount of an element with given units.
    /// @param name The name of the element
    /// @param units The units of the element amount
    auto elementAmount(std::string name, std::string units) const -> double;

    /// Return the molar amount of an element in a given phase (in units of mol).
    /// @param ielement The index of the element
    /// @param iphase The index of the phase
    auto elementAmountInPhase(Index ielement, Index iphase) const -> double;

    /// Return the molar amount of an element in a given phase (in units of mol).
    /// @param element The name of the element
    /// @param phase The name of the phase
    auto elementAmountInPhase(std::string element, std::string phase) const -> double;

    /// Return the amount of an element in a given phase with given units.
    /// @param ielement The index of the element
    /// @param iphase The index of the phase
    /// @param units The units of the element amount
    auto elementAmountInPhase(Index ielement, Index iphase, std::string units) const -> double;

    /// Return the amount of an element in a given phase with given units.
    /// @param element The name of the element
    /// @param phase The name of the phase
    /// @param units The units of the element amount
    auto elementAmountInPhase(std::string element, std::string phase, std::string units) const -> double;

    /// Return the molar amount of an element in a set of species (in units of mol).
    /// @param ielement The index of the element
    /// @param ispecies The indices of the species
    auto elementAmountInSpecies(Index ielement, const Indices& ispecies) const -> double;

    /// Return the amount of an element in a set of species with given units
    /// @param ielement The index of the element
    /// @param ispecies The indices of the species
    /// @param units The units of the element amount
    auto elementAmountInSpecies(Index ielement, const Indices& ispecies, std::string units) const -> double;

    /// Return the dual potentials of the elements (in units of J/mol).
    [[deprecated("Use state.equilibrium().elementChemicalPotentials() instead, where state is a ChemicalState object.")]]
    auto elementDualPotentials() const -> VectorConstRef;

    /// Return the molar amount of a phase (in units of mol).
    /// @param index The index of the phase
    auto phaseAmount(Index index) const -> double;

    /// Return the molar amount of a phase (in units of mol).
    /// @param name The name of the phase
    auto phaseAmount(std::string name) const -> double;

    /// Return the molar amount of a phase with given units.
    /// @param index The index of the phase
    /// @param units The units of the phase amount
    auto phaseAmount(Index index, std::string units) const -> double;

    /// Return the molar amount of a phase with given units.
    /// @param name The name of the phase
    /// @param units The units of the phase amount
    auto phaseAmount(std::string name, std::string units) const -> double;

    /// Return the stability indices of the phases with respect to chemical equilibrium.
    [[deprecated("Use phaseStabilityIndices(state) instead, where state is a ChemicalState object.")]]
    auto phaseStabilityIndices() const -> Vector;

    /// Return the chemical properties of the system.
    auto properties() const -> ChemicalProperties;

    /// Return the equilibrium properties of a calculated chemical equilibrium state.
    auto equilibrium() const -> const EquilibriumProperties&;

    /// Return the equilibrium properties of a calculated chemical equilibrium state.
    auto equilibrium() -> EquilibriumProperties&;

    /// Output the ChemicalState instance to a file.
    auto output(std::string filename) const -> void;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

/// Provide access to specific properties of a chemical equilibrium state.
class EquilibriumProperties
{
public:
    /// Construct a EquilibriumProperties instance.
    EquilibriumProperties(const ChemicalSystem& system);

    /// Construct a copy of a EquilibriumProperties instance
    EquilibriumProperties(const EquilibriumProperties& other);

    /// Destroy this EquilibriumProperties instance
    virtual ~EquilibriumProperties();

    /// Assign a EquilibriumProperties instance to this instance
    auto operator=(EquilibriumProperties other) -> EquilibriumProperties&;

    /// Set the indices of the equilibrium species partitioned as (primary, secondary).
    /// @param ips The indices of the equilibrium species ordered (primary, secondary)
    /// @param kp The number of primary species
    auto setIndicesEquilibriumSpecies(VectorXiConstRef ips, Index kp) -> void;

    /// Set the chemical potentials of the species in the equilibrium state (in units of J/mol)
    auto setSpeciesChemicalPotentials(VectorConstRef u) -> void;

    /// Set the chemical potentials of the elements in the equilibrium state (in units of J/mol)
    auto setElementChemicalPotentials(VectorConstRef y) -> void;

    /// Set the stabilities of the species in the equilibrium state (in units of J/mol)
    auto setSpeciesStabilities(VectorConstRef z) -> void;

    /// Return the number of equilibrium species.
    auto numEquilibriumSpecies() const -> Index;

    /// Return the number of primary equilibrium species.
    auto numPrimarySpecies() const -> Index;

    /// Return the number of secondary equilibrium species.
    auto numSecondarySpecies() const -> Index;

    /// Return the indices of the equilibrium species.
    auto indicesEquilibriumSpecies() const -> VectorXiConstRef;

    /// Return the indices of the primary equilibrium species.
    auto indicesPrimarySpecies() const -> VectorXiConstRef;

    /// Return the indices of the secondary equilibrium species.
    auto indicesSecondarySpecies() const -> VectorXiConstRef;

    /// Return the chemical potentials of the species (in units of J/mol).
    /// The chemical potentials of the species are the gradient of the Gibbs energy
    /// function with respect to species amounts.
    auto speciesChemicalPotentials() const -> VectorConstRef;

    /// Return the chemical potentials of the elements (in units of J/mol).
    /// The chemical potentials of the elements are the Lagrange multipliers with
    /// respect to the mass conservation constraints on the amounts of the elements.
    auto elementChemicalPotentials() const -> VectorConstRef;

    /// Return the stabilities of the species (in units of J/mol)
    /// The stabilities of the species are the slack variables with
    /// respect to the non-negative bound constraints on the amounts of the
    /// species in a chemical equilibrium calculation. They can be seen as
    /// measures of stability of a species at equilibrium, with values closer
    /// to zero meaning more stable.
    auto speciesStabilities() const -> VectorConstRef;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

/// Outputs a ChemicalState instance.
auto operator<<(std::ostream& out, const ChemicalState& state) -> std::ostream&;

/// Add two ChemicalState instances.
auto operator+(const ChemicalState& l, const ChemicalState& r) -> ChemicalState;

/// Multiply a ChemicalState instance by a scalar (from the left).
auto operator*(double scalar, const ChemicalState& state) -> ChemicalState;

/// Multiply a ChemicalState instance by a scalar (from the right).
auto operator*(const ChemicalState& state, double scalar) -> ChemicalState;

/// Return the stability indices of the phases with respect to a chemical equilibrium state.
/// The stability index of a stable phase at chemical equilibrium should be, or
/// very close to, zero. A negative stability index indicates that the
/// corresponding phase is under-stable, while a positive index indicates the
/// phase is over-stable.
auto phaseStabilityIndices(const ChemicalState& state) -> Vector;

} // namespace Reaktoro

// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021
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

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

/// Return the solubility of CO2 for given temperature, pressure, and amount of sodium chloride
auto solubilityCO2(const ChemicalSystem& system, const real& T, const real& P, double mNaCl) -> real
{
    // Define initial CO2 amount
    double n0CO2g = 10.0;
    // Define mass of water
    double mH02 = 1.0;

    // Define initial equilibrium state
    ChemicalState state(system);
    state.setTemperature(T, "celsius");
    state.setPressure(P, "bar");
    state.setSpeciesMass("H2O(l)", mH02, "kg");
    state.setSpeciesAmount("CO2(g)", n0CO2g, "mol");
    state.setSpeciesAmount("Na+", mNaCl, "mol");
    state.setSpeciesAmount("Cl-", mNaCl, "mol");

    // Define equilibrium solver and equilibrate given initial state
    EquilibriumSolver solver(system);
    solver.solve(state);

    // Fetch properties of aqueous phase from the generated state
    AqueousProps aqprops = AqueousProps(state);

    return aqprops.elementMolality("C");
}
int main()
{
    // Initialize a thermodynamic database
    SupcrtDatabase db("supcrt98.xml");

    // Create an aqueous phase
    AqueousPhase aqueousphase(speciate("H O C Na Cl"));
    aqueousphase.setActivityModel(chain(
            ActivityModelHKF(),
            ActivityModelDrummond("CO2")
    ));

    // Create a gaseous phase
    GaseousPhase gaseousphase("CO2(g)");
    gaseousphase.setActivityModel(ActivityModelPengRobinson());

    // Collecting all above-defined phases
    Phases phases(db);
    phases.add(aqueousphase);
    phases.add(gaseousphase);

    // Construct the chemical system
    ChemicalSystem system(phases);

    // Define range of temperatures
    std::vector<real> temperatures{25, 50, 100, 300};
    // Define an array with CO2 solubilities
    std::vector<real> cCO2;

    // Given pressure
    real P = 1.0;
    // Given initial amount of CO2
    double mNaCl = 10.0;

    // Evaluate solubility if the the given temperature range
    for (auto T : temperatures)
        cCO2.push_back(solubilityCO2(system, T, P, mNaCl));

    // Print the temperatures and corresponding solubilities
    std::cout << std::setw(20) << "Temperatures" << std::setw(20) << "CO2 solubilities" << std::endl;

    // Print obtained solubility for the given temperature
    for (auto i = 0; i < temperatures.size(); ++i)
        std::cout << std::setw(20) << temperatures[i] << std::setw(20) << cCO2[i] << std::endl;

    return 0;
}

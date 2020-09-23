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

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    ChemicalEditor editor;
    AqueousPhase aqueousphase = editor.addAqueousPhaseWithElements("H O Na Cl");
    aqueousphase.setActivityModelSetschenow("NaCl(aq)", 1.0);
    editor.addMineralPhase("Halite");

    ChemicalSystem system(editor);

    EquilibriumProblem problem(system);
    problem.setTemperature(30.0, "celsius");
    problem.setPressure(300.0, "bar");
    problem.add("H2O", 1, "kg");
    problem.add("NaCl", 14.0, "mol");

    ChemicalState state = equilibrate(problem);

    std::cout << state << std::endl;
}
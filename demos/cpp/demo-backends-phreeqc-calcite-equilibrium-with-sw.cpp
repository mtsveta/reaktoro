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

// This string defines a PHREEQC script problem.
// This problem was taken from the official PHREEQC example named ex1.
const std::string ex3 = R"(
TITLE Example 3, part A.--Calcite equilibrium at log Pco2 = -2.0 and 25C.        
SOLUTION 1  Pure water
        pH      7.0
        temp    25.0
EQUILIBRIUM_PHASES
        CO2(g)          -2.0    
        Calcite         0.0
SAVE solution 1
END
TITLE Example 3, part B.--Definition of seawater.
SOLUTION 2  Seawater
        units   ppm
        pH      8.22
        pe      8.451
        density 1.023
        temp    25.0
        Ca              412.3
        Mg              1291.8
        Na              10768.0
        K               399.1
        Si              4.28
        Cl              19353.0
        Alkalinity      141.682 as HCO3
        S(6)            2712.0
END
TITLE Example 3, part C.--Mix 70% groundwater, 30% seawater.
MIX 1
        1      0.7
        2      0.3
SAVE solution   3
END
TITLE Example 3, part D.--Equilibrate mixture with calcite and dolomite.
EQUILIBRIUM_PHASES 1
        Calcite         0.0
        Dolomite        0.0
USE solution 3
END
TITLE Example 3, part E.--Equilibrate mixture with calcite only.
EQUILIBRIUM_PHASES 2
        Calcite         0.0
USE solution 3
END
)";

int main()
{
    // Initialize a Phreeqc instance with the official phreeqc.dat database file
    Phreeqc phreeqc("databases/phreeqc/phreeqc.dat");

    // Execute a PHREEQC script defining a geochemical problem.
    // Here this script is actually embedded into a string named `ex1`.
    // However, `ex3` could also be a string containing the path to a script file.
    // Method execute will automatically identify when the contents are embedded in the string and
    // when the string is actually a path to a script file.
    phreeqc.execute(ex3);

    // Initialize a ChemicalSystem instance using the current state of the Phreeqc instance.
    // This will allow the use of both PHREEQC thermodynamic data and PHREEQC activity models
    // in the subsequent equilibrium calculations using Reaktoro's algorithms.
    ChemicalSystem system(phreeqc);

    // Initialize an ChemicalState instance using the current state of the Phreeqc instance.
    ChemicalState state = phreeqc.state(system);

    // Print such equilibrium state using Reaktoro format.
    std::cout << state << std::endl;
}


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

#include <PyReaktoro/PyReaktoro.hpp>

// Reaktoro includes
#include <Reaktoro/Kinetics/SmartKineticOptions.hpp>

namespace Reaktoro {

    void exportSmartKineticOptions(py::module& m)
    {
        py::class_<KineticOutputOptions>(m, "KineticOutputOptions")
                .def_readwrite("active", &KineticOutputOptions::active)
                .def_readwrite("format", &KineticOutputOptions::format)
                ;

        py::class_<SmartKineticOptions>(m, "SmartKineticOptions")
                .def(py::init<>())
                .def_readwrite("reltol", &SmartKineticOptions::reltol)
                .def_readwrite("abstol", &SmartKineticOptions::abstol)
                .def_readwrite("cutoff", &SmartKineticOptions::cutoff)
                .def_readwrite("mole_fraction_cutoff", &SmartKineticOptions::mole_fraction_cutoff)
                .def_readwrite("use_smart_equilibrium_solver", &SmartKineticOptions::use_smart_equilibrium_solver)

                .def_readwrite("equilibrium", &SmartKineticOptions::equilibrium)
                .def_readwrite("smart_equilibrium", &SmartKineticOptions::smart_equilibrium)
                .def_readwrite("learning", &SmartKineticOptions::learning)
                ;
    }

} // namespace Reaktoro

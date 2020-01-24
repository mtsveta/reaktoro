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
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Transport/ReactiveTransportSolver.hpp>
#include <Reaktoro/Transport/Mesh.hpp>

namespace Reaktoro {

void exportReactiveTransportSolver(py::module& m)
{
    py::class_<ReactiveTransportSolver>(m, "ReactiveTransportSolver")
        .def(py::init<const ChemicalSystem&>())
        .def("setMesh", &ReactiveTransportSolver::setMesh)
        .def("setVelocity", &ReactiveTransportSolver::setVelocity)
        .def("setDiffusionCoeff", &ReactiveTransportSolver::setDiffusionCoeff)
        .def("setBoundaryState", &ReactiveTransportSolver::setBoundaryState)
        .def("setTimeStep", &ReactiveTransportSolver::setTimeStep)
        .def("system", &ReactiveTransportSolver::system, py::return_value_policy::reference_internal)
        .def("output", &ReactiveTransportSolver::output)
        .def("initialize", &ReactiveTransportSolver::initialize)
        .def("step", &ReactiveTransportSolver::step)
        ;
}

} // namespace Reaktoro

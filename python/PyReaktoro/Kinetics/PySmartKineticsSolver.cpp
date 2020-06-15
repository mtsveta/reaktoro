// Reaktoro is a unified framework for modeling chemically reactive systems.
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
#include <Reaktoro/Core/ReactionSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Kinetics/SmartKineticOptions.hpp>
#include <Reaktoro/Kinetics/SmartKineticSolver.hpp>
#include <Reaktoro/Kinetics/SmartKineticResult.hpp>

namespace Reaktoro {

    void exportSmartKineticSolver(py::module& m)
    {
        auto initialize1 = static_cast<void(SmartKineticSolver::*)(ChemicalState&, double)>(&SmartKineticSolver::initialize);
        auto initialize2 = static_cast<void(SmartKineticSolver::*)(ChemicalState&, double, VectorConstRef)>(&SmartKineticSolver::initialize);

        auto solve1 = static_cast<void(SmartKineticSolver::*)(ChemicalState&, double, double, VectorConstRef)>(&SmartKineticSolver::solve);

        py::class_<SmartKineticSolver>(m, "SmartKineticSolver")
                .def(py::init<const ReactionSystem&, const Partition&>())
                .def("setOptions", &SmartKineticSolver::setOptions)
                .def("addSource", &SmartKineticSolver::addSource)
                .def("addPhaseSink", &SmartKineticSolver::addPhaseSink)
                .def("addFluidSink", &SmartKineticSolver::addFluidSink)
                .def("addSolidSink", &SmartKineticSolver::addSolidSink)
                .def("initialize", initialize1)
                .def("initialize", initialize2)
                .def("solve", solve1)
                .def("result", &SmartKineticSolver::result, py::return_value_policy::reference_internal)

                        // DEPRECATED METHODS: TO BE REMOVED
                .def("setPartition", &SmartKineticSolver::setPartition)
                ;
    }

} // namespace Reaktoro

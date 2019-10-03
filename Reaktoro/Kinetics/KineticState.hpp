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

namespace Reaktoro {

struct KineticState
{
    /// Start time of the integration t_k
    double tk;
    /// Finish time of the integration t_{k+1}
    double tk1;

    /// Initial species' amount
    Vector n0;
    /// Species' amount on the moment t_k
    Vector nk;
    /// Species' amount on the moment t_{k+1}
    Vector nk1;

    /// The partial derivatives @f$\left.\frac{\partial n}{\partial n_0}\right|@f$.
    /// These derivatives provide a measure of how much the equilibrium amounts of the species,
    /// @f$n@f$, change with an infinitesimal change in their initial condition, @f$n_0@f$.
    /// They are used when solving kinetic paths of the cells that have similar initial condition
    /// as the one already processed.
    Matrix dndn0;

    /// The right-hand size of ODE @f$\left.f(n) = \frac{\partial n}{\partial t}\right|@f$.
    Vector dndt;

    /// Jacobian of the the right-hand size f of ODE, i.e., @f$\left.J(n) = \frac{\partial f}{\partial n}\right|@f$.
    Matrix dfdn;
};

} // namespace Reaktoro
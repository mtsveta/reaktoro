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

    /// Start time of the integration
    double t0;
    /// Finish time of the integration
    double t;

    /// Initial species' amount at time t_0
    Vector u0;
    /// Species' amount on the moment t_k
    Vector u;

    /// The partial derivatives @f$\left.\frac{\partial u}{\partial u_0}\right|@f$.
    /// These derivatives provide a measure of how much @f$u@f$ changes with an infinitesimal
    /// change in their initial condition @f$u_0@f$.
    Matrix dudu0;

    /// The right-hand size of ODE @f$\left.f(u) = \frac{\partial u}{\partial t}\right|@f$.
    Vector dudt;

    /// Jacobian of the the right-hand size f of ODE, i.e., @f$\left.J(f(u)) = \frac{\partial f}{\partial u}\right|@f$.
    Matrix dfdu;
};

} // namespace Reaktoro
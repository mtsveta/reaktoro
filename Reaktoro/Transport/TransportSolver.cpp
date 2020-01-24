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

#include "TransportSolver.hpp"

// C++ includes
//#include <iomanip>
#include <iostream>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>

namespace Reaktoro {
namespace internal {

template<typename ReturnType, typename TridiagonalMatrixType>
auto row(TridiagonalMatrixType&& mat, Index index) -> ReturnType
{
    auto n = mat.size();
    auto data = mat.data();
    auto length = data.size();
    return (index == 0) ? data.segment(1, 2) : (index == n - 1) ?
        data.segment(length - 4, 2) : data.segment(3 * index, 3);
}

} // namespace internal

ChemicalField::ChemicalField(Index size, const ChemicalSystem& system)
: m_size(size),
  m_system(system),
  m_states(size, ChemicalState(system)),
  m_properties(size, ChemicalProperties(system))
{}

ChemicalField::ChemicalField(Index size, const ChemicalState& state)
: m_size(size),
  m_system(state.system()),
  m_states(size, state),
  m_properties(size, state.properties())
{}

auto ChemicalField::set(const ChemicalState& state) -> void
{
    for(auto& item : m_states)
        item = state;
}

auto ChemicalField::temperature(VectorRef values) -> void
{
    const Index len = size();
    for(Index i = 0; i < len; ++i)
        values[i] = m_states[i].temperature();
}

auto ChemicalField::pressure(VectorRef values) -> void
{
    const Index len = size();
    for(Index i = 0; i < len; ++i)
        values[i] = m_states[i].pressure();
}

auto ChemicalField::elementAmounts(VectorRef values) -> void
{
    const Index len = size();
    const Index num_elements = m_system.numElements();
    Index offset = 0;
    for(Index i = 0; i < len; ++i, offset += num_elements)
        values.segment(offset, num_elements) = m_states[i].elementAmounts();
}

auto ChemicalField::output(std::string filename, StringList quantities) -> void
{
    ChemicalOutput out(m_system);
    out.filename(filename);
    for(auto quantity : quantities)
        out.add(quantity);

}

auto TridiagonalMatrix::resize(Index size) -> void
{
    m_size = size;
    m_data.conservativeResize(size * 3);
}

auto TridiagonalMatrix::factorize() -> void
{
    const Index n = size();

    auto prev = row(0).data(); // iterator to previous row
    auto curr = row(1).data(); // iterator to current row

    for(Index i = 1; i < n; ++i, prev = curr, curr += 3)
    {
        const auto& b_prev = prev[1]; // `b` value on the previous row
        const auto& c_prev = prev[2]; // `c` value on the previous row

        auto& a_curr = curr[0]; // `a` value on the current row
        auto& b_curr = curr[1]; // `b` value on the current row

        a_curr /= b_prev; // update the a-diagonal in the tridiagonal matrix
        b_curr -= a_curr * c_prev; // update the b-diagonal in the tridiagonal matrix
    }
}

auto TridiagonalMatrix::solve(VectorRef x, VectorConstRef d) const -> void
{
    const Index n = size();

    auto curr = row(1).data(); // iterator to current row

    //-------------------------------------------------------------------------
    // Perform the forward solve with the L factor of the LU factorization
    //-------------------------------------------------------------------------
    x[0] = d[0];

    for(Index i = 1; i < n; ++i, curr += 3)
    {
        const auto& a = curr[0]; // `a` value on the current row

        x[i] = d[i] - a * x[i - 1];
    }

    curr -= 3; // step back so that curr points to the last row
    const auto& bn = curr[1]; // `b` value on the last row
    curr -= 3; // step back so that curr points to the second to last row

    //-------------------------------------------------------------------------
    // Perform the backward solve with the U factor of the LU factorization
    //-------------------------------------------------------------------------
    x[n - 1] /= bn;

    for(Index i = 2; i <= n; ++i, curr -= 3)
    {
        const auto& k = n - i; // the index of the current row
        const auto& b = curr[1]; // `b` value on the current row
        const auto& c = curr[2]; // `c` value on the current row

        x[k] = (x[k] - c * x[k + 1])/b;
    }
}

auto TridiagonalMatrix::solve(VectorRef x) const -> void
{
    solve(x, x);
}

TridiagonalMatrix::operator Matrix() const
{
    const Index n = size();
    Matrix res = zeros(n, n);
    res.row(0).head(2) = row(0).tail(2);
    for(Index i = 1; i < n - 1; ++i)
        res.row(i).segment(i - 1, 3) = row(i);
    res.row(n - 1).tail(2) = row(n - 1).head(2);
    return res;
}

Mesh::Mesh()
{}

Mesh::Mesh(Index num_cells, double xl, double xr)
{
    setDiscretization(num_cells, xl, xr);
}

auto Mesh::setDiscretization(Index num_cells, double xl, double xr) -> void
{
    Assert(xr > xl, "Could not set the discretization.",
        "The x-coordinate of the right boundary needs to be "
        "larger than that of the left boundary.");

    m_num_cells = num_cells;
    m_xl = xl;
    m_xr = xr;
    m_dx = (xr - xl) / num_cells;
    m_xcells = linspace(num_cells, xl + 0.5*m_dx, xr - 0.5*m_dx);
}

TransportSolver::TransportSolver()
{
}

//auto TransportSolver::step(VectorRef u, VectorConstRef q) -> void
//{
//    // TODO: Implement Kurganov-Tadmor method as detailed in their 2000 paper (not as in Wikipedia)
//    const double dx = mmesh.dx();
//    const double num_cells = mmesh.numCells();
//    const double alpha = velocity*dt/dx;
//    const double beta = diffusion*dt/(dx * dx);
//
//    u0 = u;
//    A.resize(num_cells);
//    phi.resize(num_cells);
//
//    // Calculate the flux limiters in the interior cells
//    for(Index i = 1; i < num_cells - 1; ++i)
//    {
//        // Calculate the variation index `r = (uP - uW)/(uE - uP)` on current cell
//        const double r = (u[i] - u[i - 1])/(u[i + 1] - u[i]);
//
//        // Calculate the flux limiter phi based on the superbee limiter (https://en.wikipedia.org/wiki/Flux_limiter)
//        phi[i] = std::max(0.0, std::max(std::min(2 * r, 1.0), std::min(r, 2.0)));
//    }
//
//    // Assemble the coefficient matrix A for the interior cells
//    for(Index icell = 1; icell < num_cells - 1; ++icell)
//    {
//        const double phiW = phi[icell - 1];
//        const double phiP = phi[icell];
//        const double aux = 1.0 + 0.5 * (phiP - phiW);
//        const double a = -beta;
//        const double b = 1 + 2*beta;
//        const double c = -beta;
//        A.row(icell) << a, b, c;
//
//        const double uW = u0[icell - 1];
//        const double uP = u0[icell];
//        u[icell] += aux*alpha * (uW - uP);
//    }
//
//    // Assemble the coefficient matrix A for the boundary cells
//    A.row(0) << 0.0, 1.0 + beta, -beta;
//    A.row(num_cells - 1) << -beta, 1.0 + beta, 0.0;
//
//    u[0] += alpha * (ul - u0[0]);
//    u[num_cells - 1] += alpha * (u0[num_cells - 2] - u0[num_cells - 1]);
//
//    u += dt * q;
//
//    A.factorize();
//    A.solve(u);
//}

auto TransportSolver::initializeFullImplicit() -> void
{
    const auto dx = mesh_.dx();
    const auto alpha = diffusion*dt/(dx * dx);
    const auto beta = velocity*dt/dx;
    const auto num_cells = mesh_.numCells();
    const auto icell0 = 0;
    const auto icelln = num_cells - 1;

    A.resize(num_cells);

    // Assemble the coefficient matrix A for the interior cells
    for(Index icell = 1; icell < icelln; ++icell)
    {
        const double a = -beta - alpha;
        const double b = 1 + beta + 2*alpha;
        const double c = -alpha;
        A.row(icell) << a, b, c;
    }

    // Assemble the coefficient matrix A for the boundary cells
    A.row(icell0) << 0.0, 1.0 + alpha + beta, -alpha;   // left boundary (flux = v * ul)
    A.row(icelln) << - beta, 1.0 + beta, 0.0;           // right boundary (free)

    // Factorize A into LU factors for future uses in method step
    A.factorize();

}
auto TransportSolver::step_implicit_fvm(Reaktoro::VectorRef u,
                                        Reaktoro::VectorConstRef q) -> void
{
    // Solving advection problem with time explicit approach
    const auto dx = mesh_.dx();
    const auto beta = velocity*dt/dx;
    const auto icell0 = 0;

    // Handle the left boundary cell
    u[icell0] += beta * ul; // left boundary (flux = v * ul)

    // Add the source contribution
    u += dt * q;

    // Solving the diffusion problem with time implicit approach
    A.solve(u);
}


auto TransportSolver::initialize() -> void
{
    const auto dx = mesh_.dx();
    const auto alpha = diffusion*dt/(dx * dx);
    const auto num_cells = mesh_.numCells();
    const auto icell0 = 0;
    const auto icelln = num_cells - 1;

    Assert(velocity * dt / dx < 0.5,
           "Could not run reactive-transport calculation reliably.",
           "The CFL number ( = velocity * dt / dx ) must be less then 0.5. "
           "Try to decrease the time step or coarsen the spatial discretization"
           "(increase the number of cells)");

    A.resize(num_cells);
    phi.resize(num_cells);

    // Assemble the coefficient matrix A for the interior cells
    for(Index icell = 1; icell < icelln; ++icell)
    {
        const double a = -alpha;
        const double b = 1 + 2*alpha;
        const double c = -alpha;
        A.row(icell) << a, b, c;
    }

    // Assemble the coefficient matrix A for the boundary cells
    // Our derivation
    // A.row(icell0) << 0.0, 1.0 + alpha, -alpha;            // forward difference approximation with first order error
    // A.row(icelln) << 0.0, 1.0, 0.0;                       // d/dx = 0 (zero flux) at the right boundary

    // Allan version ealier
    // A.row(icell0) << 0.0, 1.0 + alpha, -alpha;            // forward difference approximation with second order error
    // A.row(icelln) << -alpha, 1.0 + alpha, 0.0;            // d/dx = 0 (zero flux) at the right boundary

    // ESSS
    A.row(icell0) << 0.0, 1.0 + 4.5*alpha, -1.5*alpha;    // forward difference approximation with second order error
    A.row(icelln) << -alpha, 1.0 + alpha, 0.0;            // d/dx = 0 (zero flux) at the right boundary

    // Factorize A into LU factors for future uses in method step
    A.factorize();
}

auto TransportSolver::step(VectorRef u, VectorConstRef q) -> void
{
    // TODO: Implement Kurganov-Tadmor method as detailed in their 2000 paper (not as in Wikipedia)
    // Solving advection problem with time explicit approach
    const auto dx = mesh_.dx();
    const auto num_cells = mesh_.numCells();
    const auto beta = velocity*dt/dx;
    const auto icell0 = 0;
    const auto icelln = num_cells - 1;

    u0 = u;

    phi[0] = 2.0; //  this is very important to ensure correct flux limiting behavior for boundary cell.

    // Calculate the flux limiters in the interior cells
    for(Index icell = 1; icell < icelln; ++icell)
    {
        // Calculate the variation index `r = (uP - uW)/(uE - uP)` on current cell
        const double r = (u[icell] - u[icell - 1])/(u[icell + 1] - u[icell]);

        // Calculate the flux limiter phi based on the superbee limiter (https://en.wikipedia.org/wiki/Flux_limiter)
        phi[icell] = std::max(0.0, std::max(std::min(2 * r, 1.0), std::min(r, 2.0)));
    }

    // Compute advection contributions to u for the interior cells
    for(Index icell = 1; icell < icelln; ++icell)
    {
        const double phiW = phi[icell - 1];
        const double phiP = phi[icell];
        const double aux = 1.0 + 0.5 * (phiP - phiW);

        const double uW = u0[icell - 1];
        const double uP = u0[icell];
        u[icell] += aux*beta * (uW - uP);
    }

    // Handle the left boundary cell
    const double aux = 1 + 0.5 * phi[0];
    u[icell0] += aux * beta * (ul - u0[0]) + (3.0*diffusion*ul*dt/(dx*dx)); // prescribed amount on the wall and approximation derived by forward difference approximation with second order error
    // Handle the right boundary cell
    u[icelln] += beta * (u0[icelln - 1] - u0[icelln]); // du/dx = 0 at the right boundary

    // Add the source contribution
    u += dt * q;

    // Solving the diffusion problem with time implicit approach
    A.solve(u);
}

auto TransportSolver::step(VectorRef u) -> void
{
    //step(u, zeros(u.size()));             // flux-limiters scheme
    step_implicit_fvm(u, zeros(u.size()));  // implicit fvl scheme
}

} // namespace Reaktoro

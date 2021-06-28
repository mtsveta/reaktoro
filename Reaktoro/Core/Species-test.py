# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright (C) 2014-2021 Allan Leal
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library. If not, see <http://www.gnu.org/licenses/>.

from reaktoro import *
import pytest

def testSpecies():
    # species = Species("CO3--").withName("CO3--(aq)").withTags(["aqueous", "anion", "charged"])
    species = Species("CO3--")
    species = species.withName("CO3--(aq)")
    species = species.withTags(["aqueous", "anion", "charged"])  # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< a bad alloc error happens here
    assert species.name() == "CO3--(aq)"
    print(species.formula())
    assert species.formula() == "CO3--"
    assert species.substance() == "CO3--"
    assert species.charge() == -2
    assert species.molarMass() == pytest.approx(0.0600102972)
    assert species.aggregateState() == AggregateState.Aqueous
    assert species.elements().size() == 2
    assert species.elements().coefficient("C") == 1
    assert species.elements().coefficient("O") == 3
    assert species.tags().size() == 3
    assert "aqueous" in species.tags()
    assert "anion" in species.tags()
    assert "charged" in species.tags()
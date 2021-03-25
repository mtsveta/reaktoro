// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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

// C++ includes
#include <istream>
#include <string>

// yaml-cpp includes
#include <yaml-cpp/yaml.h>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>

namespace Reaktoro {

/// A class used to serialize/deserialize other types based on yaml format.
class yaml : public YAML::Node
{
public:
    /// Construct an object of type yaml.
    yaml();

    /// Construct an object of type yaml from a given input string.
    yaml(const char* input);

    /// Construct an object of type yaml from a given input string.
    yaml(const std::string& input);

    /// Construct an object of type yaml from a given input stream.
    yaml(std::istream& input);

    /// Construct an object of type yaml from a given Node object.
    yaml(const YAML::Node& node);

    /// Append a child node with a given value only if value is not default value.
    template<typename T, typename U = T>
    auto appendIfNotDefault(const std::string& key, const T& value, const U& defaultval = U{}) -> void {
        if(value != defaultval) (*this)[key] = value; }

    /// Return a child node with given key if found, otherwise raise an error.
    auto at(const std::string& key) const -> yaml;

    /// Return the value of the node if not empty, otherwise a given fallback value.
    template<typename T>
    auto value(const T& fallback = T{}) const -> T { return IsDefined() ? as<T>() : fallback; }

    /// Return child with given key.
    template<typename Key>
    auto operator[](const Key& key) const -> const yaml { return YAML::Node::operator[](key); }

    /// Return child with given key.
    template<typename Key>
    auto operator[](const Key& key) -> yaml { return YAML::Node::operator[](key); }

    /// Assign this yaml node with given value.
    template<typename T>
    auto operator=(const T& value) { YAML::Node::operator=(value); return *this; }

    /// Implicitly convert this yaml object into another type.
    template<typename T>
    operator T() const { return as<T>(); }

    /// Transfer the value at this yaml node to argument @p value.
    template<typename T>
    auto to(T& value) -> void
    {
        try {
            value = as<T>();
        }
        catch(...) {
            errorif(true, "Could not convert YAML node to requested value type.");
        }
    }
};

} // namespace Reaktoro

namespace YAML {

template<typename Type>
struct convert
{
    static auto encode(const Type& obj)
    {
        Reaktoro::yaml node;
        node << obj;
        return node;
    }

    static auto decode(const Node& node, Type& obj)
    {
        Reaktoro::yaml(node) >> obj;
        return true;
    }
};

} // namespace YAML

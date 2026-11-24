// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006-2026 Tiago de Paula Peixoto <tiago@skewed.de>
//
// This program is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 3 of the License, or (at your option) any
// later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#ifndef HASH_MAP_WRAP_HH
#define HASH_MAP_WRAP_HH

#include <boost/unordered/unordered_flat_set.hpp>
#include <boost/unordered/unordered_flat_map.hpp>

template<class Key,
         class Hash = boost::hash<Key>,
         class Pred = std::equal_to<Key>,
         class Allocator = std::allocator<Key>>
using gt_hash_set = boost::unordered_flat_set<Key, Hash, Pred, Allocator>;

template<class Key,
         class T,
         class Hash = boost::hash<Key>,
         class Pred = std::equal_to<Key>,
         class Allocator = std::allocator<std::pair<const Key, T>>>
using gt_hash_map = boost::unordered_flat_map<Key, T, Hash, Pred, Allocator>;

#endif // HASH_MAP_WRAP_HH

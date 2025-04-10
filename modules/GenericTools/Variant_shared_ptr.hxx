//! Copyright : see license.txt
//!
//! \brief

#pragma once

#include "../GenericMerope/StdHeaders.hxx"

namespace merope {
template<class ...A>
using variant_shared_ptr = std::variant<std::shared_ptr<A>...>;
}  // namespace  merope


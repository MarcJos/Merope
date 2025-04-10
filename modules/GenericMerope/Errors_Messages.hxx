//! Copyright : see license.txt
//!
//! \brief Error messages for Merope

#pragma once

namespace merope {


//! @brief: Merope assert. Always activated, not only on debug
#define Merope_assert(test, message)                                            \
{                                                                               \
    if (not (test)) {                                                           \
        cerr << __PRETTY_FUNCTION__ << endl;                                    \
        std::string mes_warning_macro = message;                                \
        cerr << "\n\n\n" << "Merope Warning :" << mes_warning_macro << "\n\n\n" << endl;  \
        throw runtime_error(mes_warning_macro);                                 \
    }                                                                           \
}

#define Merope_warning(test, message)                                            \
{                                                                               \
    if (not (test)) {                                                           \
        cerr << __PRETTY_FUNCTION__ << endl;                                    \
        std::string mes_warning_macro = message;                                \
        cerr << "\n\n\n" << "Merope Warning :" << mes_warning_macro << "\n\n\n" << endl;  \
}   \
}

#define Merope_error_not_done()\
Merope_assert(false, "Not programmed yet")

#define Merope_error_impossible()\
Merope_assert(false, "Impossible")

template <class...>
struct False : std::bool_constant<false> {};

//! @brief: Triggers a static error.
//! @brief: Is like static_assert(false, message)
#define Merope_static_error(MY_TYPE, message) static_assert(merope::False<MY_TYPE>{}, message)

}  // namespace  merope

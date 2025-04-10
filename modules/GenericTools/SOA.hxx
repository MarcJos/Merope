//! Copyright : see license.txt
//!
//! \briefImplements a structure of arrays
//!

#include "../../GenericMerope/StdHeaders.hxx"

template<class... T>
class SOA : public std::tuple<std::vector<T>...> {
public:
    static constexpr size_t nb_S() { return std::tuple_size_v<std::tuple<std::vector<T>...>>; }

    template<class FUNCTOR>
    void apply_on_all(FUNCTOR f) {
        apply_on_all_I<nb_S() - 1, FUNCTOR>(f);
    }

    template<class FUNCTOR>
    void apply_on_all(FUNCTOR f) const {
        apply_on_all_I<nb_S() - 1, FUNCTOR>(f);
    }

    template<class T_0, size_t I_0 = 0>
    const std::vector<T_0>& get() const {
        const auto& res = std::get<I_0>(*this);
        if constexpr (std::is_same_v<const T_0&, const typeof(res[0])&>) {
            return res;
        } else {
            return get<T_0, I_0 + 1>();
        }
    }

    template<class T_0, size_t I_0 = 0>
    std::vector<T_0>& get() {
        auto& res = std::get<I_0>(*this);
        if constexpr (std::is_same_v<const T_0&, const typeof(res[0])&>) {
            return res;
        } else {
            return get<T_0, I_0 + 1>();
        }
    }

    template<size_t I = 0>
    size_t get_total_size() const {
        if constexpr (I == nb_S()) {
            return 0.;
        } else {
            return std::get<I>(*this).size() + get_total_size<I + 1>();
        }
    }

private:
    template<size_t I, class FUNCTOR>
    void apply_on_all_I(FUNCTOR f) const {
        f(std::get<I>(*this));
        if constexpr (I > 0) {
            apply_on_all_I<I - 1, FUNCTOR>(f);
        }
    }

    template<size_t I, class FUNCTOR>
    void apply_on_all_I(FUNCTOR f) {
        f(std::get<I>(*this));
        if constexpr (I > 0) {
            apply_on_all_I<I - 1, FUNCTOR>(f);
        }
    }
};

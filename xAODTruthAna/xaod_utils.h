#ifndef XAOD_UTILS_H
#define XAOD_UTILS_H

#include <boost/throw_exception.hpp>

#define CHECK_THROW(_Expr) \
    do { \
        bool result_M_ = _Expr; \
        if (!result_M_) { \
            BOOST_THROW_EXCEPTION(std::runtime_error(std::string(#_Expr) + " failed!")); \
        } \
    } while (0)

template <typename T, typename Store>
const T* Retrieve(Store* st, const char *name)
{
    const T* result = nullptr;
    CHECK_THROW(st->retrieve(result, name));
    return result;
}

template <typename T, typename Store>
const T* Retrieve(Store* st, const std::string& name)
{
    return Retrieve<T>(st, name.c_str());
}

#endif

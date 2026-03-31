#ifndef _ASD_STD_TUPLE_CONCAT
#define _ASD_STD_TUPLE_CONCAT

#include <tuple>
#include <iostream>
#include <type_traits>

template <typename... Ts>
struct std_tuple_concat;

template <typename... Ts1, typename... Ts2, typename... Rest>
struct std_tuple_concat<std::tuple<Ts1...>, std::tuple<Ts2...>, Rest...> {
    using type = typename std_tuple_concat<std::tuple<Ts1..., Ts2...>, Rest...>::type;
};

template <typename... Ts>
struct std_tuple_concat<std::tuple<Ts...>> {
    using type = std::tuple<Ts...>;
};

template <typename... Ts>
using std_tuple_concat_Type = typename std_tuple_concat<Ts...>::type;

#endif
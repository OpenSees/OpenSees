#ifndef _ASD_UTUPLE_STORAGE
#define _ASD_UTUPLE_STORAGE


#include <iostream>
#include <tuple>
/**
 * utuple_index:
 * A helper struct to get, at compile time, the first index
 * of a type 'T' in a tuple 'Tuple'.
 */
// template <class T, class Tuple>
// struct utuple_index;

// template <class T, class... Types>
// struct utuple_index<T, std::tuple<T, Types...>> {
//     static const std::size_t value = 0;
// };

// template <class T, class U, class... Types>
// struct utuple_index<T, std::tuple<U, Types...>> {
//     static const std::size_t value = 1 + utuple_index<T, std::tuple<Types...>>::value;
// };


// utuple indexing helper struct

template <class T, class Tuple, bool NotFound = false>
struct utuple_index;

template <class T, class... Types>
struct utuple_index<T, std::tuple<T, Types...>, false> {
    static const std::size_t value = 0;
};

template <class T, class U, class... Types>
struct utuple_index<T, std::tuple<U, Types...>, false> {
    static const std::size_t value = 1 + utuple_index<T, std::tuple<Types...>>::value;
};

template <class T>
struct utuple_index<T, std::tuple<>, false> {
    static_assert(sizeof(T) == -1, "\n\n\n\n"
"===========================================================\n"
"utuple_storage.h -- Type T not found in the utuple_storage \n"
"This error means that while searching for a particular \n"
"type within a utuple_storage, the particular type was \n"
"not found. Check whether the utuple_storage was built \n"
"using the correct types and whether you are searching \n"
"for a type present within the utuple_storage, or maybe \n"
"you're not using the correct utuple_storage. \n"
"===========================================================\n\n\n\n");
};






// for calling commit and revert

// Helper function to call commit on a single element
template <class T>
void call_commit(T& t) {
    t.commit();
}

// Helper function to call revert on a single element
template <class T>
void call_revert(T& t) {
    t.revert();
}

// Recursive function to call commit on all elements in the tuple
template <std::size_t I = 0, typename... T>
typename std::enable_if<I < sizeof...(T), void>::type
call_commit_on_tuple_elements(std::tuple<T...>& tuple) {
    call_commit(std::get<I>(tuple));
    call_commit_on_tuple_elements<I + 1, T...>(tuple);
}

template <std::size_t I = 0, typename... T>
typename std::enable_if<I == sizeof...(T), void>::type
call_commit_on_tuple_elements(std::tuple<T...>&) {}

// Recursive function to call revert on all elements in the tuple
template <std::size_t I = 0, typename... T>
typename std::enable_if<I < sizeof...(T), void>::type
call_revert_on_tuple_elements(std::tuple<T...>& tuple) {
    call_revert(std::get<I>(tuple));
    call_revert_on_tuple_elements<I + 1, T...>(tuple);
}

template <std::size_t I = 0, typename... T>
typename std::enable_if<I == sizeof...(T), void>::type
call_revert_on_tuple_elements(std::tuple<T...>&) {}



// For getting parameters types

template <typename Tuple1, typename Tuple2>
struct tuple_cat_helper;

template <typename... Ts1, typename... Ts2>
struct tuple_cat_helper<std::tuple<Ts1...>, std::tuple<Ts2...>> {
    using type = std::tuple<Ts1..., Ts2...>;
};

template <typename Tuple1, typename Tuple2>
using tuple_cat_t = typename tuple_cat_helper<Tuple1, Tuple2>::type;


// For operating on utuple_storage elements

template <std::size_t I = 0, typename Function, typename... Tp>
inline typename std::enable_if<I == sizeof...(Tp), void>::type
apply_to_each_in_tuple(std::tuple<Tp...>&, Function) { }

template <std::size_t I = 0, typename Function, typename... Tp>
inline typename std::enable_if<I < sizeof...(Tp), void>::type
apply_to_each_in_tuple(std::tuple<Tp...>& t, Function f) {
    f(std::get<I>(t));
    apply_to_each_in_tuple<I + 1, Function, Tp...>(t, f);
}




/**
 * utuple_storage:
 * A heterogeneous container wrapping a std::tuple, indexable by type.
 * Only the first occurrence of a type is considered.
 */
template<typename tuple_t>
class utuple_storage {
public:
    // Sets the value of the element with type T in the tuple.
    template<class T>
    void set(const T& value) {
        std::get< utuple_index<T, tuple_t>::value >(data) = value;
    }

    // Retrieves the value of the element with type T in the tuple.
    template<class T>
    const T& get() const {
        return std::get< utuple_index<T, tuple_t>::value >(data);
    }

    // Returns the size of the tuple.
    inline static constexpr std::size_t size() {
        return std::tuple_size<tuple_t>::value;
    }

    void commit_all() {
        call_commit_on_tuple_elements(data);
    }

    void revert_all() {
        call_revert_on_tuple_elements(data);
    }

    // Prints the index, type, name, and value of all components in the utuple_storage.
    void print_components() const {
        print_components_impl<0>();
    }

    // New method to concatenate parameters_t tuples.
    auto concatenate_parameters() const {
        return concatenate_parameters_impl<0>();
    }

    //Apply a function to each element of the utuple_storage
	template<typename Function>
	void apply(Function f) {
	    apply_to_each_in_tuple(data, f);
	}

    // The wrapped std::tuple
    tuple_t data;
private:

    template <std::size_t I, typename std::enable_if<I < std::tuple_size<tuple_t>::value, int>::type = 0>
    void print_components_impl() const {
        using current_type = typename std::tuple_element<I, tuple_t>::type;
        // constexpr auto type_name = ct_type_name<current_type>::value();
        std::cout << "Index: " << I
                  // << ", Type: " << type_name
                  // << ", Name: " << get<I>().name
                  << ", Value: " << std::get<I>(data)
                  << std::endl;

        print_components_impl<I + 1>();
    }

    template <std::size_t I, typename std::enable_if<I == std::tuple_size<tuple_t>::value, int>::type = 0>
    void print_components_impl() const {
        // Base case, do nothing.
    }


    // Helper function to recursively concatenate parameters_t tuples.
    template <std::size_t I, typename std::enable_if<I < std::tuple_size<tuple_t>::value, int>::type = 0>
    auto concatenate_parameters_impl() const {
        using current_type = typename std::tuple_element<I, tuple_t>::type;
        using parameters_tuple_type = typename current_type::parameters_t;
        auto current_parameters = std::get<I>(data).template hardening_function_parameters<parameters_tuple_type>();
        auto remaining_parameters = concatenate_parameters_impl<I + 1>();
        return tuple_cat_t<parameters_tuple_type, decltype(remaining_parameters)>(std::tuple_cat(current_parameters, remaining_parameters));
    }

    template <std::size_t I, typename std::enable_if<I == std::tuple_size<tuple_t>::value, int>::type = 0>
    auto concatenate_parameters_impl() const {
        return std::tuple<>();
    }
};


namespace std {
    template<typename... Ts>
    class tuple_size<utuple_storage<Ts...>> : public std::tuple_size<std::tuple<Ts...>> {
    };
}









/**
 * utuple_concat:
 * A struct to concatenate tuples.
 */
template <typename... Ts>
struct utuple_concat;

template <typename... Ts1, typename... Ts2, typename... Rest>
struct utuple_concat<std::tuple<Ts1...>, std::tuple<Ts2...>, Rest...> {
    using type = typename utuple_concat<std::tuple<Ts1..., Ts2...>, Rest...>::type;
};

template <typename... Ts>
struct utuple_concat<std::tuple<Ts...>> {
    using type = std::tuple<Ts...>;
};

// Type alias to use the concatenated tuple type.
template <typename... Ts>
using utuple_concat_type = typename utuple_concat<Ts...>::type;


/**
 * utuple_concat_unique:
 * A struct to concatenate tuples, but make it unique.
 */
// Forward declaration of utuple_contains


// // C++17 version
// template <typename T, typename Tuple>
// struct utuple_contains;

// // Check if a tuple contains a specific type
// template <typename T, typename... Ts>
// struct utuple_contains<T, std::tuple<Ts...>> : std::bool_constant<(std::is_same_v<T, Ts> || ...)> {};


//C++14 version 
template <typename T, typename Tuple>
struct utuple_contains;

template <typename T, typename... Ts>
struct utuple_contains<T, std::tuple<T, Ts...>> : std::true_type {};

template <typename T, typename U, typename... Ts>
struct utuple_contains<T, std::tuple<U, Ts...>> : utuple_contains<T, std::tuple<Ts...>> {};

template <typename T>
struct utuple_contains<T, std::tuple<>> : std::false_type {};




// // C++17 version
// // utuple_concat_unique: A struct to concatenate tuples with unique types
// template <typename... Ts>
// struct utuple_concat_unique;

// template <typename T, typename... Ts1, typename... Ts2, typename... Rest>
// struct utuple_concat_unique<std::tuple<Ts1...>, std::tuple<T, Ts2...>, Rest...> {
//     using type = std::conditional_t<
//         utuple_contains<T, std::tuple<Ts1...>>::value,
//         typename utuple_concat_unique<std::tuple<Ts1...>, std::tuple<Ts2...>, Rest...>::type,
//         typename utuple_concat_unique<std::tuple<Ts1..., T>, std::tuple<Ts2...>, Rest...>::type>;
// };

// // Specialization for the case when the second tuple is empty
// template <typename... Ts1>
// struct utuple_concat_unique<std::tuple<Ts1...>, std::tuple<>> {
//     using type = std::tuple<Ts1...>;
// };

// template <typename... Ts>
// struct utuple_concat_unique<std::tuple<Ts...>> {
//     using type = std::tuple<Ts...>;
// };


//C++14 version
template <typename Tuple1, typename Tuple2>
struct utuple_concat_unique;

template <typename... Ts1>
struct utuple_concat_unique<std::tuple<Ts1...>, std::tuple<>> {
    using type = std::tuple<Ts1...>;
};

template <typename... Ts1, typename T, typename... Ts2>
struct utuple_concat_unique<std::tuple<Ts1...>, std::tuple<T, Ts2...>> {
    using type = typename std::conditional<
        utuple_contains<T, std::tuple<Ts1...>>::value,
        typename utuple_concat_unique<std::tuple<Ts1...>, std::tuple<Ts2...>>::type,
        typename utuple_concat_unique<std::tuple<Ts1..., T>, std::tuple<Ts2...>>::type
    >::type;
};


// Type alias to use the concatenated tuple type with unique types
template <typename... Ts>
using utuple_concat_unique_type = typename utuple_concat_unique<Ts...>::type;

#endif //_ASD_UTUPLE_STORAGE
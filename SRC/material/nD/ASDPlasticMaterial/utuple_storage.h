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
#include "std_tuple_concat.h"


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

    // The wrapped std::tuple
    tuple_t data;
    
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


    //Apply a function to each element of the utuple_storage, serves as a
    // foreach loop over stored variables
    template<typename Function>
    void apply(Function f) {
        apply_to_each_in_tuple(data, f);
    }

    //Get names of parameter types
    auto getParameterNames() const
    {
        return getParameterNamesTuple(data);
    }


    //To set a parameter by name
    void setParameterByName(const char* name, double value)  {
        setParameterByName_impl<0>(name, value);
    }

    //Get names of parameter types
    auto getInternalVariableNames() const
    {
        return getInternalVariableNamesTuple(data);
    }

    //Get an IV by name
    template<typename... Ts>
    const auto& getInternalVariableByName(const char* name) const {
        return getInternalVariableByNameImpl<0>(data, name);
    }

    //To set an internal variable by name
    int getInternalVariableSizeByName(const char* name)  const {
        return getInternalVariableSizeByName_impl<0>(name);
    }

    // Public function to get the index of an internal variable by name
    int getInternalVariableIndexByName(const char* name) const {
        return getInternalVariableIndexByNameImpl<0>(data, name);
    }

    //To set an internal variable by name
    void setInternalVariableByName(const char* name, int size, double *values)  {
        setInternalVariableByName_impl<0>(name, size, values);
    }
    
    std::string getVariableNamesAndHardeningLaws() const {
        return getVariableNamesAndHardeningLaws_impl<0>();
    }

    template<typename VarType>
    auto getTrialInternalVariable() const
    {
    	const VarType& var = this->template get<VarType> ();
        return var.trial_value;
    }

    // Assignment operator
    utuple_storage& operator=(const utuple_storage& other) {
        if (this != &other) {
            assign_from_other<0>(other.data);
        }
        return *this;
    }

private:

    // Implementations for print_components
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


    // Implementations for setParameterByName
    template <std::size_t I, typename std::enable_if<I < std::tuple_size<tuple_t>::value, int>::type = 0>
    void setParameterByName_impl(const char* name, double value)  {
        using current_type = typename std::tuple_element<I, tuple_t>::type;

        auto current_name =std::get<I>(data).getName();


        if (std::strcmp(current_name, name) == 0) {
            std::get<I>(data).value = value;
        } else {
            setParameterByName_impl<I + 1>(name, value);
        }

        setParameterByName_impl<I + 1>(name, value);
    }

    template <std::size_t I, typename std::enable_if<I == std::tuple_size<tuple_t>::value, int>::type = 0>
    void setParameterByName_impl(const char* name, double value)  {
        // Base case, do nothing.
    }

    // Implementations for getInternalVariableByNameImpl
    template <std::size_t I = 0, typename... Tp>
    typename std::enable_if<I < sizeof...(Tp), const typename std::tuple_element<I, std::tuple<Tp...>>::type&>::type
    getInternalVariableByNameImpl(const std::tuple<Tp...>& tuple, const char* name) const {
        const auto& variable = std::get<I>(tuple);
        if (strcmp(variable.getName(), name) == 0) {
            return variable;
        }
        return getInternalVariableByNameImpl<I + 1>(tuple, name);
    }

    template <std::size_t I = 0, typename... Tp>
    typename std::enable_if<I == sizeof...(Tp), const typename std::tuple_element<0, std::tuple<Tp...>>::type&>::type
    getInternalVariableByNameImpl(const std::tuple<Tp...>&, const char*) const {
        throw std::runtime_error("Variable not found");
    }

    // Implementations for getInternalVariableIndexByNameImpl
    template <std::size_t I = 0, typename... Tp>
    typename std::enable_if<I < sizeof...(Tp), int>::type
    getInternalVariableIndexByNameImpl(const std::tuple<Tp...>& tuple, const char* name) const {
        const auto& variable = std::get<I>(tuple);
        if (strcmp(variable.getName(), name) == 0) {
            return static_cast<int>(I);
        }
        return getInternalVariableIndexByNameImpl<I + 1>(tuple, name);
    }

    template <std::size_t I = 0, typename... Tp>
    typename std::enable_if<I == sizeof...(Tp), int>::type
    getInternalVariableIndexByNameImpl(const std::tuple<Tp...>&, const char*) const {
        return -1; // Return -1 if the variable is not found
    }




    // Implementations for setInternalVariableByName
    template <std::size_t I, typename std::enable_if<I < std::tuple_size<tuple_t>::value, int>::type = 0>
    int getInternalVariableSizeByName_impl(const char* name) const {
        using current_type = typename std::tuple_element<I, tuple_t>::type;

        auto current_name = std::get<I>(data).getName();

        if (std::strcmp(current_name, name) == 0) {
            return std::get<I>(data).size();
        } else {
            return getInternalVariableSizeByName_impl<I + 1>(name);
        }

        getInternalVariableSizeByName_impl<I + 1>(name);
    }

    template <std::size_t I, typename std::enable_if<I == std::tuple_size<tuple_t>::value, int>::type = 0>
    int getInternalVariableSizeByName_impl(const char* name) const  {
        // Base case, do nothing.
        return -1;
    }



    // Implementations for setInternalVariableByName
    template <std::size_t I, typename std::enable_if<I < std::tuple_size<tuple_t>::value, int>::type = 0>
    void setInternalVariableByName_impl(const char* name, int size, double * values)  {
        using current_type = typename std::tuple_element<I, tuple_t>::type;

        auto current_name =std::get<I>(data).getName();


        if (std::strcmp(current_name, name) == 0) {
            for (int i = 0; i < size; ++i)
            {
                std::get<I>(data).trial_value[i] = values[i];
                std::get<I>(data).committed_value[i] = values[i];
            }
        } else {
            setInternalVariableByName_impl<I + 1>(name, size, values);
        }

        setInternalVariableByName_impl<I + 1>(name, size, values);
    }

    template <std::size_t I, typename std::enable_if<I == std::tuple_size<tuple_t>::value, int>::type = 0>
    void setInternalVariableByName_impl(const char* name, int size, double *values)  {
        // Base case, do nothing.
    }


    template <std::size_t I, typename std::enable_if<I < std::tuple_size<tuple_t>::value, int>::type = 0>
    std::string getVariableNamesAndHardeningLaws_impl() const {
        using current_type = typename std::tuple_element<I, tuple_t>::type;

        // Get variable name
        std::string var_name = std::get<I>(data).getFullName();

        // cout << "var_name = " << var_name << endl;

        // Get hardening law type
        // std::string law_name = typeid(typename current_type::parameters_t::HardeningPolicy).name();
        // std::string law_name = current_type.;

        return var_name + ":" + getVariableNamesAndHardeningLaws_impl<I + 1>();
    }

    template <std::size_t I, typename std::enable_if<I == std::tuple_size<tuple_t>::value, int>::type = 0>
    std::string getVariableNamesAndHardeningLaws_impl() const {
        // Base case, do nothing.
        return "";
    }

    // Helper template function for assignment
    template <std::size_t I, typename std::enable_if<I < std::tuple_size<tuple_t>::value, int>::type = 0>
    void assign_from_other(const tuple_t& other) {
        std::get<I>(data) = std::get<I>(other);
        assign_from_other<I + 1>(other);
    }

    template <std::size_t I, typename std::enable_if<I == std::tuple_size<tuple_t>::value, int>::type = 0>
    void assign_from_other(const tuple_t&) {
        // Base case, do nothing.
    }
};  //end utuple_storage


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


//C++14 version 
template <typename T, typename Tuple>
struct utuple_contains;

template <typename T, typename... Ts>
struct utuple_contains<T, std::tuple<T, Ts...>> : std::true_type {};

template <typename T, typename U, typename... Ts>
struct utuple_contains<T, std::tuple<U, Ts...>> : utuple_contains<T, std::tuple<Ts...>> {};

template <typename T>
struct utuple_contains<T, std::tuple<>> : std::false_type {};




// utuple_concat_unique: A struct to concatenate tuples with unique types
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



// For parameter extraction
template <typename T>
struct ExtractNestedParameterTypes;

template <typename... Ts>
struct ExtractNestedParameterTypes<std::tuple<Ts...>> {
    using type = std_tuple_concat_Type<typename Ts::parameters_t...>;
};

template <typename T>
using ExtractNestedParameterTypes_t = typename ExtractNestedParameterTypes<T>::type;


// To get parameter names, for initialization
template<std::size_t... Is, typename... Ts>
auto getParameterNamesTuple(std::index_sequence<Is...>, const std::tuple<Ts...>& params) {
    return std::make_tuple(std::get<Is>(params).getName()...);
}

template<typename... Ts>
auto getParameterNamesTuple(const std::tuple<Ts...>& params) {
    return getParameterNamesTuple(std::make_index_sequence<sizeof...(Ts)>(), params);
}


// To get internal_variables names, for initialization
template<std::size_t... Is, typename... Ts>
auto getInternalVariableNamesTuple(std::index_sequence<Is...>, const std::tuple<Ts...>& ivs) {
    return std::make_tuple(std::get<Is>(ivs).getName()...);
}

template<typename... Ts>
auto getInternalVariableNamesTuple(const std::tuple<Ts...>& ivs) {
    return getInternalVariableNamesTuple(std::make_index_sequence<sizeof...(Ts)>(), ivs);
}




//Foreach


template<std::size_t I = 0, typename FuncT, typename... Tp>
inline typename std::enable_if<I == sizeof...(Tp), void>::type
for_each_in_tuple(std::tuple<Tp...> &, FuncT) // Unused arguments are ignored
{ }

template<std::size_t I = 0, typename FuncT, typename... Tp>
inline typename std::enable_if<I < sizeof...(Tp), void>::type
for_each_in_tuple(std::tuple<Tp...>& t, FuncT f)
{
    f(std::get<I>(t));
    for_each_in_tuple<I + 1, FuncT, Tp...>(t, f);
}



#endif //_ASD_UTUPLE_STORAGE
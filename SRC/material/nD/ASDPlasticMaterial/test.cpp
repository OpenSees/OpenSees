/******************************************************************************

                              Online C++ Compiler.
               Code, Compile, Run and Debug C++ program online.
Write your code in this editor and press "Run" button to compile and execute it.

*******************************************************************************/

#include <iostream>
#include <tuple>
#include <type_traits>
#include <array>

class NDMaterial
{
	NDMaterial(int tag): i(tag) {}
public:
	int i;
};

namespace {
	/**
	utuple_index:
	A helper function to get, at compile time, the first index
	of a type 'T' in a tuple 'Tuple'
	*/
	template <class T, class Tuple>
	struct utuple_index;
	template <class T, class... Types>
	struct utuple_index<T, std::tuple<T, Types...>> {
		static const std::size_t value = 0;
	};
	template <class T, class U, class... Types>
	struct utuple_index<T, std::tuple<U, Types...>> {
		static const std::size_t value = 1 + utuple_index<T, std::tuple<Types...>>::value;
	};
	
	/**
	utuple_storage:
	A heterogeneous container wrapping a std::tuple, indexable by type.
	only the first occurrence of a type is considered
	*/
	template<typename tuple_t>
	class utuple_storage
	{
	public:
		template<class T> 
		void set(const T& value) {
			std::get< utuple_index<T, tuple_t>::value >(data) = value;
		}
		template<class T> 
		const T& get() const {
			return std::get< utuple_index<T, tuple_t>::value >(data);
		}
		inline static constexpr std::size_t size() {
			return std::tuple_size<tuple_t>::value;
		}
	public:
		tuple_t data;
	};
	
	/**
	class to concatenate tuples 
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
    
    template <typename... Ts>
    using utuple_concat_type = typename utuple_concat<Ts...>::type;
    
    
    /**
    testing an indexed_types
    */
    
    struct EqPl {
        static constexpr const char* NAME = "EqPl";
        double value = 0.0;
        EqPl() = default;
        EqPl(double x) : value(x) {}
        inline const char* getName() const { return NAME; }
    };
    template<class T> 
    T& operator << (T& stream, const EqPl& x) {
        return stream << x.value;
    }
    
    /**
    basic  components
    */
    class VMYS {
    public:
        static constexpr const char* NAME = "VMYS";
        inline const char* getName() const { return NAME; }
        using internal_variables_t = std::tuple<double, int>;
        using parameters_t = std::tuple<yield_stress>;
    };
    
    class VMPP {
    public:
        static constexpr const char* NAME = "VMPP";
        inline const char* getName() const { return NAME; }
        using internal_variables_t = std::tuple<std::string, EqPl, std::string>;
    };
 
  
    template<class YS, class PP>
    class ASDPlastic {
    public:
        using conctat_types = utuple_concat_type<typename YS::internal_variables_t, typename PP::internal_variables_t>;
	    using storage_t = utuple_storage<conctat_types>;
	    
	    public:
	    void solve() {
	        // every func is templated and takes the storage
	        m_ys.eval(m_storage);
	    }
	    
		NDMaterial* clone(parameters_t) {
			// here we make the actual clone
		}
		
	private:
	    storage_t m_storage;
	    YS m_ys;
	    PP m_pp;
        
    };
  
  
    using proto_map_key = std::tuple<
        const char* /* yield surface*/, 
        const char* /* plastic potential*/>;
    using proto_map_value NDMaterial*;
    using proto_map_t = std::map<proto_map_key, proto_map_value>;
    inline make_prototypes() {
        proto_map_t m;
        
        { // a proto of VM + VM
	        auto iproto = new ASDPlastic<VMYS, VMPP>();
	        m[std::make_tuple(VMYS::NAME, VMPP::NAME)] = iprotp;
        }
        
        return m;
    }
    static proto_map_t TheAvailableProtos = make_prototypes();
    
 
}

int main()
{


	//using storage_t = utuple_storage<double, int, std::string, double, std::string>;
	
	using conctat_types = utuple_concat_type<VMYS::internal_variables_t, VMPP::internal_variables_t>;
	using storage_t = utuple_storage<conctat_types>;
	
	storage_t S;
	
	S.set(EqPl(3.2));
	S.set(1.2);
	std::cout << "get value: " << S.get<double>() << "\n";
	
	std::cout << "\nTuple contents:\n";
	std::cout << "[0] : " << std::get<0>(S.data) << "\n";
	std::cout << "[1] : " << std::get<1>(S.data) << "\n";
	std::cout << "[2] : " << std::get<2>(S.data) << "\n";
	std::cout << "[3] : " << std::get<3>(S.data) << " " << std::get<3>(S.data).getName() << "\n";
	std::cout << "[4] : " << std::get<4>(S.data) << "\n";
	std::cout << "SIZE = " << S.size() << "\n";
	
	
	
	return 0;
}
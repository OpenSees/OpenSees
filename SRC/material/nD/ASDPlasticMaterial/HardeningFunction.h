
#define HARDENING_FUNCTION_DEFINITION template <class EVT, class ParameterStorageType> \
    static auto f( \
        const EVT& current_value, \
        const VoigtVector& depsilon, \
        const VoigtVector& m, \
        const VoigtVector& sigma, \
        const ParameterStorageType& parameters_storage) 


// Function wrapper base class
template <typename EvolvingVariableType, class HardeningPolicy>
struct HardeningFunction {
	HARDENING_FUNCTION_DEFINITION 
    {
        return HardeningPolicy::f(current_value, depsilon, m, sigma, parameters_storage);
    }
    static constexpr const char* NAME = HardeningPolicy::NAME;
    using parameters_t = typename HardeningPolicy::parameters_t;
};


template <typename EvolvingVariableType, class HardeningPolicy>
std::ostream& operator<<(std::ostream& os, const HardeningFunction<EvolvingVariableType, HardeningPolicy>& obj) {
    os << "HardeningFunction<" << typeid(EvolvingVariableType).name() << ", " << typeid(HardeningPolicy).name() << ">";
    return os;
}



#define GET_PARAMETER_VALUE(type) parameters_storage.template get<type> ().value

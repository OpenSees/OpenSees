
//Base template struct for all internal variables
template <class EvolvingVariableType, class HardeningType, class NAMER>
struct InternalVariableType {
    static constexpr const char* NAME = NAMER::name;
    EvolvingVariableType trial_value;
    EvolvingVariableType committed_value;
    InternalVariableType() = default;
    InternalVariableType(EvolvingVariableType x) : trial_value(x), committed_value(x) {}
    void commit() {committed_value = trial_value;};
    void revert() {trial_value = committed_value;};
    const char *  getName() const { return NAME;}//(std::string(NAME) + std::string("(") + std::string(HardeningType::NAME) + std::string(")")).c_str(); }
    std::string getFullName() const 
    { 
        std::string fullname = std::string(NAME) + std::string("(") + std::string(HardeningType::NAME) + std::string(")");
        // cout << "fullname = " << fullname << endl;
        return fullname; 
    }
    int size() const {return trial_value.size();}

    template <class ParameterStorageType>
    EvolvingVariableType hardening_function(
                const VoigtVector &depsilon,
                const VoigtVector &m,
                const VoigtVector& sigma,
                const ParameterStorageType& parameters) const
    {
        return HardeningType::f(trial_value, depsilon, m, sigma, parameters);
    }

    // Assignment operator
    InternalVariableType& operator=(const InternalVariableType& other) {
        if (this != &other) { // Check for self-assignment
            trial_value = other.trial_value;
            committed_value = other.committed_value;
        }
        return *this;
    }

    using parameters_t = typename HardeningType::parameters_t;


};

//Stream operator for internal variables
template <class EvolvingVariableType, class HardeningType, class NAMER>
std::ostream& operator<<(std::ostream& os, const InternalVariableType<EvolvingVariableType, HardeningType, NAMER>& param) {
    os << param.getFullName() << 
    "\n               : [Trial: " << param.trial_value.transpose() << ", Committed: " << param.committed_value.transpose() << "]";
    os << endl;
    // os << "   HardeningType::parameters --> " << typeid(typename InternalVariableType<EvolvingVariableType, HardeningType, NAMER>::parameters_t).name() << endl;
    return os;
}

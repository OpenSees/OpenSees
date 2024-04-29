// ============================================================================
// Model Parameters base template struct
template<typename T, class NAMER>
struct ModelParameterType {

    static constexpr const char* NAME = NAMER::name;
    T value;
    ModelParameterType() = default;
    ModelParameterType(T x) : value(x) {}
    inline const char* getName() const { return NAME; }
};

//Templated stream insertion
template<typename T, class NAMER>
std::ostream& operator<<(std::ostream& os, const ModelParameterType<T, NAMER>& param) {
    os << param.NAME << " = " << param.value;
    return os;
}

#include <CFSTsteel.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Vector.h>
#include <StandardStream.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#ifdef CFSTSTEEL_AUTHOR_BASELINE
#include "OriginalCFSTsteel.h"
#endif

StandardStream testError;
OPS_Stream* opserrPtr = &testError;
extern void* OPS_CFSTsteel();

namespace {
std::vector<std::string> input;
unsigned position = 0;
int checks = 0;
void require(bool ok, const char* message) {
    ++checks;
    if (!ok) { std::cerr << "FAIL: " << message << '\n'; std::exit(1); }
}
bool near(double a, double b) {
    return std::isfinite(a) && std::isfinite(b) &&
        std::abs(a-b) <= 1.e-10 * (1.0 + std::max(std::abs(a), std::abs(b)));
}
void equal(UniaxialMaterial& a, UniaxialMaterial& b) {
    require(near(a.getStress(), b.getStress()), "stress continuation");
    require(near(a.getTangent(), b.getTangent()), "tangent continuation");
    require(near(a.getStrain(), b.getStrain()), "strain continuation");
}
std::vector<std::string> parameters(int n = 1) {

    std::vector<std::string> args = {"1", "206000", "305", "0.02", "-342.07",
        "26.58", "-74.87", "3.16", "4.07", "6.38", "3.12", "2.36", "5.82",
        "1", "1", "1", "1", "1", "116.39", "11.81", std::to_string(n)};
    for (int i = 0; i < n; ++i) {
        args.push_back(std::to_string(4467.68 / std::max(1, n)));
        args.push_back("38.39");
    }
    return args;
}
UniaxialMaterial* parse(const std::vector<std::string>& args) {
    input = args;
    position = 0;
    return static_cast<UniaxialMaterial*>(OPS_CFSTsteel());
}
class MemoryChannel : public Channel {
public:
    Vector saved;
    bool failSend = false, failReceive = false;
    char* addToProgram() override { return nullptr; }
    int setUpConnection() override { return 0; }
    int setNextAddress(const ChannelAddress&) override { return 0; }
    ChannelAddress* getLastSendersAddress() override { return nullptr; }
    int sendObj(int, MovableObject&, ChannelAddress*) override { return -1; }
    int recvObj(int, MovableObject&, FEM_ObjectBroker&, ChannelAddress*) override { return -1; }
    int sendMsg(int, int, const Message&, ChannelAddress*) override { return -1; }
    int recvMsg(int, int, Message&, ChannelAddress*) override { return -1; }
    int sendMatrix(int, int, const Matrix&, ChannelAddress*) override { return -1; }
    int recvMatrix(int, int, Matrix&, ChannelAddress*) override { return -1; }
    int sendID(int, int, const ID&, ChannelAddress*) override { return -1; }
    int recvID(int, int, ID&, ChannelAddress*) override { return -1; }
    int sendVector(int, int, const Vector& v, ChannelAddress*) override {
        if (failSend) return -1;
        saved = v; return 0;
    }
    int recvVector(int, int, Vector& v, ChannelAddress*) override {
        if (failReceive) return -1;
        v = saved; return 0;
    }
};
std::vector<double> history() {
    std::vector<double> result;
    double last = 0;
    for (double peak : {0.0002, -0.0002, 0.003, -0.004, 0.01, -0.02,
                        0.015, -0.03, 0.02, -0.04, 0.01, 0.0}) {
        for (int j = 1; j <= 100; ++j) result.push_back(last + (peak-last)*j/100.0);
        last = peak;
    }
    return result;
}
void lifecycle(int n, const char* curveOutput = nullptr) {
    std::unique_ptr<UniaxialMaterial> mat(parse(parameters(n)));
    require(bool(mat), "valid parameters accepted");
    require(near(mat->getInitialTangent(), 206000), "initial tangent");
    mat->setTrialStrain(1.e-6);
    require(near(mat->getStress(), 0.206), "elastic stress E times strain");
    mat->revertToStart();
    std::unique_ptr<UniaxialMaterial> copy;
    CFSTsteel restored;
    MemoryChannel channel;
    FEM_ObjectBroker broker;
    const auto strains = history();
    std::ofstream curve;
    if (curveOutput && n == 1) {
        curve.open(curveOutput);
        require(curve.good(), "curve output opened");
        curve << "step,strain,stress_MPa,tangent_MPa\n" << std::setprecision(17);
    }
    for (unsigned i = 0; i < strains.size(); ++i) {
        require(mat->setTrialStrain(strains[i]) == 0, "valid trial");
        require(std::isfinite(mat->getStress()) && std::isfinite(mat->getTangent()), "finite response");
        if (curve.is_open())
            curve << i+1 << ',' << strains[i] << ',' << mat->getStress() << ',' << mat->getTangent() << '\n';
        if (i % 71 == 0) {
            copy.reset(mat->getCopy());
            equal(*mat, *copy);
            require(mat->sendSelf(7, channel) == 0, "serialize trial");
            require(restored.recvSelf(7, channel, broker) == 0, "restore trial into broker object");
            equal(*mat, restored);

            copy->revertToLastCommit();
            restored.revertToLastCommit();
            equal(*copy, restored);
            copy->setTrialStrain(strains[i]);
            restored.setTrialStrain(strains[i]);
            equal(*mat, *copy);
            equal(*mat, restored);
        } else {
            copy->setTrialStrain(strains[i]);
            restored.setTrialStrain(strains[i]);
            equal(*mat, *copy);
            equal(*mat, restored);
        }
        mat->commitState(); copy->commitState(); restored.commitState();
    }

    copy.reset(mat->getCopy());
    for (double x : {-0.01, 0.02, -0.03}) mat->setTrialStrain(x);
    mat->revertToLastCommit();
    equal(*mat, *copy);
    mat->setTrialStrain(0.005); copy->setTrialStrain(0.005); equal(*mat, *copy);
    mat->revertToStart();
    std::unique_ptr<UniaxialMaterial> fresh(parse(parameters(n)));
    for (double x : strains) {
        mat->setTrialStrain(x); fresh->setTrialStrain(x); equal(*mat, *fresh);
        mat->commitState(); fresh->commitState();
    }
    channel.failSend = true;
    require(mat->sendSelf(8, channel) < 0, "send failure propagated");
    channel.failReceive = true;
    require(restored.recvSelf(8, channel, broker) < 0, "receive failure propagated");
    channel.failReceive = false;
    channel.saved(2) = 9;
    require(restored.recvSelf(7, channel, broker) < 0, "invalid state header rejected");
}
}

extern "C" int OPS_GetNumRemainingInputArgs() { return int(input.size()-position); }
extern "C" int OPS_GetIntInput(int* count, int* values) {
    if (*count < 0 || position+*count > input.size()) return -1;
    for (int i = 0; i < *count; ++i) {
        char* end = nullptr;
        long v = std::strtol(input[position++].c_str(), &end, 10);
        if (!end || *end) return -1;
        values[i] = int(v);
    }
    return 0;
}
extern "C" int OPS_GetDoubleInput(int* count, double* values) {
    if (*count < 0 || position+*count > input.size()) return -1;
    for (int i = 0; i < *count; ++i) {
        char* end = nullptr;
        values[i] = std::strtod(input[position++].c_str(), &end);
        if (!end || *end) return -1;
    }
    return 0;
}
extern "C" int OPS_Error(char* s, int) { std::cerr << s << '\n'; return 0; }

int main(int argc, char** argv) {
    require(argc <= 2, "usage: cfststeel-test [curve-output.csv]");
    for (int n : {0, 1, 2, 8}) lifecycle(n, argc == 2 ? argv[1] : nullptr);
    for (int index : {1, 2, 3, 5, 8, 13, 22}) {
        auto args = parameters(); args[index] = "0";
        require(parse(args) == nullptr, "singular parameter rejected");
    }
    for (const std::string& bad : {"-1", "9", "1.5", "not-a-number"}) {
        auto args = parameters(); args[20] = bad;
        require(parse(args) == nullptr, "invalid N rejected");
    }
    auto args = parameters(); args.push_back("1");
    require(parse(args) == nullptr, "extra arguments rejected");
    args = parameters(); args.pop_back();
    require(parse(args) == nullptr, "missing pair rejected");
    args = parameters(); args[1] = "nan";
    require(parse(args) == nullptr, "nonfinite parameter rejected");
#ifdef CFSTSTEEL_AUTHOR_BASELINE
    OriginalCFSTsteel original(1, 206000, 305, .02, -342.07, 26.58, -74.87, 3.16,
        4.07, 6.38, 3.12, 2.36, 5.82, 1, 1, 1, 1, 1, 116.39, 11.81, {4467.68}, {38.39});
    std::unique_ptr<UniaxialMaterial> current(parse(parameters()));
    for (double x : history()) {
        original.setTrialStrain(x); current->setTrialStrain(x); equal(original, *current);
        original.commitState(); current->commitState();
    }
    std::cout << "PASS: author-source comparison, 1200 cyclic strain steps\n";
#endif
    std::cout << "PASS: CFSTsteel lifecycle, N=0/1/2/8, parser and channel errors; " << checks << " checks\n";
}

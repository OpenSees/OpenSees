// Direct native trials: compare every active history/backbone field, including
// the committed state that must remain unchanged until commitState().
#include <algorithm>
#include <array>
#include <cmath>
#include <cfloat>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <UniaxialMaterial.h>
#include <Vector.h>
// Inspection only; no production layout or constitutive function is changed.
#define private public
#include <HystereticSMMaterial.h>
#undef private

namespace {
using Material = HystereticSMMaterial;
using State = std::vector<std::pair<const char*, double>>;
int checks = 0, failures = 0, stateComparisons = 0, degradedStates = 0;
void check(bool condition, const std::string& context) {
    ++checks;
    if (!condition) {
        ++failures;
        if (failures <= 12) std::cerr << "FAIL: " << context << '\n';
    }
}
State snapshot(const Material& m) {
    State s;
#define ADD(name) s.emplace_back(#name, m.name);
    ADD(Tstrain) ADD(Tstress) ADD(Ttangent) ADD(TloadIndicator)
    ADD(TrotMax) ADD(TrotMin) ADD(TrotPu) ADD(TrotNu) ADD(TenergyD)
    ADD(TanchorRotP) ADD(TanchorStressP) ADD(TanchorRotN) ADD(TanchorStressN)
    ADD(TrotMaxDuctUsed) ADD(TrotMinDuctUsed) ADD(TdamfcPeak)
    ADD(mom1p) ADD(rot1p) ADD(mom1n) ADD(rot1n)
    ADD(mom2p) ADD(mom3p) ADD(mom4p) ADD(mom5p) ADD(mom6p) ADD(mom7p)
    ADD(rot2p) ADD(rot3p) ADD(rot4p) ADD(rot5p) ADD(rot6p) ADD(rot7p)
    ADD(mom2n) ADD(mom3n) ADD(mom4n) ADD(mom5n) ADD(mom6n) ADD(mom7n)
    ADD(rot2n) ADD(rot3n) ADD(rot4n) ADD(rot5n) ADD(rot6n) ADD(rot7n)
    ADD(E1p) ADD(E2p) ADD(E3p) ADD(E4p) ADD(E5p) ADD(E6p) ADD(E7p)
    ADD(E1n) ADD(E2n) ADD(E3n) ADD(E4n) ADD(E5n) ADD(E6n) ADD(E7n)
    ADD(Eup) ADD(Eun) ADD(KeP) ADD(KeN) ADD(energyA)
    ADD(Cstrain) ADD(Cstress) ADD(Ctangent) ADD(CloadIndicator)
    ADD(CrotMax) ADD(CrotMin) ADD(CrotPu) ADD(CrotNu) ADD(CenergyD)
    ADD(CanchorRotP) ADD(CanchorStressP) ADD(CanchorRotN) ADD(CanchorStressN)
    ADD(CrotMaxDuctUsed) ADD(CrotMinDuctUsed) ADD(CdamfcPeak)
    ADD(Cmom2p) ADD(Cmom3p) ADD(Cmom4p) ADD(Cmom5p) ADD(Cmom6p) ADD(Cmom7p)
    ADD(Crot2p) ADD(Crot3p) ADD(Crot4p) ADD(Crot5p) ADD(Crot6p) ADD(Crot7p)
    ADD(Cmom2n) ADD(Cmom3n) ADD(Cmom4n) ADD(Cmom5n) ADD(Cmom6n) ADD(Cmom7n)
    ADD(Crot2n) ADD(Crot3n) ADD(Crot4n) ADD(Crot5n) ADD(Crot6n) ADD(Crot7n)
#undef ADD
    return s;
}
void same(const Material& a, const Material& b, const std::string& context) {
    ++stateComparisons;
    const auto sa = snapshot(a), sb = snapshot(b);
    for (size_t i=0; i<sa.size(); ++i)
        check(std::isfinite(sa[i].second) && sa[i].second == sb[i].second,
              context + ": " + sa[i].first);
}
struct Config { int asymmetric, pinched, damage, degradation, rotY, sevenPoints; int zeroSecondPoint=0; };
std::unique_ptr<Material> makeMaterial(Config c) {
    const int size = c.sevenPoints ? 14 : 6;
    Vector p(size), n(size), pinch(2), damage(2), deg(2), ry(2), empty, internal(48);
    const double moments[] = {100,120,60,80,70,50,45};
    const double rotations[] = {0.001,0.01,0.012,0.014,0.017,0.022,0.03};
    for (int i=0; i<size/2; ++i) {
        p(2*i) = c.zeroSecondPoint && i==1 ? 0 : c.sevenPoints ? moments[i] : (i==0 ? 100 : 120);
        p(2*i+1) = c.sevenPoints ? rotations[i] : (i==0 ? 0.001 : i==1 ? 0.01 : 0.02);
        n(2*i) = -p(2*i) * (c.asymmetric ? (i==0 ? 0.8 : 0.875) : 1);
        n(2*i+1) = -p(2*i+1) * (c.asymmetric && i==0 ? 2 : 1);
    }
    pinch(0)=c.pinched ? 0.4 : 1; pinch(1)=c.pinched ? 0.6 : 1;
    damage(0)=c.damage==2 ? -0.02 : c.damage==1 ? 0.02 : 0;
    damage(1)=c.damage ? 0.01 : 0;
    deg(0)=c.degradation==1 || c.degradation==3 ? -0.2 : c.degradation==2 ? 0.2 : 0;
    deg(1)=c.degradation==3 ? deg(0) : -deg(0);
    ry(0)=c.rotY ? 0.002 : 0; ry(1)=c.rotY ? (c.asymmetric ? -0.003 : -0.002) : 0;
    return std::make_unique<Material>(1,p,n,pinch,damage,c.damage ? 0.3 : 0,
                                      deg,ry,empty,empty,internal);
}
std::vector<Config> configurations() {
    std::vector<Config> result;
    for (int a=0;a<2;++a) for(int p=0;p<2;++p) for(int d=0;d<3;++d)
        for(int e=0;e<3;++e) for(int r=0;r<2;++r) for(int n=0;n<2;++n)
            result.push_back({a,p,d,e,r,n});
    return result;
}
void commitHistory(Material& m, const std::vector<double>& path) {
    for (double x : path) { m.setTrialStrain(x); m.commitState(); }
}
void reproducer(bool assertions) {
    for (int ry=0;ry<2;++ry) {
        Config c={0,1,2,1,ry,0};
        auto direct=makeMaterial(c), dirty=makeMaterial(c);
        commitHistory(*direct,{0.004}); commitHistory(*dirty,{0.004});
        direct->setTrialStrain(0.0041);
        dirty->setTrialStrain(0.0039); dirty->setTrialStrain(0.0041);
        if (assertions) same(*direct,*dirty,"exact reproducer trial");
        std::cout << "reproducer rotY=" << ry
            << " damage direct=" << direct->TdamfcPeak << " dirty=" << dirty->TdamfcPeak
            << " mom2n direct=" << direct->mom2n << " dirty=" << dirty->mom2n << '\n';
        direct->commitState(); dirty->commitState();
        direct->setTrialStrain(-0.005); dirty->setTrialStrain(-0.005);
        if (assertions) same(*direct,*dirty,"exact reproducer subsequent trial");
        std::cout << "later rotY=" << ry << " force direct=" << direct->getStress()
            << " dirty=" << dirty->getStress() << " tangent direct=" << direct->getTangent()
            << " dirty=" << dirty->getTangent() << '\n';
    }
}
void trialOrderTests() {
    const std::vector<std::vector<double>> prefixes={
        {},{0.004},{-0.004},{0.004,0.0039},{-0.004,-0.0039},
        {0.004,-0.004,-0.001,0.002},{-0.004,0.004,0.001,-0.002},
        {0.004,0.002,0,-0.004,-0.002,0,0.006,0.004}};
    int configId=0;
    for (auto c : configurations()) {
        int prefixId=0;
        for (const auto& prefix : prefixes) {
            auto committed=makeMaterial(c); commitHistory(*committed,prefix);
            const bool degraded=committed->mom2p!=committed->mom2p_ref || committed->mom2n!=committed->mom2n_ref;
            if (degraded) ++degradedStates;
            const double a=committed->Cstrain;
            const std::vector<double> targets={0.0041,-0.0041,0.009,-0.009,0.018,-0.018,0.025,-0.025,0.00391,-0.00391,0,
                a,std::nextafter(a,INFINITY),std::nextafter(a,-INFINITY)};
            for (double target : targets) {
                auto direct=makeMaterial(c),dirty=makeMaterial(c),rolledBack=makeMaterial(c);
                for (auto* m : {direct.get(),dirty.get(),rolledBack.get()}) commitHistory(*m,prefix);
                const std::string context="config="+std::to_string(configId)+" prefix="+std::to_string(prefixId)+" C="+std::to_string(target);
                // Reversal, monotonic, full/partial, and repeated trials from the
                // same committed state; nothing in this list is committed.
                for (double rejected : {0.0039,-0.0039,0.009,-0.009,0.001,-0.001,0.0039}) {
                    dirty->setTrialStrain(rejected); rolledBack->setTrialStrain(rejected);
                }
                rolledBack->revertToLastCommit();
                same(*committed,*rolledBack,context+" rollback retains committed degradation");
                direct->setTrialStrain(target); dirty->setTrialStrain(target); rolledBack->setTrialStrain(target);
                same(*direct,*dirty,context+" trial-order independence");
                same(*direct,*rolledBack,context+" rollback recovery");
                if (std::fabs(target-a)<DBL_EPSILON)
                    same(*committed,*direct,context+" exact/negligible return");
                direct->commitState(); dirty->commitState();
                same(*direct,*dirty,context+" committed C");
                for (double following : {-0.005,0.002,-0.003,0.006}) {
                    direct->setTrialStrain(following); dirty->setTrialStrain(following);
                    same(*direct,*dirty,context+" subsequent response");
                    direct->commitState(); dirty->commitState();
                }
            }
            ++prefixId;
        }
        ++configId;
    }
    check(degradedStates>0,"test coverage includes actually committed degraded envelopes");
}
void initializationEdgeTests() {
    for (int sign : {-1,1}) for(int ry=0;ry<2;++ry) {
        // Strength-only degradation can leave both second moments and their
        // locations unchanged while points 3-7 and their slopes change.
        Config c={0,1,2,3,ry,1,1};
        auto direct=makeMaterial(c),dirty=makeMaterial(c);
        commitHistory(*direct,{sign*0.004});commitHistory(*dirty,{sign*0.004});
        direct->setTrialStrain(sign*0.0041);
        dirty->setTrialStrain(sign*0.0039);dirty->setTrialStrain(sign*0.0041);
        same(*direct,*dirty,"zero second-point moment retains all other backbone points");
        direct->commitState();dirty->commitState();
        direct->setTrialStrain(-sign*0.015);dirty->setTrialStrain(-sign*0.015);
        same(*direct,*dirty,"later segments after rejected strength degradation");
        auto virgin=makeMaterial(c);
        dirty->revertToStart();dirty->setTrialStrain(0);
        same(*virgin,*dirty,"zero strain after reset restores the full virgin state");
        dirty->setTrialStrain(std::nextafter(0.0,INFINITY));dirty->setTrialStrain(0);
        same(*virgin,*dirty,"negligible virgin trial followed by zero");
    }
}

// Native prescribed-strain equivalents of tests/test_hystereticsm_reconnection.py.
// These exercise the same constitutive histories without an analysis interpreter.
std::unique_ptr<Material> reloadMaterial() {
    const double values[]={170.0908879721198,0.0008256150522648083,806.5664389032247,0.3299875444832722,
        68.03635518884792,0.33779606261927597,0,0.33794004766395014};
    Vector p(8),n(8),pin(2),empty,internal(48);
    for(int i=0;i<8;++i) {p(i)=values[i];n(i)=-values[i];}
    pin(0)=0.08;pin(1)=0.98;
    return std::make_unique<Material>(1,p,n,pin,empty,0.15,empty,empty,empty,empty,internal);
}
void reloadTests() {
    const double peak=0.0009391377421725002, r1=0.0008256150522648083, epsilon=1e-12;
    for (int sign : {-1,1}) {
        const std::vector<double> prefix={sign*peak,sign*0.0008689933708300856,sign*-0.0004230444166321932};
        auto m=reloadMaterial();commitHistory(*m,prefix);
        double inside=sign*(peak-epsilon);
        m->setTrialStrain(inside); double stress=m->getStress(),tangent=m->getTangent();m->commitState();
        m->setTrialStrain(sign*peak);
        double residual=m->getStress()-stress-tangent*(sign*peak-inside);
        check(std::fabs(residual)<1e-9,"retained envelope reconnects");
        std::cout << "reload sign=" << sign << " reconnectionResidual=" << residual;
        for (double step : {r1*1e-11,r1*1e-7}) {
            m=reloadMaterial();commitHistory(*m,{sign*2*r1,sign*1.9*r1});
            stress=m->getStress();tangent=m->getTangent();
            m->setTrialStrain(sign*(1.9*r1+step));
            residual=m->getStress()-stress-tangent*sign*step;
            if (step==r1*1e-11) {
                check(std::fabs(residual)<1e-9,"partial reload is continuous");
                std::cout << " partialResidual=" << residual;
            } else check(std::fabs((m->getStress()-stress)/(sign*step)-m->getTangent())<1e-3,
                         "partial reload tangent agrees with finite difference");
        }
        m=reloadMaterial();commitHistory(*m,prefix);
        const double probe=prefix.back()+sign*epsilon;
        m->setTrialStrain(probe);
        const double zero=probe-m->getStress()/m->getTangent(),target=sign*0.00055;
        const double offset=std::fabs(target-prefix.back())*1e-6;
        std::vector<double> forces,tangents;
        for (auto extra : std::vector<std::vector<double>>{{},{zero},{zero-sign*offset},{zero+sign*offset}}) {
            m=reloadMaterial();commitHistory(*m,prefix);commitHistory(*m,extra);m->setTrialStrain(target);
            forces.push_back(m->getStress());tangents.push_back(m->getTangent());
        }
        const double df=*std::max_element(forces.begin(),forces.end())-*std::min_element(forces.begin(),forces.end());
        const double dk=*std::max_element(tangents.begin(),tangents.end())-*std::min_element(tangents.begin(),tangents.end());
        check(df<1e-9 && dk<1e-7,"full reload zero crossing is subdivision independent");
        std::cout << " zeroCrossingForceSpread=" << df << " tangentSpread=" << dk << '\n';
    }
}
void acceptedTrace() {
    const std::vector<double> path={0.0005,0.001,0.002,0.004,0.0039,0.004,0.002,0,-0.001,-0.004,
        -0.0039,-0.004,-0.002,0,0.004,0.006,0.003,0.0031,0.006,-0.006,-0.003,-0.0031,-0.006,0.003,0.006,-0.006,0.006,0};
    std::vector<double> extended=path;
    extended.insert(extended.end(),{0.011,0.015,0.021,0.024,0.015,-0.011,-0.015,-0.021,-0.024,-0.015,0});
    const auto fields=snapshot(*makeMaterial({0,0,0,0,0,0}));
    std::cout << "config,step";
    for (auto f : fields) std::cout << ',' << f.first;
    std::cout << '\n';
    int id=0;
    for (auto c : configurations()) {
        auto m=makeMaterial(c);int i=0;
        for(double x : extended) {
            m->setTrialStrain(x);
            std::cout << id << ',' << i++;
            for(auto f : snapshot(*m)) std::cout << ',' << f.second;
            std::cout << '\n';m->commitState();
        }
        ++id;
    }
}
} // namespace

#ifndef HYSTERETICSM_TEST_NO_MAIN
int main(int argc,char** argv) {
    std::cout << std::setprecision(17);
    if(argc>1 && std::string(argv[1])=="--trace") {acceptedTrace();return 0;}
    if(argc>1 && std::string(argv[1])=="--reproducer") {reproducer(false);return 0;}
    reproducer(true);trialOrderTests();initializationEdgeTests();reloadTests();
    std::cout << "configurations=" << configurations().size() << " stateFields=" << snapshot(*makeMaterial({0,0,0,0,0,0})).size()
              << " stateComparisons=" << stateComparisons << " committedDegradedPrefixes=" << degradedStates
              << " checks=" << checks << " failures=" << failures << '\n';
    return failures ? 1 : 0;
}
#endif

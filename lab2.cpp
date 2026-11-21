/*
 * ============================================================================
 * README - Simulatore Integrato di Perdita di Energia (Bethe-Bloch CSDA),
 *          Spettro Cosmico e Accettanza Geometrica 3D (Telescopio 3 Strati)
 * ============================================================================
 */

#include <TGClient.h>
#include <TCanvas.h>
#include <TGComboBox.h>
#include <TGLabel.h>
#include <TGFrame.h>
#include <TGButton.h>
#include <TGTab.h>
#include <TRootEmbeddedCanvas.h>
#include <TGTextEntry.h>
#include <TGNumberEntry.h>
#include <TGLayout.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TPad.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <TAxis.h>
#include <TLine.h>
#include <TBox.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TPolyLine3D.h>
#include <TView.h>
#include <TMath.h>
#include <RQ_OBJECT.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TSystem.h>
#include <TColor.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

#define MAX_LAYERS 10
#define MAX_COMPONENTS 5

namespace {
const double kPi = 3.14159265358979323846;
const double kElectronMassMeV = 0.51099895000;
const double kMuonMassMeV = 105.6583755;
const double kBetheK = 0.307075;
const double kLn10Times2 = 4.605170185988092;
const double kEtaCut = 0.1;
const double kAvogadro = 6.02214076e23;
const int kRangeSteps = 12000;
const int kEnergyPlotPoints = 600;
const int kNGrid = 250;            // Griglia per l'integrale numerico esatto Tetto-Fondo
const double kScanRefSize = 229.0; // Riferimento scala scansione piano: 31 + 6*33 cm
}

// =========================================================================
// STRUTTURE DATI: BETHE-BLOCH, COMPOSTI & SPETTRO RASTIN
// =========================================================================
struct Medium {
    string name;
    double zOverA = 0.0, ionEV = 0.0, density = 0.0, zEffective = 0.0;
    bool hasCustomDensityParams = false;
    double cBar = 0.0, x0 = 0.0, x1 = 0.0, a = 0.0, m = 3.0, delta0 = 0.0;
};

struct IncomingPart {
    string name;
    int charge;
    double massMeV;
};

struct Results {
    double massStopping, stopping, rangeMass, rangeCm, integratedMass, tailMass;
};

struct LayerData {
    Medium mat;
    double thicknessCm, arealThickness, kineticIn, kineticOut, momentumIn, momentumOut;
};

struct CompoundInfo {
    bool valid = false;
    string name, formula;
    double density = 0.0, molarMass = 0.0, atomsPerMol = 0.0, atomDensity = 0.0;
    double zOverA = 0.0, zMean = 0.0, aMean = 0.0, ionBragg = 0.0, ionUsed = 0.0;
    bool ionOverridden = false;
    vector<int> el;
    vector<double> n, w, x;
    Medium medium;
};

// =========================================================================
// STRUTTURE DATI: MONTE CARLO GEOMETRICO (5 CATEGORIE & SCAN PIANO S)
// =========================================================================
struct TrackSample {
    double x0, y0, z0;
    double xIn, yIn, zIn;
    double xOut, yOut, zOut;
    int category; // 0=(a) Coincidenza, 4=(e) Miss
};

struct PlaneScanPoint {
    double lGen = 0.0;        // Lato piano generato [cm]
    double areaGen = 0.0;     // S = L_gen^2 [cm^2]
    long long nGen = 0;
    long long nHits = 0;
    double accPct = 0.0;      // Accettanza adimensionale A [%]
    double accPctErr = 0.0;   // Incertezza binomiale [%]
    double geomFactor = 0.0;  // G_MC [cm^2 sr]
    double geomErr = 0.0;     // Incertezza su G_MC [cm^2 sr]
    double relDiffPct = 0.0;  // Scarto % rispetto a G_Ana
};

struct MuonMCResult {
    double lxIdeal = 0.0;
    double lyIdeal = 0.0;
    double areaIdeal = 0.0;      // S = Lx * Ly [cm^2]
    double geomIdeal = 0.0;      // G_piano = S * 2*pi/(n+2) [cm^2 sr]
    double rateIdealHz = 0.0;    // nu_piano = eff * I0 * G_piano [Hz]

    double areaTop = 0.0;        // X * Y [cm^2]
    double areaSideX = 0.0;      // X * Z [cm^2]
    double areaSideY = 0.0;      // Y * Z [cm^2]
    double areaSidesTot = 0.0;   // 2*(X*Z + Y*Z) [cm^2]

    long long nEvents = 0;
    long long counts[5] = {0, 0, 0, 0, 0};
    long long countRealTot = 0;
    long long countTopEntry = 0;
    long long countSideEntry = 0;

    double acc[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    double accErr[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    double accRealTot = 0.0;
    double accRealTotErr = 0.0;
    double accTopEntry = 0.0;
    double accTopErr = 0.0;
    double accSideEntry = 0.0;
    double accSideErr = 0.0;

    double geomMC = 0.0;
    double geomMCErr = 0.0;
    double geomAna = 0.0;

    double rateMC[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    double rateMCErr[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    double rateMCRealTotHz = 0.0;
    double rateMCRealTotErrHz = 0.0;
    double rateMCTopHz = 0.0;
    double rateMCTopErrHz = 0.0;
    double rateMCSideHz = 0.0;
    double rateMCSideErrHz = 0.0;

    double meanTrackLen[4] = {0.0, 0.0, 0.0, 0.0};
    double meanTrackLenTot = 0.0;
    double meanTrackLenExp = 0.0;

    double rateAnaTopHz = 0.0;
    double rateAnaSideXHz = 0.0;
    double rateAnaSideYHz = 0.0;
    double rateAnaSidesTotHz = 0.0;
    double rateAnaTotal3DHz = 0.0;
    double rateAnaTopBotHz = 0.0;
    double rateAnaTopSideHz = 0.0;
    double rateAnaSideBotHz = 0.0;
    double rateAnaSideSideHz = 0.0;

    // Geometria Telescopio usata nel Run
    double zGen = 0.0;
    double x1=0, y1=0, z1_min=0, z1_max=0;
    double x2=0, y2=0, z2_min=0, z2_max=0;
    double x3=0, y3=0, z3_min=0, z3_max=0;

    std::vector<TrackSample> tracks3D;
    std::vector<PlaneScanPoint> planeScan;
};

vector<Medium> materials;
vector<IncomingPart> particles;

// =========================================================================
// FUNZIONI DI UTILITA' E FISICA BETHE-BLOCH / RASTIN
// =========================================================================
static string Trim(const string& s) {
    const string whitespace = " \t\r\n";
    const size_t first = s.find_first_not_of(whitespace);
    if (first == string::npos) return "";
    const size_t last = s.find_last_not_of(whitespace);
    return s.substr(first, last - first + 1);
}

static vector<string> SplitTabs(const string& line) {
    vector<string> fields;
    string field;
    istringstream iss(line);
    while (getline(iss, field, '\t')) fields.push_back(Trim(field));
    return fields;
}

static bool ParseDouble(const string& text, double& value) {
    const string s = Trim(text);
    if (s.empty()) return false;
    char* end = nullptr;
    value = strtod(s.c_str(), &end);
    return end != s.c_str() && end != nullptr && *end == '\0' && std::isfinite(value);
}

static bool ParseInt(const string& text, int& value) {
    double x = 0.0;
    if (!ParseDouble(text, x)) return false;
    const double rounded = std::round(x);
    if (std::fabs(x - rounded) > 1e-9 || rounded < -1000000 || rounded > 1000000) return false;
    value = static_cast<int>(rounded);
    return true;
}

static string G(double v, int prec = 6) {
    ostringstream os;
    os << setprecision(prec) << v;
    return os.str();
}

static void PutText(double x, double y, const char *s, double px, int align = 12) {
    TLatex *t = new TLatex(x, y, s);
    t->SetNDC(true);
    t->SetTextFont(43);
    t->SetTextSize(px);
    t->SetTextAlign(align);
    t->SetBit(kCanDelete);
    t->Draw();
}

static double PlasmaEnergyEV(const Medium& mat) {
    return 28.816 * std::sqrt(mat.density * mat.zOverA);
}

static void BuildApproxDensityParams(Medium& mat) {
    const double plasmaEV = PlasmaEnergyEV(mat);
    if (!(plasmaEV > 0.0) || !(mat.ionEV > 0.0)) {
        mat.cBar = mat.a = mat.delta0 = 0.0;
        mat.x0 = 0.2;
        mat.x1 = 3.0;
        mat.m = 3.0;
        mat.hasCustomDensityParams = false;
        return;
    }
    mat.cBar = 1.0 + 2.0 * std::log(mat.ionEV / plasmaEV);
    mat.m = 3.0;
    mat.delta0 = 0.0;
    if (mat.density > 0.05) {
        if (mat.ionEV < 100.0) {
            mat.x1 = 2.0;
            mat.x0 = (mat.cBar <= 3.681) ? 0.2 : (0.326 * mat.cBar - 1.0);
        } else {
            mat.x1 = 3.0;
            mat.x0 = (mat.cBar <= 5.215) ? 0.2 : (0.326 * mat.cBar - 1.5);
        }
    } else {
        if (mat.cBar < 10.0) { mat.x0 = 1.6; mat.x1 = 4.0; }
        else if (mat.cBar < 10.5) { mat.x0 = 1.7; mat.x1 = 4.0; }
        else if (mat.cBar < 11.0) { mat.x0 = 1.8; mat.x1 = 4.0; }
        else if (mat.cBar < 11.5) { mat.x0 = 1.9; mat.x1 = 4.0; }
        else if (mat.cBar < 12.25) { mat.x0 = 2.0; mat.x1 = 4.0; }
        else if (mat.cBar < 13.804) { mat.x0 = 2.0; mat.x1 = 5.0; }
        else { mat.x0 = 0.326 * mat.cBar - 2.5; mat.x1 = 5.0; }
    }
    mat.a = (mat.x1 > mat.x0) ? (mat.cBar - kLn10Times2 * mat.x0) / std::pow(mat.x1 - mat.x0, mat.m) : 0.0;
    mat.hasCustomDensityParams = false;
}

static double DensityCorrection(double eta, const Medium& mat) {
    if (!(eta > 0.0) || !(mat.ionEV > 0.0) || !(mat.density > 0.0) || !(mat.zOverA > 0.0)) return 0.0;
    const double x = std::log10(eta);
    double delta = 0.0;
    if (x >= mat.x1) delta = kLn10Times2 * x - mat.cBar;
    else if (x >= mat.x0 && mat.x1 > mat.x0) delta = kLn10Times2 * x - mat.cBar + mat.a * std::pow(mat.x1 - x, mat.m);
    else if (mat.delta0 > 0.0) delta = mat.delta0 * std::pow(10.0, 2.0 * (x - mat.x0));
    return std::max(0.0, delta);
}

static double ShellCorrection(double eta, double ionEV) {
    if (!(eta >= kEtaCut * (1.0 - 1.0e-12)) || !(ionEV > 0.0)) return 0.0;
    const double inv2 = 1.0 / (eta * eta), inv4 = inv2 * inv2, inv6 = inv4 * inv2;
    const double i2 = ionEV * ionEV, i3 = i2 * ionEV;
    return ((0.422377 * inv2 + 0.0304043 * inv4 - 0.00038106 * inv6) * 1.0e-6 * i2) +
           ((3.850190 * inv2 - 0.1667989 * inv4 + 0.00157955 * inv6) * 1.0e-9 * i3);
}

static double ComputeBetheBloch(double kineticMeV, double projectileMassMeV, int charge, const Medium& mat) {
    if (!(kineticMeV > 0.0) || !(projectileMassMeV > 10.0 * kElectronMassMeV) || charge == 0 ||
        !(mat.zOverA > 0.0) || !(mat.ionEV > 0.0) || !(mat.density > 0.0) || !(mat.zEffective > 0.0)) return 0.0;

    const double eta = std::sqrt(kineticMeV * (kineticMeV + 2.0 * projectileMassMeV)) / projectileMassMeV;
    const double gamma = std::sqrt(1.0 + eta * eta), beta2 = (eta * eta) / (1.0 + eta * eta);
    if (!(gamma > 1.0) || !(beta2 > 0.0 && beta2 < 1.0) || eta < kEtaCut * (1.0 - 1.0e-12)) return 0.0;

    const double massRatio = kElectronMassMeV / projectileMassMeV;
    const double tMax = (2.0 * kElectronMassMeV * eta * eta) / (1.0 + 2.0 * gamma * massRatio + massRatio * massRatio);
    const double logArg = (2.0 * kElectronMassMeV * eta * eta * tMax) / (mat.ionEV * 1e-6 * mat.ionEV * 1e-6);
    if (!(logArg > 0.0)) return 0.0;

    const double bracket = 0.5 * std::log(logArg) - beta2 - ShellCorrection(eta, mat.ionEV) / mat.zEffective - 0.5 * DensityCorrection(eta, mat);
    const double stopping = kBetheK * charge * charge * mat.zOverA / beta2 * bracket;
    return (std::isfinite(stopping) && stopping > 0.0) ? stopping : 0.0;
}

static double KineticAtEtaCut(double massMeV) {
    return massMeV * (std::sqrt(1.0 + kEtaCut * kEtaCut) - 1.0);
}

static double MomentumFromKinetic(double kineticMeV, double massMeV) {
    return std::sqrt(kineticMeV * (kineticMeV + 2.0 * massMeV));
}

static bool ComputeRangeMassCSDA(double kineticMeV, double projectileMassMeV, int charge, const Medium& mat, double& rangeMass) {
    rangeMass = 0.0;
    const double tCut = KineticAtEtaCut(projectileMassMeV);
    if (!std::isfinite(kineticMeV) || kineticMeV < tCut * (1.0 - 1.0e-12) || !(mat.density > 0.0)) return false;

    double previousT = kineticMeV, previousS = ComputeBetheBloch(previousT, projectileMassMeV, charge, mat);
    if (!(previousS > 0.0) || !std::isfinite(previousS)) return false;

    double integrated = 0.0;
    for (int i = 1; i < kRangeSteps; ++i) {
        const double currentT = kineticMeV * std::pow(tCut / kineticMeV, static_cast<double>(i) / (kRangeSteps - 1));
        const double currentS = ComputeBetheBloch(currentT, projectileMassMeV, charge, mat);
        if (!(currentS > 0.0) || !std::isfinite(currentS)) return false;
        integrated += 0.5 * (previousT - currentT) * (1.0 / previousS + 1.0 / currentS);
        previousT = currentT;
        previousS = currentS;
    }

    const double sCut = ComputeBetheBloch(tCut, projectileMassMeV, charge, mat);
    if (!(sCut > 0.0) || !std::isfinite(sCut)) return false;
    rangeMass = integrated + (tCut / sCut);
    return std::isfinite(rangeMass) && rangeMass > 0.0;
}

static bool FindKineticForRange(double targetRangeMass, double projectileMassMeV, int charge, const Medium& mat, double& kineticMeV, string& error) {
    kineticMeV = 0.0;
    const double tCut = KineticAtEtaCut(projectileMassMeV);
    double rLow = 0.0;
    if (!ComputeRangeMassCSDA(tCut, projectileMassMeV, charge, mat, rLow)) {
        error = "S(T) non valido.";
        return false;
    }
    if (targetRangeMass < rLow) {
        error = "Spessore troppo sottile.";
        return false;
    }

    double low = tCut, high = std::max(100.0, 2.0 * tCut), rHigh = 0.0;
    bool bracketed = false;
    for (int i = 0; i < 32; ++i) {
        if (!ComputeRangeMassCSDA(high, projectileMassMeV, charge, mat, rHigh)) return false;
        if (rHigh >= targetRangeMass) { bracketed = true; break; }
        high *= 2.0;
        if (high > 1.0e7) break;
    }
    if (!bracketed) {
        error = "Limite energetico superato.";
        return false;
    }

    for (int i = 0; i < 64; ++i) {
        const double mid = 0.5 * (low + high);
        double rMid = 0.0;
        if (!ComputeRangeMassCSDA(mid, projectileMassMeV, charge, mat, rMid)) return false;
        if (rMid < targetRangeMass) low = mid;
        else high = mid;
        if ((high - low) <= 1.0e-7 * std::max(1.0, mid)) break;
    }
    kineticMeV = 0.5 * (low + high);
    return true;
}

// =========================================================================
// TABELLA RAPIDA RANGE-ENERGIA CSDA (per interpolazione veloce su N tracce)
// =========================================================================
struct FastRangeTable {
    vector<double> T_grid; // MeV
    vector<double> R_grid; // g/cm^2
    double tCut = 0.0;
    double rCut = 0.0;

    void Build(double massMeV, int charge, const Medium& mat, double tMaxMeV = 2.0e6, int nPts = 500) {
        T_grid.clear();
        R_grid.clear();
        tCut = KineticAtEtaCut(massMeV);
        ComputeRangeMassCSDA(tCut, massMeV, charge, mat, rCut);
        T_grid.reserve(nPts);
        R_grid.reserve(nPts);
        for (int i = 0; i < nPts; ++i) {
            double T = tCut * std::pow(tMaxMeV / tCut, static_cast<double>(i) / (nPts - 1));
            double R = 0.0;
            if (ComputeRangeMassCSDA(T, massMeV, charge, mat, R)) {
                T_grid.push_back(T);
                R_grid.push_back(R);
            }
        }
    }

    double RangeFromT(double T) const {
        if (T <= 0.0) return 0.0;
        if (T <= tCut) return rCut * (T / tCut);
        if (T >= T_grid.back()) return R_grid.back() + (T - T_grid.back()) / 2.2;
        auto it = std::lower_bound(T_grid.begin(), T_grid.end(), T);
        size_t idx = std::distance(T_grid.begin(), it);
        if (idx == 0) return R_grid[0];
        double f = (T - T_grid[idx - 1]) / (T_grid[idx] - T_grid[idx - 1]);
        return R_grid[idx - 1] + f * (R_grid[idx] - R_grid[idx - 1]);
    }

    double TFromRange(double R) const {
        if (R <= 0.0) return 0.0;
        if (R <= rCut) return tCut * (R / rCut);
        if (R >= R_grid.back()) return T_grid.back() + (R - R_grid.back()) * 2.2;
        auto it = std::lower_bound(R_grid.begin(), R_grid.end(), R);
        size_t idx = std::distance(R_grid.begin(), it);
        if (idx == 0) return T_grid[0];
        double f = (R - R_grid[idx - 1]) / (R_grid[idx] - R_grid[idx - 1]);
        return T_grid[idx - 1] + f * (T_grid[idx] - T_grid[idx - 1]);
    }

    double EnergyLoss(double T_in, double dX) const {
        if (T_in <= 0.0 || dX <= 0.0) return 0.0;
        double R_in = RangeFromT(T_in);
        if (R_in <= dX) return T_in;
        double T_out = TFromRange(R_in - dX);
        return std::max(0.0, T_in - T_out);
    }
};

// =========================================================================
// SPETTRO MUONI COSMICI (Rastin 1984a)
// =========================================================================
static double RastinFluxP(double pGeV) {
    static const int n_tab = 18;
    static const double p_tab[n_tab] = {
        0.35, 0.40, 0.50, 0.60, 0.80, 1.0, 2.0, 5.0, 10.0,
        20.0, 50.0, 100.0, 200.0, 500.0, 1000.0, 1500.0, 2000.0, 3000.0
    };
    static const double flux_p_tab[n_tab] = {
        2.85e-3, 2.90e-3, 2.94e-3, 2.92e-3, 2.80e-3, 2.62e-3, 1.69e-3, 5.17e-4, 1.35e-4,
        2.58e-5, 2.11e-6, 2.77e-7, 3.12e-8, 1.45e-9, 1.31e-10, 3.14e-11, 1.13e-11, 2.67e-12
    };

    if (pGeV <= 0.0) return 0.0;
    if (pGeV <= p_tab[0]) return flux_p_tab[0];

    for (int i = 0; i < n_tab - 1; ++i) {
        if (pGeV >= p_tab[i] && pGeV <= p_tab[i + 1]) {
            double t = (std::log(pGeV) - std::log(p_tab[i])) / (std::log(p_tab[i + 1]) - std::log(p_tab[i]));
            return std::exp(std::log(flux_p_tab[i]) + t * (std::log(flux_p_tab[i + 1]) - std::log(flux_p_tab[i])));
        }
    }
    double slope = (std::log(flux_p_tab[n_tab - 1]) - std::log(flux_p_tab[n_tab - 2])) /
                   (std::log(p_tab[n_tab - 1]) - std::log(p_tab[n_tab - 2]));
    return flux_p_tab[n_tab - 1] * std::pow(pGeV / p_tab[n_tab - 1], slope);
}

static double IntegrateRastinMomentum(double p1GeV, double p2GeV) {
    if (!(p2GeV > p1GeV) || p2GeV <= 0.0) return 0.0;
    if (p1GeV < 0.0) p1GeV = 0.0;
    const int N = 2000;
    const double h = (p2GeV - p1GeV) / N;
    double sum = 0.5 * (RastinFluxP(p1GeV) + RastinFluxP(p2GeV));
    for (int i = 1; i < N; ++i) sum += RastinFluxP(p1GeV + i * h);
    return sum * h;
}

static double IntegralRastinAboveP(double pMinGeV) {
    const double p_ref = 0.35;
    const double int_ref = 9.13e-3;
    if (pMinGeV <= p_ref) return int_ref + IntegrateRastinMomentum(pMinGeV, p_ref);
    return std::max(0.0, int_ref - IntegrateRastinMomentum(p_ref, pMinGeV));
}

// =========================================================================
// COSTRUTTORE COMPOSTI & LETTURA FILE
// =========================================================================
static bool BuildCompound(const string& nameIn, const vector<int>& idx, const vector<double>& cnt, double rho, double ionOverride, CompoundInfo& info, string& error) {
    info = CompoundInfo();
    if (idx.empty() || idx.size() != cnt.size()) { error = "Nessun componente."; return false; }
    if (!(rho > 0.0)) { error = "Densita' non valida."; return false; }

    string name = Trim(nameIn);
    if (name.empty()) name = "Composto";

    double molar = 0.0, sumN = 0.0, sumNZ = 0.0;
    for (size_t k = 0; k < idx.size(); ++k) {
        if (idx[k] < 0 || static_cast<size_t>(idx[k]) >= materials.size() || !(cnt[k] > 0.0)) {
            error = "Componente non valido.";
            return false;
        }
        const Medium& m = materials[idx[k]];
        double A = m.zEffective / m.zOverA;
        molar += cnt[k] * A;
        sumN += cnt[k];
        sumNZ += cnt[k] * m.zEffective;
    }

    info.name = name;
    info.density = rho;
    info.molarMass = molar;
    info.atomsPerMol = sumN;
    info.atomDensity = rho * kAvogadro * sumN / molar;
    info.zOverA = sumNZ / molar;
    info.zMean = sumNZ / sumN;
    info.aMean = molar / sumN;
    info.el = idx;
    info.n = cnt;

    double lnINum = 0.0;
    ostringstream formula;
    for (size_t k = 0; k < idx.size(); ++k) {
        const Medium& m = materials[idx[k]];
        double A = m.zEffective / m.zOverA;
        double w = cnt[k] * A / molar;
        info.w.push_back(w);
        info.x.push_back(cnt[k] / sumN);
        lnINum += w * m.zOverA * std::log(m.ionEV);
        formula << m.name;
        if (std::fabs(cnt[k] - 1.0) > 1e-12) formula << G(cnt[k], 6);
    }

    info.formula = formula.str();
    info.ionBragg = std::exp(lnINum / info.zOverA);
    info.ionOverridden = (ionOverride > 0.0);
    info.ionUsed = info.ionOverridden ? ionOverride : info.ionBragg;

    Medium eff;
    eff.name = name;
    eff.zOverA = info.zOverA;
    eff.ionEV = info.ionUsed;
    eff.density = rho;
    eff.zEffective = info.zMean;
    BuildApproxDensityParams(eff);

    info.medium = eff;
    info.valid = true;
    return true;
}

static string DescribeCompound(const CompoundInfo& c) {
    ostringstream os;
    os << "\n================================================================================\n";
    os << " COSTRUTTORE COMPOSTO: " << c.name << "   (" << c.formula << ")\n";
    os << "================================================================================\n";
    os << "  Massa molare M_c     : " << G(c.molarMass, 7) << " g/mol\n";
    os << "  Unita' per molecola  : " << G(c.atomsPerMol, 6) << "\n";
    os << "  Densita' rho_c       : " << G(c.density, 6) << " g/cm^3\n";
    os << "  Densita' molecolare  : " << G(c.atomDensity, 5) << " molecole/cm^3\n\n";
    os << "  Comp       |    n_i   |   A_eff   | W_i (massa) |  Z_eff  | I_i [eV] | rho_i [g/cm^3]\n";
    os << "  ------------------------------------------------------------------------------------\n";
    for (size_t k = 0; k < c.el.size(); ++k) {
        const Medium& m = materials[c.el[k]];
        os << "  " << left << setw(10) << m.name << right
           << " | " << setw(8) << G(c.n[k], 5)
           << " | " << setw(9) << G(m.zEffective / m.zOverA, 6)
           << " | " << setw(11) << G(c.w[k], 5)
           << " | " << setw(7) << G(m.zEffective, 4)
           << " | " << setw(8) << G(m.ionEV, 5)
           << " | " << setw(14) << G(m.density, 4) << "\n";
    }
    os << "\n  Parametri efficaci (regola di Bragg nel logaritmo di Bethe):\n";
    os << "    Z/A = sum w_i (Z_i/A_i) = " << G(c.zOverA, 6) << "\n";
    os << "    I_Bragg                 = " << G(c.ionBragg, 5) << " eV\n";
    if (c.ionOverridden) os << "    I usato (override)      = " << G(c.ionUsed, 5) << " eV\n";
    os << "================================================================================\n";
    return os.str();
}

void readfile() {
    materials.clear();
    particles.clear();
    ifstream in("materials.txt");
    if (in.is_open()) {
        string line;
        while (getline(in, line)) {
            line = Trim(line);
            if (line.empty() || line[0] == '#') continue;
            vector<string> f = SplitTabs(line);
            if (f.size() < 5) continue;
            Medium m;
            m.name = f[0];
            if (!ParseDouble(f[1], m.zOverA) || !ParseDouble(f[2], m.ionEV) || !ParseDouble(f[3], m.density) || !ParseDouble(f[4], m.zEffective) ||
                !(m.zOverA > 0.0) || !(m.ionEV > 0.0) || !(m.density > 0.0) || !(m.zEffective > 0.0) || m.name.empty()) continue;
            if (f.size() >= 11) {
                double vals[6];
                bool ok = true;
                for (int j = 0; j < 6; ++j) ok = ok && ParseDouble(f[5 + j], vals[j]);
                if (ok && vals[2] > vals[1] && vals[4] > 0.0 && vals[5] >= 0.0) {
                    m.cBar = vals[0];
                    m.x0 = vals[1];
                    m.x1 = vals[2];
                    m.a = vals[3];
                    m.m = vals[4];
                    m.delta0 = vals[5];
                    m.hasCustomDensityParams = true;
                }
            }
            if (!m.hasCustomDensityParams) BuildApproxDensityParams(m);
            materials.push_back(m);
        }
    } else {
        cerr << "[WARNING] File 'materials.txt' non trovato nella directory corrente: disponibile solo materiale definito dall'utente." << endl;
    }

    ifstream in2("particles.txt");
    if (in2.is_open()) {
        string line;
        while (getline(in2, line)) {
            line = Trim(line);
            if (line.empty() || line[0] == '#') continue;
            vector<string> f = SplitTabs(line);
            if (f.size() < 3) continue;
            IncomingPart p;
            if (f[0].empty() || !ParseInt(f[1], p.charge) || !ParseDouble(f[2], p.massMeV) || p.charge == 0 || !(p.massMeV > 10.0 * kElectronMassMeV)) continue;
            p.name = f[0];
            particles.push_back(p);
        }
    } else {
        cerr << "[WARNING] File 'particles.txt' non trovato nella directory corrente: disponibile solo particella definita dall'utente." << endl;
    }
}

// =========================================================================
// CLASSE GUI UNIFICATA
// =========================================================================
class UnifiedDetectorGUI {
    RQ_OBJECT("UnifiedDetectorGUI")

private:
    TGMainFrame *fMain = nullptr;

    TGTab *fLeftTabs = nullptr;
    TGTab *fRightTabs = nullptr;

    TRootEmbeddedCanvas *fCanvasBethe = nullptr;
    TRootEmbeddedCanvas *fCanvasMC1 = nullptr;
    TRootEmbeddedCanvas *fCanvasMC2 = nullptr;

    // ---------------------------------------------------------------------
    // CONTROLLI SEZIONE BETHE-BLOCH & SPETTRO PMT
    // ---------------------------------------------------------------------
    TGVerticalFrame *leftBetheFrame = nullptr;
    TGLabel *mater = nullptr, *modeLabel = nullptr, *inputValueLabel = nullptr;
    TGNumberEntry *inputValue = nullptr, *userCharge = nullptr, *userMass = nullptr;
    TGNumberEntry *userZ = nullptr, *userA = nullptr, *userIP = nullptr, *userDens = nullptr;
    TGComboBox *inputModeList = nullptr, *matList = nullptr, *partList = nullptr;

    TGGroupFrame *userdefinedpart = nullptr;
    TGGroupFrame *userdefinedmat = nullptr;
    TGGroupFrame *multiLayerFrame = nullptr;
    bool fMultiLayerMode = false;
    TGTextButton *btnToggleMultiMode = nullptr;

    TGNumberEntry *numLayersEntry = nullptr;
    TGHorizontalFrame *layerFrames[MAX_LAYERS];
    TGComboBox *layerMatList[MAX_LAYERS];
    TGNumberEntry *layerThickness[MAX_LAYERS];

    TGLabel *answerlabel[4] = {nullptr}, *calc[4] = {nullptr};
    TGVerticalFrame *calcframe = nullptr;

    // Pannello a scomparsa: Costruttore Composti
    TGVerticalFrame *centerframe = nullptr;
    bool fBuilderVisible = false;
    TGTextEntry *fCompName = nullptr;
    TGNumberEntry *fCompDens = nullptr, *fCompI = nullptr;
    TGComboBox *fCompElem[MAX_COMPONENTS];
    TGNumberEntry *fCompCount[MAX_COMPONENTS];

    // Pannello a scomparsa: Rate Cosmici & Efficienza PMT
    TGVerticalFrame *cosmicframe = nullptr;
    bool fCosmicVisible = false;
    TGNumberEntry *fDetLx = nullptr, *fDetLy = nullptr, *fDetLz = nullptr;
    TGNumberEntry *fThresholdMeV = nullptr;
    TGLabel *fCosmicAccLabel = nullptr;
    TGLabel *fCosmicRateStopLabel = nullptr;
    TGLabel *fCosmicRatePassLabel = nullptr;
    TGLabel *fCosmicRateTotLabel = nullptr;
    TGLabel *fCosmicEffLabel = nullptr;
    double fLastPMTEfficiencyPct = 100.0;

    // ---------------------------------------------------------------------
    // CONTROLLI SEZIONE MONTE CARLO GEOMETRICO (5 CATEGORIE & SCAN S)
    // ---------------------------------------------------------------------
    TGNumberEntry *fX = nullptr;
    TGNumberEntry *fY = nullptr;
    TGNumberEntry *fZ = nullptr;
    TGNumberEntry *fLx = nullptr;
    TGNumberEntry *fLy = nullptr;
    TGNumberEntry *fNEvents = nullptr;
    TGNumberEntry *fSeed = nullptr;
    TGNumberEntry *fI0 = nullptr;
    TGNumberEntry *fExponent = nullptr;
    TGNumberEntry *fEfficiency1 = nullptr;
    TGNumberEntry *fEfficiency2 = nullptr;
    TGNumberEntry *fEfficiency3 = nullptr;
    TGLabel *fRes[16] = {nullptr};

    // ---------------------------------------------------------------------
    // CACHE DI TUTTI I GRAFICI E ISTOGRAMMI PER ESPORTAZIONE IN FILE ROOT
    // ---------------------------------------------------------------------
    TGraph *fSavedGrBethe = nullptr;
    TGraph *fSavedGrBragg = nullptr;
    TGraph *fSavedGrRastin = nullptr;
    TH1D   *fSavedHistEdepAll = nullptr;
    TH1D   *fSavedHistEdepLost = nullptr;
    TH1D   *fSavedHistEdepPass = nullptr;
    TH1D   *fSavedHistMCRates = nullptr;
    TH1D   *fSavedHistAnaRates = nullptr;
    TGraphErrors *fSavedGrScanAcc = nullptr;
    TGraphErrors *fSavedGrScanGeom = nullptr;

    void SaveAllToRoot();

    void SetBetheError(const string& message) {
        cerr << "\n[ERRORE BETHE-BLOCH] " << message << "\n" << endl;
        for (int i = 0; i < 4; ++i) calc[i]->SetText("N/A");
        calcframe->Layout();
    }

    void SetMCError(const string &message);
    bool RunMCSimulation(MuonMCResult &out, string &error);
    void DrawScene3D(TPad *pad, const MuonMCResult &r);
    void DrawSummaryTableAndHist(TPad *pad, const MuonMCResult &r);
    void DrawMCTab1(const MuonMCResult &r);
    void DrawMCTab2(const MuonMCResult &r);

public:
    UnifiedDetectorGUI(const TGWindow *p, UInt_t w, UInt_t h);
    virtual ~UnifiedDetectorGUI();

    void DrawCalc();
    void DrawMultiLayerCalc(int charge, double mass, const string& partName);
    bool ComputeMultiLayerChain(int charge, double mass, vector<LayerData>& layers, string& errorMsg);
    bool ResolveMaterial(int id, Medium& mat, string& error);
    void UpdateInputMode(Int_t id);
    void UpdateLayerVisibility();
    void OnPartSelected(Int_t id);
    void OnMatSelected(Int_t id);
    void ToggleMultiLayerMode();
    void ToggleCompoundBuilder();
    void ToggleCosmicPanel();
    void SyncDetectorHeightFromLastLayer();
    void ComputeCosmicRatesAndEff();
    void SyncPMTToMC();
    void BuildAndAddCompound();
    void AddMaterialEntry(const string &name, int id);

    void DoCalculateMC();
    void AutoScaleIdealPlane();
    void DoExit();
};

// =========================================================================
// SALVATAGGIO COMPLETO DI TUTTI I GRAFICI NEL FILE ROOT
// =========================================================================
void UnifiedDetectorGUI::SaveAllToRoot() {
    TFile rootOut("bethe_mc_results.root", "RECREATE");
    if (rootOut.IsZombie()) return;

    TDirectory *dirBethe = rootOut.mkdir("BetheBloch");
    dirBethe->cd();
    if (fSavedGrBethe) fSavedGrBethe->Write("gr_Bethe_StoppingPower", TObject::kOverwrite);
    if (fSavedGrBragg) fSavedGrBragg->Write("gr_Bethe_BraggCurve", TObject::kOverwrite);

    TDirectory *dirCosmic = rootOut.mkdir("Spettro_Rastin_PMT");
    dirCosmic->cd();
    if (fSavedGrRastin)     fSavedGrRastin->Write("gr_Rastin_Spectrum", TObject::kOverwrite);
    if (fSavedHistEdepAll)  fSavedHistEdepAll->Write("h_PMT_Edep_All", TObject::kOverwrite);
    if (fSavedHistEdepLost) fSavedHistEdepLost->Write("h_PMT_Edep_Lost", TObject::kOverwrite);
    if (fSavedHistEdepPass) fSavedHistEdepPass->Write("h_PMT_Edep_Detected", TObject::kOverwrite);

    TDirectory *dirMC = rootOut.mkdir("MonteCarlo_3D");
    dirMC->cd();
    if (fSavedHistMCRates)  fSavedHistMCRates->Write("h_MC_Rates_Categorie", TObject::kOverwrite);
    if (fSavedHistAnaRates) fSavedHistAnaRates->Write("h_Ana_Rates_Categorie", TObject::kOverwrite);
    if (fSavedGrScanAcc)    fSavedGrScanAcc->Write("gr_MC_Accettanza_vs_Lgen", TObject::kOverwrite);
    if (fSavedGrScanGeom)   fSavedGrScanGeom->Write("gr_MC_FattoreGeom_vs_Lgen", TObject::kOverwrite);

    rootOut.cd();
    if (fCanvasBethe && fCanvasBethe->GetCanvas()) fCanvasBethe->GetCanvas()->Write("Canvas_1_Bethe_Spettro_PMT", TObject::kOverwrite);
    if (fCanvasMC1 && fCanvasMC1->GetCanvas()) fCanvasMC1->GetCanvas()->Write("Canvas_2_MC_3D_Rate", TObject::kOverwrite);
    if (fCanvasMC2 && fCanvasMC2->GetCanvas()) fCanvasMC2->GetCanvas()->Write("Canvas_3_MC_Scan_Piano_S", TObject::kOverwrite);

    rootOut.Close();
}

// =========================================================================
// COSTRUTTORE INTERFACCIA GRAFICA UNIFICATA
// =========================================================================
UnifiedDetectorGUI::UnifiedDetectorGUI(const TGWindow *p, UInt_t w, UInt_t h) {
    fMain = new TGMainFrame(p, w, h, kHorizontalFrame);
    fMain->Connect("CloseWindow()", "UnifiedDetectorGUI", this, "DoExit()");
    fMain->DontCallClose();

    // =====================================================================
    // 1. PANNELLO SINISTRO CON SCHEDE
    // =====================================================================
    fLeftTabs = new TGTab(fMain, 420, 790);

    TGCompositeFrame *tabLeftBethe = fLeftTabs->AddTab("Bethe-Bloch & PMT");
    leftBetheFrame = new TGVerticalFrame(tabLeftBethe, 405, 760);
    tabLeftBethe->AddFrame(leftBetheFrame, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));

    TGLabel *incPart = new TGLabel(leftBetheFrame, "Particella incidente:");
    leftBetheFrame->AddFrame(incPart, new TGLayoutHints(kLHintsTop, 4, 0, 4, 0));
    partList = new TGComboBox(leftBetheFrame);
    partList->AddEntry("Definita dall'utente", 0);
    for (size_t i = 0; i < particles.size(); ++i) partList->AddEntry(particles[i].name.c_str(), i + 1);
    partList->Select(particles.empty() ? 0 : 1);
    partList->Resize(380, 22);
    partList->Connect("Selected(Int_t)", "UnifiedDetectorGUI", this, "OnPartSelected(Int_t)");
    leftBetheFrame->AddFrame(partList, new TGLayoutHints(kLHintsTop, 4, 4, 2, 2));

    userdefinedpart = new TGGroupFrame(leftBetheFrame, "Parametri particella:", kVerticalFrame);
    leftBetheFrame->AddFrame(userdefinedpart, new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 4, 4, 3, 2));
    {
        TGHorizontalFrame *rCh = new TGHorizontalFrame(userdefinedpart);
        rCh->AddFrame(new TGLabel(rCh, "Carica |z| [unita' di e]:"), new TGLayoutHints(kLHintsCenterY | kLHintsLeft));
        userCharge = new TGNumberEntry(rCh, 1, 9, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive, TGNumberFormat::kNELLimitMinMax, 1, 1000);
        userCharge->Resize(95, 20);
        rCh->AddFrame(userCharge, new TGLayoutHints(kLHintsCenterY | kLHintsRight));
        userdefinedpart->AddFrame(rCh, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));

        TGHorizontalFrame *rMs = new TGHorizontalFrame(userdefinedpart);
        rMs->AddFrame(new TGLabel(rMs, "Massa [MeV/c^2]:"), new TGLayoutHints(kLHintsCenterY | kLHintsLeft));
        userMass = new TGNumberEntry(rMs, kMuonMassMeV, 9, -1, TGNumberFormat::kNESRealTwo, TGNumberFormat::kNEAPositive, TGNumberFormat::kNELLimitMin, 5.12);
        userMass->Resize(115, 20);
        rMs->AddFrame(userMass, new TGLayoutHints(kLHintsCenterY | kLHintsRight));
        userdefinedpart->AddFrame(rMs, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));
    }

    btnToggleMultiMode = new TGTextButton(leftBetheFrame, "Mostra Calcolo Multi-Strato (Inverso)");
    btnToggleMultiMode->Connect("Clicked()", "UnifiedDetectorGUI", this, "ToggleMultiLayerMode()");
    leftBetheFrame->AddFrame(btnToggleMultiMode, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 4, 4, 4, 2));

    mater = new TGLabel(leftBetheFrame, "Materiale singolo:");
    leftBetheFrame->AddFrame(mater, new TGLayoutHints(kLHintsTop, 4, 0, 3, 0));
    matList = new TGComboBox(leftBetheFrame);
    matList->AddEntry("Definito dall'utente", 0);
    for (size_t i = 0; i < materials.size(); ++i) matList->AddEntry(materials[i].name.c_str(), i + 1);
    matList->Select(materials.empty() ? 0 : 1);
    matList->Resize(380, 22);
    matList->Connect("Selected(Int_t)", "UnifiedDetectorGUI", this, "OnMatSelected(Int_t)");
    leftBetheFrame->AddFrame(matList, new TGLayoutHints(kLHintsTop, 4, 4, 2, 2));

    userdefinedmat = new TGGroupFrame(leftBetheFrame, "Parametri materiale singolo:", kVerticalFrame);
    leftBetheFrame->AddFrame(userdefinedmat, new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 4, 4, 3, 2));
    {
        auto addMatRow = [](TGGroupFrame *grp, const char *txt, double val, TGNumberFormat::EStyle st, double minV) -> TGNumberEntry* {
            TGHorizontalFrame *r = new TGHorizontalFrame(grp);
            r->AddFrame(new TGLabel(r, txt), new TGLayoutHints(kLHintsCenterY | kLHintsLeft));
            TGNumberEntry *ne = new TGNumberEntry(r, val, 9, -1, st, TGNumberFormat::kNEAPositive, TGNumberFormat::kNELLimitMin, minV);
            ne->Resize(105, 20);
            r->AddFrame(ne, new TGLayoutHints(kLHintsCenterY | kLHintsRight));
            grp->AddFrame(r, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));
            return ne;
        };
        userZ = addMatRow(userdefinedmat, "Z effettivo:", 6.5, TGNumberFormat::kNESRealOne, 0.01);
        userA = addMatRow(userdefinedmat, "A efficace [g/mol]:", 13.0, TGNumberFormat::kNESRealThree, 0.01);
        userIP = addMatRow(userdefinedmat, "I medio [eV]:", 64.7, TGNumberFormat::kNESRealOne, 1.0);
        userDens = addMatRow(userdefinedmat, "Densita' [g/cm^3]:", 0.86, TGNumberFormat::kNESRealThree, 1e-6);
    }

    multiLayerFrame = new TGGroupFrame(leftBetheFrame, "Parametri strati multipli", kVerticalFrame);
    leftBetheFrame->AddFrame(multiLayerFrame, new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 4, 4, 3, 2));
    TGHorizontalFrame *hfLayers = new TGHorizontalFrame(multiLayerFrame);
    hfLayers->AddFrame(new TGLabel(hfLayers, "Numero di strati: "), new TGLayoutHints(kLHintsCenterY | kLHintsLeft));
    numLayersEntry = new TGNumberEntry(hfLayers, 2, 3, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive, TGNumberFormat::kNELLimitMinMax, 1, MAX_LAYERS);
    numLayersEntry->Connect("ValueSet(Long_t)", "UnifiedDetectorGUI", this, "UpdateLayerVisibility()");
    hfLayers->AddFrame(numLayersEntry, new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 5, 0, 0, 0));
    multiLayerFrame->AddFrame(hfLayers, new TGLayoutHints(kLHintsTop, 0, 0, 2, 4));

    for (int i = 0; i < MAX_LAYERS; ++i) {
        layerFrames[i] = new TGHorizontalFrame(multiLayerFrame, 340, 24);
        multiLayerFrame->AddFrame(layerFrames[i], new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0, 0, 1, 1));
        layerFrames[i]->AddFrame(new TGLabel(layerFrames[i], Form("S%d:", i + 1)), new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 0, 4, 0, 0));
        layerMatList[i] = new TGComboBox(layerFrames[i]);
        layerMatList[i]->AddEntry("User Def", 0);
        for (size_t m = 0; m < materials.size(); ++m) layerMatList[i]->AddEntry(materials[m].name.c_str(), m + 1);
        layerMatList[i]->Select(materials.empty() ? 0 : 1);
        layerMatList[i]->Resize(165, 20);
        layerFrames[i]->AddFrame(layerMatList[i], new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 0, 2, 0, 0));
        layerFrames[i]->AddFrame(new TGLabel(layerFrames[i], "cm"), new TGLayoutHints(kLHintsCenterY | kLHintsRight, 2, 2, 0, 0));
        layerThickness[i] = new TGNumberEntry(layerFrames[i], 20.0, 6, -1, TGNumberFormat::kNESRealOne, TGNumberFormat::kNEAPositive, TGNumberFormat::kNELLimitMinMax, 0.001, 1.0e5);
        layerThickness[i]->Resize(70, 20);
        layerThickness[i]->Connect("ValueSet(Long_t)", "UnifiedDetectorGUI", this, "SyncDetectorHeightFromLastLayer()");
        layerFrames[i]->AddFrame(layerThickness[i], new TGLayoutHints(kLHintsCenterY | kLHintsRight, 2, 2, 0, 0));
    }

    TGTextButton *compBtn = new TGTextButton(leftBetheFrame, "Mostra/Nascondi Costruttore Composti");
    compBtn->Connect("Clicked()", "UnifiedDetectorGUI", this, "ToggleCompoundBuilder()");
    leftBetheFrame->AddFrame(compBtn, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 4, 4, 4, 2));

    TGTextButton *cosmicBtn = new TGTextButton(leftBetheFrame, "Mostra/Nascondi Rate Cosmici & Efficienza PMT");
    cosmicBtn->Connect("Clicked()", "UnifiedDetectorGUI", this, "ToggleCosmicPanel()");
    leftBetheFrame->AddFrame(cosmicBtn, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 4, 4, 2, 3));

    modeLabel = new TGLabel(leftBetheFrame, "Input cinematico:");
    leftBetheFrame->AddFrame(modeLabel, new TGLayoutHints(kLHintsTop, 4, 0, 2, 0));
    inputModeList = new TGComboBox(leftBetheFrame);
    inputModeList->AddEntry("Energia cinetica T [MeV]", 0);
    inputModeList->AddEntry("Impulso p [MeV/c]", 1);
    inputModeList->Select(0);
    inputModeList->Resize(380, 22);
    inputModeList->Connect("Selected(Int_t)", "UnifiedDetectorGUI", this, "UpdateInputMode(Int_t)");
    leftBetheFrame->AddFrame(inputModeList, new TGLayoutHints(kLHintsTop, 4, 4, 2, 2));

    TGHorizontalFrame *rKinVal = new TGHorizontalFrame(leftBetheFrame);
    inputValueLabel = new TGLabel(rKinVal, "Valore T [MeV]:");
    rKinVal->AddFrame(inputValueLabel, new TGLayoutHints(kLHintsCenterY | kLHintsLeft));
    inputValue = new TGNumberEntry(rKinVal, 1000.0, 12, -1, TGNumberFormat::kNESRealTwo, TGNumberFormat::kNEAPositive, TGNumberFormat::kNELLimitMinMax, 0.001, 1.0e7);
    inputValue->Resize(130, 22);
    rKinVal->AddFrame(inputValue, new TGLayoutHints(kLHintsCenterY | kLHintsRight));
    leftBetheFrame->AddFrame(rKinVal, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

    TGHorizontalFrame *buttonframe = new TGHorizontalFrame(leftBetheFrame, 360, 36);
    leftBetheFrame->AddFrame(buttonframe, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 4, 4, 4, 2));
    TGTextButton *draw = new TGTextButton(buttonframe, "&Calcola e Disegna Bethe-Bloch");
    draw->SetBackgroundColor(TColor::RGB2Pixel(205, 235, 255));
    draw->Connect("Clicked()", "UnifiedDetectorGUI", this, "DrawCalc()");
    buttonframe->AddFrame(draw, new TGLayoutHints(kLHintsExpandX, 2, 4, 2, 2));
    TGTextButton *exitBtnBethe = new TGTextButton(buttonframe, "&Esci");
    exitBtnBethe->SetBackgroundColor(TColor::RGB2Pixel(255, 210, 210));
    exitBtnBethe->Connect("Clicked()", "UnifiedDetectorGUI", this, "DoExit()");
    buttonframe->AddFrame(exitBtnBethe, new TGLayoutHints(kLHintsRight, 2, 2, 2, 2));

    TGGroupFrame *resultframe = new TGGroupFrame(leftBetheFrame, "Risultati Bethe-Bloch / Multi-Strato", kHorizontalFrame);
    resultframe->SetTitlePos(TGGroupFrame::kCenter);
    leftBetheFrame->AddFrame(resultframe, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 4, 4, 4, 4));

    TGVerticalFrame *answerframe = new TGVerticalFrame(resultframe, 10, 10);
    resultframe->AddFrame(answerframe, new TGLayoutHints(kLHintsLeft));

    calcframe = new TGVerticalFrame(resultframe, 100, 10);
    resultframe->AddFrame(calcframe, new TGLayoutHints(kLHintsRight | kLHintsExpandX));

    const char *ansLabels[4] = {
        "Stopping iniziale [MeV/cm]: ",
        "Stopping massico [MeV cm^2/g]: ",
        "Range CSDA stimato [cm]: ",
        "Range CSDA stimato [g/cm^2]: "
    };
    for (int i = 0; i < 4; ++i) {
        answerlabel[i] = new TGLabel(answerframe, ansLabels[i]);
        answerframe->AddFrame(answerlabel[i], new TGLayoutHints(kLHintsLeft, 0, 0, 2, 2));
        calc[i] = new TGLabel(calcframe, "N/A");
        calcframe->AddFrame(calc[i], new TGLayoutHints(kLHintsRight, 0, 0, 2, 2));
    }

    // ---------------------------------------------------------------------
    // SCHEDA SINISTRA 2: MONTE CARLO GEOMETRICO
    // ---------------------------------------------------------------------
    TGCompositeFrame *tabLeftMC = fLeftTabs->AddTab("MC Coincidenza");
    TGVerticalFrame *leftMCFrame = new TGVerticalFrame(tabLeftMC, 405, 760);
    tabLeftMC->AddFrame(leftMCFrame, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));

    auto addNumberRow = [](TGGroupFrame *parent, const char *labelTxt, double initVal,
                           double minVal, double maxVal,
                           TGNumberFormat::EStyle style) -> TGNumberEntry* {
        TGHorizontalFrame *row = new TGHorizontalFrame(parent);
        row->AddFrame(new TGLabel(row, labelTxt),
                      new TGLayoutHints(kLHintsCenterY | kLHintsLeft));
        TGNumberEntry *ne = new TGNumberEntry(
            row, initVal, 9, -1, style,
            TGNumberFormat::kNEAPositive, TGNumberFormat::kNELLimitMinMax, minVal, maxVal);
        ne->Resize(92, 20);
        row->AddFrame(ne, new TGLayoutHints(kLHintsCenterY | kLHintsRight, 4, 2, 1, 1));
        parent->AddFrame(row, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));
        return ne;
    };

    TGGroupFrame *geomFrame = new TGGroupFrame(leftMCFrame, "1. Blocco SUPERIORE (il resto e' fisso nel codice)", kVerticalFrame);
    leftMCFrame->AddFrame(geomFrame, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));
    fX = addNumberRow(geomFrame, "X / larghezza [cm]:", 31.0, 0.1, 1.0e5, TGNumberFormat::kNESRealTwo);
    fY = addNumberRow(geomFrame, "Y / profondita' [cm]:", 30.0, 0.1, 1.0e5, TGNumberFormat::kNESRealTwo);
    fZ = addNumberRow(geomFrame, "Z / altezza (Top) [cm]:", 33.0, 0.1, 1.0e5, TGNumberFormat::kNESRealTwo);

    TGGroupFrame *mcFrame = new TGGroupFrame(leftMCFrame, "2. Piano Ideale di Riferimento (S = Lx * Ly) e MC", kVerticalFrame);
    leftMCFrame->AddFrame(mcFrame, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));
    fLx = addNumberRow(mcFrame, "Lx piano ideale [cm] (>> X):", 300.0, 1.0, 1.0e6, TGNumberFormat::kNESRealOne);
    fLy = addNumberRow(mcFrame, "Ly piano ideale [cm] (>> Y):", 300.0, 1.0, 1.0e6, TGNumberFormat::kNESRealOne);
    fNEvents = addNumberRow(mcFrame, "Numero eventi MC (N_gen):", 1000000, 100, 2.0e8, TGNumberFormat::kNESInteger);
    fSeed = addNumberRow(mcFrame, "Seed generatore TRandom3:", 42, 0, 9999999, TGNumberFormat::kNESInteger);

    TGTextButton *autoBtn = new TGTextButton(mcFrame, "Imposta Lx, Ly ottimali");
    autoBtn->Connect("Clicked()", "UnifiedDetectorGUI", this, "AutoScaleIdealPlane()");
    mcFrame->AddFrame(autoBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

    TGGroupFrame *fluxFrame = new TGGroupFrame(leftMCFrame, "3. Parametri Flusso Cosmico & Efficienze", kVerticalFrame);
    leftMCFrame->AddFrame(fluxFrame, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));
    fI0          = addNumberRow(fluxFrame, "I_0 verticale [muoni/(cm^2 s sr)]:", 9.13e-3, 1.0e-7, 10.0, TGNumberFormat::kNESReal);
    fExponent    = addNumberRow(fluxFrame, "Esponente angolare n (cos^n):", 2.0, 0.0, 10.0, TGNumberFormat::kNESRealTwo);
    fEfficiency1 = addNumberRow(fluxFrame, "Efficienza S1 (Top) [%]:", 100.0, 0.0, 100.0, TGNumberFormat::kNESRealTwo);
    fEfficiency2 = addNumberRow(fluxFrame, "Efficienza S2 (Mid) [%]:", 100.0, 0.0, 100.0, TGNumberFormat::kNESRealTwo);
    fEfficiency3 = addNumberRow(fluxFrame, "Efficienza S3 (Bot) [%]:", 100.0, 0.0, 100.0, TGNumberFormat::kNESRealTwo);

    TGHorizontalFrame *mcBtnRow = new TGHorizontalFrame(leftMCFrame);
    leftMCFrame->AddFrame(mcBtnRow, new TGLayoutHints(kLHintsExpandX, 4, 4, 3, 3));
    TGTextButton *runMC = new TGTextButton(mcBtnRow, "&Esegui MC 3D e Calcola Rate");
    runMC->SetBackgroundColor(TColor::RGB2Pixel(205, 235, 255));
    runMC->Connect("Clicked()", "UnifiedDetectorGUI", this, "DoCalculateMC()");
    mcBtnRow->AddFrame(runMC, new TGLayoutHints(kLHintsExpandX, 0, 4, 0, 0));

    TGTextButton *exitBtnMC = new TGTextButton(mcBtnRow, "&Esci");
    exitBtnMC->SetBackgroundColor(TColor::RGB2Pixel(255, 210, 210));
    exitBtnMC->Connect("Clicked()", "UnifiedDetectorGUI", this, "DoExit()");
    mcBtnRow->AddFrame(exitBtnMC, new TGLayoutHints(kLHintsRight, 0, 0, 0, 0));

    TGGroupFrame *resFrame = new TGGroupFrame(leftMCFrame, "Risultati Accettanza, Fattore Geometrico e Rate", kVerticalFrame);
    leftMCFrame->AddFrame(resFrame, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

    const char *labels[16] = {
        "Piano ideale S [cm^2] (v_piano): ",
        "Accettanza TOTALE A_tot [%]: ",
        "  (a) Coincidenza 3-fold [%]: ",
        "  (b) - vuoto -  [%]: ",
        "  (c) - vuoto -  [%]: ",
        "  (d) - vuoto -  [%]: ",
        "  (e) Mancati (miss) [%]: ",
        "Fattore Geom. G_MC [cm^2 sr]: ",
        "Fattore Geom. G_Ana [cm^2 sr] (solo Top): ",
        "Rate MC (a) Coincidenza 3F [Hz]: ",
        "Rate MC (b) - vuoto -  [Hz]: ",
        "Rate MC (c) - vuoto -  [Hz]: ",
        "Rate MC (d) - vuoto -  [Hz]: ",
        "Rate MC Faccia Sup. (a+b) [Hz]: ",
        "Rate MC 4 Pareti Lat. (c+d) [Hz]: ",
        "Rate MC TOTALE Reale (Coinc) [Hz]: "
    };
    for (int i = 0; i < 16; ++i) {
        TGHorizontalFrame *rr = new TGHorizontalFrame(resFrame);
        rr->AddFrame(new TGLabel(rr, labels[i]), new TGLayoutHints(kLHintsCenterY | kLHintsLeft));
        fRes[i] = new TGLabel(rr, "--");
        rr->AddFrame(fRes[i], new TGLayoutHints(kLHintsCenterY | kLHintsRight));
        resFrame->AddFrame(rr, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 0));
    }

    // =====================================================================
    // 2. PANNELLO A SCOMPARSA: COSTRUTTORE COMPOSTI
    // =====================================================================
    centerframe = new TGVerticalFrame(fMain, 290, 40);
    TGGroupFrame *compGroup = new TGGroupFrame(centerframe, "Costruisci Composto", kVerticalFrame);
    centerframe->AddFrame(compGroup, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 5, 5, 5, 5));

    TGHorizontalFrame *rName = new TGHorizontalFrame(compGroup);
    rName->AddFrame(new TGLabel(rName, "Nome: "), new TGLayoutHints(kLHintsCenterY | kLHintsLeft));
    fCompName = new TGTextEntry(rName, "Acqua (H2O)");
    fCompName->Resize(175, 22);
    rName->AddFrame(fCompName, new TGLayoutHints(kLHintsCenterY | kLHintsRight, 2, 2, 2, 2));
    compGroup->AddFrame(rName, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

    TGHorizontalFrame *rDens = new TGHorizontalFrame(compGroup);
    rDens->AddFrame(new TGLabel(rDens, "Rho [g/cm^3]: "), new TGLayoutHints(kLHintsCenterY | kLHintsLeft));
    fCompDens = new TGNumberEntry(rDens, 1.0, 7, -1, TGNumberFormat::kNESRealThree, TGNumberFormat::kNEAPositive);
    fCompDens->Resize(95, 22);
    rDens->AddFrame(fCompDens, new TGLayoutHints(kLHintsCenterY | kLHintsRight, 2, 2, 2, 2));
    compGroup->AddFrame(rDens, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

    TGHorizontalFrame *rI = new TGHorizontalFrame(compGroup);
    rI->AddFrame(new TGLabel(rI, "I [eV] (0=Bragg): "), new TGLayoutHints(kLHintsCenterY | kLHintsLeft));
    fCompI = new TGNumberEntry(rI, 75.0, 7, -1, TGNumberFormat::kNESRealOne, TGNumberFormat::kNEANonNegative);
    fCompI->Resize(95, 22);
    rI->AddFrame(fCompI, new TGLayoutHints(kLHintsCenterY | kLHintsRight, 2, 2, 2, 2));
    compGroup->AddFrame(rI, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

    compGroup->AddFrame(new TGLabel(compGroup, "Componenti (dal file materials.txt):"), new TGLayoutHints(kLHintsTop | kLHintsLeft, 2, 2, 10, 2));

    int defH = 0, defO = 0;
    for (size_t i = 0; i < materials.size(); ++i) {
        if (materials[i].name == "H" || materials[i].name == "Hydrogen") defH = i + 1;
        if (materials[i].name == "O" || materials[i].name == "Oxygen") defO = i + 1;
    }

    for (int i = 0; i < MAX_COMPONENTS; i++) {
        TGHorizontalFrame *rEl = new TGHorizontalFrame(compGroup);
        fCompElem[i] = new TGComboBox(rEl);
        fCompElem[i]->AddEntry("--", 0);
        for (size_t j = 0; j < materials.size(); ++j) fCompElem[i]->AddEntry(materials[j].name.c_str(), j + 1);
        fCompElem[i]->Resize(115, 22);
        fCompCount[i] = new TGNumberEntry(rEl, 0.0, 6, -1, TGNumberFormat::kNESRealTwo, TGNumberFormat::kNEANonNegative);
        fCompCount[i]->Resize(75, 22);

        if (i == 0) { fCompElem[i]->Select(defH); fCompCount[i]->SetNumber(2.0); }
        if (i == 1) { fCompElem[i]->Select(defO); fCompCount[i]->SetNumber(1.0); }
        if (i > 1) { fCompElem[i]->Select(0); fCompCount[i]->SetNumber(0.0); }

        rEl->AddFrame(fCompElem[i], new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 2, 5, 2, 2));
        rEl->AddFrame(fCompCount[i], new TGLayoutHints(kLHintsCenterY | kLHintsRight, 2, 2, 2, 2));
        compGroup->AddFrame(rEl, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));
    }

    TGTextButton *btnBuild = new TGTextButton(compGroup, "&Aggiungi e Stampa Info su Terminale");
    btnBuild->Connect("Clicked()", "UnifiedDetectorGUI", this, "BuildAndAddCompound()");
    compGroup->AddFrame(btnBuild, new TGLayoutHints(kLHintsTop | kLHintsCenterX, 5, 5, 15, 5));

    // =====================================================================
    // 3. PANNELLO A SCOMPARSA: RATE COSMICI & EFFICIENZA PMT 3D
    // =====================================================================
    cosmicframe = new TGVerticalFrame(fMain, 310, 40);
    TGGroupFrame *cosmicGroup = new TGGroupFrame(cosmicframe, "Rate Cosmici & Efficienza PMT", kVerticalFrame);
    cosmicframe->AddFrame(cosmicGroup, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 5, 5, 5, 5));

    cosmicGroup->AddFrame(new TGLabel(cosmicGroup, "Geometria rivelatore (ultimo strato SN):"), new TGLayoutHints(kLHintsTop | kLHintsLeft, 2, 2, 4, 2));

    TGHorizontalFrame *rLx = new TGHorizontalFrame(cosmicGroup);
    rLx->AddFrame(new TGLabel(rLx, "Larghezza Lx [cm]: "), new TGLayoutHints(kLHintsCenterY | kLHintsLeft));
    fDetLx = new TGNumberEntry(rLx, 31.0, 7, -1, TGNumberFormat::kNESRealOne, TGNumberFormat::kNEAPositive);
    fDetLx->Resize(85, 20);
    rLx->AddFrame(fDetLx, new TGLayoutHints(kLHintsCenterY | kLHintsRight));
    cosmicGroup->AddFrame(rLx, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));

    TGHorizontalFrame *rLy = new TGHorizontalFrame(cosmicGroup);
    rLy->AddFrame(new TGLabel(rLy, "Lunghezza Ly [cm]: "), new TGLayoutHints(kLHintsCenterY | kLHintsLeft));
    fDetLy = new TGNumberEntry(rLy, 30.0, 7, -1, TGNumberFormat::kNESRealOne, TGNumberFormat::kNEAPositive);
    fDetLy->Resize(85, 20);
    rLy->AddFrame(fDetLy, new TGLayoutHints(kLHintsCenterY | kLHintsRight));
    cosmicGroup->AddFrame(rLy, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));

    TGHorizontalFrame *rLz = new TGHorizontalFrame(cosmicGroup);
    rLz->AddFrame(new TGLabel(rLz, "Altezza Lz (=SN) [cm]: "), new TGLayoutHints(kLHintsCenterY | kLHintsLeft));
    fDetLz = new TGNumberEntry(rLz, 33.0, 7, -1, TGNumberFormat::kNESRealOne, TGNumberFormat::kNEAPositive);
    fDetLz->Resize(85, 20);
    rLz->AddFrame(fDetLz, new TGLayoutHints(kLHintsCenterY | kLHintsRight));
    cosmicGroup->AddFrame(rLz, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 4));

    TGHorizontalFrame *rThr = new TGHorizontalFrame(cosmicGroup);
    rThr->AddFrame(new TGLabel(rThr, "Soglia PMT [MeV]: "), new TGLayoutHints(kLHintsCenterY | kLHintsLeft));
    fThresholdMeV = new TGNumberEntry(rThr, 35.0, 7, -1, TGNumberFormat::kNESRealTwo, TGNumberFormat::kNEANonNegative, TGNumberFormat::kNELLimitMinMax, 0.0, 1.0e6);
    fThresholdMeV->Resize(90, 22);
    fThresholdMeV->Connect("ValueSet(Long_t)", "UnifiedDetectorGUI", this, "ComputeCosmicRatesAndEff()");
    rThr->AddFrame(fThresholdMeV, new TGLayoutHints(kLHintsCenterY | kLHintsRight));
    cosmicGroup->AddFrame(rThr, new TGLayoutHints(kLHintsExpandX, 2, 2, 6, 4));

    TGTextButton *btnCalcCosmic = new TGTextButton(cosmicGroup, "&Calcola Rate e Efficienza Spettro");
    btnCalcCosmic->SetBackgroundColor(TColor::RGB2Pixel(205, 235, 255));
    btnCalcCosmic->Connect("Clicked()", "UnifiedDetectorGUI", this, "ComputeCosmicRatesAndEff()");
    cosmicGroup->AddFrame(btnCalcCosmic, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 2, 6, 4));

    TGTextButton *btnSyncToMC = new TGTextButton(cosmicGroup, "Invia Geometria & Efficienza al MC 3D");
    btnSyncToMC->Connect("Clicked()", "UnifiedDetectorGUI", this, "SyncPMTToMC()");
    cosmicGroup->AddFrame(btnSyncToMC, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 2, 2, 8));

    TGGroupFrame *cosmicResGroup = new TGGroupFrame(cosmicGroup, "Risultati Spettro & Rate", kVerticalFrame);
    cosmicGroup->AddFrame(cosmicResGroup, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 2, 4, 4));

    fCosmicAccLabel      = new TGLabel(cosmicResGroup, "Accettanza geom.: N/A");
    fCosmicRateStopLabel = new TGLabel(cosmicResGroup, "Rate Stop (decad.): N/A");
    fCosmicRatePassLabel = new TGLabel(cosmicResGroup, "Rate Passanti: N/A");
    fCosmicRateTotLabel  = new TGLabel(cosmicResGroup, "Rate Totale Rivelata: N/A");
    fCosmicEffLabel      = new TGLabel(cosmicResGroup, "Efficienza PMT: N/A");

    cosmicResGroup->AddFrame(fCosmicAccLabel,      new TGLayoutHints(kLHintsTop | kLHintsLeft, 2, 2, 2, 2));
    cosmicResGroup->AddFrame(fCosmicRateStopLabel, new TGLayoutHints(kLHintsTop | kLHintsLeft, 2, 2, 2, 2));
    cosmicResGroup->AddFrame(fCosmicRatePassLabel, new TGLayoutHints(kLHintsTop | kLHintsLeft, 2, 2, 2, 2));
    cosmicResGroup->AddFrame(fCosmicRateTotLabel,  new TGLayoutHints(kLHintsTop | kLHintsLeft, 2, 2, 2, 2));
    cosmicResGroup->AddFrame(fCosmicEffLabel,      new TGLayoutHints(kLHintsTop | kLHintsLeft, 2, 2, 4, 3));

    // =====================================================================
    // 4. AREA GRAFICA DESTRA CON 3 SCHEDE
    // =====================================================================
    fRightTabs = new TGTab(fMain, 1020, 790);

    TGCompositeFrame *rt1 = fRightTabs->AddTab("1. Bethe-Bloch & Spettro PMT");
    fCanvasBethe = new TRootEmbeddedCanvas("CanvasBethe", rt1, 1010, 760);
    rt1->AddFrame(fCanvasBethe, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));

    TGCompositeFrame *rt2 = fRightTabs->AddTab("2. MC Accettanza 3D Telescopio");
    fCanvasMC1 = new TRootEmbeddedCanvas("CanvasMC1", rt2, 1010, 760);
    rt2->AddFrame(fCanvasMC1, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));

    TGCompositeFrame *rt3 = fRightTabs->AddTab("3. Scan Piano Generato S");
    fCanvasMC2 = new TRootEmbeddedCanvas("CanvasMC2", rt3, 1010, 760);
    rt3->AddFrame(fCanvasMC2, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));

    // =====================================================================
    // 5. ASSEMBLAGGIO FINALE FINESTRA
    // =====================================================================
    OnPartSelected(partList->GetSelected());
    OnMatSelected(matList->GetSelected());
    UpdateLayerVisibility();

    fMain->AddFrame(fLeftTabs,   new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 3, 3, 3, 3));
    fMain->AddFrame(centerframe, new TGLayoutHints(kLHintsLeft, 2, 2, 2, 2));
    fMain->AddFrame(cosmicframe, new TGLayoutHints(kLHintsLeft, 2, 2, 2, 2));
    fMain->AddFrame(fRightTabs,  new TGLayoutHints(kLHintsRight | kLHintsExpandX | kLHintsExpandY, 4, 4, 4, 4));

    fMain->SetWindowName("Simulatore Telescopio Cosmico - Bethe-Bloch & MC Accettanza");
    fMain->MapSubwindows();

    fMultiLayerMode = false;
    fBuilderVisible = false;
    fCosmicVisible  = false;
    leftBetheFrame->HideFrame(multiLayerFrame);
    fMain->HideFrame(centerframe);
    fMain->HideFrame(cosmicframe);

    leftBetheFrame->Layout();
    fMain->Layout();
    fMain->Resize(fMain->GetDefaultSize());
    fMain->MapWindow();

    DoCalculateMC();
    DrawCalc();
    fLeftTabs->SetTab(0);
    fRightTabs->SetTab(0);
}

UnifiedDetectorGUI::~UnifiedDetectorGUI() {
    delete fSavedGrBethe;
    delete fSavedGrBragg;
    delete fSavedGrRastin;
    delete fSavedHistEdepAll;
    delete fSavedHistEdepLost;
    delete fSavedHistEdepPass;
    delete fSavedHistMCRates;
    delete fSavedHistAnaRates;
    delete fSavedGrScanAcc;
    delete fSavedGrScanGeom;
    if (fMain) {
        fMain->Cleanup();
        delete fMain;
        fMain = nullptr;
    }
}

void UnifiedDetectorGUI::DoExit() {
    std::cout << "Chiusura dell'applicazione. Tutti i grafici sono salvati in 'bethe_mc_results.root'." << std::endl;
    gApplication->Terminate(0);
}

// =========================================================================
// METODI SEZIONE BETHE-BLOCH, COMPOSTI & SPETTRO PMT
// =========================================================================
void UnifiedDetectorGUI::OnPartSelected(Int_t id) {
    if (id == 0) {
        userCharge->SetState(kTRUE);
        userMass->SetState(kTRUE);
    } else if (id > 0 && id <= (Int_t)particles.size()) {
        userCharge->SetIntNumber(particles[id - 1].charge);
        userMass->SetNumber(particles[id - 1].massMeV);
        userCharge->SetState(kFALSE);
        userMass->SetState(kFALSE);
    }
}

void UnifiedDetectorGUI::OnMatSelected(Int_t id) {
    if (id == 0) {
        userZ->SetState(kTRUE);
        userA->SetState(kTRUE);
        userIP->SetState(kTRUE);
        userDens->SetState(kTRUE);
    } else if (id > 0 && id <= (Int_t)materials.size()) {
        const Medium& m = materials[id - 1];
        userZ->SetNumber(m.zEffective);
        userA->SetNumber(m.zEffective / m.zOverA);
        userIP->SetNumber(m.ionEV);
        userDens->SetNumber(m.density);
        userZ->SetState(kFALSE);
        userA->SetState(kFALSE);
        userIP->SetState(kFALSE);
        userDens->SetState(kFALSE);
    }
}

void UnifiedDetectorGUI::ToggleMultiLayerMode() {
    fMultiLayerMode = !fMultiLayerMode;
    if (fMultiLayerMode) {
        btnToggleMultiMode->SetText("Torna al Calcolo Singolo Materiale");
        leftBetheFrame->HideFrame(mater);
        leftBetheFrame->HideFrame(matList);
        leftBetheFrame->HideFrame(userdefinedmat);
        leftBetheFrame->HideFrame(modeLabel);
        leftBetheFrame->HideFrame(inputModeList);
        leftBetheFrame->HideFrame(inputValueLabel);
        leftBetheFrame->HideFrame(inputValue);
        leftBetheFrame->ShowFrame(multiLayerFrame);
        UpdateLayerVisibility();
        const char* labels[4] = {
            "Energia iniziale richiesta [MeV]: ",
            "Impulso iniziale richiesto [MeV/c]: ",
            "Energia uscita S1 [MeV]: ",
            "Impulso uscita S1 [MeV/c]: "
        };
        for (int i = 0; i < 4; ++i) answerlabel[i]->SetText(labels[i]);
        SyncDetectorHeightFromLastLayer();
    } else {
        btnToggleMultiMode->SetText("Mostra Calcolo Multi-Strato (Inverso)");
        leftBetheFrame->ShowFrame(mater);
        leftBetheFrame->ShowFrame(matList);
        leftBetheFrame->ShowFrame(userdefinedmat);
        leftBetheFrame->ShowFrame(modeLabel);
        leftBetheFrame->ShowFrame(inputModeList);
        leftBetheFrame->ShowFrame(inputValueLabel);
        leftBetheFrame->ShowFrame(inputValue);
        leftBetheFrame->HideFrame(multiLayerFrame);
        const char* labels[4] = {
            "Stopping iniziale [MeV/cm]: ",
            "Stopping massico [MeV cm^2/g]: ",
            "Range CSDA stimato [cm]: ",
            "Range CSDA stimato [g/cm^2]: "
        };
        for (int i = 0; i < 4; ++i) answerlabel[i]->SetText(labels[i]);
        UpdateInputMode(inputModeList->GetSelected());
    }
    leftBetheFrame->Layout();
    calcframe->Layout();
    if (fMain) fMain->Layout();
}

void UnifiedDetectorGUI::ToggleCompoundBuilder() {
    fBuilderVisible = !fBuilderVisible;
    if (fBuilderVisible) fMain->ShowFrame(centerframe);
    else fMain->HideFrame(centerframe);
    fMain->Layout();
}

void UnifiedDetectorGUI::ToggleCosmicPanel() {
    fCosmicVisible = !fCosmicVisible;
    if (fCosmicVisible) {
        SyncDetectorHeightFromLastLayer();
        fMain->ShowFrame(cosmicframe);
    } else {
        fMain->HideFrame(cosmicframe);
    }
    fMain->Layout();
}

void UnifiedDetectorGUI::SyncDetectorHeightFromLastLayer() {
    int n = numLayersEntry->GetIntNumber();
    if (n >= 1 && n <= MAX_LAYERS && fDetLz) {
        fDetLz->SetNumber(layerThickness[n - 1]->GetNumber());
    }
}

void UnifiedDetectorGUI::SyncPMTToMC() {
    fX->SetNumber(fDetLx->GetNumber());
    fY->SetNumber(fDetLy->GetNumber());
    fZ->SetNumber(fDetLz->GetNumber());
    fEfficiency1->SetNumber(fLastPMTEfficiencyPct);
    fLeftTabs->SetTab(1);
    DoCalculateMC();
}

void UnifiedDetectorGUI::BuildAndAddCompound() {
    vector<int> idx;
    vector<double> cnt;
    for (int r = 0; r < MAX_COMPONENTS; ++r) {
        int id = fCompElem[r]->GetSelected();
        double n = fCompCount[r]->GetNumber();
        if (id > 0 && n > 0.0) {
            int e = id - 1;
            bool merged = false;
            for (size_t k = 0; k < idx.size(); ++k) {
                if (idx[k] == e) { cnt[k] += n; merged = true; break; }
            }
            if (!merged) { idx.push_back(e); cnt.push_back(n); }
        }
    }

    CompoundInfo info;
    string err;
    if (!BuildCompound(fCompName->GetText(), idx, cnt, fCompDens->GetNumber(), fCompI->GetNumber(), info, err)) {
        SetBetheError("Errore costruttore composto: " + err);
        return;
    }

    cout << DescribeCompound(info) << endl;

    int pos = -1;
    for (size_t i = 0; i < materials.size(); ++i) {
        if (materials[i].name == info.medium.name) { pos = static_cast<int>(i); break; }
    }
    if (pos >= 0) {
        materials[pos] = info.medium;
        matList->Select(pos + 1);
        cout << "=> INFO: Materiale '" << info.medium.name << "' aggiornato nella lista.\n" << endl;
    } else {
        materials.push_back(info.medium);
        AddMaterialEntry(info.medium.name, materials.size());
        cout << "=> SUCCESS: Nuovo materiale '" << info.medium.name << "' aggiunto alla lista!\n" << endl;
    }
}

void UnifiedDetectorGUI::AddMaterialEntry(const string &name, int id) {
    matList->AddEntry(name.c_str(), id);
    for (int i = 0; i < MAX_LAYERS; ++i) layerMatList[i]->AddEntry(name.c_str(), id);
    for (int i = 0; i < MAX_COMPONENTS; ++i) fCompElem[i]->AddEntry(name.c_str(), id);
    matList->GetListBox()->Layout();
    for (int i = 0; i < MAX_LAYERS; ++i) layerMatList[i]->GetListBox()->Layout();
    matList->Select(id);
}

void UnifiedDetectorGUI::UpdateLayerVisibility() {
    int n = numLayersEntry->GetIntNumber();
    for (int i = 0; i < MAX_LAYERS; ++i) {
        if (i < n) multiLayerFrame->ShowFrame(layerFrames[i]);
        else multiLayerFrame->HideFrame(layerFrames[i]);
    }
    SyncDetectorHeightFromLastLayer();
    multiLayerFrame->Layout();
    if (!fMultiLayerMode) leftBetheFrame->HideFrame(multiLayerFrame);
    leftBetheFrame->Layout();
    if (fMain) fMain->Layout();
}

bool UnifiedDetectorGUI::ResolveMaterial(int matID, Medium& mat, string& error) {
    if (matID == 0) {
        mat.name = "User Definito";
        mat.zEffective = userZ->GetNumber();
        const double A = userA->GetNumber();
        mat.ionEV = userIP->GetNumber();
        mat.density = userDens->GetNumber();
        if (!(A > 0.0) || !(mat.zEffective > 0.0) || !(mat.ionEV > 0.0) || !(mat.density > 0.0) ||
            !std::isfinite(A) || !std::isfinite(mat.zEffective) || !std::isfinite(mat.ionEV) || !std::isfinite(mat.density)) {
            error = "parametri del materiale definito dall'utente non validi.";
            return false;
        }
        mat.zOverA = mat.zEffective / A;
        BuildApproxDensityParams(mat);
        return true;
    }
    if (matID > 0 && static_cast<size_t>(matID) <= materials.size()) {
        mat = materials[static_cast<size_t>(matID - 1)];
        return true;
    }
    error = "selezione del materiale non valida.";
    return false;
}

void UnifiedDetectorGUI::UpdateInputMode(Int_t id) {
    if (fMultiLayerMode) inputValueLabel->SetText("Input ignorato (multi-strato)");
    else if (id == 1) inputValueLabel->SetText("Valore p [MeV/c]:");
    else inputValueLabel->SetText("Valore T [MeV]:");
    fMain->Layout();
}

bool UnifiedDetectorGUI::ComputeMultiLayerChain(int charge, double mass, vector<LayerData>& layers, string& errorMsg) {
    int nLayers = numLayersEntry->GetIntNumber();
    if (nLayers < 1 || nLayers > MAX_LAYERS) {
        errorMsg = "numero di strati non valido.";
        return false;
    }
    layers.resize(nLayers);
    for (int i = 0; i < nLayers; ++i) {
        layers[i].thicknessCm = layerThickness[i]->GetNumber();
        if (!(layers[i].thicknessCm > 0.0) || !std::isfinite(layers[i].thicknessCm)) {
            errorMsg = Form("spessore non valido per S%d.", i + 1);
            return false;
        }
        string error;
        if (!ResolveMaterial(layerMatList[i]->GetSelected(), layers[i].mat, error)) {
            errorMsg = Form("S%d: %s", i + 1, error.c_str());
            return false;
        }
        layers[i].arealThickness = layers[i].mat.density * layers[i].thicknessCm;
    }
    double currentKineticOut = 0.0;
    for (int i = nLayers - 1; i >= 0; --i) {
        layers[i].kineticOut = currentKineticOut;
        layers[i].momentumOut = MomentumFromKinetic(currentKineticOut, mass);
        double rangeOut = 0.0;
        if (currentKineticOut > 0.0 && !ComputeRangeMassCSDA(currentKineticOut, mass, charge, layers[i].mat, rangeOut)) {
            errorMsg = Form("errore range uscita S%d.", i + 1);
            return false;
        }
        double kineticIn = 0.0;
        string error;
        if (!FindKineticForRange(rangeOut + layers[i].arealThickness, mass, charge, layers[i].mat, kineticIn, error)) {
            errorMsg = Form("S%d: %s", i + 1, error.c_str());
            return false;
        }
        layers[i].kineticIn = kineticIn;
        layers[i].momentumIn = MomentumFromKinetic(kineticIn, mass);
        currentKineticOut = kineticIn;
    }
    return true;
}

void UnifiedDetectorGUI::DrawMultiLayerCalc(int charge, double mass, const string& partName) {
    vector<LayerData> layers;
    string err;
    if (!ComputeMultiLayerChain(charge, mass, layers, err)) {
        SetBetheError(err);
        return;
    }

    int nLayers = static_cast<int>(layers.size());
    calc[0]->SetText(Form("%.7g", layers[0].kineticIn));
    calc[1]->SetText(Form("%.7g", layers[0].momentumIn));
    calc[2]->SetText(nLayers > 1 ? Form("%.7g", layers[0].kineticOut) : "0");
    calc[3]->SetText(nLayers > 1 ? Form("%.7g", layers[0].momentumOut) : "0");
    calcframe->Layout();

    double totThickCm = 0.0, totAreal = 0.0;
    for (const auto& ly : layers) {
        totThickCm += ly.thicknessCm;
        totAreal += ly.arealThickness;
    }
    std::cout << "\n================================================================================\n"
              << " CALCOLO INVERSO CSDA MULTI-STRATO (" << nLayers << " strati in successione)\n"
              << "================================================================================\n"
              << "  Particella incidente   : " << partName << "  (|z| = " << std::abs(charge)
              << ", massa = " << G(mass, 7) << " MeV/c^2)\n"
              << "  Spessore totale        : " << G(totThickCm, 6) << " cm  (" << G(totAreal, 6) << " g/cm^2)\n"
              << "--------------------------------------------------------------------------------\n"
              << "  Strato | Materiale          | L [cm]  | rho*L [g/cm^2] | T_in [MeV]  | T_out [MeV] | dE [MeV]\n"
              << "  ------------------------------------------------------------------------------\n";
    for (int i = 0; i < nLayers; ++i) {
        double dE = layers[i].kineticIn - layers[i].kineticOut;
        std::cout << "  S" << setw(2) << left << (i + 1) << right
                  << "    | " << setw(18) << left << layers[i].mat.name.substr(0, 18) << right
                  << " | " << setw(7) << fixed << setprecision(2) << layers[i].thicknessCm
                  << " | " << setw(14) << setprecision(3) << layers[i].arealThickness
                  << " | " << setw(11) << setprecision(3) << layers[i].kineticIn
                  << " | " << setw(11) << setprecision(3) << layers[i].kineticOut
                  << " | " << setw(8) << setprecision(3) << dE << "\n";
    }
    std::cout << "--------------------------------------------------------------------------------\n"
              << "  RISULTATO IN INGRESSO A S1:\n"
              << "    -> Energia cinetica minima T_in(S1) : " << G(layers[0].kineticIn, 7) << " MeV\n"
              << "    -> Impulso minimo          p_in(S1) : " << G(layers[0].momentumIn, 7) << " MeV/c\n"
              << "================================================================================\n" << std::endl;

    fRightTabs->SetTab(0);
    TCanvas *canvas = fCanvasBethe->GetCanvas();
    canvas->Clear();
    canvas->SetFillColor(10);
    canvas->cd();

    TPaveText *summary = new TPaveText(0.03, 0.03, 0.97, 0.97, "NDC");
    summary->SetBit(kCanDelete);
    summary->SetFillColor(0);
    summary->SetTextAlign(12);
    summary->SetTextSize(nLayers > 6 ? 0.022 : 0.028);
    summary->AddText(Form("Calcolo inverso CSDA: %d strati in successione", nLayers));
    summary->AddText(Form("Particella: %s; |z| = %d; massa = %.7g MeV/c^2", partName.c_str(), std::abs(charge), mass));
    summary->AddText(" ");
    for (int i = 0; i < nLayers; ++i) {
        summary->AddText(Form("S%d: %s (L = %.4g cm, rho = %.4g g/cm^3)", i + 1, layers[i].mat.name.c_str(), layers[i].thicknessCm, layers[i].mat.density));
        summary->AddText(Form("   In: T = %.7g MeV  --->  Out: T = %.7g MeV", layers[i].kineticIn, layers[i].kineticOut));
    }
    summary->AddText(" ");
    summary->AddText(Form("Energia cinetica iniziale: %.7g MeV", layers[0].kineticIn));
    summary->Draw();
    canvas->Modified();
    canvas->Update();

    SaveAllToRoot();
}

void UnifiedDetectorGUI::ComputeCosmicRatesAndEff() {
    int charge = 1;
    double mass = kMuonMassMeV;
    const int partID = partList->GetSelected();
    if (partID == 0) {
        charge = userCharge->GetIntNumber();
        mass = userMass->GetNumber();
    } else if (partID > 0 && static_cast<size_t>(partID) <= particles.size()) {
        charge = particles[partID - 1].charge;
        mass = particles[partID - 1].massMeV;
    }

    vector<LayerData> layers;
    string err;
    if (!ComputeMultiLayerChain(charge, mass, layers, err)) {
        fCosmicAccLabel->SetText("Errore Multi-Strato!");
        fCosmicEffLabel->SetText(err.c_str());
        cosmicframe->Layout();
        return;
    }

    const int nLayers = static_cast<int>(layers.size());
    const double soglia_MeV = fThresholdMeV->GetNumber();

    const double Lx = fDetLx->GetNumber();
    const double Ly = fDetLy->GetNumber();
    const double Lz = fDetLz->GetNumber();
    const double Area_Top = Lx * Ly;
    const double acceptance_Analytic = (kPi / 2.0) * Area_Top + (kPi / 4.0) * Lz * (Lx + Ly);

    vector<FastRangeTable> rangeTables(nLayers);
    for (int i = 0; i < nLayers; ++i) {
        rangeTables[i].Build(mass, charge, layers[i].mat);
    }

    double E_min_enter_vert = 0.0;
    if (nLayers > 1) {
        double kin = 0.0;
        for (int i = nLayers - 2; i >= 0; --i) {
            double rOut = (kin > 0.0) ? rangeTables[i].RangeFromT(kin) : 0.0;
            kin = rangeTables[i].TFromRange(rOut + layers[i].arealThickness);
        }
        E_min_enter_vert = kin;
    }
    const double E_max_stop_vert = layers[0].kineticIn;

    const double mGeV = mass / 1000.0;
    const double pMinGrid = 0.02;
    const double pMaxGrid = 1000.0;
    const double totalCosmicFlux = IntegralRastinAboveP(pMinGrid);

    const int N_mc = 120000;
    const double margin = 4.0 * Lz;
    const double genLx = Lx + 2.0 * margin;
    const double genLy = Ly + 2.0 * margin;
    const double Area_gen = genLx * genLy;

    const int nCdf = 1500;
    vector<double> pGrid(nCdf), cdfGrid(nCdf);
    pGrid[0] = pMinGrid;
    cdfGrid[0] = 0.0;
    for (int i = 1; i < nCdf; ++i) {
        pGrid[i] = pMinGrid * std::pow(pMaxGrid / pMinGrid, static_cast<double>(i) / (nCdf - 1));
        cdfGrid[i] = cdfGrid[i - 1] + IntegrateRastinMomentum(pGrid[i - 1], pGrid[i]);
    }
    cdfGrid[nCdf - 1] = totalCosmicFlux;

    double maxDiag = std::sqrt(Lx * Lx + Ly * Ly + Lz * Lz);
    double maxPlotE = std::max(soglia_MeV * 1.6, 2.5 * maxDiag);
    if (maxPlotE > 600.0) maxPlotE = 600.0;

    TH1D *hEdepAll  = new TH1D("hEdepAll",  "Deposito Energetico Totale nel Rivelatore (pesato su Spettro Cosmico);Energia depositata #DeltaE [MeV];Rate differenziale [Hz / MeV]", 160, 0.0, maxPlotE);
    TH1D *hEdepLost = new TH1D("hEdepLost", "Deposito Energetico Sotto Soglia PMT;Energia depositata #DeltaE [MeV];Rate differenziale [Hz / MeV]", 160, 0.0, maxPlotE);
    TH1D *hEdepPass = new TH1D("hEdepPass", "Deposito Energetico Sopra Soglia PMT;Energia depositata #DeltaE [MeV];Rate differenziale [Hz / MeV]", 160, 0.0, maxPlotE);
    hEdepAll->SetBit(kCanDelete);
    hEdepLost->SetBit(kCanDelete);
    hEdepPass->SetBit(kCanDelete);
    hEdepAll->SetDirectory(nullptr);
    hEdepLost->SetDirectory(nullptr);
    hEdepPass->SetDirectory(nullptr);

    TRandom3 rng(42);
    long N_trials = N_mc;
    double w_reach_det = 0.0, w_stop_tot  = 0.0, w_pass_tot  = 0.0, w_stop_det  = 0.0, w_pass_det  = 0.0;
    const double dGamma_per_trial = (Area_gen * (kPi / 2.0)) / static_cast<double>(N_trials);

    for (long i = 0; i < N_trials; ++i) {
        double x0 = rng.Uniform(-0.5 * genLx, 0.5 * genLx);
        double y0 = rng.Uniform(-0.5 * genLy, 0.5 * genLy);
        double z0 = Lz;

        double u = rng.Uniform(0.0, 1.0);
        double cosTheta = std::pow(u, 0.25);
        double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
        double phi = rng.Uniform(0.0, 2.0 * kPi);

        double dx = sinTheta * std::cos(phi);
        double dy = sinTheta * std::sin(phi);
        double dz = -cosTheta;

        double tMin = 0.0, tMax = 1.0e30;
        auto intersectSlab = [&](double p0, double d, double bMin, double bMax) -> bool {
            if (std::fabs(d) < 1e-12) return (p0 >= bMin && p0 <= bMax);
            double t1 = (bMin - p0) / d;
            double t2 = (bMax - p0) / d;
            if (t1 > t2) std::swap(t1, t2);
            tMin = std::max(tMin, t1);
            tMax = std::min(tMax, t2);
            return tMin <= tMax;
        };

        if (!(intersectSlab(x0, dx, -0.5 * Lx, 0.5 * Lx) &&
              intersectSlab(y0, dy, -0.5 * Ly, 0.5 * Ly) &&
              intersectSlab(z0, dz, 0.0, Lz))) continue;

        double trackLenCm = tMax - tMin;
        if (trackLenCm <= 1e-6) continue;

        double rFlux = rng.Uniform(0.0, totalCosmicFlux);
        auto it = std::lower_bound(cdfGrid.begin(), cdfGrid.end(), rFlux);
        size_t idx = std::distance(cdfGrid.begin(), it);
        double pGeV = pMaxGrid;
        if (idx == 0) pGeV = pGrid[0];
        else if (idx < cdfGrid.size()) {
            double denom = cdfGrid[idx] - cdfGrid[idx - 1];
            double frac = (denom > 0.0) ? (rFlux - cdfGrid[idx - 1]) / denom : 0.0;
            pGeV = pGrid[idx - 1] + frac * (pGrid[idx] - pGrid[idx - 1]);
        }

        double T_MeV = (std::sqrt(pGeV * pGeV + mGeV * mGeV) - mGeV) * 1000.0;

        bool stoppedAbove = false;
        double invCos = 1.0 / std::max(0.05, cosTheta);
        for (int s = 0; s < nLayers - 1; ++s) {
            double dX_above = layers[s].arealThickness * invCos;
            double dE_above = rangeTables[s].EnergyLoss(T_MeV, dX_above);
            T_MeV -= dE_above;
            if (T_MeV <= 0.0) { stoppedAbove = true; break; }
        }
        if (stoppedAbove || T_MeV <= 0.0) continue;

        w_reach_det += 1.0;

        double dX_det = layers[nLayers - 1].mat.density * trackLenCm;
        double R_in = rangeTables[nLayers - 1].RangeFromT(T_MeV);
        bool stopsInDet = (R_in <= dX_det);
        double dE_det = stopsInDet ? T_MeV : rangeTables[nLayers - 1].EnergyLoss(T_MeV, dX_det);

        double eventRateWeight = totalCosmicFlux * dGamma_per_trial;
        hEdepAll->Fill(dE_det, eventRateWeight);

        if (stopsInDet) {
            w_stop_tot += 1.0;
            if (dE_det > soglia_MeV) {
                w_stop_det += 1.0;
                hEdepPass->Fill(dE_det, eventRateWeight);
            } else {
                hEdepLost->Fill(dE_det, eventRateWeight);
            }
        } else {
            w_pass_tot += 1.0;
            if (dE_det > soglia_MeV) {
                w_pass_det += 1.0;
                hEdepPass->Fill(dE_det, eventRateWeight);
            } else {
                hEdepLost->Fill(dE_det, eventRateWeight);
            }
        }
    }

    double scaleToHz = totalCosmicFlux * dGamma_per_trial;
    double rate_stop_tot = w_stop_tot * scaleToHz;
    double rate_pass_tot = w_pass_tot * scaleToHz;
    double rate_all_tot  = w_reach_det * scaleToHz;
    double rate_stop_det = w_stop_det * scaleToHz;
    double rate_pass_det = w_pass_det * scaleToHz;
    double rate_all_det  = (w_stop_det + w_pass_det) * scaleToHz;
    double efficienza = (w_reach_det > 0.0) ? ((w_stop_det + w_pass_det) / w_reach_det) : 0.0;
    fLastPMTEfficiencyPct = efficienza * 100.0;

    double binW = hEdepAll->GetBinWidth(1);
    if (binW > 0.0) {
        hEdepAll->Scale(1.0 / binW);
        hEdepLost->Scale(1.0 / binW);
        hEdepPass->Scale(1.0 / binW);
    }

    std::cout << std::fixed << setprecision(4)
              << "\n================================================================================\n"
              << " SPETTRO COSMICI (RASTIN 1984a) & EFFICIENZA PMT 3D\n"
              << "================================================================================\n"
              << "  Rivelatore attivo (S" << nLayers << " - " << layers[nLayers - 1].mat.name << ") : "
              << setprecision(2) << Lx << " x " << Ly << " x " << Lz << " cm\n"
              << "  Accettanza geometrica analitica  : " << setprecision(2) << acceptance_Analytic << " cm^2 sr\n"
              << "  Flusso integrale Rastin (>20 MeV/c)   : " << setprecision(5) << totalCosmicFlux << " cm^-2 s^-1 sr^-1\n"
              << "  Finestra verticale stop in S" << nLayers << "    : T_in(S1) = "
              << setprecision(2) << E_min_enter_vert << " .. " << E_max_stop_vert << " MeV\n"
              << "  Soglia energetica PMT impostata  : > " << setprecision(2) << soglia_MeV << " MeV\n"
              << "--------------------------------------------------------------------------------\n"
              << "  Categoria Eventi          | Rivelati (> soglia)       | Incidenti Totali\n"
              << "  ------------------------------------------------------------------------------\n"
              << "  1) Muoni che si fermano   | " << setw(8) << setprecision(4) << rate_stop_det << " Hz (" << setw(6) << setprecision(1) << rate_stop_det * 60.0 << " /min)"
              << " | " << setw(8) << setprecision(4) << rate_stop_tot << " Hz (" << setw(6) << setprecision(1) << rate_stop_tot * 60.0 << " /min)\n"
              << "  2) Muoni passanti         | " << setw(8) << setprecision(4) << rate_pass_det << " Hz (" << setw(6) << setprecision(1) << rate_pass_det * 60.0 << " /min)"
              << " | " << setw(8) << setprecision(4) << rate_pass_tot << " Hz (" << setw(6) << setprecision(1) << rate_pass_tot * 60.0 << " /min)\n"
              << "  3) TOTALE nel rivelatore  | " << setw(8) << setprecision(4) << rate_all_det  << " Hz (" << setw(6) << setprecision(1) << rate_all_det * 60.0  << " /min)"
              << " | " << setw(8) << setprecision(4) << rate_all_tot  << " Hz (" << setw(6) << setprecision(1) << rate_all_tot * 60.0  << " /min)\n"
              << "--------------------------------------------------------------------------------\n"
              << "  => EFFICIENZA DI RIVELAZIONE PMT : " << setprecision(2) << fLastPMTEfficiencyPct << " %\n"
              << "================================================================================\n" << std::endl;

    fCosmicAccLabel->SetText(Form("Accettanza geom.: %.1f cm^2 sr", acceptance_Analytic));
    fCosmicRateStopLabel->SetText(Form("Rate Stop: %.3f / %.3f Hz", rate_stop_det, rate_stop_tot));
    fCosmicRatePassLabel->SetText(Form("Rate Pass: %.3f / %.3f Hz", rate_pass_det, rate_pass_tot));
    fCosmicRateTotLabel->SetText(Form("Rate Tot:  %.3f / %.3f Hz", rate_all_det, rate_all_tot));
    fCosmicEffLabel->SetText(Form("Efficienza PMT: %.2f %%", fLastPMTEfficiencyPct));
    cosmicframe->Layout();

    fRightTabs->SetTab(0);
    TCanvas *canvas = fCanvasBethe->GetCanvas();
    canvas->Clear();
    canvas->SetFillColor(10);

    TPad *pad1 = new TPad("pad1", "", 0.01, 0.51, 0.99, 0.99, 10);
    TPad *pad2 = new TPad("pad2", "", 0.01, 0.01, 0.99, 0.49, 10);
    pad1->SetBit(kCanDelete);
    pad2->SetBit(kCanDelete);
    pad1->SetLeftMargin(0.11);
    pad1->SetRightMargin(0.04);
    pad1->SetBottomMargin(0.12);
    pad2->SetLeftMargin(0.11);
    pad2->SetRightMargin(0.04);
    pad2->SetBottomMargin(0.12);
    pad1->Draw();
    pad2->Draw();

    pad1->cd();
    pad1->SetLogx();
    pad1->SetLogy();
    pad1->SetGrid();

    const int n_tab = 18;
    double p_tab[n_tab] = {0.35, 0.40, 0.50, 0.60, 0.80, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0, 1000.0, 1500.0, 2000.0, 3000.0};
    double flux_p_tab[n_tab] = {2.85e-3, 2.90e-3, 2.94e-3, 2.92e-3, 2.80e-3, 2.62e-3, 1.69e-3, 5.17e-4, 1.35e-4, 2.58e-5, 2.11e-6, 2.77e-7, 3.12e-8, 1.45e-9, 1.31e-10, 3.14e-11, 1.13e-11, 2.67e-12};
    const int n_pts = n_tab + 3;
    double Ek[n_pts], flux_E[n_pts];
    double p_low[3] = {0.05, 0.176, 0.265};

    for (int i = 0; i < n_pts; i++) {
        double p = (i < 3) ? p_low[i] : p_tab[i - 3];
        double f_p = (i < 3) ? flux_p_tab[0] : flux_p_tab[i - 3];
        Ek[i] = std::sqrt(p * p + mGeV * mGeV) - mGeV;
        double beta = p / (Ek[i] + mGeV);
        flux_E[i] = f_p / beta;
    }

    TGraph *gr = new TGraph(n_pts, Ek, flux_E);
    gr->SetBit(kCanDelete);
    gr->SetTitle("Spettro Differenziale Muoni Primari in Ingresso a S1 (Rastin 1984a);Energia Cinetica T_{in,S1} [GeV];dN/dE_{k} [cm^{-2} s^{-1} sr^{-1} GeV^{-1}]");
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.6);
    gr->SetMarkerColor(kAzure + 7);
    gr->SetLineColor(kAzure + 7);
    gr->SetLineWidth(2);
    gr->GetXaxis()->SetLimits(0.01, 4000.0);
    gr->GetYaxis()->SetRangeUser(1e-12, 1e-1);
    gr->GetXaxis()->SetTitleSize(0.045);
    gr->GetYaxis()->SetTitleSize(0.045);
    gr->Draw("APL");

    delete fSavedGrRastin;
    fSavedGrRastin = static_cast<TGraph*>(gr->Clone("gr_Rastin_Spectrum"));

    delete fSavedHistEdepAll;
    fSavedHistEdepAll = static_cast<TH1D*>(hEdepAll->Clone("h_PMT_Edep_All"));
    fSavedHistEdepAll->SetDirectory(nullptr);

    delete fSavedHistEdepLost;
    fSavedHistEdepLost = static_cast<TH1D*>(hEdepLost->Clone("h_PMT_Edep_Lost"));
    fSavedHistEdepLost->SetDirectory(nullptr);

    delete fSavedHistEdepPass;
    fSavedHistEdepPass = static_cast<TH1D*>(hEdepPass->Clone("h_PMT_Edep_Detected"));
    fSavedHistEdepPass->SetDirectory(nullptr);

    double E_enter_GeV = std::max(0.01, E_min_enter_vert / 1000.0);
    double E_stop_GeV  = std::max(E_enter_GeV * 1.02, E_max_stop_vert / 1000.0);
    TBox *boxStop = new TBox(E_enter_GeV, 1e-12, E_stop_GeV, 1e-1);
    boxStop->SetBit(kCanDelete);
    boxStop->SetFillColorAlpha(kOrange + 1, 0.25);
    boxStop->SetLineColor(kOrange + 2);
    boxStop->Draw("SAME");
    gr->Draw("PL SAME");

    TLegend *leg1 = new TLegend(0.42, 0.62, 0.94, 0.88);
    leg1->SetBit(kCanDelete);
    leg1->SetFillColor(0);
    leg1->SetBorderSize(1);
    leg1->SetTextSize(0.034);
    leg1->AddEntry(gr, "Spettro primario dN/dE_{k} (Rastin + plateau basse E)", "pl");
    leg1->AddEntry(boxStop, Form("Stop in S%d: T_{in,S1} = %.1f - %.1f MeV", nLayers, E_min_enter_vert, E_max_stop_vert), "f");
    if (nLayers > 1) leg1->AddEntry((TObject*)0, Form("(T_{out,S1} al limite stop = %.1f MeV)", layers[0].kineticOut), "");
    leg1->Draw();

    pad2->cd();
    pad2->SetTopMargin(0.12);
    pad2->SetGrid();

    double yMax2 = std::max(hEdepLost->GetMaximum(), hEdepPass->GetMaximum()) * 1.25;
    if (yMax2 <= 0.0) yMax2 = 1.0;

    hEdepLost->SetTitle("Deposito Energetico nel Rivelatore (pesato su Spettro Cosmico);Energia depositata #DeltaE [MeV];Rate differenziale [Hz / MeV]");
    hEdepLost->SetLineColor(kOrange + 7);
    hEdepLost->SetLineWidth(2);
    hEdepLost->SetFillColorAlpha(kOrange + 7, 0.35);
    hEdepLost->GetXaxis()->SetTitleSize(0.045);
    hEdepLost->GetYaxis()->SetTitleSize(0.045);
    hEdepLost->GetYaxis()->SetRangeUser(0.0, yMax2);
    hEdepLost->SetStats(0);
    hEdepLost->Draw("HIST");

    hEdepPass->SetLineColor(kGreen + 2);
    hEdepPass->SetLineWidth(2);
    hEdepPass->SetFillColorAlpha(kGreen + 2, 0.40);
    hEdepPass->Draw("HIST SAME");

    TLine *lineThr = new TLine(soglia_MeV, 0.0, soglia_MeV, yMax2);
    lineThr->SetBit(kCanDelete);
    lineThr->SetLineColor(kRed + 1);
    lineThr->SetLineWidth(3);
    lineThr->SetLineStyle(2);
    lineThr->Draw("SAME");

    TLegend *leg2 = new TLegend(0.56, 0.56, 0.94, 0.85);
    leg2->SetBit(kCanDelete);
    leg2->SetFillColor(0);
    leg2->SetBorderSize(1);
    leg2->SetTextSize(0.032);
    leg2->AddEntry(hEdepLost, Form("Persi (#leq %.1f MeV)", soglia_MeV), "f");
    leg2->AddEntry(hEdepPass, Form("Rilevati (> %.1f MeV)", soglia_MeV), "f");
    leg2->AddEntry(lineThr,   Form("Soglia = %.1f MeV", soglia_MeV), "l");
    leg2->Draw();

    TPaveText *pt = new TPaveText(0.56, 0.22, 0.94, 0.53, "NDC");
    pt->SetBit(kCanDelete);
    pt->SetFillColor(0);
    pt->SetBorderSize(1);
    pt->SetTextAlign(12);
    pt->SetTextSize(0.032);
    pt->AddText(Form("Accettanza (analitica singola): %.1f cm^{2} sr", acceptance_Analytic));
    pt->AddText(Form("Rate Stop: %.3f Hz (tot %.3f Hz)", rate_stop_det, rate_stop_tot));
    pt->AddText(Form("Rate Pass: %.3f Hz (tot %.3f Hz)", rate_pass_det, rate_pass_tot));
    pt->AddText(Form("Efficienza PMT: %.2f %%", fLastPMTEfficiencyPct));
    pt->Draw();

    canvas->Modified();
    canvas->Update();

    SaveAllToRoot();
}

void UnifiedDetectorGUI::DrawCalc() {
    int charge = 1;
    double mass = kMuonMassMeV;
    string partName = "Particella user";
    const int partID = partList->GetSelected();
    if (partID == 0) {
        charge = userCharge->GetIntNumber();
        mass = userMass->GetNumber();
    } else if (partID > 0 && static_cast<size_t>(partID) <= particles.size()) {
        const IncomingPart& p = particles[partID - 1];
        charge = p.charge;
        mass = p.massMeV;
        partName = p.name;
    }
    if (charge == 0 || !(mass > 10.0 * kElectronMassMeV) || !std::isfinite(mass)) {
        SetBetheError("carica/massa non valida.");
        return;
    }
    if (fMultiLayerMode) {
        DrawMultiLayerCalc(charge, mass, partName);
        return;
    }
    const double rawInput = inputValue->GetNumber();
    if (!(rawInput > 0.0) || !std::isfinite(rawInput)) {
        SetBetheError("input non valido.");
        return;
    }
    double kinetic = rawInput;
    if (inputModeList->GetSelected() == 1) kinetic = std::sqrt(rawInput * rawInput + mass * mass) - mass;
    if (!(kinetic > 0.0) || !std::isfinite(kinetic)) {
        SetBetheError("T risultante non valida.");
        return;
    }

    Medium mat;
    string materialError;
    if (!ResolveMaterial(matList->GetSelected(), mat, materialError)) {
        SetBetheError(materialError);
        return;
    }

    if (std::sqrt(kinetic * (kinetic + 2.0 * mass)) / mass < kEtaCut * (1.0 - 1.0e-12)) {
        SetBetheError("T sotto il cutoff.");
        return;
    }
    Results result{};
    result.massStopping = ComputeBetheBloch(kinetic, mass, charge, mat);
    if (!(result.massStopping > 0.0)) {
        SetBetheError("stopping iniziale nullo.");
        return;
    }
    result.stopping = result.massStopping * mat.density;

    const double tCut = KineticAtEtaCut(mass);
    vector<double> depth, stoppingPlot;
    depth.reserve(kRangeSteps + 1);
    stoppingPlot.reserve(kRangeSteps + 1);
    double massDepth = 0.0, previousT = kinetic, previousS = result.massStopping;
    depth.push_back(0.0);
    stoppingPlot.push_back(previousS);

    for (int i = 1; i < kRangeSteps; ++i) {
        const double currentT = kinetic * std::pow(tCut / kinetic, static_cast<double>(i) / (kRangeSteps - 1));
        const double currentS = ComputeBetheBloch(currentT, mass, charge, mat);
        if (!(currentS > 0.0) || !std::isfinite(currentS)) {
            SetBetheError("interruzione integrazione.");
            return;
        }
        massDepth += 0.5 * (previousT - currentT) * (1.0 / previousS + 1.0 / currentS);
        depth.push_back(massDepth);
        stoppingPlot.push_back(currentS);
        previousT = currentT;
        previousS = currentS;
    }
    result.integratedMass = massDepth;
    const double sCut = ComputeBetheBloch(tCut, mass, charge, mat);
    result.tailMass = (sCut > 0.0) ? tCut / sCut : 0.0;
    result.rangeMass = result.integratedMass + result.tailMass;
    result.rangeCm = result.rangeMass / mat.density;

    const double momentum = MomentumFromKinetic(kinetic, mass);
    const double eta = momentum / mass;
    const double gamma = std::sqrt(1.0 + eta * eta);
    const double beta = eta / gamma;
    std::cout << "\n================================================================================\n"
              << " CALCOLO BETHE-BLOCH & RANGE CSDA (Singolo Materiale)\n"
              << "================================================================================\n"
              << "  Particella incidente   : " << partName << "  (|z| = " << std::abs(charge)
              << ", massa = " << G(mass, 7) << " MeV/c^2)\n"
              << "  Materiale bersaglio    : " << mat.name << "  (Z_eff = " << G(mat.zEffective, 4)
              << ", Z/A = " << G(mat.zOverA, 5) << ", I = " << G(mat.ionEV, 5)
              << " eV, rho = " << G(mat.density, 5) << " g/cm^3)\n"
              << "  Cinematica iniziale    : T = " << G(kinetic, 7) << " MeV | p = " << G(momentum, 7)
              << " MeV/c | beta = " << G(beta, 5) << " | gamma = " << G(gamma, 5) << "\n"
              << "--------------------------------------------------------------------------------\n"
              << "  Stopping lineare (-dE/dx)       : " << G(result.stopping, 6)     << " MeV/cm\n"
              << "  Stopping massico (-dE/rho*dx)   : " << G(result.massStopping, 6) << " MeV cm^2/g\n"
              << "  Range CSDA massico (integrato)  : " << G(result.integratedMass, 6) << " g/cm^2  (coda sotto cutoff: " << G(result.tailMass, 4) << " g/cm^2)\n"
              << "  Range CSDA massico TOTALE       : " << G(result.rangeMass, 6)    << " g/cm^2\n"
              << "  Range CSDA lineare TOTALE       : " << G(result.rangeCm, 6)      << " cm\n"
              << "================================================================================\n" << std::endl;

    vector<double> ePlot, dedxPlot;
    ePlot.reserve(kEnergyPlotPoints);
    dedxPlot.reserve(kEnergyPlotPoints);
    double minS = 1e300, maxS = 0.0, ePlotMax = std::max(10000.0, kinetic * 1.5);
    for (int i = 0; i < kEnergyPlotPoints; ++i) {
        const double T = tCut * std::pow(ePlotMax / tCut, static_cast<double>(i) / (kEnergyPlotPoints - 1));
        const double S = ComputeBetheBloch(T, mass, charge, mat);
        if (S > 0.0 && std::isfinite(S)) {
            ePlot.push_back(T);
            dedxPlot.push_back(S);
            minS = std::min(minS, S);
            maxS = std::max(maxS, S);
        }
    }
    if (ePlot.size() < 2 || depth.size() < 2) return;

    fRightTabs->SetTab(0);
    TCanvas *canvas = fCanvasBethe->GetCanvas();
    canvas->Clear();
    canvas->SetFillColor(10);

    TPad *pad1 = new TPad("pad1", "", 0.02, 0.51, 0.98, 0.98, 10);
    TPad *pad2 = new TPad("pad2", "", 0.02, 0.02, 0.98, 0.49, 10);
    pad1->SetBit(kCanDelete);
    pad2->SetBit(kCanDelete);
    pad1->Draw();
    pad2->Draw();

    pad1->cd();
    pad1->SetLogx();
    pad1->SetLogy();
    pad1->SetGrid();
    TGraph *gr1 = new TGraph(static_cast<Int_t>(ePlot.size()), ePlot.data(), dedxPlot.data());
    gr1->SetBit(kCanDelete);
    gr1->SetTitle(Form("Bethe-Bloch collision stopping power: %s in %s", partName.c_str(), mat.name.c_str()));
    gr1->SetLineColor(kBlue + 1);
    gr1->SetLineWidth(2);
    gr1->Draw("AL");
    gr1->GetXaxis()->SetLimits(tCut, ePlotMax);
    gr1->GetYaxis()->SetRangeUser(std::max(1.0e-4, minS * 0.75), maxS * 1.25);
    gr1->GetXaxis()->SetTitle("T [MeV]");
    gr1->GetYaxis()->SetTitle("Mass stopping power [MeV cm^{2}/g]");

    pad2->cd();
    pad2->SetGrid();
    TGraph *gr2 = new TGraph(static_cast<Int_t>(depth.size()), depth.data(), stoppingPlot.data());
    gr2->SetBit(kCanDelete);
    gr2->SetTitle(Form("Bragg curve in %s: T_{in}=%.4g MeV; range=%.4g g/cm^{2}", mat.name.c_str(), kinetic, result.rangeMass));
    gr2->SetLineColor(kRed + 1);
    gr2->SetLineWidth(2);
    gr2->Draw("AL");
    gr2->GetXaxis()->SetLimits(0.0, std::max(1.0e-5, result.rangeMass * 1.08));
    gr2->GetYaxis()->SetRangeUser(0.0, maxS * 1.25);
    gr2->GetXaxis()->SetTitle("Depth [g/cm^{2}]");
    gr2->GetYaxis()->SetTitle("Mass stopping power [MeV cm^{2}/g]");
    TLine *rangeLine = new TLine(result.rangeMass, 0.0, result.rangeMass, maxS * 1.2);
    rangeLine->SetBit(kCanDelete);
    rangeLine->SetLineColor(kGray + 2);
    rangeLine->SetLineStyle(2);
    rangeLine->Draw("same");

    delete fSavedGrBethe;
    fSavedGrBethe = static_cast<TGraph*>(gr1->Clone("gr_Bethe_StoppingPower"));
    delete fSavedGrBragg;
    fSavedGrBragg = static_cast<TGraph*>(gr2->Clone("gr_Bethe_BraggCurve"));

    canvas->Modified();
    canvas->Update();

    calc[0]->SetText(Form("%.5g", result.stopping));
    calc[1]->SetText(Form("%.5g", result.massStopping));
    calc[2]->SetText(Form("%.5g", result.rangeCm));
    calc[3]->SetText(Form("%.5g", result.rangeMass));
    calcframe->Layout();

    SaveAllToRoot();
}

// =========================================================================
// METODI SEZIONE MONTE CARLO GEOMETRICO (TELESCOPIO A COINCIDENZA 3F)
// =========================================================================
void UnifiedDetectorGUI::AutoScaleIdealPlane() {
    const double x = fX->GetNumber();
    const double y = fY->GetNumber();
    const double z = fZ->GetNumber();
    double maxXY = std::max({x, y, 60.0}); // Teniamo conto anche degli strati da 60cm
    fLx->SetNumber(std::ceil(maxXY + 6.0 * (z + 26.7 + 2.0 + 20.5 + 2.0))); // Margine in funzione dell'altezza totale
    fLy->SetNumber(std::ceil(maxXY + 6.0 * (z + 26.7 + 2.0 + 20.5 + 2.0)));
    DoCalculateMC();
}

bool UnifiedDetectorGUI::RunMCSimulation(MuonMCResult &out, std::string &error) {
    const double x = fX->GetNumber();
    const double y = fY->GetNumber();
    const double z = fZ->GetNumber(); // Spessore blocco TOP
    const double lx = fLx->GetNumber();
    const double ly = fLy->GetNumber();
    const long long nGen = static_cast<long long>(fNEvents->GetNumber());
    const UInt_t seed = static_cast<UInt_t>(fSeed->GetNumber());
    const double i0 = fI0->GetNumber();
    const double n = fExponent->GetNumber();
    const double eff1 = fEfficiency1->GetNumber() / 100.0; // Efficienza S1 (Top)
    const double eff2 = fEfficiency2->GetNumber() / 100.0; // Efficienza S2 (Mid)
    const double eff3 = fEfficiency3->GetNumber() / 100.0; // Efficienza S3 (Bot)
    const double eff = eff1 * eff2 * eff3;                 // Efficienza totale in coincidenza 3-fold

    if (x <= 0.0 || y <= 0.0 || z <= 0.0 || lx < x || ly < y || nGen < 100) {
        error = "Parametri geometrici non validi (assicurarsi che Lx >= X, Ly >= Y e N_gen >= 100).";
        return false;
    }

    out = MuonMCResult{};
    out.lxIdeal = lx;
    out.lyIdeal = ly;
    out.areaIdeal = lx * ly;
    out.areaTop = x * y;
    out.areaSideX = x * z;
    out.areaSideY = y * z;
    out.areaSidesTot = 2.0 * (out.areaSideX + out.areaSideY);
    out.nEvents = nGen;

    const double omegaEffTop = (2.0 * kPi) / (n + 2.0);
    out.geomIdeal = out.areaIdeal * omegaEffTop;
    out.rateIdealHz = eff * i0 * out.geomIdeal;

    const double omegaEffSide = std::sqrt(kPi) * std::tgamma(0.5 * (n + 1.0)) /
                                (2.0 * std::tgamma(0.5 * (n + 4.0)));
    // La formula analitica vale *solo* per il blocco superiore, la teniamo a scopo di riferimento (non descrive il telescopio intero)
    out.geomAna = out.areaTop * omegaEffTop + out.areaSidesTot * omegaEffSide;

    out.rateAnaTopHz      = eff1 * i0 * out.areaTop * omegaEffTop;
    out.rateAnaSideXHz    = eff1 * i0 * out.areaSideX * omegaEffSide;
    out.rateAnaSideYHz    = eff1 * i0 * out.areaSideY * omegaEffSide;
    out.rateAnaSidesTotHz = 2.0 * (out.rateAnaSideXHz + out.rateAnaSideYHz);
    out.rateAnaTotal3DHz  = eff1 * i0 * out.geomAna;

    // ---------------- GEOMETRIA TELESCOPIO A 3 STRATI ----------------
    const double gap1 = 26.7; // Vuoto sotto il primo blocco
    const double x2 = 20.0, y2 = 60.0, z2 = 2.0; // Scintillatore intermedio
    const double gap2 = 20.5; // Vuoto sotto l'intermedio
    const double x3 = 20.0, y3 = 60.0, z3 = 2.0; // Scintillatore inferiore

    // Definiamo le quote assolute Z dal basso verso l'alto
    const double z3_min = 0.0;
    const double z3_max = z3_min + z3;
    const double z2_min = z3_max + gap2;
    const double z2_max = z2_min + z2;
    const double z1_min = z2_max + gap1;
    const double z1_max = z1_min + z; // z dal pannello
    const double zGen = z1_max; // Quota di generazione

    out.x1=x; out.y1=y; out.z1_min=z1_min; out.z1_max=z1_max;
    out.x2=x2; out.y2=y2; out.z2_min=z2_min; out.z2_max=z2_max;
    out.x3=x3; out.y3=y3; out.z3_min=z3_min; out.z3_max=z3_max;
    out.zGen = zGen;

    // Lambda per verificare l'intersezione di un raggio con un box arbitrario
    auto intersectBox = [&](double x0, double y0, double z_gen, double tx, double ty,
                            double X, double Y, double Zmin, double Zmax) -> bool {
        double hx = 0.5 * X, hy = 0.5 * Y;
        double u_enter = z_gen - Zmax; // Distanza lungo Z all'ingresso
        double u_exit  = z_gen - Zmin; // Distanza lungo Z all'uscita

        double uxmin = -1e30, uxmax = 1e30;
        if (std::abs(tx) < 1e-15) {
            if (x0 < -hx || x0 > hx) return false;
        } else {
            double ux1 = (-hx - x0) / tx;
            double ux2 = ( hx - x0) / tx;
            uxmin = std::min(ux1, ux2);
            uxmax = std::max(ux1, ux2);
        }

        double uymin = -1e30, uymax = 1e30;
        if (std::abs(ty) < 1e-15) {
            if (y0 < -hy || y0 > hy) return false;
        } else {
            double uy1 = (-hy - y0) / ty;
            double uy2 = ( hy - y0) / ty;
            uymin = std::min(uy1, uy2);
            uymax = std::max(uy1, uy2);
        }

        double uStart = std::max({u_enter, uxmin, uymin});
        double uEnd   = std::min({u_exit, uxmax, uymax});

        return uStart <= uEnd;
    };
    // ------------------------------------------------------------------

    TRandom3 rng(seed);
    const double hlx = 0.5 * lx;
    const double hly = 0.5 * ly;
    const double invPower = 1.0 / (n + 2.0);

    double sumTrackLen[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    out.tracks3D.reserve(180);

    for (long long ev = 0; ev < nGen; ++ev) {
        const double x0 = rng.Uniform(-hlx, hlx);
        const double y0 = rng.Uniform(-hly, hly);

        const double uRnd = rng.Uniform(1.0e-12, 1.0);
        const double cosT = std::pow(uRnd, invPower);
        const double sinT = std::sqrt(std::max(0.0, 1.0 - cosT * cosT));
        const double tanT = sinT / cosT;
        const double phi = rng.Uniform(0.0, 2.0 * kPi);

        const double tx = tanT * std::cos(phi);
        const double ty = tanT * std::sin(phi);

        // Verifica coincidenza: deve colpire tutti e tre i volumi (il telescopio)
        bool hit1 = intersectBox(x0, y0, zGen, tx, ty, x, y, z1_min, z1_max);
        bool hit2 = intersectBox(x0, y0, zGen, tx, ty, x2, y2, z2_min, z2_max);
        bool hit3 = intersectBox(x0, y0, zGen, tx, ty, x3, y3, z3_min, z3_max);

        if (hit1 && hit2 && hit3) {
            out.counts[0]++; // Mappato su "Categoria (a)" per convenienza di GUI (Coincidenza 3-Fold)
            double trackLen = (z1_max - z3_min) / cosT; // Lunghezza totale del telescopio (circa)
            sumTrackLen[0] += trackLen;

            if (out.tracks3D.size() < 180) {
                out.tracks3D.push_back({
                    x0, y0, zGen,
                    x0, y0, zGen,
                    x0 + zGen * tx, y0 + zGen * ty, 0.0,
                    0 // Green (Coincidenza)
                });
            }
        } else {
            out.counts[4]++; // Mancati (non colpisce tutte e tre)
            if (out.tracks3D.size() < 45 && out.counts[4] % 100 == 0) {
                out.tracks3D.push_back({x0, y0, zGen, x0, y0, zGen, x0 + zGen * tx, y0 + zGen * ty, 0.0, 4});
            }
        }
    }

    const double nGenD = static_cast<double>(nGen);
    for (int c = 0; c < 5; ++c) {
        const double p = static_cast<double>(out.counts[c]) / nGenD;
        out.acc[c] = p;
        out.accErr[c] = std::sqrt(std::max(0.0, p * (1.0 - p)) / nGenD);
        out.rateMC[c] = out.rateIdealHz * out.acc[c];
        out.rateMCErr[c] = out.rateIdealHz * out.accErr[c];
    }

    out.countRealTot   = out.counts[0];
    out.accRealTot    = static_cast<double>(out.countRealTot) / nGenD;
    out.accRealTotErr = std::sqrt(std::max(0.0, out.accRealTot * (1.0 - out.accRealTot)) / nGenD);

    out.geomMC    = out.accRealTot    * out.geomIdeal;
    out.geomMCErr = out.accRealTotErr * out.geomIdeal;

    out.rateMCTopHz        = out.rateMC[0];
    out.rateMCTopErrHz     = out.rateMCErr[0];
    out.rateMCSideHz       = 0.0;
    out.rateMCSideErrHz    = 0.0;
    out.rateMCRealTotHz    = out.rateMC[0];
    out.rateMCRealTotErrHz = out.rateMCErr[0];

    out.meanTrackLenTot = (out.countRealTot > 0) ? (sumTrackLen[0] / static_cast<double>(out.countRealTot)) : 0.0;
    out.meanTrackLenExp = 0.0; // Inapplicabile al telescopio a 3 strati

    // Scansione Convergenza
    const double baseWidths[12] = {40.0, 60.0, 80.0, 100.0, 130.0, 160.0, 200.0, 250.0, 300.0, 400.0, 600.0, 800.0};
    double maxDim = std::max({x, y, 60.0});
    const double scaleFactor = (maxDim + 6.0 * zGen) / kScanRefSize;
    const long long nScanEvents = std::min(nGen, 1000000LL);
    const double nScanD = static_cast<double>(nScanEvents);

    out.planeScan.clear();
    out.planeScan.reserve(12);

    for (int sIdx = 0; sIdx < 12; ++sIdx) {
        const double lGen = std::ceil(baseWidths[sIdx] * scaleFactor);
        if (lGen < maxDim) continue;
        const double halfL = 0.5 * lGen;
        const double areaS = lGen * lGen;
        long long hitsS = 0;

        for (long long ev = 0; ev < nScanEvents; ++ev) {
            const double x0 = rng.Uniform(-halfL, halfL);
            const double y0 = rng.Uniform(-halfL, halfL);
            const double cosT = std::pow(rng.Uniform(1.0e-12, 1.0), invPower);
            const double tanT = std::sqrt(std::max(0.0, 1.0 - cosT * cosT)) / cosT;
            const double phi = rng.Uniform(0.0, 2.0 * kPi);
            const double tx = tanT * std::cos(phi);
            const double ty = tanT * std::sin(phi);

            bool hit1 = intersectBox(x0, y0, zGen, tx, ty, x, y, z1_min, z1_max);
            bool hit2 = intersectBox(x0, y0, zGen, tx, ty, x2, y2, z2_min, z2_max);
            bool hit3 = intersectBox(x0, y0, zGen, tx, ty, x3, y3, z3_min, z3_max);

            if (hit1 && hit2 && hit3) hitsS++;
        }

        const double pS = static_cast<double>(hitsS) / nScanD;
        const double errS = std::sqrt(std::max(0.0, pS * (1.0 - pS)) / nScanD);
        const double gPlane = areaS * omegaEffTop;
        const double gMC = pS * gPlane;
        const double gErr = errS * gPlane;

        PlaneScanPoint pt;
        pt.lGen = lGen;
        pt.areaGen = areaS;
        pt.nGen = nScanEvents;
        pt.nHits = hitsS;
        pt.accPct = pS * 100.0;
        pt.accPctErr = errS * 100.0;
        pt.geomFactor = gMC;
        pt.geomErr = gErr;
        pt.relDiffPct = 0.0;
        out.planeScan.push_back(pt);
    }

    return true;
}

void UnifiedDetectorGUI::SetMCError(const std::string &message) {
    std::cerr << "\n[ERRORE MONTE CARLO] " << message << "\n" << std::endl;
    for (int i = 0; i < 16; ++i) fRes[i]->SetText("N/A");
    TCanvas *canvas = fCanvasMC1->GetCanvas();
    canvas->Clear();
    TPaveText *box = new TPaveText(0.08, 0.35, 0.92, 0.65, "NDC");
    box->SetBit(kCanDelete);
    box->SetFillColor(0);
    box->SetTextColor(kRed + 1);
    box->SetTextFont(43);
    box->SetTextSize(16);
    box->AddText(message.c_str());
    box->Draw();
    canvas->Modified();
    canvas->Update();
}

void UnifiedDetectorGUI::DrawScene3D(TPad *pad, const MuonMCResult &r) {
    pad->cd();
    pad->Clear();

    PutText(0.03, 0.965, "#color[1]{Vista 3D: Telescopio a Coincidenza e Tracce}", 14.5);
    PutText(0.03, 0.930, Form("A_{tot} = (%.3f #pm %.3f) %%  |  G_{MC} = %.1f #pm %.1f cm^{2}sr",
                              r.accRealTot * 100.0, r.accRealTotErr * 100.0,
                              r.geomMC, r.geomMCErr), 12.0);

    double maxDim = std::max({r.x1, r.y1, r.x2, r.y2, r.x3, r.y3});
    const double R = std::max(maxDim, r.zGen) * 1.1;
    const double zMid = 0.5 * r.zGen;
    TView *view = TView::CreateView(1);
    view->SetRange(-R, -R, zMid - R, R, R, zMid + R);

    auto drawEdge3D = [](double x1, double y1, double z1, double x2, double y2, double z2,
                         Color_t col, Width_t w, Style_t st = 1) {
        TPolyLine3D *pl = new TPolyLine3D(2);
        pl->SetBit(kCanDelete);
        pl->SetPoint(0, x1, y1, z1);
        pl->SetPoint(1, x2, y2, z2);
        pl->SetLineColor(col);
        pl->SetLineWidth(w);
        pl->SetLineStyle(st);
        pl->Draw();
    };

    auto drawBox = [&](double X, double Y, double Zmin, double Zmax, Color_t col) {
        double hx = 0.5 * X, hy = 0.5 * Y;
        double vx[4] = {-hx,  hx,  hx, -hx};
        double vy[4] = {-hy, -hy,  hy,  hy};
        for (int i = 0; i < 4; ++i) {
            int next = (i + 1) % 4;
            drawEdge3D(vx[i], vy[i], Zmax,   vx[next], vy[next], Zmax,   col, 2, 1);
            drawEdge3D(vx[i], vy[i], Zmin, vx[next], vy[next], Zmin, col, 2, 1);
            drawEdge3D(vx[i], vy[i], Zmin, vx[i],    vy[i],    Zmax,   col, 2, 1);
        }
    };

    const double pxIdeal = std::min(0.5 * r.lxIdeal, R * 0.95);
    const double pyIdeal = std::min(0.5 * r.lyIdeal, R * 0.95);
    const bool planeClipped = (0.5 * r.lxIdeal > R * 0.95) || (0.5 * r.lyIdeal > R * 0.95);
    double ix[4] = {-pxIdeal,  pxIdeal,  pxIdeal, -pxIdeal};
    double iy[4] = {-pyIdeal, -pyIdeal,  pyIdeal,  pyIdeal};
    for (int i = 0; i < 4; ++i) {
        int next = (i + 1) % 4;
        drawEdge3D(ix[i], iy[i], r.zGen, ix[next], iy[next], r.zGen, kAzure - 3, 2, 2);
    }

    // Disegna i 3 strati del telescopio
    drawBox(r.x1, r.y1, r.z1_min, r.z1_max, kBlue + 1);
    drawBox(r.x2, r.y2, r.z2_min, r.z2_max, kBlue + 1);
    drawBox(r.x3, r.y3, r.z3_min, r.z3_max, kBlue + 1);

    const Color_t catCols[5] = {kGreen + 2, kBlue + 1, kOrange + 7, kMagenta + 1, kGray + 1};
    for (const auto &tr : r.tracks3D) {
        if (tr.category == 4) {
            if (std::abs(tr.x0) > R || std::abs(tr.y0) > R) continue;
            drawEdge3D(tr.x0, tr.y0, tr.z0, tr.xOut, tr.yOut, tr.zOut, catCols[4], 1, 3);
        } else {
            if (std::abs(tr.xIn) > R || std::abs(tr.yIn) > R) continue;
            drawEdge3D(tr.xIn, tr.yIn, tr.zIn, tr.xOut, tr.yOut, tr.zOut, catCols[tr.category], 2, 1);
        }
    }

    view->ShowAxis();

    if (planeClipped) {
        PutText(0.03, 0.895, Form("#color[1]{Disegno: piano ideale troncato a %.0f#times%.0f cm (reale %.0f#times%.0f cm)}",
                                  2.0 * pxIdeal, 2.0 * pyIdeal, r.lxIdeal, r.lyIdeal), 11.5);
    }
    PutText(0.03, planeClipped ? 0.865 : 0.895, "#color[1]{Tracce illustrative (campione non in proporzione)}", 11.5);
    PutText(0.03, 0.095, "#color[418]{#bf{(a) Tripla Coincidenza}}", 12);
    PutText(0.03, 0.060, "#color[921]{(e) Mancati (miss)}   #color[857]{--- Piano ideale S (z = Z_top)}", 12);
}

void UnifiedDetectorGUI::DrawSummaryTableAndHist(TPad *pad, const MuonMCResult &r) {
    pad->cd();
    pad->Clear();

    TPad *padPlot = new TPad("padPlot", "Confronto Rate", 0.01, 0.59, 0.99, 0.99, 10);
    TPad *padText = new TPad("padText", "Tabella MC",     0.01, 0.01, 0.99, 0.58, 10);
    padPlot->SetBit(kCanDelete);
    padText->SetBit(kCanDelete);
    padPlot->SetBottomMargin(0.18);
    padPlot->SetLeftMargin(0.11);
    padPlot->SetRightMargin(0.06);
    padPlot->SetTopMargin(0.12);
    padPlot->SetGridy(1);
    padPlot->Draw();
    padText->Draw();

    padPlot->cd();
    TH1D *hMC = new TH1D("hMC", "Rate Attese Telescopio a Coincidenza;;Rate [Hz]", 3, 0.0, 3.0);
    hMC->SetBit(kCanDelete);
    hMC->SetDirectory(nullptr);
    hMC->SetStats(false);

    const char *binLabels[3] = {"Coincidenza (a)", "Analitico SINGOLO Top", "Mancati (e)"};
    hMC->GetXaxis()->SetBinLabel(1, binLabels[0]); hMC->SetBinContent(1, r.rateMCRealTotHz); hMC->SetBinError(1, r.rateMCRealTotErrHz);
    hMC->GetXaxis()->SetBinLabel(2, binLabels[1]); hMC->SetBinContent(2, r.rateAnaTotal3DHz); hMC->SetBinError(2, 0.0);
    hMC->GetXaxis()->SetBinLabel(3, binLabels[2]); hMC->SetBinContent(3, r.rateIdealHz - r.rateMCRealTotHz); hMC->SetBinError(3, 0.0);

    hMC->SetFillColor(TColor::GetColor(100, 165, 235));
    hMC->SetLineColor(kBlue + 2);
    hMC->SetBarWidth(0.42);
    hMC->SetBarOffset(0.28);
    hMC->GetYaxis()->SetRangeUser(0.0, std::max(r.rateMCRealTotHz, r.rateAnaTotal3DHz) * 1.18);
    hMC->GetXaxis()->SetLabelSize(0.055);
    hMC->GetYaxis()->SetLabelSize(0.055);
    hMC->GetYaxis()->SetTitleSize(0.060);
    hMC->GetYaxis()->SetTitleOffset(0.82);
    hMC->Draw("BAR");

    delete fSavedHistMCRates;
    fSavedHistMCRates = static_cast<TH1D*>(hMC->Clone("h_MC_Rates_Categorie"));
    fSavedHistMCRates->SetDirectory(nullptr);

    padText->cd();
    TPaveText *info = new TPaveText(0.01, 0.01, 0.99, 0.98, "NDC");
    info->SetBit(kCanDelete);
    info->SetFillColor(TColor::GetColor(250, 250, 245));
    info->SetBorderSize(1);
    info->SetTextAlign(12);
    info->SetTextFont(43);
    info->SetTextSize(11.5);
    info->SetMargin(0.025);

    info->AddText(Form("#color[4]{1) Piano Ideale S (z = %.1f cm):}", r.zGen));
    info->AddText(Form("   S = %.0f#times%.0f = %.0f cm^{2} | #nu_{piano} = #bf{%.2f Hz} | #bf{G_{MC} (telescopio) = %.2f #pm %.2f cm^{2}sr}",
                       r.lxIdeal, r.lyIdeal, r.areaIdeal, r.rateIdealHz, r.geomMC, r.geomMCErr));
    info->AddText(Form("#color[4]{2) Analisi Muoni su %lld Eventi Generati:}", r.nEvents));
    info->AddText(Form("   #bullet #bf{(a) Coincidenza a 3-Fold:} N_{a} = %lld | A_{a} = (%.4f #pm %.4f)%% | #nu_{a} = #bf{%.3f #pm %.3f Hz} (%.1f /min)",
                       r.counts[0], r.acc[0] * 100.0, r.accErr[0] * 100.0, r.rateMC[0], r.rateMCErr[0], r.rateMC[0] * 60.0));
    info->AddText(Form("   #bullet #bf{(e) Mancati (Miss):}      N_{e} = %lld | A_{e} = (%.4f #pm %.4f)%%",
                       r.counts[4], r.acc[4] * 100.0, r.accErr[4] * 100.0));
    info->AddText("#color[4]{3) Riferimento Analitico:}");
    info->AddText(Form("   La formula G_{Ana} (= %.2f cm^{2} sr) è esatta SOLO per il singolo blocco Top.", r.geomAna));
    info->AddText(Form("   Rate atteso per SINGOLO blocco Top isolato = %.2f Hz.", r.rateAnaTotal3DHz));
    info->Draw();
}

void UnifiedDetectorGUI::DrawMCTab1(const MuonMCResult &r) {
    TCanvas *canvas = fCanvasMC1->GetCanvas();
    canvas->Clear();
    canvas->SetFillColor(10);

    const double cw = (canvas->GetWw() > 0) ? canvas->GetWw() : 1010.0;
    const double ch = (canvas->GetWh() > 0) ? canvas->GetWh() : 760.0;
    const double x1 = 0.01, x2 = 0.45;
    const double sidePx = (x2 - x1) * cw;
    const double dy = std::min(0.96, sidePx / ch);
    const double y1 = 0.5 * (1.0 - dy);
    const double y2 = 0.5 * (1.0 + dy);

    TPad *pad3D = new TPad("pad3D", "Vista 3D",          x1,   y1,   x2,   y2,   10);
    TPad *pad2D = new TPad("pad2D", "Grafici e Tabella", 0.46, 0.01, 0.99, 0.99, 10);
    pad3D->SetBit(kCanDelete);
    pad2D->SetBit(kCanDelete);
    pad3D->Draw();
    pad2D->Draw();

    DrawScene3D(pad3D, r);
    DrawSummaryTableAndHist(pad2D, r);

    canvas->Modified();
    canvas->Update();
}

void UnifiedDetectorGUI::DrawMCTab2(const MuonMCResult &r) {
    TCanvas *canvas = fCanvasMC2->GetCanvas();
    canvas->Clear();
    canvas->SetFillColor(10);

    TPad *pLeft  = new TPad("pTab2Left",  "Accettanza % vs L_gen", 0.01, 0.01, 0.495, 0.99, 10);
    TPad *pRight = new TPad("pTab2Right", "G [cm^2 sr] vs L_gen",  0.505, 0.01, 0.99, 0.99, 10);
    pLeft->SetBit(kCanDelete);
    pRight->SetBit(kCanDelete);
    pLeft->SetLeftMargin(0.13);  pLeft->SetRightMargin(0.04);  pLeft->SetBottomMargin(0.12); pLeft->SetGrid();
    pRight->SetLeftMargin(0.14); pRight->SetRightMargin(0.04); pRight->SetBottomMargin(0.12); pRight->SetGrid();
    pLeft->Draw();
    pRight->Draw();

    const int nPts = static_cast<int>(r.planeScan.size());
    if (nPts == 0) return;

    TGraphErrors *grAcc = new TGraphErrors(nPts);
    grAcc->SetBit(kCanDelete);
    TGraphErrors *grGeom = new TGraphErrors(nPts);
    grGeom->SetBit(kCanDelete);

    double minG = 1e9, maxG = -1e9;
    for (int i = 0; i < nPts; ++i) {
        const auto &pt = r.planeScan[i];
        grAcc->SetPoint(i, pt.lGen, pt.accPct);
        grAcc->SetPointError(i, 0.0, pt.accPctErr);
        grGeom->SetPoint(i, pt.lGen, pt.geomFactor);
        grGeom->SetPointError(i, 0.0, pt.geomErr);

        if (pt.geomFactor - pt.geomErr < minG) minG = pt.geomFactor - pt.geomErr;
        if (pt.geomFactor + pt.geomErr > maxG) maxG = pt.geomFactor + pt.geomErr;
    }

    pLeft->cd();
    pLeft->SetLogy(1);
    grAcc->SetTitle("Accettanza Adimensionale A(S) vs Piano Generato;Lato piano ideale L_{gen} [cm]  (S = L_{gen}^{2});Accettanza geometrica A(S) [%]");
    grAcc->SetMarkerStyle(20);
    grAcc->SetMarkerSize(1.0);
    grAcc->SetMarkerColor(kBlue + 1);
    grAcc->SetLineColor(kBlue + 1);
    grAcc->SetLineWidth(2);
    grAcc->GetXaxis()->SetTitleSize(0.042);
    grAcc->GetYaxis()->SetTitleSize(0.042);
    grAcc->GetYaxis()->SetTitleOffset(1.30);
    grAcc->Draw("AP");

    PutText(0.32, 0.83, "#color[4]{Scala come A(S) #approx cost / S  #propto 1/L_{gen}^{2}}", 12.5);

    pRight->cd();
    grGeom->SetTitle("Convergenza Fattore Geometrico G(S) vs Piano Generato;Lato piano ideale L_{gen} [cm]  (S = L_{gen}^{2});Fattore geometrico G [cm^{2} sr]");
    grGeom->SetMarkerStyle(21);
    grGeom->SetMarkerSize(1.0);
    grGeom->SetMarkerColor(kAzure + 2);
    grGeom->SetLineColor(kAzure + 2);
    grGeom->SetLineWidth(2);
    const double spanG = std::max(1.0, maxG - minG);
    grGeom->GetYaxis()->SetRangeUser(std::max(0.0, minG - 0.2 * spanG), maxG + 0.2 * spanG);
    grGeom->GetXaxis()->SetTitleSize(0.042);
    grGeom->GetYaxis()->SetTitleSize(0.042);
    grGeom->GetYaxis()->SetTitleOffset(1.45);
    grGeom->Draw("AP");

    delete fSavedGrScanAcc;
    fSavedGrScanAcc = static_cast<TGraphErrors*>(grAcc->Clone("gr_MC_Accettanza_vs_Lgen"));
    delete fSavedGrScanGeom;
    fSavedGrScanGeom = static_cast<TGraphErrors*>(grGeom->Clone("gr_MC_FattoreGeom_vs_Lgen"));

    canvas->Modified();
    canvas->Update();
}

void UnifiedDetectorGUI::DoCalculateMC() {
    MuonMCResult r;
    std::string error;
    if (!RunMCSimulation(r, error)) {
        SetMCError(error);
        return;
    }

    fRes[0]->SetText(Form("%.0f (%.1f Hz)", r.areaIdeal, r.rateIdealHz));
    fRes[1]->SetText(Form("%.4f +/- %.4f %%", r.accRealTot * 100.0, r.accRealTotErr * 100.0));
    fRes[2]->SetText(Form("%.4f +/- %.4f %%", r.acc[0] * 100.0, r.accErr[0] * 100.0));
    fRes[3]->SetText(" - N/A - ");
    fRes[4]->SetText(" - N/A - ");
    fRes[5]->SetText(" - N/A - ");
    fRes[6]->SetText(Form("%.4f +/- %.4f %%", r.acc[4] * 100.0, r.accErr[4] * 100.0));
    fRes[7]->SetText(Form("%.2f +/- %.2f", r.geomMC, r.geomMCErr));
    fRes[8]->SetText(Form("%.2f", r.geomAna));
    fRes[9]->SetText(Form("%.3f +/- %.3f", r.rateMC[0], r.rateMCErr[0]));
    fRes[10]->SetText(" - N/A - ");
    fRes[11]->SetText(" - N/A - ");
    fRes[12]->SetText(" - N/A - ");
    fRes[13]->SetText(" - N/A - ");
    fRes[14]->SetText(" - N/A - ");
    fRes[15]->SetText(Form("%.3f +/- %.3f", r.rateMCRealTotHz, r.rateMCRealTotErrHz));
    fMain->Layout();

    if (fRightTabs->GetCurrent() == 0) fRightTabs->SetTab(1);
    DrawMCTab1(r);
    DrawMCTab2(r);

    const double e1 = fEfficiency1->GetNumber();
    const double e2 = fEfficiency2->GetNumber();
    const double e3 = fEfficiency3->GetNumber();
    const double eTot = (e1 * e2 * e3) / 10000.0;

    auto printReport = [&](std::ostream &os) {
        os << std::fixed
           << "\n========================================================================================\n"
           << " SIMULAZIONE MONTE CARLO TELESCOPIO (N_gen = " << r.nEvents << ")\n"
           << "========================================================================================\n"
           << "  Geometria (Z min .. Z max):\n"
           << "    Top : Z=" << r.z1_min << ".." << r.z1_max << " (" << r.x1 << "x" << r.y1 << " cm)\n"
           << "    Mid : Z=" << r.z2_min << ".." << r.z2_max << " (" << r.x2 << "x" << r.y2 << " cm)\n"
           << "    Bot : Z=" << r.z3_min << ".." << r.z3_max << " (" << r.x3 << "x" << r.y3 << " cm)\n"
           << "  Piano Ideale S (Lx x Ly)       : " << setprecision(1) << r.lxIdeal << " x " << r.lyIdeal
           << " cm  (S = " << setprecision(0) << r.areaIdeal << " cm^2, nu_piano = " << setprecision(2) << r.rateIdealHz << " Hz)\n"
           << "  Flusso & Efficienze impostate  : I0 = " << setprecision(5) << fI0->GetNumber()
           << " cm^-2 s^-1 sr^-1 | n = " << setprecision(2) << fExponent->GetNumber() << "\n"
           << "    eff_S1 (Top) = " << setprecision(2) << e1 << " % | eff_S2 (Mid) = " << e2
           << " % | eff_S3 (Bot) = " << e3 << " %  =>  eff_tot (3F) = " << eTot << " %\n"
           << "----------------------------------------------------------------------------------------\n"
           << "  Gruppo Topologico   |  Conteggi   |  Accettanza A [%]     |  Rate MC [Hz]     \n"
           << "  --------------------------------------------------------------------------------------\n"
           << "  (a) Coincidenza     | " << setw(11) << r.counts[0]
           << " | " << setw(8) << setprecision(4) << r.acc[0] * 100.0 << " +/- " << setw(6) << r.accErr[0] * 100.0
           << " | " << setw(7) << setprecision(3) << r.rateMC[0] << " +/- " << setw(6) << r.rateMCErr[0] << "\n"
           << "  (e) Mancati (miss)  | " << setw(11) << r.counts[4]
           << " | " << setw(8) << setprecision(4) << r.acc[4] * 100.0 << " +/- " << setw(6) << r.accErr[4] * 100.0
           << " |          --       \n"
           << "----------------------------------------------------------------------------------------\n"
           << "  RIEPILOGO TELESCOPIO:\n"
           << "    Fattore Geometrico G_MC      : " << setprecision(2) << r.geomMC << " +/- " << r.geomMCErr << " cm^2 sr\n"
           << "    RATE COINCIDENZA 3F          : " << setprecision(3) << r.rateMCRealTotHz << " +/- " << r.rateMCRealTotErrHz
           << " Hz   (" << setprecision(1) << r.rateMCRealTotHz * 60.0 << " /min)\n"
           << "----------------------------------------------------------------------------------------\n"
           << "  SCANSIONE CONVERGENZA IN FUNZIONE DEL PIANO GENERATO S = L_gen^2:\n"
           << "    L_gen [cm] |   S [cm^2]  |   Accettanza A [%]   |  G_MC [cm^2 sr]    \n"
           << "  --------------------------------------------------------------------------------------\n";
        for (const auto &pt : r.planeScan) {
            os << "    " << std::setw(10) << std::setprecision(0) << pt.lGen
               << " | " << std::setw(11) << std::setprecision(0) << pt.areaGen
               << " | " << std::setw(8) << std::setprecision(4) << pt.accPct << " +/- " << std::setw(6) << pt.accPctErr
               << " | " << std::setw(7) << std::setprecision(2) << pt.geomFactor << " +/- " << std::setw(5) << pt.geomErr << "\n";
        }
        os << "========================================================================================\n" << std::endl;
    };

    printReport(std::cout);
    std::ofstream reportFile("mc_results_telescopio.txt");
    if (reportFile.is_open()) printReport(reportFile);

    SaveAllToRoot();
}

void lab2() {
    readfile();
    new UnifiedDetectorGUI(gClient->GetRoot(), 1480, 820);
}

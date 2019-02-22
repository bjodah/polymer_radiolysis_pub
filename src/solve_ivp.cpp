#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <optional>
#include <unordered_map>
#include <vector>
#include <kinetgine/crn.hpp>
#include <kinetgine/json_util.hpp>

#include "clara.hpp"
#include "json.hpp"

namespace ke = kinetgine;
namespace ca = clara;
using json = nlohmann::json;

namespace {
    template<typename T, typename U>
    bool occurs(const U item, const T container){
        return std::find(std::begin(container), std::end(container), item) != std::end(container);
    }

    auto per_key(const std::vector<std::string> &keys, json j, std::optional<double> default_value) {
        std::vector<double> result;
        result.reserve(keys.size());
        for (json::iterator it = j.begin(); it != j.end(); ++it){
            if (!occurs(it.key(), keys))
                throw std::runtime_error(ke::StreamFmt() << "Unexpected key encountered: " << it.key());
        }
        for (const auto &k : keys ){
            auto it = j.find(k);
            const bool found = it != j.end();
            if (!found and !default_value)
                throw std::runtime_error(ke::StreamFmt() << "missing key in JSON object: " << k);
            result.emplace_back(found ? it->get<double>() : default_value.value());
        }
        return result;
    }
}

auto solve_ivp(
    ke::OdeSysCRN &odesys,
    double duration,
    std::string pars_s,
    std::string ic_s,
    std::string atol_s,
    std::string settings_s,
    std::string varied_path,
    int npoints)
{
    ke::SolverSettings solver_settings;
    std::vector<double> ic;
    std::vector<std::string> missing_parameter_keys;
    {   // Parameters
        std::vector<double> pars;
        json pars_j;
        try{
            pars_j = json::parse(pars_s);
        } catch (...) {
            std::cerr << __FILE__ << ":" << __LINE__ << "error: --parameters (-p) malformated" << std::endl;
            throw;
        }
        for (json::iterator it = pars_j.begin(); it != pars_j.end(); ++it){
            if (!occurs(it.key(), odesys.param_names))
                throw std::runtime_error(ke::StreamFmt() << "Unknown key encountered in parameters: " << it.key());
        }
        for (const auto &k : odesys.param_names ) {
            if (pars_j.count(k) == 0){
                missing_parameter_keys.push_back(k);
                pars.emplace_back(0.0); // filler value
            }else{
                pars.emplace_back(pars_j[k].get<double>());
            }
        }
        odesys.set_param_values(pars.data());
    }
    {   // Initial conditions
        json ic_j;
        try{
            ic_j = json::parse(ic_s);
        } catch (...) {
            std::cerr << __FILE__ << ":" << __LINE__ << "error: --initial-conditions (-c) malformated" << std::endl;
            throw;
        }
        assert(ic_j.size() == 2);
        ic = per_key(odesys.names, ic_j[0], ic_j[1].get<double>());
    }
    {   // Absolute tolerances
        if (atol_s.front() != '[' or atol_s.back() != ']') { // not json formated, assume it's a path
            std::ifstream ifs(atol_s);
            std::stringstream buffer;
            buffer << ifs.rdbuf();
            atol_s = buffer.str();
        }
        json atol_j;
        try{
            atol_j = json::parse(atol_s);
        } catch (...) {
            std::cerr << __FILE__ << ":" << __LINE__ << "error: abstol malformated" << std::endl;
            throw;
        }
        assert(atol_j.size() == 2);
        solver_settings.atol = per_key(odesys.names, atol_j[0], atol_j[1].get<double>());

    }
    { // Solver settings (other than atol)
        json j;
        try {
            j = json::parse(settings_s);
        } catch (...) {
            std::cerr << __FILE__ << ":" << __LINE__ << "error: --settings (-s) malformated" << std::endl;
            throw;
        }
#define STRINGIFY_(s) #s
#define STRINGIFY(s) STRINGIFY_(s)
#define SETTING(KEY, TYPE) \
        if (j.count(STRINGIFY(KEY))) {                           \
            solver_settings.KEY = j[STRINGIFY(KEY)].get<TYPE>(); \
            j.erase(STRINGIFY(KEY));                             \
        }
        SETTING(rtol, double);
        SETTING(mxsteps, int);
        SETTING(dx0, double);
        SETTING(dx_min, double);
        SETTING(dx_max, double);
        SETTING(method, std::string);
        SETTING(iter_type, std::string);
        SETTING(linear_solver, int);
        SETTING(maxl, int);
        SETTING(eps_lin, double);
        SETTING(nderiv,unsigned);
        SETTING(return_on_root, bool);
        SETTING(autorestart, int);
        SETTING(return_on_error, bool);
        SETTING(with_jacobian, bool);
        SETTING(with_jtimes, bool);
#undef SETTING
#undef STRINGIFY
#undef STRINGIFY_
        if (j.size())
            throw std::runtime_error("Unrecognized keys in settings");
    }
    if (varied_path.size()){
        if (duration)
            throw std::runtime_error("When varied_path is given, set duration=0");
        std::vector<double> durations;
        std::ifstream ifs(varied_path);
        json j;
        ifs >> j;
        if (j.size() != 2)
            throw std::runtime_error("Expected pair of durations/key-value");
        auto ndur = j[0].size();
        std::vector<std::string> varied_keys;
        for (json::iterator it=j[1].begin(); it != j[1].end(); ++it){
            if (it.value().size() != ndur) {
                throw std::runtime_error("Length mismath");
            }
            const auto k = it.key();
            if (occurs(k, missing_parameter_keys)){
                missing_parameter_keys.erase(missing_parameter_keys.begin() +
                                             ke::index_of_check(k, missing_parameter_keys));
            } else {
                throw std::runtime_error(ke::StreamFmt() << "Parameter key already given: " << k);
            }
            varied_keys.push_back(it.key());
        }
        if (missing_parameter_keys.size()){
            throw std::runtime_error(ke::StreamFmt() << "Missing parameter key: "
                                     << missing_parameter_keys.front());
        }
        std::vector<double> varied_values;
        for (unsigned i=0; i < ndur; ++i){
            durations.push_back(j[0][i].get<double>());
            for (json::iterator it=j[1].begin(); it != j[1].end(); ++it){
                varied_values.push_back(it.value()[i].get<double>());
            }
        }
        std::vector<int> varied_indices;
        for (const auto &k : varied_keys){
            varied_indices.push_back(ke::index_of_check(k, varied_keys));
        }
        return ke::chained_parameter_variation_predefined(
            &odesys, durations, ic, varied_values, varied_indices, npoints, solver_settings);
    } else {
        if (duration <= 0){
            throw std::runtime_error("A positive non-zero duration needed");
        }
        if (missing_parameter_keys.size()){
            throw std::runtime_error(ke::StreamFmt() << "Missing parameter key: "
                                     << missing_parameter_keys.front());
        }
        return ke::integrate_adaptive(&odesys, 0, duration, ic, solver_settings);
    }
}

int main(int argc, char **argv){
    double duration = 3600;
    int npoints = 2;
    bool verbose, clip;
    std::string serialized_data;
    std::string parameters = "{}";
    std::string ic = "";  // [{}, 0.0]
    std::string ic_path = "";
    std::string abstol = "[{}, 1e-8]";
    std::string solver_settings = "{\"rtol\": 1e-10}";
    std::string varied_path = "";
    std::string output_path = "results.txt";
    std::unique_ptr<ke::OdeSysCRN> odesys;
    {
        bool show_help = false;
        std::string input_path = "odesys_crn.dat";
        auto cli =
            ca::Help(show_help)
            | ca::Opt(input_path, input_path)["-i"]["--input"]("Path to deserialize from")
            | ca::Opt(output_path, output_path)["-o"]["--output"]("Path to save results to")
            | ca::Opt(npoints, std::to_string(npoints))["--npoints"]("Number of points per duration")
            | ca::Opt(duration, std::to_string(duration))["-d"]["--duration"]("Duration (in seconds)")
            | ca::Opt(parameters, parameters)["-p"]["--parameters"]("Values of parameters (JSON)")
            | ca::Opt(ic, ic)["-c"]["--initial-conditions"]("Initial conditions (JSON)")
            | ca::Opt(ic_path, ic_path)["--initial-conditions-path"]("Initial conditions (JSON) from file")
            | ca::Opt(abstol, abstol)["-a"]["--abstol"]("Per component absolute tolerances (JSON)")
            | ca::Opt(solver_settings, solver_settings)["-s"]["--settings"]("Solver settings (JSON)")
            | ca::Opt(varied_path, varied_path)["-v"]["--varied"]("Path to file for varied data (JSON)")
            | ca::Opt(clip)["--clip"]("Clip away negative concentrations during integration")
            | ca::Opt(verbose)["--verbose"]("Prints reactions in CRN to stdout")
            ;
        auto parsed = cli.parse( ca::Args(argc, argv) );
        if (!parsed) {
            std::cerr << "Error in command line: " << parsed.errorMessage() << std::endl;
            return 1;
        }
        if (show_help) {
            std::cout << cli;
            return 0;
        }
        std::ifstream ifs(input_path);
        cereal::BinaryInputArchive iar(ifs);
        iar(odesys);
    }
    odesys->clip_away_negative = clip;
    if (ic == "" && ic_path.size() > 0) {
        std::ifstream ifs(ic_path);
        std::stringstream buffer;
        buffer << "[" << ifs.rdbuf() << ", 0.0]";
        ic = buffer.str();
    }
    auto result = solve_ivp(*odesys, duration, parameters, ic, abstol, solver_settings, varied_path, npoints);
    auto ofh = std::ofstream(output_path);
    ofh << "#t ";
    for (auto it = odesys->names.begin();;){
        ofh << *it;
        ++it;
        if (it == odesys->names.end()){
            ofh << '\n';
            break;
        } else {
            ofh << ' ';
        }
    }
    ofh << "#";
    result->info.dump_ascii(ofh, "=", ", ");
    ofh << ", npoints=" << npoints;
    ofh << "\n";
    ofh << std::scientific << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    result->dump_ascii(ofh);
    if (verbose){
        result->info.dump_ascii(std::cout, ": ", "\n");
    }
    return !(result->info.nfo_int["success"]);
}

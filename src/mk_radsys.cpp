#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <unordered_map>
#include <vector>

#include "json.hpp"
#include <kinetgine/parsers.hpp>
#include "ensemble.hpp"
#include "clara.hpp"


using json = nlohmann::json;
namespace en = ensemble;
namespace ke = kinetgine;
namespace ca = clara;

#define N_PRINCIPALS 5
const char tokens[N_PRINCIPALS] = {'n', 'a', 'p', 'o', 'l'};

namespace {
    template<typename T>
    std::unordered_map<std::string, T> to_umap_sx(json& j){
        std::unordered_map<std::string, T> result;
        for (json::iterator it = j.begin(); it != j.end(); ++it){
            result[it.key()] = it.value().get<T>();
        }
        return result;
    }
}

void mk_odesys_crn(
    std::string output_path,
    std::string substances_path,
    std::string reactions_path,
    std::string primary_data_s,
    double n0,
    std::array<std::string, N_PRINCIPALS> napol_s,
    bool identify_quads,
    bool verbose)
{
    if (n0 < 1){
        throw std::runtime_error("n0 must be atleat 1");
    }
    std::array<std::vector<int>, N_PRINCIPALS> napol;
    for (unsigned i=0; i < N_PRINCIPALS; ++i) {
        json j;
        try {
            j = json::parse(napol_s[i]);
        } catch(...){
            std::cerr << "princial '" << tokens[i] << "' malformated" << std::endl;
            throw;
        }
        if (not j.is_array()){
            throw std::runtime_error("need an json array");
        }
        for (auto elem : j){
            napol[i].push_back(elem.get<double>());
        }
    }
    auto primary_data = en::PrimaryData();
    {
        json j;
        try {
            j = json::parse(primary_data_s);
        } catch(...) {
            std::cerr << __FILE__ << ':' << __LINE__ << ":error: --primary-data malformated" << std::endl;
            throw;
        }
#define STRINGIFY_(s) #s
#define STRINGIFY(s) STRINGIFY_(s)
#define SETTING(KEY, TYPE)                                       \
        if (j.count(STRINGIFY(KEY))) {                           \
            primary_data.KEY = j[STRINGIFY(KEY)].get<TYPE>();    \
            j.erase(STRINGIFY(KEY));                             \
        }
        SETTING(k_oxyl_H_intra_s, double);
        SETTING(k_oxyl_frag_s, double);
        SETTING(k_alkyl_O2_M_s, double);
        SETTING(k_OH_1kDa_M_s, double);
        SETTING(ratio_OH_H, double);
        SETTING(MW_kDa, double);
        SETTING(MW_monomer_g_mol, double);
        SETTING(t_half_loop_closure_intra_s_94kDa, double);
        SETTING(alkyl_inter_M_s, double);
        SETTING(peroxy_inter_M_s, double);
        SETTING(N_scal_intra, double);
        SETTING(w_scal_OH, double);
        SETTING(radrad_scal_intra, double);
        SETTING(disprop_frac, double);
        SETTING(loop_scaling_inter, double);
        SETTING(loop_scaling_intra, double);
        SETTING(max_loop_prob, bool);
#undef SETTING
#undef STRINGIFY
#undef STRINGIFY_
        if (j.size())
            throw std::runtime_error(ke::StreamFmt() << "Unrecognized keys in primary_data, first: " << j[0].get<std::string>());
    }
    std::cout << "For reference: a 1-wt% solution would have [" << tokens[0] << n0;
    for (unsigned i=1; i < N_PRINCIPALS; ++i){
        std::cout << tokens[i] << 0;
    }
    std::cout << "] / M = " << 10.0/(primary_data.MW_kDa*1e3) << std::endl;
    std::cout << "Creating colletion..." << std::flush;
    auto coll = en::mk_collection(napol);
    std::cout << " done" << std::endl;

    std::cout << "Formulation:\n";
    auto form = en::create_formulation(primary_data, n0);
    std::cout << "  Monomers per chain: " << form.monomers_per_chain << std::endl;
    std::cout << "  Loop scaling inter: " << form.loop_scaling_inter << std::endl;
    std::cout << "  Loop scaling intra: " << form.loop_scaling_intra << std::endl;
    std::cout << "  max_loop_prob:      " << form.max_loop_prob      << std::endl;

    std::cout << "Generating reactions..." << std::flush;
    auto rxns = form.generate_reactions(coll, true);
    std::cout << " done" << std::endl;

    std::cout << "Loading other reactions & substances..." << std::flush;
    ke::CRN::vec_s substances;
    {
        std::ifstream ifs(substances_path);
        std::stringstream buffer;
        buffer << ifs.rdbuf();
        substances = ke::parse_json_substances(buffer.str());
    }
    ke::CRN::vec_rN reactions;
    {
        std::ifstream ifs(reactions_path);
        std::stringstream buffer;
        buffer << ifs.rdbuf();
        auto rN_rA = ke::parse_json_reactions(buffer.str());
        reactions = rN_rA.first;
    }
    std::cout << " done" << std::endl;
    std::cout << "Generating CRN instance..." << std::flush;
    auto crn = coll.as_crn(rxns, substances, reactions, {}, identify_quads);
    std::cout << " done ("<< crn.substance_keys.size() << " substances, of which "<< crn.nquads << " are pure quadratures"<<")\n";
    if (verbose) {
        std::cout << "Number of reactions: " << crn.rN.size() + crn.rA.size() << std::endl;
        std::cout << "Number of species: " << crn.substance_keys.size() << std::endl;
        std::cout << "   of which are pure quadratures: " << crn.nquads << std::endl;
        std::cout << "\nChemical reaction network:\n";
        for (const auto& r : crn.rN){
            std::cout << "    " << r->reaction_string() << '\n';
        }
        for (const auto& r : crn.rA){
            std::cout << "    " << r->reaction_string() << '\n';
        }
        std::cout << '\n';
    }
    auto ode_settings = ke::OdeSettings{-1, -1, -1, true, 2};
    std::cout << "Creating OdeSysCRN instance..." << std::flush;
    auto odesys = ke::mk_odesys_crn(crn, ode_settings);
    const std::unordered_map<std::string, double> params = {{"doserate", 0.0}};
    std::cout << " done\n";
    odesys->set_param_values(params);
    {
        std::cout << "Saving to disk..." << std::flush;
        auto ofh = std::ofstream(output_path);
        ofh << ke::dumps(odesys);
        std::cout << " done\n";
    }
}


int main(int argc, char **argv){
    bool show_help = false;
    std::string output_path = "odesys_crn.dat";
    std::string substances_path = "substances.json";
    std::string reactions_path = "reactions.json";
    std::string primary_data = "{}";
    double n0 = 0.0;
    bool verbose = false, identify_quads = false;
    std::array<std::string, N_PRINCIPALS> napol {
        "[1,2,4]", "[0,1,2]", "[0,1,2]", "[0,1,2]", "[0,1,2]"};
    auto cli =
        ca::Help(show_help)
        | ca::Opt(output_path, output_path)["--output"]("Path to serialize to")
        | ca::Opt(reactions_path, reactions_path)["--reactions"]("Path to reactions json file (input)")
        | ca::Opt(substances_path, substances_path)["--substances"]("Path to substance json file (input)")
        | ca::Opt(n0, std::to_string(n0))["--n0"]("Number of segments from start")
        | ca::Opt(napol[0], napol[0])["-n"]["--levels-n"]("Levels for 'n'")
        | ca::Opt(napol[1], napol[1])["-a"]["--levels-a"]("Levels for 'a'")
        | ca::Opt(napol[2], napol[2])["-p"]["--levels-p"]("Levels for 'p'")
        | ca::Opt(napol[3], napol[3])["-o"]["--levels-o"]("Levels for 'o'")
        | ca::Opt(napol[4], napol[4])["-l"]["--levels-l"]("Levels for 'l'")
        | ca::Opt(primary_data, primary_data)["--primary-data"]("Primary data")
        | ca::Opt(identify_quads)["--identify-quads"]("Will treat terminal species as pure quadratures")
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
    mk_odesys_crn(output_path, substances_path, reactions_path, primary_data, n0, napol, identify_quads, verbose);
}

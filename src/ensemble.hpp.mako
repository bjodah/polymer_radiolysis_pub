// -*- coding: utf-8; compile-command: "CXX=g++-7 CXXFLAGS=-D_GLIBCXX_DEBUG make DEBUG=1" -*-
/*
This file implements a model for the aqueous radiolysis of polymers.
Copyright (C) 2019  Björn Ingvar Dahlgren

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#if defined(ENSEMBLE_VERBOSE)
#  include <iostream>
#  include <iomanip>
#endif
#include <optional>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <limits>

#include <kinetgine/crn.hpp>

#if !defined(NAMESPACE_BEGIN)
#define NAMESPACE_BEGIN(s) namespace s{
#define NAMESPACE_END(s) }
#endif

NAMESPACE_BEGIN(ensemble)

using vec_i = std::vector<int>;
using vec_s = std::vector<std::string>;
using vec_c = std::vector<char>;
using umap_ii = std::unordered_map<int, int>;
using umap_id = std::unordered_map<int, double>;
using umap_sd = std::unordered_map<std::string, double>;
using umap_si = std::unordered_map<std::string, int>;
using umap_c_vi = std::unordered_map<char, std::vector<int> >;
using umap_i_vi = std::unordered_map<int, std::vector<int> >;
using umap_c_uivi = std::unordered_map<char, umap_i_vi>;
using umap_s_uii = std::unordered_map<std::string, umap_ii>;
using vec_vec_i = std::vector<std::vector<int> >;
using vec_uii = std::vector<umap_ii>;

using o_umap_id = std::optional<umap_id>;

template<typename T, typename U>
int index_of_check(T item, U container) {
    auto it = std::find(container.begin(), container.end(), item);
    auto idx = std::distance(container.begin(), it);
    if (idx >= static_cast<int>(container.size()))
        throw std::runtime_error("Item not in container");
    return idx;
}

template<typename K, typename V, typename F>
void expel(std::unordered_map<K, V>umap, F pred)
{
    for (auto it=std::begin(umap); it != std::end(umap);)
    {
        if (pred(it->second)) {
            it = umap.erase(it);
        } else {
            ++it;
        }
    }
}

class StreamFmt
{
    std::stringstream m_s;
public:
    StreamFmt() {}
    ~StreamFmt() {}

    template <typename T>
    StreamFmt& operator << (const T& v) {
        this->m_s << v;
        return *this;
    }

    std::string str() const {
        return this->m_s.str();
    }
    operator std::string() const {
        return this->m_s.str();
    }

};

struct Collection {
    using sidx_t = unsigned short;
    using pidx_t = unsigned char;
    using lvl_val_t = int;
    using lvl_idx_t = unsigned char;
    using Reac = std::tuple<umap_ii, umap_id, double>;
    using ReacVec = std::vector<Reac>;

    sidx_t size;
    pidx_t ndim;
    vec_c principal_keys;
    vec_s other_keys;
    std::vector<int> strides;
    std::vector<std::vector<lvl_val_t> > levels_vec;
    double check_abstol = 1e-8;

    Collection(const vec_c &principal_keys, const umap_c_vi &levels,
               const umap_c_uivi &principal_composition,
               const vec_s &other_keys, const umap_s_uii &other_compositions) :
        ndim(principal_keys.size()), principal_keys(principal_keys), other_keys(other_keys)
    {
        if (principal_keys.size() > std::numeric_limits<pidx_t>::max())
            throw std::runtime_error("Too many principal keys");
        int i = 1;
        for (const auto &pk : principal_keys) {
            const auto &lvls = levels.at(pk);
            const auto sz = lvls.size();
            dim.push_back(sz);
            if (sz >= std::numeric_limits<lvl_idx_t>::max()){
                throw std::runtime_error("lvl_idx_t too restricted");
            }
            for (const auto &lvl : lvls){
                if (lvl < 0){
                    throw std::runtime_error("Negative level values do not make sense");
                }
                if (lvl > std::numeric_limits<lvl_val_t>::max() or lvl < std::numeric_limits<lvl_val_t>::min()){
                    throw std::runtime_error("lvl_val_t too restricted");
                }
            }
            {
                std::vector<lvl_val_t> lvls_;
                lvls_.reserve(sz);
                std::transform(lvls.begin(), lvls.end(), std::back_inserter(lvls_),
                               [&](int i){ return (lvl_val_t)i; });
                levels_vec.emplace_back(lvls_);
            }
            principal_composition_polys.push_back(principal_composition.at(pk));
            strides.push_back(i);
            i *= sz;
        }
        strides.push_back(i);
        int size_ = i + other_keys.size();
        if (size >= std::numeric_limits<sidx_t>::max()) {
            throw std::runtime_error("Too big system for sidx_t");
        } else {
            size = size_;
        }
        for (const auto& sk : other_keys){
            other_compositions_vec.push_back(other_compositions.at(sk));
        }
    }
    sidx_t other_index(const std::string &key) const {
        return index_of_check(key, other_keys) + strides.back();
    }
    pidx_t pki(char pk) const {
        for (int i=0; i<ndim; ++i){
            if (principal_keys[i] == pk) return i;
        }
        throw std::logic_error("unkown principal key");
    }
    lvl_val_t value(const char pk, const std::vector<lvl_idx_t>& pidxs) const
    {
        const auto pi = pki(pk);
        const auto lvl_i = pidxs[pi];
        const auto &lvls = levels_vec[pi];
        if (lvl_i >= lvls.size())
            throw std::runtime_error("level index out of range");
        return lvls[lvl_i];
    }
private:
    int principal_index(sidx_t si, pidx_t pi) const {
#if !defined(NDEBUG)
        if (pi > ndim or si > strides.back()){
            throw std::runtime_error("indices out of bounds");
        }
#endif
        const auto result = (si % strides[pi+1]) / strides[pi];
        return result;
    }
public:
    template <typename T>
    lvl_idx_t floor_val(pidx_t pi, T val) const {
        const auto lvls = levels_vec[pi];
        for (lvl_idx_t i = dim[pi]; i > 0; --i){
            if (lvls[i-1] <= val){
                return i-1;
            }
        }
        throw std::runtime_error("floor_val failed");
    }
    template <typename T>
    lvl_idx_t ceil_val(pidx_t pi, T val) const {
        const auto lvls = levels_vec[pi];
        for (lvl_idx_t i = 0; i < dim[pi]; ++i){
            if (lvls[i] >= val)
                return i;
        }
        throw std::runtime_error("ceil_val failed");
    }
    sidx_t principals_to_globidx(const std::vector<lvl_idx_t> pidxs) const {
        sidx_t result = 0;
        for (pidx_t pi=0; pi < ndim; ++pi){
#if !defined(NDEBUG)
            if (pidxs[pi] >= levels_vec[pi].size()){
                std::runtime_error("One of the principal indices too large");
            }
#endif
            result += strides[pi]*pidxs[pi];
        }
#if !defined(NDEBUG)
        if (result >= strides.back()){
            throw std::runtime_error("Bug");
        }
#endif
        return result;
    }
    void globidx_to_principals(std::vector<lvl_idx_t> &out, const sidx_t globidx) const {
#if !defined(NDEBUG)
        if (globidx >= strides.back()){
            throw std::runtime_error("Invalid globidx");
        }
        if (out.size() != ndim){
            throw std::runtime_error("globidx_to_principals out parameter of wrong size");
        }
#endif
        for (pidx_t pi=0; pi < ndim; ++pi){
            out[pi] = principal_index(globidx, pi);
        }
    }
    std::vector<lvl_idx_t> principal_indices(const std::unordered_map<char, lvl_val_t> &values) const {
        std::vector<lvl_idx_t> result;
        for (int pi=0; pi<ndim; ++pi){
            const auto &pk = principal_keys[pi];
            const auto &val = values.at(pk);
            const auto &lvls = levels_vec[pi];
            for (lvl_idx_t idx=0; idx < lvls.size(); ++idx){
                if (lvls[idx] == val){
                    result.push_back(idx);
                    goto found_index;
                }
            }
            throw std::runtime_error("Index not found");
        found_index:
            ;
        }
        return result;
    }
    std::vector<lvl_val_t> principal_values(const std::vector<lvl_idx_t> &indices) const {
        std::vector<lvl_val_t> result; result.reserve(ndim);
        for (int pi=0; pi < ndim; ++pi){
            result.push_back(levels_vec[pi][indices[pi]]);
        }
        return result;
    }
    std::unordered_map<int, int> composition(const std::vector<lvl_val_t> &values) const {
        std::unordered_map<int, int> result;
        for (pidx_t pi=0; pi<ndim; ++pi){
            for (const auto &[k, coeffs] : principal_composition_polys[pi]){
                int xn = 1;
                for (const auto &coeff : coeffs){
                    result[k] += coeff*xn;
                    xn *= values[pi];
                }
            }
        }
        return result;
    }
    std::unordered_map<int, int> composition(const sidx_t globidx) const {
        if (globidx < strides.back()) {
            std::vector<lvl_idx_t> pidxs(ndim);
            globidx_to_principals(pidxs, globidx);
            return composition(principal_values(pidxs));
        } else {
            return other_compositions_vec[globidx - strides.back()];
        }
    }
    template <typename T=lvl_val_t>
    o_umap_id balance(const double ntot_prod, const std::vector<T> &tgt) const {
        if (tgt.size() != ndim) throw std::runtime_error("incorrect size of tgt");
        std::vector<std::tuple<pidx_t, double, double>> plu;
        umap_id result;
        std::vector<lvl_idx_t> low_principals; low_principals.reserve(ndim);
        for (int pi=0; pi<ndim; ++pi){
            const auto v = tgt[pi];
            const auto lvls = levels_vec[pi];
	    if (v < lvls.front() or v > lvls.back()){
	        return std::nullopt;
	    }
            const auto l = floor_val<T>(pi, v);
            const auto u = ceil_val<T>(pi, v);
            low_principals.push_back(l);
            if (l == u){
                continue;
            } else if (u - l != 1) {
                throw std::logic_error("ceil_val & floor_val differ by more than 1");
            }
            const auto dl = v - lvls[l];
            const auto du = lvls[u] - v;
            const double span = dl + du;
            plu.emplace_back(std::tuple<pidx_t, double, double>{pi, du/span, dl/span});
        }
        const auto nmix = plu.size();
        if (nmix == 0) {
            result[principals_to_globidx(low_principals)] = ntot_prod;
        } else {
            const auto nprod = 1u << nmix;
            for (unsigned i=0; i<nprod; ++i){
                auto cur_principals = low_principals; // copy
                double c = ntot_prod;
                for (pidx_t wi=0; wi<nmix; ++wi){
                    const pidx_t pi = std::get<0>(plu[wi]);
                    const bool hi = (i >> wi & 1);
                    c *= hi ? std::get<2>(plu[wi]) : std::get<1>(plu[wi]);
                    cur_principals[pi] += hi;
                }
                result[principals_to_globidx(cur_principals)] = c;
            }
        }
        return result;
    }
    o_umap_id balance_delta(const double ntot_prod, const std::vector<lvl_idx_t> &from,
                            const std::unordered_map<char, lvl_val_t> &delta) const{
        auto tgt = principal_values(from);
        for (const auto& [k, v] : delta) {
            tgt[pki(k)] += v;
        }
        return balance<lvl_val_t>(ntot_prod, tgt);
    }
    std::string name(sidx_t index) const {
        if (index < strides.back()){
            std::string nam;
            std::vector<lvl_idx_t> indices(ndim);
            globidx_to_principals(indices, index);
            const auto values = principal_values(indices);
            for (pidx_t pi=0; pi < ndim; ++pi){
                nam += principal_keys[pi];
                nam += std::to_string(values[pi]);
            }
            return nam;
        } else {
            return other_keys[index - strides.back()];
        }
    }
    std::string reaction_string(const Reac &reac_prod_rate) const {
        const auto &[reac, prod, rate] = reac_prod_rate;
        auto coeff = [&](double v) -> std::string {
            if (v == 1) {
                return "";
            } else if (v == (int)v){
                return std::to_string((int)v) + " ";
            } else {
                return std::to_string(v) + " ";
            }
        };
        std::string s_reac;
        for (const auto &[k, v] : reac){
            s_reac += coeff(v) + name(k) + " + ";
        }
        if (s_reac.size()){
            s_reac = s_reac.substr(0, s_reac.size() - 3);
        }
        std::string s_prod;
        for (const auto &[k, v] : prod){
            s_prod += coeff(v) + name(k) + " + ";
        }
        if (s_prod.size()){
            s_prod = s_prod.substr(0, s_prod.size() - 3);
        }
        return s_reac + " -> " + s_prod + "; " + std::to_string(rate);
    }
    void check_reaction(const Reac &reac_prod_rate) const {
        const auto &[reac, prod, rate] = reac_prod_rate;
        if (rate < 0){
            throw std::runtime_error("Negative rate encountered");
        }
        std::unordered_map<int, double> net_compo;
#define ENSEMBLE_ADD_TO_NET(CONT, SIGN)                         \
        for (const auto [sp_idx, sp_amount] : CONT) {           \
            for (const auto &[ck, cm] : composition(sp_idx)) {  \
                net_compo[ck] += SIGN*cm*sp_amount;             \
            }                                                   \
        }
        ENSEMBLE_ADD_TO_NET(reac, -1);
        ENSEMBLE_ADD_TO_NET(prod, 1);
#undef ENSEMBLE_ADD_TO_NET
        for (const auto &[c_key, c_amount] : net_compo){
            if (std::abs(c_amount) > this->check_abstol) {
                std::cerr << reaction_string(reac_prod_rate) << std::endl;
                throw std::runtime_error(StreamFmt() << "Composition violation for: " << c_key);
            }
        }
    }
    void check_reactions(const ReacVec &rxns) const {
        for (const auto &reac_prod_rate : rxns){
            check_reaction(reac_prod_rate);
        }
    }
    kinetgine::CRN
    as_crn(Collection::ReacVec rxns,
           kinetgine::CRN::vec_s substances={},
           kinetgine::CRN::vec_rN elementary_reactions={},
           kinetgine::CRN::vec_rA arbitrary_reactions={},
           bool identify_quads=false
        )
    {
        kinetgine::CRN::vec_s new_substances;
        for (sidx_t glob_idx = 0; glob_idx < size; ++glob_idx){
            const auto rev_idx = size - glob_idx - 1;
            auto cur = kinetgine::substance(name(rev_idx), composition(rev_idx));
            if (std::find_if(substances.begin(), substances.end(), [&](const auto &s){ return *s == *cur; }) == substances.end()){
                new_substances.emplace_back(std::move(cur));
            }
        }
        substances.insert(substances.end(), new_substances.begin(), new_substances.end());
        elementary_reactions.reserve(elementary_reactions.size() + rxns.size());
        for (unsigned i=0; i < rxns.size(); ++i){
            const auto &[reac, prod, rate] = rxns[i];
            umap_sd net;
            umap_si act;
            for (const auto &[k, v]: reac){
                net[name(k)] -= v;
                act[name(k)] += v;
            }
            for (const auto &[k, v]: prod){
                net[name(k)] += v;
            }
            expel(net, [](double &v) -> bool { return v == 0; });
            elementary_reactions.emplace_back(kinetgine::reac_n(net, act, rate));
        }
        return kinetgine::CRN{substances, elementary_reactions, arbitrary_reactions, identify_quads};
    }
#if !defined(ENSEMBLE_TESTING)
private:
#endif
    std::vector<lvl_idx_t> dim;
    std::vector<umap_i_vi> principal_composition_polys;
    // umap_c_vi levels;

    vec_uii other_compositions_vec;
};


struct Rule {
    double ntot_prod;
    umap_si other_reac;
    umap_sd other_prod;
    std::unordered_map<char, Collection::lvl_val_t> delta;
    unsigned polymer_order=1;
};

Collection::ReacVec mkR(
    const Collection& coll,
    const Rule& rule,
    const std::vector<std::vector<Collection::lvl_idx_t>>& froms,
    const double rate)
{
    Collection::ReacVec result;
    umap_ii reac;
    umap_id prod;
    if (rule.polymer_order != froms.size()){
        throw std::runtime_error("polymer_order and number of sources inconsistent");
    }
    bool ok = true;
    for (const auto& from : froms){
        auto o_prod = coll.balance_delta(rule.ntot_prod, from, rule.delta);
        if (!o_prod){
            ok = false;
            break;
        }
        reac[coll.principals_to_globidx(from)] += 1;
        for (const auto &[k, v] : o_prod.value()){
            prod[k] += v;
        }
    }
    if (not ok){
        return {};
    }
#define ENSEMBLE_ADD_OTHER(KEY)                         \
    for (const auto &[k, v] : rule.other_ ## KEY){      \
        KEY[coll.other_index(k)] += v;                  \
    }
    ENSEMBLE_ADD_OTHER(prod);
    ENSEMBLE_ADD_OTHER(reac);
#undef ENSEMBLE_ADD_OTHER
    result.reserve(1);
    auto reac_prod_rate = std::make_tuple(reac, prod, rate);
    coll.check_reaction(reac_prod_rate);
    result.emplace_back(reac_prod_rate);
    return result;
}

std::vector<double> mul_vec(const std::vector<double> &inp, const double fact)
{
    std::vector<double> nv;
    nv.reserve(inp.size());
    std::transform(inp.begin(), inp.end(), std::back_inserter(nv),
                   [&](const double e){ return e*fact; });
    return nv;
}

template<typename From, typename To>
std::vector<To> conv_vec(const std::vector<From> &inp) {
    std::vector<To> nv;
    nv.reserve(inp.size());
    std::transform(inp.begin(), inp.end(), std::back_inserter(nv),
                   [&](const From &v) -> To { return v; });
    return nv;
}

template<typename KeyFrom, typename KeyTo, typename ValueFrom, typename ValueTo>
std::unordered_map<KeyTo, ValueTo> conv_umap(const std::unordered_map<KeyFrom, ValueFrom> &inp) {
    std::unordered_map<KeyTo, ValueTo> nu;
    for (const auto &[k, v] : inp){
        nu[k] = v;
    }
    return nu;
}

Collection::ReacVec mkFrag(
    const Collection& coll,
    const std::vector<Collection::lvl_idx_t>& from,
    const double rate)
{
    Collection::ReacVec result;
    if (from[coll.pki('o')] == 0){ // lowest oxyl count (0 presumably)
        return result;
    }
    result.reserve(1);
    umap_id prod;
    umap_ii reac {{coll.principals_to_globidx(from), 1}};
    const auto pi_n = coll.pki('n');
    const auto from_idx_n = from[pi_n];
    const auto &ns = coll.levels_vec[pi_n];
    const auto from_val_n = ns[from_idx_n];
    if (from_idx_n == 0) { // lowest number of segments
        prod = conv_umap<int, int, int, double>(reac);
        const auto idx_tiny = coll.other_index("tiny" + std::to_string(from_val_n));
        prod[idx_tiny] = 1;
    } else {
        const auto pi_o = coll.pki('o');
        const auto pi_a = coll.pki('a');
        const double span = from_val_n - ns[0];
        for (int i=0; i<from_idx_n; ++i){
            const double cur_n = ns[i];
            const double next_n = ns[i+1];
            const double w = (next_n - cur_n)/span;
            const double ratio_n1 = cur_n/from_val_n;
            const double ratio_n2 = 1.0 - ratio_n1;
            auto tgt = conv_vec<Collection::lvl_val_t, double>(coll.principal_values(from));
            tgt[pi_o] -= 1.0;
            tgt[pi_a] += 1.0;
            auto o_bal1 = coll.balance<double>(w, mul_vec(tgt, ratio_n1));
            auto o_bal2 = coll.balance<double>(w, mul_vec(tgt, ratio_n2));
            if (!o_bal1 or !o_bal2){
                return {};
            } else {
            // if (o_bal1 and o_bal2){
                for (auto &o_bal : {o_bal1, o_bal2}){
                    for (const auto &[k, v] : o_bal.value()){
                        prod[k] += v;
                    }
                }
                const auto idx_keto = coll.other_index("keto" + std::to_string(from_val_n));
                prod[idx_keto] += w;
            }
        }
    }
    auto reac_prod_rate = std::make_tuple(reac, prod, rate);
    coll.check_reaction(reac_prod_rate);
    result.emplace_back(reac_prod_rate);
    return result;
}


Collection::ReacVec mkAlkylInter(
    const Collection& coll,
    const std::vector<Collection::lvl_idx_t> &from1,
    const std::vector<Collection::lvl_idx_t> &from2,
    const double rate)
{
    Collection::ReacVec result;
    auto to = coll.principal_values(from1);
    const auto val2 = coll.principal_values(from2);
    if (to[coll.pki('a')] == 0 or val2[coll.pki('a')] == 0) {
        return {};
    }
    for (int pi=0; pi<coll.ndim; ++pi){
        to[pi] += val2[pi];
    }
    to[coll.pki('a')] -= 2;
    auto o_prod = coll.balance(1.0, to);
    if (o_prod){
        result.reserve(1);
        umap_ii reac {};
        reac[coll.principals_to_globidx(from1)] += 1;
        reac[coll.principals_to_globidx(from2)] += 1;
        auto &prod = o_prod.value();
        auto reac_prod_rate = std::make_tuple(reac, prod, rate);
        coll.check_reaction(reac_prod_rate);
        result.emplace_back(reac_prod_rate);
    }
    return result;
}

Collection::ReacVec mkOxylInterH(
    const Collection& coll,
    const std::vector<Collection::lvl_idx_t> &from1,
    const std::vector<Collection::lvl_idx_t> &from2,
    const double rate)
{
    Collection::ReacVec result;
    umap_ii reac;
    umap_id prod;
    reac[coll.principals_to_globidx(from1)] += 1;
    reac[coll.principals_to_globidx(from2)] += 1;
    if (from1 == from2){
        auto o_bal1 = coll.balance_delta(1, from1, {{'o', -1}});
        auto o_bal2 = coll.balance_delta(1, from1, {{'a', 1}});
        if (o_bal1 and o_bal2){
            prod = {
                {coll.other_index("hydroxy" + std::to_string(coll.value('n', from1))), 1}
            };
            for (auto &o_bal : {o_bal1, o_bal2}){
                for (const auto &[k, v] : o_bal.value()){
                    prod[k] += v;
                }
            }
            auto reac_prod_rate = std::make_tuple(reac, prod, rate);
            coll.check_reaction(reac_prod_rate);
            result.emplace_back(reac_prod_rate);
        }
    } else {
        for (const auto &[from_o, from_a] : {std::tuple{from1, from2}, std::tuple{from2, from1}}){
            if (coll.value('o', from_o) > 0) {
                prod = {
                    {coll.other_index("hydroxy" + std::to_string(coll.value('n', from_o))), 1}
                };
                auto o_bal1 = coll.balance_delta(1, from_o, {{'o', -1}});
                auto o_bal2 = coll.balance_delta(1, from_a, {{'a', 1}});
                if (o_bal1 and o_bal2){
                    for (auto &o_bal : {o_bal1, o_bal2}){
                        for (const auto &[k, v] : o_bal.value()){
                            prod[k] += v;
                        }
                    }
                    auto reac_prod_rate = std::make_tuple(reac, prod, rate);
                    coll.check_reaction(reac_prod_rate);
                    result.emplace_back(reac_prod_rate);
                }
            }
        }
    }
    return result;
}


struct Rates {
    double oxyl, oxyl_H_intra, oxyl_H_inter, peroxy_inter,
        alkyl_inter, alkyl_O2, radrad_intra, OH, H;
};


struct Formulation {
    double n0, monomers_per_chain, N_scal_intra, w_scal_OH, disprop_frac, radrad_scal_intra,
        loop_scaling_inter=0.0, loop_scaling_intra=0.0;
    Rates base_rates;
    bool max_loop_prob = false;
    using ps_t = std::vector<Collection::lvl_idx_t>;

    double monomers_n(const Collection& coll, const ps_t &pidxs)
    {
        return coll.value('n', pidxs)/n0*monomers_per_chain;
    }

    double inter_length_scaling(const Collection& coll, const ps_t &pidxs)
    {
        // From Fig 3, p 478, Bartoszek et al. doi:10.1002/kin.20575
        return std::pow(coll.value('n', pidxs)/n0, w_scal_OH);
    }

    double intra_length_scaling(const Collection& coll, const ps_t &pidxs) {
        return std::pow(coll.value('n', pidxs)/n0, N_scal_intra);
    }

    double intra_radical_scaling(const Collection& coll, const ps_t &pidxs, char pk) {
        return std::pow(coll.value(pk, pidxs)/2.0, radrad_scal_intra);
    }

    double loop_stiffness(const Collection& coll, const ps_t &pidxs, double exponent) {
        const auto n_loops = coll.value('l', pidxs);
        //const double loops_over_monomers = n_loops / monomers_n(coll, pidxs);
        if (n_loops == 0) {  // loops_over_monomers == 0
            return 1;
        } else {
            return std::pow(n_loops, exponent);  // loops_over_monomers
        }
    }

    double k1OH(const Collection& coll, const ps_t &pidxs){
        return monomers_n(coll, pidxs) * inter_length_scaling(coll, pidxs) * base_rates.OH;
    }
    double k1H(const Collection& coll, const ps_t &pidxs){
        return monomers_n(coll, pidxs) * inter_length_scaling(coll, pidxs) * base_rates.H;
    }
    double k2AlkylInter_(const Collection& coll, const ps_t &pidxs1, const ps_t &pidxs2, double frac){
        const double s1 = inter_length_scaling(coll, pidxs1);
        const double s2 = inter_length_scaling(coll, pidxs2);
        const double l1 = loop_stiffness(coll, pidxs1, loop_scaling_inter);
        const double l2 = loop_stiffness(coll, pidxs1, loop_scaling_inter);
        const double a1 = coll.value('a', pidxs1);
        const double a2 = coll.value('a', pidxs2);
        return s1*s2 * l1*l2 * a1*a2 * frac * base_rates.alkyl_inter;
    }
    double k2AlkylInter(const Collection& coll, const ps_t &pidxs1, const ps_t &pidxs2){
        return k2AlkylInter_(coll, pidxs1, pidxs2, 1.0 - disprop_frac);
    }
    double k1AlkylO2(const Collection& coll, const ps_t &pidxs){
        return coll.value('a', pidxs)*base_rates.alkyl_O2;
    }
    double k2PeroxyInter(const Collection& coll, const ps_t &pidxs1, const ps_t &pidxs2){
        const double s1 = inter_length_scaling(coll, pidxs1);
        const double s2 = inter_length_scaling(coll, pidxs2);
        const double l1 = loop_stiffness(coll, pidxs1, loop_scaling_inter);
        const double l2 = loop_stiffness(coll, pidxs1, loop_scaling_inter);
        const double p1 = coll.value('p', pidxs1);
        const double p2 = coll.value('p', pidxs2);
        return s1*s2 * l1*l2 * p1*p2 * base_rates.peroxy_inter;
    }
    double k1PeroxyIntra(const Collection& coll, const ps_t &pidxs){
        return intra_radical_scaling(coll, pidxs, 'p') * intra_length_scaling(coll, pidxs) *
            loop_stiffness(coll, pidxs, loop_scaling_intra) * base_rates.radrad_intra;
    }
    double loop_prob(const Collection& coll, const ps_t &pidxs){
        const double l = coll.value('l', pidxs);
        if (max_loop_prob)
            return (l == 0) ? 0.0 : 1.0;
        else
            return l/(2+l); // <|distance|> for 2 points on unit line == 1/3;
    }
    double five_prob() {
        return 3./5.;
    }
    double k1Frag(const Collection& coll, const ps_t &pidxs){
        return coll.value('o', pidxs) * (1 - loop_prob(coll, pidxs)) * (1 - five_prob()) * base_rates.oxyl;
    }
    double k1Damage(const Collection& coll, const ps_t &pidxs){
        return coll.value('o', pidxs) * five_prob() * base_rates.oxyl; // on five-ring
    }
    double k1Loop_(const Collection& coll, const ps_t &pidxs, double frac){
        auto irs = intra_radical_scaling(coll, pidxs, 'a');
        auto ils = intra_length_scaling(coll, pidxs);
        auto ls = loop_stiffness(coll, pidxs, loop_scaling_intra);
        return irs * ils * ls * frac * base_rates.radrad_intra;
    }
    double k1Loop(const Collection& coll, const ps_t &pidxs){
        return k1Loop_(coll, pidxs, 1.0 - disprop_frac);
    }
    double k1LoopOpen(const Collection& coll, const ps_t &pidxs){
        return coll.value('o', pidxs) * loop_prob(coll, pidxs) * (1 - five_prob()) * base_rates.oxyl;
    }
    double k1DispropIntra(const Collection& coll, const ps_t &pidxs){
        return k1Loop_(coll, pidxs, disprop_frac);
    }
    double k2DispropInter(const Collection& coll, const ps_t &pidxs1, const ps_t &pidxs2)
    {
        return k2AlkylInter_(coll, pidxs1, pidxs2, disprop_frac);
    }
    double k1OxylIntraH(const Collection& coll, const ps_t &pidxs)
    {
        return coll.value('n', pidxs) * loop_stiffness(coll, pidxs, loop_scaling_intra) * base_rates.oxyl_H_intra;
    }
    double k2OxylInterH(const Collection& coll, const ps_t &pidxs1, const ps_t &pidxs2)
    {
        const double n1 = coll.value('n', pidxs1);
        const double n2 = coll.value('n', pidxs2);
        const double o1 = coll.value('o', pidxs1);
        const double o2 = coll.value('o', pidxs2);
        return n0 * (o1/n1 + o2/n2) * base_rates.oxyl_H_inter;
    }
    Collection::ReacVec generate_reactions(const Collection& coll, bool check=true) {
        Collection::ReacVec reac_vec;
        auto extend = [&](const auto& nv) -> void {
            // The line below give super-linear complexity:
            //reac_vec.reserve(reac_vec.size() + std::distance(nv.begin(), nv.end()));
            reac_vec.insert(reac_vec.end(), nv.begin(), nv.end());
        };
        const Rule rOH{1.0, {{"OH", 1}}, {{"H2O", 1}}, {{'a', 1}}};
        const Rule rH {1.0, {{"H",1 }}, {{"H2", 1}}, {{'a', 1}}};
        const Rule rAlkylO2{1.0, {{"O2", 1}}, {}, {{'p', 1}, {'a', -1}}};
        const Rule rPeroxyIntra{1.0, {}, {{"O2", 1}}, {{'o', 2}, {'p', -2}}};
        Rule rDamage{1.0, {}, {}, {{'a', 1}, {'o', -1}}};
        const Rule rLoop{1.0, {}, {}, {{'a', -2}, {'l', 1}}};
        Rule rLoopOpen{1.0, {}, {}, {{'l', -1} , {'a', 1}, {'o', -1}}};
        Rule rDispropIntra{1.0, {}, {}, {{'a', -2}}};
        Rule rOxylIntraH{1.0, {}, {}, {{'a', 1}, {'o', -1}}};

        const Rule rPeroxyInter{1.0, {}, {{"O2", 1}}, {{'p', -1}, {'o', 1}}, 2};
        Rule rDispropInter{1.0, {}, {}, {{'a', -1}}, 2};

        std::vector<Collection::lvl_idx_t> pidxs(coll.ndim);
        std::vector<Collection::lvl_idx_t> pidxs2(coll.ndim);
#if defined(ENSEMBLE_VERBOSE)
        std::cout << "Will generate reactions for " << coll.size << " substances" << std::endl;
        std::cout << std::fixed << std::setprecision(3);
#endif
        for (int glob_idx = 0; glob_idx < coll.strides.back(); ++glob_idx){
#if defined(ENSEMBLE_VERBOSE)
            if (glob_idx % (coll.strides.back() / 1000 + 1) == 0 or glob_idx + 1 == coll.strides.back()){
                std::cout << "\rPercent done: " << glob_idx*1e2/(coll.strides.back()-1)
                          << " %, number of reactions: " << reac_vec.size() << std::flush;
            }
#endif
            coll.globidx_to_principals(pidxs, glob_idx);
            extend(mkR(coll, rOH, {pidxs}, k1OH(coll, pidxs)));
            extend(mkR(coll, rH, {pidxs}, k1H(coll, pidxs)));
            extend(mkR(coll, rAlkylO2, {pidxs}, k1AlkylO2(coll, pidxs)));
            extend(mkR(coll, rPeroxyIntra, {pidxs}, k1PeroxyIntra(coll, pidxs)));
            const auto cur_n1 = coll.value('n', pidxs);
            const auto cur_n1_s = std::to_string(cur_n1);
            rDamage.other_prod = {
                {"keto" + cur_n1_s, 1.0},
                {"damage_fivering" + cur_n1_s, 1.0}
            };
            extend(mkR(coll, rDamage, {pidxs}, k1Damage(coll, pidxs)));
            extend(mkR(coll, rLoop, {pidxs}, k1Loop(coll, pidxs)));
            rLoopOpen.other_prod = {
                {"keto" + cur_n1_s, 1.0},
                {"scission_loop" + cur_n1_s, 1.0}
            };
            extend(mkR(coll, rLoopOpen, {pidxs}, k1LoopOpen(coll, pidxs)));
            rDispropIntra.other_prod = {
                {"unsaturated" + cur_n1_s, 1.0}
            };
            extend(mkR(coll, rDispropIntra, {pidxs}, k1DispropIntra(coll, pidxs)));
            rOxylIntraH.other_prod = {
                {"hydroxy" + cur_n1_s, 1.0}
            };
            extend(mkR(coll, rOxylIntraH, {pidxs}, k1OxylIntraH(coll, pidxs)));

            extend(mkFrag(coll, pidxs, k1Frag(coll, pidxs)));

            for (int glob_idx2 = 0; glob_idx2 <= glob_idx; ++glob_idx2){
                coll.globidx_to_principals(pidxs2, glob_idx2);
                extend(mkAlkylInter(coll, pidxs, pidxs2, k2AlkylInter(coll, pidxs, pidxs2)));
                extend(mkR(coll, rPeroxyInter, {pidxs, pidxs2}, k2PeroxyInter(coll, pidxs, pidxs2)));
                const auto cur_n2 = coll.value('n', pidxs2);
                const auto cur_n2_s = std::to_string(cur_n2);
                rDispropInter.other_prod = {};
                rDispropInter.other_prod["unsaturated" + cur_n1_s] += 0.5;
                rDispropInter.other_prod["unsaturated" + cur_n2_s] += 0.5;
                extend(mkR(coll, rDispropInter, {pidxs, pidxs2}, k2DispropInter(coll, pidxs, pidxs2)));
                extend(mkOxylInterH(coll, pidxs, pidxs2, k2OxylInterH(coll, pidxs, pidxs2)));
            }
        }
#if defined(ENSEMBLE_VERBOSE)
        std::cout << std::endl;
#endif
        if (check){
#if defined(ENSEMBLE_VERBOSE)
            std::cout << "Running a second pass of composition violation checks..." << std::flush;
#endif
            coll.check_reactions(reac_vec);
#if defined(ENSEMBLE_VERBOSE)
            std::cout << " done" << std::endl;
#endif
        }
        return reac_vec;
    }
};

struct PrimaryData {
%for k, v in PRIMARY_DATA.items():
    double ${k}=${v};
%endfor
    bool max_loop_prob=false;
};

Formulation create_formulation(const PrimaryData& primary_data, double n0){
    double k_OH_M_s = primary_data.k_OH_1kDa_M_s *
        std::pow(primary_data.MW_kDa, primary_data.w_scal_OH);
    Formulation form;
    form.n0                 = n0;
    form.monomers_per_chain = primary_data.MW_kDa*1e3/primary_data.MW_monomer_g_mol;
    form.N_scal_intra       = primary_data.N_scal_intra;
    form.w_scal_OH          = primary_data.w_scal_OH;
    form.disprop_frac       = primary_data.disprop_frac;
    form.radrad_scal_intra  = primary_data.radrad_scal_intra;
    form.base_rates.oxyl         = primary_data.k_oxyl_frag_s;
    form.base_rates.oxyl_H_intra = primary_data.k_oxyl_H_intra_s;
    form.base_rates.oxyl_H_inter = std::sqrt(k_OH_M_s*primary_data.alkyl_inter_M_s);  // geometric mean
    form.base_rates.peroxy_inter = primary_data.peroxy_inter_M_s;
    form.base_rates.alkyl_inter  = primary_data.alkyl_inter_M_s;
    form.base_rates.alkyl_O2     = primary_data.k_alkyl_O2_M_s;
    form.base_rates.radrad_intra = std::pow(primary_data.MW_kDa/94.0, form.N_scal_intra) *
        M_LN2/primary_data.t_half_loop_closure_intra_s_94kDa;
    form.base_rates.OH           = k_OH_M_s;
    form.base_rates.H            = k_OH_M_s/primary_data.ratio_OH_H;
    form.loop_scaling_intra  = primary_data.loop_scaling_intra;
    form.loop_scaling_inter  = primary_data.loop_scaling_inter;
    form.max_loop_prob       = primary_data.max_loop_prob;
    return form;
}

std::array<std::vector<int>, 5> napol_base2(const std::array<int, 5> napol)
{
    std::array<std::vector<int>, 5> result;
    for (unsigned i=0; i<5; ++i){
        if (i) {
            result[i].push_back(0); // != 'n'
        }
        for (int j=0; j < napol[i]; ++j){
            result[i].push_back(1u << j);
        }
    }
    return result;
}

Collection mk_collection(const std::array<std::vector<int>, 5> &napol){
    const auto prim_k = std::vector<char>{'n', 'a', 'p', 'o', 'l'};
    std::unordered_map<char, std::vector<int>> lvls;
    for (unsigned i=0; i<prim_k.size(); ++i){
        if (!napol[i].size()){
            throw std::runtime_error("Need at least one level");
        }
        if (i == 0){
            if (napol[i][0] < 1)
                throw std::runtime_error("Need at least one segment");
        }
        lvls[prim_k[i]] = napol[i];
    }
    const auto comp_n = std::unordered_map<int, std::vector<int>>{{1, {2, 2}}, {6, {0, 1}}};
    const auto comp_a = std::unordered_map<int, std::vector<int>>{{1, {0, -1}}};
    const auto comp_p = std::unordered_map<int, std::vector<int>>{{1, {0, -1}}, {8, {0, 2}}};
    const auto comp_o = std::unordered_map<int, std::vector<int>>{{1, {0, -1}}, {8, {0, 1}}};
    const auto comp_l = std::unordered_map<int, std::vector<int>>{{1, {0, -2}}};

    const auto prim_comp = std::unordered_map<char, std::unordered_map<int, std::vector<int>>>{
        {'n', comp_n}, {'a', comp_a}, {'p', comp_p}, {'o', comp_o}, {'l', comp_l}
    };
    auto other_subs = std::vector<std::string>{"H2O", "OH", "O2", "H2", "H"};
    const auto comp_H2O = std::unordered_map<int, int>{{1, 2}, {8, 1}};
    const auto comp_OH = std::unordered_map<int, int>{{1, 1}, {8, 1}};
    const auto comp_O2 = std::unordered_map<int, int>{{8, 2}};
    const auto comp_H2 = std::unordered_map<int, int>{{1, 2}};
    const auto comp_H = std::unordered_map<int, int>{{1, 1}};
    auto other_comp = std::unordered_map<std::string, std::unordered_map<int, int>>{
        {"H2O", comp_H2O}, {"OH", comp_OH}, {"O2", comp_O2}, {"H2", comp_H2}, {"H", comp_H}
    };
    auto reg = [&](const std::string& name, Collection::lvl_val_t lvl_val, std::unordered_map<int, int>&& compo) -> void {
        const auto key = name + std::to_string(lvl_val);
        other_comp[key] = compo;
        other_subs.push_back(key);
    };
    for (auto nn : lvls['n']){
        reg("unsaturated", nn, {{1, -2}});
        reg("keto", nn, {{1, -2}, {8, 1}});
        reg("hydroxy", nn, {{8, 1}});
        reg("damage_fivering", nn, {{1, 2}});  // counter-balance 'keto'
        reg("scission_loop", nn, {});
        reg("tiny", nn, {});
    }
    return Collection(prim_k, lvls, prim_comp, other_subs, other_comp);
}


NAMESPACE_END(ensemble)

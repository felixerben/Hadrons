/*
 * DMixingTopD.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2025
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Matteo Di Carlo <matteo.dicarlo@cern.ch>
 * Author: Felix Erben <felix.erben@cern.ch>
 * Author: Raoul Hodgson <raoul.hodgson@desy.de>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution
 * directory.
 */

/*  END LEGAL */
#ifndef Hadrons_MContraction_DMixingTopD_hpp_
#define Hadrons_MContraction_DMixingTopD_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>
#include <Hadrons/Modules/MContraction/DMixingUtils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         DMixingTopD                                        *
 *                    (Fig. 4 (D) in arxiv:2504.16189)
 *             qCL   ┌───┐                 qUR
 *         /---->----| r |----------------->--------------\
 *        /          └───┘                                 \
 *       /          /     \          qLoop2 /-<-\           \
 *   g5 *           \     /                /     \           * g5
 *       \           \->-/ qLoop1          \     /          /
 *        \                                 ┌───┐          /
 *         \--------------<-----------------| r'|----<----/
 *                       qUL                └───┘  qCR
 *  tsrc               t1                    t2              tsnk
 *
 * Four configurations for the two weak Hamiltonians M_r M_{r'}
 * (cf. Fig. 3 in arxiv:2504.16189)
 *
 * p = +: GA x GB =   V x V + A x A
 * p = -: GA x GB = - A x V - V x A
 *
 * Contractions: [...] = tr(...)
 * rr' = 11:
 *  [qCR * g5 * qUR * GB1 * qLoop1 * GA1 * qCL * g5 * gUL * GB2 * gLoop2 * GA2]
 * rr' = 12:
 *  [qCR * g5 * qUR * GB1 * qLoop1 * GA1 * qCL * g5 * gUL * GB2]*[gLoop2 * GA2]
 * rr' = 21:
 *  [qCR * g5 * qUR * GB1 * qCL * g5 * gUL * GB2 * gLoop2 * GA2]*[qLoop1 * GA1]
 * rr' = 22:
 *  [qCR * g5 * qUR * GB1 * qCL * g5 * gUL * GB2]*[gLoop2 * GA2]*[qLoop1 * GA1]
 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class DMixingTopDPar : Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DMixingTopDPar,
                                    std::string, qULeft,
                                    std::string, qCLeft,
                                    std::string, qURight,
                                    std::string, qCRight,
                                    std::string, qLoop1,
                                    std::string, qLoop2,
                                    std::string, output);
};

template <typename FImpl>
class TDMixingTopD : public Module<DMixingTopDPar>
{
public:
    FERM_TYPE_ALIASES(FImpl, );
    class Metadata : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Metadata,
                                        int        , r,
                                        int        , s,
                                        std::string, parity,
                                        int        , eta_max);
    };
    typedef Correlator<Metadata, std::vector<Complex>> Result;

    using Parity   = typename DMixingUtils<FImpl>::Parity;
    using OpStruct = typename DMixingUtils<FImpl>::OpStruct;

    const std::array<Gamma,8> &GHs     = DMixingUtils<FImpl>::GHs;
    const std::array<Gamma,2> &parityG = DMixingUtils<FImpl>::parityG;

public:
    // constructor
    TDMixingTopD(const std::string name);
    // destructor
    virtual ~TDMixingTopD(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
    // bespoke subcontractions
    virtual SlicedPropagator contract_D_half(
        const PropagatorField &prop_c,
        const PropagatorField &prop_u_adj,
        const PropagatorField &loop);
    virtual std::vector<std::vector<Complex>> contract_D(
        const SlicedPropagator &half_if,
        const SlicedPropagator &half_fi);

    int half_idx(OpStruct r, int i, int Neta) {
        return i + Neta * r;
    }
};

MODULE_REGISTER_TMP(DMixingTopD, TDMixingTopD<FIMPL>, MContraction);

/******************************************************************************
 *                        TDMixingTopD implementation                         *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TDMixingTopD<FImpl>::TDMixingTopD(const std::string name)
    : Module<DMixingTopDPar>(name)
{
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TDMixingTopD<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().qULeft,
                                   par().qCLeft,
                                   par().qURight,
                                   par().qCRight};
    if (!par().qLoop1.empty())
        in.push_back(par().qLoop1);
    if (!par().qLoop2.empty())
        in.push_back(par().qLoop2);

    return in;
}

template <typename FImpl>
std::vector<std::string> TDMixingTopD<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

template <typename FImpl>
std::vector<std::string> TDMixingTopD<FImpl>::getOutputFiles(void)
{
    std::vector<std::string> output;

    if (!par().output.empty())
        output.push_back(resultFilename(par().output));

    return output;
}

template <typename FImpl>
typename TDMixingTopD<FImpl>::SlicedPropagator TDMixingTopD<FImpl>::contract_D_half(
    const PropagatorField &prop_c,
    const PropagatorField &prop_u_adj,
    const PropagatorField &loop)
{
    startTimer("mult");
    PropagatorField tmp = prop_u_adj * loop * prop_c;
    stopTimer("mult");

    SlicedPropagator out;
    startTimer("sliceSum");
    sliceSum(tmp, out, Tp);
    stopTimer("sliceSum");
    return out;
};

template <typename FImpl>
std::vector<std::vector<Complex>> TDMixingTopD<FImpl>::contract_D(
    const SlicedPropagator &half_l,
    const SlicedPropagator &half_r)
{
    Gamma g5(Gamma::Algebra::Gamma5);

    int Nt = env().getDim(Tdir);
    std::vector<std::vector<Complex>> corr(Nt, std::vector<Complex>(Nt));

    startTimer("trace");
    for (int t1 = 0; t1 < Nt; t1++)
    {
        for (int t2 = 0; t2 < Nt; t2++)
        {
            corr[t1][t2] = TensorRemove(trace(half_l[t1] * g5 * half_r[t2] * g5));
        }
    }
    stopTimer("trace");

    return corr;
};

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDMixingTopD<FImpl>::setup(void)
{
    int Neta1 = 1;
    int Neta2 = 1;
    if (!par().qLoop1.empty())
        Neta1 = (envGet(std::vector<PropagatorField *>, par().qLoop1)).size();
    if (!par().qLoop2.empty())
        Neta2 = (envGet(std::vector<PropagatorField *>, par().qLoop2)).size();
    
    const int Nt = env().getDim(Tdir);

    envTmp(std::vector<SlicedPropagator>, "half_l", 1, 2 * Neta1, SlicedPropagator(Nt));
    envTmp(std::vector<SlicedPropagator>, "half_r", 1, 2 * Neta2, SlicedPropagator(Nt));
    envTmp(std::vector<PropagatorField>, "GdsG1", 1, 2, PropagatorField(env().getGrid()));
    envTmp(std::vector<PropagatorField>, "GdsG2", 1, 2, PropagatorField(env().getGrid()));

    envCreate(HadronsSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDMixingTopD<FImpl>::execute(void)
{
    LOG(Message) << "Computing D-meson mixing diagram, topology D" << std::endl;
    LOG(Message) << "qULeft  : " << par().qULeft << std::endl;
    LOG(Message) << "qCLeft  : " << par().qCLeft << std::endl;
    LOG(Message) << "qURight : " << par().qURight << std::endl;
    LOG(Message) << "qCRight : " << par().qCRight << std::endl;
    bool qLoop1_empty = par().qLoop1.empty();
    bool qLoop2_empty = par().qLoop2.empty();
    if (qLoop1_empty)
        LOG(Message) << "Empty qLoop1 : (Pseudo)Scalar bilinear" << std::endl;
    else
        LOG(Message) << "qLoop1  : " << par().qLoop1 << std::endl;
    if (qLoop2_empty)
        LOG(Message) << "Empty qLoop2 : (Pseudo)Scalar bilinear" << std::endl;
    else
        LOG(Message) << "qLoop2  : " << par().qLoop2 << std::endl;

    std::vector<Result> result;

    const int Nt{env().getDim(Tdir)};
    GridCartesian *grid = envGetGrid(FermionField);

    std::vector<PropagatorField *> ql_default(1,nullptr);

    auto &qul = envGet(PropagatorField, par().qULeft);
    auto &qcl = envGet(PropagatorField, par().qCLeft);
    auto &qur = envGet(PropagatorField, par().qURight);
    auto &qcr = envGet(PropagatorField, par().qCRight);
    auto &ql1 = (qLoop1_empty)  ? ql_default
                                        : envGet(std::vector<PropagatorField *>, par().qLoop1);
    auto &ql2 = (qLoop2_empty)  ? ql_default
                                        : envGet(std::vector<PropagatorField *>, par().qLoop2);

    PropagatorField Unit(grid); Unit = 1;

    Gamma g5(Gamma::Algebra::Gamma5);
    const PropagatorField qur_adj = g5 * adj(qur) * g5;
    const PropagatorField qul_adj = g5 * adj(qul) * g5;

    bool same_loop_noise = (par().qLoop1 == par().qLoop2) and !qLoop1_empty and !qLoop2_empty;
    int Neta1 = ql1.size();
    int Neta2 = ql2.size();

    envGetTmp(std::vector<SlicedPropagator>, half_l);
    envGetTmp(std::vector<SlicedPropagator>, half_r);
    envGetTmp(std::vector<PropagatorField>, GdsG1);
    envGetTmp(std::vector<PropagatorField>, GdsG2);

    std::vector<std::vector<Complex>> tmpRes(Nt, std::vector<Complex>(Nt, 0.));
    std::vector<std::vector<Complex>> tmpSum(Nt, std::vector<Complex>(Nt, 0.));
    std::vector<std::vector<Complex>> diagSum(Nt, std::vector<Complex>(Nt, 0.));

    SlicedPropagator Lsum(Nt), Rsum(Nt);

    std::vector<OpStruct> r_vals;
    if (qLoop1_empty) r_vals = {OpStruct::Two};
    else                      r_vals = {OpStruct::One,OpStruct::Two};
    std::vector<OpStruct> s_vals;
    if (qLoop2_empty) s_vals = {OpStruct::Two};
    else                      s_vals = {OpStruct::One,OpStruct::Two};

    for (const auto p : {Parity::Pos,Parity::Neg})
    {
        for (int i = 0; i < std::max(Neta1,Neta2); i++)
        {
            if (i < Neta1) {
                startTimer("GH_cap");
                if (qLoop1_empty) {
                    GdsG1[OpStruct::One] = Unit * parityG[p];
                    GdsG1[OpStruct::Two] = Unit * parityG[p];
                } else {
                    DMixingUtils<FImpl>::GH_cap(GdsG1, *ql1[i], p);
                }
                stopTimer("GH_cap");
                for (const auto r : r_vals)
                    half_l[half_idx(r,i,Neta1)] = contract_D_half(qcl, qur_adj, GdsG1[r]);
            }
            if (i < Neta2) {
                if (same_loop_noise) {
                    GdsG2 = GdsG1;
                } else {
                    startTimer("GH_cap");
                    if (qLoop2_empty) {
                        GdsG2[OpStruct::One] = Unit * parityG[p];
                        GdsG2[OpStruct::Two] = Unit * parityG[p];
                    } else {
                        DMixingUtils<FImpl>::GH_cap(GdsG2, *ql2[i], p);
                    }
                    stopTimer("GH_cap");
                }
                for (const auto s : s_vals)
                    half_r[half_idx(s,i,Neta2)] = contract_D_half(qcr, qul_adj, GdsG2[s]);
            }
        }

        for (const auto r : r_vals)
        {
            for (const auto s : s_vals)
            {
                for (int t = 0; t < Nt; t++)
                {
                    Lsum[t] = Zero();
                    Rsum[t] = Zero();
                    if (same_loop_noise)
                        std::fill(diagSum[t].begin(), diagSum[t].end(), Complex(0.0));
                }

                for (int i = 0; i < std::max(Neta1,Neta2); i++)
                {
                    int eta_max = i + 1;
                    int i_l = std::min(i,Neta1-1);
                    int i_r = std::min(i,Neta2-1);

                    const auto &Li = half_l[half_idx(r,i_l,Neta1)];
                    const auto &Ri = half_r[half_idx(s,i_r,Neta2)];

                    // accumulate sums of Li, Ri
                    for (int t = 0; t < Nt; t++)
                    {
                        if (i < Neta1)
                            Lsum[t] += Li[t];
                        if (i < Neta2)
                            Rsum[t] += Ri[t];
                    }
                    tmpSum = contract_D(Lsum, Rsum);

                    if (same_loop_noise) // Guarantees Neta1==Neta2 so no edge cases to avoid
                    {
                        // accumulate sum of diagonal terms & remove from total
                        const auto diag = contract_D(Li, Ri);
                        for (int t1 = 0; t1 < Nt; t1++)
                            for (int t2 = 0; t2 < Nt; t2++)
                            {
                                diagSum[t1][t2] += diag[t1][t2];
                                tmpSum[t1][t2] -= diagSum[t1][t2];
                            }
                    }

                    const double norm = (!same_loop_noise)
                                            ? double(std::min(i+1,Neta1) * std::min(i+1,Neta2))
                                            : (i > 0 ? double((i + 1) * i) : 1.0);

                    for (int t1 = 0; t1 < Nt; t1++)
                        for (int t2 = 0; t2 < Nt; t2++)
                            tmpRes[t1][t2] = tmpSum[t1][t2] / norm;

                    Result res;
                    res.info.parity  = DMixingUtils<FImpl>::toString(p);
                    res.info.r       = DMixingUtils<FImpl>::toInt(r);
                    res.info.s       = DMixingUtils<FImpl>::toInt(s);
                    res.info.eta_max = eta_max;
                    
                    res.corr = tmpRes;
                    result.push_back(res);
                }
            }
        }
    }

    // save result, and hand it to environment
    saveResult(par().output, "DMixingTopD", result);
    auto &out = envGet(HadronsSerializable, getName());
    out = result;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_DMixingTopD_hpp_

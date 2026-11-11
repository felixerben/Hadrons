/*
 * DistilMeson2pt.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2026
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Felix Erben <felix.erben@cern.ch>
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
#ifndef Hadrons_MDistil_DistilMeson2pt_hpp_
#define Hadrons_MDistil_DistilMeson2pt_hpp_

#include <Hadrons/Global.hpp>
#include <sstream>
#include <algorithm>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DilutedNoise.hpp>
#include <Hadrons/NamedTensor.hpp>
#include <Hadrons/DistilMatrix.hpp>
#include <Hadrons/Modules/MDistil/DistilUtils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         DistilMeson2pt                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MDistil)

class DistilMeson2ptPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DistilMeson2ptPar,
        std::string,              output,
        std::string,              lapEigenPack,
        std::string,              leftNoise,
        std::string,              rightNoise,
        std::vector<std::string>, noisePairs,
        std::string,              leftPeramb,
        std::string,              rightPeramb,
        std::string,              leftVectorStem,
        std::string,              rightVectorStem,
        std::string,              dt,
        std::string,              gammaSrc,
        std::string,              gammaSnk,
        std::vector<std::string>, momenta,
        unsigned int,             blockSize,
        unsigned int,             cacheSize)
};

template <typename FImpl>
class TDistilMeson2pt: public Module<DistilMeson2ptPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    using Computation = DmfComputation<FImpl, HADRONS_DISTIL_TYPE, HADRONS_DISTIL_IO_TYPE>;
    using DistillationNoise = typename Computation::DistillationNoise;
    using Index = typename Computation::Index;
    using DistilVector = typename Computation::DistilVector;

    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
            std::string,          gammaSrc,
            std::string,          gammaSnk,
            std::string,          momSrc,
            std::string,          momSnk,
            unsigned int,         dt,
            std::vector<unsigned int>, timePartition,
            std::vector<Complex>, corr);
    };
    // constructor
    TDistilMeson2pt(const std::string name);
    // destructor
    virtual ~TDistilMeson2pt(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    std::vector<unsigned int> dt_;
    std::vector<std::vector<RealF>> momenta_;
    std::vector<Gamma::Algebra> gammaSrc_, gammaSnk_;
    std::map<Side, unsigned int> dilSizeLS_;
    bool isExact_ = false;
};

MODULE_REGISTER_TMP(DistilMeson2pt, TDistilMeson2pt<FIMPL>, MDistil);

/******************************************************************************
 *                 TDistilMeson2pt implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TDistilMeson2pt<FImpl>::TDistilMeson2pt(const std::string name)
: Module<DistilMeson2ptPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TDistilMeson2pt<FImpl>::getInput(void)
{
    std::vector<std::string> in = {
        par().lapEigenPack, par().leftNoise, par().rightNoise
    };
    if (!par().leftPeramb.empty())  in.push_back(par().leftPeramb);
    if (!par().rightPeramb.empty()) in.push_back(par().rightPeramb);
    return in;
}

template <typename FImpl>
std::vector<std::string> TDistilMeson2pt<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDistilMeson2pt<FImpl>::setup(void)
{
    if (par().leftPeramb.empty() == par().leftVectorStem.empty())
    {
        HADRONS_ERROR(Argument, "Specify exactly one of leftPeramb and leftVectorStem");
    }
    if (par().rightPeramb.empty() == par().rightVectorStem.empty())
    {
        HADRONS_ERROR(Argument, "Specify exactly one of rightPeramb and rightVectorStem");
    }
    if (par().blockSize == 0 || par().cacheSize == 0)
    {
        HADRONS_ERROR(Size, "blockSize and cacheSize must be positive");
    }
    if (par().output.empty())
    {
        HADRONS_ERROR(Argument, "Missing correlator output stem");
    }

    auto &nl = envGet(DistillationNoise, par().leftNoise);
    auto &nr = envGet(DistillationNoise, par().rightNoise);
    auto *g = envGetGrid(FermionField);
    const unsigned int nt = env().getDim(g->Nd() - 1);
    dilSizeLS_ = {
        {Side::left,  static_cast<unsigned int>(nl.dilutionSize(Index::l) * nl.dilutionSize(Index::s))},
        {Side::right, static_cast<unsigned int>(nr.dilutionSize(Index::l) * nr.dilutionSize(Index::s))}
    };
    if (par().blockSize > dilSizeLS_.at(Side::left) ||
        par().blockSize > dilSizeLS_.at(Side::right))
    {
        HADRONS_ERROR(Size, "blockSize exceeds a LapH-spin dilution dimension");
    }

    isExact_ = envHasDerivedType(DistillationNoise, ExactDistillationPolicy<FImpl>, par().leftNoise)
            && envHasDerivedType(DistillationNoise, ExactDistillationPolicy<FImpl>, par().rightNoise);
    if (!isExact_ && par().noisePairs.empty())
    {
        HADRONS_ERROR(Argument, "noisePairs required for stochastic distillation");
    }

    dt_ = strToVec<unsigned int>(par().dt);

    if (dt_.empty())
    {
        HADRONS_ERROR(Argument,
            "Specify at least one dt");
    }

    for (auto t : dt_)
    {
        if (t >= nl.dilutionSize(Index::t) ||
            t >= nr.dilutionSize(Index::t))
        {
            HADRONS_ERROR(Range,
                "dt outside time dilution range");
        }
    }

    gammaSrc_ = strToVec<Gamma::Algebra>(par().gammaSrc);
    gammaSnk_ = strToVec<Gamma::Algebra>(par().gammaSnk);
    if (gammaSrc_.empty() || gammaSnk_.empty())
    {
        HADRONS_ERROR(Argument, "Both gammaSrc and gammaSnk must be specified");
    }

    momenta_.clear();
    if (par().momenta.empty())
    {
        for (int px=-2; px<=2; ++px)
        for (int py=-2; py<=2; ++py)
        for (int pz=-2; pz<=2; ++pz)
            if (px*px + py*py + pz*pz <= 4)
                momenta_.push_back({RealF(px), RealF(py), RealF(pz)});
    }
    else
    {
        for (const auto &s : par().momenta)
        {
            auto p = strToVec<RealF>(s);
            if (p.size() != static_cast<size_t>(g->Nd()-1))
            {
                HADRONS_ERROR(Size, "Momentum has incorrect number of components");
            }
            momenta_.push_back(p);
        }
    }    
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDistilMeson2pt<FImpl>::execute(void)
{
    if (!isExact_)
    {
        HADRONS_ERROR(Implementation,
            "DistilMeson2pt: stochastic dilution not yet implemented");
    }

    auto *g = envGetGrid(FermionField);
    auto *g3d = envGetSliceGrid(FermionField,g->Nd()-1);

    auto &epack = envGet(
        typename DistillationNoise::LapPack,
        par().lapEigenPack);

    auto &nl = envGet(DistillationNoise,par().leftNoise);
    auto &nr = envGet(DistillationNoise,par().rightNoise);

    if (nl.getNl()!=epack.evec.size() ||
        nr.getNl()!=epack.evec.size())
    {
        HADRONS_ERROR(Size,
            "Noise and eigenpack LapH dimensions disagree");
    }

    const uint nt=env().getDim(g->Nd()-1);

    Computation computation(
        {{Side::left,"phi"},{Side::right,"phi"}},
        g,g3d,nl,nr,
        par().blockSize,
        par().cacheSize,
        nt,
        momenta_.size(),
        gammaSnk_.size(),
        true,
        vm().getTrajectory(),
        par().leftVectorStem,
        par().rightVectorStem);

    // Fourier phases, following DistilMesonFieldFixed.
    std::vector<ComplexField> phases;
    phases.reserve(momenta_.size());

    ComplexField coor(g);
    const Complex imaginary(0.0,1.0);

    for (const auto &p : momenta_)
    {
        phases.emplace_back(g);
        phases.back()=Zero();

        for (uint mu=0;mu<p.size();++mu)
        {
            LatticeCoordinate(coor,mu);
            phases.back() =
                phases.back()+
                (p[mu]/env().getDim()[mu])*coor;
        }

        phases.back() =
            exp((Real)(2*M_PI)*imaginary*phases.back());
    }

    // Each side independently accepts a perambulator
    // or a generalized-perambulator vector stem.
    std::map<Side,PerambTensor&> peramb;

    if (!par().leftPeramb.empty())
    {
        auto &p=envGet(PerambTensor,par().leftPeramb);
        peramb.emplace(Side::left,p);
    }

    if (!par().rightPeramb.empty())
    {
        auto &p=envGet(PerambTensor,par().rightPeramb);
        peramb.emplace(Side::right,p);
    }

    for (auto s : sides)
    {
        if (peramb.count(s))
        {
            const auto &available =
                peramb.at(s).MetaData.timeSources;

            for (auto dt : dt_)
            {
                if (std::find(
                        available.begin(),
                        available.end(),
                        static_cast<int>(dt))
                    == available.end())
                {
                    HADRONS_ERROR(Argument,
                        "Requested dt absent from perambulator");
                }
            }
        }
    }

    std::vector<ComplexD> corr;

    computation.execute2pt(
        corr,
        dt_,
        gammaSrc_,
        gammaSnk_,
        phases,
        epack,
        {{Side::left,0},{Side::right,0}},
        peramb);

    const size_t nm=momenta_.size();
    const size_t ng0=gammaSrc_.size();
    const size_t ng1=gammaSnk_.size();

    std::vector<Result> results;
    results.reserve(dt_.size()*nm*ng0*ng1);

    for (size_t d=0;d<dt_.size();++d)
    for (size_t m=0;m<nm;++m)
    for (size_t gs=0;gs<ng0;++gs)
    for (size_t gt=0;gt<ng1;++gt)
    {
        Result result;

        result.dt=dt_[d];
        result.timePartition =
            nl.dilutionPartition(Index::t,dt_[d]);

        {
            std::ostringstream ss;
            ss << gammaSrc_[gs];
            result.gammaSrc=ss.str();
        }

        {
            std::ostringstream ss;
            ss << gammaSnk_[gt];
            result.gammaSnk=ss.str();
        }

        {
            std::ostringstream ss;
            for (size_t mu=0;mu<momenta_[m].size();++mu)
            {
                if (mu) ss << " ";
                ss << -momenta_[m][mu];
            }
            result.momSrc=ss.str();
        }

        {
            std::ostringstream ss;
            for (size_t mu=0;mu<momenta_[m].size();++mu)
            {
                if (mu) ss << " ";
                ss << momenta_[m][mu];
            }
            result.momSnk=ss.str();
        }

        result.corr.resize(nt);

        for (uint t=0;t<nt;++t)
        {
            result.corr[t] =
                corr[((((d*nm+m)*ng0+gs)*ng1+gt)*nt+t)];
        }

        results.push_back(std::move(result));
    }

    saveResult(
        par().output,
        "distilMeson2pt",
        results);    
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_DistilMeson2pt_hpp_

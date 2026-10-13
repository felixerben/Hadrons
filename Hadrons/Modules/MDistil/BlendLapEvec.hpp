/*
 * BlendLapEvec.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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

/*  END LEGAL */
#ifndef Hadrons_MDistil_BlendLapEvec_hpp_
#define Hadrons_MDistil_BlendLapEvec_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MDistil)

/******************************************************************************
 * Concatenate selected vectors from two LapPacks (without reorthogonalizing).
 ******************************************************************************/
struct BlendLapEvecPar: Serializable
{
    GRID_SERIALIZABLE_CLASS_MEMBERS(BlendLapEvecPar,
                                    std::string, lapEvecLow,
                                    std::string, lapEvecHigh,
                                    unsigned int, nVecLow,
                                    unsigned int, nVecHigh,
                                    std::string, fileName)
};

template <typename FImpl>
class TBlendLapEvec: public Module<BlendLapEvecPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    TBlendLapEvec(const std::string name);
    virtual ~TBlendLapEvec(void) {}
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual void setup(void);
    virtual void execute(void);
};

MODULE_REGISTER_TMP(BlendLapEvec, TBlendLapEvec<FIMPL>, MDistil);

template <typename FImpl>
TBlendLapEvec<FImpl>::TBlendLapEvec(const std::string name)
: Module<BlendLapEvecPar>(name) {}

template <typename FImpl>
std::vector<std::string> TBlendLapEvec<FImpl>::getInput(void)
{
    std::vector<std::string> in{par().lapEvecLow, par().lapEvecHigh};
    return in;
}

template <typename FImpl>
std::vector<std::string> TBlendLapEvec<FImpl>::getOutput(void)
{
    std::vector<std::string> out{getName()};
    return out;
}

template <typename FImpl>
void TBlendLapEvec<FImpl>::setup(void)
{
    auto *gridHD = envGetGrid(FermionField);
    envCreate(typename DistillationNoise<FImpl>::LapPack,
              getName(), 1, par().nVecLow + par().nVecHigh, gridHD);
}

template <typename FImpl>
void TBlendLapEvec<FImpl>::execute(void)
{
    auto &low = envGet(typename DistillationNoise<FImpl>::LapPack, par().lapEvecLow);
    auto &high = envGet(typename DistillationNoise<FImpl>::LapPack, par().lapEvecHigh);
    auto &out = envGet(typename DistillationNoise<FImpl>::LapPack, getName());

    if (par().nVecLow > low.evec.size())
    {
        HADRONS_ERROR(Range, "nVecLow exceeds first input eigenpack size");
    }
    if (par().nVecHigh > high.evec.size())
    {
        HADRONS_ERROR(Range, "nVecHigh exceeds second input eigenpack size");
    }

    for (unsigned int i = 0; i < par().nVecLow; ++i)
    {
        out.evec[i] = low.evec[i];
        out.eval[i] = low.eval[i];
    }
    for (unsigned int i = 0; i < par().nVecHigh; ++i)
    {
        const unsigned int dest = par().nVecLow + i;
        out.evec[dest] = high.evec[i];
        out.eval[dest] = high.eval[i];
    }

    LOG(Message) << "Concatenated " << par().nVecLow << " vectors from '"
                 << par().lapEvecLow << "' and " << par().nVecHigh
                 << " vectors from '" << par().lapEvecHigh << "'" << std::endl;

    if (!par().fileName.empty())
    {
        out.record.solverXml = parString();
        out.record.operatorXml = "<module><id><type>MDistil::BlendLapEvec</type></id></module>";
        out.write(par().fileName + "." + std::to_string(vm().getTrajectory()), false);
    }
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_BlendLapEvec_hpp_

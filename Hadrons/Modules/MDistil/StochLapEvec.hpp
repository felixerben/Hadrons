/*
 * StochLapEvec.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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


#ifndef Hadrons_MDistil_StochLapEvec_hpp_
#define Hadrons_MDistil_StochLapEvec_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DilutedNoise.hpp>
#include <cmath>
#include <algorithm>

BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MDistil)

/******************************************************************************
 
Stochastic high-modes based on Laplacian eigenvectors

Implementation of the idea of blending

https://arxiv.org/pdf/2505.01719
 
 ******************************************************************************/

struct StochLapEvecPar: Serializable
{
    GRID_SERIALIZABLE_CLASS_MEMBERS(StochLapEvecPar,
                                    std::string, lapEigenPack,
                                    unsigned int, nVec,
                                    unsigned int, nStoch,
                                    std::string, fileName,
                                    bool, checkOrthogonality)
};

/******************************************************************************
 
 Stochastic high-modes based on Laplacian eigenvectors - Module (class) definition
 
 ******************************************************************************/

template <typename FImpl>
class TStochLapEvec: public Module<StochLapEvecPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    TStochLapEvec(const std::string name);
    virtual ~TStochLapEvec(void) {}
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual void setup(void);
    virtual void execute(void);
};

MODULE_REGISTER_TMP(StochLapEvec, TStochLapEvec<FIMPL>, MDistil);

/******************************************************************************
 *                 TStochLapEvec implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TStochLapEvec<FImpl>::TStochLapEvec(const std::string name) : Module<StochLapEvecPar>(name) {}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TStochLapEvec<FImpl>::getInput(void)
{
    std::vector<std::string> in={par().lapEigenPack};
    return in;
}


template <typename FImpl>
std::vector<std::string> TStochLapEvec<FImpl>::getOutput(void)
{
    std::vector<std::string> out{ getName() };
    return out;
}

template <typename FImpl>
void TStochLapEvec<FImpl>::setup(void)
{
    LOG(Message) << "StochLapEvec parameters: "
             << "lapEigenPack=" << par().lapEigenPack
             << ", nVec=" << par().nVec
             << ", nStoch=" << par().nStoch
             << ", fileName=" << par().fileName
             << ", checkOrthogonality=" << par().checkOrthogonality
             << std::endl;

    if (par().nStoch == 0)
    {
        HADRONS_ERROR(Range, "nStoch must be positive");
    }
    auto *gridHD = envGetGrid(FermionField);
    auto *gridLD = envGetSliceGrid(FermionField, gridHD->Nd() - 1);
    envTmp(ColourVectorField, "noise3d", 1, gridLD);
    envTmp(ColourVectorField, "work3d", 1, gridLD);
    envTmp(ColourVectorField, "low3d", 1, gridLD);
    envTmp(ColourVectorField, "basis3d", 1, gridLD);
    envTmp(ColourVectorField, "noise4d", 1, gridHD);
    envCreate(typename DistillationNoise<FImpl>::LapPack,
              getName(), 1, par().nStoch, gridHD);
}

template <typename FImpl>
void TStochLapEvec<FImpl>::execute(void)
{
    auto &low = envGet(typename DistillationNoise<FImpl>::LapPack, par().lapEigenPack);
    auto &out = envGet(typename DistillationNoise<FImpl>::LapPack, getName());
    envGetTmp(ColourVectorField, noise3d);
    envGetTmp(ColourVectorField, work3d);
    envGetTmp(ColourVectorField, low3d);
    envGetTmp(ColourVectorField, basis3d);
    envGetTmp(ColourVectorField, noise4d);

    if (par().nVec > low.evec.size())
    {
        HADRONS_ERROR(Range, "nVec exceeds the input eigenpack size");
    }

    auto *gridHD = envGetGrid(FermionField);
    const int ntLocal = gridHD->LocalDimensions()[Tdir];
    const int ntFirst = gridHD->LocalStarts()[Tdir];

    // Each four-dimensional Gaussian draw is deterministic under Hadrons rng4d().
    // The slice projection/orthogonalization is independent for each time slice.
    for (unsigned int j = 0; j < par().nStoch; ++j)
    {
        LOG(Message) << "Vec " << j << " / " << par().nStoch << std::endl;
        gaussian(this->rng4d(), noise4d);
        for (int t = 0; t < ntLocal; ++t)
        {
            //LOG(Message) << "slice " << t << " / " << ntLocal << std::endl;
            ExtractSliceLocal(noise3d, noise4d, 0, t, Tdir);
            work3d = noise3d;

            // Twice-reorthogonalized modified Gram-Schmidt against the low modes
            // and the already constructed stochastic vectors.
            for (int pass = 0; pass < 2; ++pass)
            {
                for (unsigned int i = 0; i < par().nVec; ++i)
                {
                    ExtractSliceLocal(low3d, low.evec[i], 0, t, Tdir);
                    const auto coeff = innerProduct(low3d, work3d);
                    work3d -= coeff * low3d;
                }
                for (unsigned int i = 0; i < j; ++i)
                {
                    ExtractSliceLocal(basis3d, out.evec[i], 0, t, Tdir);
                    const auto coeff = innerProduct(basis3d, work3d);
                    work3d -= coeff * basis3d;
                }
            }
            const RealD n2 = norm2(work3d);
            if (!(n2 > 0.0) || !std::isfinite(n2))
            {
                HADRONS_ERROR(Program, "Stochastic complement orthogonalization failed");
            }
            work3d *= (1.0 / std::sqrt(n2));
            InsertSliceLocal(work3d, out.evec[j], 0, t, Tdir);
        }
        // Stochastic vectors have no Laplacian eigenvalues.
        out.eval[j] = 0.0;
    }

    // Optional diagnostic: check the completed basis on every local time slice.
    // innerProduct() performs the spatial reduction over the slice grid.
    if (par().checkOrthogonality)
    {
        RealD maxLowOverlap = 0.0;
        RealD maxStochResidual = 0.0;
    
        for (int t = 0; t < ntLocal; ++t)
        {
            LOG(Message) << "checking orthogonality: slice " << t << " / " << ntLocal << std::endl;
            for (unsigned int j = 0; j < par().nStoch; ++j)
            {
                ExtractSliceLocal(work3d, out.evec[j], 0, t, Tdir);
    
                // Check orthogonality against the low Laplacian modes
                for (unsigned int i = 0; i < par().nVec; ++i)
                {
                    ExtractSliceLocal(low3d, low.evec[i], 0, t, Tdir);
    
                    const auto overlap = innerProduct(low3d, work3d);
                    const RealD overlapAbs = std::sqrt(
                        overlap.real() * overlap.real() +
                        overlap.imag() * overlap.imag());
    
                    if (overlapAbs > maxLowOverlap)
                    {
                        maxLowOverlap = overlapAbs;
                    }
                }
    
                // Check mutual orthonormality of stochastic vectors
                for (unsigned int i = 0; i <= j; ++i)
                {
                    ExtractSliceLocal(basis3d, out.evec[i], 0, t, Tdir);
    
                    const auto overlap = innerProduct(basis3d, work3d);
                    const auto residual =
                        (i == j) ? overlap - ComplexD(1.0, 0.0) : overlap;
    
                    const RealD residualAbs = std::sqrt(
                        residual.real() * residual.real() +
                        residual.imag() * residual.imag());
    
                    if (residualAbs > maxStochResidual)
                    {
                        maxStochResidual = residualAbs;
                    }
                }
            }
        }
    
        // Global maximum over all MPI ranks and time slices
        auto *grid = envGetGrid(FermionField);
        grid->GlobalMax(maxLowOverlap);
        grid->GlobalMax(maxStochResidual);
    
        LOG(Message) << "StochLapEvec orthogonality check: "
                     << "max |<low|stoch>| = " << maxLowOverlap
                     << ", max |<stoch_i|stoch_j> - delta_ij| = "
                     << maxStochResidual
                     << std::endl;
    }


    if (!par().fileName.empty())
    {
        out.record.solverXml = parString();
        out.record.operatorXml = "<module><id><type>MDistil::StochLapEvec</type></id></module>";
        out.write(par().fileName + "." + std::to_string(vm().getTrajectory()), false);
    }
    LOG(Message) << "Generated " << par().nStoch << " stochastic complement vectors"
                 << " orthogonal to " << par().nVec << " low modes" << std::endl;
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_StochLapEvec_hpp_

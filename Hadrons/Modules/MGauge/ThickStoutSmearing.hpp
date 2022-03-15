/*
 * ThickStoutSmearing.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Felix Erben      <felix.erben@ed.ac.uk>
 * Author: Antonin Portelli <antonin.portelli@me.com>
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
#ifndef Hadrons_MGauge_ThickStoutSmearing_hpp_
#define Hadrons_MGauge_ThickStoutSmearing_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                            Stout smearing                                  *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGauge)

class ThickStoutSmearingPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ThickStoutSmearingPar,
                                    std::string, gauge,
                                    unsigned int, steps,
                                    std::string, orthogDim,
                                    double, rho);
};

template <typename GImpl>
class TThickStoutSmearing: public Module<ThickStoutSmearingPar>
{
public:
    GAUGE_TYPE_ALIASES(GImpl,);
public:
    // constructor
    TThickStoutSmearing(const std::string name);
    // destructor
    virtual ~TThickStoutSmearing(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(ThickStoutSmearing, TThickStoutSmearing<GIMPL>, MGauge);

/******************************************************************************
 *                     TThickStoutSmearing implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl>
TThickStoutSmearing<GImpl>::TThickStoutSmearing(const std::string name)
: Module<ThickStoutSmearingPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl>
std::vector<std::string> TThickStoutSmearing<GImpl>::getInput(void)
{
    std::vector<std::string> in = {par().gauge};
    
    return in;
}

template <typename GImpl>
std::vector<std::string> TThickStoutSmearing<GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl>
void TThickStoutSmearing<GImpl>::setup(void)
{
    envCreateLat(GaugeField, getName());
    envTmpLat(GaugeField, "buf");
    envTmpLat(GaugeField, "Utmp");
    envTmpLat(GaugeField, "UsmrTmp");
    envTmp(Lattice<iScalar<vInteger>>, "tlat",1, envGetGrid(LatticeComplex));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl>
void TThickStoutSmearing<GImpl>::execute(void)
{
    LOG(Message) << "Smearing '" << par().gauge << "' with " << par().steps
                 << " step" << ((par().steps > 1) ? "s" : "") 
                 << " of stout smearing and rho= " << par().rho << std::endl;

    Smear_Stout<GImpl> smearer(par().rho, -1);
    auto               &U    = envGet(GaugeField, par().gauge);
    auto               &Usmr = envGet(GaugeField, getName());

    Usmr = U;
    LOG(Message) << "plaquette= " << WilsonLoops<GImpl>::avgPlaquette(U)
                 << std::endl;
    const int  Nt{env().getDim(Tdir)};
    envGetTmp(GaugeField, Utmp);
    envGetTmp(GaugeField, UsmrTmp);
    //should be input
    int t_thick=1;
    envGetTmp(Lattice<iScalar<vInteger>>, tlat);
    LatticeCoordinate(tlat, Tp);
    for (int t = 0; t < Nt; t++)
    {
        //this only works if t is sufficiently far away from the boundary
        //Utmp=where((t-t_thick <= tlat && tlat <= t+t_thick),U,0.*U);
        Utmp=where((t==tlat),U,0.*U);
        //Utmp=where((t==tlat),U,U);
        LOG(Message) << "plaquette= " << WilsonLoops<GImpl>::avgPlaquette(Utmp)
                 << std::endl;
        envGetTmp(GaugeField, buf);
        buf = Utmp;
        // smear the thick timeslice of the field using standard stout smearing
        for (unsigned int n = 0; n < par().steps; ++n)
        {
            smearer.smear(UsmrTmp, buf);
            buf = UsmrTmp;
    LOG(Message) << "plaquette= " << n << " " << WilsonLoops<GImpl>::avgPlaquette(UsmrTmp)
                 << std::endl;
        }
        // store the smeared field in the output on timeslice t
        Usmr=where((t == tlat),UsmrTmp,Usmr);
    }
    LOG(Message) << "plaquette= " << WilsonLoops<GImpl>::avgPlaquette(Usmr)
                 << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGauge_ThickStoutSmearing_hpp_

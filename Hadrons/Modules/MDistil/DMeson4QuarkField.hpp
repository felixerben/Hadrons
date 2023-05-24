#ifndef Hadrons_MDistil_DMeson4QuarkField_hpp_
#define Hadrons_MDistil_DMeson4QuarkField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/********************************************************************************
 *                         DMeson4QuarkField                                    *
 * Computes the following sub-diagram:                                          *
 *                                                                              *
 *                                                                              *
 *           ____            ___                                                *
 *          /    \          /                                                   *
 *         /      \        /                                                    *
 *        /         v1  v3                                                      *
 *  M(rho1,rho2)                                                                *
 *        \         v2  v4                                                      *
 *         \      /        \                                                    *
 *          \____/          \___                                                *
 *                                                                              *
 *                                                                              *
 *  D(t=0)           H_W(t=deltaT)                                              *
 *                                                                              *
 *                                                                              *
 *******************************************************************************/
BEGIN_MODULE_NAMESPACE(MDistil)

class DMeson4QuarkFieldPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DMeson4QuarkFieldPar,
                                    std::string,                outPath,       // file stem for the out file
                                    std::string,                DMesonField,   // M(rho,rho) meson field for the D
                                    std::string,                vectorStem1,   // charm
                                    std::string,                vectorStem2,   // SU(3) light
                                    std::string,                vectorStem3,   // SU(3) light
                                    std::string,                vectorStem4,   // SU(3) light
                                    std::vector<int>,           Dcontractions, // if this is e.g. {1, 3} then contract rho1 with v1 and rho2 with v3 
                                    std::string,                noisePol1,     // noise policy of v1 - assert compatibility with M(rho,rho) 
                                    std::string,                noisePol2,     // noise policy of v2 - assert compatibility with M(rho,rho) 
                                    std::string,                noisePol3,     // noise policy of v3
                                    std::string,                noisePol4,     // noise policy of v4
                                    std::string,                timeSources1,  // time sources used for v1 - assert compatibility with M(rho,rho) 
                                    std::string,                timeSources2,  // time sources used for v2 - assert compatibility with M(rho,rho) 
                                    std::string,                timeSources3,  // time sources used for v3
                                    std::string,                timeSources4,  // time sources used for v4
                                    std::vector<std::string>,   noisePairs,    // specifying which noise pairs to contract over - assert that they are in M(rho,rho) and compatible with v1,v2
                                    std::string,                deltaT,        // time between D and Hw
                                    std::string,                gamma,         // in Hw, we could hard-code this
                                    unsigned int,               blockSize,     // tunable parameters
                                    unsigned int,               cacheSize);
};

template <typename FImpl>
class TDMeson4QuarkField: public Module<DMeson4QuarkFieldPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    // constructor
    TDMeson4QuarkField(const std::string name);
    // destructor
    virtual ~TDMeson4QuarkField(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(DMeson4QuarkField, TDMeson4QuarkField<FIMPL>, MDistil);

/******************************************************************************
 *                 TDMeson4QuarkField implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TDMeson4QuarkField<FImpl>::TDMeson4QuarkField(const std::string name)
: Module<DMeson4QuarkFieldPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TDMeson4QuarkField<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().DMesonField, par().vectorStem1, par().vectorStem2, par().vectorStem3, par().vectorStem4};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TDMeson4QuarkField<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDMeson4QuarkField<FImpl>::setup(void)
{

    envTmpLat(FermionField, "fermion4dtmp1");
    envTmpLat(FermionField, "fermion4dtmp2");
    envTmpLat(FermionField, "fermion4dtmp3");
    envTmpLat(FermionField, "fermion4dtmp4");
    
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDMeson4QuarkField<FImpl>::execute(void)
{
    
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_DMeson4QuarkField_hpp_

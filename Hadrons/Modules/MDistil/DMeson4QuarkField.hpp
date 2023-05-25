#ifndef Hadrons_MDistil_DMeson4QuarkField_hpp_
#define Hadrons_MDistil_DMeson4QuarkField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MDistil/Base.hpp>

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
                                    unsigned int,               deltaT,            // time separation tH - tW
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
    std::vector<std::string> in = {par().vectorStem1, par().vectorStem2, par().vectorStem3, par().vectorStem4,par().noisePol1};
    
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
    GridCartesian * gridHD = envGetGrid(FermionField);
    GridCartesian * gridLD = envGetSliceGrid(FermionField,gridHD->Nd() -1);
    
    envTmp(FermionField, "fermion3dtmp1",1,gridLD);
    envTmp(FermionField, "fermion3dtmp2",1,gridLD);
    envTmp(FermionField, "fermion3dtmp3",1,gridLD);
    envTmp(PropagatorField, "prop3dtmp",1,gridLD);
    envTmp(ComplexField, "MPhiPhi",1,gridLD);

    // the 4 Distillation vectors in the 4-quark field
    envTmpLat(FermionField, "fermion4dtmp1");
    envTmpLat(FermionField, "fermion4dtmp2");
    envTmpLat(PropagatorField, "prop4dtmp");
    envTmpLat(FermionField, "fermion4dtmp3");
    envTmpLat(FermionField, "fermion4dtmp4");
    
    envTmp(ComplexField, "cplx3dtmp",1,gridLD);

    
    //envTmp(DistilMesonFieldMatrix<ComplexF>, "PerambDT",1,Nt,nVec,nDL_reduced,nNoise,nDS);
    
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDMeson4QuarkField<FImpl>::execute(void)
{

    std::map<std::string, DistilMesonFieldMatrix<ComplexD>> mf_out; 
    
    GridCartesian * gridHD = envGetGrid(FermionField);
    GridCartesian * gridLD = envGetSliceGrid(FermionField,gridHD->Nd() -1);
    const int Ntlocal{gridHD->LocalDimensions()[Tdir]};
    const int Ntfirst{gridHD->LocalStarts()[Tdir]};

    std::string mfPath = par().DMesonField;   
    LOG(Message) << "reading " << mfPath << std::endl;
    TimerArray timer;
    ContractionDistilMesonField<ComplexD,ComplexF> DMeson(mfPath,env().getDim(Tdir), timer);
    
    auto &dilNoise = envGet(DistillationNoise<FImpl>, par().noisePol1);
    int nNoise = dilNoise.size(); 
    if(nNoise>1)
    {
        HADRONS_ERROR(Range, "DMeson4QuarkField only implemented for exact distillation");
    }
    const int iNoise=0;
    int nDL = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::l);        
    int nDS = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::s);        
    int nDT = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::t);        
    int nD = nDL * nDS * nDT;
    
    int nT=env().getDim(Tdir);
    
    int deltaT = par().deltaT; 
    if(deltaT>=nT)
    {
        HADRONS_ERROR(Range, "deltaT must be smaller than nT");
    }
    
    envGetTmp(FermionField, fermion4dtmp1);
    envGetTmp(FermionField, fermion3dtmp1);
    envGetTmp(FermionField, fermion3dtmp2);
    envGetTmp(FermionField, fermion3dtmp3);
    envGetTmp(FermionField, fermion4dtmp2);
    envGetTmp(FermionField, fermion4dtmp3);
    envGetTmp(FermionField, fermion4dtmp4);
    envGetTmp(PropagatorField, prop4dtmp);
    envGetTmp(PropagatorField, prop3dtmp);
    envGetTmp(ComplexField, MPhiPhi);
    envGetTmp(ComplexField, cplx3dtmp);
    
    
    Gamma                  gX(Gamma::Algebra::GammaX);
    
    
    LOG(Message) << "WARNING: Assuming ordering s + ns*(l + nl*t) in DilutedNoise.hpp. This code will break when this changes!" << std::endl;
    int T1,T2,dk1,ds1,dk2,ds2,dSolve1,dSolve2,tD,tH;
    std::array<unsigned int, 3> index1,index2;
    for (int t = 0; t < Ntlocal; t++ )
    {
        tH = t + Ntfirst;
        tD = (tH - deltaT + nT)%nT;
        MPhiPhi=Zero();      
        for(int id1=0; id1<nDL * nDS; id1++)
        {
            // this line is where the ordering s + ns*(l + nl*t) is assumed
            index1 = dilNoise.dilutionCoordinates(id1);  
            dk1 = index1[DistillationNoise<FImpl>::Index::l];
            ds1 = index1[DistillationNoise<FImpl>::Index::s];
            auto &solve1 = envGet(std::vector<FermionField>, par().vectorStem1);
            dSolve1 = ds1 + nDS * dk1;
            fermion4dtmp1 = solve1[dSolve1];
            // this is vector 1 on timeslice tH 
            ExtractSliceLocal(fermion3dtmp1,fermion4dtmp1,0,t,Tdir);
            for(int id2=0; id2<nDL * nDS; id2++)
            {
                index2 = dilNoise.dilutionCoordinates(id2);  
                dk2 = index2[DistillationNoise<FImpl>::Index::l];
                ds2 = index2[DistillationNoise<FImpl>::Index::s];
                auto &solve2 = envGet(std::vector<FermionField>, par().vectorStem2);
                dSolve2 = ds2 + nDS * dk2;
                fermion4dtmp2 = solve2[dSolve2];
                ExtractSliceLocal(fermion3dtmp2,fermion4dtmp2,0,t,Tdir);
                fermion3dtmp3 = gX*fermion3dtmp2;
                LOG(Message) << "Meson Field " << tD << ", " << tD << " at t= " << tD << ": " << id1 << " " << id2 << " is " <<  DMeson(tD,tD,tD)(id1,id2) << std::endl;
                fermion3dtmp2 = DMeson(tD,tD,tD)(id1,id2)*fermion3dtmp3;
                prop3dtmp = outerProduct(fermion3dtmp1,fermion3dtmp2);
                //prop4dtmp = outerProduct(fermion4dtmp1,fermion4dtmp2);
                MPhiPhi += trace(prop3dtmp);
            }
        }
        for(int id1=0; id1<nDL * nDS; id1++)
        {
            // this line is where the ordering s + ns*(l + nl*t) is assumed
            index1 = dilNoise.dilutionCoordinates(id1);  
            dk1 = index1[DistillationNoise<FImpl>::Index::l];
            ds1 = index1[DistillationNoise<FImpl>::Index::s];
            auto &solve1 = envGet(std::vector<FermionField>, par().vectorStem3);
            dSolve1 = ds1 + nDS * dk1;
            fermion4dtmp1 = solve1[dSolve1];
            ExtractSliceLocal(fermion3dtmp1,fermion4dtmp1,0,t,Tdir);
            for(int id2=0; id2<nDL * nDS; id2++)
            {
                index2 = dilNoise.dilutionCoordinates(id2);  
                dk2 = index2[DistillationNoise<FImpl>::Index::l];
                ds2 = index2[DistillationNoise<FImpl>::Index::s];
                auto &solve2 = envGet(std::vector<FermionField>, par().vectorStem4);
                dSolve2 = ds2 + nDS * dk2;
                fermion4dtmp2 = solve2[dSolve2];
                ExtractSliceLocal(fermion3dtmp2,fermion4dtmp2,0,t,Tdir);
                fermion3dtmp3 = gX*fermion3dtmp2;
                LOG(Message) << "Meson Field " << tD << ", " << tD << " at t= " << tD << ": " << id1 << " " << id2 << " is " <<  DMeson(tD,tD,tD)(id1,id2) << std::endl;
                fermion3dtmp2 = fermion3dtmp3;
                prop3dtmp = outerProduct(fermion3dtmp1,fermion3dtmp2);
                cplx3dtmp = trace(prop3dtmp)*MPhiPhi;
                // MAYBE MULTIPLY BY PHASE HERE!!
            }
        }
    }
    
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_DMeson4QuarkField_hpp_

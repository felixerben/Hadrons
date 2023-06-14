#ifndef Hadrons_MDistil_DMeson4QuarkField_hpp_
#define Hadrons_MDistil_DMeson4QuarkField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/TimerArray.hpp>
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
 *  D(t=tD)           H_W(t)       tKpi                                         *
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
                                    unsigned int,               tD,            // time of D meson
                                    unsigned int,               tKpi,          // time of yet uncontracted part
                                    Gamma::Algebra,             gamma12,         // between vector1 and vector2
                                    Gamma::Algebra,             gamma34,         // between vector3 and vector4
                                    std::string,                mom,           // momentum injected into Hw
                                    unsigned int,               blockSize,     // tunable parameters
                                    unsigned int,               cacheSize);
};

template <typename FImpl>
class TDMeson4QuarkField: public Module<DMeson4QuarkFieldPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        Gamma::Algebra, gamma_snk,
                                        Gamma::Algebra, gamma_src,
                                        std::vector<Complex>, corr);
    };
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
    
    envTmpLat(FermionField,    "fermion4dtmp");
    envTmp   (FermionField,    "fermion3dtmp1" ,1, gridLD);
    envTmp   (FermionField,    "fermion3dtmp2" ,1, gridLD);
    envTmp   (FermionField,    "fermion3dtmp3" ,1, gridLD);
    envTmp   (PropagatorField, "prop3dtmp"     ,1, gridLD);
    envTmp   (ComplexField,    "MPhiPhi"       ,1, gridLD);
    envTmp   (ComplexField,    "cplx3dtmp"     ,1, gridLD);
    envTmpLat(ComplexField,    "ph");
    envTmp   (ComplexField,    "ph3d"          ,1, gridLD);
    envTmpLat(ComplexField,    "coor");

    auto &dilNoise = envGet(DistillationNoise<FImpl>, par().noisePol1);
    int nDL = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::l);        
    int nDS = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::s);     
    envTmp(Vector<HADRONS_DISTIL_IO_TYPE>, "block_buf", 1, nDL * nDS * nDL * nDS);
    envTmp(Vector<HADRONS_DISTIL_TYPE>,    "cache_buf", 1, nDL * nDS * nDL * nDS);
    
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDMeson4QuarkField<FImpl>::execute(void)
{
    // general grid setup
    GridCartesian * gridHD = envGetGrid(FermionField);
    GridCartesian * gridLD = envGetSliceGrid(FermionField,gridHD->Nd() -1);
    const int Ntlocal{gridHD->LocalDimensions()[Tdir]};
    const int Ntfirst{gridHD->LocalStarts()[Tdir]};
    int nT=env().getDim(Tdir);
    
    // block and cache to store the output in
    envGetTmp(Vector<HADRONS_DISTIL_IO_TYPE>, block_buf);
    envGetTmp(Vector<HADRONS_DISTIL_TYPE>, cache_buf);
    
    // read input D-meson field
    std::string mfPath = par().DMesonField;   
    LOG(Message) << "reading " << mfPath << std::endl;
    TimerArray timer;
    ContractionDistilMesonField<ComplexD,ComplexF> DMeson(mfPath,env().getDim(Tdir), timer);
     
    // noise class -- assert they are identical and an "exact distillation" policy
    auto &dilNoise = envGet(DistillationNoise<FImpl>, par().noisePol1);
    auto &dN2 = envGet(DistillationNoise<FImpl>, par().noisePol2);
    auto &dN3 = envGet(DistillationNoise<FImpl>, par().noisePol3);
    auto &dN4 = envGet(DistillationNoise<FImpl>, par().noisePol4);
    if(dN2.generateHash() != dilNoise.generateHash() || dN3.generateHash() != dilNoise.generateHash()|| dN4.generateHash() != dilNoise.generateHash())
    {
        HADRONS_ERROR(Implementation, "All noise policies must be identical");
    }
    int nNoise = dilNoise.size(); 
    if(nNoise>1)
    {
        HADRONS_ERROR(Implementation, "DMeson4QuarkField only implemented for exact distillation");
    }
    int nDL = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::l);        
    int nDS = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::s);        
    int nDT = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::t);        
    // other input parameters
    int tD = par().tD; 
    int tKpi = par().tKpi; 
    if(tD>=nT)
    {
        HADRONS_ERROR(Range, "tD must be smaller than nT");
    }
    if(tKpi>=nT)
    {
        HADRONS_ERROR(Range, "tKpi must be smaller than nT");
    }
    Gamma                  g12(par().gamma12);
    Gamma                  g34(par().gamma34);
    // momentum phase e^{ipx} for Hw
    Complex           i(0.0,1.0);
    std::vector<Real> p;
    p  = strToVec<Real>(par().mom);
    envGetTmp(ComplexField, coor);
    envGetTmp(ComplexField, ph);
    envGetTmp(ComplexField, ph3d);
    ph = Zero();
    for(unsigned int mu = 0; mu < env().getNd(); mu++)
    {
        LatticeCoordinate(coor, mu);
        ph = ph + (p[mu]/env().getDim(mu))*coor;
    }
    ph = exp((Real)(2*M_PI)*i*ph);

    // file name of output
    std::string outPath = par().outPath; 
    std::stringstream ss;
    ss << par().gamma12 << "_" << par().gamma34 << "_p";
    for (unsigned int mu = 0; mu < p.size(); ++mu)
            ss << p[mu] << ((mu == p.size() - 1) ? "" : "_");
    ss << ".h5";   
    outPath += "/D-Hw.tD" + std::to_string(tD) +".tKpi"+ std::to_string(tKpi) + "/" + ss.str();
    
    // Temporary objects
    envGetTmp(FermionField,    fermion4dtmp);
    envGetTmp(FermionField,    fermion3dtmp1);
    envGetTmp(FermionField,    fermion3dtmp2);
    envGetTmp(FermionField,    fermion3dtmp3);
    envGetTmp(PropagatorField, prop3dtmp);
    envGetTmp(ComplexField,    MPhiPhi);
    envGetTmp(ComplexField,    cplx3dtmp);
    
    // initialise file and metadata
    DistilMesonFieldMetadata<FImpl> md;
    for (auto pmu: p)
    {
        md.Momentum.push_back(pmu);
     }
    std::stringstream ss2;
    ss2 << par().gamma12 << "_" << par().gamma34;
    md.Operator          = ss2.str();
    md.Nt                = nT;   
    md.Nvec              = nDL;     //nvec=nDL for exact
    md.NoisePair         = {0,0};
    md.MesonFieldType    = "MD-Hw";
    md.RelativeSide      = "none";
    md.NoiseHashLeft     = "0";
    md.NoiseHashRight    = "0";
    //md.TimeDilutionLeft  = dilNoise.getMap();//[Index::t];
    //md.TimeDilutionRight = dilNoise.getMap();//[Index::t];
    //md.LapDilutionLeft   = index1[DistillationNoise<FImpl>::Index::l];
    //md.LapDilutionRight  = index1[DistillationNoise<FImpl>::Index::l];
    //md.SpinDilutionLeft  = index1[DistillationNoise<FImpl>::Index::s];
    //md.SpinDilutionRight = index1[DistillationNoise<FImpl>::Index::s];
   

    startTimer("file creation");
    makeFileDir(outPath, gridHD);
    unsigned int myRank = gridHD->ThisRank(); 
    DistilMatrixIo<HADRONS_DISTIL_IO_TYPE> matrix_io(outPath, DISTIL_MATRIX_NAME, nT, nDL * nDS, nDL * nDS);
    if(myRank==0)
    {
        matrix_io.initFile(md);
    }
    gridHD->Barrier();
    stopTimer("file creation");
    
    LOG(Message) << "WARNING: Assuming ordering s + ns*(l + nl*t) in DilutedNoise.hpp. This code will break when this changes!" << std::endl;
    // variables used in the loop structure
    int dk1,ds1,dk2,ds2,dSolve1,dSolve2,tH;
    std::array<unsigned int, 3> index1,index2;
    std::vector<TComplex>  buf;
    
    const uint i_rank =  gridHD->ThisRank();
    const uint N_ranks = gridHD->RankCount();   
 
    // loop over tH
    for (int t = 0; t < Ntlocal; t++ )
    {
        tH = t + Ntfirst;
        MPhiPhi=Zero();    
        // 3D phase e^{ipx}
        ExtractSliceLocal(ph3d,ph,0,t,Tdir);  
        startTimer("computation MPhiPhi");
        for(int id1=0; id1<nDL * nDS; id1++)
        {
            // this line is where the ordering s + ns*(l + nl*t) is assumed
            index1 = dilNoise.dilutionCoordinates(id1);  
            dk1 = index1[DistillationNoise<FImpl>::Index::l];
            ds1 = index1[DistillationNoise<FImpl>::Index::s];
            auto &solve1 = envGet(std::vector<FermionField>, par().vectorStem1);
            // full dilution assumed here, and that solve 1 comes from timeslice tD = d_{tD}
            dSolve1 = ds1 + nDS * dk1 + nDL * nDS * tD;
            fermion4dtmp = solve1[dSolve1];
            // this is vector 1 on timeslice tH 
            ExtractSliceLocal(fermion3dtmp1,fermion4dtmp,0,t,Tdir);
            for(int id2=0; id2<nDL * nDS; id2++)
            {
                index2 = dilNoise.dilutionCoordinates(id2);  
                dk2 = index2[DistillationNoise<FImpl>::Index::l];
                ds2 = index2[DistillationNoise<FImpl>::Index::s];
                auto &solve2 = envGet(std::vector<FermionField>, par().vectorStem2);
                // full dilution assumed here, and that solve 2 comes from timeslice tD = d_{tD}
                dSolve2 = ds2 + nDS * dk2 + nDL * nDS * tD; 
                fermion4dtmp = solve2[dSolve2];
                // this is vector 2 on timeslice tH 
                ExtractSliceLocal(fermion3dtmp2,fermion4dtmp,0,t,Tdir);
                fermion3dtmp3 = g12*fermion3dtmp2;
                fermion3dtmp2 = DMeson(tD,tD,tD)(id1,id2)*fermion3dtmp3;
                prop3dtmp = outerProduct(fermion3dtmp1,fermion3dtmp2);
                // this object is sum_{spin,colour,d1,d2} (DMeson[d1,d2] * vector1[d1] * gamma12 * vector2[d2]) on timeslice tH
                MPhiPhi += trace(prop3dtmp);
            }
        }
        stopTimer("computation MPhiPhi");
        /*************************************************
        FE: checked here that a sliceSum over MPhiPhi
        reproduces exactly a contraction of meson fields
        tr[ M(rho,rho; tD,tD,tD) * M(phi,phi; tD,tD,tD) ]
        *************************************************/
        startTimer("computation D-4quark");
        DistilMatrixSetIo<ComplexF> block(block_buf.data(), 1 , 1, nDL * nDS, nDL * nDS);
        for(int id1=0; id1<nDL * nDS; id1++)
        {
            // this line is where the ordering s + ns*(l + nl*t) is assumed
            index1 = dilNoise.dilutionCoordinates(id1);  
            dk1 = index1[DistillationNoise<FImpl>::Index::l];
            ds1 = index1[DistillationNoise<FImpl>::Index::s];
            auto &solve1 = envGet(std::vector<FermionField>, par().vectorStem3);
            // full dilution assumed here, and that solve 3 comes from timeslice tKpi = d_{tKpi}
            dSolve1 = ds1 + nDS * dk1 + nDL * nDS * tKpi; 
            fermion4dtmp = solve1[dSolve1];
            ExtractSliceLocal(fermion3dtmp1,fermion4dtmp,0,t,Tdir);
            for(int id2=0; id2<nDL * nDS; id2++)
            {
                // no caching for the moment - but keep this here in case anyone wants to optimise this code at some stage
                DistilMatrixSetCache<ComplexD> cache(cache_buf.data(), 1, 1, 1, 1, 1);
                index2 = dilNoise.dilutionCoordinates(id2);  
                dk2 = index2[DistillationNoise<FImpl>::Index::l];
                ds2 = index2[DistillationNoise<FImpl>::Index::s];
                auto &solve2 = envGet(std::vector<FermionField>, par().vectorStem4);
                // full dilution assumed here, and that solve 4 comes from timeslice tKpi = d_{tKpi}
                dSolve2 = ds2 + nDS * dk2 + nDL * nDS * tKpi; 
                fermion4dtmp = solve2[dSolve2];
                ExtractSliceLocal(fermion3dtmp2,fermion4dtmp,0,t,Tdir);
                fermion3dtmp3 = g34*fermion3dtmp2;          
                fermion3dtmp2 = fermion3dtmp3;
                prop3dtmp = outerProduct(fermion3dtmp1,fermion3dtmp2);
                cplx3dtmp = trace(prop3dtmp)*MPhiPhi*ph3d;
                sliceSum(cplx3dtmp,buf,Tdir);                
                cache(0,0,0,0,0)=TensorRemove(buf[0]);
                block(0,0,id1,id2) = cache(0,0,0,0,0);                
            }
        }
        stopTimer("computation D-4quark");
        startTimer("serial I/O");
        LOG(Message) << "Starting serial IO for tH = " << tH << std::endl;
        DistilMatrixSetTimeSliceIo<ComplexF> block_relative(block_buf.data(), 1, nDL * nDS, nDL * nDS);
        std::string dataset_name = std::to_string(tKpi)+"-"+std::to_string(tKpi);
        gridHD->Barrier();
        for(int iIO=0; iIO<N_ranks; iIO++)
        {
            //gridHD->Barrier();
            if(iIO==i_rank)
            {
                LOG(Message) << "Writing from rank " << i_rank << std::endl;
                matrix_io.saveBlock(block_relative, 0, 0, 0, dataset_name, 0, nDL * nDS, std::to_string(tH));
            }
            gridHD->Barrier();
        }
        //gridHD->Barrier();
        stopTimer("serial I/O");
    }
    
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_DMeson4QuarkField_hpp_

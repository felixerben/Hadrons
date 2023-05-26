#include <Hadrons/Global.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Hadrons/DiskVector.hpp>
#include <Hadrons/TimerArray.hpp>
#include <Hadrons/DistilMatrix.hpp>

using namespace Grid;
using namespace Hadrons;

#define TIME_MOD(t) (((t) + nT) % nT)
#define CLOCK() std::cout << "Contractor : " << std::setw(10) << tAr.getDTimer("total")/1e6 << " s : "

typedef std::vector<ComplexD> Corr;

namespace ContractorInputs
{
    class TrajRange : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(TrajRange,
                                        unsigned int, start,
                                        unsigned int, end,
                                        unsigned int, step);
    };

    class ContractionInput : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(ContractionInput,
                                        std::string,    outStem,
                                        std::vector<std::string>,       inputStem,
                                        std::vector<std::string>,       dmfSuffix,
                                        std::vector<std::string>,       taskList,
                                        std::string,                    diagrams,
                                        unsigned int,   nT,
                                        std::string,    tSources);
    };
} 

class ContractionMetadata : Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ContractionMetadata,
                                    unsigned int,               Traj,
                                    std::string,                OutStem,
                                    std::vector<std::string>,   DmfSuffix,
                                    std::vector<std::string>,   DmfInputFiles,
                                    unsigned int,               Nt,
                                    std::string,                TimeSources);
};
/*
struct ContractorPar
{
    ContractorInputs::TrajRange             traj;
    ContractorInputs::ContractionInput      input;
};

ContractionMetadata populateMetadata(ContractorPar& par, unsigned int traj, std::vector<std::string> taskFiles)
{
    ContractionMetadata md;
    md.Traj = traj;
    md.OutStem = par.input.outStem;
    md.DmfSuffix = par.input.dmfSuffix;
    md.Nt = par.input.nT;
    md.TimeSources = par.input.tSources;
    md.DmfInputFiles = taskFiles;
    return md;
}

std::vector<unsigned int> fetchTSources(ContractorPar& par)
{
    std::vector<unsigned int> tSources={};
    if(par.input.tSources == "0...")
        for(int t=0;t<par.input.nT;t++)   tSources.push_back(t);
    else
        tSources = strToVec<unsigned int>(par.input.tSources);
    return tSources;
}

void saveCorrelator(std::string outname, ContractorPar& par, 
    std::map<std::string,std::vector<Corr>>& result, unsigned int traj, std::vector<std::string> taskFiles)
{
    outname = par.input.outStem + outname + "." + std::to_string(traj) + ".h5";
    Hadrons::makeFileDir(outname);
    ContractionMetadata md = populateMetadata(par,traj,taskFiles);
    std::cout << "Saving correlator to " << outname << std::endl;
    ResultWriter writer(outname);
    writer.push("DistillationContraction");
    write(writer, "Metadata", md);
    writer.push("Correlators");
    std::vector<unsigned int> tSources = fetchTSources(par);
    for(auto it = result.begin(); it != result.end(); ++it)
    {
        std::string diagram = it->first;
        writer.push(diagram); //group with name of diagram
        for(int it_0=0 ; it_0<it->second.size() ; it_0++)   //looping through time sources
        {
            std::string time_source = std::to_string(tSources[it_0]);
            write(writer, time_source , it->second[it_0]);
        }
        writer.pop();
    }
}

template <typename T>
void dumpKeys(std::map<std::string, T> map)
{
    std::cout << " =============== key dump ===============" << std::endl;
    for(auto it = map.begin(); it != map.end(); ++it)
    {
        std::cout << it->first << std::endl; 
    }
    std::cout << " ==========================================" << std::endl;
}
*/

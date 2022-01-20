/*
 * Test_exact_distil.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Felix Erben <felix.erben@ed.ac.uk>
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

/**************************************************************************
 * This test showcases how to use the Hadrons xml reader to write
 * an input file for a Hadrons C++ program. An example input xml is
 * testXmlReader.xml in this directory.
**************************************************************************/
#include <typeinfo>
#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

namespace TestInputs
{
    class RunPar : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(RunPar, 
                                        std::string, runId);
    };

    class ConfigPar : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(ConfigPar,
                                        std::string,  fileStem,
                                        unsigned int, begin,
                                        unsigned int, end,
                                        unsigned int, step);
    };
    
    class OutputPar : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(OutputPar,
                                        std::string, resultDir,
                                        std::string, xmlName);
    };

    class DwfPar : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(DwfPar,
                                        double, mass,
                                        double, M5,
                                        unsigned int, Ls,
                                        double, residual);
    };
}

struct TestPar
{
    TestInputs::RunPar      runPar;
    TestInputs::ConfigPar   configPar;
    TestInputs::OutputPar   outputPar;
    TestInputs::DwfPar      dwfPar;
};

int main(int argc, char *argv[])
{

    // program expects parameter file as first argument of command line
    if (argc < 2)
    {
        std::cerr << "usage: " << argv[0] << " <parameter file>";
        std::cerr << std::endl;

        return EXIT_FAILURE;
    }
    std::string parFilename;
    parFilename = argv[1];

    TestPar testPar;
    XmlReader reader(parFilename);

    read(reader, "runPar",     testPar.runPar);
    read(reader, "configPar",  testPar.configPar);
    read(reader, "outputPar",  testPar.outputPar);
    read(reader, "dwfPar",     testPar.dwfPar);


    // initialization //////////////////////////////////////////////////////////
    Grid_init(&argc, &argv);
    HadronsLogError.Active(GridLogError.isActive());
    HadronsLogWarning.Active(GridLogWarning.isActive());
    HadronsLogMessage.Active(GridLogMessage.isActive());
    HadronsLogIterative.Active(GridLogIterative.isActive());
    HadronsLogDebug.Active(GridLogDebug.isActive());
    LOG(Message) << "Grid initialized" << std::endl;
    
    // global parameters
    Application application;
    Application::GlobalPar globalPar;
    globalPar.runId                         = testPar.runPar.runId;
    globalPar.trajCounter.start             = testPar.configPar.begin;
    globalPar.trajCounter.end               = testPar.configPar.end;
    globalPar.trajCounter.step              = testPar.configPar.step;
    globalPar.database.applicationDb        = "test-input.db";
    globalPar.database.resultDb             = "test-input-results.db";
    globalPar.database.restoreSchedule      = false;
    globalPar.database.restoreModules       = false;
    globalPar.database.restoreMemoryProfile = false;
    globalPar.database.makeStatDb           = false;
    application.setPar(globalPar);
    // you can run this test with a random field or a nersc file
    if(testPar.configPar.fileStem.empty())
    {
        application.createModule<MGauge::Random>("gauge");
    }
    else
    {
        MIO::LoadNersc::Par loadPar;
        loadPar.file = testPar.configPar.fileStem;
        application.createModule<MIO::LoadNersc>("gauge", loadPar);
    }
    // DWF actions
    MAction::DWF::Par actionPar;
    actionPar.gauge = "gauge";
    actionPar.Ls    = testPar.dwfPar.Ls; 
    actionPar.M5    = testPar.dwfPar.M5;
    actionPar.mass  = testPar.dwfPar.mass;
    actionPar.boundary = "1 1 1 -1"; // anti-periodic bcs in time
    actionPar.twist = "0. 0. 0. 0."; // no twist
    application.createModule<MAction::DWF>("dwf", actionPar);
        
    // solvers
    MSolver::RBPrecCG::Par solverPar;
    solverPar.action       = "dwf";
    solverPar.residual     = testPar.dwfPar.residual;
    solverPar.maxIteration = 10000;
    application.createModule<MSolver::RBPrecCG>("cg", solverPar);

    // execution
    application.saveParameterFile(testPar.outputPar.xmlName);
    application.run();

    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
  
    return EXIT_SUCCESS;
}

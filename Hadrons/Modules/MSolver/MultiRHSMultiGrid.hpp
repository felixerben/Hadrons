#ifndef Hadrons_MSolver_MultiRHSMultiGrid_hpp_
#define Hadrons_MSolver_MultiRHSMultiGrid_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         MultiRHSMultiGrid                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class MultiRHSMultiGridPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MultiRHSMultiGridPar,
                                    unsigned int, i);
};

template <typename FImpl>
class TMultiRHSMultiGrid: public Module<MultiRHSMultiGridPar>
{
public:
    // constructor
    TMultiRHSMultiGrid(const std::string name);
    // destructor
    virtual ~TMultiRHSMultiGrid(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(MultiRHSMultiGrid, TMultiRHSMultiGrid<FIMPL>, MSolver);

/******************************************************************************
 *                 TMultiRHSMultiGrid implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TMultiRHSMultiGrid<FImpl>::TMultiRHSMultiGrid(const std::string name)
: Module<MultiRHSMultiGridPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TMultiRHSMultiGrid<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TMultiRHSMultiGrid<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TMultiRHSMultiGrid<FImpl>::setup(void)
{
    
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TMultiRHSMultiGrid<FImpl>::execute(void)
{
    
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_MultiRHSMultiGrid_hpp_

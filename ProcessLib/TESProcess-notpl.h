#ifndef PROCESS_LIB_TESPROCESS_NOTPL_H_
#define PROCESS_LIB_TESPROCESS_NOTPL_H_


#include "MaterialsLib/adsorption/adsorption.h"


namespace ProcessLib
{

namespace TES
{


struct Materials
{
    Ads::Adsorption* adsorption;
};


class TESProcessInterface
{
public:
    Materials const& getMaterials() const {
        return _materials;
    }

    virtual ~TESProcessInterface() = default;

protected:
    Materials _materials;
};

} // namespace TES

} // namespace ProcessLib

#endif // PROCESS_LIB_TESPROCESS_NOTPL_H_

//! Copyright : see license.txt
//!
//! \brief
//

#include "../../AlgoPacking/src/StdHeaders.hxx"

#include "MesoStructure/RecurStructure.hxx"

#include "MeropeNamespace.hxx"


namespace merope {

std::function<PhaseType(PhaseType, PhaseType)> auxiMicroStructure::replace_by(
    map<PhaseType, PhaseType> pair) {
    std::function<PhaseType(PhaseType, PhaseType)> result = [pair](PhaseType ph1, PhaseType ph2) {
        if (ph2 == 1) {
            auto search = pair.find(ph1);
            if (search != pair.end()) {
                return search->second;
            }
        }
        return ph1;
        };
    return result;
}

}  // namespace merope

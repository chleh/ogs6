/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_SMALLDEFORMATIONPROCESS_FWD_H_
#define PROCESS_LIB_SMALLDEFORMATIONPROCESS_FWD_H_

#include "SmallDeformationProcess.h"
#include "NumericsConfig.h"

extern template class ProcessLib::SmallDeformationProcess<GlobalSetupType, 1>;
extern template class ProcessLib::SmallDeformationProcess<GlobalSetupType, 2>;
extern template class ProcessLib::SmallDeformationProcess<GlobalSetupType, 3>;

#endif  // PROCESS_LIB_SMALLDEFORMATIONPROCESS_FWD_H_

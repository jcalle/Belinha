/*
 *  TPZYCWillamWarnke.cpp
 *  ElastoPlasticModels
 *
 *  Created by Diogo Cecilio on 12/13/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZYCWillamWarnke.h"

int TPZYCWillamWarnke::ClassId() const{
    return Hash("TPZYCWillamWarnke");
}

template class TPZRestoreClass<TPZYCWillamWarnke>;

#pragma once

#include "samxx_const.h"
#include "samxx_utils.h"
#include "vars.h"
#include "p3_functions.hpp"
#include "p3_functions_f90.hpp"
#include "p3_f90.hpp"

void micro_p3_init();
void micro_p3_proc();
void micro_p3_diagnose();

//using namespace scream;
//using namespace scream::p3;
//using P3F        = Functions<Real, DefaultDevice>;
//P3F::P3LookupTables      lookup_tables_save;

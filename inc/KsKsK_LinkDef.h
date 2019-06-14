// Copyright (c) 2019-3-3 maxx
#ifndef INC_KSKSK_LINKDEF_H_
#define INC_KSKSK_LINKDEF_H_
#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class RooCoord+;

#pragma link C++ class RooPWAPdf+;
#pragma link C++ class MultiDalitPdf+;
// #pragma link C++ class RooPropogator+;
#pragma link C++ class PWAPlot+;
// #ifdef DOLAUROOFITSLAVE
// #pragma link C++ class LauRooFitSlave+;
// #endif
#pragma link C++ namespace LauConstants+;
#pragma link C++ namespace bes3plotstyle+;
#pragma link C++ namespace RooP4Vector+;
#pragma link C++ namespace RooSpinFactor+;
#pragma link C++ namespace RooBarrier+;
#pragma link C++ namespace LineShape+;
#pragma link C++ namespace DecayType+;
#pragma link C++ namespace Propagator+;
#pragma link C++ namespace Propagator::BW+;
#pragma link C++ namespace Propagator::RBW+;
#pragma link C++ namespace Propagator::Flatte+;
#pragma link C++ namespace Propagator::K1430_p+;
#pragma link C++ namespace Propagator::K1430_0+;
#pragma link C++ namespace Propagator::a980_0+;
#pragma link C++ namespace Propagator::a980_p+;

#endif
#endif  // INC_KSKSK_LINKDEF_H_

#include <TChain.h>


void make()
{
  TChain *ch = new TChain("ch", "analysis_dvcs_4He_t");
  ch->Add("/u/home/dupre/Out-EIChigh-M3-*.root/TOPEG");
  ch->MakeClass("analysis_dvcs_4He_t");
}


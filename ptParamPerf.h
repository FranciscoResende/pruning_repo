////////////////////////////////////////////////////////////////////////////
//                                                                        //
// ptParamPerf.h                                                               //
//                                                                        //
// Class for to be written in perfNtuple for ptParametrisation
// Created by                                                             //
// Konstantinos Ntekas (Konstantinos.Ntekas@cern.ch)                      //
// on                                                                     //
// 18.06.2020                                                             //
//                                                                        //
////////////////////////////////////////////////////////////////////////////
#ifndef PTPARAMPERF_H_INCLUDED
#define PTPARAMPERF_H_INCLUDED

// ROOT
#include "TObject.h"

// C++ std
#include <string>

namespace PERFNT {
class ptParamPerf : public TObject {
 private:
 public:
  /** Constructor */
  ptParamPerf();
  /** Destructor */
  virtual ~ptParamPerf() { ; }

  // event number
  int event;
  // SLc ID. When seeded by offline muon for a given event can be >1
  // if more than one muon or the muon overlap L/S sector, each which case 2 muons are created, one for each sector
  int iSL;
  // Offline muon quality. Needed?
  int quality;
  // Offline muon author. Needed?
  int author;
  // Offline muon is real? Needed?
  bool isReal;
  // Offline muon is combined? Needed?
  bool isCombined;
  // Offline muon is standalone? Needed?
  bool isStandalone;
  // True muon type. Needed?
  int type_truth;
  // True muon origin. Needed?
  int origin_truth;
  // Offline muon charge
  short charge;
  // True muon charge
  short charge_truth;
  // Offline muon pT
  float pt;
  // True muon pT
  float pt_truth;
  // Offline muon eta
  float eta;
  // True muon eta
  float eta_truth;
  // Offline muon phi
  float phi;
  // True muon phi
  float phi_truth;
  // Phi center of the sector [-Pi,Pi]
  float sector_phi;
  // Phi of the Roi (inner most segment)
  double roi_phi;
  // Phi of the Roi (inner most segment)
  double roi_phi_inn;
  // Phi of the Roi (inner most segment)
  double roi_phi_mid;
  // Phi of the Roi (inner most segment)
  double roi_phi_out;
  // Phi of the Roi (inner most segment)
  double roi_phi_ext;
  // Phimod = abs(phi-sector_phi)
  float phi_mod;
  // Phi ID of the sector
  short phi_index;
  // MS only muon pt, eta, phi, q
  float ms_pt;
  float ms_phi;
  float ms_eta;
  float ms_charge;

  // ComboName
  std::string comboName;
  // Candidate status from MdtMuonCandidate
  int candStatus;
  // Station has a segment
  std::vector<bool> hasSegment;
  // list of chamber names the offline muons traverses
  std::vector<std::string> chamberNames;
  // list of chamber IDs the offline muons traverses
  std::vector<short> chamberIds;
  // list of chamber iEta index
  std::vector<short> etaIdx;
  // list of chamber IDs the offline muons traverses - Use L0MDT numbering
  std::vector<short> chamberIds_L0MDT;
  // polar angle global coordinate of segment - orginal from offline segment - ie segment R & roiPhi
  std::vector<float> beta;
  // polar angle global coordinate of segment - at the RoiPhi - identical to beta for offline segment
  std::vector<float> beta_roiPhi;
  // polar angle global coordinate of segment - at the sector phi center - Used for calculated deltaBeta
  std::vector<float> beta_sectorPhi;
  // global position radial distance from offline segment at its given phi, with radial distance depending on the hit
  // content
  std::vector<float> r;
  // global position Z distance from offline segment at its given phi, with radial distance depending on the hit content
  std::vector<float> z;
  // global position radial distance extrapolated to common per sector R(barrel)/z(endcap), and projected to sector phi
  // center. Used for sagitta calc.
  std::vector<float> r_ext;
  // global position Z distance extrapolated to common per sector R(barrel)/z(endcap), and projected to sector phi
  // center . Used for sagitta calc.
  std::vector<float> z_ext;
  // global position radial distance extrapolated to common per sector R(barrel)/z(endcap), and at the phi of the
  // offline segment aka roiPhi
  std::vector<float> r_cor;
  // global position Z distance extrapolated to common per sector R(barrel)/z(endcap), and at the phi of the offline
  // segment aka roiPhi
  std::vector<float> z_cor;
  // Chamber combination ID built from the chamberIDs
  unsigned long long comboID;
  // number of MDT stations for offline muon
  int nStations;
  // sagitta stationR - sector centerPhi
  float sagitta;
  // sagitta  station R - Roi phi
  float sagitta_roiPhi;
  // deltaBeta  for 2 station, stationR - sector centerPhi
  float deltaBeta;
  // deltaBeta  for 2 station, station R - Roi phi
  float deltaBeta_roiPhi;
  // small or large sector
  bool is_small_sector;
  // has only barrel chambers?
  bool is_barrel;
  // has only endcap chambers?
  bool is_endcap;
  // has barrel+endcap chambers?
  bool is_transition_region;
  // is A side
  bool is_Aside;
  // is a two station candidate (>1 chambers)?
  bool is_two_station;
  // is a three station candidate (>2 chambers)?
  bool is_three_station;
  // is a case where the CSC segment attached to the offline muon is in a different sector than the MDT chambers
  // attached
  bool is_neighbour_csc_sector;

  // Online calculated values using the loaded parametrisation - station R, sector Phi
  // online pt
  float online_pt;
  // online pt_s
  float online_pt_s;
  // online pt_sp
  float online_pt_sp;
  // online charge
  int online_q;
  // online eta
  float online_eta;

  // online eta using roiPhi
  float online_eta_roiPhi;

  void print();

  // Initialize the vectors to fix size
  void init() {
    hasSegment.clear();
    chamberNames.clear();
    chamberIds.clear();
    etaIdx.clear();
    chamberIds_L0MDT.clear();
    beta.clear();
    beta_roiPhi.clear();
    beta_sectorPhi.clear();
    r.clear();
    z.clear();
    r_ext.clear();
    z_ext.clear();
    r_cor.clear();
    z_cor.clear();

    if (hasSegment.size() == 0) {
      hasSegment.resize(4, false);
      chamberNames.resize(4, "NA");
      chamberIds.resize(4, 99);
      etaIdx.resize(4, 10);
      chamberIds_L0MDT.resize(4, 99);
      beta.resize(4, -99);
      beta_roiPhi.resize(4, -99);
      beta_sectorPhi.resize(4, -99);
      r.resize(4, -99);
      z.resize(4, -99);
      r_ext.resize(4, -99);
      z_ext.resize(4, -99);
      r_cor.resize(4, -99);
      z_cor.resize(4, -99);
    }
  }

  void clear() {
    event        = -99;
    iSL          = -99;
    author       = -99;
    quality      = -99;
    isReal       = false;
    isCombined   = false;
    isStandalone = false;
    type_truth   = -99;
    origin_truth = -99;
    charge       = -99.;
    charge_truth = -99.;
    pt           = -99;
    pt_truth     = -99;
    eta          = -99;
    eta_truth    = -99;
    phi          = -99;
    phi_truth    = -99;
    sector_phi   = -99;
    roi_phi      = -99;
    roi_phi_inn  = -99;
    roi_phi_mid  = -99;
    roi_phi_out  = -99;
    roi_phi_ext  = -99;
    phi_mod      = -99;
    phi_index    = -99;
    ms_pt        = -99;
    ms_phi       = -99;
    ms_eta       = -99;
    ms_charge    = -99;
    comboName.clear();
    candStatus = 3;
    init();
    comboID                 = 0;
    nStations               = -99;
    sagitta                 = -99;
    sagitta_roiPhi          = -99;
    deltaBeta               = -99;
    deltaBeta_roiPhi        = -99;
    is_small_sector         = false;
    is_barrel               = false;
    is_endcap               = false;
    is_Aside                = false;
    is_transition_region    = false;
    is_two_station          = false;
    is_three_station        = false;
    is_neighbour_csc_sector = false;

    online_pt         = -99;
    online_pt_s       = -99;
    online_pt_sp      = -99;
    online_eta        = -99;
    online_eta_roiPhi = -99;
    online_q          = 0;
  }
  ClassDef(ptParamPerf, 1);
};
}  // namespace PERFNT
#endif

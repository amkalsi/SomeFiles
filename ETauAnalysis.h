// -*- C++ -*-
//
// Package:    Plots/ETauAnalysis
// Class:      ETauAnalysis
// 
/**\class ETauAnalysis ETauAnalysis.cc Plots/ETauAnalysis/plugins/ETauAnalysis.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  kaur amandeepkalsi
//         Created:  Fri, 11 Dec 2015 06:28:12 GMT
//
//


// system include files
#include <memory>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>                                                                  
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <iostream>
#include <cstring>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <map>
#include <sys/stat.h>
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"           
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "IIHEAnalysis.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/SimpleJetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/SimpleJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
//#include "RoccoR.h"
const float m_el = 0.000511 ;

namespace PU_reReco_Morind17{
	float MC_pileup_weight(int NumTrueInteraction, int scale, TString dataset){

		if (NumTrueInteraction < 0 || NumTrueInteraction > 74 ) return 1;

		double PU_data_all_noScale[75] = {
			6.54008e-06,
			2.29383e-05,
			6.32223e-05,
			8.55796e-05,
			0.000122592,
			0.000164214,
			0.000191738,
			0.000353072,
			0.000965735,
			0.00215544,
			0.00484612,
			0.00986199,
			0.0165083,
			0.0240058,
			0.0321661,
			0.0407818,
			0.0481844,
			0.0532395,
			0.0561219,
			0.0575573,
			0.058412,
			0.0588587,
			0.0583078,
			0.056491,
			0.0537587,
			0.0504445,
			0.0466722,
			0.0425747,
			0.0383287,
			0.0340574,
			0.0298202,
			0.0256705,
			0.0216922,
			0.0179859,
			0.0146378,
			0.011698,
			0.00917774,
			0.00705846,
			0.00530639,
			0.00388441,
			0.00275715,
			0.00189013,
			0.00124723,
			0.000790064,
			0.000479456,
			0.00027833,
			0.000154405,
			8.18145e-05,
			4.14101e-05,
			2.00435e-05,
			9.30691e-06,
			4.17785e-06,
			1.84606e-06,
			8.3504e-07,
			4.1498e-07,
			2.45829e-07,
			1.77914e-07,
			1.48846e-07,
			1.33924e-07,
			1.23839e-07,
			1.1526e-07,
			1.07079e-07,
			9.89863e-08,
			9.09467e-08,
			8.30145e-08,
			7.52676e-08,
			6.77837e-08,
			6.06312e-08,
			5.38664e-08,
			4.75325e-08,
			4.16595e-08,
			3.62653e-08,
			3.1356e-08,
			2.69281e-08,
			2.29691e-08
		};

		double PU_data_all_ScaleUp[75] ={
			6.37264e-06,
			1.80616e-05,
			5.98112e-05,
			7.5202e-05,
			0.000111511,
			0.000147892,
			0.000174897,
			0.000247635,
			0.000651364,
			0.00147902,
			0.00317738,
			0.00673879,
			0.0121433,
			0.0186266,
			0.0256861,
			0.0333807,
			0.0411281,
			0.0473936,
			0.0515615,
			0.0539179,
			0.0551224,
			0.0558812,
			0.0562746,
			0.0557863,
			0.0541814,
			0.0517471,
			0.0487842,
			0.0454071,
			0.0417151,
			0.0378563,
			0.0339515,
			0.0300625,
			0.0262291,
			0.0225106,
			0.018988,
			0.0157428,
			0.0128347,
			0.0102913,
			0.00811139,
			0.00627373,
			0.00474839,
			0.00350453,
			0.00251292,
			0.0017447,
			0.00116951,
			0.000755121,
			0.000468809,
			0.000279493,
			0.000159864,
			8.76833e-05,
			4.61173e-05,
			2.32764e-05,
			1.12984e-05,
			5.30207e-06,
			2.43424e-06,
			1.12208e-06,
			5.4596e-07,
			3.01607e-07,
			1.99891e-07,
			1.56737e-07,
			1.3654e-07,
			1.24901e-07,
			1.16254e-07,
			1.08552e-07,
			1.01101e-07,
			9.37102e-08,
			8.63684e-08,
			7.91252e-08,
			7.20462e-08,
			6.51962e-08,
			5.86329e-08,
			5.2404e-08,
			4.6547e-08,
			4.10889e-08,
			3.60464e-08
		};

		double PU_data_all_ScaleDown[75] ={
			6.77593e-06,
			2.92832e-05,
			6.65038e-05,
			9.76882e-05,
			0.00013671,
			0.000180566,
			0.000221755,
			0.000546682,
			0.00142169,
			0.00328014,
			0.00746986,
			0.0140533,
			0.0219725,
			0.0306384,
			0.0400658,
			0.0487539,
			0.0549158,
			0.0584733,
			0.0602103,
			0.061182,
			0.0616915,
			0.0610662,
			0.0589986,
			0.0559171,
			0.0521906,
			0.0479588,
			0.0433993,
			0.0387163,
			0.0340311,
			0.0294056,
			0.0249178,
			0.0206818,
			0.0168142,
			0.0133971,
			0.0104626,
			0.0079996,
			0.00597156,
			0.00433387,
			0.0030433,
			0.00205819,
			0.00133527,
			0.000828375,
			0.000490264,
			0.000276337,
			0.000148173,
			7.55428e-05,
			3.66294e-05,
			1.692e-05,
			7.48062e-06,
			3.20277e-06,
			1.36555e-06,
			6.1528e-07,
			3.21628e-07,
			2.09204e-07,
			1.64835e-07,
			1.44633e-07,
			1.32489e-07,
			1.22806e-07,
			1.13769e-07,
			1.04868e-07,
			9.6027e-08,
			8.73029e-08,
			7.87891e-08,
			7.05791e-08,
			6.27551e-08,
			5.53837e-08,
			4.85148e-08,
			4.2182e-08,
			3.64034e-08,
			3.11832e-08,
			2.65132e-08,
			2.23752e-08,
			1.8743e-08,
			1.55839e-08,
			1.28612e-08
		};

		double PU_data_B2F_noScale[75] ={
			1.73901e-07,
			1.45481e-05,
			5.36313e-05,
			8.93304e-05,
			0.000134637,
			0.000184277,
			0.00023659,
			0.000473075,
			0.00151712,
			0.00362251,
			0.00753906,
			0.0133722,
			0.0205964,
			0.0289851,
			0.0383967,
			0.0481983,
			0.0563721,
			0.0619701,
			0.065389,
			0.0669141,
			0.0666799,
			0.0650214,
			0.0623224,
			0.0588295,
			0.0546245,
			0.0496801,
			0.0440732,
			0.0381087,
			0.0322039,
			0.0266868,
			0.0217181,
			0.0173439,
			0.0135688,
			0.0103829,
			0.00775837,
			0.00564855,
			0.00399533,
			0.00273655,
			0.00180926,
			0.00115131,
			0.000703515,
			0.000412097,
			0.000231157,
			0.000124109,
			6.37956e-05,
			3.14347e-05,
			1.48954e-05,
			6.83995e-06,
			3.09966e-06,
			1.44331e-06,
			7.43085e-07,
			4.59278e-07,
			3.47549e-07,
			3.02845e-07,
			2.82329e-07,
			2.69486e-07,
			2.58339e-07,
			2.46917e-07,
			2.34694e-07,
			2.21652e-07,
			2.0794e-07,
			1.93758e-07,
			1.79318e-07,
			1.64827e-07,
			1.50479e-07,
			1.36447e-07,
			1.22884e-07,
			1.09919e-07,
			9.76561e-08,
			8.61735e-08,
			7.55264e-08,
			6.5747e-08,
			5.68468e-08,
			4.88192e-08,
			4.16418e-08
		};
		double PU_data_B2F_ScaleUp[75] ={
			1.08089e-07,
			1.23427e-05,
			4.69336e-05,
			7.79185e-05,
			0.000122499,
			0.000163075,
			0.000210873,
			0.000319363,
			0.000963976,
			0.00246186,
			0.00521193,
			0.00977448,
			0.0157329,
			0.0228296,
			0.0308535,
			0.0397307,
			0.0484633,
			0.0553471,
			0.0600062,
			0.062836,
			0.0640269,
			0.0636974,
			0.062143,
			0.0596839,
			0.0565276,
			0.05275,
			0.0483277,
			0.0432995,
			0.0378875,
			0.032437,
			0.0272569,
			0.0225245,
			0.0183036,
			0.0146066,
			0.0114302,
			0.00875935,
			0.00656306,
			0.00479726,
			0.0034116,
			0.00235371,
			0.00157105,
			0.00101208,
			0.000628058,
			0.000374911,
			0.000215092,
			0.000118558,
			6.27983e-05,
			3.19989e-05,
			1.57264e-05,
			7.50006e-06,
			3.5194e-06,
			1.67501e-06,
			8.56232e-07,
			5.07096e-07,
			3.6293e-07,
			3.03805e-07,
			2.77834e-07,
			2.63684e-07,
			2.53039e-07,
			2.42869e-07,
			2.32186e-07,
			2.20768e-07,
			2.08663e-07,
			1.96013e-07,
			1.8299e-07,
			1.69774e-07,
			1.56533e-07,
			1.43429e-07,
			1.30607e-07,
			1.18193e-07,
			1.06296e-07,
			9.50048e-08,
			8.43869e-08,
			7.44919e-08,
			6.53503e-08
		};
		double PU_data_B2F_ScaleDown[75] ={
			2.78703e-07,
			1.73561e-05,
			6.08057e-05,
			0.000103261,
			0.000149929,
			0.000207897,
			0.000279105,
			0.000776015,
			0.00233792,
			0.00540425,
			0.0108027,
			0.018051,
			0.0267677,
			0.0366905,
			0.0475034,
			0.057169,
			0.063952,
			0.068108,
			0.0700527,
			0.0699499,
			0.0681791,
			0.0652034,
			0.061318,
			0.0566103,
			0.0510523,
			0.0447795,
			0.0382036,
			0.0318172,
			0.0259555,
			0.0207543,
			0.0162424,
			0.012417,
			0.00925627,
			0.00671364,
			0.00472343,
			0.00321182,
			0.00210286,
			0.00132115,
			0.000794225,
			0.000455912,
			0.000249568,
			0.000130203,
			6.47576e-05,
			3.07495e-05,
			1.39952e-05,
			6.16621e-06,
			2.69456e-06,
			1.23294e-06,
			6.47737e-07,
			4.23351e-07,
			3.39038e-07,
			3.054e-07,
			2.88346e-07,
			2.75676e-07,
			2.63361e-07,
			2.50255e-07,
			2.36184e-07,
			2.21285e-07,
			2.0579e-07,
			1.89955e-07,
			1.7403e-07,
			1.58252e-07,
			1.42831e-07,
			1.27952e-07,
			1.1377e-07,
			1.00407e-07,
			8.79545e-08,
			7.64737e-08,
			6.59975e-08,
			5.65334e-08,
			4.8067e-08,
			4.05652e-08,
			3.39801e-08,
			2.82528e-08,
			2.33167e-08
		};

		double PU_data_GH_noScale[75] ={
			1.43711e-05,
			3.3259e-05,
			7.50201e-05,
			8.09657e-05,
			0.000107776,
			0.000139534,
			0.000136566,
			0.000205456,
			0.00028748,
			0.000350813,
			0.00153358,
			0.00554418,
			0.0114796,
			0.0178807,
			0.0245018,
			0.0316588,
			0.0381129,
			0.0425001,
			0.0447227,
			0.0460476,
			0.0482416,
			0.051278,
			0.0533696,
			0.0536145,
			0.0526937,
			0.0513847,
			0.0498692,
			0.0480682,
			0.0458627,
			0.0431239,
			0.0397865,
			0.035913,
			0.0316847,
			0.0273383,
			0.0231,
			0.0191393,
			0.0155526,
			0.0123748,
			0.00960817,
			0.00724636,
			0.00528331,
			0.00370823,
			0.0024971,
			0.00160925,
			0.000990756,
			0.000582034,
			0.000326015,
			0.00017404,
			8.85354e-05,
			4.29234e-05,
			1.98412e-05,
			8.75203e-06,
			3.68936e-06,
			1.48969e-06,
			5.78153e-07,
			2.16729e-07,
			7.89841e-08,
			2.82112e-08,
			9.9682e-09,
			3.51973e-09,
			1.25491e-09,
			4.56558e-10,
			1.71273e-10,
			6.68678e-11,
			2.73251e-11,
			1.16828e-11,
			5.18929e-12,
			2.36754e-12,
			1.096e-12,
			5.09461e-13,
			2.35966e-13,
			1.08341e-13,
			4.91521e-14,
			2.1992e-14,
			9.69335e-15
		};
		double PU_data_GH_ScaleUp[75] ={
			1.40786e-05,
			2.50963e-05,
			7.56518e-05,
			7.18604e-05,
			9.79943e-05,
			0.000129215,
			0.000130643,
			0.000159403,
			0.000266824,
			0.000270027,
			0.000674718,
			0.00300462,
			0.00772775,
			0.0134566,
			0.0193298,
			0.0255697,
			0.0321053,
			0.0376102,
			0.0411737,
			0.042948,
			0.044169,
			0.0462664,
			0.0490558,
			0.050992,
			0.0512954,
			0.0505135,
			0.0493456,
			0.0479996,
			0.0464234,
			0.0445225,
			0.0421864,
			0.0393349,
			0.0359783,
			0.0322331,
			0.0282848,
			0.0243332,
			0.0205493,
			0.0170494,
			0.0138926,
			0.0110957,
			0.0086568,
			0.00657046,
			0.00483147,
			0.00342967,
			0.00234352,
			0.00153815,
			0.000968239,
			0.000583933,
			0.000337166,
			0.000186316,
			9.85166e-05,
			4.98479e-05,
			2.41433e-05,
			1.12003e-05,
			4.98214e-06,
			2.12862e-06,
			8.75779e-07,
			3.48256e-07,
			1.34514e-07,
			5.07874e-08,
			1.88863e-08,
			6.9762e-09,
			2.58259e-09,
			9.66975e-10,
			3.69571e-10,
			1.45489e-10,
			5.94571e-11,
			2.53396e-11,
			1.12541e-11,
			5.17623e-12,
			2.44114e-12,
			1.16784e-12,
			5.61448e-13,
			2.69324e-13,
			1.28276e-13
		};
		double PU_data_GH_ScaleDown[75] ={
			1.47681e-05,
			4.39546e-05,
			7.35131e-05,
			9.08336e-05,
			0.000120451,
			0.000146947,
			0.00015121,
			0.000264582,
			0.000294651,
			0.000667302,
			0.00337014,
			0.00913579,
			0.0160739,
			0.0231937,
			0.0309169,
			0.0384027,
			0.0438004,
			0.0466218,
			0.0481033,
			0.0503968,
			0.0537112,
			0.0559769,
			0.0561455,
			0.0550643,
			0.0535908,
			0.0518696,
			0.0497904,
			0.0472027,
			0.0439647,
			0.0400475,
			0.0355892,
			0.0308481,
			0.0261111,
			0.0216183,
			0.0175223,
			0.013889,
			0.0107304,
			0.00803978,
			0.00580986,
			0.00402913,
			0.00267079,
			0.00168719,
			0.00101368,
			0.000578431,
			0.000313223,
			0.000160882,
			7.83722e-05,
			3.62165e-05,
			1.58857e-05,
			6.62171e-06,
			2.62824e-06,
			9.96461e-07,
			3.62569e-07,
			1.27438e-07,
			4.36388e-08,
			1.4708e-08,
			4.93493e-09,
			1.66809e-09,
			5.74938e-10,
			2.04518e-10,
			7.5916e-11,
			2.96169e-11,
			1.21442e-11,
			5.19229e-12,
			2.28453e-12,
			1.02006e-12,
			4.56865e-13,
			2.03548e-13,
			8.97288e-14,
			3.90098e-14,
			1.66947e-14,
			7.0257e-15,
			2.90561e-15,
			1.1805e-15,
			4.71041e-16
		};

		double PU_data_B2D[75] ={
			1.61496e-07,
			1.23652e-05,
			4.82966e-05,
			9.65592e-05,
			0.000150797,
			0.000204849,
			0.000274118,
			0.000655443,
			0.00226352,
			0.00539922,
			0.0105562,
			0.0177837,
			0.0272428,
			0.0374065,
			0.0467972,
			0.0557235,
			0.0638455,
			0.0701363,
			0.0740697,
			0.0755859,
			0.0747888,
			0.0717894,
			0.0669529,
			0.0608304,
			0.0538144,
			0.046114,
			0.0380152,
			0.0300213,
			0.022702,
			0.0164566,
			0.0114341,
			0.00759726,
			0.00481316,
			0.00290295,
			0.00166735,
			0.000913053,
			0.000477292,
			0.000238569,
			0.000114339,
			5.27814e-05,
			2.36331e-05,
			1.03878e-05,
			4.58664e-06,
			2.1282e-06,
			1.11697e-06,
			7.12405e-07,
			5.54821e-07,
			4.95004e-07,
			4.72643e-07,
			4.63813e-07,
			4.59037e-07,
			4.54142e-07,
			4.47847e-07,
			4.39559e-07,
			4.29e-07,
			4.16127e-07,
			4.01048e-07,
			3.83965e-07,
			3.65149e-07,
			3.44912e-07,
			3.2359e-07,
			3.01524e-07,
			2.79054e-07,
			2.56504e-07,
			2.34175e-07,
			2.12338e-07,
			1.91232e-07,
			1.71056e-07,
			1.51972e-07,
			1.34103e-07,
			1.17534e-07,
			1.02315e-07,
			8.84649e-08,
			7.59723e-08,
			6.48029e-08
		};
		double PU_data_EF[75] ={
			1.96205e-07,
			1.84728e-05,
			6.32227e-05,
			7.63337e-05,
			0.000105581,
			0.00014729,
			0.000169117,
			0.000145193,
			0.000175153,
			0.000428122,
			0.00211455,
			0.0054405,
			0.00864683,
			0.0138442,
			0.0232934,
			0.0346685,
			0.0429354,
			0.0472878,
			0.0497816,
			0.0513228,
			0.0521008,
			0.0528532,
			0.0539972,
			0.0552321,
			0.0560809,
			0.0560917,
			0.054965,
			0.0526492,
			0.0492875,
			0.0450798,
			0.0402079,
			0.0348676,
			0.0293108,
			0.0238312,
			0.0187095,
			0.0141626,
			0.0103205,
			0.00722773,
			0.00485659,
			0.00312639,
			0.00192589,
			0.00113434,
			0.000638514,
			0.000343421,
			0.000176487,
			8.66711e-05,
			4.06786e-05,
			1.82477e-05,
			7.82282e-06,
			3.20436e-06,
			1.25378e-06,
			4.68513e-07,
			1.67219e-07,
			5.70427e-08,
			1.86251e-08,
			5.83547e-09,
			1.76096e-09,
			5.14347e-10,
			1.46257e-10,
			4.07327e-11,
			1.11694e-11,
			3.02621e-12,
			8.1103e-13,
			2.14729e-13,
			5.59942e-14,
			1.43255e-14,
			3.58157e-15,
			8.71958e-16,
			2.06131e-16,
			4.7236e-17,
			1.04376e-17,
			2.23209e-18,
			4.55852e-19,
			8.80598e-20,
			1.88957e-20};

		double PU_MC[75] =  {
			1.78653e-05 ,
			2.56602e-05 ,
			5.27857e-05 ,
			8.88954e-05 ,
			0.000109362 ,
			0.000140973 ,
			0.000240998 ,
			0.00071209 ,
			0.00130121 ,
			0.00245255 ,
			0.00502589 ,
			0.00919534 ,
			0.0146697 ,
			0.0204126 ,
			0.0267586 ,
			0.0337697 ,
			0.0401478 ,
			0.0450159 ,
			0.0490577 ,
			0.0524855 ,
			0.0548159 ,
			0.0559937 ,
			0.0554468 ,
			0.0537687 ,
			0.0512055 ,
			0.0476713 ,
			0.0435312 ,
			0.0393107 ,
			0.0349812 ,
			0.0307413 ,
			0.0272425 ,
			0.0237115 ,
			0.0208329 ,
			0.0182459 ,
			0.0160712 ,
			0.0142498 ,
			0.012804 ,
			0.011571 ,
			0.010547 ,
			0.00959489 ,
			0.00891718 ,
			0.00829292 ,
			0.0076195 ,
			0.0069806 ,
			0.0062025 ,
			0.00546581 ,
			0.00484127 ,
			0.00407168 ,
			0.00337681 ,
			0.00269893 ,
			0.00212473 ,
			0.00160208 ,
			0.00117884 ,
			0.000859662 ,
			0.000569085 ,
			0.000365431 ,
			0.000243565 ,
			0.00015688 ,
			9.88128e-05 ,
			6.53783e-05 ,
			3.73924e-05 ,
			2.61382e-05 ,
			2.0307e-05 ,
			1.73032e-05 ,
			1.435e-05 ,
			1.36486e-05 ,
			1.35555e-05 ,
			1.37491e-05 ,
			1.34255e-05 ,
			1.33987e-05 ,
			1.34061e-05 ,
			1.34211e-05 ,
			1.34177e-05 ,
			1.32959e-05 ,
			1.33287e-05
		};

		if(scale ==0 && dataset=="all"){ return float(PU_data_all_noScale[NumTrueInteraction]/PU_MC[NumTrueInteraction]) ;}
		else if(scale ==1 && dataset=="all"){ return float(PU_data_all_ScaleUp[NumTrueInteraction]/PU_MC[NumTrueInteraction]) ;}
		else if(scale ==2 && dataset=="all"){ return float(PU_data_all_ScaleDown[NumTrueInteraction]/PU_MC[NumTrueInteraction]) ;}
		else if(scale ==0 && dataset=="B2F"){ return float(PU_data_B2F_noScale[NumTrueInteraction]/PU_MC[NumTrueInteraction]) ;}
		else if(scale ==1 && dataset=="B2F"){ return float(PU_data_B2F_ScaleUp[NumTrueInteraction]/PU_MC[NumTrueInteraction]) ;}
		else if(scale ==2 && dataset=="B2F"){ return float(PU_data_B2F_ScaleDown[NumTrueInteraction]/PU_MC[NumTrueInteraction]) ;}
		else if(scale ==0 && dataset=="GH"){ return float(PU_data_GH_noScale[NumTrueInteraction]/PU_MC[NumTrueInteraction]) ;}
		else if(scale ==1 && dataset=="GH"){ return float(PU_data_GH_ScaleUp[NumTrueInteraction]/PU_MC[NumTrueInteraction]) ;}
		else if(scale ==2 && dataset=="GH"){ return float(PU_data_GH_ScaleDown[NumTrueInteraction]/PU_MC[NumTrueInteraction]) ;}
		else if(scale ==0 && dataset=="B2D"){ return float(PU_data_B2D[NumTrueInteraction]/PU_MC[NumTrueInteraction]) ;}
		else if(scale ==0 && dataset=="EF"){ return float(PU_data_EF[NumTrueInteraction]/PU_MC[NumTrueInteraction]) ;}
		else {return 1;}
	}
}
//
// class declaration
//
#ifdef __MAKECINT__
#pragma link C++ class vector<TLorentzVector>+;
#endif
using namespace std;
template <typename T>
struct SortByPt
{    
	bool operator () (const T& a, const T& b) const {
		return a.first.Pt() > b.first.Pt();
	}  

};   

class ETauAnalysis : public edm::EDAnalyzer {
	public:
		explicit ETauAnalysis(const edm::ParameterSet&);
		~ETauAnalysis();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;

		//virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
		//virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
		//virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
		//virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

		// ----------member data ---------------------------
		//
		enum particleType {
			MUON = 0,
			ELECTRON = 1,
			TAU =2
		};     
		enum pairType {
			MuHad  = 0,
			EHad   = 1,
			HadHad = 2,
			MuMu   = 3,
			EE     = 4,
			EMu    = 5,
			EEPrompt = 6, // prompt Z->ee/mumu decays
			MuMuPrompt = 7,
			Other  = 8 // for e.g. h->bb
		};
		edm::Service<TFileService> fs;   
		TTree *myTree;
		float sumweight;
		bool isMuon(int type){if(type == MUON)return true; else return false;}
		bool isElectron(int type){if(type == ELECTRON)return true; else return false;}
		bool isTau(int type){if(type == TAU)return true; else return false;}
		int getPairType (int type1, int type2); // return pair type giving as input the particle types of the two composants
		bool checkBit (int word, int bitpos); // check bit "bitpos" in a word   

		bool EventSelector();
		double mtvalue ,pzeta , pzetavis, weighthis;
		int matchedgen1,matchedgen2,  mu1daumatch, Nmu1daumatch, genlevel, Ngenlevel, eventsel;
		bool OverLap05(TLorentzVector l1 , TLorentzVector l2, float conesize);
		float deltaPhi( float a, float b) ;
		float dR(float l1eta, float l1phi, float l2eta, float l2phi );
		bool muselection( int index1, float ptcut, float etacut, float isolation);
		bool tauselection(int index1, float ptcut, float etacut, float isolation , string isoWP, string muWP, string eleWP, string dmf);
		bool muonTauOverlap( int index1 );
		bool ExtraMuon(TLorentzVector, TLorentzVector );
		bool ExtraElectron(TLorentzVector , TLorentzVector , pair<double,double> &);
		bool PreselectionTauCuts( int tauindex);
		bool IsTriggerFired(int triggerbit, int triggernumber);
		int FindTriggerNumber(TString triggername);
		bool PassEleSelections( int muindex, TLorentzVector, int );
		bool PassHEEP(int eleindex);
		bool PassTauSelections( int , TLorentzVector );
		bool PassTauIsolation(int tauindex, int CR);
		bool BJets(int , int, TLorentzVector, TLorentzVector);
		bool OverlapWithMuons(int dau1index);
		bool OverlapWithElectrons(int dau1index);
		bool OverlapWithTaus(int dau1index);
		bool OverlapWithJets(int dau1index);
		bool MatchingToGenMuons(TLorentzVector , int &) ;
		bool MatchingToGenTaus(TLorentzVector , int &) ;

		float mTCalculation(float metx, float mety, float mupx, float mupy, float mupt);
		void initializePileupInfo();
		double getPileupWeight(float ntruePUInt);
		bool PassSubMuSelections( int , TLorentzVector );

		double GetCollinearMass(const TLorentzVector &tau, const TLorentzVector &mu,  const TLorentzVector MET) ;
		float PZeta(int muindex, int tauindex , float metpx, float metpy);
		float PZetaVis(int muindex, int tauindex);
		bool DiLeptonVeto(int muindex);
		bool GenSelection(TLorentzVector, TLorentzVector, int &);
		void HistoFiller(TH1F *histo, double value, double weight);
		int NumberBJets(int, int, TLorentzVector ,TLorentzVector) ;
		int ExtraMuonNumber(TLorentzVector ,TLorentzVector );
		int NumberExtraElectron(TLorentzVector , TLorentzVector);
		TFile *output_file;
		void HistoDec(const char*);
		bool GenSelection(TLorentzVector,TLorentzVector, int &maingen1, int &maingen2);
		int LeastDR(TLorentzVector, vector<int> genindx, bool, int);
		double FindMuonIDEff(double pt, double eta);
		double FindMuonISOEff(double pt, double eta);
		bool FillChain(TChain *chain, const char* inputFileList);
		TH1F *hPUmc = new TH1F("hPUmc", "hPUmc", 60, 0, 60);
		TH1F *hPUdata = new TH1F("hPUdata", "hPUdata", 60, 0, 60);
		int pairtype;
		bool gen;
		double pu_weight;
		vector <TLorentzVector> GoodMuons, GoodTaus, GoodElecs;
		vector <unsigned int> GoodMuonsIndex, GoodTausIndex, GoodElecsIndex;                                                                                                
		/*		TH1F *hEntries,*hrun,*hlumi,*hevt,*hNUP,*hdilepton_veto,*hextraelec_veto,*hextramuon_veto,*hisZmt;
				TH1F *hmu_pt,*hmu_en,*hmu_px,*hmu_py,*hmu_pz,*hmu_eta,*hmu_phi,*hmu_charge,*hmu_d0;
				TH1F *htau_pt,*htau_eta,*htau_phi,*htau_en,*htau_px,*htau_py,*htau_pz,*htau_charge,*htau_oldDM, *htau_newDM, *hZmt,*hpzetacut, *hpzeta_woCut, *hmt_woCut,*h_dilepton;

				TH1F *h_chargedIso, *h_neutral, *h_numIsoCone,*h_numIsoPhotons,*h_numIsoNeutral,*h_numCharged,*h_numPhotons,*h_numNeuHadr,*h_numCharSignal;
				TH1F *h_puCorr, *h_muonloose,*h_muontight,*h_evloose,*h_eloose,*h_emedium,*h_etight, *h_evtight;

				TH1F *htau_Loose3Hits,*htau_Medium3Hits,*htau_Tight3Hits,*htau_RawIso,*htau_oldDMwoLTraw,*htau_newDMwoLTraw,*htau_oldDMwLT,*htau_newDMwLT;

				TH1F *h_nbjets, *h_nextraElectron, *h_nextraMuon, *h_nextraElectron_after;
				TH1F *NTracks, *NTracksSignal, *NTracksIsolation;  
				*/

		ifstream file_db1; 
		char datafile[2000];               
		static const int NMAX = 6;
		static const int pair_N = 2;
		static const int  iheep_N = 2;
		static const int itauCR_N = 10;
		TH1F *h_nbjets[pair_N][iheep_N][itauCR_N];
		TH1F *hEventCounter[pair_N][iheep_N][itauCR_N], *hPV[pair_N][iheep_N][itauCR_N];    
		TH1D *h_Fill_NV[pair_N][iheep_N][itauCR_N];
		TH1F *hmu1_pt[pair_N][iheep_N][itauCR_N],  *hmu1_eta[pair_N][iheep_N][itauCR_N], *hmu1_phi[pair_N][iheep_N][itauCR_N];
		TH1F *hmu2_pt[pair_N][iheep_N][itauCR_N],  *hmu2_eta[pair_N][iheep_N][itauCR_N], *hmu2_phi[pair_N][iheep_N][itauCR_N];



		TH1F *hdimumass[pair_N][iheep_N][itauCR_N], *hdimupt[pair_N][iheep_N][itauCR_N], *hdimueta[pair_N][iheep_N][itauCR_N], *hdimuphi[pair_N][iheep_N][itauCR_N] ;

		TH1F *hmutaumass[pair_N][iheep_N][itauCR_N], *hmet[pair_N][iheep_N][itauCR_N], *h_num_jets[pair_N][iheep_N][itauCR_N], *h_num_bjets[pair_N][iheep_N][itauCR_N];

		TH1F *hLeadJetPt[pair_N][iheep_N][itauCR_N], *hSubLeadJetPt[pair_N][iheep_N][itauCR_N], *hLeadJetEta[pair_N][iheep_N][itauCR_N],*hSubLeadJetEta[pair_N][iheep_N][itauCR_N], *hLeadJetPhi[pair_N][iheep_N][itauCR_N], *hSubLeadJetPhi[pair_N][iheep_N][itauCR_N];
		TH1F *hDeltaS[pair_N][iheep_N][itauCR_N], *hDeltaPt[pair_N][iheep_N][itauCR_N],*hJetDeltaPhi[pair_N][iheep_N][itauCR_N];
		TH1D *FillCut[pair_N][iheep_N][itauCR_N][24];
		int mutaudecay, delR , mujetoverlap , taujetoverlap, mueleoverlap , taueleoverlap , mucuts , Ngoodpair , taucuts, extraElectron , mtCut , pzetA , dilepton , passDMF, extramu ,tauiso, chargeReq;
		int Nbjets,evenT,NtrigFired;
		int Nmutaudecay,NdelR, Nmujetoverlap, Ntaujetoverlap, Nmueleoverlap , Ntaueleoverlap, Nmucuts , Ntaucuts, NextraElectron, NmtCut, NpzetA, Ndilepton, NpassDMF, Nextramu, Ntauiso, NchargeReq;
		float MaxBJetPt = 670., MaxLJetPt = 1000.;
		bool isBSF= true;

		static const int tau_n_N = 7;


		TH1D          *h_Fill_Tau_N_TauFake[tau_n_N];
		TH1D          *h_Fill_Electron_N_TauFake[tau_n_N];

		TH1D          *h_EleTauCharge_TauFake[tau_n_N];
		TH1D          *h_FillmT_Tau_TauFake[tau_n_N];
		TH1D          *h_FillmT_Ele_TauFake[tau_n_N];
		TH1D          *h_Fill_DPhi_Ele_Met_TauFake[tau_n_N];
		TH1D          *h_Fill_DPhi_Tau_Met_TauFake[tau_n_N];
		TH1D          *h_Fill_DPhi_Ele_Tau_TauFake[tau_n_N];
		TH1D          *h_Fill_EleTauMass_TauFake[tau_n_N];
		TH1D          *h_Fill_TauMETMass_TauFake[tau_n_N];
		TH1D          *h_Fill_TotalMass_TauFake[tau_n_N];
		TH1D          *h_Fill_Met_TauFake[tau_n_N];
		TH1D          *h_Fill_CollMass_TauFake[tau_n_N];
		TH1D *h_Fill_PZeta_TauFake[tau_n_N];
		TH1D *h_Fill_PZetaVis_TauFake[tau_n_N];
		TH2D *h_Fill_PZeta_2D_TauFake[tau_n_N];
		TH1D *h_Fill_PZeta_Cut_TauFake[tau_n_N];
		TH1D *h_NormalizedChi2_TauFake[tau_n_N];                                                                                                                                 
		TH1D *h_CombinedQuality_chi2Local_TauFake[tau_n_N]; 
		TH1D *h_CombinedQuality_trkKink_TauFake[tau_n_N];
		TH1D *h_ValidFraction_TauFake[tau_n_N]; 
		TH1D *h_segment_compatibilty_TauFake[tau_n_N]; 
		TH1D *h_globalMuon_TauFake[tau_n_N];     
		TH1D *h_LooseMuon_TauFake[tau_n_N];  
		TH1D *h_Fill_mutaupairs_TauFake[tau_n_N];           
		TH1D *h_Fill_taupt_TauFake[tau_n_N];                           
		TH1D *h_Fill_taueta_TauFake[tau_n_N];                          
		TH1D *h_Fill_mupt_TauFake[tau_n_N];                          
		TH1D *h_Fill_mueta_TauFake[tau_n_N];                          
		TH1D *h_Fill_metphi_TauFake[tau_n_N];
		TH1D *h_Fill_mupt_SC_TauFake[tau_n_N];
		TH1D *h_Fill_mueta_SC_TauFake[tau_n_N];
		TH2D *h_Fill_mupt_mueta_TauFake[tau_n_N];
		TH2D *h_Fill_mupt_muetaSC_TauFake[tau_n_N];

		TH1D *h_Fill_tauphi_TauFake[tau_n_N];
		TH1D *h_Fill_muphi_SC_TauFake[tau_n_N];
		TH2D *h_Fill_tauphi_muphi_TauFake[tau_n_N];
		TH2D *h_Fill_tauphi_metphi_TauFake[tau_n_N];
		TH2D *h_Fill_muphi_metphi_TauFake[tau_n_N];

		///
		TH1D          *h_Fill_Tau_N_ETauFake[tau_n_N];
		TH1D          *h_Fill_Electron_N_ETauFake[tau_n_N];

		TH1D          *h_EleTauCharge_ETauFake[tau_n_N];
		TH1D          *h_FillmT_Tau_ETauFake[tau_n_N];
		TH1D          *h_FillmT_Ele_ETauFake[tau_n_N];
		TH1D          *h_Fill_DPhi_Ele_Met_ETauFake[tau_n_N];
		TH1D          *h_Fill_DPhi_Tau_Met_ETauFake[tau_n_N];
		TH1D          *h_Fill_DPhi_Ele_Tau_ETauFake[tau_n_N];
		TH1D          *h_Fill_EleTauMass_ETauFake[tau_n_N];
		TH1D          *h_Fill_TauMETMass_ETauFake[tau_n_N];
		TH1D          *h_Fill_TotalMass_ETauFake[tau_n_N];
		TH1D          *h_Fill_Met_ETauFake[tau_n_N];
		TH1D          *h_Fill_CollMass_ETauFake[tau_n_N];
		TH1D *h_Fill_PZeta_ETauFake[tau_n_N];
		TH1D *h_Fill_PZetaVis_ETauFake[tau_n_N];
		TH2D *h_Fill_PZeta_2D_ETauFake[tau_n_N];
		TH1D *h_Fill_PZeta_Cut_ETauFake[tau_n_N];
		TH1D *h_NormalizedChi2_ETauFake[tau_n_N];                                                                                                                                 
		TH1D *h_CombinedQuality_chi2Local_ETauFake[tau_n_N]; 
		TH1D *h_CombinedQuality_trkKink_ETauFake[tau_n_N];
		TH1D *h_ValidFraction_ETauFake[tau_n_N]; 
		TH1D *h_segment_compatibilty_ETauFake[tau_n_N]; 
		TH1D *h_globalMuon_ETauFake[tau_n_N];     
		TH1D *h_LooseMuon_ETauFake[tau_n_N];  
		TH1D *h_Fill_mutaupairs_ETauFake[tau_n_N];           
		TH1D *h_Fill_taupt_ETauFake[tau_n_N];                           
		TH1D *h_Fill_taueta_ETauFake[tau_n_N];                          
		TH1D *h_Fill_mupt_ETauFake[tau_n_N];                          
		TH1D *h_Fill_mueta_ETauFake[tau_n_N];                          
		TH1D *h_Fill_metphi_ETauFake[tau_n_N];
		TH1D *h_Fill_mupt_SC_ETauFake[tau_n_N];
		TH1D *h_Fill_mueta_SC_ETauFake[tau_n_N];
		TH2D *h_Fill_mupt_mueta_ETauFake[tau_n_N];
		TH2D *h_Fill_mupt_muetaSC_ETauFake[tau_n_N];

		//
		TH1D          *h_Fill_Tau_N_TauPass[tau_n_N];
		TH1D          *h_Fill_Electron_N_TauPass[tau_n_N];

		TH1D          *h_EleTauCharge_TauPass[tau_n_N];
		TH1D          *h_FillmT_Tau_TauPass[tau_n_N];
		TH1D          *h_FillmT_Ele_TauPass[tau_n_N];
		TH1D          *h_Fill_DPhi_Ele_Met_TauPass[tau_n_N];
		TH1D          *h_Fill_DPhi_Tau_Met_TauPass[tau_n_N];
		TH1D          *h_Fill_DPhi_Ele_Tau_TauPass[tau_n_N];
		TH1D          *h_Fill_EleTauMass_TauPass[tau_n_N];
		TH1D          *h_Fill_TauMETMass_TauPass[tau_n_N];
		TH1D          *h_Fill_TotalMass_TauPass[tau_n_N];
		TH1D          *h_Fill_Met_TauPass[tau_n_N];
		TH1D          *h_Fill_CollMass_TauPass[tau_n_N];
		TH1D *h_Fill_PZeta_TauPass[tau_n_N];
		TH1D *h_Fill_PZetaVis_TauPass[tau_n_N];
		TH2D *h_Fill_PZeta_2D_TauPass[tau_n_N];
		TH1D *h_Fill_PZeta_Cut_TauPass[tau_n_N];
		TH1D *h_NormalizedChi2_TauPass[tau_n_N];                                                                                                                                 
		TH1D *h_CombinedQuality_chi2Local_TauPass[tau_n_N]; 
		TH1D *h_CombinedQuality_trkKink_TauPass[tau_n_N];
		TH1D *h_ValidFraction_TauPass[tau_n_N]; 
		TH1D *h_segment_compatibilty_TauPass[tau_n_N]; 
		TH1D *h_globalMuon_TauPass[tau_n_N];     
		TH1D *h_LooseMuon_TauPass[tau_n_N];  
		TH1D *h_Fill_mutaupairs_TauPass[tau_n_N];           
		TH1D *h_Fill_taupt_TauPass[tau_n_N];                           
		TH1D *h_Fill_taueta_TauPass[tau_n_N];                          
		TH1D *h_Fill_mupt_TauPass[tau_n_N];                          
		TH1D *h_Fill_mueta_TauPass[tau_n_N];                          
		TH1D *h_Fill_metphi_TauPass[tau_n_N];
		TH1D *h_Fill_mupt_SC_TauPass[tau_n_N];
		TH1D *h_Fill_mueta_SC_TauPass[tau_n_N];
		TH2D *h_Fill_mupt_mueta_TauPass[tau_n_N];
		TH2D *h_Fill_mupt_muetaSC_TauPass[tau_n_N];

		TH1D *h_Fill_NV_TauPass[tau_n_N];
		TH1D *h_Fill_NV_TauFake[tau_n_N];
		TH1D *h_Fill_NV_TauFakeM[tau_n_N];
		TH1D *h_Fill_NV_TauFakeM_DY[tau_n_N];

		TH1D *h_Fill_NV_ETauFake[tau_n_N];

		TH1D          *h_Fill_Tau_N_TauPassM[tau_n_N];
		TH1D          *h_Fill_Electron_N_TauPassM[tau_n_N];

		TH1D          *h_EleTauCharge_TauPassM[tau_n_N];
		TH1D          *h_FillmT_Tau_TauPassM[tau_n_N];
		TH1D          *h_FillmT_Ele_TauPassM[tau_n_N];
		TH1D          *h_Fill_DPhi_Ele_Met_TauPassM[tau_n_N];
		TH1D          *h_Fill_DPhi_Tau_Met_TauPassM[tau_n_N];
		TH1D          *h_Fill_DPhi_Ele_Tau_TauPassM[tau_n_N];
		TH1D          *h_Fill_EleTauMass_TauPassM[tau_n_N];
		TH1D          *h_Fill_TauMETMass_TauPassM[tau_n_N];
		TH1D          *h_Fill_TotalMass_TauPassM[tau_n_N];
		TH1D          *h_Fill_Met_TauPassM[tau_n_N];
		TH1D          *h_Fill_CollMass_TauPassM[tau_n_N];
		TH1D *h_Fill_PZeta_TauPassM[tau_n_N];
		TH1D *h_Fill_PZetaVis_TauPassM[tau_n_N];
		TH2D *h_Fill_PZeta_2D_TauPassM[tau_n_N];
		TH1D *h_Fill_PZeta_Cut_TauPassM[tau_n_N];
		TH1D *h_NormalizedChi2_TauPassM[tau_n_N];                                                                                                                                 
		TH1D *h_CombinedQuality_chi2Local_TauPassM[tau_n_N]; 
		TH1D *h_CombinedQuality_trkKink_TauPassM[tau_n_N];
		TH1D *h_ValidFraction_TauPassM[tau_n_N]; 
		TH1D *h_segment_compatibilty_TauPassM[tau_n_N]; 
		TH1D *h_globalMuon_TauPassM[tau_n_N];     
		TH1D *h_LooseMuon_TauPassM[tau_n_N];  
		TH1D *h_Fill_mutaupairs_TauPassM[tau_n_N];           
		TH1D *h_Fill_taupt_TauPassM[tau_n_N];                           
		TH1D *h_Fill_taueta_TauPassM[tau_n_N];                          
		TH1D *h_Fill_mupt_TauPassM[tau_n_N];                          
		TH1D *h_Fill_mueta_TauPassM[tau_n_N];                          
		TH1D *h_Fill_metphi_TauPassM[tau_n_N];
		TH1D *h_Fill_mupt_SC_TauPassM[tau_n_N];
		TH1D *h_Fill_mueta_SC_TauPassM[tau_n_N];
		TH2D *h_Fill_mupt_mueta_TauPassM[tau_n_N];
		TH2D *h_Fill_mupt_muetaSC_TauPassM[tau_n_N];

		TH1D *h_Fill_NV_TauPassM[tau_n_N];


		//
		const float mu_mass = 0.10565837;                        
		TFile *fIn, *fileIn ;
		// TTree* treePtr = (TTree*) fIn->Get("HTauTauTree/HTauTauTree");
		TH1F *evCounter;
		IIHEAnalysis* tree;       
		int getTAUidNumber(TString tauIDname);
		std::string JecFileUncMC, JecFileUncData, L1Path,L2Path,L3Path,L1DATAPath,L2DATAPath,L3DATAPath,L2L3ResidualPath;


		double nEventsRaw , nEventsStored , mc_nEventsWeighted , nEventsiihe ;  
		JetCorrectorParameters *L3JetPar;
		JetCorrectorParameters *L2JetPar;
		JetCorrectorParameters *L1JetPar;

		JetCorrectorParameters *L2L3JetParDATA;
		JetCorrectorParameters *L3JetParDATA;
		JetCorrectorParameters *L2JetParDATA ;
		JetCorrectorParameters *L1JetParDATA ;

		std::vector<JetCorrectorParameters> vPar, vParData, vParDataRes;
		FactorizedJetCorrector *JetCorrector;
		FactorizedJetCorrector *JetCorrectorDATA;
		FactorizedJetCorrector *JetCorrectorDATARES;
		FactorizedJetCorrector *JetCorrectortmp;                                

		JetCorrectionUncertainty *unc ;

		JetCorrectionUncertainty *uncDATA;


		double unc_val ,ptCor_shifted;

		pair<double,double> GetElectronSF(double elept, double eleta);
		pair<double,double> GetMuonTrigger(double elept, double eleta);
		pair<double,double> GetMuonTriggerAfter(double elept, double eleta);
		double GetDYweight(double , double );                                  

		/// config parameters

		string isoTau;
		string MCHistosForPU;
		string DataHistosForPU;
		string BTagInputFile;
		string InputFile,DYWeightFile, DYSample;
		bool isdata,isSignalZ, isBkgZ;
		pair<vector<double>, vector<TLorentzVector>> ReapplyJECMC(vector <unsigned int> & indx );
		pair<vector<double>, vector<TLorentzVector>> ReapplyJECData(vector <unsigned int> & indx );
		pair<vector<double>, vector<TLorentzVector>> jecCorrAndKinematic;          
		vector <unsigned int> jecInd;      
		float qter,qter1;
		int nbjets;                                                                                                                  
		string MuonIDIsoSF ;
		double scalefactor1, scalefactor2;
		double trigwg;
		ScaleFactor * myScaleFactor ;                              
		ScaleFactor * myScaleFactor_forTrig;

		TLorentzVector FirstObj, SecondObj;
		TH1F *hmutaumasswo;
		TH1F *hPVwo;
		TH1F *hmetwo;
		TH1F *hmu1_ptwo;
		TH1F *hmu1_etawo;
		TH1F *hmu2_ptwo;
		TH1F *hmu2_etawo;
		TH1F *gen_dimuonmassCut50;
		TH1F *gen_dimuonptCut50;
		TH1F *gen_dimuonmass;
		TH1F *gen_dimuonpt;
		TH2F *ZPt_vs_DimuonMass;
		TH2F *DimuonPt_vs_DimuonMass;
		TH2F *SigDimuonPt_vs_DimuonMass;
		TH1F *Siggen_dimuonmassCut50;
		TH1F *Siggen_dimuonptCut50;
		TH1F *h_jetpt[pair_N][iheep_N][itauCR_N] ,*h_jeteta[pair_N][iheep_N][itauCR_N] ,*h_jetphi[pair_N][iheep_N][itauCR_N] , *h_jetBscore[pair_N][iheep_N][itauCR_N], *h_bjetpt[pair_N][iheep_N][itauCR_N] ,*h_bjeteta[pair_N][iheep_N][itauCR_N] ,*h_bjetphi[pair_N][iheep_N][itauCR_N];
		double mu1ptsmear, mu2ptsmear;
		TH2F *h2DPU_NVtx[pair_N][iheep_N][itauCR_N];
		float qtertmp;
		string MuonTriggerSF;
		bool dau1fired;
		bool dau2fired ;
		double GenObj2index; double GenObj1index;
		double DYSF;
		TH2F *RecoDimuonPt_vs_DimuonMass;
		double muscale;
		string filename;
		void RunSimple(const edm::Event& , const edm::EventSetup& );                     
		void FillTreeForUnfolding(const edm::Event& , const edm::EventSetup&);

		//////////////////// variables for unfolding tree //////////////////
		float recomu1pt, recomu1eta, recomu1phi, recomu2pt, recomu2eta, recomu2phi;
		float recodimupt, recodimueta, recodimuphi, recodimumass;
		float recojet1pt, recojet1eta, recojet1phi, recojet2pt, recojet2eta, recojet2phi;
		float recoDeltaS, recoDeltaPhi, recoDeltaPt;

		////////////////// gen ./////////////////////
		bool reject_event;
		float genmu1pt, genmu1eta, genmu1phi, genmu2pt, genmu2eta, genmu2phi;
		float gendimupt, gendimueta, gendimuphi, gendimumass;
		float genjet1pt, genjet1eta, genjet1phi, genjet2pt, genjet2eta, genjet2phi;

		float genDeltaS, genDeltaPhi, genDeltaPt;
		bool genZjet, recoZjet;
		double genevt;
		double pileupweight;
		double genevtwt,recoevtwt;
		vector<unsigned int> dimuindx;
		float GetEfficiency(float eta, float pt, TH2F *hist) ;
		TH2F *fhDMuMediumSF ;
		TH2F *fhDMuIsoSF;
		TH1D *h_Events_Before_Skim;
		TH1D *h_Events_After_Skim;
		TH1D *h_Events_After_GenFilter;
		TH1D *h_Fill_Mass_Gen_toChk;
		TH1D *h_Fill_Mass_Gen_toChk_before;
		double Mass;

		TH1D  *h_Fill_gsf_hadronicOverEm_E[pair_N][iheep_N][itauCR_N];
		TH1D  *h_Fill_gsf_sc_energy_E[pair_N][iheep_N][itauCR_N];
		TH2D  *h_Fill_gsf_hadronicOverEm_sc_energy_E[pair_N][iheep_N][itauCR_N];
		TH1D  *h_Fill_gsf_full5x5_sigmaIetaIeta_E[pair_N][iheep_N][itauCR_N];
		TH1D  *h_Fill_gsf_SumEt_E[pair_N][iheep_N][itauCR_N];
		TH1D  *h_Fill_TrkIso_E[pair_N][iheep_N][itauCR_N];


		TH1D  *h_Fill_gsf_hadronicOverEm_B[pair_N][iheep_N][itauCR_N];
		TH1D  *h_Fill_gsf_sc_energy_B[pair_N][iheep_N][itauCR_N];
		TH2D  *h_Fill_gsf_hadronicOverEm_sc_energy_B[pair_N][iheep_N][itauCR_N];
		TH1D  *h_Fill_gsf_full5x5_sigmaIetaIeta_B[pair_N][iheep_N][itauCR_N];
		TH1D  *h_Fill_gsf_SumEt_B[pair_N][iheep_N][itauCR_N];
		TH1D  *h_Fill_TrkIso_B[pair_N][iheep_N][itauCR_N];

		bool isOS, isSS;
		//	double TauPtBinSF_CR5(double pt);
		bool PassPreSelections(int iele, TLorentzVector ele);
		double ElePtBinSF_CR5(double elept, double eleSCeta);
		unsigned int goodHLTElecIndex1 ;
		bool matchedToGenObjetcs(TLorentzVector tau, unsigned int &);
		double MutoTauFR(TLorentzVector tau);
		TH1D          *h_Fill_Tau_N_TauFakeM[tau_n_N];
		TH1D          *h_Fill_Electron_N_TauFakeM[tau_n_N];

		TH1D          *h_EleTauCharge_TauFakeM[tau_n_N];
		TH1D          *h_FillmT_Tau_TauFakeM[tau_n_N];
		TH1D          *h_FillmT_Ele_TauFakeM[tau_n_N];
		TH1D          *h_Fill_DPhi_Ele_Met_TauFakeM[tau_n_N];
		TH1D          *h_Fill_DPhi_Tau_Met_TauFakeM[tau_n_N];
		TH1D          *h_Fill_DPhi_Ele_Tau_TauFakeM[tau_n_N];
		TH1D          *h_Fill_EleTauMass_TauFakeM[tau_n_N];
		TH1D          *h_Fill_TauMETMass_TauFakeM[tau_n_N];
		TH1D          *h_Fill_TotalMass_TauFakeM[tau_n_N];
		TH1D          *h_Fill_Met_TauFakeM[tau_n_N];
		TH1D          *h_Fill_CollMass_TauFakeM[tau_n_N];
		TH1D *h_Fill_PZeta_TauFakeM[tau_n_N];
		TH1D *h_Fill_PZetaVis_TauFakeM[tau_n_N];
		TH2D *h_Fill_PZeta_2D_TauFakeM[tau_n_N];
		TH1D *h_Fill_PZeta_Cut_TauFakeM[tau_n_N];
		TH1D *h_NormalizedChi2_TauFakeM[tau_n_N];                                                                                                                                 
		TH1D *h_CombinedQuality_chi2Local_TauFakeM[tau_n_N]; 
		TH1D *h_CombinedQuality_trkKink_TauFakeM[tau_n_N];
		TH1D *h_ValidFraction_TauFakeM[tau_n_N]; 
		TH1D *h_segment_compatibilty_TauFakeM[tau_n_N]; 
		TH1D *h_globalMuon_TauFakeM[tau_n_N];     
		TH1D *h_LooseMuon_TauFakeM[tau_n_N];  
		TH1D *h_Fill_mutaupairs_TauFakeM[tau_n_N];           
		TH1D *h_Fill_taupt_TauFakeM[tau_n_N];                           
		TH1D *h_Fill_taueta_TauFakeM[tau_n_N];                          
		TH1D *h_Fill_mupt_TauFakeM[tau_n_N];                          
		TH1D *h_Fill_mueta_TauFakeM[tau_n_N];                          
		TH1D *h_Fill_metphi_TauFakeM[tau_n_N];
		TH1D *h_Fill_mupt_SC_TauFakeM[tau_n_N];
		TH1D *h_Fill_mueta_SC_TauFakeM[tau_n_N];
		TH2D *h_Fill_mupt_mueta_TauFakeM[tau_n_N];
		TH2D *h_Fill_mupt_muetaSC_TauFakeM[tau_n_N];


		int GenTaus();
		TH1D *h_Count_Taus;
		TH1D *h_Fill_bjet_pt_TauFake[tau_n_N];
		TH1D *h_Fill_bjet_eta_TauFake[tau_n_N];
		TH1D *h_Fill_bjet_disc_TauFake[tau_n_N];

		TH1D *h_Fill_bjet_pt_TauFakeM[tau_n_N];
		TH1D *h_Fill_bjet_eta_TauFakeM[tau_n_N];
		TH1D *h_Fill_bjet_disc_TauFakeM[tau_n_N];


		////

		TH1D          *h_Fill_Tau_N_TauFakeM_DY[tau_n_N];
		TH1D          *h_Fill_Electron_N_TauFakeM_DY[tau_n_N];

		TH1D          *h_EleTauCharge_TauFakeM_DY[tau_n_N];
		TH1D          *h_FillmT_Tau_TauFakeM_DY[tau_n_N];
		TH1D          *h_FillmT_Ele_TauFakeM_DY[tau_n_N];
		TH1D          *h_Fill_DPhi_Ele_Met_TauFakeM_DY[tau_n_N];
		TH1D          *h_Fill_DPhi_Tau_Met_TauFakeM_DY[tau_n_N];
		TH1D          *h_Fill_DPhi_Ele_Tau_TauFakeM_DY[tau_n_N];
		TH1D          *h_Fill_EleTauMass_TauFakeM_DY[tau_n_N];
		TH1D          *h_Fill_TauMETMass_TauFakeM_DY[tau_n_N];
		TH1D          *h_Fill_TotalMass_TauFakeM_DY[tau_n_N];
		TH1D          *h_Fill_Met_TauFakeM_DY[tau_n_N];
		TH1D          *h_Fill_CollMass_TauFakeM_DY[tau_n_N];
		TH1D *h_Fill_PZeta_TauFakeM_DY[tau_n_N];
		TH1D *h_Fill_PZetaVis_TauFakeM_DY[tau_n_N];
		TH2D *h_Fill_PZeta_2D_TauFakeM_DY[tau_n_N];
		TH1D *h_Fill_PZeta_Cut_TauFakeM_DY[tau_n_N];
		TH1D *h_NormalizedChi2_TauFakeM_DY[tau_n_N];                                                                                                                                 
		TH1D *h_CombinedQuality_chi2Local_TauFakeM_DY[tau_n_N]; 
		TH1D *h_CombinedQuality_trkKink_TauFakeM_DY[tau_n_N];
		TH1D *h_ValidFraction_TauFakeM_DY[tau_n_N]; 
		TH1D *h_segment_compatibilty_TauFakeM_DY[tau_n_N]; 
		TH1D *h_globalMuon_TauFakeM_DY[tau_n_N];     
		TH1D *h_LooseMuon_TauFakeM_DY[tau_n_N];  
		TH1D *h_Fill_mutaupairs_TauFakeM_DY[tau_n_N];           
		TH1D *h_Fill_taupt_TauFakeM_DY[tau_n_N];                           
		TH1D *h_Fill_taueta_TauFakeM_DY[tau_n_N];                          
		TH1D *h_Fill_mupt_TauFakeM_DY[tau_n_N];                          
		TH1D *h_Fill_mueta_TauFakeM_DY[tau_n_N];                          
		TH1D *h_Fill_metphi_TauFakeM_DY[tau_n_N];
		TH1D *h_Fill_mupt_SC_TauFakeM_DY[tau_n_N];
		TH1D *h_Fill_mueta_SC_TauFakeM_DY[tau_n_N];
		TH2D *h_Fill_mupt_mueta_TauFakeM_DY[tau_n_N];
		TH2D *h_Fill_mupt_muetaSC_TauFakeM_DY[tau_n_N];

		TH1D *h_Fill_bjet_pt_TauFakeM_DY[tau_n_N];
		TH1D *h_Fill_bjet_eta_TauFakeM_DY[tau_n_N];
		TH1D *h_Fill_bjet_disc_TauFakeM_DY[tau_n_N];

		TH1D *h_Fill_bjet_pt_TauPass[tau_n_N];
		TH1D *h_Fill_bjet_eta_TauPass[tau_n_N];
		TH1D *h_Fill_bjet_disc_TauPass[tau_n_N];

		TH1D *h_Fill_bjet_pt_TauPassM[tau_n_N];
		TH1D *h_Fill_bjet_eta_TauPassM[tau_n_N];
		TH1D *h_Fill_bjet_disc_TauPassM[tau_n_N];

		TH1D *h_Fill_bjet_pt_ETauFake[tau_n_N];
		TH1D *h_Fill_bjet_eta_ETauFake[tau_n_N];
		TH1D *h_Fill_bjet_disc_ETauFake[tau_n_N];




		TH1D          *h_Fill_Tau_N_TauPassM_OS[tau_n_N];
		TH1D          *h_Fill_Electron_N_TauPassM_OS[tau_n_N];

		TH1D          *h_EleTauCharge_TauPassM_OS[tau_n_N];
		TH1D          *h_FillmT_Tau_TauPassM_OS[tau_n_N];
		TH1D          *h_FillmT_Ele_TauPassM_OS[tau_n_N];
		TH1D          *h_Fill_DPhi_Ele_Met_TauPassM_OS[tau_n_N];
		TH1D          *h_Fill_DPhi_Tau_Met_TauPassM_OS[tau_n_N];
		TH1D          *h_Fill_DPhi_Ele_Tau_TauPassM_OS[tau_n_N];
		TH1D          *h_Fill_EleTauMass_TauPassM_OS[tau_n_N];
		TH1D          *h_Fill_TauMETMass_TauPassM_OS[tau_n_N];
		TH1D          *h_Fill_TotalMass_TauPassM_OS[tau_n_N];
		TH1D          *h_Fill_Met_TauPassM_OS[tau_n_N];
		TH1D          *h_Fill_CollMass_TauPassM_OS[tau_n_N];
		TH1D *h_Fill_PZeta_TauPassM_OS[tau_n_N];
		TH1D *h_Fill_PZetaVis_TauPassM_OS[tau_n_N];
		TH2D *h_Fill_PZeta_2D_TauPassM_OS[tau_n_N];
		TH1D *h_Fill_PZeta_Cut_TauPassM_OS[tau_n_N];
		TH1D *h_NormalizedChi2_TauPassM_OS[tau_n_N];                                                                                                                                 
		TH1D *h_CombinedQuality_chi2Local_TauPassM_OS[tau_n_N]; 
		TH1D *h_CombinedQuality_trkKink_TauPassM_OS[tau_n_N];
		TH1D *h_ValidFraction_TauPassM_OS[tau_n_N]; 
		TH1D *h_segment_compatibilty_TauPassM_OS[tau_n_N]; 
		TH1D *h_globalMuon_TauPassM_OS[tau_n_N];     
		TH1D *h_LooseMuon_TauPassM_OS[tau_n_N];  
		TH1D *h_Fill_mutaupairs_TauPassM_OS[tau_n_N];           
		TH1D *h_Fill_taupt_TauPassM_OS[tau_n_N];                           
		TH1D *h_Fill_taueta_TauPassM_OS[tau_n_N];                          
		TH1D *h_Fill_mupt_TauPassM_OS[tau_n_N];                          
		TH1D *h_Fill_mueta_TauPassM_OS[tau_n_N];                          
		TH1D *h_Fill_metphi_TauPassM_OS[tau_n_N];
		TH1D *h_Fill_mupt_SC_TauPassM_OS[tau_n_N];
		TH1D *h_Fill_mueta_SC_TauPassM_OS[tau_n_N];
		TH2D *h_Fill_mupt_mueta_TauPassM_OS[tau_n_N];
		TH2D *h_Fill_mupt_muetaSC_TauPassM_OS[tau_n_N];

		TH1D *h_Fill_NV_TauPassM_OS[tau_n_N];
		TH1D *h_Fill_bjet_pt_TauPassM_OS[tau_n_N];
		TH1D *h_Fill_bjet_eta_TauPassM_OS[tau_n_N];
		TH1D *h_Fill_bjet_disc_TauPassM_OS[tau_n_N];


		///////


		TH1D          *h_Fill_Tau_N_TauPassM_SS[tau_n_N];
		TH1D          *h_Fill_Electron_N_TauPassM_SS[tau_n_N];

		TH1D          *h_EleTauCharge_TauPassM_SS[tau_n_N];
		TH1D          *h_FillmT_Tau_TauPassM_SS[tau_n_N];
		TH1D          *h_FillmT_Ele_TauPassM_SS[tau_n_N];
		TH1D          *h_Fill_DPhi_Ele_Met_TauPassM_SS[tau_n_N];
		TH1D          *h_Fill_DPhi_Tau_Met_TauPassM_SS[tau_n_N];
		TH1D          *h_Fill_DPhi_Ele_Tau_TauPassM_SS[tau_n_N];
		TH1D          *h_Fill_EleTauMass_TauPassM_SS[tau_n_N];
		TH1D          *h_Fill_TauMETMass_TauPassM_SS[tau_n_N];
		TH1D          *h_Fill_TotalMass_TauPassM_SS[tau_n_N];
		TH1D          *h_Fill_Met_TauPassM_SS[tau_n_N];
		TH1D          *h_Fill_CollMass_TauPassM_SS[tau_n_N];
		TH1D *h_Fill_PZeta_TauPassM_SS[tau_n_N];
		TH1D *h_Fill_PZetaVis_TauPassM_SS[tau_n_N];
		TH2D *h_Fill_PZeta_2D_TauPassM_SS[tau_n_N];
		TH1D *h_Fill_PZeta_Cut_TauPassM_SS[tau_n_N];
		TH1D *h_NormalizedChi2_TauPassM_SS[tau_n_N];                                                                                                                                 
		TH1D *h_CombinedQuality_chi2Local_TauPassM_SS[tau_n_N]; 
		TH1D *h_CombinedQuality_trkKink_TauPassM_SS[tau_n_N];
		TH1D *h_ValidFraction_TauPassM_SS[tau_n_N]; 
		TH1D *h_segment_compatibilty_TauPassM_SS[tau_n_N]; 
		TH1D *h_globalMuon_TauPassM_SS[tau_n_N];     
		TH1D *h_LooseMuon_TauPassM_SS[tau_n_N];  
		TH1D *h_Fill_mutaupairs_TauPassM_SS[tau_n_N];           
		TH1D *h_Fill_taupt_TauPassM_SS[tau_n_N];                           
		TH1D *h_Fill_taueta_TauPassM_SS[tau_n_N];                          
		TH1D *h_Fill_mupt_TauPassM_SS[tau_n_N];                          
		TH1D *h_Fill_mueta_TauPassM_SS[tau_n_N];                          
		TH1D *h_Fill_metphi_TauPassM_SS[tau_n_N];
		TH1D *h_Fill_mupt_SC_TauPassM_SS[tau_n_N];
		TH1D *h_Fill_mueta_SC_TauPassM_SS[tau_n_N];
		TH2D *h_Fill_mupt_mueta_TauPassM_SS[tau_n_N];
		TH2D *h_Fill_mupt_muetaSC_TauPassM_SS[tau_n_N];

		TH1D *h_Fill_NV_TauPassM_SS[tau_n_N];
		TH1D *h_Fill_bjet_pt_TauPassM_SS[tau_n_N];
		TH1D *h_Fill_bjet_eta_TauPassM_SS[tau_n_N];
		TH1D *h_Fill_bjet_disc_TauPassM_SS[tau_n_N];

		double TauPtBinSF_CR5(double taupt, double abs_taueta, double decaymode) ;
		vector<float> gen_tau_had_pt, gen_tau_had_eta, gen_tau_had_phi, gen_tau_had_energy, gen_tau_had_charge;
		double FakeRate(double taupt, double jetpt) ;
		int pairs;
		int first_index;
		int second_index;

		vector<TH1F*> hh[10][4][3];

		TH1F *h_elept;
		TH1F *h_eleeta;
		TH1F *h_elephi;
		TH1F *h_mupt;
		TH1F *h_mueta;
		TH1F *h_muphi;
		TH1F *h_emumass;
		TH1F *h_met;

		TH1F *h_taupt_num;
		TH1F *h_taueta_num;
		TH1F *h_tauphi_num;
		TH1F *h_jetpt_num;
		TH1F *h_jeteta_num;
		TH1F *h_jetphi_num;
		TH1F *h_taupt_den;
		TH1F *h_taueta_den;
		TH1F *h_tauphi_den;
		TH1F *h_jetpt_den;
		TH1F *h_jeteta_den;
		TH1F *h_jetphi_den;
		TH1F *h_taupt_failed;
		TH1F *h_taueta_failed;
		TH1F *h_tauphi_failed;
		TH1F *h_jetpt_failed;
		TH1F *h_jeteta_failed;
		TH1F *h_jetphi_failed;
		double FakeRate_SSMtLow(double taupt, double jetpt, TH2F *hist);	
		double FakeRate_mumu(double taupt, double jetpt);
		bool TauMatchedToJet(TLorentzVector taup4, unsigned int &matched_jet_indx) ;
		string FakeRateSSfile, FakeRateDYJetsfile;

		TH2F *feta_barrel_taupt_300_jetpt_150;
		TH2F *feta_endcap_taupt_300_jetpt_150;
		TH2F *feta_barrel_taupt_300_jetpt_300;
		TH2F *feta_endcap_taupt_300_jetpt_300;
		TH2F *feta_barrel_taupt_300_jetpt_1000;
		TH2F *feta_endcap_taupt_300_jetpt_1000;
		TH2F *feta_barrel_taupt_1000_jetpt_300;
		TH2F *feta_endcap_taupt_1000_jetpt_300;
		TH2F *feta_barrel_taupt_1000_jetpt_1000;
		TH2F *feta_endcap_taupt_1000_jetpt_1000;


		TFile *fake_file_DY;
		TFile *fake_file;

		TH2F *DY_feta_barrel_taupt_300_jetpt_150;
		TH2F *DY_feta_endcap_taupt_300_jetpt_150;
		TH2F *DY_feta_barrel_taupt_300_jetpt_300;
		TH2F *DY_feta_endcap_taupt_300_jetpt_300;
		TH2F *DY_feta_barrel_taupt_300_jetpt_1000;
		TH2F *DY_feta_endcap_taupt_300_jetpt_1000;
		TH2F *DY_feta_barrel_taupt_1000_jetpt_300;
		TH2F *DY_feta_endcap_taupt_1000_jetpt_300;
		TH2F *DY_feta_barrel_taupt_1000_jetpt_1000;
		TH2F *DY_feta_endcap_taupt_1000_jetpt_1000;

		bool LHEinfoMass(double &Mass) ;
		bool LHEinfoPt(double &Mass) ;

		int first_index_PP ;
		int second_index_PP ; 
		int first_index_MPP ;
		int second_index_MPP ;
		int first_index_MPP_OS ;
		int second_index_MPP_OS ;
		int first_index_MPP_SS ;
		int second_index_MPP_SS ;
		int first_index_FF;
		int second_index_FF ;
		int first_index_M;
		int second_index_M;
		double fr_weight_ss;
		int jet_index_FF;
};


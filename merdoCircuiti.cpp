#include "TGraphErrors.h"
#include "TF1.h"
#include "cmath"
#include "TFile.h"
#include <RtypesCore.h>
#include <fstream>

// Measured values

Double_t VPP {5.};
Double_t R {1.499E3};
Double_t RErr {2};
Double_t L {10.11E-3};
Double_t LErr {10E-5};
Double_t LR {38.38};
Double_t LRErr {0.03};
Double_t C {21.91E-9};
Double_t CErr {22E-11};

// Data processing functions

void correctPhiFreqResp(TGraphErrors* phiFreqResp, const TF1* phaseOffsetF) {
	for (Int_t i {0}; i < phiFreqResp->GetN(); ++i) {
		phiFreqResp->SetPoint(i, phiFreqResp->GetPointX(i), phiFreqResp->GetPointY(i) - phaseOffsetF->Eval(phiFreqResp->GetPointX(i)));
		phiFreqResp->SetPointError(i, phiFreqResp->GetErrorX(i), phaseOffsetF-);
		// need to add error on changed y value
		if (phiFreqResp->GetPointY(i) < -M_PI/2.) {
			phiFreqResp->SetPointY(i, phiFreqResp->GetPointY(i) + 2*M_PI);
		}
	}
}

void correctSins(TGraphErrors* sinusoid, Double_t deltaT, Int_t nChannel) {
	for (Int_t i {0}; i < sinusoid->GetN(); ++i) {
		sinusoid->SetPointX(i, sinusoid->GetPointX(i) + deltaT * nChannel);
	}
}

// This one provides the error, but need fitted lines to return it correctly

void addYErr(TGraphErrors* graph, TF1* yErrF) {
	for (Int_t i {0}; i < graph->GetN(); ++i) {
		graph->SetPointError(
			// xErrF->Eval(graph->GetPointX(i))
			0., // half resolution?
			yErrF->Eval(graph->GetPointX(i))
		);
	}
}

// transform two sinusoids into a lissajous shape

TGraphErrors* correlate(TGraphErrors* sin1, TGraphErrors* sin2) {
	return new TGraphErrors(sin1->GetN(), sin1->GetY(), sin2->GetY(), sin1->GetEX(), sin2->GetEY());
}

// auxiliary polarization functions

Double_t rFromPt(Double_t x, Double_t y) {
	return sqrt(x*x + y*y);
}

Double_t rErr(Double_t x, Double_t y, Double_t xErr, Double_t yErr) {
	return (std::abs(x) * xErr + std::abs(y) * yErr) / rFromPt(x, y);
}

Double_t thetaFromPt(Double_t x, Double_t y) {
	Double_t r { rFromPt(x, y) };
	Double_t cosTheta { x / r };
	Double_t sinTheta { y / r };
	Double_t tanTheta { y / x };

	bool useCos { false };
	bool useSin { false };
	bool useTan { false };

	if (std::abs(tanTheta) < std::abs(cosTheta)) {
		if (std::abs(tanTheta) < std::abs(sinTheta))
			useTan = true;
		else
			useSin = true;
	} else {
		if (std::abs(cosTheta) < std::abs(sinTheta))
			useCos = true;
		else
		  useSin = true;
	}

	if (useCos) {
		if (y >= 0.)
			return acos(cosTheta);
		else
		 	return 2*M_PI - acos(cosTheta);
	}

	if (useSin) {
		if (x >= 0.) {
			if (sinTheta >= 0.)
				return asin(sinTheta);
			else
			 	return asin(sinTheta) + 2*M_PI;
		} else {
			return M_PI - asin(sinTheta);
		}
	}

	if (useTan) {
		if (x >= 0.) {
			if (tanTheta >= 0.)
				return atan(tanTheta);
			else
			 	return atan(tanTheta) + 2*M_PI;
		} else
			return M_PI + atan(tanTheta);
	}
	
	throw std::runtime_error("Dude wtf");

}

Double_t thetaErr(Double_t x, Double_t y, Double_t xErr, Double_t yErr) {
	Double_t r { rFromPt(x, y) };
	Double_t cosTheta { x / r };
	Double_t sinTheta { y / r };
	Double_t tanTheta { y / x };

	bool useCos { false };
	bool useSin { false };
	bool useTan { false };

	if (std::abs(tanTheta) < std::abs(cosTheta)) {
		if (std::abs(tanTheta) < std::abs(sinTheta))
			useTan = true;
		else
			useSin = true;
	} else {
		if (std::abs(cosTheta) < std::abs(sinTheta))
			useCos = true;
		else
		  useSin = true;
	}

	Double_t deltaR { rErr(x, y, xErr, yErr) };

	if (useCos) {
		Double_t deltaCos { xErr/std::abs(r) + std::abs(x/(r*r)) * deltaR };
		return deltaCos/sqrt(1 - cosTheta*cosTheta);
	}
 
	if (useSin) {
		Double_t deltaSin { yErr/std::abs(r) + std::abs(y/(r*r)) * deltaR };
		return deltaSin/sqrt(1 - sinTheta*sinTheta);
	}

	if (useTan) {
		Double_t deltaTan { std::abs(y)/(x*x) * xErr + yErr / std::abs(x) };
		return deltaTan/(1 + tanTheta*tanTheta);
	}
	
	throw std::runtime_error("Dude wtf");

}

TGraphErrors* polarize(TGraphErrors* lissajous) {
	TGraphErrors* polarized = new TGraphErrors(lissajous->GetN());
	for (int i {0}; i < lissajous->GetN(); ++i) {
		polarized->SetPoint(i,
			rFromPt(lissajous->GetPointX(i), lissajous->GetPointY(i)),
			thetaFromPt(lissajous->GetPointX(i), lissajous->GetPointY(i))
		);
		polarized->SetPointError(i,
			rErr(lissajous->GetPointX(i), lissajous->GetPointY(i), lissajous->GetErrorX(i), lissajous->GetErrorY(i)),
			thetaErr(lissajous->GetPointX(i), lissajous->GetPointY(i), lissajous->GetErrorX(i), lissajous->GetErrorY(i))
		);
	}
	return polarized;
}

// Fitting functions

Double_t amplitudeFreqRespR(Double_t omega, Double_t* pars) {
	// [0] is V0
	// [1] is R
	// [2] is L
	// [3] is C
	return pars[0]*pars[1] /
		sqrt(
			pars[1]*pars[1] +
			pow((omega*pars[2] - 1/(omega * pars[3])), 2.)
		);
}

Double_t amplitudeFreqRespL(Double_t omega, Double_t* pars) {
	// [0] is V0
	// [1] is R
	// [2] is L
	// [3] is C
	return omega*pars[0]*pars[2] /
		sqrt(
			pars[1]*pars[1] +
			pow((omega*pars[2] - 1/(omega * pars[3])), 2.)
		);
}

Double_t amplitudeFreqRespC(Double_t omega, Double_t* pars) {
	// [0] is V0
	// [1] is R
	// [2] is L
	// [3] is C
	return pars[0] /
		(omega*pars[3]) /
		sqrt(
			pars[1]*pars[1] +
			pow((omega*pars[2] - 1/(omega * pars[3])), 2.)
		);
}

Double_t phaseFreqResp(Double_t omega, Double_t* pars) {
	return pars[0] + atan( (1 - omega*omega*pars[2]*pars[3]) / (omega*pars[1]*pars[3]) );
}

Double_t polarLissajous(Double_t* theta, Double_t* pars) {
	// [0] is V0
	// [1] is a 
	// [2] is phi
	
	Double_t V0 { pars[0] };
	Double_t a { pars[1] };
	Double_t phi { pars[2] };

	Double_t tanTau { 1 / tan(phi) - tan(*theta) / (a*sin(phi)) };
	Double_t cosPhi { cos(phi) };
	Double_t sinPhi { sin(phi) };
	
	return V0 * sqrt( (1 + a*a*cosPhi*cosPhi + a*a*tanTau*tanTau*sinPhi*sinPhi - 2*a*a*cosPhi*sinPhi*tanTau) / (1 + tanTau*tanTau));
}

// Data processor

void analyze(std::string dataDir) {
	
	std::vector<std::string> fileNames {
		"ampFreqRespR.txt",
		"ampFreqRespL.txt",
		"ampFreqRespC.txt",
		"phaseFreqRespR.txt",
		"phaseFreqRespL.txt",
		"phaseFreqRespC.txt",
		"ampErrCh",
		"phiErrCh",
		"sineGen",
		"sineR",
		"sineL",
		"sineC"
	};

	TGraphErrors* ampFreqRespR = new TGraphErrors((dataDir + fileNames[0]).c_str());
	TGraphErrors* ampFreqRespL = new TGraphErrors((dataDir + fileNames[1]).c_str());
	TGraphErrors* ampFreqRespC = new TGraphErrors((dataDir + fileNames[2]).c_str());

	TGraphErrors* phaseFreqRespR = new TGraphErrors((dataDir + fileNames[3]).c_str());
	TGraphErrors* phaseFreqRespL = new TGraphErrors((dataDir + fileNames[4]).c_str());
	TGraphErrors* phaseFreqRespC = new TGraphErrors((dataDir + fileNames[5]).c_str());

	// acquire amplitude errors graphs
	
	TGraphErrors* ampErrsCh[4] {};

	for (int j {0}; j < 3; ++j) {
		while (
			std::ifstream(dataDir + fileNames[6] + std::to_string(j) + ".txt").good()
		) {
			ampErrsCh[j] = new TGraphErrors(
				(dataDir + fileNames[6] + std::to_string(j) + ".txt" ).c_str()
			);
		}
	}

	/*
	// acquire phase errors graphs

	TGraphErrors* phiErrsCh[4] {};

	for (int j {0}; j < 4; ++j) {
		while (
			std::ifstream(dataDir + fileNames[7] + std::to_string(j) + ".txt").good()
		) {
			phiErrsCh[j] = new TGraphErrors(
				(dataDir + fileNames[7] + std::to_string(j) + ".txt" ).c_str()
			);
		}
	}
	*/

	// acquire phase errors graphs

	TGraphErrors* phiErrsCh[4] {};

	for (int j {0}; j < 4; ++j) {
		while (
			std::ifstream(dataDir + fileNames[7] + std::to_string(j) + ".txt").good()
		) {
			phiErrsCh[j] = new TGraphErrors(
				(dataDir + fileNames[7] + std::to_string(j) + ".txt" ).c_str()
			);
		}
	}

	// acquire phase offsets graphs
	TGraphErrors* phiOffsCh[4] {};

	for (int j {0}; j < 4; ++j) {
		while (
			std::ifstream(dataDir + fileNames[7] + std::to_string(j) + ".txt").good()
		) {
			phiOffsCh[j] = new TGraphErrors(
				(dataDir + fileNames[7] + std::to_string(j) + ".txt" ).c_str(),
				"%lg%*lg%lg"
			);
		}
	}


	// acquire sinusoids graphs

	std::vector<TGraphErrors*> sineGen {};
	std::vector<TGraphErrors*> sineR {};
	std::vector<TGraphErrors*> sineL {};
	std::vector<TGraphErrors*> sineC {};

	std::vector<TGraphErrors*> sines[4] {
		sineGen, sineR, sineL, sineC
	};

	Int_t i {0};

	for (int j {0}; j < 4; ++j) {
		i = 1;
		while (
			std::ifstream(dataDir + fileNames[8 + j] + std::to_string(i) + ".txt").good()
		) {
			sines[j].emplace_back(new TGraphErrors(
				(dataDir + fileNames[8 + j] + std::to_string(i) + ".txt" ).c_str()
				)
			);
			++i;
		}	
	}

	// Start data prep

	// Phase offset functions Fitting
	
	TF1* phiFreqRespCorrectorCh[4] {};
	
	i = 0;
	for (; i < 4; ++i) {
		phiFreqRespCorrectorCh[i] = new TF1(("correctPhiFreqRespCh" + std::to_string(i)).c_str(), "pol1", 0., 1E6);
		phiOffsCh[i]->Fit(phiFreqRespCorrectorCh[i]);
	}

	// correcting phase freq resp graph for syst phase offset
	
	correctPhiFreqResp(phaseFreqRespR, phiFreqRespCorrectorCh[1]);
	correctPhiFreqResp(phaseFreqRespL, phiFreqRespCorrectorCh[2]);
	correctPhiFreqResp(phaseFreqRespC, phiFreqRespCorrectorCh[3]);

	/*
	prob not gonna do it, hopefully negligible
	i = 0;
	for (std::vector<TGraphErrors*> sinusSet: sines) {
		for (TGraphErrors* sinusoid: sinusSet) {
			correctSins(sinusoid, AITimer, i);
			++i;
		}
	}
	*/

	// fit stdDev functions

	TF1* amplitudeErrFCh[4] {};
	TF1* phaseErrFCh[4] {};

	i = 0;
	for (; i < 4; ++i) {
		amplitudeErrFCh[i] = new TF1(("amplitudeErrFCh" + std::to_string(i)).c_str(), "pol1", 0., 1E6);
		ampErrsCh[i]->Fit(amplitudeErrFCh[i]);
		phaseErrFCh[i] = new TF1(("phaseErrFCh" + std::to_string(i)).c_str(), "pol1", 0., 1E6);
		phiErrsCh[i]->Fit(phaseErrFCh[i]);
	}

	addYErr(ampFreqRespR, amplitudeErrFCh[1]);
	addYErr(ampFreqRespL, amplitudeErrFCh[2]);
	addYErr(ampFreqRespC, amplitudeErrFCh[3]);

	addYErr(phaseFreqRespR, phaseErrFCh[1]);
	addYErr(phaseFreqRespL, phaseErrFCh[2]);
	addYErr(phaseFreqRespC, phaseErrFCh[3]);

	// coexpress sinusoids
	
	TGraphErrors* lissajousGenR1 = correlate(sines[0][0], sines[1][0]);
	TGraphErrors* lissajousGenR2 = correlate(sines[0][1], sines[1][1]);
	TGraphErrors* lissajousGenR3 = correlate(sines[0][2], sines[1][2]);

	TGraphErrors* polarLissajousGenR1 = polarize(lissajousGenR1);
	TGraphErrors* polarLissajousGenR2 = polarize(lissajousGenR2);
	TGraphErrors* polarLissajousGenR3 = polarize(lissajousGenR2);

	// fit functions, need to initialize values

	TF1* ampFreqRespRF = new TF1("ampFreqRespRF", amplitudeFreqRespR, 0., 1E6, 4);
	TF1* ampFreqRespLF = new TF1("ampFreqRespLF", amplitudeFreqRespL, 0., 1E6, 4);
	TF1* ampFreqRespCF = new TF1("ampFreqRespCF", amplitudeFreqRespC, 0., 1E6, 4);

	TF1* phaseFreqRespRF = new TF1("phaseFreqRespRF", phaseFreqResp, 0., 1E6, 4);
	TF1* phaseFreqRespLF = new TF1("phaseFreqRespLF", phaseFreqResp, 0., 1E6, 4);
	TF1* phaseFreqRespCF = new TF1("phaseFreqRespCF", phaseFreqResp, 0., 1E6, 4);

	TF1* polarLissajousGenR1F = new TF1("polarLissajousGenR1F", polarLissajous, 0., 2*M_PI, 3);
	TF1* polarLissajousGenR2F = new TF1("polarLissajousGenR2F", polarLissajous, 0., 2*M_PI, 3);
	TF1* polarLissajousGenR3F = new TF1("polarLissajousGenR3F", polarLissajous, 0., 2*M_PI, 3);

	// fitting
	
	ampFreqRespR->Fit(ampFreqRespRF);
	ampFreqRespL->Fit(ampFreqRespLF);
	ampFreqRespC->Fit(ampFreqRespCF);

	phaseFreqRespR->Fit(phaseFreqRespRF);
	phaseFreqRespL->Fit(phaseFreqRespLF);
	phaseFreqRespC->Fit(phaseFreqRespCF);

	polarLissajousGenR1->Fit(polarLissajousGenR1F);
	polarLissajousGenR2->Fit(polarLissajousGenR2F);
	polarLissajousGenR3->Fit(polarLissajousGenR3F);

	// space for cosmetic editing

	// saving EVERYTHING to TFile
	
	TFile* results = new TFile("analyzedData.ROOT", "RECREATE");

	ampFreqRespR->Write();
	ampFreqRespRF->Write();
	ampFreqRespL->Write();
	ampFreqRespLF->Write();
	ampFreqRespC->Write();
	ampFreqRespCF->Write();

	phaseFreqRespR->Write();
	ampFreqRespRF->Write();
	phaseFreqRespL->Write();
	ampFreqRespLF->Write();
	phaseFreqRespC->Write();
	ampFreqRespCF->Write();

	lissajousGenR1->Write();
	polarLissajousGenR1->Write();
	polarLissajousGenR1F->Write();
	lissajousGenR2->Write();
	polarLissajousGenR2->Write();
	polarLissajousGenR2F->Write();
	lissajousGenR3->Write();
	polarLissajousGenR3->Write();
	polarLissajousGenR3F->Write();

}

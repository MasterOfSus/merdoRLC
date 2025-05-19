#include "TGraphErrors.h"
#include "TF1.h"
#include "cmath"
#include "TFile.h"
#include <RtypesCore.h>
#include <TMath.h>
#include <fstream>

// Measured values

Double_t VPP {5.};
Double_t R {1.499E3};
Double_t RErr {2.};
Double_t L {10.11E-3};
Double_t LErr {10E-5};
Double_t LR {38.38};
Double_t LRErr {3E-2};
Double_t C {21.91E-9};
Double_t CErr {22E-11};

// Data processing functions

void correctPhiFreqResp(TGraphErrors* phiFreqResp, const TF1* phaseOffsetF) {
	for (Int_t i {0}; i < phiFreqResp->GetN(); ++i) {
		if (phiFreqResp->GetPointY(i) < -M_PI/2.) {
			phiFreqResp->SetPointY(i, phiFreqResp->GetPointY(i) + 2*M_PI);
		}
	phiFreqResp->SetPoint(i, phiFreqResp->GetPointX(i), phiFreqResp->GetPointY(i) - phaseOffsetF->Eval(phiFreqResp->GetPointX(i)));
		phiFreqResp->SetPointError(i, phiFreqResp->GetErrorX(i), 
			sqrt(
				pow((1 - phaseOffsetF->GetParameter(1))*phiFreqResp->GetErrorY(i), 2.) +
				pow(phiFreqResp->GetPointY(i)*phaseOffsetF->GetParError(1), 2.) + 
				pow(phaseOffsetF->GetParameter(0), 2.)
			)
		);
	}
}

void correctSins(TGraphErrors* sinusoid, Double_t deltaT, Int_t nChannel) {
	for (Int_t i {0}; i < sinusoid->GetN(); ++i) {
		sinusoid->SetPointX(i, sinusoid->GetPointX(i) + deltaT * nChannel);
	}
}

// This one provides the error, but needs fitted lines to return it correctly

void addYErr(TGraphErrors* graph, TF1* yErrF) {
	for (Int_t i {0}; i < graph->GetN(); ++i) {
		graph->SetPointError(
			i,
			0., // half resolution?
			yErrF->Eval(graph->GetPointX(i))
		);
	}
}

void addXYErr(TGraphErrors* graph, Double_t xErr, Double_t yErr) {
	for (Int_t i {0}; i < graph->GetN(); ++i) {
		graph->SetPointError(
			i,
			xErr, // half resolution?
			yErr
		);
	}
}
// Transform two sinusoids into a lissajous shape

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

Double_t thetaErr(Double_t x, Double_t y, Double_t xErr, Double_t yErr, Double_t cov) {
	// Double_t r { rFromPt(x, y) };

	/*
	Double_t cosTheta { x / r };
	Double_t sinTheta { y / r };
	Double_t tanTheta { y / x };
	*/

	return sqrt(y*y*xErr*xErr + x*x*yErr*yErr + x*y*cov) / pow(rFromPt(x, y), 2.);

	/*

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

	*/
}

TGraphErrors* polarize(TGraphErrors* lissajous, Double_t cov) {
	TGraphErrors* polarized = new TGraphErrors(lissajous->GetN());
	for (int i {0}; i < lissajous->GetN(); ++i) {
		polarized->SetPoint(i,
			rFromPt(lissajous->GetPointX(i), lissajous->GetPointY(i)),
			thetaFromPt(lissajous->GetPointX(i), lissajous->GetPointY(i))
		);
		polarized->SetPointError(i,
			rErr(lissajous->GetPointX(i), lissajous->GetPointY(i), lissajous->GetErrorX(i), lissajous->GetErrorY(i)),
			thetaErr(lissajous->GetPointX(i), lissajous->GetPointY(i),
							 lissajous->GetErrorX(i), lissajous->GetErrorY(i), cov));
	}
	return polarized;
}

// Fitting functions

Double_t amplitudeFreqRespR(Double_t* omega, Double_t* pars) {
	// [0] is V0
	// [1] is R
	// [2] is L
	// [3] is C
	return pars[0]*pars[1] /
		sqrt(
			pars[1]*pars[1] +
			pow((*omega*pars[2] - 1/(*omega * pars[3])), 2.)
		);
}

Double_t amplitudeFreqRespL(Double_t* omega, Double_t* pars) {
	// [0] is V0
	// [1] is R
	// [2] is L
	// [3] is C
	return *omega*pars[0]*pars[2] /
		sqrt(
			pars[1]*pars[1] +
			pow((*omega*pars[2] - 1/(*omega * pars[3])), 2.)
		);
}

Double_t amplitudeFreqRespC(Double_t *omega, Double_t* pars) {
	// [0] is V0
	// [1] is R
	// [2] is L
	// [3] is C
	return pars[0] /
		(*omega*pars[3]) /
		sqrt(
			pars[1]*pars[1] +
			pow((*omega*pars[2] - 1/(*omega * pars[3])), 2.)
		);
}

Double_t phaseFreqResp(Double_t* omega, Double_t* pars) {
	return pars[0] + atan( (1 - *omega**omega*pars[2]*pars[3]) / (*omega*pars[1]*pars[3]) );
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

	std::cout << "Started analyzing data.\n";
	
	TFile* results = new TFile("analyzedData.ROOT", "RECREATE");

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
		"sineC",
		"squareWave1Ch",
		"squareWave2Ch"
	};

	std::cout << "Started acquiring graphs.\n";

	TGraphErrors* ampFreqRespR = new TGraphErrors((dataDir + fileNames[0]).c_str());
	TGraphErrors* ampFreqRespL = new TGraphErrors((dataDir + fileNames[1]).c_str());
	TGraphErrors* ampFreqRespC = new TGraphErrors((dataDir + fileNames[2]).c_str());

	TGraphErrors* phaseFreqRespR = new TGraphErrors((dataDir + fileNames[3]).c_str());
	TGraphErrors* phaseFreqRespL = new TGraphErrors((dataDir + fileNames[4]).c_str());
	TGraphErrors* phaseFreqRespC = new TGraphErrors((dataDir + fileNames[5]).c_str());

	// acquire amplitude errors graphs
	
	TGraphErrors* ampErrsCh[4] {};

	for (int j {0}; j < 4; ++j) {
		if (
			std::ifstream(dataDir + fileNames[6] + std::to_string(j) + ".txt").good()
		) {
			ampErrsCh[j] = new TGraphErrors(
				(dataDir + fileNames[6] + std::to_string(j) + ".txt" ).c_str()
			);
			ampErrsCh[j]->Write();
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
		if (
			std::ifstream(dataDir + fileNames[7] + std::to_string(j) + ".txt").good()
		) {
			phiErrsCh[j] = new TGraphErrors(
				(dataDir + fileNames[7] + std::to_string(j) + ".txt" ).c_str()
			);
			phiErrsCh[j]->Write();
		}
	}

	// acquire phase offsets graphs
	TGraphErrors* phiOffsCh[4] {};

	for (int j {0}; j < 4; ++j) {
		if (
			std::ifstream(dataDir + fileNames[7] + std::to_string(j) + ".txt").good()
		) {
			phiOffsCh[j] = new TGraphErrors(
				(dataDir + fileNames[7] + std::to_string(j) + ".txt" ).c_str(),
				"%lg%*lg%lg"
			);
			phiOffsCh[j]->Write();
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

	// acquire square waves graphs
	
	TGraphErrors* squareWave1Ch[4] {};
	TGraphErrors* squareWave2Ch[4] {};

	i = 0;

	for (; i < 4; ++i) {
		squareWave1Ch[i] = new TGraphErrors((dataDir + fileNames[13] + std::to_string(i) + ".txt").c_str());
		squareWave2Ch[i] = new TGraphErrors((dataDir + fileNames[13] + std::to_string(i) + ".txt").c_str());
	}

	TGraphErrors* squareWaveGenR1 = new TGraphErrors(squareWave1Ch[0]->GetN(), squareWave1Ch[0]->GetY(), squareWave1Ch[1]->GetY());
	TGraphErrors* squareWaveGenR2 = new TGraphErrors(squareWave2Ch[0]->GetN(), squareWave2Ch[0]->GetY(), squareWave2Ch[1]->GetY());

	double stdDev1VCh[2] {
		TMath::StdDev(squareWaveGenR1->GetN(), squareWaveGenR1->GetX()),
		TMath::StdDev(squareWaveGenR1->GetN(), squareWaveGenR1->GetY())
	};
	Double_t cov1Ch01 {squareWaveGenR1->GetCovariance()};
	std::cout << "Found covariance for acquisition 1 of: " << cov1Ch01 << std::endl;
	Double_t stdDev2VCh[2] {
		TMath::StdDev(squareWaveGenR1->GetN(), squareWaveGenR1->GetX()),
		TMath::StdDev(squareWaveGenR1->GetN(), squareWaveGenR1->GetY())
	};
	Double_t cov2Ch01 {squareWaveGenR1->GetCovariance()};
	std::cout << "Found covariance for acquisition 2 of: " << cov2Ch01 << std::endl;
	
	// Start data prep

	std::cout << "Started data prep.\n";

	// Phase offset functions fitting
	
	TF1* phiFreqRespCorrectorCh[4] {};
	
	i = 0;
	for (; i < 4; ++i) {
		phiFreqRespCorrectorCh[i] = new TF1(("correctPhiFreqRespCh" + std::to_string(i)).c_str(), "pol1", 0., 1E6);
		phiOffsCh[i]->Fit(phiFreqRespCorrectorCh[i]);
		phiFreqRespCorrectorCh[i]->Write();
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
	std::cout << "Fitting stdDev functions.\n";

	TF1* amplitudeErrFCh[4] {};
	TF1* phaseErrFCh[4] {};

	for (i = 0; i < 4; ++i) {
		amplitudeErrFCh[i] = new TF1(("amplitudeErrFCh" + std::to_string(i)).c_str(), "pol1", 0., 1E6);
		ampErrsCh[i]->Fit(amplitudeErrFCh[i]);
		amplitudeErrFCh[i]->Write();
		phaseErrFCh[i] = new TF1(("phaseErrFCh" + std::to_string(i)).c_str(), "pol1", 0., 1E6);
		phiErrsCh[i]->Fit(phaseErrFCh[i]);
		phaseErrFCh[i]->Write();
	}

	std::cout << "Adding errors to graphs.\n";

	addYErr(ampFreqRespR, amplitudeErrFCh[1]);
	ampFreqRespR->Write();
	addYErr(ampFreqRespL, amplitudeErrFCh[2]);
	ampFreqRespL->Write();
	addYErr(ampFreqRespC, amplitudeErrFCh[3]);
	ampFreqRespC->Write();

	addYErr(phaseFreqRespR, phaseErrFCh[1]);
	phaseFreqRespR->Write();
	addYErr(phaseFreqRespL, phaseErrFCh[2]);
	phaseFreqRespL->Write();
	addYErr(phaseFreqRespC, phaseErrFCh[3]);
	phaseFreqRespC->Write();

	// coexpress sinusoids
	
	TGraphErrors* lissajousGenR1 = correlate(sines[0][0], sines[1][0]);
	addXYErr(lissajousGenR1, stdDev1VCh[0], stdDev1VCh[1]);
	lissajousGenR1->Write();
	TGraphErrors* lissajousGenR2 = correlate(sines[0][1], sines[1][1]);
	addXYErr(lissajousGenR2, stdDev1VCh[0], stdDev1VCh[1]);
	lissajousGenR2->Write();
	TGraphErrors* lissajousGenR3 = correlate(sines[0][2], sines[1][2]);
	addXYErr(lissajousGenR2, stdDev1VCh[0], stdDev1VCh[1]);
	lissajousGenR3->Write();

	TGraphErrors* polarLissajousGenR1 = polarize(lissajousGenR1, cov1Ch01);
	polarLissajousGenR1->Write();
	TGraphErrors* polarLissajousGenR2 = polarize(lissajousGenR2, cov1Ch01);
	polarLissajousGenR2->Write();
	TGraphErrors* polarLissajousGenR3 = polarize(lissajousGenR2, cov1Ch01);
	polarLissajousGenR3->Write();

	std::cout << "Started fitting functions.\n";

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
	ampFreqRespRF->Write();
	ampFreqRespL->Fit(ampFreqRespLF);
	ampFreqRespLF->Write();
	ampFreqRespC->Fit(ampFreqRespCF);
	ampFreqRespCF->Write();

	phaseFreqRespR->Fit(phaseFreqRespRF);
	phaseFreqRespRF->Write();
	phaseFreqRespL->Fit(phaseFreqRespLF);
	phaseFreqRespLF->Write();
	phaseFreqRespC->Fit(phaseFreqRespCF);
	phaseFreqRespCF->Write();

	polarLissajousGenR1->Fit(polarLissajousGenR1F);
	polarLissajousGenR1F->Write();
	polarLissajousGenR2->Fit(polarLissajousGenR2F);
	polarLissajousGenR2F->Write();
	polarLissajousGenR3->Fit(polarLissajousGenR3F);
	polarLissajousGenR2F->Write();

	std::cout << "Done!" << std::endl;

	// space for cosmetic editing

	// saving EVERYTHING to TFile
	
	results->Close();

}

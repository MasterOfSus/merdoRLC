//#include <fstream>
//#include "/home/michele/root/include/TGraph.h"
//#include "/home/michele/root/include/TFormula.h"
//#include <cmath>
//#include "/home/michele/root/include/TMath.h"
#include <Rtypes.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TFormula.h>
#include <TAxis.h>
#include <TColor.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TLegend.h>
#include <TGlobal.h>
#include <TStyle.h>
#include <TPad.h>
#include <TH2F.h>
#include <TH1F.h>
#include <fstream>
#include <functional>

// Measured values

Double_t VPP{5.};
Double_t R{1.499E3};
Double_t RErr{2.};
Double_t L{10.11E-3};
Double_t LErr{10E-5};
Double_t LR{38.38};
Double_t LRErr{3E-2};
Double_t C{21.91E-9};
Double_t CErr{22E-11};

// Data processing functions

void correctPhi2PISystErr(TGraphErrors *phiFreqResp)
{
	for (int i{0}; i < phiFreqResp->GetN(); ++i)
	{
		if (phiFreqResp->GetPointY(i) < -M_PI / 2.)
			phiFreqResp->SetPointY(i, phiFreqResp->GetPointY(i) + 2 * M_PI);
		//	phiFreqResp->SetPointY(i, phiFreqResp->GetPointY(i) - 0.5 * M_PI);
	}
}

void correctPhiFreqResp(TGraphErrors *phiFreqResp, const TF1 *phaseOffsetF)
{
	for (Int_t i{0}; i < phiFreqResp->GetN(); ++i)
	{
		phiFreqResp->SetPoint(i, phiFreqResp->GetPointX(i), phiFreqResp->GetPointY(i) - phaseOffsetF->Eval(phiFreqResp->GetPointX(i)));
		phiFreqResp->SetPointError(i, phiFreqResp->GetErrorX(i),
								   sqrt(
									   pow((1 - phaseOffsetF->GetParameter(1)) * phiFreqResp->GetErrorY(i), 2.) +
									   pow(phiFreqResp->GetPointY(i) * phaseOffsetF->GetParError(1), 2.) +
									   pow(phaseOffsetF->GetParameter(0), 2.)));
	}
}

void correctSins(TGraphErrors *sinusoid, Double_t deltaT, Int_t nChannel)
{
	for (Int_t i{0}; i < sinusoid->GetN(); ++i)
	{
		sinusoid->SetPointX(i, sinusoid->GetPointX(i) + deltaT * nChannel);
	}
}

// This one provides the error, but needs fitted lines to return it correctly

void addYErr(TGraphErrors *graph, TF1 *yErrF)
{
	for (Int_t i{0}; i < graph->GetN(); ++i)
	{
		graph->SetPointError(
			i,
			0., // half resolution?
			yErrF->Eval(graph->GetPointX(i)));
	}
}

void addXYErr(TGraphErrors *graph, Double_t xErr, Double_t yErr)
{
	for (Int_t i{0}; i < graph->GetN(); ++i)
	{
		graph->SetPointError(
			i,
			xErr, // half resolution?
			yErr);
	}
}
// Transform two sinusoids into a lissajous shape

TGraphErrors *correlate(TGraphErrors *sin1, TGraphErrors *sin2)
{
	return new TGraphErrors(sin1->GetN(), sin1->GetY(), sin2->GetY(), sin1->GetEY(), sin2->GetEY());
}

// auxiliary polarization functions

Double_t rFromPt(Double_t x, Double_t y)
{
	return sqrt(x * x + y * y);
}

Double_t rErr(Double_t x, Double_t y, Double_t xErr, Double_t yErr)
{
	return (std::abs(x) * xErr + std::abs(y) * yErr) / rFromPt(x, y);
}

Double_t thetaFromPt(Double_t x, Double_t y)
{
	Double_t r{rFromPt(x, y)};
	Double_t cosTheta{x / r};
	Double_t sinTheta{y / r};
	Double_t tanTheta{y / x};

	bool useCos{false};
	bool useSin{false};
	bool useTan{false};

	if (std::abs(tanTheta) < std::abs(cosTheta))
	{
		if (std::abs(tanTheta) < std::abs(sinTheta))
			useTan = true;
		else
			useSin = true;
	}
	else
	{
		if (std::abs(cosTheta) < std::abs(sinTheta))
			useCos = true;
		else
			useSin = true;
	}

	if (useCos)
	{
		if (y >= 0.)
			return acos(cosTheta);
		else
			return 2 * M_PI - acos(cosTheta);
	}

	if (useSin)
	{
		if (x >= 0.)
		{
			if (sinTheta >= 0.)
				return asin(sinTheta);
			else
				return asin(sinTheta) + 2 * M_PI;
		}
		else
		{
			return M_PI - asin(sinTheta);
		}
	}

	if (useTan)
	{
		if (x >= 0.)
		{
			if (tanTheta >= 0.)
				return atan(tanTheta);
			else
				return atan(tanTheta) + 2 * M_PI;
		}
		else
			return M_PI + atan(tanTheta);
	}

	throw std::runtime_error("Dude wtf");
}

Double_t thetaErr(Double_t x, Double_t y, Double_t xErr, Double_t yErr, Double_t cov)
{
	// Double_t r { rFromPt(x, y) };

	/*
	Double_t cosTheta { x / r };
	Double_t sinTheta { y / r };
	Double_t tanTheta { y / x };
	*/

	return sqrt(y * y * xErr * xErr + x * x * yErr * yErr + x * y * cov) / pow(rFromPt(x, y), 2.);

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

TGraphErrors *polarize(TGraphErrors *lissajous, Double_t cov)
{
	TGraphErrors *polarized = new TGraphErrors(lissajous->GetN());
	for (int i{0}; i < lissajous->GetN(); ++i)
	{
		polarized->SetPoint(i,
							thetaFromPt(lissajous->GetPointX(i), lissajous->GetPointY(i)),
							rFromPt(lissajous->GetPointX(i), lissajous->GetPointY(i)));
		polarized->SetPointError(i,
								 thetaErr(lissajous->GetPointX(i), lissajous->GetPointY(i),
										  lissajous->GetErrorX(i), lissajous->GetErrorY(i), cov),
								 rErr(lissajous->GetPointX(i), lissajous->GetPointY(i), lissajous->GetErrorX(i), lissajous->GetErrorY(i)));
	}
	return polarized;
}

// Fitting functions

Double_t phiErrF(Double_t *freak, Double_t *pars)
{
	if (*freak < 10000.)
		return pars[0] + *freak * pars[1] - pars[1] * 10000.;
	else
		return pars[0] + *freak * pars[2] - pars[2] * 10000.;
};

Double_t amplitudeFreqRespR(Double_t *freq, Double_t *pars)
{
	// [0] is V0
	// [1] is R
	// [2] is L
	// [3] is C
	Double_t omega{*freq * 2 * M_PI};
	return pars[0] * pars[1] /
		   sqrt(
			   pars[1] * pars[1] +
			   pow((omega * pars[2] - 1 / (omega * pars[3])), 2.));
}

Double_t amplitudeFreqRespL(Double_t *freq, Double_t *pars)
{
	// [0] is V0
	// [1] is R
	// [2] is L
	// [3] is C
	Double_t omega{*freq * 2 * M_PI};
	return omega * pars[0] * pars[2] /
		   sqrt(
			   pars[1] * pars[1] +
			   pow((omega * pars[2] - 1 / (omega * pars[3])), 2.));
}

Double_t amplitudeFreqRespC(Double_t *freq, Double_t *pars)
{
	// [0] is V0
	// [1] is R
	// [2] is L
	// [3] is C
	Double_t omega{*freq * 2 * M_PI};
	return pars[0] /
		   (omega * pars[3]) /
		   sqrt(
			   pars[1] * pars[1] +
			   pow((omega * pars[2] - 1 / (omega * pars[3])), 2.));
}

Double_t phaseFreqResp(Double_t *freq, Double_t *pars)
{
	// [0] is either pi/2 or 0 or -pi/2
	// [1] is R
	// [2] is L
	// [3] is C
	Double_t omega{*freq * 2 * M_PI};
	return (pars[0] + atan((1 / (omega * pars[3]) - omega * pars[2]) / pars[1]));
}

Double_t polarLissajous(Double_t *theta, Double_t *pars)
{
	// [0] is V0
	// [1] is a = sqrt(R^2 + (omega*L - 1/omega*C))
	// [2] is phi = arctan((1 - omega^2*L*C/omega*R*C)

	Double_t V0{pars[0]};
	Double_t a{pars[1]};
	Double_t phi{pars[2]};

	Double_t tanTau{1 / tan(phi) - tan(*theta) / (a * sin(phi))};
	Double_t cosPhi{cos(phi)};
	Double_t sinPhi{sin(phi)};

	return V0 * sqrt((1 + a * a * cosPhi * cosPhi + a * a * tanTau * tanTau * sinPhi * sinPhi - 2 * a * a * cosPhi * sinPhi * tanTau) / (1 + tanTau * tanTau));
}

// correcting offset of generator phase
void correctOffset(TGraphErrors *gen, TGraphErrors *graph)
{
	for (int i = 0; i < gen->GetN(); ++i)
	{
		graph->SetPoint(i, graph->GetPointX(i), graph->GetPointY(i) - gen->GetPointY(i));
		graph->SetPointError(i, graph->GetErrorX(i), graph->GetErrorY(i) + gen->GetErrorY(i));
	}
}
// Data processor
Double_t genAmpFreqResp(Double_t* freq, Double_t* pars) {
  Double_t parss[4] {
    pars[1], pars[2], pars[3], pars[4]
  };
  return pars[0] - amplitudeFreqRespR(freq, parss);
}

void analyze(std::string dataDir)
{

	std::cout << "Started analyzing data.\n";

	TFile *results = new TFile("analyzedData.root", "RECREATE");

	std::vector<std::string> fileNames{
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
		"squareWave2Ch"};

	std::cout << "Started acquiring graphs.\n";

	TGraphErrors *ampFreqRespR = new TGraphErrors((dataDir + fileNames[0]).c_str());
	ampFreqRespR->SetName(fileNames[0].c_str());
	ampFreqRespR->SetDrawOption("APE");
	TGraphErrors *ampFreqRespL = new TGraphErrors((dataDir + fileNames[1]).c_str());
	ampFreqRespL->SetName(fileNames[1].c_str());
	ampFreqRespL->SetDrawOption("APE");
	TGraphErrors *ampFreqRespC = new TGraphErrors((dataDir + fileNames[2]).c_str());
	ampFreqRespC->SetName(fileNames[2].c_str());
	ampFreqRespC->SetDrawOption("APE");
	TGraphErrors *ampFreqRespV = new TGraphErrors((dataDir + "ampFreqRespV.txt").c_str());
	ampFreqRespV->SetName("ampFreqRespV.txt");
	ampFreqRespV->SetDrawOption("APE");

	TGraphErrors *phaseFreqRespV = new TGraphErrors((dataDir + "phaseFreqRespV.txt").c_str());
	phaseFreqRespV->SetName("rawPhaseFreqRespV.txt");
	phaseFreqRespV->Write();
	phaseFreqRespV->SetDrawOption("APE");
	phaseFreqRespV->SetName("phaseFreqRespV.txt");

	TGraphErrors *phaseFreqRespR = new TGraphErrors((dataDir + fileNames[3]).c_str());
	phaseFreqRespR->SetName("rawPhaseFreqRespR.txt");
	phaseFreqRespR->Write();
	phaseFreqRespR->SetName(fileNames[3].c_str());
	phaseFreqRespR->SetDrawOption("APE");
	correctPhi2PISystErr(phaseFreqRespR);

	TGraphErrors *phaseFreqRespL = new TGraphErrors((dataDir + fileNames[4]).c_str());
	phaseFreqRespL->SetName("PhaseFreqRespL.txt");
	phaseFreqRespL->Write();
	phaseFreqRespL->SetName(fileNames[4].c_str());
	phaseFreqRespL->SetDrawOption("APE");
	correctPhi2PISystErr(phaseFreqRespL);
	TGraphErrors *phaseFreqRespC = new TGraphErrors((dataDir + fileNames[5]).c_str());
	phaseFreqRespC->SetName("rawPhaseFreqRespC.txt");
	phaseFreqRespC->Write();
	phaseFreqRespC->SetName(fileNames[5].c_str());
	phaseFreqRespC->SetDrawOption("APE");
	correctPhi2PISystErr(phaseFreqRespC);

	// acquire amplitude errors graphs

	TGraphErrors *ampErrsCh[4]{};

	for (int j{0}; j < 4; ++j)
	{
		if (
			std::ifstream(dataDir + fileNames[6] + std::to_string(j) + ".txt").good())
		{
			ampErrsCh[j] = new TGraphErrors(
				(dataDir + fileNames[6] + std::to_string(j) + ".txt").c_str());
			ampErrsCh[j]->SetName((fileNames[6] + std::to_string(j) + ".txt").c_str());
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

	TGraphErrors *phiOffsCh[4]{};

	for (int j{0}; j < 4; ++j)
	{
		if (
			std::ifstream(dataDir + fileNames[7] + std::to_string(j) + ".txt").good())
		{
			phiOffsCh[j] = new TGraphErrors(
				(dataDir + fileNames[7] + std::to_string(j) + ".txt").c_str());
			phiOffsCh[j]->SetName(("phiOffsCh" + std::to_string(j) + ".txt").c_str());
		}
	}

	// acquire phase error graphs
	TGraphErrors *phiErrsCh[4]{};

	for (int j{0}; j < 4; ++j)
	{
		if (
			std::ifstream(dataDir + fileNames[7] + std::to_string(j) + ".txt").good())
		{
			phiErrsCh[j] = new TGraphErrors(
				(dataDir + fileNames[7] + std::to_string(j) + ".txt").c_str(),
				"%lg%*lg%lg");
			phiErrsCh[j]->SetName((fileNames[7] + std::to_string(j) + ".txt").c_str());
			phiErrsCh[j]->Write();
		}
	}

	// acquire sinusoids graphs

	std::vector<TGraphErrors *> sineGen{};
	std::vector<TGraphErrors *> sineR{};
	std::vector<TGraphErrors *> sineL{};
	std::vector<TGraphErrors *> sineC{};

	std::vector<TGraphErrors *> sines[4]{
		sineGen, sineR, sineL, sineC};

	Int_t i{0};

	for (int j{0}; j < 4; ++j)
	{
		i = 1;
		while (
			std::ifstream(dataDir + fileNames[8 + j] + std::to_string(i) + ".txt").good())
		{
			sines[j].emplace_back(new TGraphErrors(
				(dataDir + fileNames[8 + j] + std::to_string(i) + ".txt").c_str()));
			sines[j][i - 1]->SetName((fileNames[8 + j] + std::to_string(i) + ".txt").c_str());
			sines[j][i - 1]->Write();
			sines[j][i - 1]->SetDrawOption("APE");
			++i;
		}
	}

	// acquire square waves graphs

	TGraphErrors *squareWave1Ch[4]{};
	TGraphErrors *squareWave2Ch[4]{};

	i = 0;

	for (; i < 4; ++i)
	{
		squareWave1Ch[i] = new TGraphErrors((dataDir + fileNames[12] + std::to_string(i) + ".txt").c_str());
		squareWave2Ch[i] = new TGraphErrors((dataDir + fileNames[13] + std::to_string(i) + ".txt").c_str());
		for (int j{squareWave1Ch[i]->GetN() / 4}; j >= 0; --j)
		{
			squareWave1Ch[i]->RemovePoint(j);
			squareWave2Ch[i]->RemovePoint(j);
		};
		std::cout << "\n";
	}
	squareWave1Ch[0]->SetName("squareWave1Ch0");
	squareWave1Ch[0]->Write();
	squareWave2Ch[0]->SetName("squareWave2Ch0");
	squareWave2Ch[0]->Write();
	squareWave1Ch[1]->SetName("squareWave1Ch1");
	squareWave1Ch[1]->Write();
	squareWave2Ch[1]->SetName("squareWave2Ch1");
	squareWave2Ch[1]->Write();

	TGraphErrors *squareWaveGenR1 = new TGraphErrors(squareWave1Ch[0]->GetN(), squareWave1Ch[0]->GetY(), squareWave1Ch[1]->GetY());
	squareWaveGenR1->SetName("squareWaveGenR1");
	squareWaveGenR1->Write();
	TGraphErrors *squareWaveGenR2 = new TGraphErrors(squareWave2Ch[0]->GetN(), squareWave2Ch[0]->GetY(), squareWave2Ch[1]->GetY());
	squareWaveGenR2->SetName("squareWaveGenR2");
	squareWaveGenR2->Write();

	Double_t stdDev1VCh[2]{
		TMath::StdDev(squareWaveGenR1->GetN(), squareWaveGenR1->GetX()),
		TMath::StdDev(squareWaveGenR1->GetN(), squareWaveGenR1->GetY())};
	Double_t cov1Ch01{squareWaveGenR1->GetCovariance()};
	std::cout << "Found covariance for acquisition 1 of: " << cov1Ch01 << " and voltage stdDev equal to " << stdDev1VCh[0] << std::endl;
	Double_t stdDev2VCh[2]{
		TMath::StdDev(squareWaveGenR1->GetN(), squareWaveGenR1->GetX()),
		TMath::StdDev(squareWaveGenR1->GetN(), squareWaveGenR1->GetY())};
	Double_t cov2Ch01{squareWaveGenR1->GetCovariance()};
	std::cout << "Found covariance for acquisition 2 of: " << cov2Ch01 << std::endl;

	// Start data prep

	std::cout << "Started data prep.\n";

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

	TF1 *amplitudeErrFCh[4]{};
	TF1 *phaseErrFCh[4]{};

	for (i = 0; i < 4; ++i)
	{
		amplitudeErrFCh[i] = new TF1(("amplitudeErrFCh" + std::to_string(i)).c_str(), "pol0", 0., 1E6);
		ampErrsCh[i]->Fit(amplitudeErrFCh[i]);
		amplitudeErrFCh[i]->Write();
		phaseErrFCh[i] = new TF1(("phaseErrFCh" + std::to_string(i)).c_str(), phiErrF, 0., 4E4, 3);
		phaseErrFCh[i]->SetParameters(0.00225, -0.0001, 0.0001);
		phiErrsCh[i]->Fit(phaseErrFCh[i]);
		phaseErrFCh[i]->SetNpx(1000);
		phaseErrFCh[i]->Write();
	}

	std::cout << "Adding errors to graphs.\n";

	addYErr(ampFreqRespV, amplitudeErrFCh[0]);
	ampFreqRespV->Write();
	addYErr(ampFreqRespR, amplitudeErrFCh[1]);
	ampFreqRespR->Write();
	addYErr(ampFreqRespL, amplitudeErrFCh[2]);
	ampFreqRespL->Write();
	addYErr(ampFreqRespC, amplitudeErrFCh[3]);
	ampFreqRespC->Write();

	addYErr(phaseFreqRespV, phaseErrFCh[0]);
	addYErr(phaseFreqRespR, phaseErrFCh[1]);
	addYErr(phaseFreqRespL, phaseErrFCh[2]);
	addYErr(phaseFreqRespC, phaseErrFCh[3]);

	// phase offset of generator
	correctOffset(phaseFreqRespV, phaseFreqRespR);
	correctOffset(phaseFreqRespV, phaseFreqRespL);
	correctOffset(phaseFreqRespV, phaseFreqRespC);
	// Phase offset functions fitting
	for (int i = 0; i < phaseFreqRespV->GetN(); ++i)
	{
		phaseFreqRespV->SetPoint(i, phaseFreqRespV->GetX()[i], phaseFreqRespV->GetY()[i] - M_PI / 2);
	}

	TF1 *phiFreqRespCorrectorCh[4]{};

	i = 0;

	for (TGraphErrors *g : phiOffsCh)
	{
		addYErr(g, phaseErrFCh[i]);
		g->Write();
		++i;
	}

	Double_t expSlopes[4]{0., 0.000001, 0.000002, 0.000003};

	i = 0;
	for (; i < 4; ++i)
	{
		std::cout << "Fitting phi offset function for channel " << i << std::endl;
		phiFreqRespCorrectorCh[i] = new TF1(("correctPhiFreqRespCh" + std::to_string(i)).c_str(), "pol1", 0., 4E4);
		phiFreqRespCorrectorCh[i]->SetParameter(1, expSlopes[i]);
		phiFreqRespCorrectorCh[i]->SetParameter(0, 0.);
		phiOffsCh[i]->Fit(phiFreqRespCorrectorCh[i]);
		phiFreqRespCorrectorCh[i]->Write();
	}

	// correcting phase freq resp graph for syst phase offset
	correctPhiFreqResp(phaseFreqRespV, phiFreqRespCorrectorCh[0]);
	phaseFreqRespV->Write();
	correctPhiFreqResp(phaseFreqRespR, phiFreqRespCorrectorCh[1]);
	phaseFreqRespR->Write();
	correctPhiFreqResp(phaseFreqRespL, phiFreqRespCorrectorCh[2]);
	phaseFreqRespL->Write();
	correctPhiFreqResp(phaseFreqRespC, phiFreqRespCorrectorCh[3]);
	phaseFreqRespC->Write();

	// coexpress sinusoids

	TGraphErrors *lissajousGenR1 = correlate(sines[0][0], sines[1][0]);
	addXYErr(lissajousGenR1, stdDev1VCh[0], stdDev1VCh[1]);
	lissajousGenR1->SetDrawOption("APE");
	lissajousGenR1->Write();
	TGraphErrors *lissajousGenR2 = correlate(sines[0][1], sines[1][1]);
	addXYErr(lissajousGenR2, stdDev1VCh[0], stdDev1VCh[1]);
	lissajousGenR2->SetDrawOption("APE");
	lissajousGenR2->Write();
	TGraphErrors *lissajousGenR3 = correlate(sines[0][2], sines[1][2]);
	addXYErr(lissajousGenR3, stdDev1VCh[0], stdDev1VCh[1]);
	lissajousGenR3->SetDrawOption("APE");
	lissajousGenR3->Write();

	// Drawing Lissajouses
	TCanvas *ellipses = new TCanvas();
	ellipses->SetCanvasSize(1920, 1080);
	ellipses->SetName("Lissajous Ellipses");
	ellipses->Divide(2);
	ellipses->cd(1);
	lissajousGenR2->SetTitle("Resonance Frequency");
	lissajousGenR2->GetYaxis()->SetTitle("Gen Voltage(V)");
	lissajousGenR2->GetXaxis()->SetTitle("R Voltage(V)");
	lissajousGenR2->Draw("APE");
	ellipses->cd(2);
	lissajousGenR3->SetTitle("Over Resonance Frequency (15kHz)");
	lissajousGenR3->GetYaxis()->SetTitle("Gen Voltage(V)");
	lissajousGenR3->GetXaxis()->SetTitle("R Voltage(V)");
	lissajousGenR3->Draw("APE");
	ellipses->SaveAs("lissellips.gif");
	ellipses->Write();

	TMultiGraph *polarLissajouses = new TMultiGraph();

	TGraphErrors *polarLissajousGenR1 = polarize(lissajousGenR1, cov1Ch01);
	polarLissajousGenR1->Write();
	polarLissajouses->Add(polarLissajousGenR1, "APE");
	TGraphErrors *polarLissajousGenR2 = polarize(lissajousGenR2, cov1Ch01);
	polarLissajousGenR2->Write();
	polarLissajouses->Add(polarLissajousGenR2, "APE");
	TGraphErrors *polarLissajousGenR3 = polarize(lissajousGenR3, cov1Ch01);
	polarLissajousGenR3->Write();
	polarLissajouses->Add(polarLissajousGenR3, "APE");
	polarLissajouses->Write();

	std::cout << "Started fitting functions.\n";

	// fit functions, need to initialize values

	// setting inverse lorentzian for Gen Fit
	//double_t dip = ampFreqRespV->GetMaximum() - ampFreqRespV->GetMinimum();
	//double_t gamma = R / (4 * M_PI * L);
	//double_t resfreq = 10730.;
//
	// TF1 *ampFreqRespVF = new TF1("ampFreqRespVF",
	//							 "[0] - [1]/(1 + TMath::Power((x - [2])/[3], 2))",
	//							 0., 1E6);
	//ampFreqRespVF->SetParameters(VPP/2, VPP / 2, dip, resfreq, gamma);
	TF1* ampFreqRespVF = new TF1("genAmpFreqRespF", genAmpFreqResp, 0., 1E6, 5);
	ampFreqRespVF->SetNpx(1000);
	ampFreqRespVF->SetParameters(VPP/2,VPP/2,R,L,C);
	TF1 *ampFreqRespRF = new TF1("ampFreqRespRF", amplitudeFreqRespR, 0., 1E6, 4);
	ampFreqRespRF->SetNpx(1000);
	ampFreqRespRF->SetParameters(VPP / 2., R, L, C);
	TF1 *ampFreqRespLF = new TF1("ampFreqRespLF", amplitudeFreqRespL, 1E3, 4E4, 4);
	ampFreqRespLF->SetNpx(1000);
	ampFreqRespLF->SetParameters(VPP / 2., R, L, C);
	TF1 *ampFreqRespCF = new TF1("ampFreqRespCF", amplitudeFreqRespC, 1E3, 4E4, 4);
	ampFreqRespCF->SetNpx(1000);
	ampFreqRespCF->SetParameters(VPP / 2., R, L, C);

	TF1 *phaseFreqRespVF = new TF1("phaseFreqRespVF", phaseFreqResp, 0., 1E6, 1);
	phaseFreqRespVF->SetParameters(VPP);
	phaseFreqRespVF->SetNpx(1000);
	TF1 *phaseFreqRespRF = new TF1("phaseFreqRespRF", phaseFreqResp, 0., 1E6, 4);
	phaseFreqRespRF->SetParameters(0., R, L, C);
	phaseFreqRespRF->SetNpx(1000);
	TF1 *phaseFreqRespLF = new TF1("phaseFreqRespLF", phaseFreqResp, 0., 1E6, 4);
	phaseFreqRespLF->SetParameters(0.5 * M_PI, R, L, C);
	phaseFreqRespLF->SetNpx(1000);
	TF1 *phaseFreqRespCF = new TF1("phaseFreqRespCF", phaseFreqResp, 0., 1E6, 4);
	phaseFreqRespCF->SetParameters(-0.5 * M_PI, R, L, C);
	phaseFreqRespCF->SetNpx(1000);

	TF1 *polarLissajousGenR1F = new TF1("polarLissajousGenR1F", polarLissajous, 0., 2 * M_PI, 3);
	/*	polarLissajousGenR1F->SetParameters(
			VPP/2.,
			R/sqrt(R*R + pow((2*M_PI*5E3*L - 1/(2*M_PI*5E3*C)), 2.)),
			atan((1 - 4*M_PI*M_PI*25E6*L*C)/(2*M_PI*5E3*R*C))
		);*/
	polarLissajousGenR1F->SetParameters(
		VPP / 2.,
		R / sqrt(R * R + pow((2 * M_PI * 7E3 * L - 1 / (2 * M_PI * 7E3 * C)), 2.)),
		atan((1 - 4 * M_PI * M_PI * 49E6 * L * C) / (2 * M_PI * 7E3 * R * C)));
	//polarLissajousGenR1F->FixParameter(0, VPP / 2.);
	polarLissajousGenR1F->SetNpx(1000);
	TF1 *polarLissajousGenR2F = new TF1("polarLissajousGenR2F", polarLissajous, 0., 2 * M_PI, 3);
	polarLissajousGenR2F->SetParameters(
		VPP / 2.,
		R / sqrt(R * R + pow((2 * M_PI * 10E3 * L - 1 / (2 * M_PI * 10E3 * C)), 2.)),
		atan((1 - 4 * M_PI * M_PI * 1E8 * L * C) / (2 * M_PI * 1E4 * R * C)));
	//polarLissajousGenR2F->FixParameter(0, VPP / 2.);
	polarLissajousGenR2F->SetNpx(1000);
	TF1 *polarLissajousGenR3F = new TF1("polarLissajousGenR3F", polarLissajous, 0., 2 * M_PI, 3);
	polarLissajousGenR3F->SetParameters(
		VPP / 2.,
		R / sqrt(R * R + pow((2 * M_PI * 15E3 * L - 1 / (2 * M_PI * 15E3 * C)), 2.)),
		atan((1 - 4 * M_PI * M_PI * 2.25E8 * L * C) / (2 * M_PI * 1.5E4 * R * C)));
	//polarLissajousGenR3F->FixParameter(0, VPP / 2.);
	polarLissajousGenR3F->SetNpx(1000);

	// fitting
	std::cout << "Generator amplitude frequency response fit:\n";
	ampFreqRespV->Fit(ampFreqRespVF, "NR");
	ampFreqRespVF->Write();
	std::cout << "R amplitude frequency response fit:\n";
	ampFreqRespR->Fit(ampFreqRespRF, "NR");
	ampFreqRespRF->Write();
	std::cout << "L amplitude frequency response fit:\n";
	ampFreqRespL->Fit(ampFreqRespLF, "NR");
	ampFreqRespLF->Write();
	std::cout << "C amplitude frequency response fit:\n";
	ampFreqRespC->Fit(ampFreqRespCF, "NR");
	ampFreqRespCF->Write();

	std::cout << "Generator phase frequency response fit:\n";
	phaseFreqRespV->Fit(phaseFreqRespVF, "NR");
	phaseFreqRespVF->Write();
	std::cout << "R phase frequency response fit:\n";
	phaseFreqRespR->Fit(phaseFreqRespRF, "NR");
	phaseFreqRespRF->Write();
	std::cout << "L phase frequency response fit:\n";
	phaseFreqRespL->Fit(phaseFreqRespLF, "NR", "", 5000., 40000.);
	phaseFreqRespLF->Write();
	std::cout << "C phase frequency response fit:\n";
	phaseFreqRespC->Fit(phaseFreqRespCF, "NR", "", 1000., 40000.);
	phaseFreqRespCF->Write();

	std::cout << "Polar lissajous genR1 fit:\n";
	polarLissajousGenR1->Fit(polarLissajousGenR1F);
	polarLissajousGenR1F->Write();
	std::cout << "Polar lissajous genR2 fit:\n";
	polarLissajousGenR2->Fit(polarLissajousGenR2F);
	polarLissajousGenR2F->Write();
	std::cout << "Polar lissajous genR3 fit:\n";
	polarLissajousGenR3->Fit(polarLissajousGenR3F);
	polarLissajousGenR3F->Write();

	std::cout << "Started drawing stuff." << std::endl;

	TCanvas *ampFreqRespCnvs = new TCanvas();
	ampFreqRespCnvs->SetName("ampFreqRespCnvs");
	// ampFreqRespCnvs->Divide(3);
	ampFreqRespCnvs->cd(1);
	ampFreqRespR->SetTitle("Amplitude Frequency Response");
	ampFreqRespR->GetXaxis()->SetTitle("Frequency(hz)");
	ampFreqRespR->GetYaxis()->SetTitle("Amplitude(V)");
	ampFreqRespR->SetMinimum(0.);
	ampFreqRespR->Draw("APE");
	ampFreqRespRF->SetNpx(10000);
	ampFreqRespRF->SetTitle("R Fit");
	ampFreqRespRF->Draw("SAME");
	// ampFreqRespCnvs->cd(2);
	ampFreqRespL->Draw("SAME");
	ampFreqRespLF->SetLineColor(3);
	ampFreqRespLF->SetTitle("L Fit");
	ampFreqRespLF->SetNpx(10000);
	ampFreqRespLF->Draw("SAME");
	// ampFreqRespCnvs->cd(3);
	ampFreqRespC->Draw("SAME");
	ampFreqRespCF->SetLineColor(4);
	ampFreqRespCF->SetTitle("C Fit");
	ampFreqRespCF->SetNpx(10000);
	ampFreqRespCF->Draw("SAME");
	ampFreqRespV->Draw("SAME");
	ampFreqRespVF->SetLineColor(6);
	ampFreqRespVF->SetTitle("VGen Fit");
	ampFreqRespVF->SetNpx(10000);
	ampFreqRespVF->Draw("SAME");
	ampFreqRespCnvs->Write();
	//

	TCanvas *phaseFreqRespCnvs = new TCanvas();
	phaseFreqRespCnvs->SetName("phaseFreqRespCnvs");
	phaseFreqRespCnvs->Divide(3);
	phaseFreqRespCnvs->cd(1);
	phaseFreqRespR->Draw("APE");
	phaseFreqRespRF->Draw("SAME");
	phaseFreqRespCnvs->cd(2);
	phaseFreqRespL->Draw("APE");
	phaseFreqRespLF->Draw("SAME");
	phaseFreqRespCnvs->cd(3);
	phaseFreqRespC->Draw("APE");
	phaseFreqRespCF->Draw("SAME");
	phaseFreqRespCnvs->Write();

	TCanvas *multiPhaseCnvs = new TCanvas();
	multiPhaseCnvs->SetName("multiPhaseCnvs");
	multiPhaseCnvs->Divide();
	TMultiGraph *phaseFreqRespGrphs = new TMultiGraph();
	phaseFreqRespGrphs->Add(phaseFreqRespR, "APE");
	phaseFreqRespGrphs->Add(phaseFreqRespL, "APE");
	phaseFreqRespGrphs->Add(phaseFreqRespC, "APE");
	multiPhaseCnvs->cd(1);
	phaseFreqRespGrphs->SetTitle("Phase Frequency Response");
	phaseFreqRespGrphs->GetXaxis()->SetTitle("Frequency (Hz)");
	phaseFreqRespGrphs->GetYaxis()->SetTitle("Phase (rad)");
	phaseFreqRespGrphs->SetMaximum(3.);
	phaseFreqRespGrphs->Draw("APE");
	phaseFreqRespRF->SetTitle("R Fit");
	phaseFreqRespRF->Draw("SAME");
	phaseFreqRespLF->SetTitle("L Fit");
	phaseFreqRespLF->SetLineColor(3);
	phaseFreqRespLF->SetNpx(1000);
	phaseFreqRespLF->Draw("SAME");
	phaseFreqRespCF->SetTitle("C Fit");
	phaseFreqRespCF->SetLineColor(4);
	phaseFreqRespCF->SetNpx(1000);
	phaseFreqRespCF->Draw("SAME");
	phaseFreqRespGrphs->Write();
	multiPhaseCnvs->Write();

	TCanvas *lissajousCnvs = new TCanvas();
	lissajousCnvs->SetCanvasSize(1920, 1080);
	lissajousCnvs->SetName("lissajousCnvs");
	lissajousCnvs->Divide(2);
	// lissajousCnvs->cd(1);
	// polarLissajousGenR1->Draw("APE");
	// polarLissajousGenR1F->Draw("SAME");
	lissajousCnvs->cd(1);
	polarLissajousGenR2->SetTitle("Resonance Frequency");
	polarLissajousGenR2->GetYaxis()->SetTitle("Radius (V)");
	polarLissajousGenR2->GetXaxis()->SetTitle("Angle (rad)");
	polarLissajousGenR2->Draw("APE");
	polarLissajousGenR2F->Draw("SAME");
	lissajousCnvs->cd(2);
	polarLissajousGenR3->SetTitle("Over Resonance Frequency (15kHz)");
	polarLissajousGenR3->GetYaxis()->SetTitle("Radius (V)");
	polarLissajousGenR3->GetXaxis()->SetTitle("Angle (Rad)");
	polarLissajousGenR3->Draw("APE");
	polarLissajousGenR3F->Draw("SAME");
	lissajousCnvs->Write();
	lissajousCnvs->SaveAs("polarLiss.gif");

	// space for cosmetic editing

	// saving EVERYTHING to TFile
	
	// post analysis, presentation material making
	TCanvas* multiPhaseErrGraph = new TCanvas("multiPhaseErrGraph", "Used phase error functions");
	int compColors[4] {6, 2, 3, 4};

	phaseErrFCh[0]->SetLineColor(6);
	phaseErrFCh[0]->SetTitle("Funzioni di errore utilizzate - Fase");
	phaseErrFCh[0]->GetXaxis()->SetTitle("Frequenza (Hz)");
	phaseErrFCh[0]->GetYaxis()->SetTitle("Errore associato (Rad)");
	phaseErrFCh[0]->SetNpx(1000);
	phaseErrFCh[0]->SetLineWidth(1);
	phaseErrFCh[0]->Draw();
	phiErrsCh[0]->SetMarkerStyle(3);
	phiErrsCh[0]->SetMarkerColor(6);
	phiErrsCh[0]->Draw("SAME PE");
	for (int i {1}; i < 4; ++i) {
		phaseErrFCh[i]->SetLineColor(compColors[i]);
		phaseErrFCh[i]->SetNpx(1000);
		phaseErrFCh[i]->SetLineWidth(1);
		phaseErrFCh[i]->Draw("SAME");
		phiErrsCh[i]->SetMarkerStyle(3);
		phiErrsCh[i]->SetMarkerColor(compColors[i]);
		phiErrsCh[i]->Draw("SAME PE");
	}

	TLegend* phErrGrphLegend = new TLegend(.60, .25, .90, .45, "");
	phErrGrphLegend->AddEntry(phaseErrFCh[0], "CH0 (Generatore)");
	phErrGrphLegend->AddEntry(phaseErrFCh[1], "CH1 (Resistenza)");
	phErrGrphLegend->AddEntry(phaseErrFCh[2], "CH2 (Induttanza)");
	phErrGrphLegend->AddEntry(phaseErrFCh[3], "CH3 (Condensatore)");

	phErrGrphLegend->Draw();

	multiPhaseErrGraph->Write();

	gStyle->SetLineWidth(1);
	gStyle->SetLineColor(1);

	TCanvas* multiPhaseOffsCnvs = new TCanvas("multiPhaseOffsCnvs", "Used phase offset functions");
	phiFreqRespCorrectorCh[0]->SetLineColor(6);
	phiFreqRespCorrectorCh[0]->SetTitle("Funzioni di offset utilizzate");
	phiFreqRespCorrectorCh[0]->GetXaxis()->SetTitle("Frequenza (Hz)");
	phiFreqRespCorrectorCh[0]->GetYaxis()->SetTitle("Offset (Rad)");
	phiFreqRespCorrectorCh[0]->SetNpx(1000);
	phiFreqRespCorrectorCh[0]->SetLineWidth(1);
	phiFreqRespCorrectorCh[0]->SetLineStyle(1);
	phiFreqRespCorrectorCh[0]->GetYaxis()->SetRangeUser(-0.05, 0.9);
	phiFreqRespCorrectorCh[0]->SetBit(TF1::kNotDraw);  // Optional, to prevent redraws
	phiFreqRespCorrectorCh[0]->Draw();
	phiOffsCh[0]->Draw("SAME PE");
	phiOffsCh[0]->Write();

	for (i = 1; i < 4; ++i) {
		phiFreqRespCorrectorCh[i]->SetLineColor(compColors[i]);
		phiFreqRespCorrectorCh[i]->SetLineWidth(1);
		phiFreqRespCorrectorCh[i]->SetLineStyle(1);
		phiFreqRespCorrectorCh[i]->SetNpx(1000);
		phiFreqRespCorrectorCh[i]->Draw("SAME");
		phiOffsCh[i]->Write();
		phiFreqRespCorrectorCh[i]->SetBit(TF1::kNotDraw);  // Optional, to prevent redraws
		phiOffsCh[i]->Draw("SAME PE");
	}

	TLegend* phOffsGrphLegend = new TLegend(.60, .25, .90, .45, "");
	phOffsGrphLegend->AddEntry(phiFreqRespCorrectorCh[0], "CH0 (Generatore)");
	phOffsGrphLegend->AddEntry(phiFreqRespCorrectorCh[1], "CH1 (Resistenza)");
	phOffsGrphLegend->AddEntry(phiFreqRespCorrectorCh[2], "CH2 (Induttanza)");
	phOffsGrphLegend->AddEntry(phiFreqRespCorrectorCh[3], "CH3 (Condensatore)");

	phOffsGrphLegend->Draw();

	multiPhaseOffsCnvs->Write();

	TH2F* voltagesScPl = new TH2F("voltagesScPl", "Misure di voltaggi contemporanee", 20, 2.363, 2.384, 20, 2.366, 2.378);
	TCanvas* voltagesCnvs = new TCanvas("voltagesCnvs", "Misure di voltaggi contemporanee");
	for (int i {0}; i < squareWaveGenR1->GetN(); ++i) {
		voltagesScPl->Fill(squareWaveGenR1->GetPointX(i), squareWaveGenR1->GetPointY(i));
	}
	std::cout << "Covariance: " << squareWaveGenR1->GetCovariance() << std::endl;
	voltagesScPl->Draw("LEGO");
	TH1D* voltagesSqrWvCh0 = voltagesScPl->ProjectionX();
	voltagesSqrWvCh0->SetTitle("Misure di differenza di potenziale costante - CH0");
	voltagesSqrWvCh0->SetName("Risultati del fit");
	voltagesSqrWvCh0->GetXaxis()->SetTitle("DDP (V)");
	voltagesSqrWvCh0->GetYaxis()->SetTitle("Occorrenze");
	voltagesSqrWvCh0->SetFillColor(kBlue);
	TF1* gausV0 = new TF1("gausV0", "gaus", 2.363, 2.384);
	voltagesSqrWvCh0->Fit(gausV0);
	voltagesSqrWvCh0->Write();
	TH1D* voltagesSqrWvCh1 = voltagesScPl->ProjectionY();
	voltagesSqrWvCh1->SetTitle("Misure di differenza di potenziale costante - CH1");
	voltagesSqrWvCh1->SetName("Risultati del fit");
	voltagesSqrWvCh1->GetXaxis()->SetTitle("DDP (V)");
	voltagesSqrWvCh1->GetYaxis()->SetTitle("Occorrenze");
	TF1* gausV1 = new TF1("gausV0", "gaus", 2.366, 2.378);
	voltagesSqrWvCh1->Fit(gausV1);
	voltagesSqrWvCh1->SetFillColor(kBlue);
	voltagesSqrWvCh1->Write();

	voltagesScPl->Draw();
	voltagesScPl->Write();

	results->Close();
}

void makeImages(std::string dataDir) {

	TFile* presMaterial = new TFile("pptMaterial.root", "RECREATE");

	// Section for presentation material
	
	Int_t compColors[4] {6, 2, 3, 4};
	
	TF1* expAmpFreqRespGenF = new TF1("expAmpFreqRespGenF", "[0]", 0., 4E4);
	TF1* expAmpFreqRespRF = new TF1("expAmpFreqRespRF", amplitudeFreqRespR, 0., 4E4, 4);
	TF1* expAmpFreqRespLF = new TF1("expAmpFreqRespLF", amplitudeFreqRespL, 0., 4E4, 4);
	TF1* expAmpFreqRespCF = new TF1("expAmpFreqRespCF", amplitudeFreqRespC, 0., 4E4, 4);

	TF1* expAmpFreqResps[3] = {expAmpFreqRespRF, expAmpFreqRespLF, expAmpFreqRespCF};

	expAmpFreqRespGenF->SetParameters(VPP/2.);
	expAmpFreqRespGenF->SetLineColor(compColors[0]);
	expAmpFreqRespGenF->SetTitle("Risposte in frequenza attese - Ampiezza");
	expAmpFreqRespGenF->GetXaxis()->SetTitle("Frequenza (Hz)");
	expAmpFreqRespGenF->GetYaxis()->SetTitle("Ampiezza (V)");
	expAmpFreqRespGenF->GetYaxis()->SetRangeUser(0., 2.75);
	for (int i {0}; i < 3; ++i) {
		expAmpFreqResps[i]->SetParameters(VPP/2., R, L, C);
		expAmpFreqResps[i]->SetNpx(1000);
		expAmpFreqResps[i]->SetLineColor(compColors[i + 1]);
	}
	TLegend* expAmpLegend = new TLegend(.60, .25, .90, .45, "");

	expAmpLegend->AddEntry(expAmpFreqRespGenF, "Generatore");
	expAmpLegend->AddEntry(expAmpFreqRespRF, "Resistenza");
	expAmpLegend->AddEntry(expAmpFreqRespLF, "Induttanza");
	expAmpLegend->AddEntry(expAmpFreqRespCF, "Condensatore");

	TCanvas* expAmpFreqRespsCvs = new TCanvas("expAmpFreqRespsCvs", "Piscio pene culo");

	expAmpFreqRespGenF->Draw();
	for (int i {0}; i < 3; ++i) {
		expAmpFreqResps[i]->Draw("SAME");
	}
	expAmpLegend->Draw();

	TF1* expPhaseRespGen = new TF1("expPhaseRespGen", "0", 0., 4E4);
	TF1* expPhaseRespR = new TF1("expPhaseRespR", phaseFreqResp, 0., 4E4, 4);
	TF1* expPhaseRespL = new TF1("expPhaseRespL", phaseFreqResp, 0., 4E4, 4);
	TF1* expPhaseRespC = new TF1("expPhaseRespC", phaseFreqResp, 0., 4E4, 4);

	TF1* expPhaseResps[3] {expPhaseRespR, expPhaseRespL, expPhaseRespC};

	expPhaseRespGen->SetLineColor(compColors[0]);
	expPhaseRespGen->SetTitle("Risposte in frequenza attese - Fasi");
	expPhaseRespGen->GetXaxis()->SetTitle("Freqenza (Hz)");
	expPhaseRespGen->GetYaxis()->SetTitle("Fase (Rad)");
	expPhaseRespGen->GetYaxis()->SetRangeUser(-M_PI, 4.);

	Double_t expPhaseOffs[3] {0., -M_PI/2., M_PI/2.};
	for (int i {0}; i < 3; ++i) {
		expPhaseResps[i]->SetParameters(expPhaseOffs[i], R, L, C);
		expPhaseResps[i]->SetNpx(1000);
		expPhaseResps[i]->SetLineColor(compColors[i + 1]);
	}

	TLegend* expPhaseRespsLgnd = new TLegend(.60, .65, .90, .85, "");
	expPhaseRespsLgnd->AddEntry(expPhaseRespGen, "Generatore");
	expPhaseRespsLgnd->AddEntry(expPhaseRespR, "Resistenza");
	expPhaseRespsLgnd->AddEntry(expPhaseRespL, "Induttanza");
	expPhaseRespsLgnd->AddEntry(expPhaseRespC, "Condensatore");

	TCanvas* expPhaseRespsCvs = new TCanvas("expPhaseRespsCvs", "Piscio pene cazzo");

	expPhaseRespGen->Draw();
	for (int i {0}; i < 3; ++i) {
		expPhaseResps[i]->Draw("SAME");
	}
	expPhaseRespsLgnd->Draw();

	// Amplitude error graphs
	TGraph* ampValsHz[6] = {
		new TGraph((dataDir + "rawErrors/ampVals4kHzCh3.txt").c_str()),
		new TGraph((dataDir + "rawErrors/ampVals10kHzCh3.txt").c_str()),
		new TGraph((dataDir + "rawErrors/ampVals16kHzCh3.txt").c_str()),
		new TGraph((dataDir + "rawErrors/ampVals25kHzCh3.txt").c_str()),
		new TGraph((dataDir + "rawErrors/ampVals31kHzCh3.txt").c_str()),
		new TGraph((dataDir + "rawErrors/ampVals36kHzCh3.txt").c_str()),
	};
	TGraph* ampValsGrph = new TGraph();
	ampValsGrph->SetName("ampValsGrph");
	for (int i {0}; i < 6; ++i) {
		for (int j {0}; j < ampValsHz[i]->GetN(); ++j) {
			ampValsGrph->AddPoint(ampValsHz[i]->GetPointX(j), ampValsHz[i]->GetPointY(j) - ampValsHz[i]->GetMean(2));
		}
	}
	std::cout << "AmpValsGrph pt count: " << ampValsGrph->GetN() << std::endl;
	ampValsGrph->SetTitle("Scarti misurati per l'ampiezza - Canale 3");
	ampValsGrph->GetXaxis()->SetTitle("Frequenza (Hz)");
	ampValsGrph->GetYaxis()->SetTitle("Ampiezza (scarto) (V)");
	ampValsGrph->SetMarkerStyle(3);
	ampValsGrph->SetMarkerSize(0.5);

	TH1F* ampValsHisto = new TH1F("ampValsHisto", "Valori misurati - Ampiezza - Canale 3, 25 kHz", 15, 2.4963f, 2.4970f);
	for (int i {0}; i < ampValsHz[3]->GetN(); ++i) {
		ampValsHisto->Fill(ampValsHz[3]->GetPointY(i));
	}
	std::cout << "Overflows: " << ampValsHisto->GetBinContent(ampValsHisto->GetNbinsX() + 1) << " Underflows: " << ampValsHisto->GetBinContent(0) << std::endl;
	ampValsHisto->GetXaxis()->SetTitle("Ampiezza (Hz)");
	ampValsHisto->GetYaxis()->SetTitle("Occorrenze");
	ampValsHisto->SetFillColor(kBlue);
	TF1* gaus = new TF1("gaus", "gaus", 2.4963f, 2.4970f);
	gaus->SetNpx(1000);
	TCanvas* ampValsSctrPlt = new TCanvas("ampValsSctrPlt", "Scatter plot - Ampiezza, Canale 3");
	ampValsHisto->Fit(gaus, "Q");

	ampValsGrph->Draw("AP");

	TCanvas* ampValsH = new TCanvas("ampValsH", "culo cacca");
	ampValsHisto->Draw();

	// Phase error graphs
	TGraph* phaseValsGrphHz[6] {
		new TGraph((dataDir + "rawErrors/phaseVals4kHzCh2.txt").c_str()),
		new TGraph((dataDir + "rawErrors/phaseVals10kHzCh2.txt").c_str()),
		new TGraph((dataDir + "rawErrors/phaseVals16kHzCh2.txt").c_str()),
		new TGraph((dataDir + "rawErrors/phaseVals25kHzCh2.txt").c_str()),
		new TGraph((dataDir + "rawErrors/phaseVals31kHzCh2.txt").c_str()),
		new TGraph((dataDir + "rawErrors/phaseVals36kHzCh2.txt").c_str()),
	};

	TGraph* phaseValsGrph = new TGraph();
	phaseValsGrph->SetName("phaseValsGrph");
	for (int i {}; i < 6; ++i) {
		for (int j {}; j < phaseValsGrphHz[i]->GetN(); ++j) {
			phaseValsGrph->AddPoint(phaseValsGrphHz[i]->GetPointX(j), phaseValsGrphHz[i]->GetPointY(j) - phaseValsGrphHz[i]->GetMean(2));
		}
	}
	phaseValsGrph->SetTitle("Scarti su misure di fase ripetute - Canale 2");
	phaseValsGrph->GetXaxis()->SetTitle("Frequenza (Hz)");
	phaseValsGrph->GetYaxis()->SetTitle("Fase (scarto) (Rad)");
	phaseValsGrph->SetMarkerStyle(3);
	phaseValsGrph->SetMarkerSize(0.5);
	
	TH1F* phaseValsHisto = new TH1F("phaseValsHisto", "Valori misurati - Fase - Canale 2, 16 kHz", 15, 1.395, 1.413);
	phaseValsHisto->SetFillColor(kBlue);
	for (int i {0}; i < phaseValsGrphHz[2]->GetN(); ++i) {
		phaseValsHisto->Fill(phaseValsGrphHz[2]->GetPointY(i));
	}
	TF1* moregaus = new TF1("moregaus", "gaus", 1.395, 1.4125);
	phaseValsHisto->Fit(moregaus, "Q");

	TCanvas* phaseValsScatterPlot = new TCanvas("phaseValsCnvs", "poopoopeepee");
	phaseValsGrph->Draw("AP");

	TCanvas* phaseValsHistoCnvs = new TCanvas("phaseValsHistoCnvs", "sbor");
	phaseValsHisto->Draw();

	expAmpFreqRespRF->Write();
	expAmpFreqRespLF->Write();
	expAmpFreqRespCF->Write();
	expAmpFreqRespsCvs->Write();
	expPhaseRespGen->Write();
	expPhaseRespR->Write();
	expPhaseRespL->Write();
	expPhaseRespC->Write();
	expPhaseRespsCvs->Write();
	ampValsGrph->Write();
	ampValsSctrPlt->Write();
	ampValsH->Write();
	ampValsHisto->Write();
	phaseValsGrph->Write();
	phaseValsScatterPlot->Write();
	phaseValsHistoCnvs->Write();
}

//******************************************************
//* Jettiness macro calculated with coordinate PxPyPzE
//* Author:Mi Ran Kim
//* Sungkyunkwan University.South Korea
//*  4th.March.2020
//******************************************************** 
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TFrame.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"
#include "TMath.h"
#include "TPaletteAxis.h"
#include "TROOT.h"
#include "TSpline.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include <vector>
#include "TLorentzVector.h"
#include "stdlib.h"

using namespace std;

int compare (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}

void jettiness()
{
  std::ofstream myfile;
 // myfile.open("Jettiness_NewDef_xyzE.txt");

  TFile* inputFile;
  TTree* inputTree;
  TH1F* hCounters;
  int NGenEvt;
  Float_t gen_sumWeights;
  Float_t partialSampleWeight;

  Int_t nRun;
  Long64_t nEvent;
  Int_t nLumi;
  Float_t overallEventWeight;
  Float_t xsec;

  Short_t ZZsel;
  vector<Float_t> *LepPt = 0;
  vector<Float_t> *LepEta = 0;
  vector<Float_t> *LepPhi = 0;
  vector<Float_t> *LepLepId=0;

  vector<Float_t> *ExtraLepPt = 0;
  vector<Float_t> *ExtraLepEta = 0;
  vector<Float_t> *ExtraLepPhi =0;
  vector<Float_t> *ExtraLepLepId =0;
  vector<Float_t> *ExtraLepMass =0;//NON Defined NTUPUL
  
  vector<Float_t> *PhotonPt=0;
  vector<Float_t> *PhotonEta=0;
  vector<Float_t> *PhotonPhi=0;

  Float_t ZZMass;
  vector<Float_t> *fsrPt	=0;
  vector<Float_t> *fsrEta	=0;
  vector<Float_t> *fsrPhi	=0;
  vector<Float_t> *fsrDR	=0;
  vector<Float_t> *fsrLept	=0;
  vector<Float_t> *fsrLeptID	=0;
  vector<Float_t> *fsrGenPt	=0;

  vector<Float_t> *JetPt     = 0;
  vector<Float_t> *JetEta     = 0;
  vector<Float_t> *JetPhi     = 0;
  vector<Float_t> *JetMass     = 0;


  // open input file
  //inputFile =  TFile::Open( "/afs/cern.ch/user/m/mrkim/HZZ/CMSSW_10_2_18/src/ZZAnalysis/AnalysisStep/test/ZZ4lAnalysis.root" );
  //inputFile =  TFile::Open( "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIILegacy/200205_CutBased/MC_2016/ggH125/ZZ4lAnalysis.root" );
  //inputFile =  TFile::Open( "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIILegacy/200205_CutBased/MC_2017/ggH125/ZZ4lAnalysis.root" );
  //  inputFile =  TFile::Open( "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIILegacy/200205_CutBased/MC_2018/ggH125/ZZ4lAnalysis.root" );
//  inputFile =  TFile::Open( "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIILegacy/200205_CutBased/MC_2016/ttH125/ZZ4lAnalysis.root" );
  //inputFile =  TFile::Open( "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIILegacy/200205_CutBased/MC_2017/ttH125/ZZ4lAnalysis.root" );
  inputFile =  TFile::Open( "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIILegacy/200205_CutBased/MC_2018/ttH125/ZZ4lAnalysis.root" );
  float lumi = 59.7; //fb-1

  
  hCounters = (TH1F*)inputFile->Get("ZZTree/Counters");
  NGenEvt = (Float_t)hCounters->GetBinContent(1);
  gen_sumWeights = (Float_t)hCounters->GetBinContent(40);
  partialSampleWeight = lumi * 1000 / gen_sumWeights;
  inputTree = (TTree*)inputFile->Get("ZZTree/candTree");


  inputTree->SetBranchAddress("RunNumber", &nRun);
  inputTree->SetBranchAddress("EventNumber", &nEvent);
  inputTree->SetBranchAddress("LumiNumber", &nLumi);
  inputTree->SetBranchAddress("overallEventWeight", &overallEventWeight);
  inputTree->SetBranchAddress("xsec", &xsec);
  inputTree->SetBranchAddress("ZZsel", &ZZsel);
  inputTree->SetBranchAddress("LepPt", &LepPt);
  inputTree->SetBranchAddress("LepEta", &LepEta);
  inputTree->SetBranchAddress("LepPhi", &LepPhi);
  inputTree->SetBranchAddress("LepLepId", &LepLepId);
  inputTree->SetBranchAddress("ExtraLepPt", &ExtraLepPt);
  inputTree->SetBranchAddress("ExtraLepEta", &ExtraLepEta);
  inputTree->SetBranchAddress("ExtraLepPhi", &ExtraLepPhi);
  inputTree->SetBranchAddress("ExtraLepLepId", &ExtraLepLepId);
  inputTree->SetBranchAddress("ZZMass", &ZZMass);  
  inputTree->SetBranchAddress("JetPt", &JetPt);
  inputTree->SetBranchAddress("JetEta", &JetEta);
  inputTree->SetBranchAddress("JetMass", &JetMass);
  inputTree->SetBranchAddress("JetPhi", &JetPhi);

  TH1F* HiggsMass = new TH1F("HiggsMass","Higgs Mass",100,70.,300.);
  TH1F* JetMass_H = new TH1F("JetMass_H","JetMass_H",100,0.,150.);
  //TH1F* WgtOAEW= new TH1F("WgtOAEW","WgtOAEW",100.,-100,100);  
  TH1F* WgtEW= new TH1F("WgtEW","WgtEW",100.,-10,10);  
//  TH1F* xsec= new TH1F("xsec","xsec",100.,-10,10);  

} //end multiRoot function
 /*
//---------------------------ggH-----------------------------------
 TH1F* JetN = new TH1F("JetN","ggH_JetN",20,0.,20.);//Number of Jets
 Int_t maxJN=14;//Max number of Jet in event
  TString JN[14]={"JetN0","JetN1","JetN2","JetN3","JetN4","JetN5","JetN6","JetN7","JetN8","JetN9","JetN10","JetN11","JetN12","JetN13"};//Number of Jet
  TString Jtt[14]={"Jetti_0","Jetti_1","Jetti_2","Jetti_3","Jetti_4","Jetti_5","Jetti_6","Jetti_7","Jetti_8","Jetti_9","Jetti_10","Jetti_11","Jetti_12","Jetti_13"};//Number of Jettiness
  
  TH1F* Jetti[maxJN];
  for(Int_t k=0;k<maxJN;k++){
 	Jetti[k]= new TH1F("ggH_"+Jtt[k],"ggH_"+Jtt[k],80,-2.,2.);//Number of Jets
  }

  TH1F* JetNJetti[maxJN][maxJN];
  for(Int_t jn=0;jn<maxJN;jn++){
	for(Int_t k=0;k<jn+1;k++){
  	JetNJetti[jn][k]= new TH1F("ggH_"+JN[jn]+"_"+Jtt[k],"ggH_"+JN[jn]+"_"+Jtt[k],80,-2.,2.);//Number of Jets
  	}
  }

  TH1F* Jetti2 = new TH1F("Jetti2","ggH_Jettiness",80,-2.,2.);//Number of Jets
  TH1F* Jetti1 = new TH1F("Jetti1","ggH_Jettiness",80,-2.,2.);//Number of Jets
  TH1F* Jetti0= new TH1F("Jetti0","ggH_Jettiness",320,-8.,8.);//Number of Jets
  TH1F* Jetti2_Jet3 = new TH1F("Jetti2_Jet3","ggH_Jetti2 Jet>2",80,-2.,2.);//Number of Jets
  TH1F* Jetti2_Jet2 = new TH1F("Jetti2_Jet2","ggH_Jetti2 Jet=2",80,-2.,2.);//Number of Jets
  TH1F* Jetti1_Jet2 = new TH1F("Jetti1_Jet2","ggH_Jetti1 Jet=2",80,-2.,2.);//Number of Jets
  TH1F* Jetti0_Jet2 = new TH1F("Jetti0_Jet2","ggH_Jetti0 Jet=2",80,-2.,2.);//Number of Jets
  TH1F* Jetti1_Jet1 = new TH1F("Jetti1_Jet1","ggH_Jetti1 Jet=1",80,-2.,2.);//Number of Jets
  TH1F* Jetti0_Jet1 = new TH1F("Jetti0_Jet1","ggH_Jetti0 Jet=1",80,-2.,2.);//Number of Jets
  TH1F* Jetti0_Jet0 = new TH1F("Jetti0_Jet0","ggH_Jetti0 Jet=0",80,-2.,2.);//Number of Jets
 */

 //---------------------------ttH-----------------------------------
 TH1F* JetN = new TH1F("JetN","ttH_JetN",20,0.,20.);//Number of Jets
 Int_t maxJN=18;
 TString JN[18]={"JetN0","JetN1","JetN2","JetN3","JetN4","JetN5","JetN6","JetN7","JetN8","JetN9","JetN10","JetN11","JetN12","JetN13","JetN14","JetN15","JetN16","JetN17"};//Number of Jet
  TString Jtt[18]={"Jetti_0","Jetti_1","Jetti_2","Jetti_3","Jetti_4","Jetti_5","Jetti_6","Jetti_7","Jetti_8","Jetti_9","Jetti_10","Jetti_11","Jetti_12","Jetti_13","Jetti_14","Jetti_15","Jetti_16","Jetti_17"};//Number of Jettiness
 
  TH1F* Jetti[maxJN];
  for(Int_t k=0;k<maxJN;k++){
 	Jetti[k]= new TH1F("ttH_"+Jtt[k],"ttH_"+Jtt[k],400,-10.,10.);//Number of Jets
  }

  TH1F* JetNJetti[maxJN][maxJN];
  for(Int_t jn=0;jn<maxJN;jn++){
	for(Int_t k=0;k<jn+1;k++){
  	JetNJetti[jn][k]= new TH1F("ttH_"+JN[jn]+"_"+Jtt[k],"ttH_"+JN[jn]+"_"+Jtt[k],400,-10.,10.);//Number of Jets
  	}
  }
   
  Double_t eventWeight = partialSampleWeight * xsec * overallEventWeight ;
//-------------Normal vertor of direction Z to Pt,Eta,Phi coordiatte--------------//
  TLorentzVector v4_na(0.0,0.0,1.0,1.0);//xyzE
  TLorentzVector v4_nb(0.0,0.0,-1.0,1.0);//xy-zE
  //TLorentzVector v4_na(naxyz.Pt(),naxyz.Eta(),naxyz.Phi(),naxyz.M());
  //TLorentzVector v4_nb(nbxyz.Pt(),nbxyz.Eta(),nbxyz.Phi(),nbxyz.M());
  cout<<"v4_na[0]="<<v4_na.Px()<<endl;
  cout<<"v4_na[1]="<<v4_na.Py()<<endl;
  cout<<"v4_na[2]="<<v4_na.Pz()<<endl;
  cout<<"v4_na[3]="<<v4_na.E()<<endl;
  cout<<"v4_nb[0]="<<v4_nb.Px()<<endl;
  cout<<"v4_nb[1]="<<v4_nb.Py()<<endl;
  cout<<"v4_nb[2]="<<v4_nb.Pz()<<endl;
  cout<<"v4_nb[3]="<<v4_nb.E()<<endl;

  Int_t Lep_TotN;//Total Lepton Number of Final states
  Int_t Lep_N;//
  Int_t ExtraLep_N;
  //Int_t Photon_N;
  Int_t Jet_N;
  Float_t e_mass=0.00051;//GeV/c^2
  Float_t mu_mass=0.10566;//Gev/c^2
 // Float_t JetE[Jet_N];
  
  //------------------------- Loop over tree entries--------------------------//  
  Long64_t entries = inputTree->GetEntries();
      
  myfile<<"entry "<<"ExraLep_N "<<"Jet_N  "<<"Jettiness[0]	"<<"Jettiness[1]	"<<"Jettiness[2]	"<<endl;
  for (Long64_t entry= 0; entry < entries; entry++){
  //for (Long64_t entry = 0; entry < 100; entry++) {   
      
      inputTree->GetEntry(entry);
      cout<<"------entry="<<entry<<"-----------"<<endl;
    //  myfile<<"------entry="<<entry<<"-----------"<<endl;
       	if(LepEta->size()!=4){
        cerr << "error in event " << nRun << ":" << nLumi << ":" << nEvent << "; stored leptons: "<< LepEta->size() << endl;
	continue;
	}
		     
        if( !(ZZsel >= 0) ) continue; 
 
      Lep_TotN=0;//Total Lepton Number of Final states
      Lep_N=0;//
      ExtraLep_N=0;
     //Int_t Photon_N;
      Jet_N=0;

      Jet_N=JetMass->size();//# of Jets
      ExtraLep_N=ExtraLepPt->size();//# of Extra_Lepton
      //Photon_N=Photon->size();
      Lep_N=LepPt->size();//#of Lepton* 
      Lep_TotN=Lep_N+ExtraLep_N;//Total # of Lepton in final state
      //cout<<"Lep_TotN="<< Lep_TotN<<endl;//!!BREAK OCCURE!!
      //cout<<"Lep_N="<< Lep_N<<endl;
      cout<<"ExtraLep_N="<<ExtraLep_N<<endl;
      // myfile<<"ExtraLep_N="<<ExtraLep_N<<endl;
      //cout<<"Lep_TotN="<< Lep_TotN<<endl;
      cout<<"Jet_N="<< Jet_N<<endl;
      //myfile<<"Jet_N="<<Jet_N<<endl;

  //-----------Definiton of TLonrenzvector-----------------//
     // TLorentzVector v4_Lep[Lep_TotN]; //To Add Lep_N + ExtraLep_N
      TLorentzVector v4_Lep[Lep_N];//To Add Only Lep_N
      TLorentzVector v4_ExtraLep[ExtraLep_N];
      TLorentzVector v4_LepTot;
      TLorentzVector v4_LepAdd;
      TLorentzVector v4_LepAddTot;//=LepTot+LepAdd
      TLorentzVector v4_Jet[Jet_N];
      TLorentzVector v4_JetNSum[Jet_N];//Sum # of Jet
      TLorentzVector v4_JetTot;//Sum of all Jet
      TLorentzVector v4_Photon;

//------------------------Construction of v4 for Leptons------------------
      for(Int_t i=0;i<Lep_N;i++){
		if(abs(LepLepId->at(i))==11){
		TLorentzVector v4_PEPM;//v4 Container of PtEtaPhiM valueis
		v4_PEPM.SetPtEtaPhiM(LepPt->at(i),LepEta->at(i),LepPhi->at(i),e_mass);

		v4_Lep[i].SetPxPyPzE(v4_PEPM.Px(),v4_PEPM.Py(),v4_PEPM.Pz(),v4_PEPM.E());
		}//End of if(|LepLepId)
		if(abs(LepLepId->at(i))==13){
		TLorentzVector v4_PEPM;//v4 Container of PtEtaPhiM values
		v4_PEPM.SetPtEtaPhiM(LepPt->at(i),LepEta->at(i),LepPhi->at(i),mu_mass);
		
		v4_Lep[i].SetPxPyPzE(v4_PEPM.Px(),v4_PEPM.Py(),v4_PEPM.Pz(),v4_PEPM.E());
		}//End of if(|LepLepId->at(i)
      v4_LepTot +=v4_Lep[i];//Sum of all leptons four momentum=Total Lepton momentum		
      }//End of for(Lep_N)
 
//------------------------Construction of v4 for ExtraLeptons------------------
      for(Int_t j=0;j<ExtraLep_N;j++){
		if(abs(ExtraLepLepId->at(j))==11){
			TLorentzVector v4_PEPM;//v4 Container of PtEtaPhiM values
			v4_PEPM.SetPtEtaPhiM(ExtraLepPt->at(j),ExtraLepEta->at(j),ExtraLepPhi->at(j),e_mass);

			v4_ExtraLep[j].SetPxPyPzE(v4_PEPM.Px(),v4_PEPM.Py(),v4_PEPM.Pz(),v4_PEPM.E());
		}//End of if(ExtraLepLepId==11)
		if(abs(ExtraLepLepId->at(j))==13){
			TLorentzVector v4_PEPM;//v4 Container of PtEtaPhiM valueis
			v4_PEPM.SetPtEtaPhiM(ExtraLepPt->at(j),ExtraLepEta->at(j),ExtraLepPhi->at(j),mu_mass);
			v4_ExtraLep[j].SetPxPyPzE(v4_PEPM.Px(),v4_PEPM.Py(),v4_PEPM.Pz(),v4_PEPM.E());
		}//End of if(ExtraLepLepId==13)
     // v4_LepTot +=v4_ExtraLep[j];//Sum of all leptons four momentum=Total Lepton momentum		
      v4_LepAdd +=v4_ExtraLep[j];//Sum of all leptons four momentum=Total Lepton momentum		
      cout<<"Here is extra Lepton!"<<endl;
       }//End of for(ExtraLep_N)

//-------------------------Construction of v4 for Jets------------------------------
     for(Int_t k=0;k<Jet_N;k++){		
		TLorentzVector v4_PEPM;//v4 Container of PtEtaPhiM values
		v4_PEPM.SetPtEtaPhiM(JetPt->at(k),JetEta->at(k),JetPhi->at(k),JetMass->at(k));
		v4_Jet[k].SetPxPyPzE(v4_PEPM.Px(),v4_PEPM.Py(),v4_PEPM.Pz(),v4_PEPM.E());

	//-----------The sorting Jet Energy---------------//
	//qsort(JetE,sizeof(JetE),sizeof(Float_t),compare)
	//Int_t k=0;
	//Float_t max=JetE[k];
	//for(Int_t j=0;j<Jet_N;j++){
	//		v4_Jet[j].SetPtEtaPhiM(JetPt->at(j),JetEta->at(j),JetPhi->at(j),JetMass->at(j));
//		JetE[j]=v4_Jet[j].Mag();
//			if(JetE[j]>max){
//			max=JetE[j];
//			k=i;
	//	}//End of for(Int_t j=0;j<Jet_N;j++)
	
	v4_JetNSum[k]+=v4_Jet[k];//Sume of #of Jet
	v4_JetTot+=v4_Jet[k];//Sum of all Jet
  
       	cout<<"JetNSumEnergy["<<k<<"]="<<v4_JetNSum[k].E()<<endl;
     	}//End of for(Jet_N)

//---------------Definitions for Tau(Jettiness) calculation-----------------//
	
	Int_t N=2+Jet_N;//#of Tua eliments
	Float_t Tau[Jet_N+1][Jet_N+2];//Tau elements array,
	Float_t TauSum[Jet_N+1];//Summ of mim value of Tau elements
	Float_t QQ[Jet_N+1];//Energy Normalization Factor
     	Float_t Jettiness[Jet_N+1];// Jettiness=TauSun[i]/QQ[i]
	//TLorentzVector v4_LepJetSum;//Four vector sum of JetsN(jN),leptons(l) in signal
     	TLorentzVector v4_LepJetNSum[Jet_N+1];//Four vector sum of JetsN(jN),leptons(l) in signal
     	TLorentzVector v4_LepJetNSumNor[Jet_N+1];//Four vertor sum of all Jets,Leptons in signal for calculation of QQ
     	TLorentzVector v4_qa[Jet_N+1];//Beam=Scalar products Four vector sum of j,l and the direction four vector na
      	TLorentzVector v4_qb[Jet_N+1];//Beam=Scalar products Four vector sum of j,l and the Udirection four vector nb
 
	v4_LepAddTot=v4_LepAdd+v4_LepTot;//USE IT When include additonal leptons      	
	v4_LepJetNSum[0]=v4_LepTot;//All lepton + 0 Jet 
	//v4_LepJetNSum[0]=v4_LepAddTot;//When consider extra Lepton
	for(Int_t i=1;i<Jet_N+1;i++){
		v4_LepJetNSum[i]=v4_LepTot+v4_JetNSum[i-1];//All signal leptons + NJets
		//v4_LepJetNSum[i]=v4_LepAddTot+v4_JetNSum[i-1];//All signal leptons+ Additional lepton + NJets
	}
	//Normalization of Energy, need for calculation of QQ
	v4_LepJetNSumNor[0]=v4_LepTot;//Only Jet=0
	//cout<<"v4_LepJetNSumNor[0]="<<v4_LepTot<<endl;//ERROR MESSAGY
	for(Int_t i=1;i<Jet_N+1;i++){	
		v4_LepJetNSumNor[i]=v4_LepTot+v4_JetNSum[i-1];//All signla leptons +NJets
		//v4_LepJetNSumNor[i]=v4_LepAddTot+v4_JetNSum[i-1];//All signla leptons+ Additional lepton +NJets
	}

	v4_qa[0]=(v4_LepJetNSum[0]*v4_nb)*v4_na;//Beam
        for(Int_t i=1;i<Jet_N+1;i++){
		v4_qa[i]=(v4_LepJetNSum[i]*v4_nb)*v4_na;//
	}
	v4_qb[0]=(v4_LepJetNSum[0]*v4_na)*v4_nb;//Beam
	for(Int_t i=1;i<Jet_N+1;i++){
		v4_qb[i]=(v4_LepJetNSum[i]*v4_na)*v4_nb;
	}

	QQ[0]=(v4_LepJetNSumNor[0]*v4_nb)*(v4_na*v4_LepJetNSumNor[0]);
	cout<< "Q^2[0]="<<QQ[0]<<endl;
	for(Int_t i=1;i<Jet_N+1;i++){
		QQ[i]=(v4_LepJetNSumNor[i]*v4_nb)*(v4_na*v4_LepJetNSumNor[i]);
		cout<< "Q^2["<<i<<"]="<<QQ[i]<<endl;
	}
	
	//Initialization of From TauSum[0] to TauSum[2]
	for(Int_t i=0; i<3;i++){
		TauSum[i]=0.0;
	}//for (Int_t i=0; i<3;i++)

	Jettiness[0]=1000;//Initialize
	Jettiness[1]=1000;//Initialize
	Jettiness[2]=1000;//Initialize

//-----------------(ExtraLep_N=0) Calculation of Jettineis[0] to [2]----------------//
      if(ExtraLep_N==0){
		if(Jet_N==0){ 
			Jettiness[0]=0;
			cout<<"Jettiness[0]="<<Jettiness[0]<<endl;
		//Jetti0_Jet0->Fill(Jettiness[0]);
		JetNJetti[Jet_N][0]->Fill(Jettiness[0]);
		}//End of if(Jet_N==0)
		
		else if(Jet_N==1){
			//--------Jettiniss[0]:No Jet is considerd as signal but Pk(not signals)----
			Tau[0][0]=v4_qa[0]*v4_Jet[0];
			Tau[0][1]=v4_qb[0]*v4_Jet[0];
			cout<<"Tau[0][0]="<<Tau[0][0]<<endl;
			cout<<"Tau[0][1]="<<Tau[0][1]<<endl;
			Float_t min;
			min=Tau[0][0];
			//----------Find min element of Tau-------------
				for(Int_t m=0;m<2;m++){
				if(Tau[0][m]<min) min=Tau[0][m];
				}//End of for i(Int_t m=0;m<N;m++)
			cout<<"min value[0]="<<min<<endl;
			TauSum[0] += min;//Sum of min tau elements
			cout<<"Sun of min value="<<TauSum[0]<<endl;
			
			Jettiness[0]=TauSum[0]/QQ[0];
			cout<<"Q^2["<<0<<"]="<<QQ[0]<<endl;
			//myfile<<"Jettiness["<<0<<"]="<<Jettiness[0]<<endl;
			cout<<"Jettiness[0]="<<Jettiness[0]<<endl;
			//Jetti0_Jet1->Fill(Jettiness[0],);//Jet Number  
			JetNJetti[Jet_N][0]->Fill(Jettiness[0]);
			
			//-----------Jettiness[1]-------------------
			TauSum[1]=0;//Becaus no particles can be considered as Pk
			Jettiness[1]=TauSum[1]/QQ[1];
			cout<<"Q^2["<<1<<"]="<<QQ[1]<<endl;
			//myfile<<"Jettiness["<<1<<"]="<<Jettiness[1]<<endl;
			cout<<"Jettiness[1]="<<Jettiness[1]<<endl;		
      			//Jetti1_Jet1->Fill(Jettiness[1],);//Jet Number  
			JetNJetti[Jet_N][1]->Fill(Jettiness[1]);
		}//End of else if(Jet_N==1)

		else if(Jet_N>=2){
			//-----------Jettiness[0] No Jets are considerd as signals but as Pk-----------------
			for(Int_t j=0;j<Jet_N;j++){
				Tau[0][0]=v4_qa[0]*v4_Jet[j];
				Tau[0][1]=v4_qb[0]*v4_Jet[j];
				cout<<"Tau[0][0]="<<Tau[0][0]<<endl;
				cout<<"Tau[0][1]="<<Tau[0][1]<<endl;
				Float_t min1;
				min1=Tau[0][0];
				//----------Find min element of Tau---------
					for(Int_t m=0;m<2;m++){
					if(Tau[0][m]<min1) min1=Tau[0][m];
					}//End of for i(Int_t m=0;m<N;m++)
				cout<<"min value[0]="<<min1<<endl;
				TauSum[0] += min1;//Sum of min tau elements
				cout<<"Sun of min value="<<TauSum[0]<<endl;
			}//End of for (Int_t j=0;j<Jet_N;j++)i
			Jettiness[0]=TauSum[0]/QQ[0];
			cout<<"Q^2["<<0<<"]="<<QQ[0]<<endl;
			//myfile<<"Jettiness["<<0<<"]="<<Jettiness[0]<<endl;
			cout<<"Jettiness[0]="<<Jettiness[0]<<endl;
			
		      	//Jetti0_Jet2->Fill(Jettiness[0],);//Jet Number  
			JetNJetti[Jet_N][0]->Fill(Jettiness[0]);

	
		//------------Jettiness[1]:One Jet is considerd as signal and others as Pk-------------------
			for(Int_t j=1;j<Jet_N;j++){
				Tau[1][0]=v4_qa[1]*v4_Jet[j];
				Tau[1][1]=v4_qb[1]*v4_Jet[j];
				Tau[1][2]=v4_Jet[0]*v4_Jet[j];
				cout<<"Tau[1][0]="<<Tau[1][0]<<endl;
				cout<<"Tau[1][1]="<<Tau[1][1]<<endl;
				cout<<"Tau[1][2]="<<Tau[1][2]<<endl;
				Float_t min2;
				min2=Tau[1][0];
				//----------Find min element of Tau-------------
					for(Int_t m=0;m<3;m++){
						if(Tau[1][m]<min2) min2=Tau[1][m];
					}//End of for i(Int_t m=0;m<3;m++)
				cout<<"min value[1]="<<min2<<endl;
				TauSum[1] += min2;//Sum of min tau elements
				cout<<"Sum of min value="<<TauSum[1]<<endl;
			}//for (Int_t j=1;j<Jet_N;j++)
			Jettiness[1]=TauSum[1]/QQ[1];
			cout<<"Q^2["<<1<<"]="<<QQ[1]<<endl;
			//myfile<<"Jettiness["<<1<<"]="<<TauSum[1]<<endl;
			cout<<"Jettiness[1]="<<Jettiness[1]<<endl;
      			//Jetti1_Jet2->Fill(Jettiness[1],);//Jet Number
			JetNJetti[Jet_N][1]->Fill(Jettiness[1]);

			if(Jet_N==2){
				//-----------Jettiness[2]-------------------
				TauSum[2]=0;//Becaus no particles can be considered as Pk
				Jettiness[2]=TauSum[2]/QQ[2];
				cout<<"Q^2["<<2<<"]="<<QQ[2]<<endl;
				//myfile<<"Jettiness["<<2<<"]="<<Jettiness[2]<<endl;
				cout<<"Jettiness[2]="<<Jettiness[2]<<endl;		
      				//Jetti2_Jet2->Fill(Jettiness[2],);//Jet Number
				JetNJetti[Jet_N][2]->Fill(Jettiness[2]);
      			}//End of If(JetN==2)
  
			if(Jet_N>2){

			//---Jettiness[2]:Tow Jets are considered as signals and others as Pk-----
			for(Int_t jn=2;jn<Jet_N+1;jn++){
				for(Int_t j=jn;j<Jet_N;j++){
					Tau[jn][0]=v4_qa[jn]*v4_Jet[j];//v4_qa[1]->v4_qa[jn]:Jet_N=jn +total signal lepton
					Tau[jn][1]=v4_qb[jn]*v4_Jet[j];//
					cout<<"Tau["<<jn<<"][0]="<<Tau[jn][0]<<endl;
					cout<<"Tau["<<jn<<"][1]="<<Tau[jn][1]<<endl;
					for(Int_t k=0;k<jn;k++){
					Tau[jn][k+2]=v4_Jet[k]*v4_Jet[j];
					cout<<"Tau["<<jn<<"]["<<k+2<<"]="<<Tau[jn][k+2]<<endl;
					}
				Float_t min3;
				min3=Tau[jn][0];
				//----------Find min element of Tau-------------
				for(Int_t m=0;m<jn+2;m++){
					if(Tau[jn][m]<min3) min3=Tau[jn][m];
				}//End of for i(Int_t m=0;m<jn+2;m++)
				cout<<"min value["<<jn<<"]="<<min3<<endl;
				TauSum[jn] += min3;//Sum of min tau elements
			 	cout<<"The sums of min values["<<jn<<"]="<<TauSum[jn]<<endl;
				}//End of for(Int_t j=jn;j<Jet_N;j++)
			if(jn==Jet_N) TauSum[jn]=0.;
			Jettiness[jn]=TauSum[jn]/QQ[jn];
			cout<<"Q^2["<<jn<<"]="<<QQ[jn]<<endl;
			//myfile<<"Jettiness["<<2<<"]="<<Jettiness[2]<<endl;
			cout<<"Jettiness["<<jn<<"]="<<Jettiness[jn]<<endl;
      			
	   		//Jetti2->Fill(Jettiness[jn],);//FUNZIONA  
			JetNJetti[Jet_N][jn]->Fill(Jettiness[jn]);//*** Break *** segmentation violatio
			}//End of for(Int_r jn=2;jn<Jet_N+1;jn++)
  			
		}//End of if(Jet_N>2)
	}//End of else if(JetN>=2)		
    }//End of (ExtraLep_N==0)


//-----------------(ExtraLep_N!=0) Calculation of Jettineis[0] to [2]----------------//
      if(ExtraLep_N!=0){
		if(Jet_N==0){ //---Jettiness[0]:v4_ExtraLep is considered as Pk--------
			for(Int_t j=0;j<ExtraLep_N;j++){
			Tau[0][0]=v4_qa[0]*v4_ExtraLep[j];
			Tau[0][1]=v4_qb[0]*v4_ExtraLep[j];
			cout<<"Tau[0][0]="<<Tau[0][0]<<endl;
			cout<<"Tau[0][1]="<<Tau[0][1]<<endl;
			Float_t min;
			min=Tau[0][0];
			//----------Find min element of Tau-------------
				for(Int_t m=0;m<2;m++){
					if(Tau[0][m]<min) min=Tau[0][m];
				}//End of for i(Int_t m=0;m<N;m++)
			cout<<"min value[0]="<<min<<endl;
			TauSum[0] += min;//Sum of min tau elements
			cout<<"The sums of min values="<<TauSum[0]<<endl;
			}//End of for(Int_t j=0;j<ExtraLep_N;j++)
		Jettiness[0]=TauSum[0]/QQ[0];
		cout<<"Q^2["<<0<<"]="<<QQ[0]<<endl;
		//myfile<<"Jettiness["<<0<<"]="<<Jettiness[0]<<endl;
		cout<<"Jettiness[0]="<<Jettiness[0]<<endl;
      		//Jetti0_Jet0->Fill(Jettiness[0],);//Jet Number  
		JetNJetti[Jet_N][0]->Fill(Jettiness[0]);
		}//End of if(Jet_N==0)
		
		else if(Jet_N==1){
			//---------Jettiness[0]:v4_Jet[] is considered as Pk----
			Tau[0][0]=v4_qa[0]*v4_Jet[0];
			Tau[0][1]=v4_qb[0]*v4_Jet[0];
			cout<<"Tau[0][0]="<<Tau[0][0]<<endl;
			cout<<"Tau[0][1]="<<Tau[0][1]<<endl;
			Float_t min0;
			min0=Tau[0][0];
			//----------Find min element of Tau-------------
				for(Int_t m=0;m<2;m++){
				if(Tau[0][m]<min0) min0=Tau[0][m];
				}//End of for i(Int_t m=0;m<N;m++)
			cout<<"min value[0]="<<min0<<endl;
			TauSum[0] += min0;//Sum of min tau elements
			cout<<"The sums of min values="<<TauSum[0]<<endl;

			//--------Jettiness[0]:v4_ExtraLep are considered as Pk-------
			for(Int_t j=0;j<ExtraLep_N;j++){
				Tau[0][0]=v4_qa[0]*v4_ExtraLep[0];
				Tau[0][1]=v4_qb[0]*v4_ExtraLep[0];
				cout<<"Tau[0][0]="<<Tau[0][0]<<endl;
				cout<<"Tau[0][1]="<<Tau[0][1]<<endl;
				Float_t min1;
				min1=Tau[0][0];
			//----------Find min element of Tau-------------
					for(Int_t m=0;m<2;m++){
					if(Tau[0][m]<min1) min1=Tau[0][m];
					}//End of for i(Int_t m=0;m<N;m++)
			cout<<"min value[0]="<<min1<<endl;
			TauSum[0] += min1;//Sum of min tau elements
			cout<<"The sums of min values="<<TauSum[0]<<endl;
			}//End of for (Int_t j=0;j<ExtraLep_N;j++)

			Jettiness[0]=TauSum[0]/QQ[0];
			cout<<"Q^2["<<0<<"]="<<QQ[0]<<endl;
			//myfile<<"Jettiness["<<0<<"]="<<Jettiness[0]<<endl;
			cout<<"Jettiness[0]="<<Jettiness[0]<<endl;
      			//Jetti0_Jet1->Fill(Jettiness[0],);//Jet Number  
			JetNJetti[Jet_N][0]->Fill(Jettiness[0]);
			
			//-------Jettiness[1]:v4_Jet is considered as signal and v4_ExtraLep[] as Pk-------
			for(Int_t j=0;j<ExtraLep_N;j++){
				Tau[1][0]=v4_qa[1]*v4_ExtraLep[j];
				Tau[1][1]=v4_qb[1]*v4_ExtraLep[j];
				Tau[1][2]=v4_Jet[0]*v4_ExtraLep[j];
				cout<<"Tau[1][0]="<<Tau[1][0]<<endl;
				cout<<"Tau[1][1]="<<Tau[1][1]<<endl;
				cout<<"Tau[1][2]="<<Tau[1][2]<<endl;
				Float_t min2;
				min2=Tau[1][0];
			//----------Find min element of Tau-------------
					for(Int_t m=0;m<3;m++){
					if(Tau[1][m]<min2) min2=Tau[1][m];
					}//End of for i(Int_t m=0;m<N;m++)
				cout<<"min value[1]="<<min2<<endl;
				TauSum[1] += min2;//Sum of min tau elements
			 	cout<<"The sums of min values="<<TauSum[1]<<endl;
			}//End of for (Int_t j=0;j<ExtraLep_N;j++)
			Jettiness[1]=TauSum[1]/QQ[1];
			cout<<"Q^2["<<1<<"]="<<QQ[1]<<endl;
			//myfile<<"Jettiness["<<1<<"]="<<Jettiness[1]<<endl;
			cout<<"Jettiness[1]="<<Jettiness[1]<<endl;
      			//Jetti1_Jet1->Fill(Jettiness[1],);//Jet Number  
			JetNJetti[Jet_N][1]->Fill(Jettiness[1]);
		}//End of else if(Jet_N==1)

		else if(Jet_N>=2){
			//-----------Jettiness[0]:v4_Jet are considered as Pk-----------------
			for(Int_t j=0;j<Jet_N;j++){
				Tau[0][0]=v4_qa[0]*v4_Jet[j];
				Tau[0][1]=v4_qb[0]*v4_Jet[j];
				cout<<"Tau[0][0]="<<Tau[0][0]<<endl;
				cout<<"Tau[0][1]="<<Tau[0][1]<<endl;
				Float_t min1;
				min1=Tau[0][0];
				//----------Find min element of Tau---------
					for(Int_t m=0;m<2;m++){
					if(Tau[0][m]<min1) min1=Tau[0][m];
					}//End of for i(Int_t m=0;m<N;m++)
				cout<<"min value[0]="<<min1<<endl;
				TauSum[0] += min1;//Sum of min tau elements
			 	cout<<"The sums of min values="<<TauSum[0]<<endl;
			}//End of for (Int_t j=0;j<Jet_N;j++)
			
			//-----------Jettiness[0] v4_extraLep are considred as Pk-----------------
			for(Int_t j=0;j<ExtraLep_N;j++){ 
				Tau[0][0]=v4_qa[0]*v4_ExtraLep[j];
				Tau[0][1]=v4_qb[0]*v4_ExtraLep[j];
				cout<<"Tau[0][0]="<<Tau[0][0]<<endl;
				cout<<"Tau[0][1]="<<Tau[0][1]<<endl;
				Float_t min2;
				min2=Tau[0][0];
				//----------Find min element of Tau---------
					for(Int_t m=0;m<2;m++){
					if(Tau[0][m]<min2) min2=Tau[0][m];
					}//End of for i(Int_t m=0;m<N;m++)
				cout<<"min value[0]="<<min2<<endl;
				TauSum[0] += min2;//Sum of min tau elements
			 	cout<<"The sums of min values="<<TauSum[0]<<endl;
			}//End of for(Int_t j=0;j<ExtraLep_N;j++)
			Jettiness[0]=TauSum[0]/QQ[0];
			cout<<"Q^2["<<0<<"]="<<QQ[0]<<endl;
			//myfile<<"Jettiness["<<0<<"]="<<Jettiness[0]<<endl;
			cout<<"Jettiness[0]="<<Jettiness[0]<<endl;
			
      			//Jetti0_Jet2->Fill(Jettiness[0],);//Jet Number 
			JetNJetti[Jet_N][0]->Fill(Jettiness[0]);
	 
			//------------Jettiness[1]:One Jet as signal, lest jet as Pk----------
			for(Int_t j=1;j<Jet_N;j++){
				Tau[1][0]=v4_qa[1]*v4_Jet[j];
				Tau[1][1]=v4_qb[1]*v4_Jet[j];
				Tau[1][2]=v4_Jet[0]*v4_Jet[j];
				cout<<"Tau[1][0]="<<Tau[1][0]<<endl;
				cout<<"Tau[1][1]="<<Tau[1][1]<<endl;
				cout<<"Tau[1][2]="<<Tau[1][2]<<endl;
				Float_t min3;
				min3=Tau[1][0];
				//----------Find min element of Tau-------------
					for(Int_t m=0;m<3;m++){
					if(Tau[1][m]<min3) min3=Tau[1][m];
					}//End of for i(Int_t m=0;m<3;m++)
				cout<<"min value[0]="<<min3<<endl;
				TauSum[1] += min3;//Sum of min tau elements
			 	cout<<"The sums of min values="<<TauSum[1]<<endl;
			}//End of for(Int_t j=1;j<Jet_N;j++)i

			//------------Jettiness[1]:One Jet as signal, v4_ExtraLep as Pk-------------------
			for(Int_t j=0;j<ExtraLep_N;j++){
				Tau[1][0]=v4_qa[1]*v4_ExtraLep[j];
				Tau[1][1]=v4_qb[1]*v4_ExtraLep[j];
				Tau[1][2]=v4_Jet[0]*v4_ExtraLep[j];
				cout<<"Tau[1][0]="<<Tau[1][0]<<endl;
				cout<<"Tau[1][1]="<<Tau[1][1]<<endl;
				cout<<"Tau[1][2]="<<Tau[1][2]<<endl;
				Float_t min4;
				min4=Tau[1][0];
				//----------Find min element of Tau-------------
					for(Int_t m=0;m<3;m++){
					if(Tau[1][m]<min4) min4=Tau[1][m];
					}//End of for i(Int_t m=0;m<3;m++)
				cout<<"min value[1]="<<min4<<endl;
				TauSum[1] += min4;//Sum of min tau elements
			 	cout<<"The sums of min values="<<TauSum[1]<<endl;
				}//End of for  (Int_t j=1;j<ExtraLep_N;j++)

			Jettiness[1]=TauSum[1]/QQ[1];
			cout<<"Q^2["<<1<<"]="<<QQ[1]<<endl;
			//myfile<<"Jettiness["<<1<<"]="<<Jettiness[1]<<endl;
			cout<<"Jettiness[1]="<<Jettiness[1]<<endl;
      			//Jetti1_Jet2->Fill(Jettiness[1],);//Jet Number  
			JetNJetti[Jet_N][1]->Fill(Jettiness[1]);
			
			if(Jet_N==2){
				//-------------Jettiness[2]-----------------
				for(Int_t j=0;j<ExtraLep_N;j++){
					Tau[2][0]=v4_qa[2]*v4_ExtraLep[j];
					Tau[2][1]=v4_qb[2]*v4_ExtraLep[j];
					Tau[2][2]=v4_Jet[0]*v4_ExtraLep[j];
					Tau[2][3]=v4_Jet[1]*v4_ExtraLep[j];
					cout<<"Tau[2][0]="<<Tau[2][0]<<endl;
					cout<<"Tau[2][1]="<<Tau[2][1]<<endl;
					cout<<"Tau[2][2]="<<Tau[2][2]<<endl;
					cout<<"Tau[2][3]="<<Tau[2][3]<<endl;
					Float_t min3;
					min3=Tau[2][0];
					//----------Find min element of Tau-------------
						for(Int_t m=0;m<4;m++){
							if(Tau[2][m]<min3) min3=Tau[2][m];
						}//End of for i(Int_t m=0;m<4;m++)
					cout<<"min value[2]="<<min3<<endl;
					TauSum[2] += min3;//Sum of min tau elements
			 		cout<<"The sums of min values="<<TauSum[2]<<endl;
				}//End of for(Int_t j=0;j<Extra_N;j++)
				
				Jettiness[2]=TauSum[2]/QQ[2];
				cout<<"Q^2["<<2<<"]="<<QQ[2]<<endl;
				//myfile<<"Jettiness["<<2<<"]="<<Jettiness[2]<<endl;
				cout<<"Jettiness[2]="<<Jettiness[2]<<endl;
				//Jetti2_Jet2->Fill(Jettiness[2],);	
				JetNJetti[Jet_N][2]->Fill(Jettiness[2]);
			}//Endo of(Jet_N==2)
			
			if(Jet_N>2){
				//---------Jettiness[2]:Two jets are considred as signals and others as Pk-----
				for(Int_t jn=2;jn<Jet_N+1;jn++){
					for(Int_t j=jn;j<Jet_N;j++){
						Tau[jn][0]=v4_qa[jn]*v4_Jet[j];
						Tau[jn][1]=v4_qb[jn]*v4_Jet[j];
						cout<<"Tau["<<jn<<"][0]="<<Tau[jn][0]<<endl;
						cout<<"Tau["<<jn<<"][1]="<<Tau[jn][1]<<endl;
						for(Int_t k=0;k<jn;k++){
							Tau[jn][k+2]=v4_Jet[k]*v4_Jet[j];
							cout<<"Tau["<<jn<<"]["<<k+2<<"]="<<Tau[jn][k+2]<<endl;
						}//End of for(Int_t j=jn;j<Jet_N;j++)
	
					Float_t min3;
					min3=Tau[jn][0];
					//----------Find min element of Tau-------------
					for(Int_t m=0;m<jn+2;m++){
						if(Tau[jn][m]<min3) min3=Tau[jn][m];
					}//End of for i(Int_t m=0;m<jn+2;m++)
					cout<<"min value["<<jn<<"]="<<min3<<endl;
					TauSum[jn] += min3;//Sum of min tau elements
		 			if(jn==Jet_N) TauSum[jn]=0.;
			 		cout<<"The sums of min values["<<jn<<"]="<<TauSum[jn]<<endl;
				}//End of for(Int_t j=2;j<Jet_N;j++)
					
				//-------Jettiness[2]:Two jets are considered as signals and ExtraLep as Pk-----
				for(Int_t j=jn;j<ExtraLep_N;j++){
					Tau[jn][0]=v4_qa[jn]*v4_ExtraLep[j];
					Tau[jn][1]=v4_qb[jn]*v4_ExtraLep[j];
					cout<<"Tau["<<jn<<"][0]="<<Tau[jn][0]<<endl;
					cout<<"Tau["<<jn<<"][1]="<<Tau[jn][1]<<endl;
					for(Int_t k=0;k<jn;k++){
						Tau[jn][k+2]=v4_Jet[k]*v4_ExtraLep[j];
						cout<<"Tau["<<jn<<"]["<<k+2<<"]="<<Tau[jn][k+2]<<endl;
					}//End of for(Int_t k=0;k<jn;k++)

				Float_t min3;
				min3=Tau[jn][0];
				//----------Find min element of Tau-------------
					for(Int_t m=0;m<jn+2;m++){
						if(Tau[jn][m]<min3) min3=Tau[jn][m];
					}//End of for i(Int_t m=0;m<jn+2;m++)
				cout<<"min value["<<jn<<"]="<<min3<<endl;
				TauSum[jn] += min3;//Sum of min tau elements
			 	cout<<"The sums of min values="<<TauSum[jn]<<endl;
				}//End of for(Int_t j=0;j<Extra_N;j++)
			
			Jettiness[jn]=TauSum[jn]/QQ[jn];
			cout<<"Q^2["<<jn<<"]="<<QQ[jn]<<endl;
			//myfile<<"Jettiness["<<2<<"]="<<Jettiness[2]<<endl;
			cout<<"Jettiness["<<jn<<"]="<<Jettiness[jn]<<endl;
      			JetNJetti[Jet_N][jn]->Fill(Jettiness[jn]);//Jet Number
      			
			}//End of for(jn=2;jn<JetN;jn++)  
		}//End of If(Jet_N>2)		
	}//End of else if(Jet_N>=2)
     }//End of (ExtraLep_N!==0)

      //myfile<< entry<<"	"<<ExtraLep_N <<" "<<Jet_N <<"  "<<Jettiness[0] <<"		"<< Jettiness[1]<<"		"<< Jettiness[2]<<endl;

      //Int_t JetM_size=JetMass->size();
     //Int_t LepPt_size=LepPt->size();

      	//Jetti0->Fill(Jettiness[0],);//Jet Number  
      	//Jetti1->Fill(Jettiness[1],);//Jet Number  
      	//Jetti2->Fill(Jettiness[2],);//Jet Number  
        JetN->Fill(Jet_N);//Jet Number 
        WgtEW->Fill(eventWeight);
	cout<<"xsec="<<xsec<<endl;
	cout<<"overAllEventWeight="<<overallEventWeight<<end;   
	for(Int_t j=0;j<Jet_N+1;j++){
		Jetti[j]->Fill(Jettiness[j]);
	}//for(Int_t j=0;j<Jet_N+1;j++)
 
      //HiggsMass->Fill(ZZMass,);
      //for(UInt_t i=0;i<JetMass->size();i++){
	//JetMass_H->Fill(JetMass->at(i),); 
       //cout<<"JetM_size="<< JetM_size<<endl;
     //}//Endof for(UInt_t i=0;i<JetMass->size();i++
    } //end of Loop entry
    //myfile.close();

    
 // TFile *f= new TFile ("ggH3W_JettiN.root","CREATE");//W=weight, 1:eventWeight, 2:partialSampleWeight,3:gen_sumWeights
  //TFile *f= new TFile ("ttHW_JettiN.root","CREATE");//W=weight, 1:eventWeight, 2:partialSampleWeight,3:gen_sumWeights
   // TFile *f= new TFile("Weight.root","CREATE");
  TCanvas *c1= new TCanvas();  
  WgtEW->Draw();
  WgtEW->SaveAs("eventWeight.png");
  //JetN->Write();
 // Jetti2->Write();
  //Jetti1->Write();
  //Jetti0->Write();
  //Jetti2_Jet3->Write();
  //Jetti2_Jet2->Write();
  //Jetti1_Jet2->Write();
  //Jetti0_Jet2->Write();
  //Jetti1_Jet1->Write();
  //Jetti0_Jet1->Write();
  //Jetti0_Jet0->Write();

  for(Int_t j=0;j<maxJN;j++){
	Jetti[j]->Write();
  }//End of for(Int_t j=0;j<maxJN;j++)

  for(Int_t j=0;j<maxJN;j++){
	for(Int_t k=0;k<j+1;k++){
		JetNJetti[j][k]->Write();
        }//End of for(Int_t k=0;k<j+1;k++)
  }//for(Int_t j=0;j<maxJN;j++)
 
  //f->Write();
  //f->Close();

} //end histo funtion Void

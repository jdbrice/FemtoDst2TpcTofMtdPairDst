#ifndef PAIR_DST_MAKER_H
#define PAIR_DST_MAKER_H


#include "TreeAnalyzer.h"
#include "XmlRange.h"

#include "vendor/loguru.h"

#include "TNamed.h"
#include "TTree.h"


#include "FemtoDstFormat/BranchReader.h"
#include "FemtoDstFormat/BranchWriter.h"
#include "FemtoDstFormat/TClonesArrayReader.h"
#include "FemtoDstFormat/FemtoEvent.h"

#include "Filters/MuonMVAFilter.h"
#include "Filters/TrackFilter.h"
#include "Filters/LowPtMuonFilter.h"

#include "PairDstFormat/FemtoPair.h"


template <>
TString XmlConfig::get<TString>( string path ) const {
	TString r( getString( path ) );
	return r;
}

template <>
TString XmlConfig::get<TString>( string path, TString dv ) const {
	if ( !exists( path ) )
		return dv;
	TString r( getString( path ) );
	return r;
}

class PairDstMaker : public TreeAnalyzer
{
protected:

	FemtoEvent *_event;

	BranchReader<FemtoEvent> _rEvent;
	TClonesArrayReader<FemtoTrack> _rTracks;
	TClonesArrayReader<FemtoTrackHelix> _rHelices;
	TClonesArrayReader<FemtoBTofPidTraits> _rBTofPid;
	TClonesArrayReader<FemtoMtdPidTraits> _rMtdPid;

	TrackFilter _lowPtTrackFilter;
	LowPtMuonFilter _lowFilter;

	vector<FemtoTrackProxy> pos_mtd;
	vector<FemtoTrackProxy> neg_mtd;
	vector<FemtoTrackProxy> pos_tof;
	vector<FemtoTrackProxy> neg_tof;

	//

	// BranchReader<FemtoEvent> _fer;
	// TClonesArrayReader<FemtoTrack> _ftr;
	// TClonesArrayReader<FemtoMtdPidTraits> _fmtdr;

	vector <MuonMVAFilter> mlps;
	vector <MuonMVAFilter> bdts;
	HistoBins bins_pt_mva;


	TTree * _pairDst = nullptr;
	BranchWriter<FemtoPair> _fpw;
	FemtoPair _pair; 

	int cenBinMin; //inclusive
	int cenBinMax; //inclusive

public:

	const int DEBUG = 1;
	PairDstMaker() {}
	~PairDstMaker() {}

	virtual void initialize(){
		TreeAnalyzer::initialize();

		this->_rEvent.setup   ( chain, "Event" );
		this->_rTracks.setup  ( chain, "Tracks" );
		this->_rHelices.setup ( chain, "Helices" );
		this->_rBTofPid.setup ( chain, "BTofPidTraits" );
		this->_rMtdPid.setup  ( chain, "MtdPidTraits" );

		_lowPtTrackFilter.load( config, nodePath + ".LowPtTrackFilter" );
		_lowFilter.load( config, nodePath + ".LowPtMuonFilter" );
		

		bins_pt_mva.load( config, nodePath + ".MuonMVAFilter.ptBins" );
		string mlp_template_str = config.getString( nodePath + ".MuonMVAFilter.WeightsMLP" );
		string bdt_template_str = config.getString( nodePath + ".MuonMVAFilter.WeightsBDT" );

		/* LOAD THE MLPs */
		// use this first one to setup the sinlgeton reader and load the variables
		MuonMVAFilter m_setup;
		m_setup.loadVars( config, nodePath + ".MuonMVAFilter" );

		if ( mlp_template_str.find( "%" ) != string::npos ) {
			for ( size_t i = 0; i < bins_pt_mva.nBins(); i++ ){
				MuonMVAFilter m;
				TString wf = TString::Format( mlp_template_str.c_str(), i );
				TString n = TString::Format( "mlp_%zu", i );
				m.load( string( wf ), string( n ) );
				LOG_F( INFO, "Loading %s from %s", n.Data(), wf.Data() );
				mlps.push_back( m );
			}
		} else {
			MuonMVAFilter m;
			m.load( mlp_template_str, "DNN" );
			mlps.push_back( m );
		}

		book->cd();

		this->_pairDst = new TTree( "PairDst", "" );
		this->_fpw.createBranch( this->_pairDst, "Pairs" );

		cenBinMin = config.getInt( nodePath + ".CenBin:min", -1 );
		cenBinMax = config.getInt( nodePath + ".CenBin:max", -1 );

	}
protected:

	void addPair(FemtoEvent * event, FemtoTrackProxy &p1, FemtoTrackProxy &p2){
		
		if ( p1._track->mId == p2._track->mId) return;

		this->_pair.reset();
		this->_pair.mVertexZ              = event->mPrimaryVertex_mX3;
		this->_pair.mDeltaVertexZ         = event->mPrimaryVertex_mX3 - event->mWeight;
		this->_pair.mGRefMult             = event->mGRefMult;



		if ( nullptr != p1._mtdPid ){
			this->_pair.d1_mDeltaY            = p1._mtdPid->mDeltaY;
			this->_pair.d1_mDeltaZ            = p1._mtdPid->mDeltaZ;
			this->_pair.d1_mDeltaTimeOfFlight = p1._mtdPid->mDeltaTimeOfFlight;
			this->_pair.d1_mMatchFlag         = p1._mtdPid->mMatchFlag;
			this->_pair.d1_mTriggerFlag       = p1._mtdPid->mTriggerFlag;
			this->_pair.d1_mCell              = p1._mtdPid->cell();
			this->_pair.d1_mModule            = p1._mtdPid->module();
			this->_pair.d1_mBackleg           = p1._mtdPid->backleg();
		} else {
			// store the TPC TOF muon info instead
			this->_pair.d1_mDeltaTimeOfFlight = _lowFilter.zb( p1, "mu" );
			this->_pair.d1_mDeltaY            = _lowFilter.zbCorr( p1, "mu" );
			this->_pair.d1_mDeltaZ            = _lowFilter.zbCorr( p1, "pi" );
			this->_pair.d1_mCell              = -1;
		}
		
		this->_pair.d1_mPt                = p1._track->mPt;
		this->_pair.d1_mEta               = p1._track->mEta;
		this->_pair.d1_mPhi               = p1._track->mPhi;
		this->_pair.d1_mId                = p1._track->mId;
		this->_pair.d1_mNHitsFit          = p1._track->mNHitsFit;
		this->_pair.d1_mNHitsMax          = p1._track->mNHitsMax;
		this->_pair.d1_mNHitsDedx         = p1._track->mNHitsDedx;
		this->_pair.d1_mNSigmaPion        = p1._track->nSigmaPion();
		this->_pair.d1_mDCA               = p1._track->gDCA();
		this->_pair.d1_mPid               = p1._pid;

		LOG_F( 9, "p2._track = %p", p2._track );
		LOG_F( 9, "p2._mtdPid = %p", p2._mtdPid );
		LOG_F( 9, "p2._btofPid = %p", p2._btofPid );
		// return;

		if ( nullptr != p2._mtdPid ){
			LOG_F( 9, "MTD info" );
			this->_pair.d2_mDeltaY            = p2._mtdPid->mDeltaY;
			this->_pair.d2_mDeltaZ            = p2._mtdPid->mDeltaZ;
			this->_pair.d2_mDeltaTimeOfFlight = p2._mtdPid->mDeltaTimeOfFlight;
			this->_pair.d2_mMatchFlag         = p2._mtdPid->mMatchFlag;
			this->_pair.d2_mTriggerFlag       = p2._mtdPid->mTriggerFlag;
			this->_pair.d2_mCell              = p2._mtdPid->cell();
			this->_pair.d2_mModule            = p2._mtdPid->module();
			this->_pair.d2_mBackleg           = p2._mtdPid->backleg();
		} else {
			LOG_F( 9, "BTOF info" );
			// store the TPC TOF muon info instead
			this->_pair.d2_mDeltaTimeOfFlight = _lowFilter.zb( p2, "mu" );
			this->_pair.d2_mDeltaY            = _lowFilter.zbCorr( p2, "mu" );
			this->_pair.d2_mDeltaZ            = _lowFilter.zbCorr( p2, "pi" );
			this->_pair.d2_mCell              = -1;
		}
		
		this->_pair.d2_mPt                = p2._track->mPt;
		this->_pair.d2_mEta               = p2._track->mEta;
		this->_pair.d2_mPhi               = p2._track->mPhi;
		this->_pair.d2_mId                = p2._track->mId;
		this->_pair.d2_mNHitsFit          = p2._track->mNHitsFit;
		this->_pair.d2_mNHitsMax          = p2._track->mNHitsMax;
		this->_pair.d2_mNHitsDedx         = p2._track->mNHitsDedx;
		this->_pair.d2_mNSigmaPion        = p2._track->nSigmaPion();
		this->_pair.d2_mDCA               = p2._track->gDCA();
		this->_pair.d2_mPid               = p2._pid;


		TLorentzVector lv1, lv2, lv;

		lv = p1._track->lv( 0.1056583745 ) + p2._track->lv( 0.1056583745 ); 
		
		this->_pair.mChargeSum = p1._track->charge() + p2._track->charge();
		this->_pair.mMass      = lv.M();
		this->_pair.mPt        = lv.Pt();
		this->_pair.mEta       = lv.PseudoRapidity();
		this->_pair.mPhi       = lv.Phi();
		this->_pair.mRapidity  = lv.Rapidity();

		this->_fpw.set( this->_pair );
		this->_pairDst->Fill();

	}

	virtual void preEventLoop(){
		TreeAnalyzer::preEventLoop();
		book->cd();
	}

	virtual void makePairs( FemtoEvent* event, vector<FemtoTrackProxy> &col1, vector<FemtoTrackProxy> &col2 ){
		TLorentzVector lv1, lv2, lv;
		for ( FemtoTrackProxy& _proxy1 : col1 ){
			for ( FemtoTrackProxy& _proxy2 : col2 ){
				
				// make the MTD track always first
				if ( nullptr != _proxy1._mtdPid )
					addPair( event, _proxy1, _proxy2 );
				else if ( nullptr != _proxy2._mtdPid )
					addPair( event, _proxy2, _proxy1 );
			} // loop col1
		} // loop col2
	}


	virtual void analyzeEvent(){
		FemtoEvent * event = this->_rEvent.get();


		if ( cenBinMin >= 0 && cenBinMax >= 0 ){
			if ( event->mBin16 < cenBinMin )
				return;
			if ( event->mBin16 > cenBinMax )
				return;
		}

		// if ( fabs(event->mPrimaryVertex_mX3 - event->mWeight) < 200 )
		// 	return;

		size_t nTracks = this->_rTracks.N();

		size_t nTOF = 0;
		size_t nMTD = 0;

		pos_mtd.clear();
		neg_mtd.clear();
		pos_tof.clear();
		neg_tof.clear();
		
		FemtoTrackProxy p1;
		for ( size_t i = 0; i < nTracks; i++ ){
			p1.reset();
			p1.assemble( i, this->_rTracks, this->_rHelices, this->_rBTofPid );
			p1.setMtdPidTraits( this->_rMtdPid );

			if ( p1._track->mPt < 0.01 ) continue;

			int charge = p1._track->charge();

			if ( p1._track->mBTofPidTraitsIndex >= 0 && _lowPtTrackFilter.pass( p1 ) ){
				
				double p = p1._track->mPt * cosh( p1._track->mEta );
				// double zb = _lowFilter.zb( p1, "mu" );
				// double zbMu = _lowFilter.zbCorr( p1, "mu" );
				// double zbPi = _lowFilter.zbCorr( p1, "pi" );

				// book->fill( "zb_p", p, zb );

				// if ( _lowFilter.pass( p1 ) ){
					// book->fill( "zb_p_signal", p,zb );
					if ( charge > 0 )
						pos_tof.push_back( p1 );
					else 
						neg_tof.push_back( p1 );
					nTOF++;
				// }

				// mutually exclusive with MTD muons, so if we get here then dont need to check MTD  muons
				continue;
			} // low p muon


			// if we get here it was not a LowP Muon Candidate
			if ( nullptr == p1._mtdPid  ) continue;

			int ipt1 = bins_pt_mva.findBin( p1._track->mPt );
			if ( ipt1 < 0 ){
				p1._pid = -999;
				continue;
			} else {
				MuonMVAFilter::fillVars( p1 );
				p1._pid = mlps[ipt1].evaluate( p1 );

				if ( charge > 0 )
					pos_mtd.push_back( p1 );
				else 
					neg_mtd.push_back( p1 );
			}
		}// loop track i

		// if ( (pos_mtd.size() + neg_tof.size()) > 1 || (neg_mtd.size() + pos_tof.size()) > 1  ){
		// 	LOG_F( INFO, ">>pos_mtd.size() = %lu", pos_mtd.size() );
		// 	LOG_F( INFO, "neg_mtd.size() = %lu", neg_mtd.size() );
			
		// 	LOG_F( INFO, "pos_tof.size() = %lu", pos_tof.size() );
		// 	LOG_F( INFO, "<<neg_tof.size() = %lu", neg_tof.size() );
		// }
		

		// make pairs from the lists
		makePairs( event, pos_mtd, neg_tof );
		makePairs( event, neg_mtd, pos_tof );
		
		makePairs( event, pos_mtd, pos_tof );
		makePairs( event, neg_mtd, neg_tof );
	}


	virtual void postMake(){
		TreeAnalyzer::postMake();

		this->_pairDst->Write();
	}

};


#endif

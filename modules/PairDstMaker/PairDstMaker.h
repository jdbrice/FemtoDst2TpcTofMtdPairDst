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

#include "Filters/MuonMLPFilter.h"

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

	BranchReader<FemtoEvent> _fer;
	TClonesArrayReader<FemtoTrack> _ftr;
	TClonesArrayReader<FemtoMtdPidTraits> _fmtdr;

	MuonMLPFilter _mlp;
	TTree * _pairDst = nullptr;
	BranchWriter<FemtoPair> _fpw;
	FemtoPair _pair; 

public:

	const int DEBUG = 1;
	PairDstMaker() {}
	~PairDstMaker() {}

	virtual void initialize(){
		TreeAnalyzer::initialize();

		this->_fer.setup( this->chain, "Event" );
		this->_ftr.setup( this->chain, "Tracks" );
		this->_fmtdr.setup( this->chain, "MtdPidTraits" );

		_mlp.load( config, nodePath + ".MuonMLPFilter" );

		book->cd();

		this->_pairDst = new TTree( "PairDst", "" );
		this->_fpw.createBranch( this->_pairDst, "Pairs" );

	}
protected:

	void addPair(FemtoEvent * event, FemtoTrackProxy &p1, FemtoTrackProxy &p2){
		
		if ( p1._track->mId == p2._track->mId) return;

		this->_pair.reset();
		this->_pair.mVertexZ              = event->mPrimaryVertex_mX3;
		this->_pair.mGRefMult             = event->mGRefMult;

		this->_pair.d1_mDeltaY            = p1._mtdPid->mDeltaY;
		this->_pair.d1_mDeltaZ            = p1._mtdPid->mDeltaZ;
		this->_pair.d1_mDeltaTimeOfFlight = p1._mtdPid->mDeltaTimeOfFlight;
		this->_pair.d1_mMatchFlag         = p1._mtdPid->mMatchFlag;
		this->_pair.d1_mTriggerFlag       = p1._mtdPid->mTriggerFlag;
		this->_pair.d1_mCell              = p1._mtdPid->cell();
		this->_pair.d1_mModule            = p1._mtdPid->module();
		this->_pair.d1_mBackleg           = p1._mtdPid->backleg();

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


		this->_pair.d2_mDeltaY            = p2._mtdPid->mDeltaY;
		this->_pair.d2_mDeltaZ            = p2._mtdPid->mDeltaZ;
		this->_pair.d2_mDeltaTimeOfFlight = p2._mtdPid->mDeltaTimeOfFlight;
		this->_pair.d2_mMatchFlag         = p2._mtdPid->mMatchFlag;
		this->_pair.d2_mTriggerFlag       = p2._mtdPid->mTriggerFlag;
		this->_pair.d2_mCell              = p2._mtdPid->cell();
		this->_pair.d2_mModule            = p2._mtdPid->module();
		this->_pair.d2_mBackleg           = p2._mtdPid->backleg();

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

		// if ( lv.M() < 0.22 ){
		// 	LOG_F(INFO, "track1.pt, eta, phi = %f, %f, %f", p1._track->mPt, p1._track->mEta, p1._track->mPhi );
		// 	LOG_F(INFO, "track2.pt, eta, phi = %f, %f, %f", p2._track->mPt, p2._track->mEta, p2._track->mPhi );
		// }
		
		this->_pair.mChargeSum = p1._track->charge() + p2._track->charge();
		this->_pair.mMass      = lv.M();
		this->_pair.mPt        = lv.Pt();
		this->_pair.mEta       = lv.PseudoRapidity();
		this->_pair.mPhi       = lv.Phi();
		this->_pair.mRapidity  = lv.Rapidity();

		this->_fpw.set( this->_pair );
		this->_pairDst->Fill();

	}



	virtual void analyzeEvent(){
		FemtoEvent * event = this->_fer.get();

		

		size_t nTracks = this->_ftr.N();
		
		for ( size_t i = 0; i < nTracks; i++ ){
			FemtoTrackProxy p1;
			p1._track = this->_ftr.get( i );
			p1.setMtdPidTraits( this->_fmtdr );

			if ( nullptr == p1._mtdPid  ) continue;
			if ( p1._track->mPt < 0.01 ) continue;

			p1._pid = this->_mlp.evaluate( p1 );


			for ( size_t j = 0; j < nTracks; j++ ){
				if ( i == j ) continue;

				FemtoTrackProxy p2;
				p2._track = this->_ftr.get( j );
				p2.setMtdPidTraits( this->_fmtdr );

				if ( nullptr == p2._mtdPid  ) continue;
				if ( p2._track->mPt < 0.01 ) continue;

				p2._pid = this->_mlp.evaluate( p2 );

				addPair( event, p1, p2 );
			} // loop track j
		}// loop track i
	}


	virtual void postMake(){
		TreeAnalyzer::postMake();

		this->_pairDst->Write();
	}

};


#endif
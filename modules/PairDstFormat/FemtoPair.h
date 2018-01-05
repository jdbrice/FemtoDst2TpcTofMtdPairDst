#ifndef FEMTO_PAIR_H
#define FEMTO_PAIR_H

#include "FemtoDstFormat/FemtoDstBranch.h"

class FemtoPair : public FemtoDstBranch
{
public:
	virtual const char* classname() const { return "FemtoPair"; } 
	FemtoPair() {
		
	}
	~FemtoPair() {}

	void reset(){
		mVertexZ              = -999;
		mDeltaVertexZ         = -999;
		mGRefMult             = 0;
		mChargeSum            = -3;
		mPt                   = 0.0;
		mEta                  = 0.0;
		mPhi                  = 0.0;
		mMass                 = 0.0;
		mRapidity             = 0.0;
		d1_mDeltaY            = 0.0;
		d1_mDeltaZ            = 0.0;
		d1_mDeltaTimeOfFlight = 0.0;
		d1_mMatchFlag         = 0.0;
		d1_mTriggerFlag       = 0.0;
		d1_mCell              = 0.0;
		d1_mModule            = 0.0;
		d1_mBackleg           = 0.0;
		d1_mPt                = 0.0;
		d1_mEta               = 0.0;
		d1_mPhi               = 0.0;
		d1_mId                = 0.0;
		d1_mNHitsFit          = 0.0;
		d1_mNHitsMax          = 0.0;
		d1_mNHitsDedx         = 0.0;
		d1_mNSigmaPion        = 0.0;
		d1_mDCA               = 0.0;
		d1_mPid               = 0.0;
		d2_mDeltaY            = 0.0;
		d2_mDeltaZ            = 0.0;
		d2_mDeltaTimeOfFlight = 0.0;
		d2_mMatchFlag         = 0.0;
		d2_mTriggerFlag       = 0.0;
		d2_mCell              = 0.0;
		d2_mModule            = 0.0;
		d2_mBackleg           = 0.0;
		d2_mPt                = 0.0;
		d2_mEta               = 0.0;
		d2_mPhi               = 0.0;
		d2_mId                = 0.0;
		d2_mNHitsFit          = 0.0;
		d2_mNHitsMax          = 0.0;
		d2_mNHitsDedx         = 0.0;
		d2_mNSigmaPion        = 0.0;
		d2_mDCA               = 0.0;
		d2_mPid               = 0.0;
	}
	void copy( FemtoPair * that ){
		this->mVertexZ              = that->mVertexZ;
		this->mDeltaVertexZ         = that->mDeltaVertexZ;
		this->mGRefMult             = that->mGRefMult;
		this->mChargeSum            = that->mChargeSum;
		this->mPt                   = that->mPt;
		this->mEta                  = that->mEta;
		this->mPhi                  = that->mPhi;
		this->mMass                 = that->mMass;
		this->mRapidity             = that->mRapidity;
		this->d1_mDeltaY            = that->d1_mDeltaY;
		this->d1_mDeltaZ            = that->d1_mDeltaZ;
		this->d1_mDeltaTimeOfFlight = that->d1_mDeltaTimeOfFlight;
		this->d1_mMatchFlag         = that->d1_mMatchFlag;
		this->d1_mTriggerFlag       = that->d1_mTriggerFlag;
		this->d1_mCell              = that->d1_mCell;
		this->d1_mModule            = that->d1_mModule;
		this->d1_mBackleg           = that->d1_mBackleg;
		this->d1_mPt                = that->d1_mPt;
		this->d1_mEta               = that->d1_mEta;
		this->d1_mPhi               = that->d1_mPhi;
		this->d1_mId                = that->d1_mId;
		this->d1_mNHitsFit          = that->d1_mNHitsFit;
		this->d1_mNHitsMax          = that->d1_mNHitsMax;
		this->d1_mNHitsDedx         = that->d1_mNHitsDedx;
		this->d1_mNSigmaPion        = that->d1_mNSigmaPion;
		this->d1_mDCA               = that->d1_mDCA;
		this->d1_mPid               = that->d1_mPid;
		this->d2_mDeltaY            = that->d2_mDeltaY;
		this->d2_mDeltaZ            = that->d2_mDeltaZ;
		this->d2_mDeltaTimeOfFlight = that->d2_mDeltaTimeOfFlight;
		this->d2_mMatchFlag         = that->d2_mMatchFlag;
		this->d2_mTriggerFlag       = that->d2_mTriggerFlag;
		this->d2_mCell              = that->d2_mCell;
		this->d2_mModule            = that->d2_mModule;
		this->d2_mBackleg           = that->d2_mBackleg;
		this->d2_mPt                = that->d2_mPt;
		this->d2_mEta               = that->d2_mEta;
		this->d2_mPhi               = that->d2_mPhi;
		this->d2_mId                = that->d2_mId;
		this->d2_mNHitsFit          = that->d2_mNHitsFit;
		this->d2_mNHitsMax          = that->d2_mNHitsMax;
		this->d2_mNHitsDedx         = that->d2_mNHitsDedx;
		this->d2_mNSigmaPion        = that->d2_mNSigmaPion;
		this->d2_mDCA               = that->d2_mDCA;
		this->d2_mPid               = that->d2_mPid;
	}

	// Event
	Float_t 	mVertexZ;				// Event VertexZ
	Float_t 	mDeltaVertexZ;			// delta vz
	UShort_t 	mGRefMult;				// global RefMult
	

	// Pair
	Char_t    mChargeSum;
	Float_t   mPt;
	Float_t   mEta;
	Float_t   mPhi;
	Float_t   mMass;
	Float_t   mRapidity;

	// MTD Stuff
	Float_t   d1_mDeltaY;					// DeltaY between matched track-hit pair
	Float_t   d1_mDeltaZ;					// DeltaZ between matched track-hit pair
	Float_t   d1_mDeltaTimeOfFlight;		// Difference between measured and expected time-of-flight
	Char_t    d1_mMatchFlag;				// Matching flag indicating multiple matches
	Char_t    d1_mTriggerFlag;
	Char_t    d1_mCell;
	Char_t    d1_mModule;
	Char_t    d1_mBackleg;

	// Track Stuff
	Float_t 	d1_mPt;					// primary track px
	Float_t 	d1_mEta;				// primary track py
	Float_t 	d1_mPhi;				// primary track pz
	UShort_t 	d1_mId;					// track Id
	Char_t 		d1_mNHitsFit;			// q*nHitsFit - TPC
	Char_t 		d1_mNHitsMax;			// nHitsMax - TPC
	UChar_t 	d1_mNHitsDedx;			// nHitsDedx - TPC
	Float_t 	d1_mNSigmaPion;			// nsigmaPi * 100
	Float_t 	d1_mDCA;
	Float_t 	d1_mPid;

	// MTD Stuff
	Float_t   d2_mDeltaY;					// DeltaY between matched track-hit pair
	Float_t   d2_mDeltaZ;					// DeltaZ between matched track-hit pair
	Float_t   d2_mDeltaTimeOfFlight;		// Difference between measured and expected time-of-flight
	Char_t    d2_mMatchFlag;				// Matching flag indicating multiple matches
	Char_t    d2_mTriggerFlag;
	Char_t    d2_mCell;
	Char_t    d2_mModule;
	Char_t    d2_mBackleg;

	// Track Stuff
	Float_t 	d2_mPt;					// primary track px
	Float_t 	d2_mEta;				// primary track py
	Float_t 	d2_mPhi;				// primary track pz
	UShort_t 	d2_mId;					// track Id
	Char_t 		d2_mNHitsFit;			// q*nHitsFit - TPC
	Char_t 		d2_mNHitsMax;			// nHitsMax - TPC
	UChar_t 	d2_mNHitsDedx;			// nHitsDedx - TPC
	Float_t 	d2_mNSigmaPion;			// nsigmaPi * 100
	Float_t 	d2_mDCA;
	Float_t 	d2_mPid;

	ClassDef( FemtoPair, 2 )
};


#endif

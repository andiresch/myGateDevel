/*----------------------
  Copyright (C): OpenGATE Collaboration

  This software is distributed under the terms
  of the GNU Lesser General  Public Licence (LGPL)
  See GATE/LICENSE.txt for further details
  ----------------------*/

/*
  \brief Class GateDoseActor :
  \brief
*/

// gate
#include "GateDoseActor.hh"
#include "GateMiscFunctions.hh"

// g4
#include <G4EmCalculator.hh>
#include <G4VoxelLimits.hh>
#include <G4NistManager.hh>
#include <G4PhysicalConstants.hh>
#include <G4Gamma.hh>
#include <G4Proton.hh>
#include <G4Positron.hh>
#include <G4Deuteron.hh>
#include <G4Electron.hh>

//-----------------------------------------------------------------------------
GateDoseActor::GateDoseActor(G4String name, G4int depth):
  GateVImageActor(name,depth) {
  GateDebugMessageInc("Actor",4,"GateDoseActor() -- begin\n");

  mCurrentEvent=-1;
  mIsEdepImageEnabled = false;
  mIsLastHitEventImageEnabled = false;
  mIsEdepSquaredImageEnabled = false;
  mIsEdepUncertaintyImageEnabled = false;
  mIsDoseImageEnabled = true;
  mIsDoseSquaredImageEnabled = false;
  mIsDoseUncertaintyImageEnabled = false;
  mIsDoseToWaterImageEnabled = false;
  mIsDoseToWaterSquaredImageEnabled = false;
  mIsDoseToWaterUncertaintyImageEnabled = false;
  mIsNumberOfHitsImageEnabled = false;
  mIsDoseNormalisationEnabled = false;
  mIsDoseToWaterNormalisationEnabled = false;
  mDoseAlgorithmType = "VolumeWeighting";
  mImportMassImage = "";
  mExportMassImage = "";

  mVolumeFilter = "";
  mMaterialFilter = "";

  mDose2WaterWarningFlag = true;

  pMessenger = new GateDoseActorMessenger(this);
  GateDebugMessageDec("Actor",4,"GateDoseActor() -- end\n");
  emcalc = new G4EmCalculator;
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
/// Destructor
GateDoseActor::~GateDoseActor()  {
  delete pMessenger;
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void GateDoseActor::EnableDoseNormalisationToMax(bool b) {
  mIsDoseNormalisationEnabled = b;
  mDoseImage.SetNormalizeToMax(b);
  mDoseImage.SetScaleFactor(1.0);
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void GateDoseActor::EnableDoseNormalisationToIntegral(bool b) {
  mIsDoseNormalisationEnabled = b;
  mDoseImage.SetNormalizeToIntegral(b);
  mDoseImage.SetScaleFactor(1.0);
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
/// Construct
void GateDoseActor::Construct() {
  GateDebugMessageInc("Actor", 4, "GateDoseActor -- Construct - begin\n");
  GateVImageActor::Construct();

  // Find G4_WATER. This it needed here because we will used this
  // material for dedx computation for DoseToWater.
  G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");

  // Record the stepHitType
  mUserStepHitType = mStepHitType;

  // Enable callbacks
  EnableBeginOfRunAction(true);
  EnableBeginOfEventAction(true);
  EnablePreUserTrackingAction(true);
  EnableUserSteppingAction(true);

  // Check if at least one image is enabled
  if (!mIsEdepImageEnabled &&
      !mIsDoseImageEnabled &&
      !mIsDoseToWaterImageEnabled &&
      !mIsNumberOfHitsImageEnabled &&
      mExportMassImage=="")  {
    GateError("The DoseActor " << GetObjectName()
              << " does not have any image enabled ...\n Please select at least one ('enableEdep true' for example)");
  }

  // Output Filename
  mEdepFilename = G4String(removeExtension(mSaveFilename))+"-Edep."+G4String(getExtension(mSaveFilename));
  mDoseFilename = G4String(removeExtension(mSaveFilename))+"-Dose."+G4String(getExtension(mSaveFilename));
  mDoseToWaterFilename = G4String(removeExtension(mSaveFilename))+"-DoseToWater."+G4String(getExtension(mSaveFilename));
  mNbOfHitsFilename = G4String(removeExtension(mSaveFilename))+"-NbOfHits."+G4String(getExtension(mSaveFilename));

  // Set origin, transform, flag
  SetOriginTransformAndFlagToImage(mEdepImage);
  SetOriginTransformAndFlagToImage(mDoseImage);
  SetOriginTransformAndFlagToImage(mNumberOfHitsImage);
  SetOriginTransformAndFlagToImage(mLastHitEventImage);
  SetOriginTransformAndFlagToImage(mDoseToWaterImage);
  SetOriginTransformAndFlagToImage(mMassImage);

  // Resize and allocate images
  if (mIsEdepSquaredImageEnabled || mIsEdepUncertaintyImageEnabled ||
      mIsDoseSquaredImageEnabled || mIsDoseUncertaintyImageEnabled ||
      mIsDoseToWaterSquaredImageEnabled || mIsDoseToWaterUncertaintyImageEnabled) {
    mLastHitEventImage.SetResolutionAndHalfSize(mResolution, mHalfSize, mPosition);
    mLastHitEventImage.Allocate();
    mIsLastHitEventImageEnabled = true;
  }
  if (mIsEdepImageEnabled) {
    //  mEdepImage.SetLastHitEventImage(&mLastHitEventImage);
    mEdepImage.EnableSquaredImage(mIsEdepSquaredImageEnabled);
    mEdepImage.EnableUncertaintyImage(mIsEdepUncertaintyImageEnabled);
    // Force the computation of squared image if uncertainty is enabled
    if (mIsEdepUncertaintyImageEnabled) mEdepImage.EnableSquaredImage(true);
    mEdepImage.SetResolutionAndHalfSize(mResolution, mHalfSize, mPosition);
    mEdepImage.Allocate();
    mEdepImage.SetFilename(mEdepFilename);
  }
  if (mIsDoseImageEnabled) {
    // mDoseImage.SetLastHitEventImage(&mLastHitEventImage);
    mDoseImage.EnableSquaredImage(mIsDoseSquaredImageEnabled);
    mDoseImage.EnableUncertaintyImage(mIsDoseUncertaintyImageEnabled);
    mDoseImage.SetResolutionAndHalfSize(mResolution, mHalfSize, mPosition);
    // Force the computation of squared image if uncertainty is enabled
    if (mIsDoseUncertaintyImageEnabled) mDoseImage.EnableSquaredImage(true);

    // DD(mDoseImage.GetVoxelVolume());
    //mDoseImage.SetScaleFactor(1e12/mDoseImage.GetVoxelVolume());
    mDoseImage.Allocate();
    mDoseImage.SetFilename(mDoseFilename);
  }
  if (mIsDoseToWaterImageEnabled) {
    mDoseToWaterImage.EnableSquaredImage(mIsDoseToWaterSquaredImageEnabled);
    mDoseToWaterImage.EnableUncertaintyImage(mIsDoseToWaterUncertaintyImageEnabled);
    // Force the computation of squared image if uncertainty is enabled
    if (mIsDoseToWaterUncertaintyImageEnabled) mDoseToWaterImage.EnableSquaredImage(true);
    mDoseToWaterImage.SetResolutionAndHalfSize(mResolution, mHalfSize, mPosition);
    mDoseToWaterImage.Allocate();
    mDoseToWaterImage.SetFilename(mDoseToWaterFilename);
  }
  if (mIsNumberOfHitsImageEnabled) {
    mNumberOfHitsImage.SetResolutionAndHalfSize(mResolution, mHalfSize, mPosition);
    mNumberOfHitsImage.Allocate();
  }

  if (mIsDoseImageEnabled &&
      (mExportMassImage != "" || mDoseAlgorithmType == "MassWeighting" ||
       mVolumeFilter != ""    || mMaterialFilter != "")) {
    mVoxelizedMass.SetMaterialFilter(mMaterialFilter);
    mVoxelizedMass.SetVolumeFilter(mVolumeFilter);
    mVoxelizedMass.SetExternalMassImage(mImportMassImage);
    mVoxelizedMass.Initialize(mVolumeName, &mDoseImage.GetValueImage());

    if (mExportMassImage != "") {
      mMassImage.SetResolutionAndHalfSize(mResolution, mHalfSize, mPosition);
      mMassImage.Allocate();

      mVoxelizedMass.UpdateImage(&mMassImage);

      mMassImage.Write(mExportMassImage);
    }
  }

  if (mExportMassImage!="" && mImportMassImage!="")
    GateWarning("Exported mass image will be the same as the imported one.");

  if (mDoseAlgorithmType != "MassWeighting") {
    mDoseAlgorithmType = "VolumeWeighting";

    if (mImportMassImage != "")
      GateWarning("importMassImage command is only compatible with MassWeighting algorithm. Ignored. ");
  }

  // Print information
  GateMessage("Actor", 1,
              "Dose DoseActor    = '" << GetObjectName() << "'\n" <<
              "\tDose image        = " << mIsDoseImageEnabled << Gateendl <<
              "\tDose squared      = " << mIsDoseSquaredImageEnabled << Gateendl <<
              "\tDose uncertainty  = " << mIsDoseUncertaintyImageEnabled << Gateendl <<
              "\tDose to water image        = " << mIsDoseToWaterImageEnabled << Gateendl <<
              "\tDose to water squared      = " << mIsDoseToWaterSquaredImageEnabled << Gateendl <<
              "\tDose to water uncertainty  = " << mIsDoseToWaterUncertaintyImageEnabled << Gateendl <<
              "\tEdep image        = " << mIsEdepImageEnabled << Gateendl <<
              "\tEdep squared      = " << mIsEdepSquaredImageEnabled << Gateendl <<
              "\tEdep uncertainty  = " << mIsEdepUncertaintyImageEnabled << Gateendl <<
              "\tNumber of hit     = " << mIsNumberOfHitsImageEnabled << Gateendl <<
              "\t     (last hit)   = " << mIsLastHitEventImageEnabled << Gateendl <<
              "\tDose algorithm    = " << mDoseAlgorithmType << Gateendl <<
              "\tMass image (import) = " << mImportMassImage << Gateendl <<
              "\tMass image (export) = " << mExportMassImage << Gateendl <<
              "\tEdepFilename      = " << mEdepFilename << Gateendl <<
              "\tDoseFilename      = " << mDoseFilename << Gateendl <<
              "\tNb Hits filename  = " << mNbOfHitsFilename << Gateendl);

  ResetData();
  GateMessageDec("Actor", 4, "GateDoseActor -- Construct - end\n");
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
/// Save data
void GateDoseActor::SaveData() {
  GateVActor::SaveData(); // (not needed because done into GateImageWithStatistic)

  if (mIsEdepImageEnabled) mEdepImage.SaveData(mCurrentEvent+1);
  if (mIsDoseImageEnabled) {
    if (mIsDoseNormalisationEnabled)
      mDoseImage.SaveData(mCurrentEvent+1, true);
    else
      mDoseImage.SaveData(mCurrentEvent+1, false);
  }

  if (mIsDoseToWaterImageEnabled) {
    if (mIsDoseToWaterNormalisationEnabled)
      mDoseToWaterImage.SaveData(mCurrentEvent+1, true);
    else
      mDoseToWaterImage.SaveData(mCurrentEvent+1, false);
  }

  if (mIsLastHitEventImageEnabled) {
    mLastHitEventImage.Fill(-1); // reset
  }

  if (mIsNumberOfHitsImageEnabled) {
    mNumberOfHitsImage.Write(mNbOfHitsFilename);
  }
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void GateDoseActor::ResetData() {
  if (mIsLastHitEventImageEnabled) mLastHitEventImage.Fill(-1);
  if (mIsEdepImageEnabled) mEdepImage.Reset();
  if (mIsDoseImageEnabled) mDoseImage.Reset();
  if (mIsDoseToWaterImageEnabled) mDoseToWaterImage.Reset();
  if (mIsNumberOfHitsImageEnabled) mNumberOfHitsImage.Fill(0);
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void GateDoseActor::BeginOfRunAction(const G4Run * r) {
  GateVActor::BeginOfRunAction(r);
  GateDebugMessage("Actor", 3, "GateDoseActor -- Begin of Run\n");
  mDose2WaterWarningFlag = true;
  // ResetData(); // Do no reset here !! (when multiple run);
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// Callback at each event
void GateDoseActor::BeginOfEventAction(const G4Event * e) {
  GateVActor::BeginOfEventAction(e);
  mCurrentEvent++;
  GateDebugMessage("Actor", 3, "GateDoseActor -- Begin of Event: "<< mCurrentEvent << Gateendl);
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void GateDoseActor::UserPreTrackActionInVoxel(const int /*index*/, const G4Track* track)
{
  if(track->GetDefinition() == G4Gamma::Gamma()) { mStepHitType = PostStepHitType; }
  else { mStepHitType = mUserStepHitType; }
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void GateDoseActor::UserSteppingActionInVoxel(const int index, const G4Step* step) {
  GateDebugMessageInc("Actor", 4, "GateDoseActor -- UserSteppingActionInVoxel - begin\n");
  GateDebugMessageInc("Actor", 4, "enedepo = " << step->GetTotalEnergyDeposit() << Gateendl);
  GateDebugMessageInc("Actor", 4, "weight = " <<  step->GetTrack()->GetWeight() << Gateendl);
  const double weight = step->GetTrack()->GetWeight();
  const double edep = step->GetTotalEnergyDeposit()*weight;//*step->GetTrack()->GetWeight();

  // if no energy is deposited or energy is deposited outside image => do nothing
  if (edep == 0) {
    GateDebugMessage("Actor", 5, "edep == 0 : do nothing\n");
    GateDebugMessageDec("Actor", 4, "GateDoseActor -- UserSteppingActionInVoxel -- end\n");
    return;
  }
  if (index < 0) {
    GateDebugMessage("Actor", 5, "index < 0 : do nothing\n");
    GateDebugMessageDec("Actor", 4, "GateDoseActor -- UserSteppingActionInVoxel -- end\n");
    return;
  }

  if (mVolumeFilter != "" && mVolumeFilter+"_phys" != step->GetPreStepPoint()->GetPhysicalVolume()->GetName())
    return;

  if (mMaterialFilter != "" && mMaterialFilter != step->GetPreStepPoint()->GetMaterial()->GetName())
    return;

  // compute sameEvent
  // sameEvent is false the first time some energy is deposited for each primary particle
  bool sameEvent=true;
  if (mIsLastHitEventImageEnabled) {
    GateDebugMessage("Actor", 2,  "GateDoseActor -- UserSteppingActionInVoxel: Last event in index = " << mLastHitEventImage.GetValue(index) << Gateendl);
    if (mCurrentEvent != mLastHitEventImage.GetValue(index)) {
      sameEvent = false;
      mLastHitEventImage.SetValue(index, mCurrentEvent);
    }
  }

  //---------------------------------------------------------------------------------
  // Volume weighting
  double density = step->GetPreStepPoint()->GetMaterial()->GetDensity();
  //---------------------------------------------------------------------------------

  //---------------------------------------------------------------------------------
  // Mass weighting OR filter
  if (mDoseAlgorithmType == "MassWeighting" || mMaterialFilter != "" || mVolumeFilter != "")
    density = mVoxelizedMass.GetDoselMass(index)/mDoseImage.GetVoxelVolume();
  //---------------------------------------------------------------------------------

  if (mMaterialFilter != "") {
    GateDebugMessage("Actor", 3,  "GateDoseActor -- UserSteppingActionInVoxel: material filter debug = " << Gateendl
                     << " material name        = " << step->GetPreStepPoint()->GetMaterial()->GetName() << Gateendl
                     << " density              = " << G4BestUnit(mVoxelizedMass.GetPartialMassWithMatName(index)/mVoxelizedMass.GetPartialVolumeWithMatName(index), "Volumic Mass") << Gateendl
                     << " dosel cubic volume   = " << G4BestUnit(mDoseImage.GetVoxelVolume(), "Volume") << Gateendl
                     << " partial cubic volume = " << G4BestUnit(mVoxelizedMass.GetPartialVolumeWithMatName(index), "Volume") << Gateendl);
  }

  if (mVolumeFilter != "") {
    GateDebugMessage("Actor", 3,  "GateDoseActor -- UserSteppingActionInVoxel: volume filter debug = " << Gateendl
                     << " volume name          = " << step->GetPreStepPoint()->GetPhysicalVolume()->GetName() << Gateendl
                     << " Dose scored inside volume filtered volume !" << Gateendl);
  }

  double dose=0.;
  if (mIsDoseImageEnabled) {
    // ------------------------------------
    // Convert deposited energy into Gray
    dose = edep/density/mDoseImage.GetVoxelVolume()/gray;
    // ------------------------------------
  //G4double edepA = step->GetTotalEnergyDeposit()*weight;//*step->GetTrack()->GetWeight();

    // ------------------------------------
    //G4double doseA;
    //G4double densityA = step->GetPreStepPoint()->GetMaterial()->GetDensity();
    //double eee = edepA/(eV);
    //double ddd = densityA/(kg/m3);
    //double vvv = mDoseImage.GetVoxelVolume()/m3;
    //doseA = edepA/(eV)/densityA/(kg/m3)/mDoseImage.GetVoxelVolume()/m3;
    //doseA = eee/ddd/vvv;
	  //G4cout<<"density orig: "<< densityA<<G4endl;
	  //G4cout<<"density: "<< densityA/(kg/m3)<<" kg/m3"<<G4endl;
	  //G4cout<<"volume orig: "<< mDoseImage.GetVoxelVolume()<<" "<<G4endl;
	  //G4cout<<"volume: "<< mDoseImage.GetVoxelVolume()/cm3<<" cm3"<<G4endl;
	  //G4cout<<"volume: "<< mDoseImage.GetVoxelVolume()/m3<<" m3"<<G4endl;
	  //G4cout<<"Edep orig: "<< edep<<G4endl;
	  //G4cout<<"Edep: "<<edep/eV<<" eV"<<G4endl;
	  //G4cout<<"Edep: "<<edep/keV<<" keV"<<G4endl;
	  //G4cout<<"Edep: "<<edep/MeV<<" MeV"<<G4endl;
	  //G4cout<<"dose orig: "<< dose<<G4endl;
	  //G4cout<<"doseA orig: "<< doseA<<G4endl;
	  //G4cout<<"doseA orig*qe: "<< doseA*1.60218e-19<<G4endl;
	  //G4cout<<"doseA : "<< doseA/(eV/kg)<<"eV/kg"<<G4endl;
	  //G4cout<<"doseA : "<< doseA/(MeV/kg)<<"MeV/kg"<<G4endl;
	  //G4cout<<"doseA : "<< doseA/(eV/kg)*1.60218e-19<<"J/kg"<<G4endl;
	  //G4cout<<"gray f : "<< 1/gray<<" unit"<<G4endl<<G4endl ;
    GateDebugMessage("Actor", 2,  "GateDoseActor -- UserSteppingActionInVoxel:\tdose = "
                     << G4BestUnit(dose, "Dose")
                     << " rho = "
                     << G4BestUnit(density, "Volumic Mass")<< Gateendl );
  }

  double doseToWater = 0;
  if (mIsDoseToWaterImageEnabled)
    {
      // to get nuclear inelastic cross-section, see "geant4.9.4.p01/examples/extended/hadronic/Hadr00/"
      // #include "G4HadronicProcessStore.hh"
      // G4HadronicProcessStore* store = G4HadronicProcessStore::Instance();
      // store->GetInelasticCrossSectionPerAtom(particle,e,elm);

      double cut = DBL_MAX;
      cut=1;

      G4String material = step->GetPreStepPoint()->GetMaterial()->GetName();
      double Energy = step->GetPreStepPoint()->GetKineticEnergy();
      G4String PartName = step->GetTrack()->GetDefinition()->GetParticleName();
      //double DEDX=0, DEDX_Water=0;

      // Dose to water: it could be possible to make this process more
      // generic by choosing any material in place of water
      double Volume = mDoseToWaterImage.GetVoxelVolume();

  //========================= AR 3.2.2017 ============================

	 const G4Material* materialTest = step->GetPreStepPoint()->GetMaterial();//->GetName();
     double energy1 = step->GetPreStepPoint()->GetKineticEnergy();
     double energy2 = step->GetPostStepPoint()->GetKineticEnergy();
     double energy=(energy1+energy2)/2;
     const G4ParticleDefinition* partname = step->GetTrack()->GetDefinition();//->GetParticleName();

	//G4double dedxAR = emcalc->ComputeElectronicDEDX(energyAR, partnameAR, materialAR);
	G4double DEDX = 1;
	G4double DEDX_Water = 1;
	// if not a photon get DEDX
	if ( partname->GetPDGEncoding() != 22  )
	{       
		DEDX = emcalc->ComputeTotalDEDX(energy,partname,materialTest, cut);
        DEDX_Water = emcalc->ComputeTotalDEDX(energy, PartName, "G4_WATER", cut);
		//G4cout<< "PDG code 22: " << PartName<<G4endl;
	}
	doseToWater=edep/density/Volume/gray*(DEDX_Water/1.)/(DEDX/(density*e_SI));
	
	if (DEDX_Water == 0 || DEDX == 0) {
		GateWarning("DEDX = 0 in doseToWater, Edep ommited");
		G4cout<<"PartName: "<< PartName<<" Edep: "<<edep/gray<<G4endl; 
		doseToWater = 0.0;
	}
	
    //if ( partnameAR->GetPDGEncoding() != 2212 && PartName!= "e-" &&PartName!= "e+" )
    //{
		//G4cout<< "name: "<<PartName<< " pdg code:" <<  partnameAR->GetPDGEncoding() <<G4endl;
	//}
  //========================= AR 3.2.2017 ============================

      // Other particles should be taken into account (Helium etc), but bug ? FIXME
      //if (PartName== "proton" || PartName== "e-" || PartName== "e+" || PartName== "deuteron") {
        ////if (PartName != "O16[0.0]" && PartName != "alpha" && PartName != "Be7[0.0]" && PartName != "C12[0.0]" ){

        //DEDX = emcalc->ComputeTotalDEDX(Energy, PartName, material, cut);
        //DEDX_Water = emcalc->ComputeTotalDEDX(Energy, PartName, "G4_WATER", cut);

        //doseToWater=edep/density/Volume/gray*(DEDX_Water/1.)/(DEDX/(density*e_SI));
      //}
      //else {
        //DEDX = emcalc->ComputeTotalDEDX(100, "proton", material, cut);
        //DEDX_Water = emcalc->ComputeTotalDEDX(100, "proton", "G4_WATER", cut);
        //doseToWater=edep/density/Volume/gray*(DEDX_Water/1.)/(DEDX/(density*e_SI));
      //}

      //GateDebugMessage("Actor", 2,  "GateDoseActor -- UserSteppingActionInVoxel:\tdose to water = "
                       //<< G4BestUnit(doseToWater, "Dose to water")
                       //<< " rho = "
                       //<< G4BestUnit(density, "Volumic Mass")<< Gateendl );
    //}
//=========================================================================================================
// ===== GATE develop version
//=========================================================================================================

      //G4Material * material = step->GetPreStepPoint()->GetMaterial();
      //static G4Material * water = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
      //double energy = step->GetPreStepPoint()->GetKineticEnergy();
      //double DEDX=0, DEDX_Water=0;

      //// Dose to water: it could be possible to make this process more
      //// generic by choosing any material in place of water
      //double volume = mDoseToWaterImage.GetVoxelVolume();

      //// Get current particle
      //const G4ParticleDefinition * p = step->GetTrack()->GetParticleDefinition();
      //if (p == G4Proton::Proton() or
          //p == G4Electron::Electron() or
          //p == G4Positron::Positron() or
          //p == G4Deuteron::Deuteron() or
          //p == G4Gamma::Gamma()) {
        //// For Gamma, we consider the DEDX from Electron
        //if (p == G4Gamma::Gamma()) p = G4Electron::Electron();
        //DEDX = emcalc->ComputeTotalDEDX(energy, p, material, cut);
        //DEDX_Water = emcalc->ComputeTotalDEDX(energy, p, water, cut);
        //doseToWater = edep/density/volume/gray*(DEDX_Water/1.0)/(DEDX/(density*e_SI));
        //if (DEDX_Water == 0 or DEDX == 0) doseToWater = 0.0; // to avoid inf or NaN
      //}
      //else {
        //if (mDose2WaterWarningFlag) {
          //GateMessage("Actor", 0, "WARNING: DoseToWater with a particle which is not proton/electron/positron/gamma/deuteron: results could be wrong." << G4endl);
          //mDose2WaterWarningFlag = false;
        //}
      //}

      //GateDebugMessage("Actor", 2,  "GateDoseActor -- UserSteppingActionInVoxel:\tdose to water = "
                       //<< G4BestUnit(doseToWater, "Dose to water")
                       //<< " rho = "
                       //<< G4BestUnit(density, "Volumic Mass")<< Gateendl );
    //}
//=========================================================================================================

  if (mIsEdepImageEnabled) {
    GateDebugMessage("Actor", 2, "GateDoseActor -- UserSteppingActionInVoxel:\tedep = " << G4BestUnit(edep, "Energy") << Gateendl);
  }

  if (mIsDoseImageEnabled)
    {
      if (mIsDoseUncertaintyImageEnabled || mIsDoseSquaredImageEnabled)
        {
          if (sameEvent) mDoseImage.AddTempValue(index, dose);
          else mDoseImage.AddValueAndUpdate(index, dose);
        }
      else mDoseImage.AddValue(index, dose);
    }

  if (mIsDoseToWaterImageEnabled)
    {
      if (mIsDoseToWaterUncertaintyImageEnabled || mIsDoseToWaterSquaredImageEnabled)
        {
          if (sameEvent) mDoseToWaterImage.AddTempValue(index, doseToWater);
          else mDoseToWaterImage.AddValueAndUpdate(index, doseToWater);
        }
      else mDoseToWaterImage.AddValue(index, doseToWater);
    }

  if (mIsEdepImageEnabled)
    {
      if (mIsEdepUncertaintyImageEnabled || mIsEdepSquaredImageEnabled)
        {
          if (sameEvent) mEdepImage.AddTempValue(index, edep);
          else mEdepImage.AddValueAndUpdate(index, edep);
        }
      else
        {
          mEdepImage.AddValue(index, edep);

        }
    }

  if (mIsNumberOfHitsImageEnabled) mNumberOfHitsImage.AddValue(index, weight);

  GateDebugMessageDec("Actor", 4, "GateDoseActor -- UserSteppingActionInVoxel -- end\n");
}
}
//-----------------------------------------------------------------------------

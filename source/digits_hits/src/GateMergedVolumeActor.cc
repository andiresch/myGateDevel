/*----------------------
  Copyright (C): OpenGATE Collaboration

  This software is distributed under the terms
  of the GNU Lesser General  Public Licence (LGPL)
  See GATE/LICENSE.txt for further details
  ----------------------*/

#include "GateMergedVolumeActor.hh"
#include "GateMergedVolumeActorMessenger.hh"
#include "GateActorMessenger.hh"
#include "GateDetectorConstruction.hh"
#include "G4TransportationManager.hh"
#include "G4SteppingManager.hh"
#include "G4EventManager.hh"

GateMergedVolumeActor::GateMergedVolumeActor(G4String name, G4int depth)
: GateVActor(name,depth)
{
  GateMessage("Actor",4,"GateMergedVolumeActor() -- begin\n");
  pActorMessenger = new GateActorMessenger(this);
  pMergedVolumeActorMessenger = new GateMergedVolumeActorMessenger(this);
  GateMessage("Actor",4,"GateMergedVolumeActor() -- end\n");
}

GateMergedVolumeActor::~GateMergedVolumeActor()
{
  GateMessage("Actor",4,"~GateMergedVolumeActor() -- begin\n");
  delete pMergedVolumeActorMessenger;
  delete pActorMessenger;
  GateMessage("Actor",4,"~GateMergedVolumeActor() -- end\n");
}

void GateMergedVolumeActor::Construct()
{
  GateVActor::Construct();

  // Enable callbacks
  EnableBeginOfRunAction( false );
  EnableBeginOfEventAction( false );
  EnablePreUserTrackingAction( false );
  EnableUserSteppingAction( true );

  // Store the pointer to solid to merge
  for( std::vector<G4String>::size_type i = 0; i < mVolToMerge.size(); ++i )
  {
    GateVVolume *vol = GateObjectStore::GetInstance()->FindVolumeCreator( mVolToMerge[i] );
    mSolidVolToMerge.push_back( vol->GetLogicalVolume()->GetSolid() );
    mPhysicalVolToMerge.push_back( vol->GetPhysicalVolume() );
  }

  ResetData();
}

void GateMergedVolumeActor::ListOfVolumesToMerge( G4String& vol )
{
  // Read volume to merge separated by a comma and store the solid
  // The priority of the volume is from the first declared volume to the last
  std::istringstream iss( vol );
  while( iss.good() )
  {
    G4String volumeName;
    getline( iss, volumeName, ',' );
    volumeName.erase(
      remove_if( volumeName.begin(), volumeName.end(), isspace ),
      volumeName.end()
    );

    // Check if the volume exists
    GateObjectStore::GetInstance()->FindVolumeCreator( volumeName );
    // Store the name
    mVolToMerge.push_back( volumeName );
  }
}

/*G4bool GateMergedVolumeActor::ProcessHits(G4Step * step , G4TouchableHistory* th)
{
  G4cout << "GateMergedVolumeActor::ProcessHits" << G4endl;
  return true;
}*/

G4bool GateMergedVolumeActor::ProcessHits(G4Step * step , G4TouchableHistory* )
{
  // Get the stack and stepping manager
  static G4EventManager* theEventMgr = G4EventManager::GetEventManager();
  G4SteppingManager* theStepMgr = theEventMgr->GetTrackingManager()->GetSteppingManager();
  G4StackManager* theStackMgr = theEventMgr->GetStackManager();

  // Get information about pre step (position, momentum, step length)
  G4ThreeVector const preStepPos = step->GetPreStepPoint()->GetPosition();
  G4ThreeVector const preStepDir = step->GetPreStepPoint()->GetMomentumDirection();
  G4double const stepLength = step->GetStepLength();

  // Get information about post step (position)
  G4ThreeVector const postStepPos = step->GetPostStepPoint()->GetPosition();
  G4ThreeVector const postStepDir = step->GetPostStepPoint()->GetMomentumDirection();

  // Loop over the volume(s) to merge and priority to the first volume
  for( std::vector<G4VSolid*>::size_type i = 0; i < mSolidVolToMerge.size(); ++i )
  {
    //Get the tolerance of the volume
    G4double const tolerance = mSolidVolToMerge[ i ]->GetTolerance();

    // Get the coordinates in the world space
    G4AffineTransform transform( mPhysicalVolToMerge[ i ]->GetRotation(), mPhysicalVolToMerge[ i ]->GetTranslation() );
    transform.Invert();
    // Add a small distance (tolerance) to position to avoid boundary problems
    //G4ThreeVector const preStepTolerance = preStepPos + tolerance * preStepDir;
    G4ThreeVector const positionVolRef = transform.TransformPoint( preStepPos );
    G4ThreeVector const dirVolRef = transform.TransformAxis( preStepDir );

    // Check the distance to in for each solid
    G4double const distanceToInVolumeToMerge = mSolidVolToMerge[ i ]->DistanceToIn( positionVolRef, dirVolRef );

    // Check if the particle is already inside the merged volume
    EInside const insideVol = mSolidVolToMerge[ i ]->Inside( positionVolRef );

    /*G4cout << "-----------------" << G4endl;
    G4cout << "*** IN VOLUME ***: " << i << G4endl;
    mSolidVolToMerge[ i ]->DumpInfo();
    G4cout << "Pre step position  : " << G4BestUnit( preStepPos, "Length" ) << G4endl;
    G4cout << "Post step position : " << G4BestUnit( postStepPos, "Length" ) << G4endl;
    G4cout << "Point Pre after transformation: " << G4BestUnit( positionVolRef, "Length" ) << G4endl;
    G4cout << "Dir Pre after transformation: " << G4BestUnit( dirVolRef, "Length" ) << G4endl;
    G4cout << "Step length: " << G4BestUnit( stepLength, "Length" ) << G4endl;
    G4cout << "Distance to in volume to merge: " << G4BestUnit( distanceToInVolumeToMerge, "Length" ) << G4endl;*/

    if( distanceToInVolumeToMerge <= stepLength && insideVol == kOutside )
    {
      // Compute the coordinates on the merged volume
      G4ThreeVector const intersectionPointOnVolume = preStepPos + ( distanceToInVolumeToMerge + tolerance ) * preStepDir;

      // Get some important track information before to kill it
      G4double const trackEnergy = theStepMgr->GetStep()->GetPreStepPoint()->GetTotalEnergy(); // Take energy of the pre step
      G4double const trackTime = theStepMgr->GetTrack()->GetGlobalTime();
      G4ParticleDefinition const* particleDefinition = theStepMgr->GetTrack()->GetParticleDefinition();
      G4int const newTrackID = theStepMgr->GetTrack()->GetTrackID();
      G4int const newParentID = theStepMgr->GetTrack()->GetParentID();

      // Kill the track and set the position of the new track on the merged volume
      theStepMgr->GetTrack()->SetWeight( 0 );
      theStepMgr->GetTrack()->SetKineticEnergy( 0.0 );
      //theStepMgr->GetTrack()->SetTrackStatus( fStopAndKill );
      theStepMgr->GetTrack()->SetTrackStatus( fKillTrackAndSecondaries );

      // Create a new particle
      G4DynamicParticle* dynParticle = new G4DynamicParticle( particleDefinition, preStepDir, trackEnergy );

      // Create a new track
      G4Track* newTrack = new G4Track( dynParticle, trackTime, intersectionPointOnVolume );
      newTrack->SetTrackID( newTrackID );
      newTrack->SetParentID( newParentID );

      // Push the new track to the stack
      theStackMgr->PushOneTrack( newTrack );

      /*G4cout << "***************************" << G4endl;
      G4cout << "The particle is entering in the merged volume number: " << i << G4endl;
      G4cout << "Coordinates of intersection: " << G4BestUnit( intersectionPointOnVolume, "Length" ) << G4endl;*/
    }

    //G4cout << "Post Step is inside?: " << mSolidVolToMerge[ i ]->Inside( postStepPos ) << G4endl;
    //G4cout << "Distance to out pre step: " << G4BestUnit( distanceToOutVolumeToMerge, "Length" ) << G4endl;
    // Check if the track is already in the merged volume
    /*if( mSolidVolToMerge[ i ]->Inside( preStepPos ) == kInside )
    {
      G4cout << "Distance to out volume to merge: " << G4BestUnit( distanceToOutVolumeToMerge, "Length" ) << G4endl;
      // Check if post step in on the surface
      if( mSolidVolToMerge[ i ]->Inside( postStepPos ) == kSurface )
      {
        G4cout << "Particle on boundary!!!" << G4endl;
      }
    }*/
  }
  return true;
}

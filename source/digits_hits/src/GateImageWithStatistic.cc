/*----------------------
  Copyright (C): OpenGATE Collaboration

  This software is distributed under the terms
  of the GNU Lesser General  Public Licence (LGPL)
  See GATE/LICENSE.txt for further details
  ----------------------*/



#ifndef GATEIMAGEWITHSTATISTIC_CC
#define GATEIMAGEWITHSTATISTIC_CC

#include "GateImageWithStatistic.hh"
#include "GateMessageManager.hh"
#include "GateMiscFunctions.hh"

//-----------------------------------------------------------------------------
/// Constructor
GateImageWithStatistic::GateImageWithStatistic()  {
  mIsSquaredImageEnabled = false;
  mIsUncertaintyImageEnabled = false;
  mIsValuesMustBeScaled = false;
  mOverWriteFilesFlag = true;
  mNormalizedToMax = false;
  mNormalizedToIntegral = false;
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
/// Destructor
GateImageWithStatistic::~GateImageWithStatistic()  {
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void GateImageWithStatistic::SetOrigin(G4ThreeVector o) {
  mValueImage.SetOrigin(o);
  mSquaredImage.SetOrigin(o);
  mTempImage.SetOrigin(o);
  mUncertaintyImage.SetOrigin(o);
  mScaledValueImage.SetOrigin(o);
  mScaledSquaredImage.SetOrigin(o);
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void GateImageWithStatistic::SetTransformMatrix(const G4RotationMatrix & m) {
  mValueImage.SetTransformMatrix(m);
  mSquaredImage.SetTransformMatrix(m);
  mTempImage.SetTransformMatrix(m);
  mUncertaintyImage.SetTransformMatrix(m);
  mScaledValueImage.SetTransformMatrix(m);
  mScaledSquaredImage.SetTransformMatrix(m);
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void GateImageWithStatistic::SetScaleFactor(double s) {
  mScaleFactor = s;
  mIsValuesMustBeScaled = true;
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void GateImageWithStatistic::SetResolutionAndHalfSize(const G4ThreeVector & resolution,
						      const G4ThreeVector & halfSize,
						      const G4ThreeVector & position)  {

  mValueImage.SetResolutionAndHalfSize(resolution, halfSize, position);
  if (mIsUncertaintyImageEnabled) {
    mUncertaintyImage.SetResolutionAndHalfSize(resolution, halfSize, position);
    if (!mIsSquaredImageEnabled) {
      mSquaredImage.SetResolutionAndHalfSize(resolution, halfSize, position);
      mTempImage.SetResolutionAndHalfSize(resolution, halfSize, position);
      mScaledSquaredImage.SetResolutionAndHalfSize(resolution, halfSize, position);
    }
  }
  if (mIsSquaredImageEnabled) {
    mSquaredImage.SetResolutionAndHalfSize(resolution, halfSize, position);
    mTempImage.SetResolutionAndHalfSize(resolution, halfSize, position);
    mScaledSquaredImage.SetResolutionAndHalfSize(resolution, halfSize, position);
  }

  mScaledValueImage.SetResolutionAndHalfSize(resolution, halfSize, position);
}
//-----------------------------------------------------------------------------

// A.Resch
//-----------------------------------------------------------------------------
void GateImageWithStatistic::SetResolutionAndHalfSizeCylinder(const G4ThreeVector & resolution,
						      const G4ThreeVector & halfSize, const G4ThreeVector & position)  {
  
  mValueImage.SetResolutionAndHalfSizeCylinder(resolution, halfSize, position);
  if (mIsUncertaintyImageEnabled) {
    mUncertaintyImage.SetResolutionAndHalfSizeCylinder(resolution, halfSize, position);
    if (!mIsSquaredImageEnabled) {
      mSquaredImage.SetResolutionAndHalfSizeCylinder(resolution, halfSize, position);
      mTempImage.SetResolutionAndHalfSizeCylinder(resolution, halfSize, position);
      mScaledSquaredImage.SetResolutionAndHalfSizeCylinder(resolution, halfSize, position);
    }
  }
  if (mIsSquaredImageEnabled) {
    mSquaredImage.SetResolutionAndHalfSizeCylinder(resolution, halfSize, position);
    mTempImage.SetResolutionAndHalfSizeCylinder(resolution, halfSize, position);
    mScaledSquaredImage.SetResolutionAndHalfSizeCylinder(resolution, halfSize, position);
  }

  mScaledValueImage.SetResolutionAndHalfSizeCylinder(resolution, halfSize, position);
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void GateImageWithStatistic::SetResolutionAndHalfSizeCylinder(const G4ThreeVector & resolution,
						      const G4ThreeVector & halfSize)  {
  G4ThreeVector mPosition = G4ThreeVector(0.0, 0.0, 0.0);
   SetResolutionAndHalfSizeCylinder(resolution, halfSize, mPosition);
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void GateImageWithStatistic::SetResolutionAndHalfSize(const G4ThreeVector & resolution,
						      const G4ThreeVector & halfSize)  {
  G4ThreeVector mPosition = G4ThreeVector(0.0, 0.0, 0.0);
  SetResolutionAndHalfSize(resolution, halfSize, mPosition);
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void GateImageWithStatistic::Allocate() {
  mValueImage.Allocate();
  if (mIsUncertaintyImageEnabled) {
    mUncertaintyImage.Allocate();
    if (!mIsSquaredImageEnabled) {
      mSquaredImage.Allocate();
      mTempImage.Allocate();
      if (mIsValuesMustBeScaled) mScaledSquaredImage.Allocate();
    }
  }
  if (mIsSquaredImageEnabled) {
    mSquaredImage.Allocate();
    mTempImage.Allocate();
    if (mIsValuesMustBeScaled) mScaledSquaredImage.Allocate();
  }
  if (mIsValuesMustBeScaled) mScaledValueImage.Allocate();
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void GateImageWithStatistic::Reset(double val) {
  mValueImage.Fill(val);
  if (mIsUncertaintyImageEnabled) {
    mUncertaintyImage.Fill(0.0);
    if (!mIsSquaredImageEnabled) {
      mSquaredImage.Fill(val*val);
      mTempImage.Fill(0.0);
      if (mIsValuesMustBeScaled) mScaledSquaredImage.Fill(0.0);
    }
  }
  if (mIsSquaredImageEnabled) {
    mSquaredImage.Fill(val*val);
    mTempImage.Fill(0.0);
    if (mIsValuesMustBeScaled) mScaledSquaredImage.Fill(0.0);
  }
  if (mIsValuesMustBeScaled) mScaledValueImage.Fill(0.0);
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void GateImageWithStatistic::Fill(double value) {
  mValueImage.Fill(value);
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
double GateImageWithStatistic::GetValue(const int index) {
  return mValueImage.GetValue(index);
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void GateImageWithStatistic::SetValue(const int index, double value) {
  mValueImage.SetValue(index, value);
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void GateImageWithStatistic::AddValue(const int index, double value) {
  GateDebugMessage("Actor", 2, "AddValue index=" << index << " value=" << value << Gateendl);
  mValueImage.AddValue(index, value);
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void GateImageWithStatistic::AddTempValue(const int index, double value) {
  GateDebugMessage("Actor", 2, "AddTempValue index=" << index << " value=" << value << Gateendl);
  mTempImage.AddValue(index, value);
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void GateImageWithStatistic::AddValueAndUpdate(const int index, double value) {

  GateDebugMessageInc("Actor", 2, "AddValue and update -- start: "<<mTempImage.GetSize() << Gateendl);
  double tmp = mTempImage.GetValue(index);
  mValueImage.AddValue(index, tmp);
  if (mIsSquaredImageEnabled || mIsUncertaintyImageEnabled) mSquaredImage.AddValue(index, tmp*tmp);
  mTempImage.SetValue(index, value);
  GateDebugMessageDec("Actor", 2, "AddValue and update -- end"<< Gateendl);
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void GateImageWithStatistic::SetFilename(G4String f) {
  mFilename = f;
  mSquaredFilename = G4String(removeExtension(f))+"-Squared."+G4String(getExtension(f));
  mUncertaintyFilename = G4String(removeExtension(f))+"-Uncertainty."+G4String(getExtension(f));

  mInitialFilename = mFilename;
  mSquaredInitialFilename = mSquaredFilename;
  mUncertaintyInitialFilename = mUncertaintyFilename;
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void GateImageWithStatistic::SaveData(int numberOfEvents, bool normalise) {

  // Filename
  if (!mOverWriteFilesFlag) {
    mFilename = GetSaveCurrentFilename(mInitialFilename);
    mSquaredFilename = GetSaveCurrentFilename(mSquaredInitialFilename);
    mUncertaintyFilename = GetSaveCurrentFilename(mUncertaintyInitialFilename);
  }

  static double factor=1.0;
  if (mIsSquaredImageEnabled || mIsUncertaintyImageEnabled) {UpdateImage();}

  if (mIsValuesMustBeScaled == true) {
    //DD(mScaleFactor);
    factor = mScaleFactor;
  }

  if (normalise) {
    mIsValuesMustBeScaled = true;
    double sum = 0.0;
    double max = 0.0;
    GateImageDouble::const_iterator pi = mValueImage.begin();
    GateImageDouble::const_iterator pe = mValueImage.end();
    while (pi != pe) {
      if (*pi > max) max = *pi;
      sum += *pi*factor;
      ++pi;
    }
    if (mNormalizedToMax) SetScaleFactor(factor*1.0/max);
    if (mNormalizedToIntegral) SetScaleFactor(factor*1.0/sum);
  }

  GateMessage("Actor", 2, "Save " << mFilename << " with scaling = "
	      << mScaleFactor << "(" << mIsValuesMustBeScaled << ")\n");

  if (!mIsValuesMustBeScaled) mValueImage.Write(mFilename);
  else {
    GateImageDouble::iterator po = mScaledValueImage.begin();
    GateImageDouble::iterator pi = mValueImage.begin();
    GateImageDouble::const_iterator pe = mValueImage.end();
    while (pi != pe) {
      *po = (*pi)*mScaleFactor;
      ++pi;
      ++po;
    }
    mScaledValueImage.Write(mFilename);
    SetScaleFactor(factor);
  }

  if (mIsSquaredImageEnabled) {
    UpdateSquaredImage();
    if (!mIsUncertaintyImageEnabled) { // only write if square enable and no uncertainty
      if(!mIsValuesMustBeScaled)  mSquaredImage.Write(mSquaredFilename);
      else mScaledSquaredImage.Write(mSquaredFilename);
    }
  }
  if (mIsUncertaintyImageEnabled) {
    if (!mIsSquaredImageEnabled) UpdateSquaredImage();
    UpdateUncertaintyImage(numberOfEvents);
    mUncertaintyImage.Write(mUncertaintyFilename);
    mSquaredImage.Write(mSquaredFilename); // force output of squared dose for grid
  }

}
//-----------------------------------------------------------------------------


////-----------------------------------------------------------------------------
//void GateImageWithStatistic::SaveData(bool normaliseCylindricalArea) {

  //// Filename
  ////if (!mOverWriteFilesFlag) {
    ////mFilename = GetSaveCurrentFilename(mInitialFilename);
    ////mSquaredFilename = GetSaveCurrentFilename(mSquaredInitialFilename);
    ////mUncertaintyFilename = GetSaveCurrentFilename(mUncertaintyInitialFilename);
  ////}

  ////static double factor=1.0;
  ////if (mIsSquaredImageEnabled || mIsUncertaintyImageEnabled) {UpdateImage();}

  ////if (mIsValuesMustBeScaled == true) {
    //////DD(mScaleFactor);
    ////factor = mScaleFactor;
  ////}

  ////if (normalise) {
    ////mIsValuesMustBeScaled = true;
    ////double sum = 0.0;
    ////double max = 0.0;
    ////GateImageDouble::const_iterator pi = mValueImage.begin();
    ////GateImageDouble::const_iterator pe = mValueImage.end();
    ////while (pi != pe) {
      ////if (*pi > max) max = *pi;
      ////sum += *pi*factor;
      ////++pi;
    ////}
    ////if (mNormalizedToMax) SetScaleFactor(factor*1.0/max);
    ////if (mNormalizedToIntegral) SetScaleFactor(factor*1.0/sum);
  ////}

  ////GateMessage("Actor", 2, "Save " << mFilename << " with scaling = "
	      ////<< mScaleFactor << "(" << mIsValuesMustBeScaled << ")\n");
  //G4ThreeVector  voxSize = mValueImage.GetVoxelSize();
  //G4ThreeVector  resolut = mValueImage.GetResolution();
  //if (!normaliseCylindricalArea) mValueImage.Write(mFilename);
  //else {
	  //int n = 1;
    ////GateImageDouble::iterator po = mScaledValueImage.begin();
    //GateImageDouble::iterator pi = mValueImage.begin();
    //GateImageDouble::const_iterator pe = mValueImage.end();
    //while (pi != pe) {
      //*pi = (M_PI)*voxSize.x()*voxSize.x()*(2*n-1);
      //n++;
      //if (n == resolut.x()) { n =1; }
      //++pi;
    //}
    //mScaledValueImage.Write(mFilename);
    ////SetScaleFactor(factor);
  //}


//}
////-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void GateImageWithStatistic::UpdateImage() {
  GateImageDouble::iterator pi = mValueImage.begin();
  GateImageDouble::iterator pt = mTempImage.begin();
  GateImageDouble::const_iterator pe = mValueImage.end();
  while (pi != pe) {
    *pi += (*pt);
    ++pt;
    ++pi;
  }
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void GateImageWithStatistic::UpdateSquaredImage() {
  GateImageDouble::iterator pi = mSquaredImage.begin();
  GateImageDouble::iterator pt = mTempImage.begin();
  GateImageDouble::const_iterator pe = mSquaredImage.end();
  GateImageDouble::iterator po;

  if(mIsValuesMustBeScaled) po = mScaledSquaredImage.begin();

  double fact = mScaleFactor*mScaleFactor;
  while (pi != pe) {
    *pi += (*pt)*(*pt);
    if(mIsValuesMustBeScaled)  {  *po = (*pi)*fact;++po;}
    *pt = 0;
    ++pt;
    ++pi;
  }
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void GateImageWithStatistic::UpdateUncertaintyImage(int numberOfEvents)
{

  GateImageDouble::iterator po = mUncertaintyImage.begin();
  GateImageDouble::iterator pi;
  GateImageDouble::iterator pii;
  GateImageDouble::const_iterator pe;

  if(mIsValuesMustBeScaled)  pi = mScaledValueImage.begin();
  else pi = mValueImage.begin();

  if(mIsValuesMustBeScaled) pii = mScaledSquaredImage.begin();
  else pii = mSquaredImage.begin();

  if(mIsValuesMustBeScaled)  pe = mScaledValueImage.end();
  else pe = mValueImage.end();

  int N = numberOfEvents;

  while (pi != pe) {
    double squared = (*pii);
    double mean = (*pi);


    // Ma2002 p1679 : relative statistical uncertainty
    /*	if (mean != 0.0)
     *po = sqrt( (N*squared - mean*mean) / ((N-1)*(mean*mean)) );
     else *po = 1;*/


    // Chetty2006 p1250 : relative statistical uncertainty
    // exactly same than Ma2002
    if (mean != 0.0 && N != 1 && squared != 0.0){
      *po = sqrt( (1.0/(N-1))*(squared/N - pow(mean/N, 2)))/(mean/N);
    }
    else *po = 1;

    /*
    // Ma2002 p1679 : relative statistical uncertainty (estimation)
    if (mean != 0.0)
    *po = sqrt( squared/(mean*mean) );
    else *po = 1;
    */

    /*
    // Walters2002 p2745 : statistical uncertainty
    if (mean != 0.0) {
    *po = sqrt((1.0/((double)N-1.0)) *
    (squared/(double)N - pow(mean/(double)N, 2)));
    }
    else *po = 1.0;
    */
    ++po;
    ++pi;
    ++pii;
  }
}
//-----------------------------------------------------------------------------

#endif /* end #define GATEIMAGEWITHSTATISTIC_CC */

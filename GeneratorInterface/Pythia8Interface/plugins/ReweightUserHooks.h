#include "Pythia8/Pythia.h"

class PtHatReweightUserHook : public Pythia8::UserHooks
{
  public:
    PtHatReweightUserHook(double _pt = 15, double _power = 4.5) :
      pt(_pt), power(_power) {}
    virtual ~PtHatReweightUserHook() {}

    virtual bool canBiasSelection() { return true; }

    virtual double biasSelectionBy(const Pythia8::SigmaProcess* sigmaProcessPtr,
                      const Pythia8::PhaseSpace* phaseSpacePtr, bool inEvent)
    {
      //the variable selBias of the base class should be used;
      if ((sigmaProcessPtr->nFinal() == 2)) {
        selBias = pow(phaseSpacePtr->pTHat() / pt, power);
        return selBias;
      }
      selBias = 1.;
      return selBias;
    }

  private:
    double pt, power;
};

class RapReweightUserHook : public Pythia8::UserHooks
{
  public:
    RapReweightUserHook() {}
    virtual ~RapReweightUserHook() {}

    virtual bool canBiasSelection() { return true; }

    virtual double biasSelectionBy(const Pythia8::SigmaProcess* sigmaProcessPtr,
                      const Pythia8::PhaseSpace* phaseSpacePtr, bool inEvent)
    {
      //the variable selBias of the base class should be used;
      if ((sigmaProcessPtr->nFinal() == 2)) {
        double x1 = phaseSpacePtr->x1();
        double x2 = phaseSpacePtr->x2();
        double yLab = 0.5*log(x1/x2);
        double yCM = 0.5*log( phaseSpacePtr->tHat() / phaseSpacePtr->uHat() );
        double pTHat = phaseSpacePtr->pTHat();
        double sigmaLab = 15.44/pow(pTHat,0.0253) - 12.56; // empirical parametrization
        double sigmaCM  = 5.45/pow(pTHat+64.84,0.34);      // empirical parametrization
        selBias = exp( pow(fabs(yLab),2.)/(2*sigmaLab*sigmaLab) + pow(fabs(yCM),2.)/(2*sigmaCM*sigmaCM) ); // empirical reweighting function
        return selBias;
      }
      selBias = 1.;
      return selBias;
    }
};

class PtHatRapReweightUserHook : public Pythia8::UserHooks
{
  public:
    PtHatRapReweightUserHook(double _pt = 15, double _power = 4.5) :
      pt(_pt), power(_power) {}
    virtual ~PtHatRapReweightUserHook() {}

    virtual bool canBiasSelection() { return true; }

    virtual double biasSelectionBy(const Pythia8::SigmaProcess* sigmaProcessPtr,
                      const Pythia8::PhaseSpace* phaseSpacePtr, bool inEvent)
    {
      //the variable selBias of the base class should be used;
      if ((sigmaProcessPtr->nFinal() == 2)) {
        double x1 = phaseSpacePtr->x1();
        double x2 = phaseSpacePtr->x2();
        double yLab = 0.5*log(x1/x2);
        double yCM = 0.5*log( phaseSpacePtr->tHat() / phaseSpacePtr->uHat() );
        double pTHat = phaseSpacePtr->pTHat();
        double sigmaLab = 15.44/pow(pTHat,0.0253) - 12.56; // empirical parametrization
        double sigmaCM  = 5.45/pow(pTHat+64.84,0.34);      // empirical parametrization
        selBias = pow(pTHat / pt, power) * exp( pow(fabs(yLab),2.)/(2*sigmaLab*sigmaLab) + pow(fabs(yCM),2.)/(2*sigmaCM*sigmaCM) ); // empirical reweighting function
        return selBias;
      }
      selBias = 1.;
      return selBias;
    }

  private:
    double pt, power;
};

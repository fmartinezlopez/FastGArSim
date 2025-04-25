#include "PhysicsList.hh"

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4HadronPhysicsFTFP_BERT.hh"
#include "G4HadronPhysicsQGSP_BERT.hh"
#include "G4HadronPhysicsQGSP_BIC.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4IonPhysics.hh"
#include "G4ProcessManager.hh"
#include "G4StepLimiter.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

PhysicsList::PhysicsList()
: G4VModularPhysicsList(),
  fEmPhysicsList(nullptr),
  fHadronPhysicsList(nullptr),
  fDecayPhysicsList(nullptr),
  fMessenger(nullptr),
  fEmName("emstandard_opt4"),
  fHadronName("FTFP_BERT"),
  fCutValue(1.0*mm)
{
    // Define messenger for physics list configuration
    fMessenger = new G4GenericMessenger(this, "/physics/", "Physics List Commands");
    
    // Define command for EM physics variant
    fMessenger->DeclareProperty("emModel", fEmName, 
                                "Set EM physics model (emstandard_opt0-4, emlivermore, empenelope)");
    
    // Define command for hadronic physics variant
    fMessenger->DeclareProperty("hadronicModel", fHadronName,
                                "Set hadronic physics model (FTFP_BERT, QGSP_BERT, QGSP_BIC, QBBC)");
    
    // Define command for cut value
    fMessenger->DeclarePropertyWithUnit("cutValue", "mm", fCutValue,
                                        "Set production threshold");
    
    // Initialize with default physics
    fDecayPhysicsList = new G4DecayPhysics();
    
    // Set EM physics list
    AddPhysicsList(fEmName);
    
    // Set hadronic physics list
    AddPhysicsList(fHadronName);
    
    // Set default cut value
    SetVerboseLevel(1);
}

PhysicsList::~PhysicsList()
{
    delete fMessenger;
    delete fDecayPhysicsList;
    delete fEmPhysicsList;
    delete fHadronPhysicsList;
}

void PhysicsList::ConstructParticle()
{
    // Construct all particles
    fDecayPhysicsList->ConstructParticle();
}

void PhysicsList::ConstructProcess()
{
    // Add standard transportation
    AddTransportation();
    
    // Construct decay processes
    fDecayPhysicsList->ConstructProcess();
    
    // Construct EM processes
    if(fEmPhysicsList) {
        fEmPhysicsList->ConstructProcess();
    }
    
    // Construct hadronic processes
    if(fHadronPhysicsList) {
        fHadronPhysicsList->ConstructProcess();
    }

    // Enable step limitation    
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle;
    G4ProcessManager* pmanager;
    
    auto particleIterator = particleTable->GetIterator();
    particleIterator->reset();
    
    while ((*particleIterator)()) {
        particle = particleIterator->value();
        pmanager = particle->GetProcessManager();
        if (pmanager) {
            pmanager->AddDiscreteProcess(new G4StepLimiter());
        }
    }

}

void PhysicsList::SetCuts()
{
    // Set production threshold for all particles
    SetCutValue(fCutValue, "gamma");
    SetCutValue(fCutValue, "e-");
    SetCutValue(fCutValue, "e+");
    SetCutValue(fCutValue, "proton");

    G4cout << "\n==============================================" << G4endl;
    G4cout << "PHYSICS LIST: Production threshold definitely set to " 
           << G4BestUnit(fCutValue, "Length") << G4endl;
    G4cout << "==============================================" << G4endl;
}

void PhysicsList::AddPhysicsList(const G4String& name)
{
    // Configure EM physics
    if(name == "emstandard_opt0") {
        delete fEmPhysicsList;
        fEmPhysicsList = new G4EmStandardPhysics();
        fEmName = name;
    } else if(name == "emstandard_opt1") {
        delete fEmPhysicsList;
        fEmPhysicsList = new G4EmStandardPhysics_option1();
        fEmName = name;
    } else if(name == "emstandard_opt2") {
        delete fEmPhysicsList;
        fEmPhysicsList = new G4EmStandardPhysics_option2();
        fEmName = name;
    } else if(name == "emstandard_opt3") {
        delete fEmPhysicsList;
        fEmPhysicsList = new G4EmStandardPhysics_option3();
        fEmName = name;
    } else if(name == "emstandard_opt4") {
        delete fEmPhysicsList;
        fEmPhysicsList = new G4EmStandardPhysics_option4();
        fEmName = name;
    } else if(name == "emlivermore") {
        delete fEmPhysicsList;
        fEmPhysicsList = new G4EmLivermorePhysics();
        fEmName = name;
    } else if(name == "empenelope") {
        delete fEmPhysicsList;
        fEmPhysicsList = new G4EmPenelopePhysics();
        fEmName = name;
    }
    
    // Configure hadronic physics
    else if(name == "FTFP_BERT") {
        delete fHadronPhysicsList;
        fHadronPhysicsList = new G4HadronPhysicsFTFP_BERT();
        fHadronName = name;
    } else if(name == "QGSP_BERT") {
        delete fHadronPhysicsList;
        fHadronPhysicsList = new G4HadronPhysicsQGSP_BERT();
        fHadronName = name;
    } else if(name == "QGSP_BIC") {
        delete fHadronPhysicsList;
        fHadronPhysicsList = new G4HadronPhysicsQGSP_BIC();
        fHadronName = name;
    } else if(name == "QBBC") {
        delete fHadronPhysicsList;
        fHadronPhysicsList = new G4HadronInelasticQBBC();
        fHadronName = name;
    } else {
        G4cout << "Physics list " << name << " is not defined." << G4endl;
    }
    
    // Inform the user
    G4cout << "Physics list set to " << name << G4endl;
}
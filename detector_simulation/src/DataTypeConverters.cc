#include "DataTypeConverters.hh"
#include "G4SystemOfUnits.hh"

// Convert Geant4 TPCHit to ROOT TPCHit
root::TPCHit ConvertTPCHit(const TPCHit& g4Hit) {
    return root::TPCHit(
        g4Hit.GetPosition().x()/cm,
        g4Hit.GetPosition().y()/cm,
        g4Hit.GetPosition().z()/cm,
        g4Hit.GetEnergyDeposit()/MeV,
        g4Hit.GetStepSize()/cm
    );
}

// Convert Geant4 ECalHit to ROOT ECalHit
root::ECalHit ConvertECalHit(const ECalHit& g4Hit) {
    return root::ECalHit(
        g4Hit.GetPosition().x()/cm,
        g4Hit.GetPosition().y()/cm,
        g4Hit.GetPosition().z()/cm,
        g4Hit.GetEnergyDeposit()/MeV
    );
}

// Convert Geant4 MuIDHit to ROOT MuIDHit
root::MuIDHit ConvertMuIDHit(const MuIDHit& g4Hit) {
    return root::MuIDHit(
        g4Hit.GetPosition().x()/cm,
        g4Hit.GetPosition().y()/cm,
        g4Hit.GetPosition().z()/cm,
        g4Hit.GetEnergyDeposit()/MeV
    );
}

// Convert Geant4 TrajectoryPoint to ROOT TrajectoryPoint
root::TrajectoryPoint ConvertTrajectoryPoint(const TrajectoryPoint& g4Point) {
    return root::TrajectoryPoint(
        g4Point.position.x()/cm,
        g4Point.position.y()/cm,
        g4Point.position.z()/cm,
        g4Point.time/ns
    );
}

// Convert Geant4 Trajectory to ROOT Trajectory
root::Trajectory ConvertTrajectory(const Trajectory& g4Trajectory) {
    root::Trajectory rootTrajectory;
    
    // Convert volume names
    rootTrajectory.start_vol_name = g4Trajectory.GetStartVolName().c_str();
    rootTrajectory.stop_vol_name = g4Trajectory.GetStopVolName().c_str();
    
    // Convert trajectory points
    const auto& g4Points = g4Trajectory.GetPoints();
    for (const auto& g4Point : g4Points) {
        root::TrajectoryPoint rootPoint = ConvertTrajectoryPoint(g4Point);
        rootTrajectory.points.push_back(rootPoint);
    }
    
    return rootTrajectory;
}

// Convert Geant4 Particle to ROOT Particle
root::Particle ConvertParticle(const Particle& g4Particle) {
    root::Particle rootParticle(
        g4Particle.GetTrackID(),
        g4Particle.GetPDGCode(),
        g4Particle.GetCreatorProcess().c_str(),
        g4Particle.GetMotherID()
    );
    
    // Convert trajectory
    rootParticle.trajectory = ConvertTrajectory(g4Particle.GetTrajectory());
    
    // Convert TPC hits
    for (const auto& g4Hit : g4Particle.GetTPCHits()) {
        rootParticle.tpcHits.push_back(ConvertTPCHit(g4Hit));
    }
    
    // Convert ECal hits
    for (const auto& g4Hit : g4Particle.GetECalHits()) {
        rootParticle.ecalHits.push_back(ConvertECalHit(g4Hit));
    }
    
    // Convert MuID hits
    for (const auto& g4Hit : g4Particle.GetMuIDHits()) {
        rootParticle.muidHits.push_back(ConvertMuIDHit(g4Hit));
    }
    
    // Convert secondary hits
    for (const auto& g4Hit : g4Particle.GetSecondaryTPCHits()) {
        rootParticle.sec_tpcHits.push_back(ConvertTPCHit(g4Hit));
    }
    
    for (const auto& g4Hit : g4Particle.GetSecondaryECalHits()) {
        rootParticle.sec_ecalHits.push_back(ConvertECalHit(g4Hit));
    }
    
    for (const auto& g4Hit : g4Particle.GetSecondaryMuIDHits()) {
        rootParticle.sec_muidHits.push_back(ConvertMuIDHit(g4Hit));
    }
    
    return rootParticle;
}

// Convert Geant4 Event to ROOT Event
root::Event ConvertEvent(const Event& g4Event) {
    root::Event rootEvent(g4Event.GetEventID());
    
    // Convert particles
    for (const auto& g4Particle : g4Event.GetParticles()) {
        rootEvent.particles.push_back(ConvertParticle(g4Particle));
    }
    
    return rootEvent;
}
sc = satelliteScenario;
load("timetableSatelliteTrajectory.mat","positionTT","velocityTT");
sat = satellite(sc,positionTT,velocityTT,"CoordinateFrame","ecef");
play(sc);
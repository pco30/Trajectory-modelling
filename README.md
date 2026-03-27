# Trajectory-modelling

MATLAB-based trajectory simulation project with two model pipelines:
- a ballistic/parachute trajectory model with 2D and 3D visualization
- a thrust-assisted orbit-style ascent model

## Project overview

This repository simulates vehicle motion through Earth's atmosphere using:
- RK4 numerical integration
- gravity and aerodynamic drag
- atmosphere interpolation (density, pressure, temperature)
- shooting methods (secant-style) to match target conditions

Outputs are primarily plots and animations generated in MATLAB.

## Repository structure

```text
Trajectory-modelling/
|- README.md
|- 2D and 3D models/
|  |- EarthProjection.m
|  |- ShootingMethod.m
|  |- ivpSolver.m
|  |- stepRungeKutta.m
|  |- stateDeriv.m
|  |- atmosEarth.m
|- Orbit/
   |- OrbitShootingMethod.m
   |- ivpSolver.m
   |- stepRungeKutta.m
   |- stateDeriv.m
   |- thrust.m
   |- atmosEarth.m
```

## Requirements

- MATLAB (recommended)
- No external package installation required

## Quick start

Open MATLAB in the repository root, then run one of the paths below.

### 1) 2D/3D ballistic trajectory

```matlab
addpath('2D and 3D models');
EarthProjection(200000, 2500, 0.1);
```

Parameters in `EarthProjection(apogee, v, dt)`:
- `apogee`: target apogee (m)
- `v`: initial speed (m/s)
- `dt`: integration time step (s)

### 2) Orbit pipeline

```matlab
addpath('Orbit');
OrbitShootingMethod(2000000);
```

Parameter in `OrbitShootingMethod(alt)`:
- `alt`: target altitude (m)

## How it works (high level)

### 2D and 3D models

1. `EarthProjection` calls `ShootingMethod` to solve for launch angle.
2. `ShootingMethod` repeatedly runs `ivpSolver` and updates angle via secant logic.
3. `ivpSolver` integrates the ODE using RK4 (`stepRungeKutta` -> `stateDeriv`).
4. `stateDeriv` uses atmospheric data from `atmosEarth` and applies drag/gravity.
5. Results are rendered in 2D and projected on a 3D Earth surface.

### Orbit

1. `OrbitShootingMethod` solves for a thrust-start condition.
2. `ivpSolver` propagates state and mass across pre-thrust/burn/post-burn phases.
3. `stateDeriv` combines gravity, drag, and `thrust`.
4. Atmospheric properties are provided by `atmosEarth`.

## Notes and known issues

- There is no automated test suite yet; validation is currently plot/trajectory based.
- The orbit shooting path may require debugging/tuning for convergence in some cases.
- `EarthProjection` fetches an Earth texture image from a remote URL at runtime.

## Development suggestions

- Keep simulation constants grouped and documented for easier calibration.
- Add lightweight regression scripts to compare key outputs after changes.
- Consider reducing duplicated logic between `Orbit/` and `2D and 3D models/`.

# flow_properties

The flow along a trajectory: from the `tecplot.dat` a run has written, the Mach, Reynolds and
Knudsen numbers, the pressure and the dynamic pressure, the stagnation-point heat flux and the
heat load, in one Tecplot file. **It does not run the solver**, and needs numpy and PyYAML only.

```console
cd validation/pathfinder
python3 ../../src_helper/flow_properties/flow_properties.py output_result/tecplot.dat \
        --nose-radius 0.6625 -o tecplot_output.dat
```

It writes the same quantities as `src_helper/compute_flowfield_fromTrajectory` of Tacode v1,
whose `tecplot_output.dat` it replaces, one row per point of the trajectory:

| Column | Meaning |
|---|---|
| `Time[s]`, `Alti[km]`, `Dens[kg/m3]`, `Temp[K]` | copied from the run |
| `Pres[Pa]` | `ρ R T`, `R = R_u/M` |
| `Velair[m/s]` | the air-relative speed: `VelairAbs` when the run had a wind, `VelplAbs` otherwise |
| `Mach` | `V/sqrt(γ R T)` |
| `Re` | `ρ V L/μ`, `μ` from the Sutherland law |
| `Kn` | **copied from the run** (built by the solver from the species of the atmosphere table and `satellite.characteristic_length`) |
| `Qdyn[Pa]` | `ρ V²/2` |
| `HeatFlux[W/m2]` | stagnation-point convective heat flux, `C ρ^a V^b/sqrt(Rn)` |
| `HeatLoad[J/m2]` | its time integral (trapezoidal) |

## Where the numbers come from

- **`γ` and `M`** (`--specific-heat-ratio`, `--molecular-weight`, g/mol): given on the command
  line or in the settings file; otherwise from `atmosphere.specific_heat_ratio` and
  `atmosphere.molecular_weight` of the run's control file (`--control`, default `config.yml`),
  which are the values the solver uses with a Mach table — so the `Mach` column is the solver's
  own when the run has one; otherwise the default of `--gas`.
- **`L`** of the Reynolds number (`--length`): otherwise `satellite.characteristic_length`.
- **`--gas`** (`co2`, the default, or `air`) gives the Sutherland law, the Sutton-Graves
  coefficient, and `γ`, `M` when nothing else does:

  | | `γ` | `M` g/mol | Sutherland `μ0`, `T0`, `S` | Sutton-Graves `k` kg^0.5/m |
  |---|---|---|---|---|
  | `co2` | 1.290723 | 43.49 | 1.37e-5 Pa s, 293 K, 240 K | 1.9027e-4 |
  | `air` | 1.4 | 28.9644 | 1.716e-5 Pa s, 273.15 K, 110.4 K | 1.7415e-4 |

  The `co2` values are those of the v1 tool (CO2 95.57 %, N2 2.70 %, Ar 1.60 %, O2 0.13 %,
  mixed with 7, 5, 3 and 5 degrees of freedom); the `air` ones are the standard atmosphere and
  White, *Viscous Fluid Flow*.
- **`--heating`**: `sutton-graves` (the default; `C = k`, `a = 0.5`, `b = 3`), `tauber-mars`
  (`1.35e-4`, `0.5`, `3.04`: Tauber, Bowles and Yang 1989, whose `1.35e-8` is in W/cm²), or
  `custom` with `--heat-flux-law C a b`. `--nose-radius` is `Rn` in metres. **Without it the
  tool uses 1 m and says so**, on the screen and in the header of the file it writes (`NOT
  GIVEN: a placeholder`): 1 m is no more than a placeholder, and the heat flux scales as
  `1/sqrt(Rn)`.

## Differences from the v1 tool

- The v1 tool took the heat-flux coefficient of Tauber as `1.35e-7`, which gives kW/m², and
  labelled the column W/m². This one writes W/m².
- The v1 header had no comma before `"HeatFlux[W/m2]"` and no name for the heat-load column it
  wrote. Here the header names every column, in the form Tacode's own output uses, so that
  `src_helper/animate_trajectory/tecplot_reader.py` and Tecplot both read it.
- The v1 tool rebuilt the Knudsen number with a single molecular diameter of 4e-10 m; this one
  uses the solver's.
- The v1 tool used the planet-relative speed; this one uses the air-relative speed.

## Settings file

Section `flow_properties` of `config_helper.yml` (see `src_helper/general/README.md`), for example:

```yaml
flow_properties:
  filename: output_result/tecplot.dat
  output: tecplot_output.dat
  nose_radius: 0.22
  heating: tauber-mars
```

`test/test_flow_properties.py` checks the formulas against the solver's Mach number, the
Sutton-Graves expression, the Sutherland law and the heat load of a constant flux, the order in
which the settings are taken, and that the file it writes is read back.

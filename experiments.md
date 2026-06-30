# GIAMIP Experiments

Summary of the GIAMIP experiment protocol. There are two categories of experiments:

- **Tier 1 (mandatory)** — `Exp01`–`Exp05`: a common, relatively controlled set of comparisons characterising inter-model differences from ice history, Earth structure, rheology, initial topography, and numerical implementation.
- **Tier 2 (optional)** — `Exp06`–`Exp12`: additional configurations and sensitivities. `Exp06`–`Exp07` repeat `Exp01`/`Exp02` over the shorter 80 ka–future interval to test sensitivity to loading-history length.

Each experiment is driven by a different **global input ice model**, which determines the forcing file passed to the checker via `--forcing-filepath` (e.g. `Exp05` uses PaleoMIST). The model- and full-level checkers select this file automatically per experiment via the `EXPERIMENT_FORCING` mapping in `giamip_compliance_checker.py`. Optional experiments may submit multiple members using the `ExpXX_mYY` naming convention (`m01` = baseline).

| Tier | Experiment | Global input ice model (forcing) | Viscosity structure | Elastic structure | Earth rheology | Initial polar topography¹ | Period |
|------|-----------|----------------------------------|---------------------|-------------------|----------------|---------------------------|--------|
| 1 (mandatory) | `Exp01` | ICE-6G_D | VM5a | PREM | Maxwell | Bedmap2 (ANT) + BedMachineV5 (Greenland) | 122 ka – 1000 yrs into future |
| 1 (mandatory) | `Exp02` | ICE-7G_NA | VM7 | PREM | Open | Bedmap2 (ANT) + BedMachineV5 (Greenland) | 122 ka – 1000 yrs into future |
| 1 (mandatory) | `Exp03` | GLAC3b (Profile 1) | Within suggested range | Open | Open | Bedmap3 (ANT) + BedMachineV5 (Greenland) | 122 ka – 1000 yrs into future |
| 1 (mandatory) | `Exp04` | GLAC3b (Profile 2) | Within suggested range | Open | Open | Bedmap3 (ANT) + BedMachineV5 (Greenland) | 122 ka – 1000 yrs into future |
| 1 (mandatory) | `Exp05` | PaleoMIST (Version a1) | Gowan et al. (2021) | PREM | Maxwell | BedMachineV3 (ANT) + BedMachineV5 (Greenland) | 80 ka – 1000 yrs into future |
| 2 (optional) | `Exp06` | ICE-6G_D | VM5a | PREM | Maxwell | Bedmap2 (ANT) + BedMachineV5 (Greenland) | 80 ka – 1000 yrs into future |
| 2 (optional) | `Exp07` | ICE-7G_NA | VM7 | PREM | Same as Exp02 | Bedmap2 (ANT) + BedMachineV5 (Greenland) | 80 ka – 1000 yrs into future |
| 2 (optional) | `Exp08` | PaleoMIST | Open | Open | Open | Open | 80 ka – 1000 yrs into future |
| 2 (optional) | `Exp09` | GLAC3b (Profile 1) | Open | Open | Open | Open | 122 ka – 1000 yrs into future |
| 2 (optional) | `Exp10` | GLAC3b (Profile 2) | Open | Open | Open | Open | 122 ka – 1000 yrs into future |
| 2 (optional) | `Exp11` | GLAC3b (Profile 3) | Open | Open | Open | Open | 122 ka – 1000 yrs into future |
| 2 (optional) | `Exp12` | ICE-6G_D | Open | Open | Open | Open | 122 ka – 1000 yrs into future |

¹ Initial polar topography as listed; GEBCO 2024 is used for the rest of the globe.

**Notes**
- *Viscosity / Elastic / Rheology*: prescribed profiles (`VM5a`, `VM7`, `Gowan et al.`) follow the specified profiles in the protocol; "Within suggested range" allows group-selected parameters within the prescribed ranges; "Open" allows free choice. All choices must be documented in the model report.
- *Period*: negative time stamps are past, positive are future, with `0` = present day (2000 CE).

*Source: GIAMIP Experiment Protocol, v3.0 (June 10, 2026), "GIAMIP Experiment Protocol → Experiments" table.*

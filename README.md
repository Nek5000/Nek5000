# Notes: Important Solver updates

## b15364f

- Pressure solver updates. Significantly improves robustness at high CFL

## 0670afa

- constrainTLSR routine updated to prevent TLSR update if the sign switches
- Earlier technique was to flip the sign of computed TLSR

### b2859c8

- Auto dt and number of time steps calculation for CLSR and TLSR
- Based on tests performed on damBreak and circVortex problems

### 7a5b510

- TLSR now does not use the tanh profile for sign function. It was creating boundary related problems and spurious currents at interface.
- There are still some issues at the boundary for extreme cases (see clsTest3D case). 

# Log of Algorithm Testing

## Model setup

Variable naming for model/data specifications:
- `p`: the order or AR model, i.e. AR(p)
- `T`: the length of each time series
- `nGroup`: the number of groups/models of the data
- `nLim`: the number of time series per group, if being an integer; or the
  lower and upper bound of the number of time series per group, if being a
  2-dim vector.
- `tol`: tolerance for EM in numerical computation (default: `1e-8`)

## k-ARs

### Record 1

```
p = 5;
T = 250;
nGroup = 160;
nLim = 100;
```

Results:
```
> elapsed time: 129.210996 seconds
> 546 time series misclustered
> 29 groups missed
```

### Record 2

```
p = 5;
T = 250;
nGroup = 160;
nLim = 150;
```

Results:

```
> Elapsed time is 213.942913 seconds.
> 6685 time series misclustered
> 27 groups missed
```

## mix-ARs

### Record 2

```
p = 5;
T = 250;
nGroup = 160;
nLim = 150;
```

Results:

To compile test problem, `dtst2c.f`

```console
gfortran -o myprogram dtst2c.f dsrc2c.f jcg.f jsi.f rscg.f rssi.f sor.f ssorcg.f ssorsi.f -lblas
```

Additional compilation notes

+ Rename `IRAND` to `IRANDBLH` (due to name conflict)
+ Changed line 1 to `PROGRAM ITPTST` (removed args)

# Welcome to HistrogramsFit.jl's documentation

The problem of fitting curves to histograms is ubiquitous in High Energy Physics (HEP)
and usually it involves three steps:

1. Determining the "best fit" parameters of a curve;
2. Determining the errors on the parameters;
3. Judging the goodness of the fit. 

This simple Julia module takes an $nD$ histogram and a data distribution model and
creates a theoretically sound chi square statistics that can be used to perform

1. Point estimation;
2. Confidence interval estimation;
3. Goodness-of-fit testing.

```
!!! warning
        This is still a work in progress and it is in pre-alpha stage.
        Expect issues and bug. If you want to contribute, please feel free
        to open a issue/start a discussion on GitHub! Thank you!
```

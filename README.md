# cpppPrototipe
a prototipe R package for cppp estimation.


### Two main functions
`calcPPP()` function that calculates the posterior predictive p-value (ppp).

Arguments:
- `data`: the data used to fit the model (in nimble this information is stored in the model object, so we do not neet to explictly pass the object).
- `samplerPosterior`: a function that samples from the posterior distribution of the model.
- `samplesPostPred`: a function that samples from the posterior predictive distribution of the model.
- `discrepancy` : a function that simulates from the discrepancy

# fitness-assay-python

This repository holds code for the _bulk fitness assay_ described in this 
[2016 Cell Paper](http://dx.doi.org/10.1016/j.cell.2016.08.002). In short it infers fitnesses from barcode abundance data, 
specifically designed for serial growth/dilution experiments.

## Experimental scope

The methods implemented here are for inferring the (relative) fitness of 
individuals in a population undergoing serial growth and 
dilution. The methods should also work for inferring fitness in other situations, but have not been tested/optimized for
those conditions. It is important to note that fitness is not an absolute quantity, especially in experiments with
near constant population sizes; instead, it is differences in fitness between organisms which can be measured.

This code is primarily for use in analyzing data from _bulk fitness assays_, where mutants are competed against an
ancestor or known reference strain in order to estimate their fitnesses relative to the ancestor. 
[DNA barcoding technology](http://www.nature.com/nature/journal/v519/n7542/abs/nature14279.html), among other tools,
can be used to measure abundances of different lineages throughout the course of an experiment/assay. This abundance
data can then be translated to information about the relative growth rates of different organisms. This particular
method measures the fitness effects of pre-existing variation and is optimized for competetive assays with
a short duration (so no further evolution occurs during the experiment).

This method is designed around the assumption that
the ancestral/reference strain is the largest component of the pool of organisms throughout
the experiment. It can either be given a list of putatively neutral (relative to the reference) strains, and it will
further filter that list to the "true" neutrals. If no list is given, then it attempts to isolate the majority phenotype
and mark those as the reference strain.

The assumptions of the model tend to break down as mutant populations become larger, due to the advent
of frequency dependent selection (that is, ecology). In cases where the model assumptions break down, the methods can still
be used for quantification of different types; however, the biological meaning of the reported "fitness" may be different,
and the variability across experiments may be much larger than suggested by the reported errors.

## Basics of inference algorithm

This code infers the fitnesses of barcoded individuals from the abundance data. For each pair of timepoints, the fitness
increase _s_ per timepoint is estimated as the log ratio of frequencies at adjacent timepoints, compared
to the reference type. The fitnesses from each
pair of timepoints is assigned a noise value. This noise model assumes Gaussian noise (which was theoretically
motivated and fit well the original data). The code then takes an appropriate weighted average of all measurements
for a single barcode to come up with an overall estimate of the fitness (and its error).

A few notes on the noise model: there are two contributions to variation in fitness measurements.
One part, the "counting noise" comes from genetic drift and sampling error, and the
other part an impossible to remove "multiplicative noise" which is some unknown function of the DNA extraction
and sequencing
pipeline. The algorithm uses within replicate data to estimate the counting noise coefficient, and uses
between replicate comparisons to estimate the multiplicative noise.

The algorithm also attempts to identify cells of the reference (often ancestral) type. There are two ways for
this inference to proceed. If there are known reference barcodes (or a family with high probability of being reference),
the algorithm can use them as the reference type (after cleaning up any adaptive or maladaptive outliers). Otherwise, it 
assumes that most lineages have the reference phenotype, and tries to isolate them by iteratively estimating the
mean and variance of the largest peak in fitness values from the data. These "neutral" lineages are then used
as a reference to measure fitness against, and to estimate the scale of the counting noise. Around
~300 of these lineages are needed in order to get good estimates of the mean fitness and noise coefficients.
(Note: the algorithm tries to center the distribution of the reference strain at 0 fitness, but does not force
it to do so. If the data is good, the centering is close to perfect; however small deviations always exist but are
not meaninful.)

Note that all fitnesses are returned as fitnesses per cycle. This notation is appropriate since the fitness gain
in serial
growth-dilution experiments often comes from dynamics at the beginning and end of the cycle, rather than being
evenly distributed across the cycle (as would be the case for an increased exponential growth rate).

For more details, please see the supplement of the above paper. It is recommended that the user familiarize themselves
with the basic ideas in the method to understand how to interpret results, debug, and choose
parameters. Details available upon request.

## Usage

The main method to be called is `inferFitness`. This function does fitness inference for one experiment across
all replicates at once. It requires the following arguments:

* `barcodes`: _N_ x 1 list of all barcodes. These need to be unique, sortable identifiers in any format (numbers, strings, etc.).
* `cycleTimes`: 1 x _q_ list of cycle times. For example, if the assay was run for three complete cycles and samples were taken before the experiment and after each cycle, except the first, the list would be: [0,2,3]
* `allReads`: Python dictionary of length r, where key is replicate name, and value is N x q matrix of reads.

It the following optional arguments:

* `outputFolder`: Folder to save fitness data output. Saved in tab separated column formats.
* `experimentName`: Name of experiment, used in naming saved fitness files.
* `neutralBarcodes`: List of putatively neutral barcodes. Useful when ancestral barcodes are less abundant, or when fitness distribution is multimodal.
* `multNoiseThresh`: Read threshold used to compute multiplicative noise coefficient.
* `zCutoff`: Cutoff used to establish neutrals in Gaussian model.
* `multNoiseBase`: Default value of multiplicative noise, used if only 1 replicate exists (or multiplicative noise can't be computed).
* `lowCoverageThresh`: Minimum number of reads for a timepoint to be considered.
* `sparsityThresh`: Maximum fraction of lineages which can have 0 reads.

The code returns the following:
* `repFitnessData`: dictionary with same keys as allReads. values are dictionaries, with key value pairs:
  * `neutralBarcodes`: _N<sub>neut</sub>_ x 1 list of barcodes which were assumed to be neutral.
  * `timePointsUsed`: 1 x _q_ list of timepoints used in inference
  * `multNoiseParams`: 1 x _q_-1 list of multiplicative noise parameters
  * `kappas`: 1 x _q_-1 list of additive noise parameters
  * `meanFitnesses`: 1 x _q_-1 list of mean fitnesses
  * `allTimeFitness`: _N_ x _q_-1 list of fitnesses
  * `allTimeErrors`: _N_ x _q_-1 list of errors
  * `aveFitness`: _N_ x 1 list of overall fitness estimates
  * `aveError`: _N_ x 1 list of overall error estimates

In addition, the code prints out a simple consistency check for the noise model at each timepoint for each replicate.
Briefly, this checks if the mixed noise model holds; if timepoints are not consistent,
either `multNoiseThresh` is too low, or `multNoiseBase` should be set to zero. Please see the supplement of
[the original work](http://dx.doi.org/10.1016/j.cell.2016.08.002) for details.

Finally, the function saves files to `outputFolder`, labelled by the experiment name, replicate, and the type of data stored.

## Common questions/errors

### I get a lot of warnings at runtime. Is that normal?

Division by zero warnings are to be expected; the come for lineages which have 0 reads at some timepoints. Depending on
your version of Python, there may be a few `FutureWarnings` as well. Other warnings may be indicative of flawed behavior.

The code should work for both Python 2 and Python 3; however, it is more regularly tested on Python 3. Sometimes
warnings/errors can be explained by the differences between the two.

### I have inhomogenous cycles/want fitness in different units. How can I do that?

The `cycleTimes` parameter can take time in whichever units you desire. If you want per generation fitnesses, then
use the sequence of total generations since first passage. If you have longer or shorter cycles, then
`cycleTimes` should reflect that fact.

### My fitness/error values are all `NaN`. What happened?

If you only have a few `NaN` values in `allTimeFitness`, then those are likely just lineages with poor coverage. If
all values at any timepoint are `NaN`, then something has gone wrong with the estimation of the mean fitness; check the
lineages that are labelled as neutral at input and output to make sure those classes behave correctly.

If the values of `allTimeErrors` are `NaN`, there are two possibilities. One is that the estimation of the additive
noise parameters `kappas` went awry, which is an issue with the neutral types. The other is that estimation of
the `multNoiseParams` is somehow flawed; check the behavior of lineages with average number of reads larger than
`multNoiseThresh`.

Note that if the errors are all `NaN`, the averaged fitnesses will be as well. This is because the averages are
inverse variance weighted.

### The code has suggested that some of my timepoints are inconsistent. What does that mean?

As per [this paper](http://dx.doi.org/10.1016/j.cell.2016.08.002), this code uses a mixed noise model to analyze the data.
At lower read numbers (below 10<sup>3</sup>), the variation in read numbers seems to be well approximated by
an _additive noise model_, where each cell (or read) contributes an independent amount to the errors. These
give errors in fitness which scale as as _r_ <sup>-1/2</sup>, where _r_ is the typical number of reads for that
lineage. This noise comes from drift, sampling, and other error processes.

However, high frequency lineages have noise in fitness which does not decrease with
increasing read depth. This _multiplicative noise_ has unknown origin; it has yet to be studied systematically, but
there are some indications that it may come from biases in sequencing.

The algorithm uses the neutral lineages to estimate the additive noise (by averaging _within_ replicates) and
high frequency lineages to estimate the multiplicative noise (by averaging _between_ replicates). Ideally,
one would like to use lineages whose additive noise is much less than their multiplicative noise in order to estimate
the multiplicative noise.

This gives us a simple consistency check; if the reported multiplicative noise is larger than the additive noise of the
lineages used to estimate it, the noise model makes sense. This gives us the consistency condition
`multNoiseParam > kappa/multNoiseThresh`, where `kappa` is the additive noise parameter. The code computes this
inequality, and reports timepoints which fail it.

If you have inconsistent timepoints, they may not be estimating the multiplicative noise well. You can then try
increasing `multNoiseThresh` to see if this solves the problem. If you wish to ignore multiplicative noise all
together, then set `multNoiseThresh` higher than the max of your data, and `multNoiseBase` to zero. Note that this
may underestimate errors for high frequency lineages.


## Dependencies

Requires Python 2.7.x or 3.5.x or newer. Also requires the [SciPy](https://www.scipy.org/) stack.

## Acknowledgments

Methods developed by Atish Agarwala under the supervision of Daniel S. Fisher. Python implementation by Atish Agarwala.
Additional feedback and support from the Sherlock and Petrov Labs at Stanford.

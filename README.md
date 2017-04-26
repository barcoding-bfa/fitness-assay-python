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

* `barcodes`: N x 1 list of all barcodes. These need to be unique, sortable identifiers in any format (numbers, strings, etc.).
* `cycleTimes`: 1 x q list of cycle times. For example, if the assay was run for three complete cycles and samples were taken before the experiment and after each cycle, except the first, the list would be: [0,2,3]
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
  * `timePointsUsed`: 1 x q list of timepoints used in inference
  * `multNoiseParams`: 1 x q-1 list of multiplicative noise parameters
  * `meanFitnesses`: 1 x q-1 list of mean fitnesses
  * `allTimeFitness`: N x q-1 list of fitnesses
  * `allTimeErrors`: N x q-1 list of errors
  * `aveFitness`: N x 1 list of overall fitness estimates
  * `aveError`: N x 1 list of overall error estimates

In addition, it saves files to `outputFolder`, labelled by the experiment name, replicate, and the type of data stored.

## Dependencies

Requires Python 2.7.x or 3.5.x or newer. Also requires the [SciPy](https://www.scipy.org/) stack.

## Acknowledgments

Methods developed by Atish Agarwala under the supervision of Daniel S. Fisher. Python implementation by Atish Agarwala.
Additional feedback and support from the Sherlock and Petrov Labs at Stanford.

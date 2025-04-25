# tosnipper

## Description
This function takes in a genotype file and converts it to a SNIPPER 3.5 analysis-ready file. A metadata of the samples is required.

## Function
```
tosnipper(input, references, target.pop = TRUE, population.name = population, markers = snps)
```

## Parameters
@param *input* is the input file where the samples are in the first column and the allelic genotypes are in the succeeding columns. See **sampledata.csv** as an example.

@param *references* is the metadata of the samples included in the input file. The sample codes **should be the same** in both the input and reference file.

@param *target.pop* should be indicated TRUE if there is a target population to be classified. Set to FALSE if all populations will be used as training data.

@param *population.name* should be indicated if target.pop is TRUE.

@param *markers* is the number of markers in the dataset.

## Examples
```
tosnipper("excelfile.xlsx", "references.csv", target.pop = TRUE, population.name = "Filipino", markers = 40)
tosnipper("excelfile.csv", "metadata.xlsx", target.pop = FALSE, markers = 32)
```

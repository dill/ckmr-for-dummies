---
title: Kin pair matching for CKMR
author: David L Miller

---

# Dave doesn't know genetics

Some definitions because I don't know about genetic

- **allele**: varying bit of the gene
- **locus**/**loci**: location (of a region?) on the chromosome
- **Medelian exclusion**: are a given parent-offspring pair compatible at all to be matched, some are not (e.g., `AA`-`AB`, yes but `BB`-`AA`, no).


# Doing CKMR

There are three stages to CKMR:

1. *genotyping*: actually getting the biological samples and coding them into something that the computer can read. This is done by "elves" according to MVB -- soem biological thing that takes the samples, cuts up the DNA and gets interesting sections (alleles/loci/SNPs).
2. *kin pair finding*: working out whether a given pair of animals are related by looking at whether there are matches in alleles. "Kin" can refer to a number of different relationships, only some of which are "interesting".
3. *close-kin mark-recapture*: actually doing the CKMR, thinking about population dynamics etc.

# But what are "genotypes"?

6 possible genotypes per individual per locus: `AA`, `AB`, `BB`, `OO`, `AO`, `BO`.

# What about doing the pair finding?

There are (at least) 4 different pairings that are interesting and some others that are not:

- *duplicates* (DUP): if there was some problem with the samples at stage 1 above or if the same individual gets sampled more than once.
- *parent-offspring* (POP): what you'd expect.
- *full-sibling* (FSP): relatively unlikely, depending on the population size and mating strategy.
- *half-siblings* (HSP): if you have 1 parent in common.
- *cousins*: different demographics, probably want to exclude?
- *unrelated pairs* (UP): animals that are not related, don't want to do anything with them, but do want to know who they are!


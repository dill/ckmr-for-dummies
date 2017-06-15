---
title: Kin pair matching for CKMR
author: David L Miller
mathjax-url: https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js

---

# Notational notes

- The words "fish" and "animal" (and probably, later, "tuna") are used interchangably.
- CKMR is close-kin mark-recapture.

# Dave doesn't know genetics

Some definitions because I don't know about genetics

- **locus**/**loci**: contiguous set of `G`s, `A`s, `T`s, `C`s, in a given location along the chromosome. *Think*: regular expression; need to match a pre/post-amble thing and the bit in the middle has some variation (at just one position?). The trivially largest loci is the whole genome. The smallest one is a bit more tricky.
- **allele**: a varying (at the population level) locus. There may be multiple alleles in the population but each individual will carry at most 2 of them.
- **Medelian exclusion**: are a given parent-offspring pair compatible at all to be matched, some are not (e.g., `AA`-`AB`, yes because you can get the `A` from either parent but `BB`-`AA`, no since you only have `B`s, so where could you get the A from?).
- **linkage**: issue that there is some correlation along the genome between alleles, so if they are physically close, you are more likely to inherit them.


# Doing CKMR

There are three stages to CKMR:

1. *genotyping*: actually getting the biological samples and coding them into something that the computer can read. This is done by "elves" according to MVB -- soem biological thing that takes the samples, cuts up the DNA and gets interesting sections (alleles/loci/SNPs).
2. *kin pair finding*: working out whether a given pair of animals are related by looking at whether there are matches in alleles. "Kin" can refer to a number of different relationships, only some of which are "interesting".
3. *close-kin mark-recapture*: actually doing the CKMR, thinking about population dynamics etc.

# But what *are* "genotypes"?

6 possible genotypes per individual per locus: `AA`, `AB`, `BB`, `OO`, `AO`, `BO`. These vary between the pairs ("copies"). Here `A` and `B` represent whole loci (e.g., "`GCCCTACCTA`") and `O` represents `NULL`, where the code could not be read.

# What about doing the pair finding?

There are (at least) 4 different pairings that are interesting and some others that are not:

- *duplicates* (DUP): if there was some problem with the samples at stage 1 above or if the same individual gets sampled more than once.
- *parent-offspring* (POP): what you'd expect.
- *full-sibling* (FSP): relatively unlikely, depending on the population size and mating strategy.
- *half-siblings* (HSP): if you have 1 parent in common.
- *cousins*: different demographics, probably want to exclude?
- *unrelated pairs* (UP): animals that are not related, don't want to do anything with them, but do want to know who they are!

## Data

Before talking about kin finding, how is the data formatted? At different stages, the data are formatted differently, roughly speaking:

1. Fresh from the sequencer we have a matrix with columns as the fish and rows as the alleles (for a given loci there might be many alleles). Entries are the frequency (number of "reads") of the allele in the given fish.
2. Process to get presence absence of the alleles (`A`, `B` or `O`). This is a bit complicated and involves processing multiple alleles down to the most useful two (plus `O`).
3. Estimate the population levels of each allele.
4. Format to fish as rows and columns as loci. Entries are then then the genotype of that loci for that fish (e.g, `AB`).

## Getting kin pairs

Basic idea here is calculate a statistic between a pair of fish to evaluate the relationship (more on that below), then plot a histogram of the pair scores across the samples. Humps/separations in the histogram show particular relationships. This will get more obvious below.

We remove the pairings *in order*.

### Duplicates

To find duplicates we calculate the number of loci that are the same in a fishpair -- actually number of exclusions. The number of exclusions for a fish with itself is 0, but there is some error in this. In any case the duplicates pile themselves up near zero, the rest are in a big lump further down the horizontal axis (this includes all the other relationships).

### Full sibling pairs

Not many of these, but need to remove them.

### Parent-offspring pairs

POPs are found using *weighted pseudo-exclusions*, WPSEX. We can re-group the genotypes to have, e.g., `AO` and `AA` be in the same group `AAO` if we're not sure. In that case there are some pseudo-exclusions to think about. If we had `AAO` and `BBO` then `AO`-`BB`, `AA`-`BO` and `AA`-`BB` are not allowed (for PO or OP relations) but we don't know what these are (e.g., if `O` are rare), so we don't exclude `AAO` `BBO` in that case (hence "pseudo").

### Half-siblings

HSPs require a different metric again, pseudo-log-odds: PLODs. We calculate this with something like


$$
\text{LOD}_\text{HSP, UP} = \log \frac{\mathbb{P}\left[\tilde{g}_1, \tilde{g}_2 \vert \text{HSP} \right]}{\mathbb{P}\left[\tilde{g}_1, \tilde{g}_2 \vert \text{UP}\right]}
$$

for genotypes $g_1$, $g_2$, recorded as $\tilde{g}_1$, $\tilde{g}_2$ (strictly speaking there should be a subscript $l$ on these for the given loci, but that's dropped for notation simplicity). We can write the above about any kin relationship combination, here we separate HSP from UP though. The pseudo comes from suming over the loci, giving us a score for an pairing.

Can calculate $\mathbb{P}\left[\tilde{g}_1, \tilde{g}_2 \vert \text{HSP} \right]$ as:

$$
\mathbb{P}\left[\tilde{g}_1, \tilde{g}_2 \big\vert \text{HSP} \right] = \sum_{g_1, g_2} \left\{ \mathbb{P}\left[\tilde{g}_1 \big\vert g_1 \right] \mathbb{P}\left[\tilde{g}_2 \big\vert g_2 \right] \mathbb{P}\left[\tilde{g}_1, \tilde{g}_2 \big\vert \text{HSP} \right]\right\}
$$

where "HSP" could be whatever kinship relationship we are interested in (thinking more generally here), we can call this $k_{12}$, say. In that case we can calculate:


$$
\mathbb{P}\left[\tilde{g}_1, \tilde{g}_2 \big\vert \text{HSP} \right] = \kappa_0  \mathbb{P}\left[g_1 \right] \mathbb{P}\left[g_2 \right] + \kappa_1 \mathbb{P}\left[g_1, g_2 \big\vert \text{1 allele shared} \right] + \kappa_2 \mathbb{P}\left[g_1 \right] \mathbb{I}\left[g_1 = g_2 \right]
$$
Thinking simply about the above equation, we're looking at the probabilities that 0, 1 and 2 alleles are shared. Also worth thinking here that *all* inheritence is from parent-offspring pairs (somewhere down the line; this seems obvious but wasn't immediately to me, esp. given equation 5.1 of Bravington, Skaug and Anderson, 2016). These are identical by descent (*idb*) In the above:
$$
\kappa_m = \mathbb{P} \left[ m \text{ shared alleles} \big\vert k_{ij} \right]
$$
where appropriate values can be found in Table 1 of Bravington, Skaug and Anderson, 2016; and:
$$
\mathbb{P}\left[g_i \right]  = \pi(g_{i,a}) \pi(g_{i,b}) \left(1 + \mathbb{I}\left[g_{i,a} \neq g_{i,b} \right] \right)
$$
where $\pi(g_{i,a})$ is the population frequency of $g_{i,a}$ (the above expression assumes "Hardy-Weinberg equilibrium without linkage disequilibrium".

With all of that calculated the PLOD is then:
$$
\text{PLOD}_\text{HSP, UP} = \sum_l \log \frac{\mathbb{P}\left[\tilde{g}_{l,1}, \tilde{g}_{l,2} \vert \text{HSP} \right]}{\mathbb{P}\left[\tilde{g}_{l,1}, \tilde{g}_{l,2} \vert \text{UP}\right]}
$$
where the summation is over loci, indexed by $l$.



# References

Bravington, M.V., Skaug, H.J. and Anderson, E.C. (2016) Close-Kin Mark-Recapture. Statistical Science.



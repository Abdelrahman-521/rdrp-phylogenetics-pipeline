# ClustalW vs T-Coffee: what changes in the RdRp tree?

38 viral RdRp domains, aligned two ways, each alignment turned into a maximum-likelihood tree with IQ-TREE (ModelFinder + 1,000 ultrafast bootstraps).

| | ClustalW | T-Coffee |
|---|---|---|
| Alignment columns | 634 | 697 |
| Gaps | 21.3% | 28.4% |
| Mean pairwise identity | 24.4% | 25.7% |
| Parsimony-informative sites | 525 | 504 |
| Best-fit model (BIC) | LG+F+R5 | LG+F+R5 |
| Tree log-likelihood | -36569.124 | -36172.7734 |
| Mean support, clades both trees share | 95.0 | 92.5 |
| Mean support, clades only this tree has | 61.2 | 60.7 |

**Robinson-Foulds distance: 12** out of a possible 70 (17%). The trees agree on 29 of 35 clades, and 7 of 38 viruses have a different closest relative (red links in the tanglegram).

Log-likelihoods are not comparable across the two alignments (different columns), so the useful signals are the topology difference and where the support is weak.

## Clades only one tree has (strongest first)

| Tree | Support | Clade |
|---|---|---|
| ClustalW | 92 | CuPV-1, HplV-35 |
| ClustalW | 77 | AaIFV, BrBV, CuPV-1, DWV, EoV, HplV-34 (+9 more) |
| ClustalW | 71 | KV, VDV-1 |
| ClustalW | 67 | CPMV, EMCV, FMDV, HRV, HplV-81, LniV-1 (+5 more) |
| ClustalW | 32 | CuPV-1, EoV, HplV-34, HplV-35, HuAV_1, PnV (+1 more) |
| T-Coffee | 90 | AaIFV, BrBV, CuPV-1, DWV, HplV-34, HplV-35 (+7 more) |
| T-Coffee | 85 | CuPV-1, HplV-34, HplV-35, HuAV_1, IFV, SBV |
| T-Coffee | 78 | DWV, KV |
| T-Coffee | 45 | CuPV-1, HplV-34 |
| T-Coffee | 45 | EMCV, EoV, FMDV, HRV, HplV-81, LniV-1 (+5 more) |

## Viruses whose closest relatives change

| Virus | ClustalW sister group | T-Coffee sister group |
|---|---|---|
| CuPV-1 | HplV-35 | HplV-34 |
| DWV | KV, VDV-1 | KV |
| HplV-34 | CuPV-1, HplV-35 | CuPV-1 |
| HplV-35 | CuPV-1 | CuPV-1, HplV-34 |
| IFV | AaIFV, BrBV, CuPV-1, DWV (+10 more) | CuPV-1, HplV-34, HplV-35, HuAV_1 (+1 more) |
| KV | VDV-1 | DWV |
| VDV-1 | KV | DWV, KV |

![Tanglegram](comparison_tanglegram.png)

Ryan Smith – mosquito malaria single cell rnaseq, finish it up using monocle. Doesn’t care if it is monocle, any software can be used. He wants to identify progenitor cells.
He has a hunch on which cells are precursor, and he thinks that 1/8 groups are junk, as they are likely doublets.
He is sending me expression tables, metadata tables, and data for individual clusters.
How do cells from a previous scRNA seq analysis.  How does it fit with ours? Wants to compare the 24 cells from the previous study (PNAS).  Simple as here are our groupings and clusters, and want to see how these 24 cells group with their clusters. How do their 24 cells cluster with their monocle cell groups.  
They have other data that negates the two cell type group's hypothesis in PNAS. Do not reinvent the wheel, do it simple and as basic as possible, a crude analysis.  He wants to say that "these 24 cells group with these clusters", and some kind of visual to show that.

Because the data from Ryan's project is log normalized, and the data could be questionably different, we will perform 2 analyses.  1.  Monocle clustering with Ryan's data. 2. Monocle clustering of Ryan's data + 24 external samples.

We were also tasked with performing a monocle pseudotime analysis with Ryan's data.

May 05/04-07?/20
1.  Use lineage 5 as root for just the smith data.
2.  Remove all non-common genes and see how smith/severo cluster –expect those to be dispersed around, not clustering to a single spot.
3.  Do severo data by themselves
4.  Find genes with expression that define the populations.

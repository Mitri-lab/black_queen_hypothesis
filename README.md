# black_queen_hypothesis

This repository hosts figures and code for testing the black queen hypothesis.  
So far I tested if I see more deletions in co-evolved samples in PacBio data.  
This analysis was carried out using the [deletion_detection](https://github.com/nahanoo/deletion_detection).  

Double checking the identified PacBio deletions I realized that some deletions found in the Illumina data were probably due to low coverage areas rather than actual deletions.  
I suspected that this could be due tho high GC areas in the *Microbacterium saperdae* which already has a very high average GC content.  
In order to test this I developed [gc_bias](https://github.com/nahanoo/gc_bias).
This analysis indeed showed that for *Microbacterium saperdae* there are areas in the genome with no coverage correlating with high GC content in such areas.
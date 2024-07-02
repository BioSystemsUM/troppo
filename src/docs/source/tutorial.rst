Tutorials
==============================

A typical workflow follows two main steps. The first is to attribute a score to each reaction of the model, in accordance with the omics data imputed. The second is to use the scores and apply an integration method to select a subset of reactions to build the final model.

Integration scoring methods implemented in *Troppo* are:

- continuous: `ContinuousScoreIntegrationStrategy`
- threshold: `ThresholdSelectionIntegrationStrategy`
- default_core: `DefaultCoreIntegrationStrategy`
- adjusted_score: `AdjustedScoreIntegrationStrategy`
- custom: `CustomSelectionIntegrationStrategy`

Omics integration methods implemented in *Troppo* are:

- gimme: `GIMME`
- tinit: `tINIT`
- fastcore: `GIMME`
- imat: `IMAT`
- swiftcore: `SWIFTCORE`
- corda: `CORDA`

Note that the appropriate integration scoring method can differ between integration algorithms. For instance, for *GIMME* a continuous scoring method can be used, while for `fastcore` a threshold scoring method is more adequate.

Moreover, gene-level thresholding can be applied to the omics data before integration. This can be done using the `GeneLevelThresholding` class.

.. toctree::
   :maxdepth: 2

   tutorial_gimme
   tutorial_batch_run
   tutorial_task_eval
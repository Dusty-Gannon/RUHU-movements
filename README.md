# Pairing automated mark-recapture and social network models to explore the effects of landscape configuration on hummingbird foraging patterns

D. G. Gannon, A. S. Hadley, and S. J. K. Frey

## Summary

This repository documents the analysis my coauthors and I conducted to test for effects of landscape features on Rufous Hummingbird (*Selasphorus rufus*) foraging movement patterns. We aimed to test two hypotheses of how forest encroachment could reduce the funtional connectivity of the landscape by influencing hummingbird foraging patterns.

* **Barrier hypothesis**: If hummingbirds fly low, foraging at flowers in meadows where resources are abundant and avoid flying into the forest, forests may act as a barrier to hummingbird movement. If true, we would predict few movements among resources that are separated from one-another by closed-canopy forest as well as those located inside the forest.

* **Resource discovery hypothesis**: Hummingbirds are capable fliers, so they could alternatively fly over the canopy to forage in other meadows. If true, we would predict frequent movements to and from all food resources except those inside the forest since hummingbirds would be less likely to see the resources when flying over the canopy. 

We (mostly my coauthors, Drs. [Sarah J. K. Frey](http://sarahjkfrey.com/) and [Adam S. Hadley](https://www.forestbiodiversity.org/adam-hadley)) established four arrays of five hummingbird feeders on four peaks along Frizzel Ridge in the H. J. Andrews Experimental Forest. We established a center feeder in a large, central alpine meadow and four satellite feeders c.a. 250m from the center. The satellite feeders were positioned such that at least one was in the open and connected to the center feeder by open habitat, one was in the open but separated from the center by coniferous forest canopy, and one was placed under coniferous forest canopy. We (again, mostly my coauthors) implanted 163 Rufous Hummingbirds with subcutaneous transponders that could be read at recording stations fitted to the feeders. 

We utilized a class of generalized regression models developed for social network data often called "sender-receiver" or "round-robin" models (see [Warner et al. 1979](https://www.researchgate.net/publication/232567388_A_new_round_robin_analysis_of_variance_for_social_interaction_data) for a statistical description and [Fletcher et al. 2011](https://www-pnas-org.ezproxy.proxy.library.oregonstate.edu/content/108/48/19282) for an application to landscape ecology). A full model description for our model can be found in the [appendix](https://github.com/Dusty-Gannon/RUHU-movements/blob/main/Gannon_et_al_RSBL_ESM_appendixS1.pdf) of our manuscript.

## Repository contents

### [Data](https://github.com/Dusty-Gannon/RUHU-movements/tree/main/Data) folder

* `forest.tif`: A raster layer that shows forested and non-forested regions based on digitization throughout the study area.

* `hja_hummingbird_captures_rfid.csv`: Data collected on hummingbirds that were caught and implanted with Passive Integrate Transponders

* `hja_reader_locations_centercomb.txt`: Data on feeder locations and characteristics of the locations (tab-separated data).

* `total_movements_by_year.tsv`: Data on the total number of movements between each pair of feeders for each year the feeders were out on the landscape. These data are summaries of the raw data that are used to fit the sender-receiver model (see [code for fitting the model here](https://github.com/Dusty-Gannon/RUHU-movements/blob/main/Analysis/Bilinear_model_bird_movement.md)). Raw data will be made available ASAP on the [H. J. Andrews LTER data portal](http://andlter.forestry.oregonstate.edu/data/abstract.aspx?dbcode=SA028). Code to get from the raw data to this dataset can be found in the [Data prep folder](https://github.com/Dusty-Gannon/RUHU-movements/tree/main/Data_prep).

### [Data_prep](https://github.com/Dusty-Gannon/RUHU-movements/tree/main/Data_prep) folder

* `raw_reads_to_visits.Rmd`: Documents our process of removing records that all belong to the same visit to the feeder. This is the first step to reproducing our results if one were to download the raw data from the [H.J. Andrews Data Portal](https://andrewsforest.oregonstate.edu/data) (database code SA028).

* `hummer_movements.(R)md`: Documents code used to get from "visit data" (records of birds separated by at least 30s) to the movement data (`total_movements_by_year.tsv`). This is the second step to reproducting our results.

* `list_of_rfid_tests_and_glitches.R`: An R file containing a list of reader glitches and test tags that should be removed from raw data.

### [Analysis](https://github.com/Dusty-Gannon/RUHU-movements/tree/main/Analysis)

* `Bilinear_model_figure_1.Rmd`: This code can be used to reproduce a map of the feeders on the landscape with lines connecting each that are weighted by the frequency of hummingbird movements among locations.

* `Bilinear_model_bird_movement.(R)md`: This file documents the model fitting procedure using the `rstan` package (see [here](https://mc-stan.org/users/interfaces/rstan)). It also contains a brief model description, though a more detailed on the model can be found in the [appendix](https://github.com/Dusty-Gannon/RUHU-movements/blob/main/Gannon_et_al_RSBL_ESM_appendixS1.pdf) of our manuscript.

* `Competition1_bird_use_of_feeders.Rmd`: This file documents a supplementary analysis we did in response to a reviewer's question. Here, we looked deeper into how many birds were using each feeder and if there was evidence that some birds monopolized the use of certain feeders.

* `Competition2_movement_from_busy_feeders.Rmd`: This file documents a supplementary analysis we did in response to a reviewer's question. Here, we asked whether birds were more likely to move to another feeder after visiting a busy feeder (a feeder being used by many birds).

* `Competition3_birds_notin_movement_data.Rmd`: This file documents a supplementary analysis we did in response to a reviewer's question. Here, we summarized data on the birds that were not recorded moving among feeders but were recorded at feeders on at least 5 separate days. We asked whether these birds may have monopolized certain feeders, forcing others to move to other feeders.

* `summaries.(R)md`: Demographic and overall data summaries. 






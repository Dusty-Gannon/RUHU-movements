# Pairing automated mark-recapture and social network models to explore the effects of landscape configuration on hummingbird foraging patterns

D. G. Gannon, A. S. Hadley, and S. J. K. Frey

## Summary

This repository documents the analysis my coauthors and I conducted to test for a effects of landscape features on Rufous Hummingbird (*Selasphorus rufus*) foraging movement patterns. We aimed to test two hypotheses of how forest encroachment could reduce the funtional connectivity of the landscape by influencing hummingbird foraging patterns.

* **Barrier hypothesis**: If hummingbirds fly low, foraging at flowers in meadows where resources are abundant and avoid flying into the forest, forests may act as a barrier to hummingbird movement. If true, we would predict few movements to feeders that are diconnected from others by closed-canopy forest and those placed inside the forest.

* **Resource discovery hypothesis**: Hummingbirds are capable fliers, so they could alternatively fly over the canopy to forage in other meadows. If true, we would predict frequent movements to and from all feeders except those placed inside the forest since hummingbirds would be less likely to see the feeders when flying over the canopy. 

We (mostly my coauthors, Drs. [Sarah J. K. Frey](http://sarahjkfrey.com/) and [Adam S. Hadley](https://www.forestbiodiversity.org/adam-hadley)) established four arrays of five hummingbird feeders on four peaks along Frizzel Ridge in the H. J. Andrews Experimental Forest. We established a center feeder in a large, central alpine meadow and four satellite feeders c.a. 250m from the center. The satellite feeders were positioned such that at least one was in the open and connected to the center feeder by open habitat, one was in the open but separated from the center by coniferous forest canopy, and one was placed under coniferous forest canopy. We (again, mostly my coauthors) implanted 163 Rufous Hummingbirds with subcutaneous transponders that could be read at recording stations fitted to the feeders. 

We utilized a class of generalized regression models developed for social network data often called "sender-receiver" or "round-robin" models (see [Warner et al. 1979](https://www.researchgate.net/publication/232567388_A_new_round_robin_analysis_of_variance_for_social_interaction_data) for a statistical description and [Fletcher et al. 2011](https://www-pnas-org.ezproxy.proxy.library.oregonstate.edu/content/108/48/19282) for an application to landscape ecology). A full model description for our model can be found in the [appendix](https://github.com/Dusty-Gannon/RUHU-movements/blob/main/Gannon_et_al_RSBL_ESM_appendixS1.pdf) of our manuscript.

## Repository contents

### [Data](https://github.com/Dusty-Gannon/RUHU-movements/tree/main/Data) folder

* `forest.tif`: A raster layer that shows forested and non-forested regions based on digitization throughout the study area.

* `hja_hummingbird_captures_rfid.csv`: Data collected on hummingbirds that were caught and implanted with Passive Integrate Transponders

* `hja_reader_locations_centercomb.txt`: Data on feeder locations and characteristics of the locations (tab-separated data).

* `total_movements_by_year.tsv`: Data on the total number of movements between each pair of feeders for each year the feeders were out on the landscape. These data are summaries of the raw data that are used to fit the sender-receiver model (see [code for fitting the model here](https://github.com/Dusty-Gannon/RUHU-movements/blob/main/Analysis/Bilinear_model_bird_movement.md)). Raw data will be made available ASAP on the [H. J. Andrews LTER data portal](http://andlter.forestry.oregonstate.edu/data/abstract.aspx?dbcode=SA028). Code to get from the raw data to this dataset can be found in the [Data prep folder](https://github.com/Dusty-Gannon/RUHU-movements/tree/main/Data_prep).

### [Data_prep](https://github.com/Dusty-Gannon/RUHU-movements/tree/main/Data_prep) folder

* `hummer_movements.md`: Documents code used to get from "visit data" (records of birds separated by at least 30s) to the movement data (`total_movements_by_year.tsv`).

* `list_of_rfid_tests_and_glitches.R`: An R file containing a list of reader glitches and test tags that should be removed from raw data.




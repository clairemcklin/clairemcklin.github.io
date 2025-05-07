# Data-Intensive Ecology Spring 2025

## Abstract
<p>Plant-microbe relationships play an important role in maintaining plant health and ecosystem resilience, but global climate change threatens to disrupt these vital ecosystem interactions. This study examines how soil microbial biodiversity affects plant biodiversity using R to visualize data from the National Ecological Observatory Network (NEON). Based on existing literature, I predict that plant biodiversity will increase as soil microbe biodiversity increases. I will measure biodiversity using the Shannon Diversity Index and species richness values. The Soil Microbe Community Composition, Plant Presence and Percent Cover, and Soil Physical and Chemical Properties datasets from NEON will be used to investigate the relationship between plant and microbe biodiversity. Soil temperature from nineteen sites across different climate zones in the United States will be compared to the plant and microbe biodiversity data to show differences in diversity and species richness among different climates and habitat types. Plant and microbe species reliant on interspecies symbiotic relationships may be under harm by rising global temperatures; differences in climate niches between interdependent species may cause increased rates of extinction as global temperatures increase. However, the relationship between plants and soil microbes may also be harnessed to restore ecosystems degraded by natural or anthropogenic disasters through sustainable agricultural practices.</p>

## Graphs

### Comparing Total Microbe and Plant Species Richness
<img src="/rich_scatter.png" alt="Relationship Between Plant and Microbe Richness" title="Relationship Between Plant and Microbe Richness" width="470"/>
<img src="/Rplot47.png" alt="Relationship Between Plant and Microbe Richness" title="Relationship Between Plant and Microbe Richness" width="470"/>


### Comparing Total Microbe and Plant Shannon Diversity Indices by Soil Properties
<img src="/boxplot_temp.png" alt="Shannon Diversity Index by Soil Temperature" title="Shannon Diversity Index by Soil Temperature" width="470"/>
<img src="/boxplot_moisture.png" alt="Shannon Diversity Index by Soil Moisture" title="Shannon Diversity Index by Soil Moisture" width="470"/>

### Comparing Fungi and Plant Species Richness


## Code

[Link To Biodiversity Code Used in this Project](https://github.com/clairemcklin/clairemcklin.github.io/blob/b0df39e89076a31a77f2a6da017f4422db911f20/biodiversity.R)


## Conclusion

<ul>
  <li>There is a positive correlation between plant and soil microbe biodiversity that is partially caused by soil moisture and temperature.</li>
  <li>Lower soil temperatures and higher moisture contents generally showed a higher correlation between plant and microbe/fungi biodiversity.</li>
  <li>Future directions: analyze correlations between plant and microbe biodiversity in sites before and after ecosystem restoration or degradation</li>
</ul>

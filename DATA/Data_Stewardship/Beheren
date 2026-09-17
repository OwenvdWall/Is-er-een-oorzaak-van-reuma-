Beheer van onderzoeksgegevens

Goed beheer van onderzoeksgegevens (Research Data Management) waarborgt dat
onderzoek overzichtelijk, controleerbaar en reproduceerbaar blijft. Dit project
is gestructureerd opgeslagen in een openbare GitHub-repository om de transparantie
en het versiebeheer te optimaliseren.

Mappenstructuur en bestandsnaamgeving
De repository is logisch onderverdeeld om data, scripts en resultaten
strikt gescheiden te houden:
- DATA: Opgesplitst in de submappen RUWE DATA en GEFILTERDE DATA om de
  originele bronbestanden intact te laten.
- RESULTATEN: Bevat alle gegenereerde visualisaties met herkenbare en
  beschrijvende namen (zoals GO_Enrichment.png en PCA&Vulcanoplot.png).
- Scripts en Documentatie: Het R-script en de centrale README.md staan
  direct in de root voor snelle toegang.

Reproduceerbaarheid en documentatie
Om de analyse volledig reproduceerbaar te maken, documenteert de README.md
de exacte softwareomgeving en statistische drempelwaarden:
- Softwareomgeving: Uitgevoerd in R 4.5.2 met specifieke package-versies voor
  dataverwerking (Rsubread 2.24.0, Rsamtools 2.26.0, DESeq2 1.50.2) en
  visualisatie (ggplot2 4.0.3, EnhancedVolcano 1.28.2).
- Analysecriteria: De pijplijn loopt van RNA-seq-verwerking tot differentiële
  genexpressie en functionele verrijking (GO, KEGG en Pathview). Genenselectie
  vond plaats op basis van een adjusted p-value < 0,05 en een absolute
  log2 fold change > 1.

Versiebeheer en privacy
- Versiebeheer: Via GitHub worden alle codewijzigingen nauwkeurig bijgehouden.
  Dit maakt het mogelijk om fouten snel te traceren en zo nodig terug te
  grijpen naar eerdere versies.
- Privacy en Data-ethiek: Om de privacy van patiënten te beschermen, bevat
  deze repository geen direct herleidbare medische gegevens. Er is uitsluitend
  gewerkt met een reeds gepubliceerde, openbare en geanonymiseerde dataset van
  Platzer et al. (2019), bestaande uit synoviumbiopten van vier RA-patiënten
  en vier controles.

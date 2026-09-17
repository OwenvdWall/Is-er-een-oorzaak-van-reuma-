Tijdens dit transcriptomics-project is GitHub gebruikt als
platform voor het beheren van de onderzoeksgegevens, de
R-scripts en de bijbehorende documentatie. Een gestructureerd
databeheer zorgt ervoor dat gegevens overzichtelijk blijven,
analyses reproduceerbaar zijn en keuzes traceerbaar blijven.

Mappenstructuur en bestandsbeheer
Binnen de repository is een vaste mappenstructuur gebruikt om
alle onderdelen overzichtelijk te organiseren:
- DATA: Bevat de submappen RUWE DATA (de oorspronkelijke dataset
  van Platzer et al., 2019) en GEFILTERDE DATA (de bewerkte
  datasets uit R). Dit houdt bronbestanden strikt gescheiden.
- R-SCRIPT: Staat in de rootmap en bevat de workflow voor read
  mapping (Rsubread), genkwantificatie (featureCounts),
  differentiële expressie (DESeq2) en pathway-analyses.
- RESULTATEN: Bevat alle visualisaties met herkenbare namen,
  zoals GO_Enrichment.png en hsa04660.pathview.png.

Documentatie en reproduceerbaarheid
In het README.md-bestand is de achtergrond van het onderzoek naar
reumatoïde artritis beschreven, samen met de analysemethode,
resultaten en conclusie. Hierin staan ook de exacte softwareversies
(zoals R 4.5.2) en de criteria (adjusted p-value < 0,05 en absolute
log2 fold change > 1) vastgelegd. Dit zorgt ervoor dat een andere
onderzoeker de resultaten exact kan reproduceren.

Versiebeheer en privacy
- Versiebeheer: Via commits zijn alle wijzigingen in de code en
  documentatie nauwkeurig bijgehouden. Dit maakt het mogelijk om
  bij fouten snel terug te grijpen naar eerdere versies.
- Privacy: De repository bevat uitsluitend de openbare, geanonymiseerde
  data van 8 specifieke sample ID's (SRR4785819 t/m SRR4785988),
  waardoor patiëntgegevens volledig beschermd blijven.

Door deze inrichting blijft de volledige RNA-seq workflow transparant,
veilig en direct bruikbaar voor derden.

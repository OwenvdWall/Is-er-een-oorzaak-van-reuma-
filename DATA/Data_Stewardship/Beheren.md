Tijdens dit transcriptomics-project is GitHub gebruikt als
platform voor het beheren van de onderzoeksgegevens, de
R-scripts en de bijbehorende documentatie. Een gestructureerd
databeheer zorgt ervoor dat gegevens overzichtelijk blijven,
analyses reproduceerbaar zijn en keuzes traceerbaar blijven.

Binnen de repository is een vaste mappenstructuur gebruikt om
alle onderdelen overzichtelijk te organiseren. De hoofdmap
DATA bevat de submappen RUWE DATA voor de oorspronkelijke dataset
van Platzer et al. (2019) en GEFILTERDE DATA voor de bewerkte
datasets uit R, zodat bronbestanden strikt gescheiden blijven.
Het bestand R-SCRIPT staat direct in de rootmap en bevat de
volledige workflow voor read mapping met Rsubread, genkwantificatie
via featureCounts, differentiële expressie met DESeq2 en de
pathway-analyses. Alle gegenereerde visualisaties zijn opgeslagen
in de map RESULTATEN en hebben herkenbare namen gekregen, zoals
GO_Enrichment.png en hsa04660.pathview.png.

In het README.md-bestand is de achtergrond van het onderzoek naar
reumatoïde artritis beschreven, samen met de analysemethode,
de resultaten en de uiteindelijke conclusie. Hierin staan ook
de exacte softwareversies zoals R 4.5.2 en de statistische
criteria van een adjusted p-value < 0,05 en een absolute
log2 fold change > 1 vastgelegd. Dit zorgt ervoor dat een andere
onderzoeker de resultaten in de toekomst exact kan reproduceren.

Het versiebeheer is ingericht via commits, waarmee alle
wijzigingen in de code en documentatie nauwkeurig zijn bijgehouden.
Dit maakt het mogelijk om bij eventuele fouten in het script
snel terug te grijpen naar eerdere, werkende versies. Wat betreft
de privacy bevat de repository uitsluitend de openbare, volledig
geanonymiseerde data van de acht specifieke sample ID's van
SRR4785819 tot en met SRR4785988, waardoor medische patiëntgegevens
optimaal beschermd blijven. Door deze zorgvuldige inrichting
blijft de volledige RNA-seq workflow transparant, veilig en
direct bruikbaar voor wie het onderzoek wil bekijken.

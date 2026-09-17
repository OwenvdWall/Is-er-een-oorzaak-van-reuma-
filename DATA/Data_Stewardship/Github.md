
In dit project werk ik als data steward en beheer ik de gegevens
via een openbare GitHub-repository. De pagina is zo ingericht dat
alles overzichtelijk, veilig en makkelijk te volgen is voor wie
het onderzoek wil bekijken.

Binnen de repository heb ik een vaste mappenstructuur gebruikt.
De hoofdmap DATA heeft twee submappen: RUWE DATA voor de originele
dataset van Platzer et al. (2019) en GEFILTERDE DATA voor de
bestanden uit R. Zo blijven de bronbestanden netjes apart staan.
Het bestand R-SCRIPT staat direct in de hoofdmap en bevat de hele
code voor het uitlijnen (Rsubread), tellen (featureCounts), de
expressieanalyse (DESeq2) en de pathways. Alle plaatjes staan in de
map RESULTATEN en hebben logische namen zoals GO_Enrichment.png en
hsa04660.pathview.png. Hierdoor kan een ander alles snel vinden.

De pagina helpt ook om het onderzoek reproduceerbaar te maken.
In het README.md-bestand staat de inleiding over reuma, de methode,
de resultaten en de conclusie. Ook worden daar de exacte versies
genoemd (zoals R 4.5.2) en de instellingen die ik heb gebruikt
(adjusted p-value < 0,05 en log2 fold change > 1). Zo kan een ander
onderzoek precies zo worden nagedaan in de toekomst.

Het versiebeheer doe ik via commits, waarmee alle aanpassingen
in de code en tekst netjes worden opgeslagen. Als er een foutje in
mijn script sluipt, kan ik snel terug naar een eerdere versie.
Voor de privacy staan er geen persoonlijke patiëntgegevens op de
pagina. Ik gebruik alleen de openbare, geanonymiseerde data van
de 8 sample ID's (SRR4785819 t/m SRR4785988). De data is dus veilig
en het project voldoet aan de privacyregels.

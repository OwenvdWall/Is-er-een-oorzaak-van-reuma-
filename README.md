# Van RNA-seq tot Pathview: brede opregulatie van de T-Cel receptor cascade bij Reumatoïde Artritis
## Inleiding

Reumatoïde artritis (RA) is een chronische auto-immuunziekte die wordt gekenmerkt door ontstekingen van de gewrichten. Deze ontstekingsreactie leidt tot schade aan kraakbeen en bot en kan zich uitbreiden naar omliggende weefsels. RA ontstaat door een interactie tussen genetische factoren en omgevingsfactoren, zoals roken (Gabriel, 2001). De exacte oorzaak is nog niet volledig bekend, maar onderzoek toont aan dat zowel genetische aanleg als ontregeling van het immuunsysteem een belangrijke rol spelen.

RA treft voornamelijk synoviale gewrichten, die bestaan uit een gewrichtskapsel, synoviaal membraan en synoviaal vocht. Klinisch presenteert de aandoening zich meestal symmetrisch, waarbij de ontsteking vaak begint in de kleine gewrichten van handen en voeten en zich vervolgens uitbreidt naar grotere gewrichten (Radu & Bungau, 2021).

Om dieper inzicht te verkrijgen in de onderliggende mechanismen van reumatoïde artritis is in dit project een transcriptomicsanalyse uitgevoerd op basis van RNA-seq-data. Het doel van deze analyse is het identificeren van verschillen in genexpressie tussen gezonde controles en RA-patiënten. Hiervoor is gebruikgemaakt van genexpressiedata uit verschillende studies, die diverse klinische stadia van RA omvatten, waaronder zowel vroege als latere fasen van de ziekte (Platzer et al., 2019). Door vergelijking met gezonde monsters kunnen specifieke genen en biologische pathways worden geïdentificeerd die bijdragen aan de pathogenese van RA. 



## Methode

Voor dit onderzoek is gebruikgemaakt van een gepubliceerde RNA-seq-dataset van Platzer et al. (2019), afkomstig uit de NCBI Sequence Read Archive (SRA). De dataset bestaat uit synoviumbiopten van vier patiënten met reumatoïde artritis (RA) en vier gezonde controles. De klinische kenmerken van de gebruikte samples zijn weergegeven in Tabel 1. De toegepaste bio-informatische workflow is schematisch weergegeven in Figuur 1.

##### Tabel 1: Overzicht van de klinische kenmerken van de gebruikte RNA-seq dataset (Platzer et al., 2019).
| Sample ID | Leeftijd | Geslacht | Conditie |
| :--- | :--- | :--- | :--- |
| **SRR4785819** | 31 | Vrouw | Controle |
| **SRR4785820** | 15 | Vrouw | Controle |
| **SRR4785828** | 31 | Vrouw | Controle |
| **SRR4785831** | 42 | Vrouw | Controle |
| **SRR4785979** | 54 | Vrouw | Reumatoïde artritis |
| **SRR4785980** | 66 | Vrouw | Reumatoïde artritis |
| **SRR4785986** | 60 | Vrouw | Reumatoïde artritis |
| **SRR4785988** | 59 | Vrouw | Reumatoïde artritis |
<br>

<p align="center">
  <img src="https://github.com/user-attachments/assets/1b22e28e-56ac-415a-a54d-bc6d37b3a9b0" width="450" style="image-rendering: -webkit-optimize-contrast;" alt="RNA-seq-analyse stroomschema">
  <br>
  <em>Figuur 1. RNA-seq-analyse stroomschema. Workflow toegepast voor de verwerking van de rauwe RNA-seq-data tot aan de pathway-visualisatie.</em>
</p>

De RNA-seq-data werden geanalyseerd in R (versie 4.5.2). Reads werden uitgelijnd op het humane referentiegenoom met Rsubread (versie 2.24.0; Liao et al., 2019), waarna de BAM-bestanden werden verwerkt met Rsamtools (versie 2.26.0; Morgan et al., 2024). Met featureCounts (Liao et al., 2014) werd een count matrix gegenereerd. Differentieel geëxpresseerde genen werden geïdentificeerd met DESeq2 (versie 1.50.2; Love et al., 2014). Visualisaties werden gemaakt met ggplot2 (versie 4.0.3; Wickham, 2016) en EnhancedVolcano (versie 1.28.2; Blighe et al., 2024).

Voor de GO- en KEGG-enrichmentanalyse werden genen geselecteerd op basis van een adjusted p-value < 0,05 en een absolute log2 fold change > 1. Hierbij werd de Benjamini-Hochberg-correctie toegepast om te corrigeren voor multiple testing. De geselecteerde genen werden geannoteerd met org.Hs.eg.db (versie 3.22.0; Carlson, 2024) en geanalyseerd met clusterProfiler (versie 4.18.4; Wu et al., 2021). De meest relevante pathway werd vervolgens gevisualiseerd met pathview (versie 1.50.0; Luo & Brouwer, 2013).




## Resultaten

### PCA plot & Volcano plot

Om de algehele datastructuur te verkennen en individuele verschillen in genexpressie te identificeren, zijn er visualisaties gemaakt. De PCA-plot (Figuur 2A) toont een scherpe scheiding tussen controles en reuma-samples. Omdat de eerste twee hoofdcomponenten samen 85% van de variatie verklaren (PC1: 74%, PC2: 11%), is direct zichtbaar dat de aanwezigheid van de ziekte de oorzaak is van de biologische verschillen. De Volcano plot (Figuur 2B) brengt dit op gen niveau in kaart en toont de actievere (rechts) en minder actieve genen (links). 

<p align="center">
  <img src="https://github.com/user-attachments/assets/a82d384c-c4dd-4fc5-8db4-aeec02d1c033" width="100%" style="image-rendering: -webkit-optimize-contrast;" alt="PCA en Volcano Plot Resultaten">
  <br>
 
  <em>Figuur 2. Resultaten van de genexpressie-analyse. (A) De PCA-plot laat een duidelijke scheiding zien tussen de gezonde controles (normal) en de reuma-samples (Rheuma). De eerste as (PC1) verklaart hierbij 74% van de verschillen. (B) De Volcano plot laat zien welke genen significant verschillen, uitgezet tegen de log2 fold change. De belangrijkste genen voor reuma zijn apart gelabeld met lijntjes, zodat de namen goed leesbaar zijn.</em>
</p>





### GO en KEGG Analyse

De Gene Ontology (GO) enrichment-analyse (Figuur 3) toont aan dat de top 10 biologische processen gedomineerd wordt door lymfocyten- en leukocytenactivatie, wat de immuuninfiltratie in het synoviale weefsel weerspiegelt. De KEGG-pathwayanalyse laat een significante verrijking zien van signaalroutes zoals cytokine-receptorinteracties (Figuur 4). De     T- cell receptor signaling pathway (hsa04660) is specifiek geselecteerd vanwege de sleutelrol van T-cellen bij het verloop van RA.

<p align="center">
  <img src="https://github.com/user-attachments/assets/1f826903-223a-4b6d-bbe2-d7ae246f3045" width="550" style="image-rendering: -webkit-optimize-contrast;" alt="Gene Ontology Enrichment Analyse">
  <br>
  <em>Figuur 3. Gene Ontology Enrichment. Dotplot van de top 10 significante biologische processen (GO), gerangschikt op GeneRatio.</em>
</p>

<p align="center">
  <img src="https://github.com/user-attachments/assets/c92b6b0a-9073-4357-bbaf-6631bf6cc375" width="550" style="image-rendering: -webkit-optimize-contrast;" alt="KEGG Pathway Enrichment Analyse">
  <br>
  <em>Figuur 4. KEGG Pathway Enrichment. Dotplot van de meest verrijkte moleculaire signaalroutes binnen de dataset.</em>
</p>

### PATHVIEW 

De Pathview-visualisatie (Figuur 5) onthult een breedschalige opregulatie (rode boxen) binnen deze cascade, zoals de MAP-kinasecascade. Dit sluit aan bij de literatuur (Radu & Bungau, 2021), waarin staat dat aanhoudende immuun signalering de pro-inflammatoire status in stand houdt.

<p align="center">
  <img src="https://github.com/user-attachments/assets/10ec3509-f01f-4929-b23d-6d7bc684da68" width="650" style="image-rendering: -webkit-optimize-contrast;" alt="Pathview T Cell Receptor Signaling Pathway">
  <br>
  <em>Figuur 5. Pathview-visualisatie van de 'T cell receptor signaling pathway' (hsa04660). Rode boxen induceren geüpreguleerde genen in de RA-monsters; groene boxen geven gedownreguleerde genen aan.</em>
</p>

## Conclusie 

De transcriptomicsanalyse bevestigt ingrijpende herprogrammering in het synoviale weefsel bij RA. De PCA-scheiding, immuun-gerelateerde GO/KEGG-verrijking en opregulatie van de T-celreceptorcascade bewijzen in lijn met Radu & Bungau (2021) dat autoreactieve T-cellen de drijvende kracht vormen achter de gewrichtsontsteking. Dit biedt waardevolle mechanistische aanknopingspunten voor gerichte immunotherapieën bij RA-patiënten.

## AI-verklaring

 Dit verslag is door de auteur zelf geschreven. AI is gebruikt om de figuren biologisch te interpreteren, waarna de tekst zelfstandig is geformuleerd en taalkundige controle heeft plaatsgevonden.


## Bronnen

Platzer, A., Nussbaumer, T., Karonitsch, T., Smolen, J. S., & Aletaha, D. (2019). Analysis of gene expression in rheumatoid arthritis and related conditions offers insights into sex-bias, gene biotypes and co-expression patterns. PLoS ONE, 14(7), e0219698. https://doi.org/10.1371/journal.pone.0219698

Radu, A., & Bungau, S. G. (2021). Management of Rheumatoid Arthritis: An Overview. Cells, 10(11), 2857. https://doi.org/10.3390/cells10112857

## Packages

 Rsubread & featureCounts:
  Liao, Y., Smyth, G. K., & Shi, W. (2019). The R package Rsubread is 
  easier, faster, cheaper and better for alignment and quantification 
  of RNA sequencing reads. Nucleic Acids Research, 47(8), e47.

DESeq2:
  Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of 
  fold change and dispersion for RNA-seq data with DESeq2. 
  Genome Biology, 15(12), 550.

clusterProfiler:
  Wu, T., Hu, E., Xu, S., Chen, M., Guo, P., Dai, Z., ... & Yu, G. 
  (2021). clusterProfiler 4.0: A universal enrichment tool for 
  interpreting omics data. The Innovation, 2(3), 100141.

EnhancedVolcano:
  Blighe, K., Rana, S., & Lewis, M. (2020). EnhancedVolcano: 
  Publication-ready volcano plots with enhanced colouring and labeling. 
  R package version 1.7.14. https://github.com/kevinblighe/EnhancedVolcano

Pathview:
  Luo, W., & Brouwer, C. (2013). Pathview: an R/Bioconductor package 
  for pathway-based data integration and visualization. 
  Bioinformatics, 29(14), 1830-1831.

ggplot2:
  Wickham, H. (2016). ggplot2: Elegant Graphics for Data Analysis. 
  Springer-Verlag New York.


Van RNA-seq tot Pathview: Brede Opregulatie van de T-Cel Receptor Cascade bij Reumatoïde Artritis
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

in totaal werden er 29.407 genen geanalyseerd. de PCA analyse laat een heldere scheiding zien tussen de gezonde groep. de eerste twee hoofd de eerste twee hoofdcomponenten verklaren namelijk 74% en 11% van alle variatie. Zie (Figuur 2A) uit verder onderzoek bleek dat 4572 genen een significant verschil in expressie hadden (padj <= 0.05 en |log2 fold change| >= 1). De Volcano plot in figuur 2B brengt deze genen in beeld. hierop is te zien dat genen opgereguleerd zijn of juis neergereguleerd zijn.

<p align="center">
  <img src="https://github.com/user-attachments/assets/a82d384c-c4dd-4fc5-8db4-aeec02d1c033" width="100%" style="image-rendering: -webkit-optimize-contrast;" alt="PCA en Volcano Plot Resultaten">
  <br>
 
  <em>Figuur 2. Resultaten van de genexpressie-analyse. (A) De PCA-plot laat een duidelijke scheiding zien tussen de gezonde controles (normal) en de reuma-samples (Rheuma). De eerste as (PC1) verklaart hierbij 74% van de verschillen. (B) De Volcano plot laat zien welke genen significant verschillen, uitgezet tegen de log2 fold change. De belangrijkste genen voor reuma zijn apart gelabeld met lijntjes, zodat de namen goed leesbaar zijn.</em>
</p>





### GO ANALYSE 
<p align="center">
  <img src="https://github.com/user-attachments/assets/1f826903-223a-4b6d-bbe2-d7ae246f3045" width="550" style="image-rendering: -webkit-optimize-contrast;" alt="Gene Ontology Enrichment Analyse">
  <br>
  <em>Figuur 3. Gene Ontology Enrichment. Dotplot van de top 10 significante biologische processen (GO), gerangschikt op GeneRatio.</em>
</p>



### KEGG ANALYSE
<p align="center">
  <img src="https://github.com/user-attachments/assets/c92b6b0a-9073-4357-bbaf-6631bf6cc375" width="550" style="image-rendering: -webkit-optimize-contrast;" alt="KEGG Pathway Enrichment Analyse">
  <br>
  <em>Figuur 4. KEGG Pathway Enrichment. Dotplot van de meest verrijkte moleculaire signaalroutes binnen de dataset.</em>
</p>




### PATHVIEW 
<p align="center">
  <img src="https://github.com/user-attachments/assets/10ec3509-f01f-4929-b23d-6d7bc684da68" width="650" style="image-rendering: -webkit-optimize-contrast;" alt="Pathview T Cell Receptor Signaling Pathway">
  <br>
  <em>Figuur 5. Pathview-visualisatie van de 'T cell receptor signaling pathway' (hsa04660). Rode boxen induceren geüpreguleerde genen in de RA-monsters; groene boxen geven gedownreguleerde genen aan.</em>
</p>





## Conclusie 


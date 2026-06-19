# Is er een oorzaak voor reuma? Tijdelijke titel
## Inleiding

Reumatoïde artritis (RA) is een chronische auto-immuunziekte die wordt gekenmerkt door ontstekingen van de gewrichten. Deze ontstekingsreactie leidt tot schade aan kraakbeen en bot en kan zich uitbreiden naar omliggende weefsels. RA ontstaat door een interactie tussen genetische factoren en omgevingsfactoren, zoals roken (Gabriel, 2001). De exacte oorzaak is nog niet volledig bekend, maar onderzoek toont aan dat zowel genetische aanleg als ontregeling van het immuunsysteem een belangrijke rol spelen.

RA treft voornamelijk synoviale gewrichten, die bestaan uit een gewrichtskapsel, synoviaal membraan en synoviaal vocht. Klinisch presenteert de aandoening zich meestal symmetrisch, waarbij de ontsteking vaak begint in de kleine gewrichten van handen en voeten en zich vervolgens uitbreidt naar grotere gewrichten (Radu & Bungau, 2021).

Om dieper inzicht te verkrijgen in de onderliggende mechanismen van reumatoïde artritis is in dit project een transcriptomicsanalyse uitgevoerd op basis van RNA-seq-data. Het doel van deze analyse is het identificeren van verschillen in genexpressie tussen gezonde controles en RA-patiënten. Hiervoor is gebruikgemaakt van genexpressiedata uit verschillende studies, die diverse klinische stadia van RA omvatten, waaronder zowel vroege als latere fasen van de ziekte (Platzer et al., 2019). Door vergelijking met gezonde monsters kunnen specifieke genen en biologische pathways worden geïdentificeerd die bijdragen aan de pathogenese van RA. 



## Methode
Voor dit onderzoek is gebruikgemaakt van een gepubliceerde RNA-seq dataset van Platzer et al. (2019)bestaande uit synoviumbiopten van vier RA-patiënten  en vier gezonde controle. De RA-groep betreft patiënten met een gevestigde diagnose (>12 maanden0 en ACPA-positiviteit, terwijl de controlegroep ACPA-negatief is. De ruwe data zijn beschikbaar via NCBI Sequence Read Archive (SRA) een geïdentificeerd met SRR-accessienummers. Deze dataset maakt een vergelijking mogelijk tussen gezonde en aangedane weefsels gericht op het identificeren van differentieel tot expressie komende genen en betrokken biologische pathways bij RA. Zie (Figuur 1)


<p align="center">
  src=<img width="342" height="2181" alt="image" src="https://github.com/user-attachments/assets/1484c261-433c-47e8-aa98-1f17e10517e3" />
 <br>
</p>
<p align="left">
  <p>
<em><b>Figuur 1.</b> RNA-seq-analyse stroomschema. Overzicht van de bio-informatische workflow die is toegepast voor de analyse van RNA-seq-data van gezonde controles en patiënten met reumatoïde artritis. De workflow start met RNA-seq-data en het humane referentiegenoom, gevolgd door genome indexing, read alignment, BAM-verwerking en het opstellen van een count matrix. Vervolgens werd met behulp van metadata een differentiële genexpressieanalyse uitgevoerd met DESeq2. De resultaten werden gevisualiseerd met een Volcano plot en een PCA-plot. Op basis van de significante differentieel tot expressie komende genen werden een Gene Ontology (GO)-analyse en een KEGG-pathwayanalyse uitgevoerd. De geselecteerde pathway werd ten slotte gevisualiseerd met Pathview. Het stroomschema is gemaakt met Microsoft Visio.</em>
</p>

De analyse werd uitgevoerd in R (versie 4.5.2). BiocManager (1.30.27) werd gebruikt voor het installeren van Bioconductor-packages. Het referentiegenoom werd geïndexeerd en de RNA-seq-data werd uitgelijnd met Rsubread (2.24.0), waarna de BAM-bestanden werden gesorteerd en geïndexeerd met Rsamtools (2.26.0). Met featureCounts (Rsubread) werd een count matrix opgesteld. De differentiële genexpressieanalyse werd uitgevoerd met DESeq2 (1.50.2). De resultaten werden gevisualiseerd met EnhancedVolcano (1.28.2), ggplot2 (4.0.3) en enrichplot (1.30.5). Voor de Gene Ontology- en KEGG-analyse werd clusterProfiler (4.18.4) gebruikt, in combinatie met org.Hs.eg.db (3.22.0) voor genannotatie. De geselecteerde KEGG-pathway werd gevisualiseerd met pathview (1.50.0). dplyr (1.2.0) werd gebruikt voor het filteren en bewerken van de resultaten.


## Resultaten

### PCA plot & Volcano plot

in totaal werden er 29.407 genen geanalyseerd. de PCA analyse laat een heldere scheiding zien tussen de gezonde groep. de eerste twee hoofd de eerste twee hoofdcomponenten verklaren namelijk 74% en 11% van alle variatie. Zie (Figuur 2A) uit verder onderzoek bleek dat 4572 genen een significant verschil in expressie hadden (padj <= 0.05 en |log2 fold change| >= 1). De Volcano plot in figuur 2B brengt deze genen in beeld. hierop is te zien dat genen opgereguleerd zijn of juis neergereguleerd zijn.

<p align="center">
  <img width="342" height="2181" alt="image" src="<img width="1280" height="720" alt="image" src=<img width="1280" height="720" alt="image" src="https://github.com/user-attachments/assets/1fe016d0-b6e3-4e82-a2b8-bfbf910dbd4c" />
 />
" />
 <br>
</p>
<p align="left">
  <p>
<em><b>Figuur 2.</b>





### GO ANALYSE 


### KEGG ANALYSE



### PATHVIEW 


## Conclusie 


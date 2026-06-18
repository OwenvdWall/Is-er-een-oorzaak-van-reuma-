# Is er een oorzaak voor reuma?
## Inleiding

Reumatoïde artritis (RA) is een chronische auto-immuunziekte die wordt gekenmerkt door ontstekingen van de gewrichten. Deze ontstekingsreactie leidt tot schade aan kraakbeen en bot en kan zich uitbreiden naar omliggende weefsels. RA ontstaat door een interactie tussen genetische factoren en omgevingsfactoren, zoals roken (Gabriel, 2001). De exacte oorzaak is nog niet volledig bekend, maar onderzoek toont aan dat zowel genetische aanleg als ontregeling van het immuunsysteem een belangrijke rol spelen.

RA treft voornamelijk synoviale gewrichten, die bestaan uit een gewrichtskapsel, synoviaal membraan en synoviaal vocht. Klinisch presenteert de aandoening zich meestal symmetrisch, waarbij de ontsteking vaak begint in de kleine gewrichten van handen en voeten en zich vervolgens uitbreidt naar grotere gewrichten (Radu & Bungau, 2021).

Om dieper inzicht te verkrijgen in de onderliggende mechanismen van reumatoïde artritis is in dit project een transcriptomicsanalyse uitgevoerd op basis van RNA-seq-data. Het doel van deze analyse is het identificeren van verschillen in genexpressie tussen gezonde controles en RA-patiënten. Hiervoor is gebruikgemaakt van genexpressiedata uit verschillende studies, die diverse klinische stadia van RA omvatten, waaronder zowel vroege als latere fasen van de ziekte (Platzer et al., 2019). Door vergelijking met gezonde monsters kunnen specifieke genen en biologische pathways worden geïdentificeerd die bijdragen aan de pathogenese van RA.



## Methode
Voor dit onderzoek is gebruikgemaakt van een gepubliceerde RNA-seq dataset van Platzer et al. (2019)bestaande uit synoviumbiopten van vier RA-patiënten  en vier gezonde controle. De RA-groep betreft patiënten met een gevestigde diagnose (>12 maanden0 en ACPA-positiviteit, terwijl de controlegroep ACPA-negatief is. De ruwe data zijn beschikbaar via NCBI Sequence Read Archive (SRA) een geïdentificeerd met SRR-accessienummers. Deze dataset maakt een vergelijking mogelijk tussen gezonde en aangedane weefsels gericht op het identificeren van differentieel tot expressie komende genen en betrokken biologische pathways bij RA.

<p align="center">
  <img width="700" height="700" alt="Image" src="https://github.com/user-attachments/files/29107849/Methode_Reuma.pdf" />
  <br>
</p>
<p align="left">
  <em><b>Figuur 1. </b>RNA-seq analysepipeline. Overzicht van de bioinformatische workflow voor de verwerking van RNA-seq data, beginnend bij FASTQ-reads en referentiegenoom, gevolgd door read alignment, BAM-verwerking en genkwantificatie. De resulterende count matrix werd gebruikt voor differentiële expressieanalyse met DESeq2. Significant tot expressie komende genen werden vervolgens geanalyseerd via GO- en KEGG-enrichment en gevisualiseerd met Pathview. Gemaakt met Microsoft Visio</em>
</p>



## Resultaten

### Vulcanoplot



### PCA plot



### GO ANALYSE 


### KEGG ANALYSE



### PATHVIEW (T-CELL REceptor signaling)


## Conclusie 


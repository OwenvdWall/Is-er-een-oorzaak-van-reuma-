# Is er een oorzaak voor reuma?
## Inleiding

Reumatoïde artritis (RA) is een chronische auto-immuunziekte die wordt gekenmerkt door ontstekingen van de gewrichten. Deze ontstekingsreactie leidt tot schade aan kraakbeen en bot en kan zich uitbreiden naar omliggende weefsels. RA ontstaat door een interactie tussen genetische factoren en omgevingsfactoren, zoals roken (Gabriel, 2001). De exacte oorzaak is nog niet volledig bekend, maar onderzoek toont aan dat zowel genetische aanleg als ontregeling van het immuunsysteem een belangrijke rol spelen.

RA treft voornamelijk synoviale gewrichten, die bestaan uit een gewrichtskapsel, synoviaal membraan en synoviaal vocht. Klinisch presenteert de aandoening zich meestal symmetrisch, waarbij de ontsteking vaak begint in de kleine gewrichten van handen en voeten en zich vervolgens uitbreidt naar grotere gewrichten (Radu & Bungau, 2021).

Om dieper inzicht te verkrijgen in de onderliggende mechanismen van reumatoïde artritis is in dit project een transcriptomicsanalyse uitgevoerd op basis van RNA-seq-data. Het doel van deze analyse is het identificeren van verschillen in genexpressie tussen gezonde controles en RA-patiënten. Hiervoor is gebruikgemaakt van genexpressiedata uit verschillende studies, die diverse klinische stadia van RA omvatten, waaronder zowel vroege als latere fasen van de ziekte (Platzer et al., 2019). Door vergelijking met gezonde monsters kunnen specifieke genen en biologische pathways worden geïdentificeerd die bijdragen aan de pathogenese van RA.



## Methode
Voor dit onderzoek is gebruikgemaakt van een gepubliceerde RNA-seq dataset van Platzer et al. (2019)bestaande uit synoviumbiopten van vier RA-patiënten  en vier gezonde controle. De RA-groep betreft patiënten met een gevestigde diagnose (>12 maanden0 en ACPA-positiviteit, terwijl de controlegroep ACPA-negatief is. De ruwe data zijn beschikbaar via NCBI Sequence Read Archive (SRA) een geïdentificeerd met SRR-accessienummers. Deze dataset maakt een vergelijking mogelijk tussen gezonde en aangedane weefsels gericht op het identificeren van differentieel tot expressie komende genen en betrokken biologische pathways bij RA.


<p align="center">
  <img width="342" height="2181" alt="image" src="https://github.com/user-attachments/assets/ae96a182-2b78-478b-8a95-8cc8089f7384" />
 <br>
</p>
<p align="left">
  <p>
<em><b>Figuur 1.</b> RNA-seq-analyse stroomschema. Overzicht van de bio-informatische workflow die is toegepast voor de analyse van RNA-seq-data van gezonde controles en patiënten met reumatoïde artritis. De workflow start met RNA-seq-data en het humane referentiegenoom, gevolgd door genome indexing, read alignment, BAM-verwerking en het opstellen van een count matrix. Vervolgens werd met behulp van metadata een differentiële genexpressieanalyse uitgevoerd met DESeq2. De resultaten werden gevisualiseerd met een Volcano plot en een PCA-plot. Op basis van de significante differentieel tot expressie komende genen werden een Gene Ontology (GO)-analyse en een KEGG-pathwayanalyse uitgevoerd. De geselecteerde pathway werd ten slotte gevisualiseerd met Pathview. Het stroomschema is gemaakt met Microsoft Visio.</em>
</p>



## Resultaten

### Vulcanoplot



### PCA plot



### GO ANALYSE 


### KEGG ANALYSE



### PATHVIEW (T-CELL REceptor signaling)


## Conclusie 


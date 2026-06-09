# Codi RStudio implementat en el Treball Final de Grau 
Repositori que conté tots els scripts en R utilitzats per al Treball de Fi de Grau "Els efectes de la sequera en la comunitat d'amfibis de la XPN de la DIBA". Inclou els scripts dels models TRIM, d’anàlisi descriptiva, models DIM i generació de heatmaps.
L'autora d'aquest TFG és l'alumna Martina Buisán Rodríguez, estudiant del Grau de Biologia de la Universitat de Barcelona. 

A continuació es presenta una descripció individual de cada script, indicant-ne l'objectiu, l'estructura general i el tipus de resultat que genera.

# TRIMS ADULTS  [🔗](./TRIMS_ADULTS.R)
Aquest script conté tot el flux de treball necessari per ajustar els models TRIM per a les espècies amb dades suficients. Inclou:
- Preparació de les taules de comptatge per espècie (site × year × season).
- Ajust del model TRIM amb sobre-dispersió i efectes de lloc i any.
- Automatització del procés per a totes les espècies seleccionades.
- Extracció de resultats: índexs formals i escalats, totals imputats, heatmaps i gràfics de tendència.
- Guardat dels models en format .rds per facilitar consultes posteriors.

Aquest script és el nucli de l’anàlisi de tendències poblacionals adultes i reprodueix exactament els resultats presentats al cos del treball.

# Índex descriptiu ADULTS  [](./index_descriptiu_adults.R)
Aquest script inclou els procediments necessaris per calcular els índexs descriptius anuals d’abundància d’adults per a les espècies que no compleixen els requisits mínims per ajustar un model TRIM. Inclou:
- Preparació de dades: lectura, filtratge i agregació per any i lloc.
- Càlcul de l’índex simple: mitjana d’individus per lloc mostrejat.
- Generació de gràfics: sèries temporals per espècie amb punts i línia.

Aquest script proporciona una alternativa robusta i interpretable quan el volum de dades no permet aplicar models de tendència més complexos.

# DIM EVIDÈNCIES DE REPRODUCCIÓ  [](./DIM_evidencies_reproduccio.R)
Aquest és el script més extens i complex, ja que implementa tot el procés d’ajust dels models dinàmics d’ocupació (DIM) per les poblacions reproductores. Inclou:
- Neteja i imputació de Variables (ECELS, percentatge d’ompliment, hidroperíode, etc.).
- Construcció de matrius de detecció (primàries i secundàries).
- Creació de l’objecte unmarkedMultFrame amb Variables de lloc, d’observació i anuals.
- Ajust del model base i del model complet.
- Estructura del bucle de selecció de models amb criteris estrictes de validesa estadística.
- Extracció dels 10 millors models segons AIC per espècie.
- Matriu de correlacions per avaluar col·linealitat entre variables.

Aquest script reprodueix íntegrament els resultats dels DIM presentats al treball i garanteix la seva replicabilitat.

# Heatmap EVIDÈNCIES DE REPRODUCCIÓ  [](./heatmaps_evidencies_reproduccio.R)
Aquest script genera els heatmaps de presència/absència d’evidències de reproducció per espècie i any. Inclou:
- Construcció de matrius binàries a partir de les dades de camp.
- Ordenació temporal i espacial per facilitar la lectura visual.
- Generació de heatmaps amb codificació de colors per presència/absència.
- Exportació opcional de les figures per incloure-les als annexos.

Aquest script és essencial per visualitzar patrons temporals i espacials de les evidències de reproducció, i complementa els resultats dels DIM.

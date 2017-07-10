# Cieľ

Vývoj algoritmu ktorý by identifikoval objekty a extrahoval ich pozície z astronomických CCD snímok, ktoré obsahujú pozorovania hviezdneho 
pozadia spolu s vesmírnym odpadom ako nefukčné družice, nosné rakety a úlomky satelitov.

# Anotácia

Počas astronomických pozorovaní sa získavajú snímky nočnej oblohy, prevažne jej kokrétnej časti, ktoré sa ukladajú do tzv. Flexible Image Transport System (FITS) formátu. Tieto snímky obsahujú signál rôzneho charakteru od šumu spôsobeného elektrickým prúdom a vyčítavaním snímky zo CCD kamier, cez pozadie oblohy až po skutočné objekty ako hviezdy alebo objekty slnečnej sústavy (asteroidy, kométy, vesmírny odpad, atď.). Každý pixel FITS snímky je charakterisktický svojou pozíciou na CCD kamere (x,y) a intenzitou. Tieto údaje sa využívajú na výpočet polohy objektu na CCD snímke a na jeho súhrnú intenzitu. Na typických astronomických snímkach sa hviezdy javia ako statické body, ktoré možno popísať tzv. rozptýlovou funkciou (z ang. Point Spread Function). To neplatí v prípade, keď sa uskutočnia pozorovania vesmírneho odpadu, ktorý sa pohybuje relatívne rýchlo vzhľadom k hviezdnemu pozadiu. V tomto prípade sa objekty javia ako predlžené čiary a nie ako body. Ak sa počas pozorovaní ďalekohľad pohybuje za objektom vesmírneho odpadu nastáva situácia, že všetky hviezdy sa javia ako predlžené čiary s rovnakou dĺžkou a smerom, zatiaľ čo snímaný objekt sa javí ako bod. Úlohou študenta/-ky bude naštudovať si literatúru venujúcu sa spracovaniu astronomických FITS snímok, ktoré obsahujú objekty vesmírneho odpadu. Následne študent/-ka navrhne najvhodnejší, alebo aj vlastný algoritmus na segmentáciu snímok, ktorý následne naprogramuje a otestuje. Počas segmentácie sa identifikujú všetky objekty na snímke a pre každý taký objekt sa vyextrahuje jeho pozícia na CCD snímke (x,y) ako aj súhrná intenzita. Testovanie algoritmu bude uskutočnené na reálnych snímkach na ktorých sa nachádza hviezdne pozadie ako aj vesmírny odpad. Výsledky sa porovnajú s predpoveďami pozícii vesmírneho odpadu, ktoré budú študentovi dodané spolu s reálnymi snímkami získanými ďalekohľadmi na Astronomickom a geofyzikálnom observatóriu v Modre, FMFI UK.
## diffusion map analysis of Hausgarten pikoplankton data

# Zusammenfassung der Methode:
  
# 1) Daten als Tabelle: Wir wollen eine Tabelle haben die die einzelnen Messungen als Zeilen enth?lt und die 
# verschiedenen Parameter die gemessen werden als Zeilen. Der Methode ist es relativ egal wie diese Tabelle konstruiert wird. Im prinzip eignen sich sogar Bilddateien als eingabedaten, Der Farbwert es ersten pixels ist dann der erste Parameter, der Farbwert des zweiten pixels der zweite parameter und so weiter. Wenn ihr die Daten in anderer Form habt sollte das auch funktionieren.
# 2) Standardtisierung: Wenn die Parameter vergleichbar sind, wie im Beispiel der Bild datei ist hier nicht zu tun, 
# wenn die Parameter sehr verschieden sind wie bei unserern Zensusdaten ist es gut die Spalten so zu normieren das jede 
# spalte den Mittelwert 0 und Varianz 1 bekommt. 
# 3) Lokale Metrik: Jetzt definieren wir ein Mass f?r die ?hnlichkeit von zwei Zeilen D(x,y). Das kann zum Beispiel 
# einfach der Euklidische  Abstand im Parameterraum sein. Aus dem Abstandsma? machen wir ein ?hnlichkeitsma? mit 
# S(x,y)=1/D(x,y). Diese ?hnlichkeit wird f?r all Paare von Zeilen. F?r einen Datensatz mit N Zeilen erhalten wir 
# so eine NxN Matrix von ?hnlichkeiten, S. 
# 4) Thresholding: Wir Vertrauen der ?hnlichkeit nur wenn sie hinreichend gro? ist. In einigen Systemen weiss man 
# wie gro?, gro? genug ist. Es funktioniert aber auch die einfache regel "Eine ?hnlichkeit zwischen zwei Messungen 
# ist signifikant wenn die ?hnlichkeit unter den 10 gr??ten ?hnlichkeiten f?r mindestens eine der Messungen ist". 
# Alle zu kleinen/nicht signifikanten Eintr?ge in der S-Matrix werden null gesetzt.
# 5) Laplace Matrix:  Aus der S-Matrix machen wir eine Laplace Matrix L. Die Formel dazu ist 
# L[i,j]= - S[i,j] f?r i,j verschieden und  L[i,i] = sum_j S[i,j]. Die Originalmethode verwendet die normierte 
# Laplcematrix, in der Zeile i nochmal durch L[i,i] geteilt wird. (Also L[i,i] wird 1) - kann man machen die einfache 
# Laplace matrix funktioniert in meiner Erfahrung genau so gut.  
# 6) Eigenvektoren: Jetzt die Eigenvektoren und eigenwerte der Laplacematrix berechnen. Die eigenvektoren die am 
# n?chsten an der null sind am wichtigsten sie sind die Hauptkomponenten des Datensatzes. Jedes eigenwert eigenvector 
# paar representiert eine Variable die die Methode detektiert hat. Der entsprechende Eigenvector ordnet den einzelnen 
# Messungen/Zeilen der Ausgangstabllen einen Wert in der entsprechenden Variable zu. 



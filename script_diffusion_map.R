## diffusion map analysis of Hausgarten pikoplankton data

rm(list=ls())
# Zusammenfassung der Methode:
  
# 1) Daten als Tabelle: Wir wollen eine Tabelle haben die die einzelnen Messungen als Zeilen enth?lt und die 
# verschiedenen Parameter die gemessen werden als Zeilen. Der Methode ist es relativ egal wie diese Tabelle konstruiert wird. Im prinzip eignen sich sogar Bilddateien als eingabedaten, Der Farbwert es ersten pixels ist dann der erste Parameter, der Farbwert des zweiten pixels der zweite parameter und so weiter. Wenn ihr die Daten in anderer Form habt sollte das auch funktionieren.

# produce matrix to be replaced by data: -------------------

set.seed(423)
vec <- runif(440)

mat <- matrix(vec, ncol=11)

#-----------------------------------------------------------

# 2) Standardisierung: Wenn die Parameter vergleichbar sind, wie im Beispiel der Bild datei ist hier nicht zu tun, 
# wenn die Parameter sehr verschieden sind wie bei unserern Zensusdaten ist es gut die Spalten so zu normieren das jede 
# spalte den Mittelwert 0 und Varianz 1 bekommt. 

# standardize matrix columns to mean = 0 and sd = 1 --------

for(i in 1:ncol(mat)){
  mat[, i] <- (mat[, i] - mean(mat[, i]))/sd(mat[, i]) 
}

# ----------------------------------------------------------

# 3) Lokale Metrik: Jetzt definieren wir ein Mass f?r die ?hnlichkeit von zwei Zeilen D(x,y). Das kann zum Beispiel 
# einfach der Euklidische  Abstand im Parameterraum sein. Aus dem Abstandsma? machen wir ein ?hnlichkeitsma? mit 
# S(x,y)=1/D(x,y). Diese ?hnlichkeit wird f?r all Paare von Zeilen. F?r einen Datensatz mit N Zeilen erhalten wir 
# so eine NxN Matrix von ?hnlichkeiten, S. 

# calculating euclidean distance for matrix mat ------------

eucdists <- matrix(data = NA, nrow = nrow(mat), ncol = nrow(mat))

for (i in 1:nrow(mat)){
  for (j in 1:nrow(mat)){
    eucdists[i,j] <- dist(rbind(mat[i, ], mat[j, ]), method = "euclidean")
  }
}

# calculating degree of similarity
simil <- 1/eucdists
diag(simil) <- 0

#------------------------------------------------------------

# 4) Thresholding: Wir Vertrauen der ?hnlichkeit nur wenn sie hinreichend gro? ist. In einigen Systemen weiss man 
# wie gro?, gro? genug ist. Es funktioniert aber auch die einfache regel "Eine ?hnlichkeit zwischen zwei Messungen 
# ist signifikant wenn die ?hnlichkeit unter den 10 gr??ten ?hnlichkeiten f?r mindestens eine der Messungen ist". 
# Alle zu kleinen/nicht signifikanten Eintr?ge in der S-Matrix werden null gesetzt.

# reducing similarity matrix --------------------------------
simil_red <- matrix(data = NA, nrow = nrow(simil), ncol = nrow(simil))

for (i in 1:nrow(simil)){
  for (j in 1:nrow(simil)){
    if (simil[i,j] < sort(simil[i, ], decreasing = T)[10]){
      simil_red[i,j] <- 0
    }else{
      simil_red[i,j] <- simil[i,j]
    }
  }
}

#-----------------------------------------------------------

simil_red <- matrix(data = NA, nrow = nrow(simil), ncol = nrow(simil))

for (i in 1:nrow(simil)){
  for (j in 1:nrow(simil)){
    if (simil[i,j] < sort(simil[i, ], decreasing = T)[10]){
      simil_red[i,j] <- 0
    }else{
      simil_red[i,j] <- simil[i,j]
    }
  }
}

# 5) Laplace Matrix:  Aus der S-Matrix machen wir eine Laplace Matrix L. Die Formel dazu ist 
# L[i,j]= - S[i,j] f?r i,j verschieden und  L[i,i] = sum_j S[i,j]. Die Originalmethode verwendet die normierte 
# Laplcematrix, in der Zeile i nochmal durch L[i,i] geteilt wird. (Also L[i,i] wird 1) - kann man machen die einfache 
# Laplace matrix funktioniert in meiner Erfahrung genau so gut. 

# calculate laplace-matrix----------------------------------
lap <- matrix(data = NA, nrow = nrow(simil_red), ncol = nrow(simil_red))

for (i in 1:nrow(simil_red)){
  for (j in 1:nrow(simil_red)){
    if (i == j){
<<<<<<< HEAD
      lap[i, j] <- 1
    }else{
      lap[i, j] <- -simil_red[i, j]/sum(simil_red[i,])
=======
      lap[i, j] <- sum(simil_red[, j])
    }else{
      lap[i, j] <- -simil_red[i, j]
>>>>>>> 0c983f18b58b23d880f12190445846236e5a59dd
    }
  }
}

#-----------------------------------------------------------

# 6) Eigenvektoren: Jetzt die Eigenvektoren und eigenwerte der Laplacematrix berechnen. Die eigenvektoren die am 
# n?chsten an der null sind am wichtigsten sie sind die Hauptkomponenten des Datensatzes. Jedes eigenwert eigenvector 
# paar representiert eine Variable die die Methode detektiert hat. Der entsprechende Eigenvector ordnet den einzelnen 
# Messungen/Zeilen der Ausgangstabllen einen Wert in der entsprechenden Variable zu. 

eigen(lap, symmetric = FALSE)
x
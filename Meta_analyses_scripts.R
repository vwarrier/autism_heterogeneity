###Meta_analyses scripts

Betassc = 0.0580253
SEssc = 0.022
wssc = 1/(SEssc*SEssc)


Betaagre = 0.04491
SEagre = 0.05391
wagre = 1/(SEagre*SEagre)


Betameta = (wagre*Betaagre + wssc*Betassc)/(wssc + wagre)
SEmeta = sqrt(1/(wssc + wagre))
Zmeta = Betameta/SEmeta
Pmeta =  2*pnorm(-abs(Zmeta))

Pmeta

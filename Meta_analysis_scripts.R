data1 = read.delim("clipboard")
data1[is.na(data1)] <- 0
data1$w = 1/(data1$SE*data1$SE)
data1$w.1 = 1/(data1$SE.1*data1$SE.1)


data1$Betameta = (data1$Mean*data1$w + data1$Mean.1*data1$w.1 )/(data1$w.1 + data1$w )
data1$SEmeta = sqrt(1/(data1$w.1 + data1$w))
data1$Zmeta = data1$Betameta/data1$SEmeta
data1$Pmeta = 2*pnorm(-abs(data1$Zmeta))
data1$P_BY = p.adjust(data1$Pmeta, method = "BY", n = length(data1$Pmeta))
data1$P_bonferroni = p.adjust(data1$Pmeta, method = "bonferroni", n = length(data1$Pmeta))

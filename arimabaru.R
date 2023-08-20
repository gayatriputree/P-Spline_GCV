library(fpp)
library(dplyr)
library(car)

#input data kurs
Data_Kurs=read.csv("arima.csv", header = TRUE, sep = ";" )
Data_Kurs
min(Data_Kurs$Kurs)
max(Data_Kurs$Kurs)

#Sorting Data Kurs
Kurs_bulan_tahun = Data_Kurs %>% select(Month,Year,Kurs) %>% group_by(Month,Year) %>% summarise(Kurs)
print(Kurs_bulan_tahun, n=42)
Kurs_sorting=arrange(select(Kurs_bulan_tahun,Month,Year,Kurs),Year)
print(Kurs_sorting, n=42)

#Data Time Series Kurs
Kurs_ts=ts(Kurs_sorting$Kurs,start = c(2020,1),end=c(2023,6), frequency = 12)
Kurs_ts

#Plot Time Series Data Kurs
library(TSstudio)
plot(Kurs_ts,xlab = "Waktu", ylab = "Kurs", col="#29648A",lwd=2.5, ylim=c(13700,16000))
ts_plot(Kurs_ts, title = "Harga Nilai Tukar Rupiah (IDR) terhadap Dollar Amerika (USD) Periode 2020 - 2023", 
        Xtitle = "Waktu", Ytitle = "Kurs", line.mode = "lines")

#Uji Stasioneritas
library(tseries)
adf.test(Kurs_ts)
BoxCox.lambda(Kurs_ts)


#Differencing Data Kurs
Kurs_diff=diff(Kurs_ts, differences = 1)
Kurs_diff
adf.test(Kurs_diff)
BoxCox.lambda(Kurs_diff)
plot(Kurs_diff,xlab = "Waktu", ylab = "Kurs", col="#29648A",lwd=2)
abline(h=0, col="#29648A", lwd=1)


acf(Kurs_diff, lag.max = 20)
pacf(Kurs_diff, lag.max = 20)

#ARIMA
Model=auto.arima(Kurs_ts)
Model
Modelarima=arima(Kurs_ts,order=c(2,1,0))
Modelarima

summary(Modelarima)
lmtest::coeftest(Modelarima)
Modelforecast=forecast(Modelarima, h=12)
Modelforecast
plot(forecast(Modelforecast))

#Uji Asumsi Model ARIMA
#White Noise
Box.test(Modelarima$residuals, type = "Ljung-Box")
#Normalitas
shapiro.test(Modelarima$residuals)

#plot data aktual dan ARIMA
A=Kurs_ts
A
B=fitted(Modelarima)
B
plot(A,type = "l", col = "darkblue",lwd=1,xlab="Waktu", ylab = "Kurs")
lines(B, type = "l" ,col = "green",lwd=1)
legend("top", inset=0.05,legend=c("Nilai Kurs Aktual","Nilai Kurs Prediksi ARIMA"), 
       cex=0.7,lty = c(1,1),col=c("darkblue","green"), box.col="white", lwd = 2,
       bg="floralwhite",text.font = 4)

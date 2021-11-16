library(ggplot2)
data(msleep)

tabla = msleep %>% group_by(vore) %>% summarize_at("sleep_total", .funs = list(Media = mean, Desvio = sd)) %>% dplyr::filter(!is.na(vore)) %>% arrange(desc(Media))

sueño=data.frame(medias = tapply(msleep$sleep_total,msleep$vore, mean),
                 desvios = tapply(msleep$sleep_total,msleep$vore, sd))
orden = order(sueño$medias, decreasing=T)
sueño = sueño[orden, ]

## clase 4
## ggplot2
data("mtcars")
ggplot(mtcars, aes(x = wt, y=mpg)) + geom_point() + theme_classic()
# en aes siempre x e y, donde x es la predictora e y es la var respuesta
ggplot(mtcars, aes(x = wt, y=mpg)) + geom_point() + theme_linedraw()
ggplot(mtcars, aes(x = wt, y=mpg)) + geom_point() + theme_minimal()
# siempre q hago referencia a una columna lo tengo q hacer dentro de aes( )
# 
ggplot(mtcars, aes(x = wt, y=mpg)) + geom_point(aes(color=cyl)) + theme_classic()

data("ChickWeight")
ggplot(ChickWeight, aes(x = Time, y = weight))+geom_point()
# no se distingue cuale esta con q dieta, solo q el peso sube con el paso del tiempo
ggplot(ChickWeight, aes(x = Time, y = weight))+geom_point(color = "red")

# Como quiero q cambie el color de acuerdo a otra variable (Diet) y eso es una columna va dentro de "aes"
ggplot(ChickWeight, aes(x = Time, y = weight))+geom_point(aes(color = Diet))
ggplot(ChickWeight, aes(x = Time, y = weight))+geom_point(aes(color = Diet, size= weight))
ggplot(ChickWeight, aes(x = Time, y = weight))+geom_point(aes(color = Diet, shape = Diet, alpha = weight))

# puedo hacer que el tamaño del punto dependa del valor de una variable continua con size
ggplot(mtcars, aes(x = wt, y=mpg)) + geom_point(aes(color=cyl, size = hp)) + theme_classic()

ggplot(iris, aes(x=Species, y = Petal.Length)) + geom_boxplot()
ggplot(iris, aes(x=Species, y = Petal.Length)) + geom_boxplot(aes(color = Species))
ggplot(iris, aes(x=Species, y = Petal.Length)) + geom_boxplot(aes(fill = Species))

ggplot(iris, aes(x=Species, y = Petal.Length)) + geom_point(aes(color = Species))
# para evitar puntos encimados, geom_jitter arma como una franja q corresponde a cada nivel del factor y dentro de esa franja los mueve aleatoriamente a izq o derecha
ggplot(iris, aes(x=Species, y = Petal.Length)) + geom_jitter(aes(color = Species))

ggplot(iris, aes(x=Species, y = Petal.Length)) + geom_violin()
# en este se ve la distribucion si cortamos cada "violin" al medio con una linea vertical. Por ej se ve q setosa tiene distribucion normal (o eso pareciera)
ggplot(iris, aes(x=Species, y = Petal.Length)) + geom_violin()+coord_flip()

# los geom_... se pueden combinar, IMPORTA EL ORDEN!
ggplot(iris, aes(x= Species, y= Sepal.Width)) + geom_violin() + geom_jitter(aes(color=Species))
ggplot(iris, aes(x= Species, y= Sepal.Width)) + geom_jitter(aes(color=Species)) + geom_violin() ## puso ultima la capa de violin con rellenos solido entonces tapa el primer geom
ggplot(iris, aes(x= Species, y= Sepal.Width)) + geom_jitter(aes(color=Species))+ geom_violin(alpha=0.3) 
ggplot(iris, aes(x= Species, y= Sepal.Width)) + geom_jitter(aes(color=Species))+ geom_violin(alpha=0) 

data("diamonds")
# carat son los kilates
# tiene como un millon de puntos.....
ggplot(diamonds, aes(x=carat, y=price))+geom_point()
# geom_hex funciona como un histograma, se cuenta cuantos puntos hay en cada corte (cada hexagono) y el color representa la cantidad de puntos que cayeron dentro
ggplot(diamonds, aes(x=carat, y=price))+geom_hex()
# la mayoria son baratos de pocos kilates
ggplot(diamonds, aes(x=carat, y=price))+geom_smooth()
ggplot(ChickWeight, aes(x = Time, y = weight)) + geom_smooth()
# error muy bajo
# agrego los puntos y veo mas la realidad, los puntos deberian agregarse al final
ggplot(ChickWeight, aes(x = Time, y = weight)) + geom_smooth() + geom_point() 

ggplot(ChickWeight, aes(x = Time, y = weight)) + geom_smooth(aes(color=Diet)) + geom_point(alpha=0.1) 

ggplot(ChickWeight, aes(x = Time, y = weight)) + geom_smooth(aes(fill=Diet)) + geom_point(aes(color=Diet)) 
# geom_smooth usa LOES o LOWESS (una regresion local) para generar la curva de tendencia. Eso se puede cambiar con el argumento method
ggplot(ChickWeight, aes(x = Time, y = weight)) + geom_smooth(aes(fill=Diet), method="lm")
## OJO AHI HIZO 4 MODELOS LINEALES, UNO PARA CADA NIVEL DEL FACTOR. 
## UNO HARIA UN MODELO LINEAL EN TODO CASO CON TIME EXPLICADO POR DIETA PERO CONSIDERANDO TODOS LOS NIVELES EN UN SOLO MODELO.
# 
ggplot(ChickWeight, aes(x = Time, y = weight)) + geom_smooth(method="lm")

ggplot(mtcars, aes(x=wt, y=mpg))+geom_smooth(method="lm")+geom_point()


# facet_wrap. Esta función nos permite generar el gráfico deseado al agregar como argumento dentro de la función el simbolo ~ seguido del nombre de la variable a utilizar para separar los gráficos, 
ggplot(ChickWeight, aes(x = Time, y = weight)) + geom_point(aes(color = Diet)) + 
  geom_line(aes(color = Diet, group = Chick)) + facet_wrap(~Diet)

githubURL <- ("https://raw.githubusercontent.com/derek-corcoran-barrios/derek-corcoran-barrios.github.io/master/Clase4/TempHum.rds")
download.file(githubURL, "TempHum.rds", method = "curl")
library(tidyverse)
TempHum <- read_rds("TempHum.rds") %>% mutate(Mes = as.numeric(Mes))
# filtro Punta Arenas
PA <- TempHum %>% filter(Ciudad_localidad == "Punta Arenas")

ggplot(PA, aes(x = Mes, y = Temperatura)) + geom_point()
ggplot(PA, aes(x = Mes, y = Temperatura)) + geom_point() + stat_smooth(method = "lm")
ggplot(PA, aes(x = Mes, y = Temperatura)) + geom_point() + stat_smooth(method = "lm", formula = y ~ x + I(x^2))



San <- TempHum %>% filter(Ciudad_localidad == "Quinta Normal") %>% 
  pivot_longer(cols = c(Temperatura, Humedad), names_to = "Unidad", 
               values_to = "medida")

ggplot(San, aes(x = Mes, y = medida)) + geom_point() + 
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), 
              aes(fill = Unidad, color = Unidad))
# las q se sacan son q tenian NA
# 
Algunos = TempHum %>% filter(Ciudad_localidad %in% c("Arica", "Rapa Nui", "La Serena", "Valparaíso", "Quinta Normal", "Concepción", "Valdivia", "Punta Arenas"))

# Algunos = Algunos %>% gather(key = Unidad, value = medida, Temperatura, Humedad)
ggplot(Algunos, aes(x=Mes, y =  Temperatura)) + stat_smooth(method="lm", formula =y ~ x + I(x^2), aes(fill= Ciudad_localidad) ) + geom_point(aes(color=Ciudad_localidad))
ggplot(Algunos, aes(x=Mes, y =  Temperatura)) + stat_smooth(method="lm", formula =y ~ x + I(x^2)) + geom_point() + facet_wrap(~Ciudad_localidad, ncol=2)
 
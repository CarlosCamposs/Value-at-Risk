# Cálculo del VaR y tVaR usando cópulas

## Descripción

En este proyecto se plantea la hipótesis de que existe una dependencia entre los rendimientos de las acciones de FEMSA y WALMEX y se quiere cuantificar algunas medidas de riesgo a través de un modelo matemático que considere dicha dependencia. Mediante el uso de cópulas se modela la dependencia de los rendimientos de ambas acciones y se calculan medidas como el Value at Risk (VaR) y Tail Value at Risk (tVaR) para conocer el valor de pérdida máxima que se puede incurrir al ser poseedor de las acciones de ambas compañías en un horizonte de tiempo de un día. 

Finalmente, los resultados obtenidos empleando cópulas se comparan con otros métodos matemáticos para visualizar la similitud que existe y comprobar que el uso de cópulas representa una alternativa viable.

<img src="https://raw.githubusercontent.com/CarlosCamposs/Value-at-Risk/master/imagen/VaR.jpg",height="375",alt="VaR">

## Documentación
El proyecto se desarrolla en los siguientes archivos

- [Archivo PDF](https://github.com/CarlosCamposs/Value-at-Risk/blob/main/(PDF)%20Calculo%20del%20VaR%20y%20tVaR%20usando%20copulas.pdf). Es el documento final donde se extiende la explicación del proyecto así como la justificación y resultados del mismo.

- [Script de R](https://github.com/CarlosCamposs/Value-at-Risk/blob/main/(R)%20Codigo.R). Corresponde al código realizado en RStudio para la elaboración del proyecto 

## Comentarios finales

Se decidió realizar el análisis sobre los rendimientos de dos acciones por la facilidad que se tiene al trabajar con cópulas bivariadas, pero sin duda podría extenderse a considerar un mayor número de acciones. La programación en R para el cálculo de las medidas de riesgo no representó problema alguno, sin embargo al momento de realizar muchas simulaciones se tuvo demoras en la obtención de resultados, posiblemente se podría optimizar el código para que tarde menos al realizar un gran número de simulaciones.

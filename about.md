---
title: <font size="4"> `HDMT` Shiny App </font>
author: <font size="2"> James Dai,  Kevin Wang </font>
date: " "
output:
  html_document
runtime: shiny
---
<br>

Mediation analysis is of rising interest in clinical trials and epidemiology.
Suppose data are collected for an exposures X, a candidate mediators M and a continuous outcome Y. The classical mediation test entails the following two linear regression models for (X,M,Y),

$$\mathbb{E} (M|X) = \alpha_{0} + \alpha_{} X $$
$$\mathbb{E} (Y|M,X) = \beta_{0} + \beta_{1} M_j + \beta X$$  

A mediation relationship stipulates that two conditions have to be met simultaneously: first, the mediator shall be associated with the exposure ($\alpha\neq 0$); second, the mediator has to be associated with the outcome conditional on the exposure ($\beta \neq 0$). This constitutes a composite null hypothesis
$$
H_0: \alpha =0 ,   \beta_{}=0,  
$$
 Properly controlling the type I error rate remains 
a challenge due to this composite null hypothesis, particularly when testing for mediation in high-dimensional studies. 


In this Shiny App package we implement the multiple-testing procedure that was published in Dai et al. (2020) JASA, to accurately controls the family-wise error rate (FWER) and the false discovery rate (FDR) for testing 
high-dimensional mediation composite null hypotheses. 



The `HDMT` Shiny App takes the following user-defined inputs

  - Input dataset: either user input .csv file or a demo dataset `meQTLdata` as one of two examples in the paper. The input dataset shall contain two columns: the first is p-value for testing $\alpha=0$; the second is p-value for testing $\beta=0$. Each row is a molecular feature that is tested for mediation.
  - A multiple-testing error measure: either FDR or FWER.
  - A signficant threshold from 0.05 to 0.2. 
  - A multiple-testing correction method: either `Exact` or `Approximate`. See the JASA paper for details.   
  


After inputting the dataset and the parameters, click `Calculate` to start execution of the algorithm. The main result panel will show a q-q plot, the number of signficant tests based on the user-defined threshold and, if there is signifiance, the `download` button for the result dataset.



#### Reference: 
  - James Dai, Janet Stanford, Michael LeBlanc. A multiple-testing procedure for high-dimensional mediation hypotheses. Journal of the American Statistical Association, 2020. DOI: 10.1080/01621459.2020.1765785.

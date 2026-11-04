# Replicação do modelo semiestrutural de pequeno porte do BCB

O trabalho busca replicar o modelo agregado de pequeno porte do Banco Central do Brasil, na versão do boxe "Atualização dos modelos semiestruturais de pequeno porte" (Relatório de Inflação, junho de 2024).

## Como rodar

1. Instale os pacotes:

```r
install.packages(c("KFAS", "dplyr", "lubridate", "readxl", "ggplot2",
                   "tidyr", "patchwork", "adaptMCMC", "coda"))
```

2. Defina a pasta `dados/` como diretório de trabalho e rode os scripts "1_modelo", "2_lacunas" e "3_irf_modelo" na mesma sessão do R.

- **1_modelo.R**: monta os dados e o modelo, estima o modo e a posteriori (MCMC) e gera as tabelas e gráficos de parâmetros, hiato do produto e taxa de juros neutra.
- **2_lacunas.R**: mostra de onde vêm os valores calibrados e mede quanto cada escolha própria pesa nos resultados, reestimando o modelo com cada alternativa.
- **3_irf_modelo.R**: simula as respostas da inflação a choques de juros, câmbio e demanda, compara-as com as do boxe e gera gráficos.

Os scripts `coleta_*.R` e `consolida_dados.R` só são necessários para refazer os dados, que já estão em `dados/`.

# Replicação do modelo semiestrutural de pequeno porte do BCB

O trabalho busca replicar o modelo agregado de pequeno porte do Banco Central do Brasil, na versão do boxe "Atualização dos modelos semiestruturais de pequeno porte" (Relatório de Inflação, junho de 2024).

## Como rodar

1. Instale os pacotes:

```r
install.packages(c("KFAS", "dplyr", "lubridate", "readxl", "ggplot2",
                   "tidyr", "patchwork", "adaptMCMC", "coda"))
```

2. Defina a pasta `dados/` como diretório de trabalho e rode os scripts "1_modelo", "2_lacunas" e "3_irf_modelo" na mesma sessão do R.

Os scripts `coleta_*.R` e `consolida_dados.R` só são necessários para refazer os dados, que já estão em `dados/`.

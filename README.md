# Replicação do modelo semiestrutural de pequeno porte do BCB

O trabalho busca replicar o modelo agregado de pequeno porte do Banco Central do Brasil, na versão do boxe "Atualização dos modelos semiestruturais de pequeno porte" (Relatório de Inflação, junho de 2024).

## Como rodar

1. Instale os pacotes:

```r
install.packages(c("KFAS", "dplyr", "lubridate", "readxl", "ggplot2",
                   "tidyr", "patchwork", "adaptMCMC", "coda"))
```

2. Defina a pasta `dados/` como diretório de trabalho e rode os scripts "1_modelo", "2_lacunas" e "3_irf_modelo" na mesma sessão do R.

- **1_modelo.R**: lê a base pronta (`dados_modelo.csv`), monta o modelo, estima a moda e a posteriori (MCMC) e gera as tabelas e gráficos de parâmetros, hiato do produto e taxa de juros neutra.
- **2_lacunas.R**: mostra de onde vêm os valores calibrados e mede quanto cada escolha própria pesa nos resultados, reestimando o modelo com cada alternativa.
- **3_irf_modelo.R**: simula as respostas da inflação a choques de juros, câmbio e demanda, compara com as do boxe do BCB e gera gráficos.

Os scripts `coleta_*.R` e `consolida_dados.R` só são necessários para refazer os dados, que já estão em `dados/`.

## Refazer os dados

Rode as coletas, em qualquer ordem, e depois o `consolida_dados.R`. Cada coleta grava um CSV trimestral e mostra gráficos de conferência. Pacotes adicionais:

```r
install.packages(c("rbcb", "sidrar", "seasonal", "mFilter", "jsonlite",
                   "readr", "purrr", "fredr"))
```

A coleta 3 exige uma chave da API do FRED na variável de ambiente `FRED_API_KEY`. Os arquivos `ri202406anp.xlsx` (anexo estatístico do RI de jun/2024) e `dados_novos.xlsx` (CDS e hiato mundial) devem estar na pasta.

- **coleta_1.R** → `dados1.csv`: metas de inflação do CMN, meta para a taxa Selic (SGS 432), IPCA cheio e de preços livres (SGS 433 e 11428), câmbio R$/US$ (SGS 3698) e IC-Br agregado (SGS 27574).
- **coleta_2a.R** → `dados2a.csv`: PIB trimestral (IBGE, SIDRA 1620), expectativas Focus para o IPCA em 12 meses, para a Selic por reunião do Copom e para a Selic de fim de ano, e o calendário do Copom. Usa o cenário de juro real do anexo do RI (Gráfico 2.2.2) para calcular a tendência do juro real.
- **coleta_2b.R** → `dados2b.csv`: taxa de desocupação da PNAD Contínua (SGS 24369), Nuci da indústria de transformação (SGS 28561), estoque de empregos formais do Caged (SGS 28763) e resultado primário do governo federal e BC (SGS 5783).
- **coleta_2c.R** → `dados2c.csv`: IPCA, preços livres e administrados (SGS 433, 11428 e 4449) e expectativas Focus para os administrados e para o IPCA em 12 meses.
- **coleta_3.R** → `dados3.csv`: Fed funds (FRED), CDS de cinco anos do Brasil e hiato do produto mundial de Moraes (2025), ambos de `dados_novos.xlsx`, e o Oceanic Niño Index (NOAA/CPC).
- **consolida_dados.R** → `dados_modelo.csv`: junta as coletas e aplica as transformações do modelo (anualização, termos climáticos, metas, defasagens), de 2003T2 a 2023T4.

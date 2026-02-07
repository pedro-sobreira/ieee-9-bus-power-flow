# Calculo de Fluxo de Potencia - Sistema IEEE 9 Barras

Este repositorio contem um programa em Python para realizar o calculo de fluxo de potencia no sistema de teste IEEE 9 barras, utilizando o metodo iterativo de Newton-Raphson.

## Descricao

O fluxo de potencia e uma ferramenta fundamental na analise de sistemas de energia eletrica, permitindo determinar as tensoes, angulos de fase, potencias ativas e reativas em todas as barras do sistema, bem como os fluxos de potencia nas linhas de transmissao. Este programa implementa o metodo de Newton-Raphson, conhecido por sua robustez e rapida convergencia, para resolver as equacoes nao lineares do fluxo de potencia.

O sistema IEEE 9 barras e um sistema de teste padrao da industria, composto por 9 barras, 3 geradores e 3 cargas, frequentemente utilizado para estudos e validacao de algoritmos de fluxo de potencia.

## Funcionalidades

- Calculo de tensoes (magnitude e angulo) em todas as barras.
- Determinacao das potencias ativas e reativas geradas/consumidas em cada barra.
- Relatorio detalhado do fluxo de potencia ativa e reativa em cada linha de transmissao (ida e volta).
- Calculo das perdas de potencia ativa em cada linha.
- Exibicao clara e formatada dos resultados.

## Como Usar

### Requisitos

Certifique-se de ter o Python 3 e a biblioteca NumPy instalados em seu ambiente.

```bash
pip install numpy
```

### Execucao

1. Clone este repositorio para sua maquina local:
   ```bash
   git clone https://github.com/pedro-sobreira/ieee-9-bus-power-flow.git
   cd ieee-9-bus-power-flow
   ```

2. Execute o programa Python:
   ```bash
   python3 power_flow_ieee9.py
   ```

O programa ira imprimir os resultados diretamente no console.

## Saida do Programa

A saida do programa e dividida em duas secoes principais:

### 1. Relatorio de Tensoes e Potencias nos Barramentos

Esta secao apresenta a magnitude da tensao (V em pu), o angulo de fase (em graus), e as potencias ativas (P_gen em MW) e reativas (Q_gen em MVAr) injetadas em cada barra apos a convergencia do fluxo de potencia.

```
RELATORIO DE TENSOES E POTENCIAS NOS BARRAMENTOS
--------------------------------------------------------------------------------
Barra   V (pu)     Angulo (deg)    P_gen (MW)   Q_gen (MVAr)
--------------------------------------------------------------------------------
1       1.0400     0.0000          71.6410      27.0459     
...
```

### 2. Relatorio Detalhado de Fluxo de Potencia por Conexao

Esta secao fornece um detalhamento do fluxo de potencia para cada barra, mostrando as potencias ativas e reativas que fluem de uma barra para suas barras conectadas. Tambem inclui a geracao e carga locais para um balanco de potencia na barra.

```
RELATORIO DETALHADO DE FLUXO DE POTENCIA POR CONEXAO
================================================================================
BARRA 1:
  Conectada a   P (MW)          Q (MVAr)       
  ---------------------------------------------
  Barra 4        71.6410         27.0459        
  GERACAO       71.6410         27.0459        
  ---------------------------------------------
  TOTAL LIQUIDO 71.6410         27.0459         (Deve ser ~0)
================================================================================
BARRA 4:
  Conectada a   P (MW)          Q (MVAr)       
  ---------------------------------------------
  Barra 1        -71.6410        -23.9231       
  Barra 5        40.9373         22.8931        
  Barra 6        30.7037         1.0300         
  ---------------------------------------------
  TOTAL LIQUIDO -0.0000         0.0000          (Deve ser ~0)
...
```

## Dados do Sistema IEEE 9 Barras

Os parametros do sistema IEEE 9 barras (dados de barra e linha) estao codificados diretamente no arquivo `power_flow_ieee9.py`. Eles representam um sistema padrao para estudos de fluxo de potencia.

## Contribuicao

Contribuicoes sao bem-vindas! Sinta-se a vontade para abrir issues ou pull requests para melhorias, correcoes de bugs ou novas funcionalidades.

## Licenca

Este projeto esta licenciado sob a licenca MIT. Veja o arquivo `LICENSE` para mais detalhes. (Nota: O arquivo LICENSE nao esta incluido neste repositorio, mas e uma boa pratica adiciona-lo.)

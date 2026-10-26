<h1 align="center">Engrenamento hipoide / Simulador tribológico de contato lubrificado</h1>

<p align="center">
  <em>Contato de Hertz, tensões subsuperficiais, filme EHL, térmica, desgaste, fadiga e pitting em um único modelo,<br>auditado contra soluções analíticas e executável no celular.</em>
</p>

<p align="center">
  <img alt="Python 3" src="https://img.shields.io/badge/Python-3-3776AB?logo=python&logoColor=white">
  <img alt="Sem dependências" src="https://img.shields.io/badge/Depend%C3%AAncias-apenas%20biblioteca%20padr%C3%A3o-2E7D32">
  <img alt="Tkinter" src="https://img.shields.io/badge/Interface-Tkinter-555555">
  <img alt="Gráficos em Canvas" src="https://img.shields.io/badge/Gr%C3%A1ficos-Canvas%20puro-1976D2">
  <img alt="Alvo: Pydroid 3" src="https://img.shields.io/badge/Alvo-Pydroid%203%20(Android)-3DDC84?logo=android&logoColor=white">
</p>

<p align="center">
  <img src="graficos.png" alt="Painel de nove gráficos: pressão de contato, tensões em profundidade, filme, Λ, temperatura, integridade, curva S-N, margens por força e evolução por ciclos" width="900">
</p>
<p align="center"><sub>Os nove gráficos do simulador para o caso de referência. A figura foi renderizada a partir dos mesmos comandos de desenho que o Canvas do Tkinter recebe; não é uma captura de tela do aparelho.</sub></p>

---

Este repositório documenta o **modelo**, as **decisões** e a **verificação**. O programa (Python + Tkinter, roda no celular) é o veículo.

O ponto de partida foi um simulador anterior de cerca de 20 mil linhas. Ele executava sem erro e imprimia números plausíveis, mas tinha uma fórmula incorreta no campo de tensões, um fator de carga calculado de duas formas diferentes, funções chamadas e nunca definidas, e muito código que nada usava. O trabalho foi auditar, corrigir, enxugar e validar contra soluções conhecidas, registrando o que ainda **não** está validado.

## Sumário

1. [O problema](#1-o-problema)
2. [Referências](#2-referências)
3. [O trabalho, em quatro etapas](#3-o-trabalho-em-quatro-etapas)
4. [Resultados](#4-resultados)
5. [A interface](#5-a-interface)
6. [Premissas e limites conhecidos](#6-premissas-e-limites-conhecidos)
7. [Execução correta não é lógica correta](#7-execução-correta-não-é-lógica-correta)
8. [Estrutura do repositório](#8-estrutura-do-repositório)
9. [Como executar](#9-como-executar)
10. [Citação e licença](#10-citação-e-licença)

## 1. O problema

Num contato lubrificado entre dois corpos (esfera contra plano, por exemplo), a carga produz uma pressão de Hertz, essa pressão gera um campo de tensões abaixo da superfície, o filme de óleo separa total ou parcialmente as superfícies e o atrito aquece o contato. Cada um desses efeitos limita a vida do componente por um mecanismo diferente:

- **fadiga subsuperficial**, governada pela tensão de von Mises máxima e por uma curva S-N;
- **pitting**, pela pressão média contra a tensão limite de contato;
- **filme insuficiente**, pelo parâmetro Λ = h_min/Rq;
- **temperatura**, por uma elevação bulk e uma flash;
- **desgaste**, pela fração de carga que as asperezas suportam e pela lei de Archard.

O simulador calcula esses efeitos em cadeia a partir de um conjunto de entradas e responde qual limite é atingido primeiro e em quantos ciclos.

## 2. Referências

O modelo usa fórmulas fechadas e correlações da literatura. Nenhuma constante empírica foi ajustada neste trabalho.

| Parte do modelo | Fonte |
|---|---|
| Contato de Hertz e campo de Boussinesq | Johnson, K. L. *Contact Mechanics*. Cambridge University Press, 1985. |
| Espessura mínima e central do filme (h_min, h_c) | Hamrock, B. J.; Dowson, D. Isothermal elastohydrodynamic lubrication of point contacts, Parte III. *ASME J. Lubrication Technology*, 1977. |
| Razão de carga das asperezas (L_a, Eq. 23) | Masjedi, M.; Khonsari, M. M. On the effect of surface roughness in point-contact EHL: formulas for film thickness and asperity load. *Tribology International* 82 (2015) 228–244. |
| Desgaste | Archard, J. F. Contact and rubbing of flat surfaces. *J. Applied Physics* 24 (1953) 981–988. |
| Dano acumulado de fadiga | Miner, M. A. Cumulative damage in fatigue. *J. Applied Mechanics* 12 (1945) A159–A164. |
| Critérios de auditoria do código | Al-Khalili, Z.; Howell, N.; Klakow, D. Evaluating intermediate reasoning of code-assisted large language models for mathematics. *Proc. GEM² 2025*, ACL, pp. 741–758. |

O glossário embutido no programa (87 verbetes) registra a fórmula exata usada para cada grandeza.

## 3. O trabalho, em quatro etapas

### 3.1 Auditoria do simulador anterior

| Achado | Consequência |
|---|---|
| 20.055 linhas, 580 KB; 227 funções inalcançáveis a partir da interface e do caso de referência (10.846 linhas) | Código morto de versões anteriores; impossível de verificar por inteiro |
| Tensões laterais σx, σy e τxy de Boussinesq com fórmula incorreta | σVM máxima cerca de 12% acima do correto no caso de referência; profundidade do máximo e vida por fadiga erradas |
| Fator de carga do pitting: S_H³ na derivação, S_H² no relatório, e S_H,min ignorado | O mesmo conceito aparecia com dois valores (11,7 e 40,2) |
| 4 funções chamadas e nunca definidas na página de histórico | Os botões "adicionar ponto" e "avaliar histórico" falhavam |
| Campos de R1, R2, uₑ, v_slide e frequência desativados na interface, ou lidos e ignorados pelo código | Entradas que pareciam editáveis e não alteravam o cálculo |
| Cálculo executado no mesmo fluxo da interface | Rolagem e barra congeladas durante a simulação |
| Variáveis atribuídas e não usadas, `import` duplicado, versão inconsistente entre docstring e relatório | Ruído que esconde erros reais |

### 3.2 Correções matemáticas

**Tensões de Boussinesq.** O termo `(1 − z/R)·x²/r²` de σx deve ser `(1 − z/R)·(x² − y²)/r²` (e simétrico em σy). Em τxy, o primeiro termo leva coeficiente 2, não 1. A forma correta foi derivada de σr e σθ de Johnson e conferida contra a solução analítica no eixo. O ramo r = 0 do código já usava o limite correto, o que denunciava a inconsistência.

**Pitting.** O fator de carga até o limite passou a ser (S_H/S_H,min)³ nos dois lugares, coerente com p_m ∝ F^(1/3) e com a margem macro.

### 3.3 Reengenharia

```mermaid
flowchart LR
    I["Entradas<br/>34 campos, SI"] --> H["Hertz<br/>E', R', a, p0, pm"]
    H --> S["Tensões<br/>Boussinesq, σVM, cache"]
    S --> E["Filme EHL<br/>h_min, h_c, Λ"]
    E --> T["Atrito e térmica<br/>τ, q, T_bulk, T_flash"]
    T --> W["Desgaste<br/>L_a, Archard"]
    W --> F["Fadiga e pitting<br/>S-N, Miner, S_H"]
    F --> U["Integridade e vida<br/>U, N_lim, mecanismo"]
```

- **Tamanho:** de 20.055 linhas e 580 KB para um arquivo de 3.303 linhas e 121 KB, que já inclui os nove gráficos, o glossário e a gravação em TXT.
- **Desempenho:** o campo de tensões ocupava quase todo o tempo. A integração de Boussinesq foi reescrita sem dicionário nem chamada de função por ponto, mantendo a ordem das operações (4,6 s para 1,9 s, resultado idêntico bit a bit). Em seguida, um cache do campo normalizado σ/p0 reaproveita o resultado quando só F, uₑ, v_slide, T₀ etc. mudam.
- **Interface:** página única em cartões, cálculo em segundo plano com barra de progresso, gráficos que aparecem à medida que ficam prontos, rolagem por toque e gravação do relatório em TXT na pasta Download.
- **Entradas:** todos os valores são editáveis.

### 3.4 Validação e verificação

| Verificação | Resultado |
|---|---|
| E', R', a, p0 e p_m contra as fórmulas fechadas de Hertz | Iguais (erro relativo < 1e-9) |
| ∫ p dA contra F | 1,000002 |
| σz no eixo contra −p0/(1 + z²/a²) | Erro ≤ 0,15% |
| σVM no eixo (ν = 0,3) contra a solução analítica | 0,545 contra 0,5427 (z/a = 0,25); 0,6202 contra 0,6196 (0,5); 0,4712 contra 0,471 (1,0) |
| Máximo de σVM no campo | 0,621·p0 em z/a = 0,5; analítico 0,620·p0 em z/a ≈ 0,48 |
| Antes da correção | σVM no eixo 0,69·p0 em z/a = 0,5 (cerca de 12% acima); σx = −0,109·p0 contra −0,180·p0 analítico |
| Limpeza e otimização contra o original | Relatório idêntico e casos com força, atrito, velocidade, temperatura e ν alterados, todos idênticos bit a bit |
| Cache do campo | Desvio ≤ 2e-14·p0 (arredondamento); relatório idêntico |
| Estrutura do código | Sem nomes indefinidos e sem função morta (verificado por AST) |

A limpeza e a otimização foram verificadas contra o original **antes** das correções matemáticas, para que as duas mudanças não se confundissem.

## 4. Resultados

Caso de referência: F = 1000 N, contato quase plano (R1 = 10⁶ mm, R2 = 72,5 mm), propriedades de aço nos dois corpos (E = 210 GPa, ν = 0,3).

| Grandeza | Valor |
|---|---:|
| Semi-eixo de contato a = b | 0,7782 mm |
| Pressão máxima p0 / média p_m | 788,5 / 525,7 MPa |
| Tensão de von Mises máxima | 489,6 MPa (U_σ = 0,272), em z = 0,389 mm |
| Filme h_min / h_c | 4,250 / 7,002 µm |
| Λ = h_min/Rq | 10,02 |
| Temperatura total | 62,6 °C |
| U_global | 0,741, governado pela temperatura |
| Vida limite | 5,15 × 10¹⁰ ciclos, por fadiga |

Efeito das correções matemáticas sobre o mesmo caso:

| Grandeza | Antes | Depois |
|---|---:|---:|
| σVM máxima | 550 MPa | 490 MPa |
| Profundidade do máximo | 0,195 mm | 0,389 mm |
| U_σ | 0,306 | 0,272 |
| Vida por fadiga | 1,95 × 10¹⁰ ciclos | 5,15 × 10¹⁰ ciclos |

> [!IMPORTANT]
> As entradas de lubrificação, térmica, desgaste e limites do caso de referência são **estimadas**, e o relatório marca `[ESTIMADO]` nas principais. Os números acima mostram o comportamento do modelo, não um componente real. O arquivo [`exemplo_resultado.txt`](exemplo_resultado.txt) é o relatório completo gravado pelo programa para esse caso, com as tabelas dos nove gráficos.

## 5. A interface

Uma página rolável, em cartões: entradas agrupadas (contato e materiais, lubrificação e cinemática, térmica e atrito, desgaste, limites, ciclos/fadiga/pitting), execução com barra de progresso, resultados principais, nove gráficos e o relatório técnico.

| # | Gráfico | Varredura |
|:-:|---|---|
| 1 | Pressão de contato p(r) | raio de −a a +a |
| 2 | Tensões × profundidade no eixo (σz, σr, σVM, τmax) | z de 0,1a a 3a |
| 3 | Espessura de filme h_min e h_c × velocidade | uₑ de 0,2 a 2 vezes |
| 4 | Λ × velocidade | uₑ de 0,2 a 2 vezes |
| 5 | Temperatura × velocidade de deslizamento | v_slide de 0,1 a 2 vezes |
| 6 | Integridade (índices U e dano) | estado atual |
| 7 | Curva S-N com o ponto de operação | σ em torno do ponto de operação |
| 8 | Margens × força (U_σ, U_global, S_H) | F de 0,2 a 3 vezes |
| 9 | Evolução por ciclos (D, U_desgaste, U_global) | N de 0 a 1,2·N_lim |

O programa também traz 87 verbetes técnicos com as fórmulas efetivamente usadas e a gravação do relatório completo em `Download/Tribologico_resultado_V781.txt`.

## 6. Premissas e limites conhecidos

Os limites abaixo fazem parte do resultado.

**Entradas e escopo físico**

- Várias entradas do caso de referência são estimadas (uₑ, η₀, αp, μ, T₀ e os limites). Substitua pelos dados do seu caso.
- A tensão máxima de von Mises vem de uma malha (9 × 9 em x–y, 12 níveis em z); a posição do máximo em z tem passo de 0,25a.
- O campo de tensões usa só a pressão normal. A tração superficial do atrito não entra; com μ alto, as tensões reais perto da superfície são maiores.

**Correlações e simplificações**

- As correlações de h_min, h_c e L_a são empíricas. O modelo **não verifica** se as entradas estão dentro das faixas em que elas foram obtidas.
- A térmica é de forma fechada (flash simplificada e bulk concentrada). Não há acoplamento térmico–EHL.
- A taxa de desgaste é constante no tempo, sem evolução do perfil de contato.
- O coeficiente térmico β é uma entrada que não entra em nenhuma equação.
- Dois critérios de Λ coexistem: U_Λ = 1/Λ (limite em Λ = 1) e Λ mín = h_min,lim/Rq. Eles não são o mesmo limite.
- U_hmin e U_pitting são mostrados, mas não entram no U_global.

**O que foi e o que não foi validado**

- Hertz e Boussinesq foram conferidos contra soluções analíticas. As correlações de EHL, a térmica, o desgaste e a fadiga **não** foram validados contra dados experimentais.
- Toda a verificação automatizada da interface usou *stubs* do Tkinter, nunca Tk real. *(Autor: acrescentar aqui a avaliação visual feita no aparelho.)*
- Desempenho medido só em PC: primeira execução do caso em cerca de 2 s; execuções seguintes, sem mudar ν₁, em milissegundos; gráficos em cerca de 0,2 s. O tempo no celular aparece em "Tempo computacional" no painel de resultados.

## 7. Execução correta não é lógica correta

Este trabalho aplicou ao código os critérios do artigo de Al-Khalili, Howell e Klakow: um programa gerado com assistência de IA pode produzir uma resposta numericamente plausível por um caminho de lógica **não** fundamentada, e métricas de execução não revelam isso. O artigo propõe avaliar a solidez do raciocínio dentro do programa (seis classes, de "conceitual" a "sem lógica") e características estruturais como variáveis definidas e não usadas, usadas e não definidas, número de chamadas e complexidade.

No simulador, o exemplo concreto é o campo de tensões. O programa rodava, o relatório saía completo e os valores tinham ordem de grandeza razoável, mas σx, σy e τxy estavam errados, e só a comparação com a solução analítica no eixo mostrou isso. Daí a regra seguida no trabalho: **cada fórmula crítica é conferida contra uma solução conhecida, e o que não foi conferido está declarado na seção 6.**

O modelo foi desenvolvido com assistência de IA. Os critérios do artigo foram usados como roteiro de auditoria (código morto, variáveis soltas, nomes indefinidos, complexidade, valores fixos sem origem) e não como uma reprodução do método do artigo.

## 8. Estrutura do repositório

```text
.
├── README.md
├── modelo_tribologico.py     # aplicação completa (modelo + interface + gráficos + glossário)
├── graficos.png              # figura exibida no topo deste README
└── exemplo_resultado.txt     # relatório gravado pelo programa para o caso de referência
```

## 9. Como executar

**Requisitos:** Python 3 e Tkinter. Não usa NumPy nem matplotlib. Alvo: Pydroid 3 (Android).

```bash
python modelo_tribologico.py        # interface gráfica
```

No Pydroid 3, abra o arquivo e execute. Ajuste os campos, toque em **EXECUTAR** e acompanhe a barra de progresso. **GRAVAR TXT** salva o relatório e as tabelas dos gráficos em `Download`. **VERBETES TÉCNICOS** abre o glossário.

O caso também pode ser calculado sem interface:

```python
from modelo_tribologico import run_v776_reference_case
resultado, relatorio = run_v776_reference_case()
print(relatorio)
```

## 10. Citação e licença

Se usar este trabalho, cite também as fontes da [seção 2](#2-referências).

> [!IMPORTANT]
> Este é um modelo de estudo. Não é uma ferramenta validada para projeto de componentes.


# Estudo geométrico de engrenagens espirais e hipóides — origem do TRIBOSOFT

## ![Engrenagem hipóide](SBP_2236-1.jpg)

Este trabalho reúne a investigação geométrica e cinemática de engrenagens cônicas espirais e hipóides e o desenvolvimento matemático que antecedeu o simulador tribológico. O software final é uma etapa posterior de um percurso mais amplo: primeiro foram estudadas as superfícies, a geração das curvas e o movimento relativo; depois, essa base foi organizada em um modelo geral de contato e transformada em código.

## 1. Estudo das engrenagens espirais e hipóides

![Desenho geométrico da engrenagem](20260922_123620.jpg)

A investigação aborda a espiral de Arquimedes, sua representação plana e cônica, os ângulos associados à curva e as relações geométricas relevantes para a geração de dentes. Também examina o percurso de contato entre flancos conjugados, distinguindo a trajetória do ponto de contato da linha de contato instantânea.

![Simulação de corte de engrenagem espiral](<Eloid spiral bevel gear cutting sim(360P)_2.gif>)

A comparação entre engrenagens espirais de eixos intersectantes e engrenagens hipóides considera o *offset* entre os eixos e suas consequências para a cinemática e para o deslizamento. O estudo também diferencia o contato geométrico ideal do contato sob carga, no qual a deformação elástica modifica a distribuição de pressão e a região efetiva de contato.

Essa etapa fornece a base conceitual para descrever o contato como um fenômeno geométrico e cinemático, e não apenas como uma característica isolada do perfil do dente.

## 2. Da geometria ao modelo matemático de contato

O desenvolvimento avança para uma estrutura geral aplicável a superfícies engrenadas. A sequência considerada é:

**superfícies parametrizadas → conjugação geométrica e cinemática → linha de contato → movimento relativo → deformação → pressão → microtopologia → comportamento tribológico.**

Essa formulação separa a geometria das superfícies dos fenômenos físicos que surgem quando o contato é carregado e lubrificado. Ela estabelece a ligação entre posição e movimento dos flancos, condições de contato e grandezas necessárias à análise tribológica.

## 3. O software inicial

A primeira implementação computacional dessa estrutura foi concebida como uma ferramenta genérica: não impõe uma geometria física única, mas recebe expressões paramétricas para as duas superfícies. A versão inicial documentada prevê derivadas numéricas, cálculo das normais, estimativa da velocidade relativa, testes de coincidência geométrica e da condição cinemática **v_rel · n = 0**, busca de candidatos à linha de contato e construção de uma representação do volume tribológico em torno desses candidatos.

Portanto, o software inicial não deve ser confundido com o simulador tribológico integrado posterior. Sua função é estrutural: transformar a descrição matemática das superfícies e do movimento em uma base computacional para identificar e caracterizar o contato.

## 4. Evolução para o simulador tribológico

Sobre essa linha de desenvolvimento foram organizadas formulações para contato de Hertz, tensões superficiais e subsuperficiais, lubrificação elastohidrodinâmica (EHL), espessura de filme, atrito, aquecimento, desgaste, fadiga e mecanismos de *pitting*. Essas etapas culminaram no simulador integrado, que reúne entradas, cálculos e resultados em uma única aplicação.

O código chegou a aproximadamente 20 mil linhas durante sua evolução. Esse tamanho registra a trajetória do software, mas não constitui, por si só, uma medida da validade matemática ou da necessidade de cada trecho. A avaliação do código e a verificação das fórmulas são atividades distintas da apresentação do estudo geométrico que lhe deu origem.

## Escopo e ressalva

Este README apresenta a etapa anterior ao simulador final: o estudo das engrenagens espirais e hipóides, a passagem para uma formulação geral do contato e a concepção do software inicial. Não substitui a documentação específica do modelo tribológico nem afirma que todas as formulações do desenvolvimento estejam experimentalmente validadas.

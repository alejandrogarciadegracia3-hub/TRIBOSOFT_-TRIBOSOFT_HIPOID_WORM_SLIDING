# TRIBO_SOFT — do estudo geométrico de engrenagens hipóides ao simulador de contato lubrificado

<p align="center">
  <em>Geometria e cinemática do engrenamento, contato de Hertz, tensões subsuperficiais, filme EHL, atrito, temperatura, desgaste e fadiga,<br>
  verificados contra soluções analíticas, com o que ainda não está validado declarado, e executáveis no celular.</em>
</p>

<p align="center">
  <img alt="Python 3" src="https://img.shields.io/badge/Python-3-3776AB?logo=python&logoColor=white">
  <img alt="Sem dependências" src="https://img.shields.io/badge/Depend%C3%AAncias-apenas%20biblioteca%20padr%C3%A3o-2E7D32">
  <img alt="Tkinter" src="https://img.shields.io/badge/Interface-Tkinter-555555">
  <img alt="Gráficos em Canvas" src="https://img.shields.io/badge/Gr%C3%A1ficos-Canvas%20puro-1976D2">
  <img alt="Alvo: Pydroid 3" src="https://img.shields.io/badge/Alvo-Pydroid%203%20(Android)-3DDC84?logo=android&logoColor=white">
</p>

---

Este repositório conta uma trajetória em quatro fases, e a conta inteira importa: o programa atual só faz sentido sabendo de onde veio e por que cada decisão foi tomada.

1. **O estudo geométrico**: espiral de Arquimedes, engrenagens espirais e hipóides, linha de contato, cinemática.
2. **O MAIN (`MAIN_TRIBOSOFT`)**: a formulação matemática geral do contato tribológico, em 39 seções mais os módulos E1–E5.
3. **O simulador integrado V7.81 e sua auditoria**: cerca de 20 mil linhas que rodavam e imprimiam números plausíveis, até que a conferência com soluções analíticas mostrou fórmulas erradas.
4. **A linha TRIBOCONTATO**: um modelo de um contato elíptico aço–aço lubrificado, construído sobre o MAIN, com validação embutida e, na v11, os módulos que o MAIN previa na teoria mas o código simplificava.

## Sumário

1. [Fase 1 — O estudo geométrico](#fase-1--o-estudo-geométrico)
2. [Fase 2 — O MAIN](#fase-2--o-main)
3. [Fase 3 — O simulador V7.81 e a auditoria](#fase-3--o-simulador-v781-e-a-auditoria)
4. [Fase 4 — TRIBOCONTATO, da v5 à v11](#fase-4--tribocontato-da-v5-à-v11)
5. [O que a v11 muda](#o-que-a-v11-muda)
6. [Como o modelo é verificado](#como-o-modelo-é-verificado)
7. [Premissas e limites](#premissas-e-limites)
8. [Execução correta não é lógica correta](#execução-correta-não-é-lógica-correta)
9. [Como executar](#como-executar)
10. [Estrutura do repositório](#estrutura-do-repositório)
11. [Referências](#referências)
12. [Citação e aviso](#citação-e-aviso)

---

## Fase 1 — O estudo geométrico

O ponto de partida não foi um programa, foi geometria. Foram estudadas a espiral de Arquimedes (plana e cônica), os ângulos associados à curva, a geração dos dentes de engrenagens cônicas espirais e hipóides, e o percurso de contato entre flancos conjugados, distinguindo a **trajetória do ponto de contato** da **linha de contato instantânea**.

A comparação entre espirais de eixos concorrentes e hipóides leva em conta o *offset* entre os eixos e suas consequências para a cinemática e para o deslizamento. O estudo também separa o contato geométrico ideal do contato sob carga, em que a deformação elástica altera a pressão e a região efetiva de contato.

Disso saiu a primeira estrutura de pensamento do projeto:

**superfícies parametrizadas → conjugação geométrica e cinemática → linha de contato → movimento relativo → deformação → pressão → microtopologia → comportamento tribológico.**

A primeira ferramenta computacional foi genérica: recebia expressões paramétricas das duas superfícies, calculava derivadas numéricas e normais, estimava a velocidade relativa, testava a coincidência geométrica e a condição cinemática **v_rel · n = 0**, buscava candidatos à linha de contato e montava o volume tribológico em torno deles. Ela não impunha uma geometria física única e não deve ser confundida com o simulador tribológico posterior.

## Fase 2 — O MAIN

O **MAIN** (`MAIN_TRIBOSOFT.docx`) é o texto matemático que organiza o trabalho. Ele não reduz o comportamento tribológico a um único critério de falha: trata o contato como sistema acoplado, em que carga, geometria, velocidade, materiais, lubrificação, temperatura e superfície se influenciam. Cobre contato hertziano, tensões subsuperficiais, EHL, espessura de filme, razão de filme Λ, atrito, calor, temperatura de flash, desgaste, micropitting, macropitting, scuffing, critérios de integridade e a interação entre os mecanismos.

Dois pontos do MAIN pesam em tudo o que veio depois:

- Ele **já descreve** o acoplamento dinâmico (por exemplo, o desgaste que modifica a geometria real, que modifica a pressão, que modifica o desgaste), mas a formulação computacional operacional (seções 22 e 25) reduz essa cadeia a um fluxo linear e estático.
- A tabela-mestra traz propriedades validadas do aço AISI 8620, mas **não** traz atrito, coeficiente de Archard nem limites de desgaste e temperatura, e a curva S-N de contato tem apenas dois pontos (800 MPa e 700 MPa).

Por isso a qualidade do modelo varia por módulo, e o MAIN foi lido como declaração de limites, não só como fonte de fórmulas.

## Fase 3 — O simulador V7.81 e a auditoria

A primeira implementação integrada chegou a **20.055 linhas e 580 KB**. Rodava sem erro e imprimia valores plausíveis. A auditoria mostrou o custo disso:

| Achado | Consequência |
|---|---|
| 227 funções inalcançáveis a partir da interface e do caso de referência | Código morto de versões anteriores, impossível de verificar por inteiro |
| Tensões σx, σy e τxy de Boussinesq com fórmula incorreta | σVM máxima cerca de 12% acima do correto; profundidade do máximo e vida por fadiga erradas |
| Fator de carga do pitting calculado de duas formas (S_H³ e S_H²) | O mesmo conceito com dois valores |
| 4 funções chamadas e nunca definidas | Botões da página de histórico falhavam |
| Campos da interface lidos e ignorados | Entradas que pareciam editáveis e não alteravam o cálculo |
| Cálculo no mesmo fluxo da interface | Rolagem congelada durante a simulação |

A correção das tensões foi derivada de σr e σθ de Johnson e conferida contra a solução analítica no eixo. Efeito no caso de referência: σVM máxima de 550 para 490 MPa, profundidade do máximo de 0,195 para 0,389 mm, vida por fadiga de 1,95·10¹⁰ para 5,15·10¹⁰ ciclos. A limpeza e a otimização (para 3.303 linhas e 121 KB) foram verificadas contra o original **antes** das correções matemáticas, para que as duas mudanças não se confundissem.

A lição que ficou, e que organiza a fase seguinte: **cada fórmula crítica é conferida contra uma solução conhecida, e o que não foi conferido é declarado.**

## Fase 4 — TRIBOCONTATO, da v5 à v11

A linha `TRIBOCONTATO*.py` é um programa separado do V7.81. Modela **um** contato elíptico aço–aço com deslizamento rotativo, pensado para o hipoide de alta redução, em Pydroid 3 e Tkinter, só com a biblioteca padrão. A cadeia é:

```mermaid
flowchart LR
    I["Entradas<br/>(SI)"] --> H["Hertz elíptico exato<br/>a, b, p0"]
    H --> F["Filme EHL<br/>hmin, hc, Λ"]
    F --> A["Atrito e térmica"]
    H --> S["Tensões<br/>subsuperficiais"]
    S --> Y["Historico do tensor<br/>em ponto material"]
    A --> W["Desgaste e vida"]
    Y --> W
    W --> U["Utilizações e veredito"]
```

O que cada marco acrescentou, conforme o próprio código registra:

| Versão | Marco |
|---|---|
| **v5** | Valida entradas, converte rpm e horas em ciclos por dente, separa desgaste máximo de desgaste médio na faixa. |
| **v9 → v10** | Integral de Boussinesq com malha graduada de Gauss–Legendre e transformação angular que remove a singularidade da borda. Corrige uma divergência identificada da v9 (σz > p0 em z = 0,1c para k = 7) e a mantém preservada só para comparação. Tensor completo (6 componentes) em ponto qualquer e **histórico temporal** em ponto material fixo, ligado ao Miner do MAIN 13.22. |
| **Validação externa** | Comparação com dados publicados de hipoides de alta redução (Nakagawa et al., 2024) e de pontos do ciclo de engrenamento (Paouris et al., 2016), com a divisão calibra/valida do MAIN 27.20. Calibrações são registradas e **não aplicadas ao cálculo**. |
| **v11** | Fadiga multiaxial por profundidade, atrito misto, cota de temperatura e desgaste incremental com realimentação (seção seguinte). |

O resultado nunca é apresentado como previsão absoluta de vida: o veredito validado usa só a tensão de von Mises; filme, fadiga, temperatura e desgaste aparecem como indicativos com nível de confiança explícito.

## O que a v11 muda

A v11 ataca as quatro simplificações do v10 que o MAIN já descrevia como acopladas:

| Simplificação do v10 | Na v11 | Nível |
|---|---|---|
| Desgaste linear sobre elipse fixa, w_N = N·w₁ | Desgaste incremental: o perfil w(y, N) altera a folga, o contato é recalculado em grade, a pressão cai e o desgaste desacelera (Euler modificado, passo controlado por estabilidade) | solver verificado; resultado físico [N] |
| Atrito de Coulomb constante (μ = 0,05) | μ misto: asperezas (forma de Greenwood–Williamson) mais fluido (Eyring) | limites e integrais verificados; parâmetros [N] |
| Blok sem solução para Pe < 10 | Blok mantido; para Pe < 10 informa uma **cota superior estacionária**, não uma temperatura estimada | identidades verificadas |
| S-N de dois pontos sobre a pressão superficial | S-N reescrita em Φ(z) = τ_amp + λ·p_hid,max, com dano D(z) por profundidade | identidades verificadas; λ e H(z) pendentes |

Exemplo (K reduzido a 2·10⁻¹¹ para ficar no regime em que a malha converge): o desgaste no centro é 41,5 µm contra 118,5 µm do cálculo linear estático; a pressão máxima cai de 843 para 456 MPa; a área de contato cresce 1,71 vezes. É o efeito de conformidade que o v10 não podia mostrar.

Os módulos novos são **opt-in** (botões próprios) quando custam mais que o cálculo padrão. A formulação completa está em `MAIN_NOVO_v11.docx`.

## Como o modelo é verificado

Cada teste carrega um nível, e o relatório os exibe:

| Nível | Significado |
|:-:|---|
| **V** | Verificação matemática: identidade ou solução fechada |
| **R** | Referência analítica publicada, comparada numericamente |
| **C** | Calibrado em dado externo |
| **E** | Comparação com dado externo |
| **N** | Não validado |

Três conjuntos de testes executam pela linha de comando:

| Comando | Resultado (ambiente de desenvolvimento) |
|---|---|
| `--validar` | 69 de 69 (V 32 · R 10 · E 27) |
| `--validar-hist` | 47 de 47 |
| `--validar-v11` | 37 de 37 |

Vale ler esses números com ceticismo, e o repositório diz como:

- Parte dos testes são **identidades por construção** (por exemplo, o dano na profundidade crítica igual ao dano em pH depois da calibração). Elas verificam fiação e unidades, não física.
- Os testes com conteúdo físico independente incluem: núcleo de Love contra quadratura; contato em grade contra Hertz exato (erro de 0,3% na pressão máxima); hidrostática no eixo contra Johnson; escalamento de Φ; integral do fluido contra grade cartesiana; limites de Eyring; convergência do desgaste com refino de malha e de passo.
- **As tolerâncias da v11 foram fixadas depois de prototipar o solver**, e as de `--validar` até 21/10/2026 sem registro prévio. Pelo critério 28.26 do MAIN, devem ser tratadas como posteriores ao resultado. Nenhuma tolerância deve ser alterada para um teste passar.
- Não foi usado nenhum dado externo novo na v11.

## Premissas e limites

Fazem parte do resultado:

- **Escopo:** um contato elíptico, não multicontato. Atrito não entra no campo de tensões (só pressão normal).
- **Filme:** hmin vem da fórmula do MAIN (Hamrock–Dowson) ou, no arrastamento angular, de uma variante calibrada nos nove pontos do Paouris (domínio restrito: u_e de 1,5 a 8,7 m/s, F de 230 a 10.600 N, k de 7 a 30).
- **Parâmetros que o MAIN não fornece**, com valores padrão que são **ilustrativos**: atrito, K de Archard, μ das asperezas (0,12), τ₀ de Eyring (5 MPa) e λ da fadiga (0,3).
- **Fadiga:** a S-N tem dois pontos; fora de 700–800 MPa é extrapolação. O perfil de dureza H(z) do aço cementado não está incluído. A margem G_sub do MAIN continua sem R_sub nem Φ_sub.
- **Desgaste incremental:** w depende só de y (ângulo de arrastamento ignorado), N igual nos dois corpos, K e H constantes. Com o K padrão (2·10⁻⁹) o desgaste é da ordem de 10 mm; o laço para em 60 passos e o relatório marca *INCOMPLETO*, porque nesse regime o semiespaço elástico já não se sustenta. A pressão máxima é sensível à malha nas bordas de sulco; o relatório traz também uma medida mais estável (p0eq).
- **Temperatura:** a cota estacionária é um teto, muito conservador em deslizamento rápido, não uma estimativa. Não há acoplamento térmico–EHL, nem realimentação do desgaste para o filme.
- **Interface:** a verificação automatizada nunca usou Tk real. A avaliação visual no aparelho e o tempo de execução no celular ainda precisam ser registrados pelo autor.

## Execução correta não é lógica correta

Um programa gerado com assistência de IA pode produzir uma resposta numericamente plausível por um caminho de lógica não fundamentado, e métricas de execução não revelam isso (Al-Khalili, Howell e Klakow, 2025). No V7.81 o exemplo concreto foi o campo de tensões: o programa rodava, o relatório saía completo e a ordem de grandeza era razoável, mas σx, σy e τxy estavam errados, e só a comparação com a solução analítica mostrou.

O modelo foi desenvolvido com assistência de IA. Os critérios do artigo serviram de roteiro de auditoria (código morto, variáveis soltas, nomes indefinidos, complexidade, valores fixos sem origem), não de reprodução do método. A mesma disciplina vale para a v11: as simplificações são declaradas, as premissas estáticas são nomeadas como tais, e o que é apenas identidade por construção é dito como identidade.

## Como executar

**Requisitos:** Python 3 e Tkinter. Sem NumPy nem matplotlib. Alvo: Pydroid 3 (Android).

```bash
python TRIBOCONTATOGRACIA11.py                 # interface gráfica
python TRIBOCONTATOGRACIA11.py --validar       # validação do núcleo, com dados externos
python TRIBOCONTATOGRACIA11.py --validar-hist  # testes do histórico de tensões
python TRIBOCONTATOGRACIA11.py --validar-v11   # testes dos módulos da v11 (cerca de 8 s num PC)
python TRIBOCONTATOGRACIA11.py --historico     # relatório com histórico e S-N em Φ
```

Na interface: **CALCULAR** roda o cálculo padrão; **HISTÓRICO DE TENSÕES** e **DESGASTE INCREMENTAL** são mais lentos e opcionais; **Mostrar módulos não calibrados** exibe temperatura, atrito, desgaste, vida e rendimento. **SALVAR TXT** grava o relatório na pasta Download.

## Estrutura do repositório

```text
.
├── README.md
├── TRIBOCONTATOGRACIA11.py   # programa atual (modelo + validação + interface Tkinter)
├── MAIN_NOVO_v11.docx        # formulação operacional da v11, com níveis e premissas
├── MAIN_TRIBOSOFT.docx       # MAIN original (39 seções + E1–E5)
└── TRIBOCONTATOGRACIA10.py   # versão anterior, preservada
```

Os estudos geométricos (espiral de Arquimedes, cinemática inversa do cutter) e o simulador V7.81 pertencem às fases anteriores e podem ficar em pastas próprias.

## Referências

| Parte do modelo | Fonte |
|---|---|
| Contato de Hertz, campo de Boussinesq, contato de linha | Johnson, K. L. *Contact Mechanics*. Cambridge University Press, 1985. |
| Filme EHL (h_min, h_c) | Hamrock, B. J.; Dowson, D. *ASME J. Lubrication Technology*, 1977. |
| Filme com arrastamento angular e dados de validação | Paouris, L. et al. *Proc. IMechE Part C* 230 (7–8), 2016. |
| Hipoides de alta redução (dados de validação) | Nakagawa et al. *Trans. JSME* 90 (939), 2024. |
| Relações geométricas e de força do hipoide | Stadtfeld, H. J. *Gear Technology*, ago. 2011. |
| Razão de carga das asperezas | Masjedi, M.; Khonsari, M. M. *Tribology International* 82 (2015) 228–244. |
| Temperatura de flash | Blok (fonte móvel, partição por temperaturas iguais), conforme MAIN seção 11 e 20. |
| Desgaste | Archard, J. F. *J. Applied Physics* 24 (1953) 981–988. |
| Dano acumulado | Miner, M. A. *J. Applied Mechanics* 12 (1945) A159–A164. |
| Auditoria de código gerado com IA | Al-Khalili, Z.; Howell, N.; Klakow, D. *Proc. GEM² 2025*, ACL, pp. 741–758. |

## Citação e aviso

Se usar este trabalho, cite também as fontes acima.

> [!IMPORTANT]
> Este é um modelo de estudo. Não é uma ferramenta validada para projeto de componentes. Os números absolutos de vida e de desgaste devem ser lidos como tendências, não como previsões, até que os parâmetros marcados [N] sejam calibrados com dados do par real.

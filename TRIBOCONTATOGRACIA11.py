# -*- coding: utf-8 -*-
"""
TRIBOCONTATOGRACIA11 - contato lubrificado aco-aco, deslizamento rotativo (versao evolutiva)
Pydroid 3 / Tkinter, so biblioteca padrao. Base: MAIN_TRIBOSOFT (secoes 9-11, 14, tabela-mestra).
Cadeia: Hertz eliptico exato -> filme EHL (Hamrock-Dowson) -> Coulomb -> flash (Blok)
        -> tensoes subsuperficiais -> desgaste (Archard) e vida (S-N) -> utilizacoes.
v5: valida entradas, converte rpm/horas em ciclos por dente e separa desgaste maximo de desgaste medio.
v10: tensor completo em (x,y,z) e historico temporal em ponto material fixo (historico, campo_tensor), integrado ao Miner do MAIN 13.22.
v11: (3) fadiga multiaxial Phi(z) = tau_amp + lambda p_hid,max com D(z) por profundidade; (2) atrito misto (Eyring + asperezas) e cota estacionaria de temperatura;
     (1) desgaste incremental com realimentacao geometria -> pressao (contato em grade, Heun). Premissas estaticas listadas em cada funcao; sem calibracao [N].
Uso: modo grafico (padrao) | python TRIBOCONTATOGRACIA11.py --validar | --validar-hist | --validar-v11 | --historico
"""
import math, os, sys
PI = math.pi
DEG = PI / 180

# ---------------------------------------------------------------- PARAMETROS
# (chave, rotulo, unidade, padrao, fator->SI)
CAMPOS = [
 ("# GEOMETRIA (x = direcao de arrastamento; raio > 0 convexo, < 0 concavo)",),
 ("R1x", "Corpo 1  R1x", "mm", 50.0, 1e-3), ("R1y", "Corpo 1  R1y", "mm", 1000000.0, 1e-3),
 ("R2x", "Corpo 2  R2x", "mm", 60.0, 1e-3), ("R2y", "Corpo 2  R2y", "mm", 200.0, 1e-3),
 ("# CARGA E CINEMATICA (u1,u2 > 0, mesmo sentido)",),
 ("F", "Forca normal F", "N", 1000.0, 1.0),
 ("u1", "Velocidade superficial u1", "m/s", 8.0, 1.0), ("u2", "Velocidade superficial u2", "m/s", 0.3, 1.0),
 ("th", "Angulo de arrastamento (a partir de x)", "graus", 0.0, DEG),
 ("N", "Ciclos no ponto de superficie N", "passagens", 3175e6, 1.0),
 ("# ACO (AISI 8620 cementado)",),
 ("E", "Modulo de Young E", "GPa", 205.0, 1e9), ("nu", "Poisson nu", "-", 0.29, 1.0),
 ("ks", "Condutividade k", "W/mK", 46.6, 1.0), ("cps", "Calor especifico cp", "J/kgK", 475.0, 1.0),
 ("rhos", "Densidade rho", "kg/m3", 7850.0, 1.0),
 ("Ra1", "Rugosidade Ra1", "um", 0.5, 1e-6), ("Ra2", "Rugosidade Ra2", "um", 0.3, 1e-6),
 ("# LUBRIFICANTE (SAE 80W-90 GL-5)",),
 ("nu40", "Viscosidade cinematica 40C", "cSt", 148.0, 1.0), ("nu100", "Viscosidade cinematica 100C", "cSt", 15.5, 1.0),
 ("rho", "Densidade do oleo", "kg/m3", 890.0, 1.0), ("alp", "Coef. pressao-viscosidade", "1e-8/Pa", 2.0, 1e-8),
 ("Tb", "Temperatura de bulk", "C", 75.0, 1.0),
 ("# ATRITO / DESGASTE / LIMITES (nao constam na tabela-mestra do MAIN)",),
 ("mu", "Coeficiente de atrito mu", "-", 0.05, 1.0), ("K", "Coef. de Archard K", "-", 2e-9, 1.0),
 ("H", "Dureza H", "GPa", 7.4, 1e9),
 ("Tlim", "Temperatura limite", "C", 150.0, 1.0), ("hlim", "Filme minimo limite", "um", 0.15, 1e-6),
 ("wlim", "Desgaste limite", "um", 20.0, 1e-6), ("slim", "Tensao von Mises limite", "GPa", 1.5, 1e9),
 ("# ATRITO MISTO E FADIGA MULTIAXIAL v11 (nao constam na tabela-mestra; premissas estaticas)",),
 ("mua", "Atrito limite das asperezas mu_asp", "-", 0.12, 1.0), ("tau0", "Tensao de Eyring tau0", "MPa", 5.0, 1e6),
 ("lam", "Coef. hidrostatico lambda", "-", 0.3, 1.0),
]
PAD = {c[0]: c[3] * c[4] for c in CAMPOS if len(c) > 1}
# Configuracao fixa dos testes (nao muda quando o usuario troca os padroes da ENTRADA)
PAD_REF = dict(PAD, R1x=22e-3, R1y=40e-3, R2x=115e-3, R2y=-200e-3, u1=8.0, u2=5.0)
# (chave, rotulo, unidade, padrao ilustrativo, fator->SI)
CAMPOS_HRH = [
 ("Am", "Distancia conica media Am", "mm", 60.0, 1e-3), ("a", "Offset hipoide a", "mm", 50.0, 1e-3),
 ("bg", "Espiral da coroa beta_g", "graus", 35.0, DEG), ("al", "Angulo de pressao alfa", "graus", 20.0, DEG),
 ("gm", "Angulo de passo do pinhao gama", "graus", 20.0, DEG), ("T", "Torque no pinhao T", "N.m", 10.0, 1.0),
 ("n1", "Rotacao do pinhao", "rpm", 1500.0, 1.0), ("nc", "Pares que dividem a carga", "-", 1.0, 1.0),
 ("bp0", "beta pinhao medido (0 = Stadtfeld)", "graus", 0.0, DEG),
]
# Nakagawa et al. (Trans. JSME 90(939), 2024) Tabela 1: z1, i, Am_p, Am_g [mm], gama_p, delta_g, beta_p, beta_g [graus],
# razao de contato teorica, rendimento da tabela [%], razao de contato real a 100 N.m (lida da Fig. 22, +-0,03)
NAK = {"3/60": (3, 20, 100.73, 94.65, 4.40, 85.15, 59.99, 35.00, 3.84, 75.24, 2.92),
       "6/60": (6, 10, 140.90, 94.95, 6.25, 83.13, 59.75, 35.00, 3.89, 82.65, 2.87),
       "12/60": (12, 5, 117.60, 97.46, 13.45, 75.33, 52.93, 29.00, 3.54, 93.51, 2.48)}
NAK_A, NAK_AL, NAK_T = 40.0, 20.0, 100.0
# Caso A. Gracia 2/60 (dados do usuario, sem saida de referencia): Am, offset [mm], espirais, passo [graus], torque [N.m], rpm, pares
# Dados do usuario: dp = diametros primitivos [mm] (tratados como diametros no ponto de contato); Ae = geratriz do cone primitivo da coroa [mm].
# dp_g/2 = 180 mm <= Ae = 183,4 mm: coerente para um cone primitivo (r = A sen delta).
AGRACIA = dict(dp_p=62.0, dp_g=360.0, Ae=183.4, z1=2, z2=60, a=137.77, bg=32.0, al=20.0, gm=8.9, T=10.0, n1=2500.0, nc=5.0, bp0=79.0)
AM_G_FRAC = 0.857   # Am_coroa/Ae (Nakagawa: 0,857 nos 3 projetos; Paouris: 0,843)
# Geometria da coroa e do pinhao do par 13/36 (Paouris 2016 Tab. 1): distancia conica media da coroa = externa - largura/2
PAO_GEO = dict(Amg=(95.598 - 29.999 / 2) * 1e-3, a=24.0e-3, bg=27.601, dg=59.653, bp=45.989)
# Paouris et al. (Proc IMechE C 230(7-8), 2016), Tab. 3 e Fig. 8: W [N], U, V [m/s], Rzx, Rzy [m], hc, hmin numericos [um],
# pico de pressao de Hertz lido da Fig. 7 [GPa, +-0,02]. Oleo eta0=0,135 Pa.s, alfa=1e-8 /Pa; aco E=210 GPa, nu=0,33; largura da coroa 29,999 mm.
PAO = ((652.14, 2.240, 3.655, 0.0170, 1.198, 0.899, 0.742, 0.62), (1778.11, 2.289, 3.447, 0.0179, 1.222, 0.821, 0.675, 0.85),
       (3016.08, 2.265, 3.069, 0.0183, 1.280, 0.760, 0.624, 0.98), (4679.21, 2.286, 2.922, 0.0191, 1.299, 0.763, 0.634, 1.13),
       (5142.83, 2.326, 2.790, 0.0196, 1.302, 0.790, 0.664, 1.16), (5323.86, 2.346, 2.690, 0.0199, 1.310, 0.796, 0.670, 1.17),
       (5005.89, 2.344, 2.656, 0.0200, 1.320, 0.805, 0.677, 1.15), (4198.77, 2.340, 2.628, 0.0201, 1.330, 0.825, 0.694, 1.05),
       (463.45, 2.448, 2.202, 0.0224, 1.329, 1.148, 0.953, 0.48))   # offset [mm], pressao [graus], torque na coroa [N.m]

# envelope MAIN (tabela-mestra): chave -> (min, max, rotulo)
ENV = {"ph": (500e6, 3000e6, "pH 500-3000 MPa"), "rred": (5e-3, 150e-3, "rho_red 5-150 mm"),
       "ue": (0, 13, "ue <= 13 m/s"), "vs": (0, 12, "v_desl <= 12 m/s"),
       "Ra": (0.3e-6, 0.5e-6, "Ra 0,3-0,5 um"), "nu40": (128, 168, "nu40 128-168 cSt"),
       "nu100": (14, 17, "nu100 14-17 cSt"), "alp": (1.3e-8, 2.8e-8, "alpha 1,3-2,8e-8 /Pa"),
       "Tb": (57, 90, "Tbulk 57-90 C"), "k": (1, 8, "k 1-8 (dominio Hamrock-Dowson)"),
       "Pe": (10, 1e99, "Pe >= 10 (dominio Blok)")}
# Criterios de aceitacao (MAIN 28.26). ATENCAO: ate 21/10/2026 as tolerancias abaixo foram fixadas sem registro previo:
# tratar como posteriores ao resultado. Dado novo exige inserir a tolerancia aqui ANTES de comparar. Identidades usam 1e-9 a 1e-12 no proprio teste.
TOL = {"eta_main": 0.01, "E_main": 0.03, "hertz_desl": 0.01, "k_HB": 0.02, "johnson_sz": 0.005, "johnson_sr": 0.01,
       "johnson_tau": 0.03, "johnson_z": 0.08, "blok": 0.005, "ch_cont": 0.05, "ch_hd": 0.20, "eta_pot": 0.01, "H_HRC": 0.03,
       "hrh_bp": 0.04, "hrh_Fn": 0.005, "hrh_vg": 0.10, "pao_p0": 0.06, "pao_geom": 0.05, "pao_h": 0.20, "tau_lin": 0.01,
       "pao_tau": 0.05, "stress_k15": 0.05, "bp_cal": 0.02, "pao_split": 0.15}
# --- v10: tolerancias dos testes do historico de tensoes. REGISTRADAS EM 25/10/2026 ANTES DA PRIMEIRA EXECUCAO (MAIN 28.26).
# Unidades: "em p0" = erro absoluto dividido por p0; "rel" = erro relativo. Nao alterar para fazer um teste passar:
# alteracao posterior exige registro da data e do motivo no relatorio.
TOL.update({"h_pt": 1e-6,     # equilibrio local do nucleo de Boussinesq (derivadas centrais), normalizado
            "h_id": 1e-10,    # identidades algebricas (traco, invariantes, Miner), rel
            "h_sim": 1e-4,    # simetrias do campo numerico (quadratura graduada nao e simetrica), em p0
            "h_circ": 2e-3,   # circular no eixo vs Johnson, em p0
            "h_axis": 5e-3,   # eixo do campo completo vs tensoes() existente, em p0
            "h_surf": 3e-2,   # condicoes de superficie z=0,01c (sz=p, tau=0), em p0
            "h_ext": 3e-2,    # superficie externa circular vs (1-2nu)p0a^2/(3r^2), z=0,005a, rel
            "h_eq": 2e-2,     # equilibrio local do campo numerico (dif. finitas), normalizado pelo maior gradiente
            "h_eqF": 1e-2,    # equilibrio de forca: integral de sz no plano = F, rel
            "h_conv": 1e-3,   # refinamento de malha, de amostragem no tempo e de janela, rel
            "h_line": 3e-2,   # limite de linha (k=100) vs Johnson 4.45 e vs tau_zx_linha, em p0 / rel
            "h_vmax": 1e-2,   # max global de von Mises do volume vs max no eixo, rel
            "h_dir": 1e-3,    # troca de eixos (arrastamento em y) vs elipse trocada, rel
            "h_phi": 2e-3})   # S-N em phi: p_eq vs pH com k fixo (escala F^(1/3)), rel. REGISTRADA EM 25/10/2026 ANTES DA EXECUCAO
# --- v11: tolerancias dos testes --validar-v11. REGISTRADAS EM 26/10/2026 APOS PROTOTIPAR o solver de grade (escolha de n=10 e dos criterios de passo
# feita com Hertz exato e convergencia a vista): tratar como POSTERIORES ao resultado (MAIN 28.26). Nao alterar para fazer um teste passar.
TOL.update({"v_k": 1e-4, "v_p0": 1e-2, "v_dl": 5e-3, "v_eq": 1e-9, "v_w5": 5e-2, "v_conv": 3e-2, "v_f32": 1e-6, "v_fint": 2e-2, "v_G": 2e-2, "v_Phi": 1e-4})
NIVEL_TXT = "Niveis: [V] validado | [C] calibrado em dado externo | [R] formula do MAIN/literatura, sem teste numerico aqui | [N] nao validado"
SOBRE_TXT = ("Modelo de UM contato eliptico aco-aco lubrificado (deslizamento rotativo), com o MAIN como base e artigos como apoio. "
             "Nucleo validado: Hertz exato, tensoes subsuperficiais, viscosidade, forca do HRH, veredito por von Mises. "
             "Filme (hmin, Lambda): formula do MAIN [R] ou calibrada no Paouris [C]. "
             "Temperatura, atrito, desgaste, vida e rendimento: sem dado, ocultos por padrao. Nao cobre multicontato.")
# S-N: dois pontos do MAIN (secao 9 da tabela de valores gerais)
SN = ((800e6, 2.75e6), (700e6, 6.73e6))
SN_M = math.log(SN[1][1] / SN[0][1]) / math.log(SN[0][0] / SN[1][0])
N_LONGA = 5e7


def avaliar_fadiga_SN(pH, N):
    """Avalia a curva S-N existente no MAIN sem adicionar pontos ou constantes.

    A lei de potencia usa exclusivamente os dois pontos documentados. A vida
    estimada e a pressao admissivel inversa sao calculadas para entradas
    positivas, com classificacao explicita de interpolacao/extrapolacao.
    O marco MAIN de 5e7 ciclos e uma referencia de vida longa, nao um terceiro
    ponto experimental da curva.
    """
    pH = float(pH); N = float(N)
    if not math.isfinite(pH) or not math.isfinite(N):
        raise ValueError("pH e N devem ser finitos.")
    if pH <= 0 or N < 0:
        raise ValueError("pH deve ser > 0 e N deve ser >= 0.")

    Nf = SN[0][1] * (SN[0][0] / pH) ** SN_M
    N_dentro = SN[0][1] <= N <= SN[1][1]
    pH_adm = SN[0][0] * (SN[0][1] / N) ** (1.0 / SN_M) if N > 0 else None
    pH_dentro = pH_adm is not None and SN[1][0] <= pH_adm <= SN[0][0]
    # Sensibilidade analitica da lei Nf = C*(pH_ref/pH)^m.
    # A elasticidade e exata: +1% em pH implica aproximadamente -m% em Nf.
    dNf_dpH = -SN_M * Nf / pH
    return {"Nf_estimado": Nf, "N_dentro_SN": N_dentro,
            "pH_adm_N": pH_adm, "pH_adm_dentro_SN": pH_dentro,
            "faixa_N": "interpolacao" if N_dentro else "extrapolacao",
            "faixa_pH": "interpolacao" if SN[1][0] <= pH <= SN[0][0] else "extrapolacao",
            "margem_pH_N": pH_adm / pH - 1.0 if pH_adm is not None else None,
            "elasticidade_Nf_pH": -SN_M,
            "dNf_dpH": dNf_dpH,
            "vida_sobre_operacao": Nf / N if N > 0 else math.inf}

# ---------------------------------------------------------------- NUCLEO
def kelip(m):
    """Integrais elipticas completas K(m), E(m) (AGM), m < 1."""
    a, b, s, p = 1.0, math.sqrt(1.0 - m), 0.5 * m, 0.5
    for _ in range(40):
        a, b, c = (a + b) / 2, math.sqrt(a * b), (a - b) / 2
        p *= 2
        s += p * c * c
        if abs(c) < 1e-16:
            break
    K = PI / (2 * a)
    return K, K * (1 - s)

def hertz(F, Rx, Ry, Es):
    """Hertz eliptico exato. Devolve a (eixo x), b (eixo y), k=b/a, p0, pm."""
    inv = max(Ry / Rx, Rx / Ry)
    kg = 1.0
    if inv > 1 + 1e-9:
        lo, hi = 1.0, 1e4
        for _ in range(200):
            kg = math.sqrt(lo * hi)
            K, E = kelip(1 - 1 / kg ** 2)
            if (kg ** 2 * E - K) / (K - E) < inv:
                lo = kg
            else:
                hi = kg
        kg = math.sqrt(lo * hi)
    E = kelip(1 - 1 / kg ** 2)[1]
    Rr = 1 / (1 / Rx + 1 / Ry)
    longo = (6 * kg ** 2 * E * F * Rr / (PI * 2 * Es)) ** (1 / 3)   # E' = 2E*
    curto = longo / kg
    a, b = (curto, longo) if Ry >= Rx else (longo, curto)
    return a, b, b / a, 3 * F / (2 * PI * a * b), F / (PI * a * b)

def walther(nu40, nu100, T):
    """ASTM D341: nu(T) em cSt, T em C."""
    x40, x100 = math.log10(313.15), math.log10(373.15)
    y40 = math.log10(math.log10(nu40 + 0.7)); y100 = math.log10(math.log10(nu100 + 0.7))
    B = (y40 - y100) / (x100 - x40)
    A = y40 + B * x40
    return 10 ** (10 ** (A - B * math.log10(T + 273.15))) - 0.7

def filme(F, Rx, k, Ep, eta0, ue, alp):
    """Hamrock-Dowson (MAIN 10.8-10.9)."""
    U = eta0 * ue / (Ep * Rx); G = alp * Ep; W = F / (Ep * Rx * Rx)
    hmin = Rx * 3.63 * U ** 0.68 * G ** 0.49 * W ** -0.073 * (1 - math.exp(-0.68 * k))
    hc = Rx * 2.69 * U ** 0.67 * G ** 0.53 * W ** -0.067 * (1 - 0.61 * math.exp(-0.73 * k))
    return hmin, hc, U, G, W


def filme_ch(F, Re, Rt, Ep, eta0, ue, alp):
    """Chittenden et al. (1985), via tribonet.org: arrastamento no eixo menor (Re<=Rt) ou maior (Re>=Rt)."""
    U = eta0 * ue / (Ep * Re); G = alp * Ep; W = F / (Ep * Re * Re)
    c = U ** 0.68 * G ** 0.49 * W ** -0.073
    if Re <= Rt:
        r = (Rt / Re) ** (2 / 3.0)
        hm, hc = 3.64 * c * (1 - math.exp(-0.67 * r)), 4.30 * c * (1 - math.exp(-1.30 * r))
    else:
        r = Rt / Re
        hm, hc = 3.00 * c * (1 - math.exp(-0.96 * r)), 3.06 * c * (1 - math.exp(-3.36 * r))
    return Re * hm, Re * hc, U, G, W

def hrh(Am, a, bg, al, gm, T, n1, nc, mu, bp0=0.0):
    """Hipoide de alta reducao. Stadtfeld (Gear Technology, ago/2011): beta_p = beta_g + atan(a/Am),
    lead = 90 - beta_p, componentes de forca. Cinematica de linha de contato (derivacao propria, eixos a 90)."""
    bp = bp0 if bp0 > 0 else bg + math.atan(a / Am); lam = PI / 2 - bp
    rm = Am * math.sin(gm); vp = 2 * PI * n1 / 60 * rm
    vg = vp * math.tan(lam); vs = vp / math.cos(lam)
    Ft = T / rm; Fn = Ft / (math.cos(al) * math.cos(bp))
    cb, ca = math.cos(bp) * math.cos(al), math.cos(al)
    Fy = -Ft * (math.sin(gm) * math.sin(bp) * ca + math.cos(gm) * math.sin(al)) / cb
    Fz = -Ft * (math.cos(gm) * math.sin(bp) * ca - math.sin(gm) * math.sin(al)) / cb
    rho = math.atan(mu / math.cos(al))
    return dict(bp=bp, lam=lam, psi=2 * lam, rm=rm, vp=vp, vg=vg, vs=vs, ue=vs / 2,
                u1=vp * math.cos(lam), u2=vp * math.sin(lam) ** 2 / math.cos(lam),
                Ft=Ft, Fn=Fn, Fnc=Fn / nc, Fx=-Ft, Fy=Fy, Fz=Fz, eta=math.tan(lam) / math.tan(lam + rho))

def filme_ang(F, Rx, Ry, th, Ep, eta0, ue, alp):
    """hmin com arrastamento angular (th a partir do eixo x): Paouris et al. 2016, Eq. (34).
    Re, Rs = raios do paraboloide na direcao de arrastamento e transversal; ue = modulo da velocidade de arrastamento."""
    Re = 1 / (math.cos(th) ** 2 / Rx + math.sin(th) ** 2 / Ry)
    Rs = 1 / (math.sin(th) ** 2 / Rx + math.cos(th) ** 2 / Ry)
    U = eta0 * ue / (Ep * Re); G = alp * Ep; W = F / (Ep * Re * Re)
    hm = Re * 3.68 * U ** 0.68 * G ** 0.49 * W ** -0.073 * (1 - math.exp(-0.67 * (Rs / Re) ** (2 / 3.0)))
    return hm, U, G, W

_TZX = []
def tau_zx_linha(n=300):
    """max |tau_xz|/p0 do contato de linha de Hertz (deformacao plana, Johnson), por integracao numerica."""
    if not _TZX:
        h = 2.0 / n; ws = [(-1 + (i + .5) * h, math.sqrt(1 - (-1 + (i + .5) * h) ** 2)) for i in range(n)]
        best = 0.0
        for ix in range(0, 31):
            for iz in range(2, 31):
                x, z = ix / 20.0, iz / 20.0
                v = abs(2 * z * z / PI * sum(p * (x - w) / ((x - w) ** 2 + z * z) ** 2 * h for w, p in ws))
                best = max(best, v)
        _TZX.append(best)
    return _TZX[0]

def _j(X, n=400):
    """int_{-1}^{X} sqrt(1-s^2)/sqrt(X-s) ds  (s = X - w^2, Simpson)."""
    Wm = math.sqrt(X + 1); h = Wm / n; t = 0.0
    for i in range(n + 1):
        w = i * h
        f = math.sqrt(max(1 - (X - w * w) ** 2, 0.0))
        t += f * (1 if i in (0, n) else (4 if i % 2 else 2))
    return 2 * t * h / 3

_JMAX = []
def jmax():
    if not _JMAX:
        _JMAX.append(max(_j(-1 + 2 * i / 200) for i in range(1, 201)))
    return _JMAX[0]

def flash(mu, p0, vs, a, u1, u2, ks, rho, cp):
    """Flash de Blok: fonte movel, fluxo Coulomb semi-eliptico, particao por temperaturas iguais."""
    e = math.sqrt(ks * rho * cp)
    q0 = mu * p0 * vs
    S = e * (math.sqrt(u1) + math.sqrt(u2))
    Tf = q0 * math.sqrt(a) * jmax() / (math.sqrt(PI) * S)
    f1 = e * math.sqrt(u1) / S
    return Tf, q0, f1

def tensoes_v9(a, b, p0, nu, nz=40, nt=60, nf=16):
    """[v9, preservado so para comparacao; tensoes() agora usa campo_tensor] Tensoes no eixo z (Boussinesq integrado sobre p = p0*sqrt(1-rho^2) na elipse).
    Divergencia identificada na v10: quadratura uniforme insuficiente em z pequeno (sz > p0 em z = 0,1c para k=7; erro 5,7e-2 p0, cai com refino).
    Devolve z, sx, sy, sz (compressao positiva), svm, tau."""
    c = min(a, b)
    zs = [c * 0.1 * (30.0 ** (i / (nz - 1))) for i in range(nz)]
    dt, dp = PI / 2 / nt, PI / 2 / nf
    pts = []
    for i in range(nt):
        th = (i + .5) * dt; rho = math.sin(th)
        for j in range(nf):
            ph = (j + .5) * dp
            x, y = a * rho * math.cos(ph), b * rho * math.sin(ph)
            dP = 4 * p0 * a * b * math.cos(th) ** 2 * math.sin(th) * dt * dp
            r2 = x * x + y * y
            pts.append((r2, x * x / r2, dP))
    SX, SY, SZ, SV, TA = [], [], [], [], []
    for z in zs:
        sx = sy = sz = 0.0
        for r2, c2, dP in pts:
            R = math.sqrt(r2 + z * z); R3 = R ** 3
            szt = 3 * z ** 3 / (2 * PI * R3 * R * R)
            sr = ((1 - 2 * nu) / (R * (R + z)) - 3 * r2 * z / (R3 * R * R)) / (2 * PI)
            st = (1 - 2 * nu) * (z / R3 - 1 / (R * (R + z))) / (2 * PI)
            sz += szt * dP
            sx += (sr * c2 + st * (1 - c2)) * dP
            sy += (sr * (1 - c2) + st * c2) * dP
        sx, sy = -sx, -sy                      # compressao positiva
        SX.append(sx); SY.append(sy); SZ.append(sz)
        SV.append(math.sqrt(.5 * ((sx - sy) ** 2 + (sy - sz) ** 2 + (sz - sx) ** 2)))
        TA.append(.5 * (max(sx, sy, sz) - min(sx, sy, sz)))
    return zs, SX, SY, SZ, SV, TA

# ---------------------------------------------------------------- v10: CAMPO DE TENSAO COMPLETO E HISTORICO TEMPORAL
# Convencao: compressao positiva. x = direcao de arrastamento (semi-eixo a), y = transversal (semi-eixo b), z = profundidade (para dentro).
# Tensor (sx, sy, sz, txy, txz, tyz). Mesmo metodo de tensoes(): Boussinesq integrado sobre p = p0*sqrt(1-(x/a)^2-(y/b)^2),
# agora com o tensor completo em ponto qualquer (x,y,z) e malha graduada em torno do ponto de campo (tensoes() so cobre x=y=0).
_GLC = {}
def _gl(n):
    """Nos e pesos de Gauss-Legendre em [-1,1] (Newton sobre P_n)."""
    if n not in _GLC:
        xs, ws = [], []
        for i in range(n):
            x = math.cos(PI * (i + 0.75) / (n + 0.5))
            for _ in range(100):
                p0, p1 = 1.0, x
                for k in range(2, n + 1):
                    p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
                dx = p1 / (n * (x * p1 - p0) / (x * x - 1.0)); x -= dx
                if abs(dx) < 1e-15: break
            p0, p1 = 1.0, x
            for k in range(2, n + 1):
                p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            dp = n * (x * p1 - p0) / (x * x - 1.0)
            xs.append(x); ws.append(2.0 / ((1.0 - x * x) * dp * dp))
        _GLC[n] = (xs, ws)
    return _GLC[n]

def _quebras(c0, s, lo, hi, q):
    """Extremos de paineis em [lo,hi], graduados geometricamente (razao q) em torno de c0, primeira escala s."""
    pts = {lo, hi}
    if lo < c0 < hi: pts.add(c0)
    d = max(s, 1e-12 * (hi - lo))
    while True:
        for v in (c0 - d, c0 + d):
            if lo < v < hi: pts.add(v)
        if c0 - d <= lo and c0 + d >= hi: break
        d *= q
    return sorted(pts)

def _nos(bp, n, L):
    """Gauss-Legendre por painel na variavel angular u = asin(coord/L) (remove a singularidade sqrt do contorno)."""
    xg, wg = _gl(n)
    ang = [math.asin(max(-1.0, min(1.0, v / L))) for v in bp]
    A, W = [], []
    for i in range(len(bp) - 1):
        h = ang[i + 1] - ang[i]
        if h <= 1e-15: continue
        m = 0.5 * (ang[i] + ang[i + 1])
        for xk, wk in zip(xg, wg):
            A.append(0.5 * h * xk + m); W.append(0.5 * h * wk)
    return A, W

def boussinesq(dx, dy, z, nu):
    """Tensor (sx, sy, sz, txy, txz, tyz) por unidade de forca normal compressiva P=1 na superficie, ponto a (dx, dy) do ponto de carga, profundidade z>0.
    Compressao positiva. Obtido de (sr, st, sz, trz) de Boussinesq e rotacionado: sx = sr c^2 + st s^2 etc.; txz = trz c, tyz = trz s."""
    r2 = dx * dx + dy * dy; R2 = r2 + z * z; R = math.sqrt(R2); R3 = R * R2; R5 = R3 * R2
    f = 1.0 / (R * (R + z)); k1 = 1.0 - 2.0 * nu
    sr = k1 * f - 3.0 * r2 * z / R5
    st = k1 * (z / R3 - f)
    sz = -3.0 * z ** 3 / R5
    if r2 > 1e-300:
        c2, s2, sc = dx * dx / r2, dy * dy / r2, dx * dy / r2
        g = -3.0 * z * z / R5                                          # trz/r (o r cancela): txz = g*dx, tyz = g*dy
        txz, tyz = g * dx, g * dy
    else:
        c2, s2, sc, txz, tyz = 1.0, 0.0, 0.0, 0.0, 0.0
    k = -1.0 / (2.0 * PI)                                              # tracao positiva -> compressao positiva
    return (k * (sr * c2 + st * s2), k * (sr * s2 + st * c2), k * sz, k * (sr - st) * sc, k * txz, k * tyz)

def campo_tensor(a, b, p0, nu, x, y, z, n=5, q=3.0):
    """Tensor completo (sx, sy, sz, txy, txz, tyz) [compressao positiva] no ponto (x,y,z), z>0, do contato hertziano eliptico.
    Variaveis (chi, psi): xi = a sin chi, eta = b cos chi sin psi; dP = p0 a b cos^3 chi cos^2 psi dchi dpsi (suave).
    Malha graduada (razao q, n nos Gauss por painel) em torno de (x, y)."""
    if z <= 0.0: raise ValueError("z deve ser > 0")
    S = [0.0] * 6
    chi, wc = _nos(_quebras(x, z, -a, a, q), n, a)
    for ch, wch in zip(chi, wc):
        cc = math.cos(ch); em = b * cc
        if em <= 0.0: continue
        dx = x - a * math.sin(ch)
        psi, wp = _nos(_quebras(y, math.hypot(z, dx), -em, em, q), n, em)
        base = p0 * a * b * cc ** 3 * wch
        for ps, wps in zip(psi, wp):
            cp = math.cos(ps)
            t = boussinesq(dx, y - em * math.sin(ps), z, nu)
            w = base * cp * cp * wps
            for i in range(6): S[i] += t[i] * w
    return tuple(S)

def tensoes(a, b, p0, nu, nz=40, n=5, q=3.0):
    """Tensoes no eixo z (x=y=0), mesmas profundidades e mesma saida do v9, agora por campo_tensor (malha graduada).
    Devolve z, sx, sy, sz (compressao positiva), svm, tau."""
    c = min(a, b)
    zs = [c * 0.1 * (30.0 ** (i / (nz - 1))) for i in range(nz)]
    SX, SY, SZ, SV, TA = [], [], [], [], []
    for z in zs:
        sx, sy, sz = campo_tensor(a, b, p0, nu, 0.0, 0.0, z, n, q)[:3]
        SX.append(sx); SY.append(sy); SZ.append(sz)
        SV.append(math.sqrt(.5 * ((sx - sy) ** 2 + (sy - sz) ** 2 + (sz - sx) ** 2)))
        TA.append(.5 * (max(sx, sy, sz) - min(sx, sy, sz)))
    return zs, SX, SY, SZ, SV, TA

def principais(sx, sy, sz, txy, txz, tyz):
    """Autovalores do tensor simetrico (forma trigonometrica fechada), s1 >= s2 >= s3."""
    q = (sx + sy + sz) / 3.0
    p2 = (sx - q) ** 2 + (sy - q) ** 2 + (sz - q) ** 2 + 2.0 * (txy ** 2 + txz ** 2 + tyz ** 2)
    p = math.sqrt(p2 / 6.0)
    if p < 1e-300: return q, q, q
    bxx, byy, bzz, bxy, bxz, byz = (sx - q) / p, (sy - q) / p, (sz - q) / p, txy / p, txz / p, tyz / p
    det = bxx * (byy * bzz - byz * byz) - bxy * (bxy * bzz - byz * bxz) + bxz * (bxy * byz - byy * bxz)
    ph = math.acos(max(-1.0, min(1.0, det / 2.0))) / 3.0
    e1 = q + 2.0 * p * math.cos(ph); e3 = q + 2.0 * p * math.cos(ph + 2.0 * PI / 3.0)
    return e1, 3.0 * q - e1 - e3, e3

def vm_comp(sx, sy, sz, txy, txz, tyz):
    """von Mises a partir das 6 componentes (MAIN 9.9)."""
    return math.sqrt(0.5 * ((sx - sy) ** 2 + (sy - sz) ** 2 + (sz - sx) ** 2) + 3.0 * (txy ** 2 + txz ** 2 + tyz ** 2))

def historico(R, z0, w0=0.0, nx=41, jan=6.0, n=5, q=3.0, kap=2.2):
    """Historico temporal do tensor em UM ponto material fixo (profundidade z0, deslocamento lateral w0 perpendicular a trajetoria),
    durante a passagem do contato. Quase-estatico, so pressao normal (atrito mu nao entra: [N] no modelo).
    Trajetoria no referencial do contato: direcao e = (cos th, sin th) [th = R['th'], de x]; s = coordenada ao longo de e, de -L a +L,
    com L = jan*max(l, z0), l = semi-corda da elipse na direcao e. Amostragem s = L sinh(kap u)/sinh(kap), u uniforme (mais densa no contato).
    Tempo do ponto material do corpo i: t = (s + L)/u_i."""
    a, b, p0, nu = R["a"], R["b"], R["p0"], R["nu"]
    th = R.get("th", 0.0); ct, st = math.cos(th), math.sin(th)
    ell = 1.0 / math.sqrt(ct * ct / (a * a) + st * st / (b * b))
    L = jan * max(ell, z0)
    S, T = [], []
    for j in range(nx):
        u = -1.0 + 2.0 * j / (nx - 1)
        S.append(L * math.sinh(kap * u) / math.sinh(kap))
    comp = [[] for _ in range(6)]; svm, tsz, tmax, pr, ph = [], [], [], [], []
    for s in S:
        x, y = s * ct - w0 * st, s * st + w0 * ct
        t = campo_tensor(a, b, p0, nu, x, y, z0, n, q)
        for i in range(6): comp[i].append(t[i])
        ph.append((t[0] + t[1] + t[2]) / 3.0)
        svm.append(vm_comp(*t)); e = principais(*t); pr.append(e); tmax.append(0.5 * (e[0] - e[2]))
        tsz.append(t[4] * ct + t[5] * st)
    H = dict(s=S, L=L, ell=ell, z=z0, w0=w0, th=th, sx=comp[0], sy=comp[1], sz=comp[2], txy=comp[3], txz=comp[4], tyz=comp[5],
             svm=svm, tsz=tsz, tmax=tmax, pr=pr, ph=ph)
    # Extremos refinados (secao aurea no intervalo vizinho ao extremo amostrado): picos de tsz sao estreitos (largura ~ z) e a
    # amostragem uniforme os subestimava (ate 2% em z = 0,2c com nx = 81). Sem refino o resultado dependeria de nx.
    def f_(s_, k):
        t_ = campo_tensor(a, b, p0, nu, s_ * ct - w0 * st, s_ * st + w0 * ct, z0, n, q)
        return (t_[4] * ct + t_[5] * st) if k == "tsz" else vm_comp(*t_)
    def ext_(vals, k, sinal):
        i = max(range(nx), key=lambda m: sinal * vals[m])
        if i in (0, nx - 1): return vals[i]
        lo, hi = S[i - 1], S[i + 1]; g = 0.5 * (math.sqrt(5.0) - 1.0)
        x1, x2 = hi - g * (hi - lo), lo + g * (hi - lo); f1, f2 = sinal * f_(x1, k), sinal * f_(x2, k)
        for _ in range(22):
            if f1 > f2: hi, x2, f2 = x2, x1, f1; x1 = hi - g * (hi - lo); f1 = sinal * f_(x1, k)
            else: lo, x1, f1 = x1, x2, f2; x2 = lo + g * (hi - lo); f2 = sinal * f_(x2, k)
        return sinal * max(f1, f2, sinal * vals[i])
    H["ext"] = dict(tsz_max=ext_(tsz, "tsz", 1), tsz_min=ext_(tsz, "tsz", -1), svm_max=ext_(svm, "svm", 1))
    for k in ("u1", "u2"):
        if R.get(k, 0) > 0: H["t" + k[1]] = [(s + L) / R[k] for s in S]
    return H

def resumo_hist(H):
    """Medidas do ciclo H_c (MAIN 13.4). tsz = cisalhamento ortogonal (plano trajetoria-profundidade), amplitude = (max-min)/2."""
    tz = H["tsz"]; vm = H["svm"]
    mx = max(max(abs(v) for v in H[k]) for k in ("sx", "sy", "sz", "txy", "txz", "tyz"))
    cauda = max(max(abs(H[k][0]), abs(H[k][-1])) for k in ("sx", "sy", "sz", "txy", "txz", "tyz")) / mx if mx > 0 else 0.0
    e = H["ext"]
    return dict(z=H["z"], svm_max=e["svm_max"], svm_min=min(vm), tsz_amp=0.5 * (e["tsz_max"] - e["tsz_min"]), tsz_max=e["tsz_max"], tsz_min=e["tsz_min"],
                tresca_max=max(H["tmax"]), phyd_max=max(H["ph"]), i_vm=vm.index(max(vm)), s_vm=H["s"][vm.index(max(vm))], cauda=cauda)

def zs_hist(R, nz=14):
    """Profundidades da varredura: mesma faixa de tensoes() (0,1c a 3c, geometrica), c = min(a,b)."""
    c = min(R["a"], R["b"])
    return [c * 0.1 * (30.0 ** (i / (nz - 1))) for i in range(nz)]

def varrer_z(R, w0=0.0, zs=None, **kw):
    """Resumo do historico em cada profundidade. Devolve lista de resumos e indices criticos (por amplitude de tsz e por von Mises)."""
    zs = zs or zs_hist(R); V = [resumo_hist(historico(R, z, w0, **kw)) for z in zs]
    return V, max(range(len(V)), key=lambda i: V[i]["tsz_amp"]), max(range(len(V)), key=lambda i: V[i]["svm_max"])

def dano_miner(niveis):
    """MAIN 13.22: D = sum n_i/N_i, com N_i da curva S-N do MAIN (avaliar_fadiga_SN, dois pontos documentados). niveis = [(n_i, pH_i)]."""
    D = 0.0
    for ni, pHi in niveis:
        D += ni / avaliar_fadiga_SN(pHi, ni)["Nf_estimado"]
    return D

def pressao_equivalente(Phi, rho_ref):
    """S-N em Phi = amplitude de tau ortogonal do historico. Os dois pontos do MAIN (800 MPa/2,75e6; 700 MPa/6,73e6) valem para o contato de
    referencia, cuja razao e rho_ref = Phi_ref/pH_ref. Um Phi qualquer corresponde a p_eq = Phi/rho_ref, e Nf = Nf_SN(p_eq) (mesma lei do MAIN)."""
    return Phi / rho_ref

def dano_miner_phi(niveis, rho_ref):
    """MAIN 13.22 com a S-N em Phi: niveis = [(n_i, Phi_i)]; D = sum n_i/Nf(Phi_i/rho_ref)."""
    return dano_miner([(ni, pressao_equivalente(Phi, rho_ref)) for ni, Phi in niveis])

# ---------------------------------------------------------------- v11: ATRITO MISTO (MAIN 19, 25.11), COTA ESTACIONARIA, DESGASTE INCREMENTAL (MAIN 14.38, 22.16, 25.13, 25.20)
# Premissas ESTATICAS (valores fixos durante o calculo, nao evoluem): tau0, mu_asp, lambda, K, H, N igual nos dois corpos, filme hf = hc (ou hmin sem hc), eta(p) de Barus em Tb.
def _f32(t, n=240):
    """E[(X-t)+^(3/2)], X ~ N(0,1) (asperezas gaussianas, Greenwood-Williamson). s = t + u^2 remove a singularidade; Simpson em u de 0 a 3,4."""
    um = 3.4; h = um / n; tot = 0.0
    for i in range(n + 1):
        u = i * h; s = t + u * u
        f = 2.0 * u ** 4 * math.exp(-0.5 * s * s) / math.sqrt(2.0 * PI)
        tot += f * (1 if i in (0, n) else (4 if i % 2 else 2))
    return tot * h / 3.0

def alfa_asp(Lam):
    """Fracao de carga nas asperezas: forma de GW (gaussiana), normalizada a 1 em Lambda = 0. Forma [N]: a magnitude absoluta nao e calibrada."""
    return _f32(max(Lam, 0.0)) / _f32(0.0)

def tau_eyring(p, hf, vs, eta0, alp, tau0):
    """Cisalhamento do fluido: tau = tau0 asinh(eta gamma_ponto / tau0), eta = eta0 exp(alp p) (Barus em Tb), gamma_ponto = vs/hf."""
    eta = eta0 * math.exp(min(alp * p, 700.0))
    return tau0 * math.asinh(eta * (vs / hf) / tau0)

def atrito_misto(R, P, n=80):
    """mu medio ponderado pela carga: mu = alfa mu_asp + (1-alfa) (int tau_f dA)/F, sobre p = p0 sqrt(1-rho^2).
    int tau_f dA = 2 pi a b int_0^1 tau_f(p0 s) s ds (s = sqrt(1-rho^2)). Devolve dict; todas as saidas [N]."""
    hf = R["hc"] if R.get("hc") else R["hmin"]
    al = alfa_asp(R["Lam"]); h = 1.0 / n; tot = 0.0
    for i in range(n + 1):
        s = i * h
        tot += tau_eyring(R["p0"] * s, hf, R["vs"], R["eta0"], P["alp"], P["tau0"]) * s * (1 if i in (0, n) else (4 if i % 2 else 2))
    tint = 2.0 * PI * R["a"] * R["b"] * tot * h / 3.0
    muf = tint / R["F"]; mud = al * P["mua"] + (1.0 - al) * muf
    return dict(mu_din=mud, alfa_asp=al, mu_fluido=muf, hf=hf)

def G_estac(a, b, n=360):
    """G = (1/8) int_0^{2pi} rho(phi) dphi, rho = semi-corda da elipse na direcao phi (centro: T = q0 G/k por corpo, semiespaco, fonte q0 sqrt(1-rho^2)). Circulo: pi a/4."""
    h = 0.5 * PI / n; s = 0.0
    for i in range(n):
        ph = (i + .5) * h
        s += 1.0 / math.sqrt(math.cos(ph) ** 2 / (a * a) + math.sin(ph) ** 2 / (b * b))
    return 4.0 * s * h / 8.0

def flash_estac(mu, p0, vs, a, b, ks):
    """Limite estacionario no centro (conducao pura, dois corpos iguais, particao 1/2): T = q0 G / (2 k). Cota superior da fonte movel com o mesmo fluxo
    (nucleo de Green movel = nucleo estacionario x exp(-v(R+xi)/2kappa) <= 1). Substitui 'n/a' quando Pe < 10; nao e a temperatura de Blok."""
    return mu * p0 * vs * G_estac(a, b) / (2.0 * ks)

def _Iq(X, Y):
    """Primitiva dupla de 1/r, impar em X e em Y."""
    if X == 0.0 or Y == 0.0: return 0.0
    ax, ay = abs(X), abs(Y)
    v = ax * math.asinh(ay / ax) + ay * math.asinh(ax / ay)
    return v if X * Y > 0 else -v

def _rect(x1, x2, y1, y2):
    """int int dxi deta / r sobre [x1,x2]x[y1,y2], ponto de campo na origem."""
    return _Iq(x2, y2) - _Iq(x1, y2) - _Iq(x2, y1) + _Iq(x1, y1)

def matriz_influencia(nx, ny, dx, dy, Es):
    """Deslocamento normal combinado no centro da celula (i,k) por pressao unitaria na celula (j,l), quadrante x,y>0 com as 3 imagens espelhadas
    (simetria do desgaste). u = (1/(pi E*)) int p/r dA (mesma convencao de _u_hertz)."""
    c = 1.0 / (PI * Es)
    K = [[c * _rect((i - .5) * dx, (i + .5) * dx, (j - .5) * dy, (j + .5) * dy) for j in range(2 * ny + 1)] for i in range(2 * nx + 1)]
    n = nx * ny; A = [[0.0] * n for _ in range(n)]
    for i in range(nx):
        for k in range(ny):
            Ar = A[i * ny + k]
            for j in range(nx):
                for l in range(ny):
                    Ar[j * ny + l] = K[abs(i - j)][abs(k - l)] + K[i + j + 1][abs(k - l)] + K[abs(i - j)][k + l + 1] + K[i + j + 1][k + l + 1]
    return A

def _chol(M):
    m = len(M); L = [[0.0] * m for _ in range(m)]
    for i in range(m):
        Li = L[i]
        for j in range(i + 1):
            s = M[i][j]; Lj = L[j]
            for q in range(j): s -= Li[q] * Lj[q]
            if i == j:
                if s <= 0.0: raise ValueError("matriz de influencia nao positiva")
                Li[j] = math.sqrt(s)
            else: Li[j] = s / Lj[j]
    return L

def _chsolve(L, b):
    m = len(L); y = [0.0] * m
    for i in range(m):
        s = b[i]; Li = L[i]
        for q in range(i): s -= Li[q] * y[q]
        y[i] = s / Li[i]
    x = [0.0] * m
    for i in range(m - 1, -1, -1):
        s = y[i]
        for q in range(i + 1, m): s -= L[q][i] * x[q]
        x[i] = s / L[i][i]
    return x

def contato_grade(g, A, dA4, F, S0):
    """Contato sem atrito por conjunto ativo: A p = delta - g nas celulas com p>0, p >= 0, nao-penetracao fora, 4 dA sum p = F. g = folga por celula.
    Devolve p, delta, iteracoes. Cholesky de A_SS resolvido para (1) e (-g); delta pela restricao de forca."""
    n = len(g); S = sorted(S0)
    for it in range(80):
        if not S: S = [min(range(n), key=lambda i: g[i])]
        m = len(S); L = _chol([[A[i][j] for j in S] for i in S])
        z1 = _chsolve(L, [1.0] * m); z2 = _chsolve(L, [-g[i] for i in S])
        delta = (F / dA4 - sum(z2)) / sum(z1)
        pS = [z2[q] + delta * z1[q] for q in range(m)]
        if min(pS) < 0.0:
            S = [S[q] for q in range(m) if pS[q] >= 0.0]; continue
        Sset = set(S); tol = 1e-9 * abs(delta); viol = []
        for i in range(n):
            if i not in Sset and sum(A[i][S[q]] * pS[q] for q in range(m)) + g[i] - delta < -tol: viol.append(i)
        if not viol:
            p = [0.0] * n
            for q in range(m): p[S[q]] = pS[q]
            return p, delta, it
        S = sorted(Sset | set(viol))
    raise RuntimeError("contato por conjunto ativo nao convergiu")

def _interp(xs, vs, x, fora=None):
    """Interpolacao linear; abaixo do primeiro no = primeiro valor; acima do ultimo = fora (ou o ultimo valor se None)."""
    if x <= xs[0]: return vs[0]
    if x >= xs[-1]: return vs[-1] if fora is None else fora
    lo, hi = 0, len(xs) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if xs[mid] <= x: lo = mid
        else: hi = mid
    t = (x - xs[lo]) / (xs[hi] - xs[lo]); return vs[lo] + t * (vs[hi] - vs[lo])

def desgaste_incremental(R, n=10, gam=0.5, passos_max=60, borda=0.02, fdom=1.35, M=200, tol_p=0.35):
    """MAIN 14.38/22.16/25.13/25.20: desgaste -> geometria real -> pressao -> desgaste. Folga g(x,y) = x^2/2Rx + y^2/2Ry + w(y) - delta (Hertz em grade, quadrante).
    Archard por passagem (ambos os corpos somados): dw = K vs/H (1/u1 + 1/u2) int p dx. ESTATICO assumido: w depende so de y (trilha ao longo de x; th ignorado,
    como no calcular) e N e a mesma contagem de passagens nos dois corpos. Integracao: Euler modificado (Heun) com w em tabela fina em y (independente da grade);
    passo limitado por dw <= gam (w_max + 0,1 delta0) e por variacao relativa da taxa <= tol_p (estabilidade). Dominio da grade: x livre, y adaptativo.
    p_max e a pressao numerica normalizada pelo Hertz exato do contato novo (erro de grade ~0,3%); em bordas de sulco p_max e sensivel a malha [N]."""
    Rx, Ry, Es, F = R["Rx"], R["Ry"], R["Es"], R["F"]
    a0, b0, p00 = R["a"], R["b"], R["p0"]
    Kc = R["K"] * R["vs"] / R["H"] * (1.0 / R["u1"] + 1.0 / R["u2"]); Ntot = R["N"]
    T = dict(Ymax=5.0 * b0, wt=[0.0] * (M + 1)); S = dict(Lx=fdom * a0, Ly=fdom * b0, ex=a0, ey=b0); cache = {}
    def grade(Lx, Ly):
        if (Lx, Ly) not in cache:
            if len(cache) > 8: cache.clear()
            cache[(Lx, Ly)] = (Lx / n, Ly / n, matriz_influencia(n, n, Lx / n, Ly / n, Es))
        return cache[(Lx, Ly)]
    def ynos(): return [T["Ymax"] * j / M for j in range(M + 1)]
    def estender(Ly):
        while Ly > T["Ymax"]:
            ys = ynos(); nov = 2.0 * T["Ymax"]
            T["wt"] = [_interp(ys, T["wt"], nov * j / M, 0.0) for j in range(M + 1)]; T["Ymax"] = nov
    def extensao(p, dx, dy):
        pm = max(p); ix = iy = 0
        for i in range(n):
            for k in range(n):
                if p[i * n + k] > 0.01 * pm: ix = max(ix, i); iy = max(iy, k)
        return (ix + 1) * dx, (iy + 1) * dy, pm
    def resolver(wt):
        for _ in range(14):
            Lx, Ly = S["Lx"], S["Ly"]; estender(Ly); dx, dy, A = grade(Lx, Ly); ys = ynos()
            wv = [_interp(ys, wt if wt is not None else T["wt"], (k + .5) * dy) for k in range(n)]
            g = [((i + .5) * dx) ** 2 / (2 * Rx) + ((k + .5) * dy) ** 2 / (2 * Ry) + wv[k] for i in range(n) for k in range(n)]
            ex, ey = S["ex"], S["ey"]
            S0 = [i * n + k for i in range(n) for k in range(n) if ((i + .5) * dx / ex) ** 2 + ((k + .5) * dy / ey) ** 2 < 1.0] or [0]
            p, d = contato_grade(g, A, 4 * dx * dy, F, S0)[:2]
            exn, eyn, pm = extensao(p, dx, dy); S["ex"], S["ey"] = exn, eyn; mud = False
            if max(p[(n - 1) * n + k] for k in range(n)) > borda * pm: S["Lx"] = Lx * 1.5; mud = True
            elif exn < 0.55 * Lx: S["Lx"] = fdom * exn; mud = True
            if max(p[i * n + n - 1] for i in range(n)) > borda * pm: S["Ly"] = Ly * 1.5; mud = True
            elif eyn < 0.55 * Ly: S["Ly"] = fdom * eyn; mud = True
            if not mud: break
        return p, d, dx, dy
    def taxas(p, dx, dy):
        rc = [Kc * 2 * dx * sum(p[i * n + k] for i in range(n)) for k in range(n)]
        yc = [(k + .5) * dy for k in range(n)]; Ly = n * dy
        return [(_interp(yc, rc, y) if y <= Ly else 0.0) for y in ynos()], max(rc)
    hist = []
    p, d, dx, dy = resolver(None); pnum0, delta0 = max(p), d
    def registrar(N, p, dx, dy, d):
        ex, ey, pm = extensao(p, dx, dy)
        hist.append(dict(N=N, w0=T["wt"][0], pmax=p00 * pm / pnum0, p0eq=4.0 / 3.0 * sum(v * v for v in p) / sum(p),
                         area=4 * dx * dy * sum(1 for v in p if v > 0), xa=ex, yb=ey, delta=d))
    registrar(0.0, p, dx, dy, d)
    N = 0.0; cnt = 0; rej = 0
    while N < Ntot * (1 - 1e-12) and cnt < passos_max:
        cnt += 1
        r1, rmax = taxas(p, dx, dy)
        dN = min(Ntot - N, gam * (max(T["wt"]) + 0.1 * delta0) / rmax) if rmax > 0 else Ntot - N
        w_old = list(T["wt"])
        for _ in range(14):
            p2, d2, dx2, dy2 = resolver([w_old[j] + dN * r1[j] for j in range(M + 1)])
            r2, rm2 = taxas(p2, dx2, dy2)
            if rmax <= 0 or max(abs(r2[j] - r1[j]) for j in range(M + 1)) / rmax <= tol_p: break
            dN *= 0.5; rej += 1
        T["wt"] = [w_old[j] + 0.5 * dN * (r1[j] + r2[j]) for j in range(M + 1)]
        N += dN
        p, d, dx, dy = resolver(None); registrar(N, p, dx, dy, d)
    D = 0.0
    for h0, h1 in zip(hist[:-1], hist[1:]):
        D += (h1["N"] - h0["N"]) * 0.5 * (1.0 / avaliar_fadiga_SN(h0["pmax"], 1.0)["Nf_estimado"] + 1.0 / avaliar_fadiga_SN(h1["pmax"], 1.0)["Nf_estimado"])
    Ly = n * dy; ys = [y for y in ynos() if y <= Ly]; sh = max(1.0 / R["u1"], 1.0 / R["u2"]) / (1.0 / R["u1"] + 1.0 / R["u2"])
    wlin = Ntot * (R["w1"] + R["w2"])
    return dict(hist=hist, passos=cnt, rejeitados=rej, n=n, completo=N >= Ntot * (1 - 1e-9), N_alcancado=N, w0=T["wt"][0], wN_din=T["wt"][0] * sh,
                w_lin_soma=wlin, razao_w=(T["wt"][0] / wlin if wlin > 0 else 1.0), D_inc=D, delta0=delta0,
                px=[((i + .5) * dx, p[i * n]) for i in range(n)], py=[((k + .5) * dy, p[k]) for k in range(n)],
                prof=[(y, T["wt"][j]) for j, y in enumerate(ynos()) if y <= Ly], pnum0=pnum0)

def fadiga_historico(R, w0=0.0, nz=10, rho_ref=None, **kw):
    """v11 (MAIN 13.22, 9.12): Phi(z) = tau_amp(z) + lambda p_hid,max(z) em cada profundidade (tipo Crossland/Dang Van; tau_amp = amplitude de tau ortogonal do
    historico, p_hid,max = max do (sx+sy+sz)/3 no ciclo). S-N do MAIN reescrita em Phi (p_eq = Phi/rho_ref, mesma lei de dois pontos). D(z) = N/Nf(p_eq(z)),
    Miner de um nivel por profundidade; D_phi = max_z D(z). ESTATICO assumido: lambda constante (entrada, nao consta no MAIN); rho_ref calibrado no contato nominal
    (com isso D no nominal nao depende de lambda); dureza H(z) do perfil cementado NAO entra (o MAIN nao fornece o perfil): Nf igual em todas as profundidades [N].
    G_sub = R_sub - Phi_sub (MAIN 9.12, 13.16) segue sem R_sub nem Phi_sub: nao calculado."""
    V, ic, iv = varrer_z(R, w0, zs_hist(R, nz), **kw)
    lam = R.get("lam", 0.0)
    Phiz = [v["tsz_amp"] + lam * v["phyd_max"] for v in V]
    ic = max(range(len(V)), key=lambda i: Phiz[i]); Phi = Phiz[ic]
    cal = rho_ref is None                       # sem contato de referencia dado: calibra no proprio contato nominal (premissa estatica, ver motivo)
    if cal: rho_ref = Phi / R["p0"]
    peq = pressao_equivalente(Phi, rho_ref); sn = avaliar_fadiga_SN(peq, R["N"])
    Dz = [R["N"] / avaliar_fadiga_SN(pressao_equivalente(max(f, 1e-30), rho_ref), R["N"])["Nf_estimado"] for f in Phiz]
    return dict(V=V, i_tau=ic, i_vm=iv, D_main=dano_miner([(R["N"], R["p0"])]), Phi=Phi, Phiz=Phiz, Dz=Dz, lam=lam, zs=[v["z"] for v in V],
                rho_ref=rho_ref, rho_calibrado_no_nominal=cal,
                p_eq=peq, Nf_phi=sn["Nf_estimado"], D_phi=R["N"] / sn["Nf_estimado"], faixa_phi=sn["faixa_pH"], G_sub=None, Phi_sub=None,
                motivo="S-N do MAIN reescrita em Phi(z) = tau_amp + lambda p_hid,max (lambda de entrada, nao do MAIN). O MAIN nao informa a geometria dos ensaios de 800/700 MPa: "
                       "premissa estatica = esses ensaios correspondem ao contato nominal (rho_ref = Phi/p0 do proprio caso). Sem H(z) do perfil cementado. G_sub = R_sub - Phi_sub segue sem R_sub nem Phi_sub no MAIN.")

def calcular(P, tens=True, hist=False, inc=False):
    R = dict(P)
    R["Rx"] = 1 / (1 / P["R1x"] + 1 / P["R2x"]); R["Ry"] = 1 / (1 / P["R1y"] + 1 / P["R2y"])
    R["rred"] = 1 / (1 / R["Rx"] + 1 / R["Ry"])
    R["Es"] = P["E"] / (2 * (1 - P["nu"] ** 2)); R["Ep"] = 2 * R["Es"]
    a, b, k, p0, pm = hertz(P["F"], R["Rx"], R["Ry"], R["Es"])
    R.update(a=a, b=b, k=k, p0=p0, pm=pm)
    R["ue"] = (P["u1"] + P["u2"]) / 2; R["vs"] = abs(P["u1"] - P["u2"])
    R["nuTb"] = walther(P["nu40"], P["nu100"], P["Tb"])
    R["eta0"] = R["nuTb"] * 1e-6 * P["rho"]
    if abs(P.get("th", 0.0)) > 1e-9 and R["Rx"] > 0 and R["Ry"] > 0:
        hmin, U, G, W = filme_ang(P["F"], R["Rx"], R["Ry"], P["th"], R["Ep"], R["eta0"], R["ue"], P["alp"])
        hc = None; R["mh"] = "Paouris Eq.34, arrastamento angular"
    elif k >= 1:
        hmin, hc, U, G, W = filme(P["F"], R["Rx"], k, R["Ep"], R["eta0"], R["ue"], P["alp"]); R["mh"] = "Hamrock-Dowson"
    else:
        hmin, hc, U, G, W = filme_ch(P["F"], R["Rx"], R["Ry"], R["Ep"], R["eta0"], R["ue"], P["alp"]); R["mh"] = "Chittenden, eixo maior"
    R.update(hmin=hmin, hc=hc, U=U, G=G, W=W)
    R["Rq"] = 1.25 * math.hypot(P["Ra1"], P["Ra2"]); R["Lam"] = hmin / R["Rq"]
    R["Tf"], R["q0"], R["f1"] = flash(P["mu"], p0, R["vs"], a, P["u1"], P["u2"], P["ks"], P["rhos"], P["cps"])
    kap = P["ks"] / (P["rhos"] * P["cps"])
    R["Pe"] = min(P["u1"], P["u2"]) * a / kap
    if R["Pe"] < ENV["Pe"][0]:
        R["Tf"] = None                       # Blok nao se aplica (Pe < 10): sem numero
    R["Tc"] = P["Tb"] + R["Tf"] if R["Tf"] is not None else None
    R.update(atrito_misto(R, P))                 # v11: mu dinamico (misto) [N]; mu de entrada segue valendo para Tf e Tc acima
    R["Tf_din"] = flash(R["mu_din"], p0, R["vs"], a, P["u1"], P["u2"], P["ks"], P["rhos"], P["cps"])[0]
    if R["Pe"] < ENV["Pe"][0]: R["Tf_din"] = None
    R["Tc_din"] = P["Tb"] + R["Tf_din"] if R["Tf_din"] is not None else None
    R["Tup"] = flash_estac(P["mu"], p0, R["vs"], a, b, P["ks"]); R["Tup_din"] = flash_estac(R["mu_din"], p0, R["vs"], a, b, P["ks"])
    # Archard local: integra p(x)=p0*sqrt(1-(x/a)^2) ao longo do eixo x.
    # Integral de -a a +a = pi*p0*a/2. N significa passagens do ponto de
    # superficie pelo contato, não número total de dentes engrenando no sistema.
    # O mesmo N só representa ambos os corpos se for a contagem de passagens
    # do ponto material avaliado em cada corpo.
    integral_p_dx = 0.5 * PI * p0 * a
    # Profundidade no centro da faixa (y=0): estimativa conservadora maxima.
    R["w1"] = P["K"] * R["vs"] * integral_p_dx / (P["u1"] * P["H"])
    R["w2"] = P["K"] * R["vs"] * integral_p_dx / (P["u2"] * P["H"])
    # Media sobre a largura da elipse: integral de (1-y^2/b^2) = 2/3 do centro.
    R["wmean1"] = (2.0 / 3.0) * R["w1"]
    R["wmean2"] = (2.0 / 3.0) * R["w2"]
    R["wN1"] = P["N"] * R["w1"]
    R["wN2"] = P["N"] * R["w2"]
    R["wN"] = max(R["wN1"], R["wN2"])
    R["wNmean"] = P["N"] * max(R["wmean1"], R["wmean2"])
    R["Nw"] = P["wlim"] / max(R["w1"], R["w2"]) if max(R["w1"], R["w2"]) > 0 else math.inf
    R.update(avaliar_fadiga_SN(p0, P["N"]))
    R["Nf"] = R["Nf_estimado"]
    R["D"] = P["N"] / R["Nf"]
    R["u_T"] = (R["Tc"] - P["Tb"]) / (P["Tlim"] - P["Tb"]) if R["Tc"] is not None else None
    R["u_h"] = P["hlim"] / hmin; R["u_w"] = R["wN"] / P["wlim"]; R["u_D"] = R["D"]
    if tens:
        zs, sx, sy, sz, sv, ta = tensoes(a, b, p0, P["nu"])
        i = sv.index(max(sv)); j = ta.index(max(ta))
        R.update(zs=zs, sx=sx, sy=sy, sz=sz, sv=sv, ta=ta, svm=sv[i], zsv=zs[i], tmax=ta[j], ztm=zs[j])
        R["u_s"] = R["svm"] / P["slim"]
    if hist:                                    # opt-in: o historico custa ~35x o calculo padrao (1,5 s contra 0,04 s no PC)
        R["hist"] = fadiga_historico(R)
        R["Phi"], R["p_eq"], R["Nf_phi"], R["D_phi"] = R["hist"]["Phi"], R["hist"]["p_eq"], R["hist"]["Nf_phi"], R["hist"]["D_phi"]
    if inc:                                     # opt-in v11: desgaste incremental com realimentacao (contato em grade)
        R["inc"] = desgaste_incremental(R)
    return R

def avisos(R):
    out = []
    def f(ok, msg):
        if not ok: out.append(msg)
    for ch, v in (("ph", R["p0"]), ("rred", R["rred"]), ("ue", R["ue"]), ("vs", R["vs"]),
                  ("nu40", R["nu40"]), ("nu100", R["nu100"]), ("alp", R["alp"]), ("Tb", R["Tb"]),
                  ("Pe", R["Pe"])):
        f(ENV[ch][0] <= v <= ENV[ch][1], "fora do envelope MAIN: " + ENV[ch][2])
    if abs(R.get("th", 0.0)) > 1e-9:
        out.append("arrastamento angular: hmin calibrado no Paouris 2016 (so vale com ue 1,5-8,7 m/s, F 230-10600 N, k 7-30); hc nao validado; flash segue o eixo x")
    if R["Tf"] is None:
        out.append("flash de Blok n/a: Pe < 10 (conducao domina); temperatura sem valor")
    if R["k"] < 1:
        out.append("k<1: arrastamento no eixo maior da elipse; hmin por Chittenden (fora do dominio de Hamrock-Dowson)")
    elif R["k"] > 8:
        out.append("fora do envelope MAIN: " + ENV["k"][2])
    out.append("v11: atrito misto, cota estacionaria e desgaste incremental usam premissas estaticas (tau0, mu_asp, lambda, K, H, N igual nos dois corpos, w so em y); sem calibracao [N].")
    out.append("v5: desgaste Archard e vida S-N sao estimativas condicionais; K e curva S-N nao calibrados para o par real. Fora de 700-800 MPa, Nf e extrapolado. N e passagens do ponto material.")
    for ch in ("Ra1", "Ra2"):
        f(ENV["Ra"][0] - 1e-12 <= R[ch] <= ENV["Ra"][1] + 1e-12, "fora do envelope MAIN: " + ENV["Ra"][2] + " (" + ch + ")")
    return out

def incerteza(P, n=150):
    """Monte Carlo (MAIN 29.21) sobre as faixas da tabela-mestra: nu40/nu100 (totalmente correlacionados), alfa, Tb, Ra1, Ra2."""
    import random
    rg = random.Random(1); o = {"n": n, "hmin": [], "Lam": [], "Tc": []}
    for _ in range(n):
        u = rg.random()
        Q = dict(P, nu40=ENV["nu40"][0] + u * (ENV["nu40"][1] - ENV["nu40"][0]),
                 nu100=ENV["nu100"][0] + u * (ENV["nu100"][1] - ENV["nu100"][0]), alp=rg.uniform(*ENV["alp"][:2]),
                 Tb=rg.uniform(*ENV["Tb"][:2]), Ra1=rg.uniform(*ENV["Ra"][:2]), Ra2=rg.uniform(*ENV["Ra"][:2]))
        r = calcular(Q, False)
        o["hmin"].append(r["hmin"] * 1e6); o["Lam"].append(r["Lam"])
        if r["Tc"] is not None: o["Tc"].append(r["Tc"])
    return o

def nivel_h(R):
    """hmin e Lambda: C com angulo de arrastamento dentro do dominio dos 9 pontos do Paouris ampliado 2x
    (ue 1,5-8,7 m/s; F 230-10600 N; k 7-30); N com angulo fora desse dominio; R (formula do MAIN) com arrastamento colinear."""
    kk = max(R["k"], 1 / R["k"])
    if abs(R.get("th", 0.0)) > 1e-9:
        return "C" if (1.5 <= R["ue"] <= 8.7 and 230 <= R["F"] <= 10600 and 7 <= kk <= 30) else "N"
    return "R"

def faixa(x):
    v = sorted(x); m = len(v) - 1
    return v[int(0.05 * m)], v[int(0.5 * m)], v[int(0.95 * m)]

def _pao_modelo(row):
    W, U, V, Rx, Ry, hcn, hmn, pf = row
    Q = dict(PAD_REF); Q.update(E=210e9, nu=0.33, nu40=0.135 / 890e-6, nu100=0.135 / 890e-6, rho=890.0, alp=1e-8, Tb=40.0, R2x=1e12, R2y=1e12,
                           F=W, R1x=Rx, R1y=Ry, u1=math.hypot(U, V), u2=math.hypot(U, V), th=math.atan2(V, U))
    return calcular(Q, False)

def calibracoes():
    """K_BP: razao beta_p(tabela)/beta_p(Stadtfeld) media dos 3 pares do Nakagawa. C_h: fator de hmin ajustado nos pontos impares do Paouris."""
    ks = []
    for k, (z1, ii, Amp, Amg, gp, dg, bpt, bgt, eps, eta, ar) in NAK.items():
        bpe = hrh(Amg * 1e-3, NAK_A * 1e-3, bgt * DEG, NAK_AL * DEG, dg * DEG, 1.0, 60.0, 1.0, 0.05)["bp"]
        ks.append(bpt * DEG / bpe)
    odd = [PAO[i] for i in (0, 2, 4, 6, 8)]
    ch = sum(r[6] / (_pao_modelo(r)["hmin"] * 1e6) for r in odd) / len(odd)
    return sum(ks) / len(ks), ch

# ---------------------------------------------------------------- VALIDACAO
def avisos_f(R, exp):
    return [a for a in avisos(R) if exp or not ("Pe" in a or "flash" in a)]

def _u_hertz(x0, y0, a, b, p0, Es, nr=200, nf=72):
    """Deslocamento normal combinado u = (1/(pi E*)) int p/r dA (polar centrado no ponto)."""
    tot = 0.0
    for j in range(nf):
        ph = (j + .5) * 2 * PI / nf; c, s = math.cos(ph), math.sin(ph)
        A = c * c / a ** 2 + s * s / b ** 2; B = x0 * c / a ** 2 + y0 * s / b ** 2
        C = x0 ** 2 / a ** 2 + y0 ** 2 / b ** 2 - 1
        rm = (-B + math.sqrt(B * B - A * C)) / A
        sj = 0.0
        for i in range(nr):
            r = (i + .5) / nr * rm
            q = 1 - ((x0 + r * c) / a) ** 2 - ((y0 + r * s) / b) ** 2
            sj += p0 * math.sqrt(max(q, 0.0))
        tot += sj * rm / nr * 2 * PI / nf
    return tot / (PI * Es)

def validar():
    L = []
    def t(nome, val, ref, tol, rel=True):
        tol = TOL[tol] if isinstance(tol, str) else tol
        err = abs(val - ref) / abs(ref) if rel else abs(val - ref)
        L.append((nome, val, ref, err, err <= tol, None, tol))
    P = dict(PAD_REF); R = calcular(P)
    # 1 lubrificante
    t("Walther reproduz nu40", walther(148, 15.5, 40), 148, 1e-9)
    t("Walther reproduz nu100", walther(148, 15.5, 100), 15.5, 1e-9)
    t("eta40(128 cSt) vs MAIN 0,114", 128e-6 * 890, 0.114, "eta_main")
    t("eta40(168 cSt) vs MAIN 0,150", 168e-6 * 890, 0.150, "eta_main")
    t("eta100(14 cSt) vs MAIN 0,0125", 14e-6 * 890, 0.0125, "eta_main")
    t("eta100(17 cSt) vs MAIN 0,015", 17e-6 * 890, 0.015, "eta_main")
    # 2 material
    t("E* = E/(2(1-nu^2)) vs MAIN E'~115 GPa", R["Es"], 115e9, "E_main")
    # 3 Hertz circular (Johnson): a=(3FR/4E*)^(1/3), p0=(6FE*^2/(pi^3 R^2))^(1/3)
    Rc, F, Es = 0.02, 200.0, R["Es"]
    a, b, k, p0, pm = hertz(F, Rc, Rc, Es)             # Rx=Ry=Rc (esfera equivalente)
    t("Hertz circular: a", a, (3 * F * Rc / (4 * Es)) ** (1 / 3), 1e-9)
    t("Hertz circular: p0", p0, (6 * F * Es ** 2 / (PI ** 3 * Rc ** 2)) ** (1 / 3), 1e-9)
    # 4 Hertz eliptico: relacao exata e Hamrock-Brewe
    a, b, k, p0, pm = hertz(R["F"], R["Rx"], R["Ry"], Es)
    Kk, Ee = kelip(1 - 1 / k ** 2)
    t("Elipse: (k^2E-K)/(K-E) = Ry/Rx", (k * k * Ee - Kk) / (Kk - Ee), R["Ry"] / R["Rx"], 1e-9)
    t("Elipse: k exato vs Hamrock-Brewe 1,0339(Ry/Rx)^0,636", k, 1.0339 * (R["Ry"] / R["Rx"]) ** 0.636, "k_HB")
    t("Equilibrio: F = int p dA = (2/3) pi a b p0", 2 * PI * a * b * p0 / 3, R["F"], 1e-12)
    # 5 condicao de Hertz independente: u(0,0)-u(x,0) = x^2/(2Rx) ; u(0,0)-u(0,y) = y^2/(2Ry)
    u0 = _u_hertz(0, 0, a, b, p0, Es)
    t("Hertz (deslocamento) em x=a/2", u0 - _u_hertz(a / 2, 0, a, b, p0, Es), (a / 2) ** 2 / (2 * R["Rx"]), "hertz_desl")
    t("Hertz (deslocamento) em y=b/2", u0 - _u_hertz(0, b / 2, a, b, p0, Es), (b / 2) ** 2 / (2 * R["Ry"]), "hertz_desl")
    # 6 tensoes: circular vs Johnson (nu=0,3)
    nu, a1, p1 = 0.3, 1e-3, 1e9
    zs, sx, sy, sz, sv, ta = tensoes(a1, a1, p1, nu)
    i = min(range(len(zs)), key=lambda q: abs(zs[q] - 0.5 * a1))
    z = zs[i]
    t("sigma_z(z) circular vs p0/(1+z^2/a^2)", sz[i], p1 / (1 + z * z / a1 ** 2), "johnson_sz")
    sr = p1 * ((1 + nu) * (1 - z / a1 * math.atan(a1 / z)) - 1 / (2 * (1 + z * z / a1 ** 2)))
    t("sigma_r(z) circular vs Johnson", sx[i], sr, "johnson_sr")
    t("tau_max circular vs 0,31 p0", max(ta), 0.31 * p1, "johnson_tau")
    t("z(tau_max) circular vs 0,48 a", zs[ta.index(max(ta))], 0.48 * a1, "johnson_z")
    # 7 filme: expoentes HD e limite
    h1 = filme(R["F"], R["Rx"], k, R["Ep"], R["eta0"], R["ue"], P["alp"])[0]
    h2 = filme(R["F"], R["Rx"], k, R["Ep"], R["eta0"], 2 * R["ue"], P["alp"])[0]
    t("hmin ~ ue^0,68", h2 / h1, 2 ** 0.68, 1e-9)
    h3 = filme(2 * R["F"], R["Rx"], k, R["Ep"], R["eta0"], R["ue"], P["alp"])[0]
    t("hmin ~ F^-0,073", h3 / h1, 2 ** -0.073, 1e-9)
    # 8 flash: constante de Blok 1,11 (MAIN 11.6) e particao
    t("Blok: 2*jmax*sqrt(2)/pi^1,5 vs 1,11", 2 * jmax() * math.sqrt(2) / PI ** 1.5, 1.11, "blok")
    ap = 0.5e-3; u1, u2 = 9.0, 4.0
    Tf, q0, f1 = flash(0.05, 1e9, 5.0, ap, u1, u2, 46.6, 7850, 475)
    e = math.sqrt(46.6 * 7850 * 475)
    Tf1 = f1 * q0 * math.sqrt(ap) * jmax() / (math.sqrt(PI) * e * math.sqrt(u1))
    Tf2 = (1 - f1) * q0 * math.sqrt(ap) * jmax() / (math.sqrt(PI) * e * math.sqrt(u2))
    t("Particao de Blok: T1 = T2", Tf1, Tf2, 1e-9)
    t("Particao de Blok: T1 = Tflash", Tf1, Tf, 1e-9)
    # 9 S-N reproduz os dois pontos do MAIN
    t("S-N: N(800 MPa)", SN[0][1] * (SN[0][0] / 800e6) ** SN_M, 2.75e6, 1e-9)
    t("S-N: N(700 MPa)", SN[0][1] * (SN[0][0] / 700e6) ** SN_M, 6.73e6, 1e-9)
    t("S-N inversa: N=2,75e6 retorna 800 MPa", avaliar_fadiga_SN(750e6, 2.75e6)["pH_adm_N"], 800e6, 1e-9)
    t("S-N inversa em 5e7 e rotulada extrapolacao", 1.0 if avaliar_fadiga_SN(750e6, 5e7)["faixa_N"] == "extrapolacao" else 0.0, 1.0, 1e-12)
    t("S-N inversa em 5e7 pela lei existente", avaliar_fadiga_SN(750e6, 5e7)["pH_adm_N"], 518.9796859942684e6, 1e-8)
    t("S-N sensibilidade: elasticidade dlnNf/dlnpH = -m", avaliar_fadiga_SN(750e6, 5e6)["elasticidade_Nf_pH"], -SN_M, 1e-12)
    t("S-N faixa de pressao identifica extrapolacao independente de N", 1.0 if avaliar_fadiga_SN(476.1e6, 5e6)["faixa_pH"] == "extrapolacao" else 0.0, 1.0, 1e-12)
    t("S-N derivada analitica concorda com diferenca central", avaliar_fadiga_SN(750e6, 5e6)["dNf_dpH"], (avaliar_fadiga_SN(750.001e6, 5e6)["Nf_estimado"]-avaliar_fadiga_SN(749.999e6, 5e6)["Nf_estimado"])/(2e3), 2e-5)
    # 9b conversao de ciclos por dente e media transversal de desgaste
    t("Ciclos/dente pinhao 2/60, 3000 rpm, 1 h", ciclos_por_dente(2, 60, 3000, 1, "pinhao"), 180000, 1e-12)
    t("Ciclos/dente coroa 2/60, 3000 rpm, 1 h", ciclos_por_dente(2, 60, 3000, 1, "coroa"), 6000, 1e-12)
    t("Desgaste medio transversal = 2/3 do maximo", R["wmean1"] / R["w1"], 2 / 3, 1e-12)
    # 10 potencia dissipada: int mu p vs dA = mu F vs
    t("Potencia: (2pi/3) a b p0 mu vs = mu F vs", 2 * PI * R["a"] * R["b"] * R["p0"] / 3 * 0.05 * 3.0, 0.05 * R["F"] * 3.0, 1e-12)
    # 11 filme: Chittenden (continuidade entre eixos e acordo com HD)
    hA = filme_ch(R["F"], 0.02, 0.02, R["Ep"], R["eta0"], R["ue"], P["alp"])[0]
    hB = filme_ch(R["F"], 0.02, 0.02 * (1 + 1e-9), R["Ep"], R["eta0"], R["ue"], P["alp"])[0]
    hC = filme_ch(R["F"], 0.02 * (1 + 1e-9), 0.02, R["Ep"], R["eta0"], R["ue"], P["alp"])[0]
    t("Chittenden: continuidade eixo menor/maior (Rt=Re)", hB, hC, "ch_cont")
    hD = filme_ch(R["F"], R["Rx"], R["Ry"], R["Ep"], R["eta0"], R["ue"], P["alp"])[0]
    t("Chittenden vs Hamrock-Dowson (caso padrao)", hD, h1, "ch_hd")
    # 12 HRH: relacoes de Stadtfeld e cinematica de linha de contato
    A = hrh(0.06, 0.05, 35 * DEG, 20 * DEG, 20 * DEG, 10.0, 1500.0, 1.0, 0.05)
    t("Stadtfeld: tan(beta_p - beta_g) = a/Am", math.tan(A["bp"] - 35 * DEG), 0.05 / 0.06, 1e-12)
    t("Stadtfeld: offset 0 -> beta_p = beta_g", hrh(0.06, 0.0, 35 * DEG, 20 * DEG, 20 * DEG, 10.0, 1500.0, 1.0, 0.05)["bp"], 35 * DEG, 1e-12)
    Fv = math.sqrt(A["Fx"] ** 2 + A["Fy"] ** 2 + A["Fz"] ** 2)
    t("Stadtfeld: |Fx,Fy,Fz| = Ft/(cos a cos b)", Fv, A["Ft"] / (math.cos(20 * DEG) * math.cos(A["bp"])), 1e-12)
    l, vp = A["lam"], A["vp"]
    u1v = (vp * math.cos(l), vp * math.sin(l)); u2v = (-A["vg"] * math.sin(l), A["vg"] * math.cos(l))
    t("Linha de contato: componente normal a linha igual nos dois corpos", u1v[1], u2v[1], 1e-12)
    t("Linha de contato: v_desl = vp/cos(lead)", u1v[0] - u2v[0], A["vs"], 1e-12)
    ex, ey = (u1v[0] + u2v[0]) / 2, (u1v[1] + u2v[1]) / 2
    t("Linha de contato: angulo de arrastamento = 2*lead", math.atan2(ey, ex), 2 * l, 1e-12)
    t("Linha de contato: u_e = v_desl/2", math.hypot(ex, ey), A["vs"] / 2, 1e-12)
    ee = (math.cos(2 * l), math.sin(2 * l))
    t("Projecao de u1 no arrastamento = vp cos(lead)", u1v[0] * ee[0] + u1v[1] * ee[1], A["u1"], 1e-12)
    t("Projecao de u2 no arrastamento = vp sin2/cos", u2v[0] * ee[0] + u2v[1] * ee[1], A["u2"], 1e-12)
    # 13 rendimento: parafuso exato e balanco de potencia (lead grande, mu pequeno)
    B = hrh(0.06, 0.0, 70 * DEG, 20 * DEG, 20 * DEG, 10.0, 1500.0, 1.0, 0.01)
    lb = B["lam"]; ca = math.cos(20 * DEG)
    t("Rendimento: tan(l)/tan(l+rho) = (cos a - mu tan l)/(cos a + mu cot l)", B["eta"], (ca - 0.01 * math.tan(lb)) / (ca + 0.01 / math.tan(lb)), 1e-12)
    t("Rendimento vs balanco de potencia 1 - mu/(cos a sin l cos l)", B["eta"], 1 - 0.01 / (ca * math.sin(lb) * math.cos(lb)), "eta_pot")
    # 14 dureza: Stadtfeld indica 62 HRC no pinhao (HRC->HV de memoria, ASTM E140: ~746 HV)
    t("H=7,4 GPa vs 62 HRC (~746 HV x 9,807e-3)  [conversao de memoria]", 7.4, 746 * 9.807e-3, "H_HRC")
    # 15 Nakagawa et al. (2024), Tabela 1: hipoides de alta reducao 3/60, 6/60, 12/60 (offset 40 mm)
    for k, (z1, ii, Amp, Amg, gp, dg, bpt, bgt, eps, eta, ar) in NAK.items():
        Amp_, Amg_ = Amp * 1e-3, Amg * 1e-3
        t("Nakagawa %s: beta_p Stadtfeld (Am da coroa) vs tabela" % k,
          hrh(Amg_, NAK_A * 1e-3, bgt * DEG, NAK_AL * DEG, dg * DEG, 1.0, 60.0, 1.0, 0.05)["bp"], bpt * DEG, "hrh_bp")
        Hg = hrh(Amg_, 0.0, bgt * DEG, NAK_AL * DEG, dg * DEG, NAK_T, 60.0, 1.0, 0.05)
        Hp = hrh(Amp_, 0.0, bpt * DEG, NAK_AL * DEG, gp * DEG, NAK_T / ii, 60.0, 1.0, 0.05)
        t("Nakagawa %s: Fn pelo pinhao = Fn pela coroa" % k, Hp["Fn"], Hg["Fn"], "hrh_Fn")
        Hv = hrh(Amp_, NAK_A * 1e-3, bgt * DEG, NAK_AL * DEG, gp * DEG, NAK_T / ii, 60.0, 1.0, 0.05, bpt * DEG)
        vg_real = (60.0 * 2 * PI / 60 / ii) * Amg_ * math.sin(dg * DEG)
        t("Nakagawa %s: v_coroa = vp tan(lead) vs cinematica da tabela" % k, Hv["vg"], vg_real, "hrh_vg")
    Hg = hrh(NAK["3/60"][3] * 1e-3, 0.0, NAK["3/60"][7] * DEG, NAK_AL * DEG, NAK["3/60"][5] * DEG, NAK_T, 60.0, 1.0, 0.05)
    t("Nakagawa Eq.(2) F = T/(A_G sin d_G) = Ft de Stadtfeld", Hg["Ft"], NAK_T / (NAK["3/60"][3] * 1e-3 * math.sin(NAK["3/60"][5] * DEG)), 1e-12)
    # 16 Paouris et al. (2016): 9 pontos do ciclo de engrenamento (Tab. 3, Fig. 7 e 8), ponta a ponta por calcular()
    Q0 = dict(PAD_REF); Q0.update(E=210e9, nu=0.33, nu40=0.135 / 890e-6 * 1.0, nu100=0.135 / 890e-6, rho=890.0, alp=1e-8, Tb=40.0,
                              R2x=1e12, R2y=1e12)
    pior = {"p": (0, 0, 0), "h": (0, 0, 0)}; rmin, rmax = 9.0, 0.0
    for i, (W, U, V, Rx, Ry, hcn, hmn, pf) in enumerate(PAO, 1):
        Q = dict(Q0, F=W, R1x=Rx, R1y=Ry, u1=math.hypot(U, V), u2=math.hypot(U, V), th=math.atan2(V, U))
        Rr = calcular(Q, False)
        for key, v, ref in (("p", Rr["p0"] / 1e9, pf), ("h", Rr["hmin"] * 1e6, hmn)):
            e_ = abs(v - ref) / ref
            if e_ >= pior[key][0]: pior[key] = (e_, v, ref, i)
        r_ = 2 * Rr["b"] / 29.999e-3; rmin, rmax = min(rmin, r_), max(rmax, r_)
    t("Paouris: p0 de Hertz nos 9 pontos vs Fig. 7 (pior P%d)" % pior["p"][3], pior["p"][1], pior["p"][2], "pao_p0")
    t("Paouris: eixo maior / largura da coroa, minimo vs 0,16", rmin, 0.16, "pao_geom")
    t("Paouris: eixo maior / largura da coroa, maximo vs 0,38", rmax, 0.38, "pao_geom")
    t("Paouris: hmin nos 9 pontos vs EHL numerico (pior P%d)" % pior["h"][3], pior["h"][1], pior["h"][2], "pao_h")
    t("Hertz de linha: max |tau_zx|/p0 vs 0,25 (sigma_e = 0,5 no artigo)", tau_zx_linha(), 0.25, "tau_lin")
    t("Paouris P5: amplitude de tau_zx/p0 do EHL vs Hertz de linha", (0.227 + 0.2556) / 2, tau_zx_linha(), "pao_tau")
    # 17 integrador eliptico em k=15 contra o limite de linha de Johnson (nu=0,3)
    _z, _x, _y, _s, _v, tak = tensoes(0.3e-3, 4.5e-3, 1e9, 0.3)
    t("Tensao eliptica k=15: tau_max/p0 vs limite de linha 0,30", max(tak) / 1e9, 0.30, "stress_k15")
    # 18 calibracao com divisao calibra/valida (MAIN 27.20): K_BP no Nakagawa -> valida na geometria do Paouris; C_h nos pontos impares -> valida nos pares
    kbp, ch = calibracoes(); G = PAO_GEO
    bpp = hrh(G["Amg"], G["a"], G["bg"] * DEG, NAK_AL * DEG, G["dg"] * DEG, 1.0, 60.0, 1.0, 0.05)["bp"]
    t("Paouris geometria: beta_p de Stadtfeld sem calibrar vs 45,989 graus", bpp, G["bp"] * DEG, "hrh_bp")
    t("Paouris geometria: beta_p com K=%.4f calibrado no Nakagawa" % kbp, bpp * kbp, G["bp"] * DEG, "bp_cal")
    pior_s = (0, 0, 0, 0)
    for idx in (1, 3, 5, 7):
        rw = PAO[idx]; v_ = ch * _pao_modelo(rw)["hmin"] * 1e6; e_ = abs(v_ - rw[6]) / rw[6]
        if e_ >= pior_s[0]: pior_s = (e_, v_, rw[6], idx + 1)
    t("Paouris split: hmin dos pares com fator %.3f ajustado nos impares (pior P%d)" % (ch, pior_s[3]), pior_s[1], pior_s[2], "pao_split")
    V_IDX = {1, 2, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 50, 57}
    R_IDX = {3, 4, 5, 6, 7, 11, 21, 27, 28, 40}
    for i, x in enumerate(L, 1):
        L[i - 1] = x[:5] + ("V" if i in V_IDX else ("R" if i in R_IDX else "E"),) + x[6:]
    return L

# ---------------------------------------------------------------- v10: VALIDACAO DO HISTORICO DE TENSOES
def _linha_johnson(a, p0, x, z):
    """Contato de linha de Hertz (Johnson, Contact Mechanics, eq. 4.45), deformacao plana: (sx, sz, txz) em compressao positiva."""
    A = a * a - x * x + z * z; D = math.sqrt(A * A + 4.0 * x * x * z * z)
    m = math.sqrt(0.5 * (D + A)); nn = math.copysign(math.sqrt(max(0.5 * (D - A), 0.0)), x if x != 0 else 1.0)
    k = p0 / a; mn = m * m + nn * nn
    return (k * (m * (1 + (z * z + nn * nn) / mn) - 2 * z), k * m * (1 - (z * z + nn * nn) / mn), k * nn * (m * m - z * z) / mn)

def _lcg(seed):
    s = [seed]
    def f():
        s[0] = (1103515245 * s[0] + 12345) % 2147483648
        return s[0] / 2147483648.0 - 0.5
    return f

def validar_historico(rapido=False):
    """Testes do historico de tensoes. Cada item: (nome, valor, referencia, erro, ok, nivel, tol). Niveis: V identidade/solucao fechada, R referencia analitica publicada."""
    L = []
    def t(nome, val, ref, tol, nivel="V", rel=True, esc=None):
        tl = TOL[tol] if isinstance(tol, str) else tol
        err = abs(val - ref) / abs(ref) if rel else abs(val - ref) / (esc if esc else 1.0)
        L.append((nome, val, ref, err, err <= tl, nivel, tl))
    nu = 0.29
    # 1 nucleo de Boussinesq: equilibrio local, traco e simetrias
    p = (0.7, -0.4, 0.9); h = 1e-5
    def d(k, ax):
        pp = [list(p), list(p)]; pp[0][ax] += h; pp[1][ax] -= h
        return (boussinesq(*pp[0], nu)[k] - boussinesq(*pp[1], nu)[k]) / (2 * h)
    SX, SY, SZ, TXY, TXZ, TYZ = range(6)
    esc = abs(d(SZ, 2)) + abs(d(SX, 0))
    for nome, e in (("x", d(SX, 0) + d(TXY, 1) + d(TXZ, 2)), ("y", d(TXY, 0) + d(SY, 1) + d(TYZ, 2)), ("z", d(TXZ, 0) + d(TYZ, 1) + d(SZ, 2))):
        t("Boussinesq: equilibrio local, direcao %s (residuo/gradiente)" % nome, e / esc, 0.0, "h_pt", rel=False)
    bs = boussinesq(*p, nu); Rr = math.sqrt(p[0] ** 2 + p[1] ** 2 + p[2] ** 2)
    t("Boussinesq: traco sx+sy+sz = (1+nu) z/(pi R^3)", bs[0] + bs[1] + bs[2], (1 + nu) * p[2] / (PI * Rr ** 3), "h_id")
    t("Boussinesq: sz = 3 z^3/(2 pi R^5) (compressao +)", bs[2], 3 * p[2] ** 3 / (2 * PI * Rr ** 5), "h_id")
    # 2 invariantes e principais
    rg = _lcg(7); worst = 0.0; worst_vm = 0.0
    for _ in range(200):
        c = [rg() * 2e9 for _ in range(6)]; e = principais(*c)
        I1 = c[0] + c[1] + c[2]; I2 = c[0] * c[1] + c[1] * c[2] + c[2] * c[0] - c[3] ** 2 - c[4] ** 2 - c[5] ** 2
        I3 = c[0] * (c[1] * c[2] - c[5] ** 2) - c[3] * (c[3] * c[2] - c[5] * c[4]) + c[4] * (c[3] * c[5] - c[1] * c[4])
        sc = max(abs(v) for v in e) + 1e-30
        worst = max(worst, abs(sum(e) - I1) / sc, abs(e[0] * e[1] + e[1] * e[2] + e[2] * e[0] - I2) / sc ** 2, abs(e[0] * e[1] * e[2] - I3) / sc ** 3)
        vme = math.sqrt(0.5 * ((e[0] - e[1]) ** 2 + (e[1] - e[2]) ** 2 + (e[2] - e[0]) ** 2))
        worst_vm = max(worst_vm, abs(vme - vm_comp(*c)) / sc)
    t("Principais: erro maximo de I1, I2, I3 em 200 tensores aleatorios", worst, 0.0, "h_id", rel=False)
    t("von Mises dos autovalores = von Mises das 6 componentes (200 tensores)", worst_vm, 0.0, "h_id", rel=False)
    # 3 circular no eixo vs Johnson, e superficie externa exata
    a1, p1 = 1e-3, 1e9; w = 0.0; wv = 0.0
    for zr in (0.25, 0.5, 1.0, 2.0):
        z = zr * a1; T = campo_tensor(a1, a1, p1, nu, 0.0, 0.0, z)
        szj = p1 / (1 + zr * zr); srj = p1 * ((1 + nu) * (1 - zr * math.atan(1 / zr)) - 1 / (2 * (1 + zr * zr)))
        w = max(w, abs(T[2] - szj) / p1, abs(T[0] - srj) / p1, abs(T[1] - srj) / p1)
        wv = max(wv, abs(vm_comp(*T) - abs(szj - srj)) / p1)
    t("Circular, eixo: sz, sx, sy vs Johnson (4 profundidades)", w, 0.0, "h_circ", "R", rel=False)
    t("Circular, eixo: von Mises (6 componentes) vs |sz - sr| de Johnson", wv, 0.0, "h_circ", "R", rel=False)
    for rr in (2.0, 3.0):
        T = campo_tensor(a1, a1, p1, nu, rr * a1, 0.0, 0.005 * a1)
        ref = (1 - 2 * nu) * p1 / (3 * rr * rr)
        t("Circular, superficie externa r=%ga: radial = -(1-2nu)p0/(3 (r/a)^2) (tracao)" % rr, T[0], -ref, "h_ext", "R")
        t("Circular, superficie externa r=%ga: circunferencial = +(1-2nu)p0/(3 (r/a)^2)" % rr, T[1], ref, "h_ext", "R")
    # 4 elipse k=7: superficie livre e simetrias
    a2, b2, p2 = 1e-3, 7e-3, 1e9; c2 = a2
    w = 0.0
    for (x, y) in ((0.0, 0.0), (0.3 * a2, 0.2 * b2), (0.6 * a2, -0.4 * b2)):
        T = campo_tensor(a2, b2, p2, nu, x, y, 0.01 * c2)
        pp = p2 * math.sqrt(max(0.0, 1 - (x / a2) ** 2 - (y / b2) ** 2))
        w = max(w, abs(T[2] - pp) / p2, abs(T[4]) / p2, abs(T[5]) / p2)
    t("Elipse k=7, z=0,01c: sz = p(x,y), txz = tyz = 0 (3 pontos)", w, 0.0, "h_surf", rel=False)
    x, y, z = 0.37 * a2, 0.21 * b2, 0.4 * c2; T0 = campo_tensor(a2, b2, p2, nu, x, y, z)
    par = {"x": (1, 1, 1, -1, -1, 1), "y": (1, 1, 1, -1, 1, -1)}      # paridade de (sx,sy,sz,txy,txz,tyz) sob x->-x e y->-y
    w = 0.0
    for eixo, (xx, yy) in (("x", (-x, y)), ("y", (x, -y))):
        T = campo_tensor(a2, b2, p2, nu, xx, yy, z)
        w = max(w, max(abs(T[i] - par[eixo][i] * T0[i]) for i in range(6)) / p2)
    Tb = campo_tensor(a2, b2, p2, nu, -x, -y, z)
    w = max(w, max(abs(Tb[i] - par["x"][i] * par["y"][i] * T0[i]) for i in range(6)) / p2)
    t("Elipse k=7: paridades sx,sy,sz,tyz pares em x; txy,txz impares em x (e analogo em y)", w, 0.0, "h_sim", rel=False)
    # 5 equilibrio local do campo numerico (dif. finitas), 3 pontos
    hh = 0.04 * c2; w = 0.0
    for (x, y, z) in ((0.4 * a2, 0.5 * a2, 0.5 * c2), (1.0 * a2, 0.0, 0.8 * c2), (2.5 * a2, 1.0 * a2, 0.6 * c2)):
        def F(dx, dy, dz): return campo_tensor(a2, b2, p2, nu, x + dx, y + dy, z + dz)
        g = {}
        for ax, (dx, dy, dz) in enumerate(((hh, 0, 0), (0, hh, 0), (0, 0, hh))):
            fp, fm = F(dx, dy, dz), F(-dx, -dy, -dz); g[ax] = [(fp[i] - fm[i]) / (2 * hh) for i in range(6)]
        e1 = g[0][SX] + g[1][TXY] + g[2][TXZ]; e2 = g[0][TXY] + g[1][SY] + g[2][TYZ]; e3 = g[0][TXZ] + g[1][TYZ] + g[2][SZ]
        gm = max(abs(v) for ax in g for v in g[ax]); w = max(w, abs(e1) / gm, abs(e2) / gm, abs(e3) / gm)
    t("Elipse k=7: equilibrio local d(sij)/dxj = 0 por dif. finitas (3 pontos, residuo/gradiente)", w, 0.0, "h_eq", rel=False)
    # 6 equilibrio de forca: integral de sz no plano z = F (circular e elipse k=3)
    for (aa, bb, nome) in ((1e-3, 1e-3, "circular"), (1e-3, 3e-3, "elipse k=3")):
        F0 = 2 * PI * aa * bb * p2 / 3; z = max(aa, bb); Rm = 15 * bb; nr, nph = 15, (1 if aa == bb else 10)
        xg, wg = _gl(5); tot = 0.0
        for jp in range(nph):
            ph = (jp + 0.5) * (PI / 2) / nph if nph > 1 else 0.0
            for ip in range(nr):
                for xk, wk in zip(xg, wg):
                    r = (ip + 0.5 * (xk + 1)) * Rm / nr
                    sz_ = campo_tensor(aa, bb, p2, nu, r * math.cos(ph), r * math.sin(ph), z)[2]
                    tot += sz_ * r * (Rm / nr) * 0.5 * wk * (2 * PI / nph)
        t("Equilibrio de forca %s: integral de sz dxdy em z=%.0f mm = F" % (nome, z * 1e3), tot, F0, "h_eqF")
    # 7 compatibilidade com o tensoes() do v9 (eixo z), causa da divergencia e limite fisico
    zs, sx_, sy_, sz_, sv_, ta_ = tensoes(a2, b2, p2, nu)
    z9 = tensoes_v9(a2, b2, p2, nu); ws = 0.0
    def dif(Z9, idx):
        return max(max(abs(campo_tensor(a2, b2, p2, nu, 0.0, 0.0, zs[i])[k] - Z9[1 + k][i]) / p2 for k in range(3)) for i in idx)
    for i in range(0, 40, 4):
        T = campo_tensor(a2, b2, p2, nu, 0.0, 0.0, zs[i]); ws = max(ws, abs(T[3]) / p2, abs(T[4]) / p2, abs(T[5]) / p2)
    t("Eixo x=y=0: txy = txz = tyz = 0", ws, 0.0, "h_sim", rel=False)
    prof = [i for i in range(40) if zs[i] >= 0.5 * c2]
    t("Compatibilidade com o v9 (60x16) onde ele convergiu (z >= 0,5c): sx, sy, sz", dif(z9, prof), 0.0, "h_axis", rel=False)
    idx10 = list(range(0, 40, 4)); e1 = dif(z9, idx10)
    e2 = dif(tensoes_v9(a2, b2, p2, nu, 40, 120, 64), idx10); e3 = dif(tensoes_v9(a2, b2, p2, nu, 40, 240, 256), idx10)
    t("Causa da divergencia do v9 em z < 0,5c: o erro cai com a quadratura (60x16 > 120x64 > 240x256) [1 = sim]", 1.0 if e1 > e2 > e3 else 0.0, 1.0, 1e-12)
    t("v9 com quadratura 240x256 vs campo completo no eixo (10 profundidades)", e3, 0.0, "h_axis", rel=False)
    mz = max(sz_) / p2
    L.append(("Limite fisico: sz no eixo <= p0 em todas as profundidades (tensoes corrigido; v9 dava %.4f p0 em z=0,1c)" % (max(z9[3]) / p2),
              mz, 1.0, max(0.0, mz - 1.0), mz <= 1.0 + 1e-9, "V", 1e-9))
    # 8 maximo global de von Mises (MAIN 9.10) vs maximo no eixo (o que tensoes() calcula)
    ii = list(range(6, 40, 3)); gmax = 0.0; pos = None
    for xk in (0.0, 0.5, 1.0, 1.5):
        for yk in (0.0, 0.1, 0.25, 0.5, 0.75):
            for i in ii:
                v = vm_comp(*campo_tensor(a2, b2, p2, nu, xk * a2, yk * b2, zs[i]))
                if v > gmax: gmax, pos = v, (xk, yk, zs[i] / c2)
    amax = max(sv_[i] for i in ii)
    t("Maximo de von Mises na grade x,y,z (MAIN 9.10) vs maximo no eixo de tensoes() (pos. x/a=%g, y/b=%g, z/c=%.2f)" % pos, gmax, amax, "h_vmax")
    # 9 convergencia: malha, amostragem no tempo e janela
    Rk = dict(a=a2, b=b2, p0=p2, nu=nu, th=0.0, u1=8.0, u2=5.0); zc = 0.5 * c2
    def met(**kw):
        r = resumo_hist(historico(Rk, zc, **kw)); return r["svm_max"], r["tsz_amp"]
    base, fino = met(), met(n=8, q=2.0)
    t("Convergencia de malha (n=5,q=3 vs n=8,q=2): von Mises max do ciclo", base[0], fino[0], "h_conv")
    t("Convergencia de malha (n=5,q=3 vs n=8,q=2): amplitude de tau ortogonal", base[1], fino[1], "h_conv")
    ee = [abs(met(n=n_, q=q_)[1] - fino[1]) / fino[1] for (n_, q_) in ((4, 3.0), (5, 3.0), (6, 2.0))]
    t("Convergencia monotona do erro (n,q)=(4,3) >= (5,3) >= (6,2) [1 = sim]", 1.0 if ee[0] >= ee[1] >= ee[2] else 0.0, 1.0, 1e-12)
    n161 = met(nx=161); n41 = met(nx=41)
    t("Amostragem no tempo com extremos refinados (nx=41 vs 161): amplitude de tau ortogonal", n41[1], n161[1], "h_conv")
    t("Amostragem no tempo (nx=81 vs 161): von Mises max", base[0], n161[0], "h_conv")
    t("Amostragem no tempo (nx=81 vs 161): amplitude de tau ortogonal", base[1], n161[1], "h_conv")
    j12 = met(jan=12.0)
    t("Janela do ciclo (6 vs 12 semi-cordas): von Mises max", base[0], j12[0], "h_conv")
    t("Janela do ciclo (6 vs 12 semi-cordas): amplitude de tau ortogonal", base[1], j12[1], "h_conv")
    # 10 limite de linha (k=100) vs Johnson 4.45, e a transcricao de Johnson satisfaz equilibrio
    al, bl = 1e-3, 100e-3; w = 0.0
    for xk in (0.0, 0.5, 1.0, 1.5, 2.0):
        for zk in (0.25, 0.5, 1.0):
            T = campo_tensor(al, bl, p2, nu, xk * al, 0.0, zk * al); J = _linha_johnson(al, p2, xk * al, zk * al)
            w = max(w, abs(T[0] - J[0]) / p2, abs(T[2] - J[1]) / p2, abs(T[4] - J[2]) / p2)
    t("Limite de linha k=100: sx, sz, txz vs Johnson 4.45 (15 pontos)", w, 0.0, "h_line", "R", rel=False)
    hj = 1e-6 * al; xj, zj = 0.6 * al, 0.7 * al
    f = lambda xx, zz: _linha_johnson(al, p2, xx, zz)
    dsx = (f(xj + hj, zj)[0] - f(xj - hj, zj)[0]) / (2 * hj); dtx = (f(xj + hj, zj)[2] - f(xj - hj, zj)[2]) / (2 * hj)
    dtz = (f(xj, zj + hj)[2] - f(xj, zj - hj)[2]) / (2 * hj); dsz = (f(xj, zj + hj)[1] - f(xj, zj - hj)[1]) / (2 * hj)
    gj = max(abs(dsx), abs(dtx), abs(dtz), abs(dsz))
    t("Johnson 4.45 transcrito: equilibrio d sx/dx + d txz/dz = 0 (residuo/gradiente)", (dsx + dtz) / gj, 0.0, "h_pt", rel=False)
    t("Johnson 4.45 transcrito: equilibrio d txz/dx + d sz/dz = 0 (residuo/gradiente)", (dtx + dsz) / gj, 0.0, "h_pt", rel=False)
    # 11 amplitude de tau ortogonal no limite de linha vs tau_zx_linha() existente
    Rl = dict(a=al, b=bl, p0=p2, nu=nu, th=0.0, u1=8.0, u2=5.0)
    V, ic, iv = varrer_z(Rl, 0.0, zs_hist(Rl, 25))
    t("Limite de linha: amplitude de tau ortogonal / p0 (historico, k=100) vs tau_zx_linha() existente", V[ic]["tsz_amp"] / p2, tau_zx_linha(), "h_line")
    # 12 troca de eixos: arrastamento em y numa elipse (a,b) = arrastamento em x na elipse (b,a)
    Ra = dict(a=a2, b=b2, p0=p2, nu=nu, th=PI / 2, u1=8.0, u2=5.0); Rb = dict(a=b2, b=a2, p0=p2, nu=nu, th=0.0, u1=8.0, u2=5.0)
    A = resumo_hist(historico(Ra, zc)); B = resumo_hist(historico(Rb, zc))
    t("Troca de eixos (th=90 graus vs elipse trocada): von Mises max", A["svm_max"], B["svm_max"], "h_dir")
    t("Troca de eixos (th=90 graus vs elipse trocada): amplitude de tau ortogonal", A["tsz_amp"], B["tsz_amp"], "h_dir")
    # 13 fadiga (MAIN 13.22) e eixo de tempo
    Rc = calcular(PAD_REF, hist=True)
    t("Miner de um nivel = D = N/Nf do calcular()", dano_miner([(Rc["N"], Rc["p0"])]), Rc["D"], "h_id")
    n1_, n2_, pa_, pb_ = 1e5, 2e6, 750e6, 700e6
    ref2 = n1_ / (SN[0][1] * (SN[0][0] / pa_) ** SN_M) + n2_ / (SN[0][1] * (SN[0][0] / pb_) ** SN_M)
    t("Miner de dois niveis = n1/N1 + n2/N2 pela S-N do MAIN", dano_miner([(n1_, pa_), (n2_, pb_)]), ref2, "h_id")
    t("Miner: 6,73e6 ciclos a 700 MPa (ponto do MAIN) da D = 1", dano_miner([(6.73e6, 700e6)]), 1.0, "h_id")
    H = historico(Rk, zc, nx=21)
    t("Eixo de tempo: (t_fim - t_ini) u1 = 2L", (H["t1"][-1] - H["t1"][0]) * 8.0, 2 * H["L"], "h_id")
    # 14 S-N reescrita em Phi (amplitude de tau ortogonal), contato eliptico
    t("S-N em Phi, contato nominal: D_phi = D em pH (por construcao, um nivel)", Rc["D_phi"], Rc["D"], "h_id")
    t("S-N em Phi: p_eq = pH no contato de calibracao", Rc["p_eq"], Rc["p0"], "h_id")
    fh = lambda RR: fadiga_historico(dict(RR, N=Rc["N"]), nz=10)
    f1 = fh(dict(Rk, F=1.0)); f2 = fh(dict(Rk, a=2 * a2, b=2 * b2, p0=2 * p2))
    t("Phi/p0 invariante a escala do contato (a,b,p0) -> (2a,2b,2p0), mesmo k", f2["rho_ref"], f1["rho_ref"], "h_phi")
    rho = f1["rho_ref"]; Phi1, Phi2 = 0.8 * f1["Phi"], 1.1 * f1["Phi"]; n1_, n2_ = 3e5, 4e5
    ref2 = n1_ / avaliar_fadiga_SN(Phi1 / rho, n1_)["Nf_estimado"] + n2_ / avaliar_fadiga_SN(Phi2 / rho, n2_)["Nf_estimado"]
    t("Miner em Phi, dois niveis = n1/Nf(Phi1/rho) + n2/Nf(Phi2/rho)", dano_miner_phi([(n1_, Phi1), (n2_, Phi2)], rho), ref2, "h_id")
    fl = fadiga_historico(dict(Rl, N=Rc["N"]), nz=10, rho_ref=rho)
    t("Contato de referencia dado: p_eq/pH = rho_contato/rho_ref (k=100 contra k=7)", fl["p_eq"] / p2, fl["Phi"] / p2 / rho, "h_phi")
    return L

def validar_v11():
    """Testes v11 (atrito misto, cota estacionaria, Phi(z) com hidrostatica, D(z), contato em grade e desgaste incremental). Mesmo formato de validar_historico."""
    L = []
    def t(nome, val, ref, tol, nivel="V", rel=True):
        tl = TOL[tol] if isinstance(tol, str) else tol
        err = abs(val - ref) / abs(ref) if rel else abs(val - ref)
        L.append((nome, val, ref, err, err <= tl, nivel, tl))
    def tmenor(nome, val, lim):
        L.append((nome, val, lim, max(0.0, val - lim), val < lim, "V", 0.0))
    # 1 nucleo de influencia e contato em grade
    c = 1.7e-3
    t("Nucleo: int dA/r no quadrado 2c x 2c = 8 c asinh(1)", _rect(-c, c, -c, c), 8 * c * math.asinh(1.0), "h_id")
    x1, x2, y1, y2 = 0.5e-3, 1.5e-3, 0.2e-3, 1.0e-3; nn = 400; hx, hy = (x2 - x1) / nn, (y2 - y1) / nn
    bf = sum(1.0 / math.hypot(x1 + (i + .5) * hx, y1 + (j + .5) * hy) for i in range(nn) for j in range(nn)) * hx * hy
    t("Nucleo: retangulo deslocado vs quadratura de ponto medio 400x400", _rect(x1, x2, y1, y2), bf, "v_k")
    P = dict(PAD_REF); R = calcular(P, False); n = 10; dx, dy = 1.35 * R["a"] / n, 1.35 * R["b"] / n
    A = matriz_influencia(n, n, dx, dy, R["Es"])
    g = [((i + .5) * dx) ** 2 / (2 * R["Rx"]) + ((k + .5) * dy) ** 2 / (2 * R["Ry"]) for i in range(n) for k in range(n)]
    S0 = [i * n + k for i in range(n) for k in range(n) if ((i + .5) * dx / R["a"]) ** 2 + ((k + .5) * dy / R["b"]) ** 2 < 1.0]
    p, d, it = contato_grade(g, A, 4 * dx * dy, R["F"], S0)
    t("Grade n=10, sem desgaste: p_max vs p0 de Hertz exato", max(p), R["p0"], "v_p0", "R")
    t("Grade n=10, sem desgaste: delta vs u(0,0) de Hertz (_u_hertz)", d, _u_hertz(0, 0, R["a"], R["b"], R["p0"], R["Es"]), "v_dl", "R")
    t("Grade: equilibrio 4 dx dy sum p = F", 4 * dx * dy * sum(p), R["F"], "v_eq")
    t("Grade: p >= 0 em todas as celulas (violacao/p0)", max(0.0, -min(p)) / R["p0"], 0.0, 0.0, rel=False)
    u = [sum(A[i][j] * p[j] for j in range(n * n)) for i in range(n * n)]
    t("Grade: compatibilidade u + g = delta nas celulas em contato (residuo/delta)", max(abs(u[i] + g[i] - d) for i in range(n * n) if p[i] > 0) / d, 0.0, "v_eq", rel=False)
    t("Grade: sem penetracao fora do contato (violacao/delta)", max(0.0, -min(u[i] + g[i] - d for i in range(n * n) if p[i] == 0.0)) / d, 0.0, "v_eq", rel=False)
    # 2 desgaste incremental
    I0 = desgaste_incremental(dict(R, K=0.0))
    t("K = 0: w = 0", I0["w0"], 0.0, 0.0, rel=False)
    t("K = 0: p0eq inalterada", I0["hist"][-1]["p0eq"], I0["hist"][0]["p0eq"], "h_id")
    ws = R["w1"] + R["w2"]; Np = 0.02 * d / ws
    I1 = desgaste_incremental(dict(R, N=Np))
    t("Desgaste pequeno (2% do delta): w(0) incremental vs N (w1+w2) linear [Archard integrado]", I1["w0"], Np * ws, "v_w5")
    Nm = 100e-6 / ws; Rm = dict(R, N=Nm)
    Ia = desgaste_incremental(Rm, n=10, tol_p=0.35); Ib = desgaste_incremental(Rm, n=12, tol_p=0.2)
    t("Convergencia (n=10, tol 0,35 vs n=12, tol 0,2): w(0) apos N com w linear = 100 um", Ia["w0"], Ib["w0"], "v_conv")
    tmenor("Realimentacao: w(0) incremental / w linear estatico < 1 (conformidade desacelera o desgaste)", Ia["razao_w"], 1.0)
    tmenor("Realimentacao: p0eq final / p0eq inicial < 1", Ia["hist"][-1]["p0eq"] / Ia["hist"][0]["p0eq"], 1.0)
    mono = all(b["w0"] >= a_["w0"] for a_, b in zip(Ia["hist"][:-1], Ia["hist"][1:]))
    t("w(0) nao decresce com N [1 = sim]", 1.0 if mono else 0.0, 1.0, 1e-12)
    t("Passos terminam no N pedido [1 = sim]", 1.0 if Ia["completo"] else 0.0, 1.0, 1e-12)
    # 3 atrito misto
    t("alfa_asp(0) = 1", alfa_asp(0.0), 1.0, "h_id")
    t("E[(X)+^(3/2)] = 2^(1/4) Gamma(5/4)/sqrt(2 pi)", _f32(0.0), 2 ** 0.25 * math.gamma(1.25) / math.sqrt(2 * PI), "v_f32")
    t("alfa_asp decresce com Lambda (0 a 4, 40 passos) [1 = sim]", 1.0 if all(alfa_asp(0.1 * i) > alfa_asp(0.1 * (i + 1)) for i in range(40)) else 0.0, 1.0, 1e-12)
    tmenor("alfa_asp(4) < 1e-3", alfa_asp(4.0), 1e-3)
    t("Eyring: limite newtoniano tau = eta gamma_ponto", tau_eyring(0.0, 1e-6, 1.0, 5e-6, 0.0, 5e6), 5.0, "h_id")
    t("Eyring: alto cisalhamento tau = tau0 ln(2x)", tau_eyring(0.0, 1e-6, 1.0, 5e6, 0.0, 5e6), 5e6 * math.log(2e6), "h_id")
    t("mu dinamico com Lambda = 0 (alfa = 1) = mu_asp", atrito_misto(dict(R, Lam=0.0), P)["mu_din"], P["mua"], "h_id")
    am = atrito_misto(dict(R, Lam=60.0), P)
    t("mu dinamico com Lambda = 60 (alfa = 0) = mu do fluido", am["mu_din"], am["mu_fluido"], "h_id")
    a_, b_ = R["a"], R["b"]; hf = am["hf"]; nb = 200; tot = 0.0
    for i in range(nb):
        for j in range(nb):
            xx, yy = -a_ + (i + .5) * 2 * a_ / nb, -b_ + (j + .5) * 2 * b_ / nb; r2 = (xx / a_) ** 2 + (yy / b_) ** 2
            if r2 < 1.0: tot += tau_eyring(R["p0"] * math.sqrt(1.0 - r2), hf, R["vs"], R["eta0"], P["alp"], P["tau0"])
    tot *= (2 * a_ / nb) * (2 * b_ / nb)
    t("int tau_f dA (s = sqrt(1-rho^2), Simpson 1D) vs grade cartesiana 200x200 na elipse", am["mu_fluido"] * R["F"], tot, "v_fint", "R")
    # 4 cota estacionaria
    t("G circulo = pi a / 4", G_estac(1e-3, 1e-3), PI * 1e-3 / 4, "h_id")
    t("G(2a,2b) = 2 G(a,b)", G_estac(2e-3, 6e-3), 2 * G_estac(1e-3, 3e-3), "h_id")
    t("G(a,b) = G(b,a)", G_estac(1e-3, 3e-3), G_estac(3e-3, 1e-3), "h_id")
    ae, be = 1e-3, 3e-3; nb = 200; tot = 0.0
    for i in range(nb):
        for j in range(nb):
            xx, yy = -ae + (i + .5) * 2 * ae / nb, -be + (j + .5) * 2 * be / nb; r2 = (xx / ae) ** 2 + (yy / be) ** 2
            if r2 < 1.0: tot += math.sqrt(1.0 - r2) / math.hypot(xx, yy)
    tot *= (2 * ae / nb) * (2 * be / nb) / (2 * PI)
    t("G elipse k=3 vs integral dupla cartesiana 200x200 (1/2pi) int sqrt(1-rho^2)/r dA", G_estac(ae, be), tot, "v_G", "R")
    t("flash_estac circular = mu p0 vs (pi a/4)/(2 k)", flash_estac(0.05, 1e9, 3.0, 1e-3, 1e-3, 46.6), 0.05 * 1e9 * 3.0 * (PI * 1e-3 / 4) / (2 * 46.6), "h_id")
    # 5 fadiga multiaxial: Phi(z) = tau_amp + lambda p_hid,max e D(z)
    nu = 0.29; a2, b2, p2 = 1e-3, 7e-3, 1e9
    Rk = dict(a=a2, b=b2, p0=p2, nu=nu, th=0.0, u1=8.0, u2=5.0, N=1e7)
    f0 = fadiga_historico(dict(Rk, lam=0.0), nz=8); f3 = fadiga_historico(dict(Rk, lam=0.3), nz=8); V = f3["V"]
    t("lambda = 0: Phi = max_z amplitude de tau ortogonal (v10)", f0["Phi"], max(v["tsz_amp"] for v in f0["V"]), "h_id")
    t("Phi(z_i) = tau_amp + lambda p_hid,max em 8 profundidades (erro/Phi)", max(abs(f3["Phiz"][i] - (V[i]["tsz_amp"] + 0.3 * V[i]["phyd_max"])) for i in range(8)) / f3["Phi"], 0.0, "h_id", rel=False)
    a1, p1, zr = 1e-3, 1e9, 0.5; T = campo_tensor(a1, a1, p1, nu, 0.0, 0.0, zr * a1)
    szj = p1 / (1 + zr * zr); srj = p1 * ((1 + nu) * (1 - zr * math.atan(1 / zr)) - 1 / (2 * (1 + zr * zr)))
    t("Hidrostatica circular no eixo (z=0,5a): (sx+sy+sz)/3 vs Johnson (erro/p0)", abs((T[0] + T[1] + T[2]) / 3 - (szj + 2 * srj) / 3) / p1, 0.0, "h_circ", "R", rel=False)
    t("D(z) maximo = D_phi", max(f3["Dz"]), f3["D_phi"], "h_id")
    t("D na profundidade critica = D em pH (calibrado no nominal)", f3["Dz"][f3["i_tau"]], f3["D_main"], "h_id")
    t("D(z) >= 0 em todas as profundidades (menor valor)", min(0.0, min(f3["Dz"])), 0.0, 0.0, rel=False)
    f4 = fadiga_historico(dict(Rk, a=2 * a2, b=2 * b2, p0=2 * p2, lam=0.3), nz=8)
    t("Phi max escala com p0 sob (a,b,p0) -> (2a,2b,2p0), lambda = 0,3", f4["Phi"], 2 * f3["Phi"], "v_Phi")
    return L

def comparar_nakagawa():
    """Comparacoes SEM aprovacao: o modelo atual nao reproduz estes dados."""
    out = ["COMPARACOES SEM APROVACAO"]
    Q0 = dict(PAD_REF); Q0.update(E=210e9, nu=0.33, nu40=0.135 / 890e-6, nu100=0.135 / 890e-6, rho=890.0, alp=1e-8, Tb=40.0, R2x=1e12, R2y=1e12)
    e0, e1, hcs = [], [], []
    for (W, U, V, Rx, Ry, hcn, hmn, pf) in PAO:
        m = math.hypot(U, V)
        Rr = calcular(dict(Q0, F=W, R1x=Rx, R1y=Ry, u1=U, u2=U), False)          # arrastamento so em x, velocidade U
        e0.append(100 * (Rr["hmin"] * 1e6 - hmn) / hmn)
        tt = math.atan2(V, U)
        Re = 1 / (math.cos(tt) ** 2 / Rx + math.sin(tt) ** 2 / Ry); Rs = 1 / (math.sin(tt) ** 2 / Rx + math.cos(tt) ** 2 / Ry)
        hcs.append(100 * (filme_ch(W, Re, Rs, Rr["Ep"], 0.135, m, 1e-8)[1] * 1e6 - hcn) / hcn)
    out.append("Paouris: ignorar o angulo (arrastamento so em x) superestima hmin em %.0f a %.0f%%; por isso existe o campo de angulo." % (min(e0), max(e0)))
    out.append("Paouris: filme central hc por Chittenden angular erra %.0f a %.0f%% vs EHL numerico: hc nao validado no ramo angular." % (min(hcs), max(hcs)))
    out.append("Paouris: pontos 1 e 9 (carga baixa) tem amplitude de tau_zx 0,216 e 0,188 contra 0,25 de Hertz: Hertz superestima 14 e 25%.")
    out.append("Mohammadpour 2017: Hertz a meio do ciclo bate (1,14 GPa vs ~1,2 do grafico), mas a Tabela 4 (U=11,64 m/s) da hmin=1,36 um contra ~0,2-0,3 um do artigo: nao reconciliei, nao usei.")
    A = AGRACIA
    lam = PI / 2 - A["bp0"] * DEG; dgs = 79 * DEG
    kin = (A["dp_g"] / (A["z2"] / A["z1"])) / A["dp_p"]
    out.append("A.Gracia 2/60: fechamento cinematico v_coroa/v_pinhao = (dp_g/razao)/dp_p = %.4f vs tan(avanco 11 graus) = %.4f (%+.1f%%); Nakagawa deu +5 a +10%%." % (kin, math.tan(lam), 100 * (kin - math.tan(lam)) / math.tan(lam)))
    out.append("A.Gracia 2/60: raio primitivo da coroa %.0f mm <= geratriz %.0f mm (coerente: r = A sen delta, delta = %.1f graus)." % (A["dp_g"] / 2, A["Ae"], math.asin(A["dp_g"] / 2 / A["Ae"]) / DEG))
    for nome, amg in (("Am_coroa = 0,857*Ae = %.1f mm" % (AM_G_FRAC * A["Ae"]), AM_G_FRAC * A["Ae"] * 1e-3),
                      ("Am_coroa = (dp_g/2)/sen(delta) = %.1f mm" % (A["dp_g"] / 2 / math.sin(dgs)), A["dp_g"] / 2e3 / math.sin(dgs))):
        bpe = hrh(amg, A["a"] * 1e-3, A["bg"] * DEG, A["al"] * DEG, A["gm"] * DEG, 1.0, 60.0, 1.0, 0.05)["bp"] / DEG
        out.append("A.Gracia 2/60: beta do pinhao pela relacao de Stadtfeld com %s = %.2f vs 79 medido (%+.1f%%)." % (nome, bpe, 100 * (bpe - A["bp0"]) / A["bp0"]))
    out.append("Nenhum dos artigos traz velocidade de deslizamento: flash de Blok segue sem validacao de temperatura de hipoide.")
    for k, (z1, ii, Amp, Amg, gp, dg, bpt, bgt, eps, eta, ar) in NAK.items():
        lam = (90 - bpt) * DEG
        e5 = hrh(Amp * 1e-3, 0.0, bpt * DEG, NAK_AL * DEG, gp * DEG, 1.0, 60.0, 1.0, 0.05)["eta"]
        mui = math.tan(math.atan(math.tan(lam) / (eta / 100)) - lam) * math.cos(NAK_AL * DEG)
        out.append("%s: rendimento tabela %.1f%% | parafuso mu=0,05 %.1f%% | mu implicito %.3f" % (k, eta, 100 * e5, mui))
    out.append("A tabela tem rendimentos diferentes com lead quase igual (3/60 e 6/60): o parafuso nao basta.")
    out.append("Temperaturas medidas (ate 6,4; 5,2; 2,4 C) sao pos-contato a 60 rpm, Pe << 10: nao comparaveis ao flash de Blok.")
    return out

def texto_validacao(L):
    nv = {c: sum(1 for x in L if x[5] == c) for c in "VRE"}
    s = ["VALIDACAO E VERIFICACAO",
         "V = verificacao matematica (identidade ou solucao fechada) | R = referencia publicada ou valor do MAIN | E = dado externo.",
         "MAIN 28.26 pede tolerancia antes do resultado: ate 21/10/2026 elas foram fixadas sem registro previo.",
         "hmin angular (Paouris) e CALIBRADO: variante escolhida nos mesmos 9 pontos (MAIN 28.17 nao satisfeito)."]
    for n, v, r, e, ok, c, tl in L:
        s.append("%s [%s] %-52s %10.4g vs %-10.4g (erro %.2g%%, tol %.2g%%)" % ("OK " if ok else "ERR", c, n, v, r, 100 * e, 100 * tl))
    n_ok = sum(1 for x in L if x[4])
    s.append("Resultado: %d/%d aprovados (V %d | R %d | E %d). So E e comparacao com dado externo." % (n_ok, len(L), nv["V"], nv["R"], nv["E"]))
    kbp, ch = calibracoes()
    s += ["", "CALIBRACOES REGISTRADAS (divisao calibra/valida, MAIN 27.20)",
          "K_BP = %.4f (ajustado nos 3 pares do Nakagawa; validado na geometria do Paouris); NAO aplicado ao calculo" % kbp,
          "C_h = %.3f (ajustado nos pontos impares do Paouris; validado nos pares); NAO aplicado ao calculo (|C_h-1| pequeno)" % ch,
          "Dominio do hmin calibrado: ue 1,5-8,7 m/s; F 230-10600 N; k 7-30; arrastamento angular.", ""] + comparar_nakagawa()
    return "\n".join(s)

# ---------------------------------------------------------------- RELATORIO
def relatorio(R, U=None, exp=False):
    g = lambda x: "%.4g" % x
    reg = "filme limite/rugoso (L<1)" if R["Lam"] < 1 else ("misto (1<=L<3)" if R["Lam"] < 3 else "filme pleno (L>=3)")
    nf_tag = "C" if SN[1][0] <= R["p0"] <= SN[0][0] else "N"; nh = nivel_h(R)
    hcs = ("%s um" % g(R["hc"] * 1e6)) if R["hc"] is not None else "n/v"
    s = [NIVEL_TXT, "Escopo: " + SOBRE_TXT, "=== CONTATO ===",
         "Rx=%s mm  Ry=%s mm  rho_red=%s mm" % (g(R["Rx"] * 1e3), g(R["Ry"] * 1e3), g(R["rred"] * 1e3)),
         "E*=%s GPa  E'(EHL)=%s GPa [V]" % (g(R["Es"] / 1e9), g(R["Ep"] / 1e9)),
         "a=%s mm (x)  b=%s mm (y)  k=b/a=%s [V]" % (g(R["a"] * 1e3), g(R["b"] * 1e3), g(R["k"])),
         "p0=pH=%s MPa  pm=%s MPa [V]" % (g(R["p0"] / 1e6), g(R["pm"] / 1e6)),
         "ue=%s m/s  v_desl=%s m/s" % (g(R["ue"]), g(R["vs"])),
         "=== FILME EHL ===",
         "nu(Tb)=%s cSt  eta0=%s Pa.s [V]" % (g(R["nuTb"]), g(R["eta0"])),
         "U=%s  G=%s  W=%s" % (g(R["U"]), g(R["G"]), g(R["W"])),
         "hmin=%s um [%s]  hc=%s  [%s]" % (g(R["hmin"] * 1e6), nh, hcs, R["mh"]),
         "Rq=%s um  Lambda=%s [%s]  -> %s" % (g(R["Rq"] * 1e6), g(R["Lam"]), nh, reg)]
    if exp:
        s += ["=== ATRITO E TEMPERATURA (mu e entrada) [N] ===",
              "q0=mu p0 v_desl=%s MW/m2" % g(R["q0"] / 1e6), "fracao do calor no corpo 1=%s" % g(R["f1"]),
              ("T_flash=n/a (Pe=%s < 10: Blok nao se aplica)" % g(R["Pe"])) if R["Tf"] is None else
              "T_flash=%s K  T_contato=%s C  Pe=%s" % (g(R["Tf"]), g(R["Tc"]), g(R["Pe"])),
              "=== DESGASTE E VIDA ===",
              "desgaste max./passagem no centro: corpo1=%s nm  corpo2=%s nm [N]" % (g(R["w1"] * 1e9), g(R["w2"] * 1e9)),
              "desgaste medio/passagem na faixa: corpo1=%s nm  corpo2=%s nm [N]" % (g(R["wmean1"] * 1e9), g(R["wmean2"] * 1e9)),
              "N=%s passagens do ponto material; max=%s um; medio=%s um [N]" % (g(R["N"]), g(R["wN"] * 1e6), g(R["wNmean"] * 1e6)),
              "N ate limite (conservador, linear, geometria fixa)=%s [N]" % g(R["Nw"]),
              "Nf(S-N MAIN)=%s ciclos; D=N/Nf=%s [%s]" % (g(R["Nf"]), g(R["D"]), nf_tag),
              "Sensibilidade S-N: +1%% em pH implica ~%.2f%% em Nf (lei de potencia); vida/N=%s" % (100 * R["elasticidade_Nf_pH"], g(R["vida_sobre_operacao"])),
              "Pressao admissivel para N=%s MPa (%s); margem pH=%s%%" % (g(R["pH_adm_N"] / 1e6) if R["pH_adm_N"] is not None else "n/a", R["faixa_N"], g(100 * R["margem_pH_N"]) if R["margem_pH_N"] is not None else "n/a"),
              "pH admissivel dentro de 700-800 MPa: %s" % ("sim" if R["pH_adm_dentro_SN"] else "nao; extrapolacao"),
              "MAIN: 5e7 ciclos e referencia de vida longa, nao terceiro ponto S-N.",
              "Curva S-N do MAIN em pH; com o historico (hist) ela e reescrita em Phi = amplitude de tau ortogonal (ver secao HISTORICO)."]
        s += ["=== v11: ATRITO MISTO (Eyring + asperezas) E COTA ESTACIONARIA [N] ===",
              "mu dinamico=%s = alfa_asp %s x mu_asp %s + (1-alfa) x mu_fluido %s; mu de entrada=%s; filme do fluido hf=%s um" % (g(R["mu_din"]), g(R["alfa_asp"]), g(R["mua"]), g(R["mu_fluido"]), g(R["mu"]), g(R["hf"] * 1e6)),
              ("T_flash (mu dinamico)=n/a (Pe=%s < 10)" % g(R["Pe"])) if R["Tf_din"] is None else "T_flash (mu dinamico)=%s K  T_contato=%s C" % (g(R["Tf_din"]), g(R["Tc_din"])),
              "T estacionaria (cota superior, conducao pura): mu entrada=%s K | mu dinamico=%s K. Nao e a temperatura de Blok; vale como teto quando Pe < 10." % (g(R["Tup"]), g(R["Tup_din"]))]
        if R.get("inc"):
            I = R["inc"]; h0, h1 = I["hist"][0], I["hist"][-1]
            s += ["=== v11: DESGASTE INCREMENTAL COM REALIMENTACAO (grade %dx%d por quadrante, Heun) [N] ===" % (I["n"], I["n"]),
                  "N alcancado=%s de %s%s; passos=%d (rejeitados por estabilidade: %d)" % (g(I["N_alcancado"]), g(R["N"]), "" if I["completo"] else " [INCOMPLETO: limite de passos]", I["passos"], I["rejeitados"]),
                  "w(0) incremental (soma dos 2 corpos)=%s um | linear estatico N(w1+w2)=%s um | razao=%s" % (g(I["w0"] * 1e6), g(I["w_lin_soma"] * 1e6), g(I["razao_w"])),
                  "desgaste maximo do corpo mais desgastado: incremental=%s um | linear estatico=%s um" % (g(I["wN_din"] * 1e6), g(R["wN"] * 1e6)),
                  "p_max: %s -> %s MPa | p0eq (4/3 p_ef): %s -> %s MPa | area de contato x%s | semi-largura em y: %s -> %s mm" % (g(h0["pmax"] / 1e6), g(h1["pmax"] / 1e6), g(h0["p0eq"] / 1e6), g(h1["p0eq"] / 1e6), g(h1["area"] / h0["area"]), g(h0["yb"] * 1e3), g(h1["yb"] * 1e3)),
                  "Miner com p_max(t) na S-N do MAIN: D=%s (estatico: %s) [N; p_max e sensivel a malha em bordas de sulco]" % (g(I["D_inc"]), g(R["D"])),
                  "Premissas estaticas: w depende so de y (th ignorado); N igual nos dois corpos; K, H constantes; semiespaco elastico (invalido se w >> a)."]
    else:
        s.append("Modulos sem calibracao (temperatura, atrito, desgaste, vida, rendimento) ocultos: ative em ENTRADA.")
    if "svm" in R:
        s += ["=== TENSOES (eixo z) [V] ===",
              "sigma_VM max=%s MPa em z=%s mm" % (g(R["svm"] / 1e6), g(R["zsv"] * 1e3)),
              "tau_max=%s MPa em z=%s mm" % (g(R["tmax"] / 1e6), g(R["ztm"] * 1e3))]
        s += ["=== VEREDITO VALIDADO (so tensao de von Mises) ===",
              "utilizacao %s -> %s" % (g(R["u_s"]), "INTEGRO" if R["u_s"] <= 1 else "FORA DO LIMITE")]
    if "hist" in R:
        Hh = R["hist"]; Vh = Hh["V"]; c_ = Vh[Hh["i_tau"]]; v_ = Vh[Hh["i_vm"]]
        s += ["=== HISTORICO DE TENSOES, ponto material fixo em y=0 (so pressao normal; mu nao entra) ===",
              "Ciclo H_c: janela de +-6 max(semi-corda, z) em torno do centro; compressao positiva; nucleo testado em --validar-hist",
              "tau ortogonal: amplitude max=%s MPa (%s p0) em z=%s mm" % (g(c_["tsz_amp"] / 1e6), g(c_["tsz_amp"] / R["p0"]), g(c_["z"] * 1e3)),
              "von Mises no ciclo: max=%s MPa em z=%s mm; faixa max-min=%s MPa" % (g(v_["svm_max"] / 1e6), g(v_["z"] * 1e3), g((v_["svm_max"] - v_["svm_min"]) / 1e6)),
              "Tresca max no ciclo=%s MPa (z=%s mm) | resolucao de z: %d profundidades (varredura discreta)" % (g(max(x["tresca_max"] for x in Vh) / 1e6), g(Vh[max(range(len(Vh)), key=lambda i: Vh[i]["tresca_max"])]["z"] * 1e3), len(Vh)),
              "S-N em Phi (MAIN 13.22): Phi=%s MPa; rho_ref=Phi/p0=%s%s; p_eq=%s MPa (%s); Nf=%s ciclos; D=N/Nf=%s [N]" % (g(Hh["Phi"] / 1e6), g(Hh["rho_ref"]), " (calibrado no contato nominal)" if Hh["rho_calibrado_no_nominal"] else "", g(Hh["p_eq"] / 1e6), Hh["faixa_phi"], g(Hh["Nf_phi"]), g(Hh["D_phi"])),
              "D em pH = %s (igual ao D em Phi no contato nominal, por construcao; muda quando o ciclo tem varios contatos/geometrias)" % g(Hh["D_main"]),
              "Premissa: " + Hh["motivo"],
              "Von Mises, Tresca: diagnostico [N]; nao entram em D nem no veredito.",
              "Phi(z)=tau_amp+lambda p_hid,max, lambda=%s [N]: Phi max=%s MPa em z=%s mm" % (g(Hh["lam"]), g(Hh["Phi"] / 1e6), g(Hh["zs"][Hh["i_tau"]] * 1e3)),
              "D(z): max=%s em z=%s mm; menor=%s em z=%s mm (%d profundidades). Sem H(z) cementado: Nf igual em todas as profundidades [N]" % (g(max(Hh["Dz"])), g(Hh["zs"][Hh["Dz"].index(max(Hh["Dz"]))] * 1e3), g(min(Hh["Dz"])), g(Hh["zs"][Hh["Dz"].index(min(Hh["Dz"]))] * 1e3), len(Hh["Dz"]))]
    ind = [("Filme minimo", R["u_h"], nh)]
    if exp or nf_tag == "C": ind.append(("Fadiga (Miner)", R["u_D"], nf_tag))
    if exp: ind += [("Temperatura", R["u_T"], "N"), ("Desgaste", R["u_w"], "N")]
    s.append("=== INDICATIVOS (nao entram no veredito) ===")
    s += ["%-16s %s [%s]" % (n, "n/a" if u is None else g(u), tg) for n, u, tg in ind]
    if U:
        s.append("=== FAIXAS DE INCERTEZA (MAIN s.29: nu40, nu100, alfa, Tb, Ra; n=%d; P5 / P50 / P95) ===" % U["n"])
        for nome, ch, un in (("hmin", "hmin", "um"), ("Lambda", "Lam", ""), ("T_contato", "Tc", "C")):
            if U[ch] and (exp or ch != "Tc"):
                lo, md, hi = faixa(U[ch]); s.append("%s: %s / %s / %s %s" % (nome, g(lo), g(md), g(hi), un))
        s.append("P(Lambda<1) = %.0f%%" % (100.0 * sum(1 for x in U["Lam"] if x < 1) / len(U["Lam"])))
        s.append("Nao inclui a incerteza do modelo de filme (ate 40%, Paouris 2015). Entradas amostradas de forma independente (MAIN diz que sao acopladas).")
    av = avisos_f(R, exp)
    if av:
        s.append("=== AVISOS ==="); s += av
    return "\n".join(s)

# ---------------------------------------------------------------- TEMA / GRAFICOS
BG, CARD, HDR, TXT, MUT = "#0b1220", "#16213a", "#0f1a30", "#e2e8f0", "#94a3b8"
ACC, OK, WARN, BAD, PURP, GRID = "#38bdf8", "#4ade80", "#fbbf24", "#f87171", "#c084fc", "#2b3a55"
ESC = 1.0                                   # escala de tela (px por dp), ajustada na interface

def fnt(n, b=False):
    return ("Helvetica", -int(n * ESC), "bold") if b else ("Helvetica", -int(n * ESC))

def plot(cv, titulo, xl, yl, series, hl=(), zero=True):
    cv.delete("all")
    W, H = int(cv["width"]), int(cv["height"])
    Lm, Rm, Tm, Bm = int(46 * ESC), int(12 * ESC), int(34 * ESC), int(44 * ESC)
    xs = [x for s in series for x in s[0]]
    ys = [y for s in series for y in s[1]] + [h[0] for h in hl]
    x0, x1 = min(xs), max(xs); y0, y1 = (min(0.0, min(ys)) if zero else min(ys) - 0.2), max(ys) * 1.08
    if x1 <= x0: x1 = x0 + 1
    if y1 <= y0: y1 = y0 + 1
    px = lambda x: Lm + (x - x0) / (x1 - x0) * (W - Lm - Rm)
    py = lambda y: H - Bm - (y - y0) / (y1 - y0) * (H - Tm - Bm)
    cv.create_text(W / 2, Tm / 2, text=titulo, fill=TXT, font=fnt(12, True))
    cv.create_text(W / 2, H - int(10 * ESC), text=xl, fill=MUT, font=fnt(10))
    cv.create_text(Lm, Tm - int(8 * ESC), text=yl, anchor="w", fill=MUT, font=fnt(10))
    for i in range(6):
        xv = x0 + (x1 - x0) * i / 5; yv = y0 + (y1 - y0) * i / 5
        cv.create_line(px(xv), Tm, px(xv), H - Bm, fill=GRID)
        cv.create_line(Lm, py(yv), W - Rm, py(yv), fill=GRID)
        cv.create_text(px(xv), H - Bm + int(12 * ESC), text="%.3g" % xv, fill=MUT, font=fnt(9))
        cv.create_text(Lm - int(4 * ESC), py(yv), text="%.3g" % yv, anchor="e", fill=MUT, font=fnt(9))
    cv.create_rectangle(Lm, Tm, W - Rm, H - Bm, outline="#475569")
    for v, cor, nome in hl:
        cv.create_line(Lm, py(v), W - Rm, py(v), fill=cor, dash=(5, 4), width=2)
        cv.create_text(W - Rm - int(4 * ESC), py(v) - int(9 * ESC), text=nome, anchor="e", fill=cor, font=fnt(10))
    for n, (xv, yv, cor, nome) in enumerate(series):
        cv.create_line(*[c for p in zip(map(px, xv), map(py, yv)) for c in p], fill=cor, width=max(2, int(1.5 * ESC)))
        y = Tm + int(12 * ESC) + int(17 * ESC) * n
        cv.create_line(Lm + int(8 * ESC), y, Lm + int(26 * ESC), y, fill=cor, width=max(2, int(1.5 * ESC)))
        cv.create_text(Lm + int(32 * ESC), y, text=nome, anchor="w", fill=cor, font=fnt(10))

def graficos(P, R, cvs, exp=False):
    n = 61
    xa = [-R["a"] + 2 * R["a"] * i / (n - 1) for i in range(n)]
    yb = [-R["b"] + 2 * R["b"] * i / (n - 1) for i in range(n)]
    px_ = [R["p0"] / 1e9 * math.sqrt(max(1 - (x / R["a"]) ** 2, 0)) for x in xa]
    py_ = [R["p0"] / 1e9 * math.sqrt(max(1 - (y / R["b"]) ** 2, 0)) for y in yb]
    plot(cvs[0], "Pressão de Hertz [V]", "posição (mm)", "p (GPa)",
         [([x * 1e3 for x in xa], px_, BAD, "p(x,0)"), ([y * 1e3 for y in yb], py_, ACC, "p(0,y)")])
    z = [q * 1e3 for q in R["zs"]]
    plot(cvs[1], "Tensões no eixo z [V]", "profundidade (mm)", "tensão (GPa)",
         [(z, [v / 1e9 for v in R["sv"]], BAD, "von Mises"), (z, [v / 1e9 for v in R["ta"]], OK, "τ máx"),
          (z, [v / 1e9 for v in R["sz"]], ACC, "σ z")], [(P["slim"] / 1e9, PURP, "limite VM")])
    xs, hs = [], []
    for i in range(15):
        s = 0.2 + 1.4 * i / 14
        r = calcular(dict(P, u1=P["u1"] * s, u2=P["u2"] * s), False)
        xs.append(r["ue"]); hs.append(r["hmin"] * 1e6)
    plot(cvs[2], "Filme mínimo × velocidade [indicativo]", "ue (m/s)  [escala u1,u2]", "hmin (µm)",
         [(xs, hs, ACC, "hmin")], [(R["Rq"] * 1e6, BAD, "Λ=1"), (P["hlim"] * 1e6, PURP, "hlim")])
    fs, tc = [], []
    for i in range(15):
        s = 0.3 + 1.4 * i / 14
        fs.append(P["F"] * s); tc.append(calcular(dict(P, F=P["F"] * s), False)["Tc"])
    if not exp or any(v is None for v in tc):
        cvs[3].delete("all")
        cvs[3].create_text(int(cvs[3]["width"]) / 2, int(cvs[3]["height"]) / 2, width=int(cvs[3]["width"]) - 40,
                           text=("Módulo de temperatura não calibrado: oculto (ative em ENTRADA)" if not exp else "Temperatura n/a: Pe < 10, Blok não se aplica [N]"),
                           fill=WARN, font=fnt(12))
    else:
        plot(cvs[3], "Temperatura de contato × carga [N]", "F (N)", "Tc (°C)",
             [(fs, tc, BAD, "Tc")], [(P["Tlim"], PURP, "Tlim")])
    fs, uT, uh, uw, uv, uD = [], [], [], [], [], []
    for i in range(15):
        s = 0.3 + 1.4 * i / 14
        r = calcular(dict(P, F=P["F"] * s), False)
        fs.append(P["F"] * s); uT.append(r["u_T"]); uh.append(r["u_h"]); uw.append(r["u_w"])
        uv.append(R["svm"] / P["slim"] * s ** (1 / 3.0)); uD.append(r["u_D"])
    ser = [(fs, uv, WARN, "von Mises [V]"), (fs, uh, ACC, "filme [%s]" % nivel_h(R))]
    if exp:
        ser += [(fs, uw, OK, "desgaste [N]"), (fs, uD, PURP, "fadiga")]
        if None not in uT: ser.append((fs, uT, BAD, "temperatura [N]"))
    plot(cvs[4], "Utilizações × carga (só von Mises é validado)", "F (N)", "utilização", ser, [(1.0, TXT, "limite")])
    ps = [0.4e9 + 2.1e9 * i / 40 for i in range(41)]
    if exp or SN[1][0] <= R["p0"] <= SN[0][0]:
        plot(cvs[5], "Curva S-N do MAIN (calibrada em 700-800 MPa) [C]", "pH (GPa)", "log10 N",
             [([p / 1e9 for p in ps], [math.log10(SN[0][1] * (SN[0][0] / p) ** SN_M) for p in ps], BAD, "S-N (700-800 MPa)")],
             [(math.log10(P["N"]), WARN, "N operação"), (math.log10(R["Nf"]), ACC, "Nf")], zero=False)
    else:
        cvs[5].delete("all")
        cvs[5].create_text(int(cvs[5]["width"]) / 2, int(cvs[5]["height"]) / 2, width=int(cvs[5]["width"]) - 40,
                           text="Vida: pH fora de 700-800 MPa (faixa calibrada no MAIN). Oculta (ative em ENTRADA).", fill=WARN, font=fnt(12))

def graficos_hrh(H, cvs, exp=False):
    amax = H["Am"] * math.tan(max(0.01, 85 * DEG - H["bg"]))
    ao = [amax * i / 24 for i in range(25)]
    L = [hrh(H["Am"], a, H["bg"], H["al"], H["gm"], H["T"], H["n1"], H["nc"], H["mu"]) for a in ao]
    x = [a * 1e3 for a in ao]
    plot(cvs[0], "Ângulos × offset [C]", "offset a (mm)", "graus",
         [(x, [r["bp"] / DEG for r in L], BAD, "β espiral pinhão"), (x, [r["lam"] / DEG for r in L], ACC, "avanço (lead)"),
          (x, [r["psi"] / DEG for r in L], OK, "ângulo de arrastamento 2·lead")])
    plot(cvs[1], "Velocidades × offset [C]", "offset a (mm)", "m/s",
         [(x, [r["vs"] for r in L], BAD, "v deslizamento"), (x, [r["vg"] for r in L], ACC, "v coroa"),
          (x, [r["ue"] for r in L], OK, "u arrastamento")])
    plot(cvs[2], "Forças × offset [V]", "offset a (mm)", "N",
         [(x, [r["Fn"] for r in L], BAD, "Fn normal"), (x, [r["Ft"] for r in L], ACC, "Ft tangencial"),
          (x, [abs(r["Fy"]) for r in L], OK, "|Fy|"), (x, [abs(r["Fz"]) for r in L], PURP, "|Fz|")])
    ser = []
    for f, cor in ((0.5, OK), (1.0, ACC), (1.5, BAD)):
        rr = [hrh(H["Am"], a, H["bg"], H["al"], H["gm"], H["T"], H["n1"], H["nc"], H["mu"] * f) for a in ao]
        ser.append((x, [100 * r["eta"] for r in rr], cor, "μ × %.1f" % f))
    if exp:
        plot(cvs[3], "Rendimento do parafuso × offset [N]", "offset a (mm)", "η (%)", ser)
    else:
        cvs[3].delete("all")
        cvs[3].create_text(int(cvs[3]["width"]) / 2, int(cvs[3]["height"]) / 2, width=int(cvs[3]["width"]) - 40,
                           text="Rendimento não calibrado: oculto (ative em ENTRADA)", fill=WARN, font=fnt(12))

def graficos_v11(R, cvs):
    """Graficos v11 (canvases 6 a 8): desgaste incremental x linear, pressao x N, D(z). Mensagem quando o modulo opt-in nao foi executado."""
    def msg(cv, t):
        cv.delete("all"); cv.create_text(int(cv["width"]) / 2, int(cv["height"]) / 2, width=int(cv["width"]) - 40, text=t, fill=WARN, font=fnt(12))
    I = R.get("inc")
    if I:
        h = I["hist"]; xs = [x["N"] for x in h]
        plot(cvs[6], "Desgaste no centro × N [N]", "N (passagens)", "w (µm)",
             [(xs, [x["w0"] * 1e6 for x in h], BAD, "incremental (realimentado)"), (xs, [x["N"] * (R["w1"] + R["w2"]) * 1e6 for x in h], ACC, "linear estático")])
        plot(cvs[7], "Pressão × N [N]", "N (passagens)", "p (GPa)",
             [(xs, [x["pmax"] / 1e9 for x in h], WARN, "p_max (sensível a malha)"), (xs, [x["p0eq"] / 1e9 for x in h], OK, "p0eq = 4/3 p_ef")], [(R["p0"] / 1e9, PURP, "Hertz inicial")])
    else:
        msg(cvs[6], "Desgaste incremental não calculado: use DESGASTE INCREMENTAL em ENTRADA"); msg(cvs[7], "Pressão × N: idem")
    if R.get("hist"):
        Hh = R["hist"]; plot(cvs[8], "Dano D(z) por profundidade [N]", "profundidade (mm)", "D", [([z * 1e3 for z in Hh["zs"]], Hh["Dz"], BAD, "D(z)")], [(1.0, TXT, "D = 1")])
    else:
        msg(cvs[8], "D(z) não calculado: use HISTÓRICO DE TENSÕES em ENTRADA")

# ---------------------------------------------------------------- CICLOS POR DENTE
def ciclos_por_dente(z1, z2, rpm_pinhao, horas, corpo="pinhao"):
    """Ciclos de engrenamento por dente: um evento por revolucao daquele elemento."""
    vals = (z1, z2, rpm_pinhao, horas)
    if not all(math.isfinite(float(v)) for v in vals):
        raise ValueError("Entradas de ciclos devem ser numeros finitos.")
    if z1 <= 0 or z2 <= 0 or rpm_pinhao < 0 or horas < 0:
        raise ValueError("Z1/Z2 devem ser > 0; rpm e horas devem ser >= 0.")
    rpm = rpm_pinhao if corpo == "pinhao" else rpm_pinhao * z1 / z2
    return rpm * 60.0 * horas

# ---------------------------------------------------------------- INTERFACE
def interface():
    import tkinter as tk
    global ESC
    root = tk.Tk(); root.title("TRIBOCONTATOGRACIA11"); root.configure(bg=BG)
    sw = root.winfo_screenwidth()
    ESC = max(1.0, sw / 360.0)
    pd = int(8 * ESC); wc = sw - 4 * pd
    ent, est, pag, abas = {}, {}, {}, {}
    expv = tk.BooleanVar(value=False)

    # ---- cabecalho e abas
    hd = tk.Frame(root, bg=HDR); hd.pack(fill="x")
    tk.Label(hd, text="TRIBOCONTATO", bg=HDR, fg=ACC, font=fnt(20, True)).pack(pady=(pd, 0))
    tk.Label(hd, text="contato aço–aço hipoide · alta redução", bg=HDR, fg=MUT, font=fnt(11)).pack(pady=(0, pd))
    tb = tk.Frame(root, bg=BG); tb.pack(fill="x", pady=(pd // 2, 0))
    body = tk.Frame(root, bg=BG); body.pack(fill="both", expand=True)

    def nova(nome, rot):
        host = tk.Frame(body, bg=BG)
        cv = tk.Canvas(host, bg=BG, highlightthickness=0)
        sb = tk.Scrollbar(host, orient="vertical", command=cv.yview, width=int(16 * ESC))
        cv.configure(yscrollcommand=sb.set); sb.pack(side="right", fill="y"); cv.pack(side="left", fill="both", expand=True)
        inn = tk.Frame(cv, bg=BG); jn = cv.create_window((0, 0), window=inn, anchor="nw")
        inn.bind("<Configure>", lambda e: cv.configure(scrollregion=cv.bbox("all")))
        cv.bind("<Configure>", lambda e: cv.itemconfigure(jn, width=e.width))
        pag[nome] = (host, cv, inn)
        b = tk.Button(tb, text=rot, command=lambda: ir(nome), font=fnt(10, True), relief="flat", bd=0,
                      bg=CARD, fg=TXT, activebackground=ACC, pady=int(7 * ESC))
        b.pack(side="left", fill="x", expand=True, padx=1); abas[nome] = b
        return inn

    def ir(nome):
        for h, c, i in pag.values(): h.pack_forget()
        pag[nome][0].pack(fill="both", expand=True); est["cv"] = pag[nome][1]
        for n, b in abas.items():
            b.configure(bg=ACC if n == nome else CARD, fg=BG if n == nome else TXT)
        if nome == "VAL" and "val" not in est: validar_ui()

    # arrastar com o dedo (nao em campos/botoes)
    est["y"] = 0
    def _p(e): est["y"] = e.y_root
    def _m(e):
        if e.widget.winfo_class() in ("Entry", "Button", "Scrollbar", "Text"): return
        cv = est.get("cv")
        if cv is None: return
        dy = est["y"] - e.y_root; est["y"] = e.y_root
        cv.yview_scroll(int(dy / max(1, 6 * ESC)), "units")
    root.bind_all("<ButtonPress-1>", _p, add="+"); root.bind_all("<B1-Motion>", _m, add="+")
    root.bind_all("<MouseWheel>", lambda e: est["cv"].yview_scroll(-1 if e.delta > 0 else 1, "units"), add="+")

    def cartao(pai, titulo=None, sub=None):
        f = tk.Frame(pai, bg=CARD, padx=int(10 * ESC), pady=int(8 * ESC)); f.pack(fill="x", padx=pd, pady=pd // 2)
        if titulo: tk.Label(f, text=titulo, bg=CARD, fg=ACC, font=fnt(13, True), anchor="w").pack(fill="x")
        if sub: tk.Label(f, text=sub, bg=CARD, fg=MUT, font=fnt(10), anchor="w", justify="left", wraplength=wc - 6 * pd).pack(fill="x", pady=(0, 4))
        return f

    def botao(pai, txt, cmd, cor=ACC, fg=BG):
        b = tk.Button(pai, text=txt, command=cmd, bg=cor, fg=fg, activebackground=cor, font=fnt(15, True),
                      relief="flat", bd=0, pady=int(9 * ESC)); b.pack(fill="x", padx=pd, pady=pd // 2 + 2); return b

    # ---- paginas
    pE = nova("ENT", "ENTR."); pR = nova("RES", "RESULT."); pG = nova("GRA", "GRÁF.")
    pH = nova("HRH", "HRH"); pV = nova("VAL", "VALID."); pT = nova("TXT", "TXT")

    # ENTRADA
    cs = cartao(pE, "O que é o TRIBOCONTATO", SOBRE_TXT.replace("-", "–") if False else SOBRE_TXT)
    tk.Checkbutton(cs, text="Mostrar módulos não calibrados (temperatura, atrito, desgaste, vida, rendimento)", variable=expv,
                   command=lambda: (calc(False), calc_h()), bg=CARD, fg=TXT, selectcolor=BG, activebackground=CARD,
                   activeforeground=TXT, font=fnt(11), anchor="w", justify="left", wraplength=wc - 8 * pd).pack(fill="x", pady=(4, 0))
    card, g = None, 0
    for c in CAMPOS:
        if len(c) == 1:
            raw = c[0][2:]; tit = raw.split(" (")[0]; sub = raw[len(tit):].strip(" ()")
            card = cartao(pE, tit, sub or None); card.columnconfigure(0, weight=1); g = 0
            continue
        cf = tk.Frame(card, bg=CARD); cf.pack(fill="x", pady=2); cf.columnconfigure(0, weight=1)
        tk.Label(cf, text=c[1], bg=CARD, fg=TXT, font=fnt(12), anchor="w").grid(row=0, column=0, sticky="w")
        e = tk.Entry(cf, width=9, justify="right", font=fnt(14), bg=BG, fg=TXT, insertbackground=TXT, relief="flat",
                     highlightthickness=1, highlightbackground="#334155", highlightcolor=ACC)
        e.insert(0, "%g" % c[3]); e.grid(row=0, column=1, padx=int(6 * ESC), ipady=int(3 * ESC))
        tk.Label(cf, text=c[2], bg=CARD, fg=MUT, font=fnt(11), width=7, anchor="w").grid(row=0, column=2)
        ent[c[0]] = e

    # Conversor explicito de rpm/tempo para ciclos por dente, sem confundir ciclos do pinhao e da coroa.
    cc = cartao(pE, "CONVERSOR DE CICLOS POR DENTE", "Um ciclo por dente a cada revolucao do respectivo elemento. O resultado pode preencher N automaticamente.")
    cvars = {}
    for key, label, default in (("z1", "Dentes do pinhao Z1", "2"), ("z2", "Dentes da coroa Z2", "60"),
                                ("rpm", "Rotacao do pinhao (rpm)", "3000"), ("horas", "Tempo de operacao (h)", "1")):
        rr = tk.Frame(cc, bg=CARD); rr.pack(fill="x", pady=2)
        tk.Label(rr, text=label, bg=CARD, fg=TXT, font=fnt(11), anchor="w").pack(side="left", fill="x", expand=True)
        ee = tk.Entry(rr, width=10, justify="right", font=fnt(13), bg=BG, fg=TXT, insertbackground=TXT, relief="flat")
        ee.insert(0, default); ee.pack(side="right"); cvars[key] = ee
    corpo_ciclo = tk.StringVar(value="pinhao")
    rb = tk.Frame(cc, bg=CARD); rb.pack(fill="x", pady=3)
    tk.Radiobutton(rb, text="Ciclo por dente do pinhao", variable=corpo_ciclo, value="pinhao", bg=CARD, fg=TXT, selectcolor=BG,
                   activebackground=CARD, font=fnt(10)).pack(anchor="w")
    tk.Radiobutton(rb, text="Ciclo por dente da coroa", variable=corpo_ciclo, value="coroa", bg=CARD, fg=TXT, selectcolor=BG,
                   activebackground=CARD, font=fnt(10)).pack(anchor="w")
    ciclo_msg = tk.Label(cc, text="2/60 a 3000 rpm: pinhao 180.000 ciclos/h; coroa 6.000 ciclos/h.", bg=CARD, fg=MUT,
                         font=fnt(10), anchor="w", justify="left", wraplength=wc - 8 * pd)
    ciclo_msg.pack(fill="x", pady=3)
    def aplicar_ciclos():
        try:
            z1 = float(cvars["z1"].get().replace(",", ".")); z2 = float(cvars["z2"].get().replace(",", "."))
            rpm = float(cvars["rpm"].get().replace(",", ".")); horas = float(cvars["horas"].get().replace(",", "."))
            ncy = ciclos_por_dente(z1, z2, rpm, horas, corpo_ciclo.get())
            ent["N"].delete(0, "end"); ent["N"].insert(0, "%.8g" % ncy)
            ciclo_msg.configure(text="N preenchido: %.8g ciclos por dente do %s. Confira se esse e o corpo cuja vida deseja avaliar." %
                                (ncy, "pinhao" if corpo_ciclo.get() == "pinhao" else "coroa"), fg=OK)
        except Exception as ex:
            ciclo_msg.configure(text="Erro no conversor: " + str(ex), fg=BAD)
    tk.Button(cc, text="CALCULAR CICLOS E PREENCHER N", command=aplicar_ciclos, bg=ACC, fg=BG,
              font=fnt(11, True), relief="flat", bd=0, pady=int(7 * ESC)).pack(fill="x", pady=4)

    def padroes():
        for c in CAMPOS:
            if len(c) > 1:
                ent[c[0]].delete(0, "end"); ent[c[0]].insert(0, "%g" % c[3])
    botao(pE, "CALCULAR", lambda: calc(True))
    botao(pE, "HISTÓRICO DE TENSÕES + S-N em Φ (mais lento)", lambda: calc_hist(), CARD, ACC)
    botao(pE, "DESGASTE INCREMENTAL (realimentação, lento)", lambda: calc_inc(), CARD, ACC)
    lb_h = tk.Label(pE, text="", bg=BG, fg=OK, font=fnt(10), wraplength=wc, justify="left"); lb_h.pack()
    botao(pE, "Restaurar padrões", padroes, CARD, TXT)

    def ler():
        P = {}
        for c in CAMPOS:
            if len(c) > 1:
                P[c[0]] = float(ent[c[0]].get().replace(",", ".")) * c[4]
        for ch in ("u1", "u2", "F", "E", "H", "wlim", "hlim", "Tlim"):
            if P[ch] <= 0: raise ValueError(ch + " deve ser > 0")
        for ch in ("R1x", "R1y", "R2x", "R2y"):
            if P[ch] == 0: raise ValueError(ch + " nao pode ser zero")
        for ch, val in P.items():
            if not math.isfinite(val): raise ValueError(ch + " deve ser finito")
        if not (-1.0 < P["nu"] < 0.5): raise ValueError("Poisson nu deve estar entre -1 e 0,5")
        if P["N"] < 0: raise ValueError("N deve ser >= 0")
        if P["K"] < 0: raise ValueError("K de Archard deve ser >= 0")
        if P["mu"] < 0: raise ValueError("mu deve ser >= 0")
        if P["mua"] < 0 or P["lam"] < 0: raise ValueError("mu_asp e lambda devem ser >= 0")
        if P["tau0"] <= 0: raise ValueError("tau0 deve ser > 0")
        if P["nu40"] <= 0 or P["nu100"] <= 0 or P["rho"] <= 0: raise ValueError("Viscosidades e densidade devem ser > 0")
        Rx = 1.0 / (1.0 / P["R1x"] + 1.0 / P["R2x"])
        Ry = 1.0 / (1.0 / P["R1y"] + 1.0 / P["R2y"])
        if Rx <= 0 or Ry <= 0: raise ValueError("Raios equivalentes Rx e Ry devem ser positivos para este modelo de contato")
        return P

    # GRAFICOS
    wg = wc - 2 * int(10 * ESC) - 2 * pd
    cvs = []
    for _ in range(9):
        cd = cartao(pG)
        cvs.append(tk.Canvas(cd, width=wg, height=int(wg * 0.8), bg=BG, highlightthickness=0)); cvs[-1].pack()

    # RELATORIO
    cdT = cartao(pT)
    txt = tk.Text(cdT, width=max(20, int(wg / (6 * ESC))), height=48, font=("Courier", -int(10 * ESC)), wrap="word",
                  bg=BG, fg=TXT, relief="flat", insertbackground=TXT)
    txt.pack(fill="both")
    def atualizar_txt():
        txt.delete("1.0", "end"); txt.insert("end", est.get("s", "") + "\n\n" + est.get("h", "") + "\n\n" + est.get("v", ""))
    def salvar():
        pasta = "/storage/emulated/0/Download"
        if not os.path.isdir(pasta): pasta = os.getcwd()
        cam = os.path.join(pasta, "TRIBOCONTATO_resultado.txt")
        with open(cam, "w", encoding="utf-8") as f:
            f.write(est.get("s", "") + "\n\n" + est.get("h", "") + "\n\n" + est.get("v", ""))
        lb.configure(text="Salvo em " + cam)
    botao(pT, "SALVAR TXT", salvar)
    lb = tk.Label(pT, text="", bg=BG, fg=OK, font=fnt(10), wraplength=wc); lb.pack()

    def limpar(f):
        for w in f.winfo_children(): w.destroy()
    # RESULTADO
    def limpar_(f):
        for w in f.winfo_children(): w.destroy()
    def cor_u(u): return OK if u < 0.7 else (WARN if u <= 1 else BAD)

    def mostrar(R, U=None, exp=False):
        limpar(pR)
        nf_tag = "C" if SN[1][0] <= R["p0"] <= SN[0][0] else "N"; nh = nivel_h(R)
        okk = R["u_s"] <= 1
        bn = tk.Frame(pR, bg=OK if okk else BAD, padx=pd, pady=pd); bn.pack(fill="x", padx=pd, pady=pd // 2)
        tk.Label(bn, text="ÍNTEGRO (tensão)" if okk else "FORA DO LIMITE (tensão)", bg=bn["bg"], fg=BG, font=fnt(20, True)).pack()
        tk.Label(bn, text="veredito só com o validado: von Mises, utilização %.2f" % R["u_s"], bg=bn["bg"], fg=BG, font=fnt(12)).pack()
        lam = R["Lam"]
        reg = "limite" if lam < 1 else ("misto" if lam < 3 else "filme pleno")
        flash_ok = R["Tf"] is not None
        cf_ = TXT if nh in ("C", "R") else MUT
        tiles = [("Pressão pH", "%.0f" % (R["p0"] / 1e6), "MPa · V", TXT), ("Elipse k = b/a", "%.2f" % R["k"], "V", TXT),
                 ("Semi-eixos a × b", "%.3f × %.3f" % (R["a"] * 1e3, R["b"] * 1e3), "mm · V", TXT),
                 ("von Mises máx", "%.0f" % (R["svm"] / 1e6), "MPa · V", cor_u(R["u_s"])),
                 ("τ máx", "%.0f" % (R["tmax"] / 1e6), "MPa · V", TXT),
                 ("Viscosidade η0", "%.4f" % R["eta0"], "Pa·s · V", TXT),
                 ("Filme mínimo", "%.3f" % (R["hmin"] * 1e6), "µm · " + nh, cf_),
                 ("Λ = hmin/Rq", "%.2f" % lam, reg + " · " + nh, cf_)]
        if exp:
            tiles += [("Temp. contato", "%.1f" % R["Tc"] if flash_ok else "n/a", "°C · N" if flash_ok else "Pe<10 · N", MUT),
                      ("Flash de Blok", "%.1f" % R["Tf"] if flash_ok else "n/a", "K · N" if flash_ok else "Pe<10 · N", MUT),
                      ("Vida Nf", "%.2e" % R["Nf"], "ciclos · " + nf_tag, MUT),
                      ("pH admissível para N", "%.0f" % (R["pH_adm_N"] / 1e6) if R["pH_adm_N"] is not None else "n/a", "MPa · " + R["faixa_N"], MUT),
                      ("Desgaste em N", "%.3f" % (R["wN"] * 1e6), "µm · N", MUT),
                      ("N até limite desgaste", "%.2e" % R["Nw"], "passagens · N", MUT),
                      ("μ dinâmico (misto)", "%.3f" % R["mu_din"], "α asp. %.2f · N" % R["alfa_asp"], MUT),
                      ("T estac. (teto)", "%.0f" % R["Tup_din"], "K · N", MUT)]
        gr = tk.Frame(pR, bg=BG); gr.pack(fill="x", padx=pd // 2)
        gr.columnconfigure(0, weight=1, uniform="t"); gr.columnconfigure(1, weight=1, uniform="t")
        for i, (t, v, u, cr) in enumerate(tiles):
            tl = tk.Frame(gr, bg=CARD, padx=int(8 * ESC), pady=int(7 * ESC))
            tl.grid(row=i // 2, column=i % 2, sticky="nsew", padx=pd // 2, pady=pd // 2)
            tk.Label(tl, text=t, bg=CARD, fg=MUT, font=fnt(10), anchor="w").pack(fill="x")
            tk.Label(tl, text=v, bg=CARD, fg=cr, font=fnt(18, True), anchor="w").pack(fill="x")
            tk.Label(tl, text=u or " ", bg=CARD, fg=MUT, font=fnt(10), anchor="w").pack(fill="x")
        uts = [("von Mises [V]", R["u_s"], cor_u(R["u_s"])), ("Filme [%s]" % nh, R["u_h"], "#64748b")]
        if exp or nf_tag == "C": uts.append(("Fadiga [%s]" % nf_tag, R["u_D"], "#64748b"))
        if exp: uts += [("Temperatura [N]", R["u_T"], "#64748b"), ("Desgaste [N]", R["u_w"], "#64748b")]
        cd = cartao(pR, "Utilizações", "só a barra colorida (von Mises) decide · cinzas são indicativas · linha branca = limite")
        wb = wg; hb = int(30 * ESC) * len(uts) + int(8 * ESC)
        cb = tk.Canvas(cd, width=wb, height=hb, bg=CARD, highlightthickness=0); cb.pack()
        x0, x1 = int(wb * 0.34), wb - int(8 * ESC); sc = (x1 - x0) / 1.5
        for i, (n, u, cor) in enumerate(uts):
            y = int(30 * ESC) * i + int(15 * ESC)
            cb.create_text(0, y, text=n, anchor="w", fill=TXT, font=fnt(11))
            cb.create_rectangle(x0, y - int(8 * ESC), x1, y + int(8 * ESC), fill=BG, outline="")
            if u is not None:
                cb.create_rectangle(x0, y - int(8 * ESC), x0 + min(u, 1.5) * sc, y + int(8 * ESC), fill=cor, outline="")
            cb.create_text(x1, y - int(14 * ESC), text="n/a" if u is None else "%.2f" % u, anchor="e", fill=MUT, font=fnt(9))
        cb.create_line(x0 + sc, 0, x0 + sc, hb, fill="#ffffff", width=2)
        if U:
            ci = cartao(pR, "Faixas de incerteza", "MAIN §29 · ν40, ν100, αp, Tb, Ra amostrados nas faixas da tabela-mestra · P5 / P50 / P95")
            linhas = []
            for nome, ch, un, fm in (("hmin", "hmin", "µm", "%.3f"), ("Λ", "Lam", "", "%.2f"), ("Temp. contato", "Tc", "°C", "%.1f")):
                if U[ch] and (exp or ch != "Tc"):
                    lo, md, hi = faixa(U[ch]); linhas.append("%s: %s / %s / %s %s" % (nome, fm % lo, fm % md, fm % hi, un))
            linhas.append("P(Λ < 1) = %.0f%%" % (100.0 * sum(1 for x in U["Lam"] if x < 1) / len(U["Lam"])))
            linhas.append("Não inclui a incerteza do modelo de filme (até 40%, Paouris 2015). O MAIN diz que as variáveis são acopladas; aqui são amostradas de forma independente.")
            for ln in linhas:
                tk.Label(ci, text=ln, bg=CARD, fg=TXT, font=fnt(11), anchor="w", justify="left", wraplength=wg).pack(fill="x", pady=1)
        av = avisos_f(R, exp)
        if av:
            ca = cartao(pR, "Avisos")
            for a in av:
                tk.Label(ca, text="• " + a, bg=CARD, fg=WARN, font=fnt(11), anchor="w", justify="left", wraplength=wg).pack(fill="x")

    def calc(vai):
        try:
            P = ler(); R = calcular(P); U = incerteza(P)
            est["s"] = relatorio(R, U, expv.get()); mostrar(R, U, expv.get()); graficos(P, R, cvs, expv.get()); graficos_v11(R, cvs); atualizar_txt()
            if vai: ir("RES")
        except Exception as ex:
            limpar(pR)
            tk.Label(pR, text="ERRO: %s" % ex, bg=BG, fg=BAD, font=fnt(13), wraplength=wc, justify="left").pack(padx=pd, pady=pd)
            if vai: ir("RES")


    def calc_hist():
        try:
            P = ler(); lb_h.configure(text="Calculando histórico...", fg=OK); root.update_idletasks()
            R = calcular(P, hist=True); est["s"] = relatorio(R, None, expv.get()); graficos(P, R, cvs, expv.get()); graficos_v11(R, cvs); atualizar_txt()
            lb_h.configure(text="Histórico calculado: ver aba TXT.", fg=OK); ir("TXT")
        except Exception as ex:
            lb_h.configure(text="ERRO: %s" % ex, fg=BAD)

    def calc_inc():
        try:
            P = ler(); lb_h.configure(text="Calculando desgaste incremental...", fg=OK); root.update_idletasks()
            R = calcular(P, inc=True); est["s"] = relatorio(R, None, expv.get()); graficos(P, R, cvs, expv.get()); graficos_v11(R, cvs); atualizar_txt()
            I = R["inc"]
            lb_h.configure(text="Desgaste incremental: w(0)=%.3g um (linear estático %.3g um), %d passos%s. Ver GRÁF. e TXT." % (I["w0"] * 1e6, I["w_lin_soma"] * 1e6, I["passos"], "" if I["completo"] else " [INCOMPLETO: limite de passos]"),
                           fg=OK if I["completo"] else WARN); ir("GRA")
        except Exception as ex:
            lb_h.configure(text="ERRO: %s" % ex, fg=BAD)

    # HRH
    ceh = {}
    ch = cartao(pH, "HRH · offset alto", "Relações de Stadtfeld (Gear Technology, ago/2011) e cinemática de linha de contato. Valores padrão são ilustrativos, não do artigo.")
    for c in CAMPOS_HRH:
        cf = tk.Frame(ch, bg=CARD); cf.pack(fill="x", pady=2); cf.columnconfigure(0, weight=1)
        tk.Label(cf, text=c[1], bg=CARD, fg=TXT, font=fnt(12), anchor="w").grid(row=0, column=0, sticky="w")
        e = tk.Entry(cf, width=9, justify="right", font=fnt(14), bg=BG, fg=TXT, insertbackground=TXT, relief="flat",
                     highlightthickness=1, highlightbackground="#334155", highlightcolor=ACC)
        e.insert(0, "%g" % c[3]); e.grid(row=0, column=1, padx=int(6 * ESC), ipady=int(3 * ESC))
        tk.Label(cf, text=c[2], bg=CARD, fg=MUT, font=fnt(11), width=7, anchor="w").grid(row=0, column=2)
        ceh[c[0]] = e
    def ler_h():
        H = {c[0]: float(ceh[c[0]].get().replace(",", ".")) * c[4] for c in CAMPOS_HRH}
        H["mu"] = float(ent["mu"].get().replace(",", ".")) * 1.0
        if H["Am"] <= 0 or H["nc"] <= 0 or (H["bp0"] if H["bp0"] > 0 else H["bg"] + math.atan(H["a"] / H["Am"])) >= PI / 2:
            raise ValueError("beta_p deve ser < 90 graus")
        return H
    hres = tk.Frame(pH, bg=BG); hres_ref = {}
    wgh = wg
    hcvs = []
    def calc_h():
        try:
            H = ler_h(); r = hrh(H["Am"], H["a"], H["bg"], H["al"], H["gm"], H["T"], H["n1"], H["nc"], H["mu"], H["bp0"])
            hres_ref["r"] = r
            limpar(hres)
            it = [("β espiral pinhão", "%.2f" % (r["bp"] / DEG), "° · C"), ("Avanço (lead)", "%.2f" % (r["lam"] / DEG), "° · C"),
                  ("Arrastamento 2·lead", "%.2f" % (r["psi"] / DEG), "° · C"), ("v pinhão", "%.2f" % r["vp"], "m/s · V"),
                  ("v coroa", "%.3f" % r["vg"], "m/s · C"), ("v deslizamento", "%.2f" % r["vs"], "m/s · C"),
                  ("u arrastamento", "%.2f" % r["ue"], "m/s · C"), ("Ft", "%.0f" % r["Ft"], "N · V"),
                  ("Fn normal", "%.0f" % r["Fn"], "N · V"), ("Fn por par", "%.0f" % r["Fnc"], "N · C"),
                  ("u1 / u2 no arrastamento", "%.2f / %.3f" % (r["u1"], r["u2"]), "m/s · C")]
            if expv.get(): it.insert(10, ("Rendimento", "%.1f" % (100 * r["eta"]), "% · N"))
            gr = tk.Frame(hres, bg=BG); gr.pack(fill="x", padx=pd // 2)
            gr.columnconfigure(0, weight=1, uniform="h"); gr.columnconfigure(1, weight=1, uniform="h")
            for i, (t_, v_, u_) in enumerate(it):
                tl = tk.Frame(gr, bg=CARD, padx=int(8 * ESC), pady=int(7 * ESC))
                tl.grid(row=i // 2, column=i % 2, sticky="nsew", padx=pd // 2, pady=pd // 2)
                tk.Label(tl, text=t_, bg=CARD, fg=MUT, font=fnt(10), anchor="w").pack(fill="x")
                tk.Label(tl, text=v_, bg=CARD, fg=TXT, font=fnt(18, True), anchor="w").pack(fill="x")
                tk.Label(tl, text=u_, bg=CARD, fg=MUT, font=fnt(10), anchor="w").pack(fill="x")
            graficos_hrh(H, hcvs, expv.get())
            est["h"] = "HRH: " + ", ".join("%s=%s%s" % (a, b, c) for a, b, c in it)
            atualizar_txt()
        except Exception as ex:
            limpar(hres); tk.Label(hres, text="ERRO: %s" % ex, bg=BG, fg=BAD, font=fnt(13), wraplength=wc).pack()
    def enviar():
        r = hres_ref.get("r")
        if r is None: calc_h(); r = hres_ref.get("r")
        if r is None: return
        for ch_, v_ in (("F", r["Fnc"]), ("u1", r["u1"]), ("u2", r["u2"]), ("th", (PI / 2 - r["psi"]) / DEG)):
            ent[ch_].delete(0, "end"); ent[ch_].insert(0, "%.4g" % v_)
        lbh.configure(text="F, u1, u2 e ângulo de arrastamento enviados à ENTRADA. Raios de curvatura não vêm do HRH: ajuste-os.")
    def preset(k):
        z1, ii, Amp, Amg, gp, dg, bpt, bgt, eps, eta, ar = NAK[k]
        for ch_, v_ in (("Am", Amp), ("a", NAK_A), ("bg", bgt), ("al", NAK_AL), ("gm", gp), ("T", NAK_T / ii),
                        ("n1", 60.0), ("nc", ar), ("bp0", bpt)):
            ceh[ch_].delete(0, "end"); ceh[ch_].insert(0, "%g" % v_)
        calc_h()
        lbh.configure(text="Caso Nakagawa %s: Tabela 1, 100 N.m na coroa, pinhão a 60 rpm; pares = razão de contato real lida da Fig. 22." % k)
    def preset_ag():
        A = AGRACIA
        am = (A["dp_p"] / 2) / math.sin(A["gm"] * DEG)
        for ch_, v_ in (("Am", am), ("a", A["a"]), ("bg", A["bg"]), ("al", A["al"]), ("gm", A["gm"]), ("T", A["T"]),
                        ("n1", A["n1"]), ("nc", A["nc"]), ("bp0", A["bp0"])):
            ceh[ch_].delete(0, "end"); ceh[ch_].insert(0, "%g" % round(v_, 3))
        calc_h()
        lbh.configure(text="Caso A. Gracia 2/60: Am do pinhão = (dp 62 / 2) / sen 8,9° = %.1f mm, com dp = diâmetro no ponto de contato. Geratriz da coroa 183,4 mm. Pressão 20° assumida; torque no pinhão." % am)
    pr = tk.Frame(pH, bg=BG); pr.pack(fill="x", padx=pd)
    for k_ in NAK:
        tk.Button(pr, text="Nakagawa " + k_, command=lambda k=k_: preset(k), bg=CARD, fg=ACC, activebackground=ACC,
                  font=fnt(12, True), relief="flat", bd=0, pady=int(7 * ESC)).pack(side="left", fill="x", expand=True, padx=2)
    botao(pH, "A.GRACIA 2 / 60", preset_ag, CARD, ACC)
    botao(pH, "CALCULAR HRH", calc_h)
    botao(pH, "ENVIAR F, u1, u2 AO CONTATO", enviar, CARD, TXT)
    lbh = tk.Label(pH, text="", bg=BG, fg=OK, font=fnt(10), wraplength=wc, justify="left"); lbh.pack()
    hres.pack(fill="x")
    for _ in range(4):
        cdh = cartao(pH)
        hcvs.append(tk.Canvas(cdh, width=wg, height=int(wg * 0.8), bg=BG, highlightthickness=0)); hcvs[-1].pack()

    # VALIDACAO
    cvl = tk.Frame(pV, bg=BG); cvl.pack(fill="x")
    def validar_ui():
        L = validar(); est["val"] = 1; est["v"] = texto_validacao(L); atualizar_txt(); limpar(cvl)
        n_ok = sum(1 for x in L if x[4]); nv = {c: sum(1 for x in L if x[5] == c) for c in "VRE"}
        bn = tk.Frame(cvl, bg=OK if n_ok == len(L) else BAD, pady=pd); bn.pack(fill="x", padx=pd, pady=pd // 2)
        tk.Label(bn, text="%d / %d aprovados" % (n_ok, len(L)), bg=bn["bg"], fg=BG, font=fnt(18, True)).pack()
        tk.Label(bn, text="V verificação %d · R referência %d · E dado externo %d" % (nv["V"], nv["R"], nv["E"]), bg=bn["bg"], fg=BG, font=fnt(11)).pack()
        cr = cartao(cvl, "Como ler", "Só E compara com dado externo. Tolerâncias até 21/10/2026 foram fixadas sem registro prévio (MAIN 28.26): tratar como posteriores. hmin angular é calibrado nos mesmos 9 pontos (MAIN 28.17 não satisfeito).")
        for n, v, r, e, ok, c, tl in L:
            rw = tk.Frame(cvl, bg=CARD, padx=int(8 * ESC), pady=int(5 * ESC)); rw.pack(fill="x", padx=pd, pady=1)
            tk.Label(rw, text="✔" if ok else "✘", bg=CARD, fg=OK if ok else BAD, font=fnt(16, True), width=2).pack(side="left")
            tk.Label(rw, text=c, bg=CARD, fg=ACC if c == "E" else MUT, font=fnt(12, True), width=2).pack(side="left")
            tk.Label(rw, text="%s\n%.4g  vs  %.4g   (erro %.2g%% · tol %.2g%%)" % (n, v, r, 100 * e, 100 * tl), bg=CARD, fg=TXT, font=fnt(10),
                     anchor="w", justify="left", wraplength=wg - int(60 * ESC)).pack(side="left", fill="x")
        kbp, ch = calibracoes()
        cc = cartao(cvl, "Calibrações registradas", "divisão calibra/valida (MAIN 27.20)")
        for ln in ("K_BP = %.4f: ajustado nos 3 pares do Nakagawa, validado na geometria do Paouris; não aplicado." % kbp,
                   "C_h = %.3f: ajustado nos pontos ímpares do Paouris, validado nos pares; não aplicado." % ch,
                   "Domínio do hmin calibrado: ue 1,5–8,7 m/s · F 230–10600 N · k 7–30 · arrastamento angular."):
            tk.Label(cc, text=ln, bg=CARD, fg=TXT, font=fnt(10), anchor="w", justify="left", wraplength=wg).pack(fill="x", pady=1)
        cmn = cartao(cvl, "Comparações sem aprovação")
        for ln in comparar_nakagawa()[1:]:
            tk.Label(cmn, text=ln, bg=CARD, fg=WARN, font=fnt(10), anchor="w", justify="left", wraplength=wg).pack(fill="x", pady=1)
    pV.pack_forget()
    botao(pV, "EXECUTAR VALIDAÇÃO", validar_ui)
    def validar_v11_ui():
        L = validar_v11(); n_ok = sum(1 for x in L if x[4])
        lbv.configure(text="v11: %d / %d aprovados (detalhe na aba TXT)" % (n_ok, len(L)), fg=OK if n_ok == len(L) else BAD)
        est["v"] = est.get("v", "") + "\n\nVALIDACAO v11\n" + "\n".join("%s [%s] %s: %.4g vs %.4g (erro %.3g, tol %.3g)" % ("OK " if x[4] else "ERR", x[5], x[0], x[1], x[2], x[3], x[6]) for x in L)
        atualizar_txt()
    botao(pV, "VALIDAÇÃO v11 (lenta)", validar_v11_ui, CARD, ACC)
    lbv = tk.Label(pV, text="", bg=BG, fg=OK, font=fnt(10), wraplength=wc, justify="left"); lbv.pack()
    cvl.pack_forget(); cvl.pack(fill="x")

    calc(False); calc_h(); ir("ENT")
    root.mainloop()

if __name__ == "__main__":
    if "--validar-v11" in sys.argv:
        import time as _t; _t0 = _t.time(); _L = validar_v11()
        for _n, _v, _r, _e, _ok, _c, _tl in _L:
            print("%s [%s] %-100s %11.4g vs %-11.4g (erro %.3g, tol %.3g)" % ("OK " if _ok else "ERR", _c, _n, _v, _r, _e, _tl))
        print("Resultado: %d/%d aprovados (%.0f s)" % (sum(1 for x in _L if x[4]), len(_L), _t.time() - _t0)); sys.exit(0 if all(x[4] for x in _L) else 1)
    if "--validar-hist" in sys.argv:
        import time as _t; _t0 = _t.time(); _L = validar_historico()
        for _n, _v, _r, _e, _ok, _c, _tl in _L:
            print("%s [%s] %-96s %11.4g vs %-11.4g (erro %.3g, tol %.3g)" % ("OK " if _ok else "ERR", _c, _n, _v, _r, _e, _tl))
        print("Resultado: %d/%d aprovados (%.0f s)" % (sum(1 for x in _L if x[4]), len(_L), _t.time() - _t0)); sys.exit(0 if all(x[4] for x in _L) else 1)
    if "--historico" in sys.argv:
        _P = dict(PAD); print(relatorio(calcular(_P, hist=True), None, True)); sys.exit(0)
    if "--validar" in sys.argv:
        print(texto_validacao(validar())); print(); _P = dict(PAD); print(relatorio(calcular(_P), incerteza(_P), "--completo" in sys.argv))
    else:
        interface()

# -*- coding: utf-8 -*-
"""
TRIBOCONTATO - contato lubrificado aco-aco, deslizamento rotativo (hipoide de alta reducao)
Pydroid 3 / Tkinter, so biblioteca padrao. Base: MAIN_TRIBOSOFT (secoes 9-11, 14, tabela-mestra).
Cadeia: Hertz eliptico exato -> filme EHL (Hamrock-Dowson) -> Coulomb -> flash (Blok)
        -> tensoes subsuperficiais -> desgaste (Archard) e vida (S-N) -> utilizacoes.
Uso: modo grafico (padrao) ou  python TRIBOCONTATO.py --validar
"""
import math, os, sys
PI = math.pi

# ---------------------------------------------------------------- PARAMETROS
# (chave, rotulo, unidade, padrao, fator->SI)
CAMPOS = [
 ("# GEOMETRIA (x = direcao de arrastamento; raio > 0 convexo, < 0 concavo)",),
 ("R1x", "Corpo 1  R1x", "mm", 22.0, 1e-3), ("R1y", "Corpo 1  R1y", "mm", 40.0, 1e-3),
 ("R2x", "Corpo 2  R2x", "mm", 115.0, 1e-3), ("R2y", "Corpo 2  R2y", "mm", -200.0, 1e-3),
 ("# CARGA E CINEMATICA (u1,u2 > 0, mesmo sentido)",),
 ("F", "Forca normal F", "N", 150.0, 1.0),
 ("u1", "Velocidade superficial u1", "m/s", 8.0, 1.0), ("u2", "Velocidade superficial u2", "m/s", 5.0, 1.0),
 ("N", "Ciclos de operacao N", "-", 1e6, 1.0),
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
]
PAD = {c[0]: c[3] * c[4] for c in CAMPOS if len(c) > 1}

# envelope MAIN (tabela-mestra): chave -> (min, max, rotulo)
ENV = {"ph": (500e6, 3000e6, "pH 500-3000 MPa"), "rred": (5e-3, 150e-3, "rho_red 5-150 mm"),
       "ue": (0, 13, "ue <= 13 m/s"), "vs": (0, 12, "v_desl <= 12 m/s"),
       "Ra": (0.3e-6, 0.5e-6, "Ra 0,3-0,5 um"), "nu40": (128, 168, "nu40 128-168 cSt"),
       "nu100": (14, 17, "nu100 14-17 cSt"), "alp": (1.3e-8, 2.8e-8, "alpha 1,3-2,8e-8 /Pa"),
       "Tb": (57, 90, "Tbulk 57-90 C"), "k": (1, 8, "k 1-8 (dominio Hamrock-Dowson)"),
       "Pe": (10, 1e99, "Pe >= 10 (dominio Blok)")}
# S-N: dois pontos do MAIN (secao 9 da tabela de valores gerais)
SN = ((800e6, 2.75e6), (700e6, 6.73e6))
SN_M = math.log(SN[1][1] / SN[0][1]) / math.log(SN[0][0] / SN[1][0])
N_LONGA = 5e7

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

def tensoes(a, b, p0, nu, nz=40, nt=60, nf=16):
    """Tensoes no eixo z (Boussinesq integrado sobre p = p0*sqrt(1-rho^2) na elipse).
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

def calcular(P, tens=True):
    R = dict(P)
    R["Rx"] = 1 / (1 / P["R1x"] + 1 / P["R2x"]); R["Ry"] = 1 / (1 / P["R1y"] + 1 / P["R2y"])
    R["rred"] = 1 / (1 / R["Rx"] + 1 / R["Ry"])
    R["Es"] = P["E"] / (2 * (1 - P["nu"] ** 2)); R["Ep"] = 2 * R["Es"]
    a, b, k, p0, pm = hertz(P["F"], R["Rx"], R["Ry"], R["Es"])
    R.update(a=a, b=b, k=k, p0=p0, pm=pm)
    R["ue"] = (P["u1"] + P["u2"]) / 2; R["vs"] = abs(P["u1"] - P["u2"])
    R["nuTb"] = walther(P["nu40"], P["nu100"], P["Tb"])
    R["eta0"] = R["nuTb"] * 1e-6 * P["rho"]
    hmin, hc, U, G, W = filme(P["F"], R["Rx"], k, R["Ep"], R["eta0"], R["ue"], P["alp"])
    R.update(hmin=hmin, hc=hc, U=U, G=G, W=W)
    R["Rq"] = 1.25 * math.hypot(P["Ra1"], P["Ra2"]); R["Lam"] = hmin / R["Rq"]
    R["Tf"], R["q0"], R["f1"] = flash(P["mu"], p0, R["vs"], a, P["u1"], P["u2"], P["ks"], P["rhos"], P["cps"])
    R["Tc"] = P["Tb"] + R["Tf"]
    kap = P["ks"] / (P["rhos"] * P["cps"])
    R["Pe"] = min(P["u1"], P["u2"]) * a / kap
    R["w1"] = P["K"] * pm * R["vs"] * 2 * a / P["u1"] / P["H"]       # por passagem
    R["w2"] = P["K"] * pm * R["vs"] * 2 * a / P["u2"] / P["H"]
    R["wN"] = P["N"] * max(R["w1"], R["w2"])
    R["Nf"] = SN[0][1] * (SN[0][0] / p0) ** SN_M
    R["D"] = P["N"] / R["Nf"]
    R["u_T"] = (R["Tc"] - P["Tb"]) / (P["Tlim"] - P["Tb"])
    R["u_h"] = P["hlim"] / hmin; R["u_w"] = R["wN"] / P["wlim"]; R["u_D"] = R["D"]
    if tens:
        zs, sx, sy, sz, sv, ta = tensoes(a, b, p0, P["nu"])
        i = sv.index(max(sv)); j = ta.index(max(ta))
        R.update(zs=zs, sx=sx, sy=sy, sz=sz, sv=sv, ta=ta, svm=sv[i], zsv=zs[i], tmax=ta[j], ztm=zs[j])
        R["u_s"] = R["svm"] / P["slim"]
    return R

def avisos(R):
    out = []
    def f(ok, msg):
        if not ok: out.append(msg)
    for ch, v in (("ph", R["p0"]), ("rred", R["rred"]), ("ue", R["ue"]), ("vs", R["vs"]),
                  ("nu40", R["nu40"]), ("nu100", R["nu100"]), ("alp", R["alp"]), ("Tb", R["Tb"]),
                  ("k", max(R["k"], 1 / R["k"])), ("Pe", R["Pe"])):
        f(ENV[ch][0] <= v <= ENV[ch][1], "fora do envelope MAIN: " + ENV[ch][2])
    for ch in ("Ra1", "Ra2"):
        f(ENV["Ra"][0] - 1e-12 <= R[ch] <= ENV["Ra"][1] + 1e-12, "fora do envelope MAIN: " + ENV["Ra"][2] + " (" + ch + ")")
    return out

# ---------------------------------------------------------------- VALIDACAO
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
        err = abs(val - ref) / abs(ref) if rel else abs(val - ref)
        L.append((nome, val, ref, err, err <= tol))
    P = dict(PAD); R = calcular(P)
    # 1 lubrificante
    t("Walther reproduz nu40", walther(148, 15.5, 40), 148, 1e-9)
    t("Walther reproduz nu100", walther(148, 15.5, 100), 15.5, 1e-9)
    t("eta40(128 cSt) vs MAIN 0,114", 128e-6 * 890, 0.114, 0.01)
    t("eta40(168 cSt) vs MAIN 0,150", 168e-6 * 890, 0.150, 0.01)
    t("eta100(14 cSt) vs MAIN 0,0125", 14e-6 * 890, 0.0125, 0.01)
    t("eta100(17 cSt) vs MAIN 0,015", 17e-6 * 890, 0.015, 0.01)
    # 2 material
    t("E* = E/(2(1-nu^2)) vs MAIN E'~115 GPa", R["Es"], 115e9, 0.03)
    # 3 Hertz circular (Johnson): a=(3FR/4E*)^(1/3), p0=(6FE*^2/(pi^3 R^2))^(1/3)
    Rc, F, Es = 0.02, 200.0, R["Es"]
    a, b, k, p0, pm = hertz(F, Rc, Rc, Es)             # Rx=Ry=Rc (esfera equivalente)
    t("Hertz circular: a", a, (3 * F * Rc / (4 * Es)) ** (1 / 3), 1e-9)
    t("Hertz circular: p0", p0, (6 * F * Es ** 2 / (PI ** 3 * Rc ** 2)) ** (1 / 3), 1e-9)
    # 4 Hertz eliptico: relacao exata e Hamrock-Brewe
    a, b, k, p0, pm = hertz(R["F"], R["Rx"], R["Ry"], Es)
    Kk, Ee = kelip(1 - 1 / k ** 2)
    t("Elipse: (k^2E-K)/(K-E) = Ry/Rx", (k * k * Ee - Kk) / (Kk - Ee), R["Ry"] / R["Rx"], 1e-9)
    t("Elipse: k exato vs Hamrock-Brewe 1,0339(Ry/Rx)^0,636", k, 1.0339 * (R["Ry"] / R["Rx"]) ** 0.636, 0.02)
    t("Equilibrio: F = int p dA = (2/3) pi a b p0", 2 * PI * a * b * p0 / 3, R["F"], 1e-12)
    # 5 condicao de Hertz independente: u(0,0)-u(x,0) = x^2/(2Rx) ; u(0,0)-u(0,y) = y^2/(2Ry)
    u0 = _u_hertz(0, 0, a, b, p0, Es)
    t("Hertz (deslocamento) em x=a/2", u0 - _u_hertz(a / 2, 0, a, b, p0, Es), (a / 2) ** 2 / (2 * R["Rx"]), 0.01)
    t("Hertz (deslocamento) em y=b/2", u0 - _u_hertz(0, b / 2, a, b, p0, Es), (b / 2) ** 2 / (2 * R["Ry"]), 0.01)
    # 6 tensoes: circular vs Johnson (nu=0,3)
    nu, a1, p1 = 0.3, 1e-3, 1e9
    zs, sx, sy, sz, sv, ta = tensoes(a1, a1, p1, nu)
    i = min(range(len(zs)), key=lambda q: abs(zs[q] - 0.5 * a1))
    z = zs[i]
    t("sigma_z(z) circular vs p0/(1+z^2/a^2)", sz[i], p1 / (1 + z * z / a1 ** 2), 0.005)
    sr = p1 * ((1 + nu) * (1 - z / a1 * math.atan(a1 / z)) - 1 / (2 * (1 + z * z / a1 ** 2)))
    t("sigma_r(z) circular vs Johnson", sx[i], sr, 0.01)
    t("tau_max circular vs 0,31 p0", max(ta), 0.31 * p1, 0.03)
    t("z(tau_max) circular vs 0,48 a", zs[ta.index(max(ta))], 0.48 * a1, 0.08)
    # 7 filme: expoentes HD e limite
    h1 = filme(R["F"], R["Rx"], k, R["Ep"], R["eta0"], R["ue"], P["alp"])[0]
    h2 = filme(R["F"], R["Rx"], k, R["Ep"], R["eta0"], 2 * R["ue"], P["alp"])[0]
    t("hmin ~ ue^0,68", h2 / h1, 2 ** 0.68, 1e-9)
    h3 = filme(2 * R["F"], R["Rx"], k, R["Ep"], R["eta0"], R["ue"], P["alp"])[0]
    t("hmin ~ F^-0,073", h3 / h1, 2 ** -0.073, 1e-9)
    # 8 flash: constante de Blok 1,11 (MAIN 11.6) e particao
    t("Blok: 2*jmax*sqrt(2)/pi^1,5 vs 1,11", 2 * jmax() * math.sqrt(2) / PI ** 1.5, 1.11, 0.005)
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
    # 10 potencia dissipada: int mu p vs dA = mu F vs
    t("Potencia: (2pi/3) a b p0 mu vs = mu F vs", 2 * PI * R["a"] * R["b"] * R["p0"] / 3 * 0.05 * 3.0, 0.05 * R["F"] * 3.0, 1e-12)
    return L

def texto_validacao(L):
    s = ["VALIDACAO MATEMATICA / MAIN"]
    for n, v, r, e, ok in L:
        s.append("%s %-52s %10.4g vs %-10.4g (%.2g%%)" % ("OK " if ok else "ERR", n, v, r, 100 * e))
    n_ok = sum(1 for x in L if x[4])
    s.append("Resultado: %d/%d testes aprovados" % (n_ok, len(L)))
    return "\n".join(s)

# ---------------------------------------------------------------- RELATORIO
def relatorio(R):
    g = lambda x: "%.4g" % x
    reg = "filme limite/rugoso (L<1)" if R["Lam"] < 1 else ("misto (1<=L<3)" if R["Lam"] < 3 else "filme pleno (L>=3)")
    uts = [("Temperatura", R["u_T"]), ("Filme minimo", R["u_h"]), ("Desgaste", R["u_w"]),
           ("von Mises", R.get("u_s", 0.0)), ("Fadiga (Miner)", R["u_D"])]
    crit = max(uts, key=lambda q: q[1])
    s = ["=== CONTATO ===",
         "Rx=%s mm  Ry=%s mm  rho_red=%s mm" % (g(R["Rx"] * 1e3), g(R["Ry"] * 1e3), g(R["rred"] * 1e3)),
         "E*=%s GPa  E'(EHL)=%s GPa" % (g(R["Es"] / 1e9), g(R["Ep"] / 1e9)),
         "a=%s mm (x)  b=%s mm (y)  k=b/a=%s" % (g(R["a"] * 1e3), g(R["b"] * 1e3), g(R["k"])),
         "p0=pH=%s MPa  pm=%s MPa" % (g(R["p0"] / 1e6), g(R["pm"] / 1e6)),
         "ue=%s m/s  v_desl=%s m/s" % (g(R["ue"]), g(R["vs"])),
         "=== FILME EHL ===",
         "nu(Tb)=%s cSt  eta0=%s Pa.s" % (g(R["nuTb"]), g(R["eta0"])),
         "U=%s  G=%s  W=%s" % (g(R["U"]), g(R["G"]), g(R["W"])),
         "hmin=%s um  hc=%s um" % (g(R["hmin"] * 1e6), g(R["hc"] * 1e6)),
         "Rq=%s um  Lambda=%s  -> %s" % (g(R["Rq"] * 1e6), g(R["Lam"]), reg),
         "=== ATRITO E TEMPERATURA ===",
         "q0=mu p0 v_desl=%s MW/m2" % g(R["q0"] / 1e6),
         "fracao do calor no corpo 1=%s" % g(R["f1"]),
         "T_flash=%s K  T_contato=%s C  Pe=%s" % (g(R["Tf"]), g(R["Tc"]), g(R["Pe"])),
         "=== DESGASTE E VIDA ===",
         "desgaste/passagem: corpo1=%s nm  corpo2=%s nm" % (g(R["w1"] * 1e9), g(R["w2"] * 1e9)),
         "desgaste em N=%s ciclos: %s um" % (g(R["N"]), g(R["wN"] * 1e6)),
         "Nf(S-N MAIN)=%s ciclos  D=N/Nf=%s%s" % (g(R["Nf"]), g(R["D"]), "  (vida longa)" if R["Nf"] > N_LONGA else "")]
    if "svm" in R:
        s += ["=== TENSOES (eixo z) ===",
              "sigma_VM max=%s MPa em z=%s mm" % (g(R["svm"] / 1e6), g(R["zsv"] * 1e3)),
              "tau_max=%s MPa em z=%s mm" % (g(R["tmax"] / 1e6), g(R["ztm"] * 1e3))]
    s.append("=== UTILIZACOES (<=1 admissivel) ===")
    s += ["%-16s %s" % (n, g(v)) for n, v in uts]
    s.append("Critico: %s (%s) -> %s" % (crit[0], g(crit[1]), "INTEGRO" if crit[1] <= 1 else "FORA DO LIMITE"))
    av = avisos(R)
    if av:
        s.append("=== AVISOS ==="); s += av
    return "\n".join(s)

# ---------------------------------------------------------------- TEMA / GRAFICOS
BG, CARD, HDR, TXT, MUT = "#0b1220", "#16213a", "#0f1a30", "#e2e8f0", "#94a3b8"
ACC, OK, WARN, BAD, PURP, GRID = "#38bdf8", "#4ade80", "#fbbf24", "#f87171", "#c084fc", "#2b3a55"
ESC = 1.0                                   # escala de tela (px por dp), ajustada na interface

def fnt(n, b=False):
    return ("Helvetica", -int(n * ESC), "bold") if b else ("Helvetica", -int(n * ESC))

def plot(cv, titulo, xl, yl, series, hl=()):
    cv.delete("all")
    W, H = int(cv["width"]), int(cv["height"])
    Lm, Rm, Tm, Bm = int(46 * ESC), int(12 * ESC), int(34 * ESC), int(44 * ESC)
    xs = [x for s in series for x in s[0]]
    ys = [y for s in series for y in s[1]] + [h[0] for h in hl]
    x0, x1 = min(xs), max(xs); y0, y1 = min(0.0, min(ys)), max(ys) * 1.08
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

def graficos(P, R, cvs):
    n = 61
    xa = [-R["a"] + 2 * R["a"] * i / (n - 1) for i in range(n)]
    yb = [-R["b"] + 2 * R["b"] * i / (n - 1) for i in range(n)]
    px_ = [R["p0"] / 1e9 * math.sqrt(max(1 - (x / R["a"]) ** 2, 0)) for x in xa]
    py_ = [R["p0"] / 1e9 * math.sqrt(max(1 - (y / R["b"]) ** 2, 0)) for y in yb]
    plot(cvs[0], "Pressão de Hertz", "posição (mm)", "p (GPa)",
         [([x * 1e3 for x in xa], px_, BAD, "p(x,0)"), ([y * 1e3 for y in yb], py_, ACC, "p(0,y)")])
    z = [q * 1e3 for q in R["zs"]]
    plot(cvs[1], "Tensões no eixo z", "profundidade (mm)", "tensão (GPa)",
         [(z, [v / 1e9 for v in R["sv"]], BAD, "von Mises"), (z, [v / 1e9 for v in R["ta"]], OK, "τ máx"),
          (z, [v / 1e9 for v in R["sz"]], ACC, "σ z")], [(P["slim"] / 1e9, PURP, "limite VM")])
    xs, hs = [], []
    for i in range(15):
        s = 0.2 + 1.4 * i / 14
        r = calcular(dict(P, u1=P["u1"] * s, u2=P["u2"] * s), False)
        xs.append(r["ue"]); hs.append(r["hmin"] * 1e6)
    plot(cvs[2], "Filme mínimo × velocidade", "ue (m/s)  [escala u1,u2]", "hmin (µm)",
         [(xs, hs, ACC, "hmin")], [(R["Rq"] * 1e6, BAD, "Λ=1"), (P["hlim"] * 1e6, PURP, "hlim")])
    fs, tc = [], []
    for i in range(15):
        s = 0.3 + 1.4 * i / 14
        fs.append(P["F"] * s); tc.append(calcular(dict(P, F=P["F"] * s), False)["Tc"])
    plot(cvs[3], "Temperatura de contato × carga", "F (N)", "Tc (°C)",
         [(fs, tc, BAD, "Tc")], [(P["Tlim"], PURP, "Tlim")])

# ---------------------------------------------------------------- INTERFACE
def interface():
    import tkinter as tk
    global ESC
    root = tk.Tk(); root.title("TRIBOCONTATO"); root.configure(bg=BG)
    sw = root.winfo_screenwidth()
    ESC = max(1.0, sw / 360.0)
    pd = int(8 * ESC); wc = sw - 4 * pd
    ent, est, pag, abas = {}, {}, {}, {}

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
        b = tk.Button(tb, text=rot, command=lambda: ir(nome), font=fnt(11, True), relief="flat", bd=0,
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
    pE = nova("ENT", "ENTRADA"); pR = nova("RES", "RESULTADO"); pG = nova("GRA", "GRÁFICOS")
    pV = nova("VAL", "VALIDAÇÃO"); pT = nova("TXT", "RELATÓRIO")

    # ENTRADA
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

    def padroes():
        for c in CAMPOS:
            if len(c) > 1:
                ent[c[0]].delete(0, "end"); ent[c[0]].insert(0, "%g" % c[3])
    botao(pE, "CALCULAR", lambda: calc(True))
    botao(pE, "Restaurar padrões", padroes, CARD, TXT)

    def ler():
        P = {}
        for c in CAMPOS:
            if len(c) > 1:
                P[c[0]] = float(ent[c[0]].get().replace(",", ".")) * c[4]
        for ch in ("u1", "u2", "F", "E", "H"):
            if P[ch] <= 0: raise ValueError(ch + " deve ser > 0")
        return P

    # GRAFICOS
    wg = wc - 2 * int(10 * ESC) - 2 * pd
    cvs = []
    for _ in range(4):
        cd = cartao(pG)
        cvs.append(tk.Canvas(cd, width=wg, height=int(wg * 0.8), bg=BG, highlightthickness=0)); cvs[-1].pack()

    # RELATORIO
    cdT = cartao(pT)
    txt = tk.Text(cdT, width=max(20, int(wg / (6 * ESC))), height=48, font=("Courier", -int(10 * ESC)), wrap="word",
                  bg=BG, fg=TXT, relief="flat", insertbackground=TXT)
    txt.pack(fill="both")
    def atualizar_txt():
        txt.delete("1.0", "end"); txt.insert("end", est.get("s", "") + "\n\n" + est.get("v", ""))
    def salvar():
        pasta = "/storage/emulated/0/Download"
        if not os.path.isdir(pasta): pasta = os.getcwd()
        cam = os.path.join(pasta, "TRIBOCONTATO_resultado.txt")
        with open(cam, "w", encoding="utf-8") as f:
            f.write(est.get("s", "") + "\n\n" + est.get("v", ""))
        lb.configure(text="Salvo em " + cam)
    botao(pT, "SALVAR TXT", salvar)
    lb = tk.Label(pT, text="", bg=BG, fg=OK, font=fnt(10), wraplength=wc); lb.pack()

    # RESULTADO
    def limpar(f):
        for w in f.winfo_children(): w.destroy()
    def cor_u(u): return OK if u < 0.7 else (WARN if u <= 1 else BAD)

    def mostrar(R):
        limpar(pR)
        uts = [("Temperatura", R["u_T"]), ("Filme mínimo", R["u_h"]), ("Desgaste", R["u_w"]),
               ("von Mises", R["u_s"]), ("Fadiga (Miner)", R["u_D"])]
        crit = max(uts, key=lambda q: q[1]); okk = crit[1] <= 1
        bn = tk.Frame(pR, bg=OK if okk else BAD, padx=pd, pady=pd); bn.pack(fill="x", padx=pd, pady=pd // 2)
        tk.Label(bn, text="ÍNTEGRO" if okk else "FORA DO LIMITE", bg=bn["bg"], fg=BG, font=fnt(20, True)).pack()
        tk.Label(bn, text="crítico: %s (utilização %.2f)" % (crit[0], crit[1]), bg=bn["bg"], fg=BG, font=fnt(12)).pack()
        lam = R["Lam"]
        reg = "limite" if lam < 1 else ("misto" if lam < 3 else "filme pleno")
        tiles = [("Pressão pH", "%.0f" % (R["p0"] / 1e6), "MPa", TXT), ("Elipse k = b/a", "%.2f" % R["k"], "", TXT),
                 ("Semi-eixos a × b", "%.3f × %.3f" % (R["a"] * 1e3, R["b"] * 1e3), "mm", TXT),
                 ("Filme mínimo", "%.3f" % (R["hmin"] * 1e6), "µm", cor_u(R["u_h"])),
                 ("Λ = hmin/Rq", "%.2f" % lam, reg, OK if lam >= 3 else (WARN if lam >= 1 else BAD)),
                 ("Temp. contato", "%.1f" % R["Tc"], "°C", cor_u(R["u_T"])),
                 ("Flash de Blok", "%.1f" % R["Tf"], "K", TXT), ("Viscosidade η0", "%.4f" % R["eta0"], "Pa·s", TXT),
                 ("Vida Nf", "%.2e" % R["Nf"], "ciclos", cor_u(R["u_D"])),
                 ("Desgaste em N", "%.3f" % (R["wN"] * 1e6), "µm", cor_u(R["u_w"])),
                 ("von Mises máx", "%.0f" % (R["svm"] / 1e6), "MPa", cor_u(R["u_s"])),
                 ("τ máx", "%.0f" % (R["tmax"] / 1e6), "MPa", TXT)]
        gr = tk.Frame(pR, bg=BG); gr.pack(fill="x", padx=pd // 2)
        gr.columnconfigure(0, weight=1, uniform="t"); gr.columnconfigure(1, weight=1, uniform="t")
        for i, (t, v, u, cr) in enumerate(tiles):
            tl = tk.Frame(gr, bg=CARD, padx=int(8 * ESC), pady=int(7 * ESC))
            tl.grid(row=i // 2, column=i % 2, sticky="nsew", padx=pd // 2, pady=pd // 2)
            tk.Label(tl, text=t, bg=CARD, fg=MUT, font=fnt(10), anchor="w").pack(fill="x")
            tk.Label(tl, text=v, bg=CARD, fg=cr, font=fnt(18, True), anchor="w").pack(fill="x")
            tk.Label(tl, text=u or " ", bg=CARD, fg=MUT, font=fnt(10), anchor="w").pack(fill="x")
        cd = cartao(pR, "Utilizações", "≤ 1 admissível · linha branca = limite")
        wb = wg; hb = int(30 * ESC) * len(uts) + int(8 * ESC)
        cb = tk.Canvas(cd, width=wb, height=hb, bg=CARD, highlightthickness=0); cb.pack()
        x0, x1 = int(wb * 0.34), wb - int(8 * ESC); sc = (x1 - x0) / 1.5
        for i, (n, u) in enumerate(uts):
            y = int(30 * ESC) * i + int(15 * ESC)
            cb.create_text(0, y, text=n, anchor="w", fill=TXT, font=fnt(11))
            cb.create_rectangle(x0, y - int(8 * ESC), x1, y + int(8 * ESC), fill=BG, outline="")
            cb.create_rectangle(x0, y - int(8 * ESC), x0 + min(u, 1.5) * sc, y + int(8 * ESC), fill=cor_u(u), outline="")
            cb.create_text(x1, y - int(14 * ESC), text="%.2f" % u, anchor="e", fill=MUT, font=fnt(9))
        cb.create_line(x0 + sc, 0, x0 + sc, hb, fill="#ffffff", width=2)
        av = avisos(R)
        if av:
            ca = cartao(pR, "Avisos")
            for a in av:
                tk.Label(ca, text="• " + a, bg=CARD, fg=WARN, font=fnt(11), anchor="w", justify="left", wraplength=wg).pack(fill="x")

    def calc(vai):
        try:
            P = ler(); R = calcular(P)
            est["s"] = relatorio(R); mostrar(R); graficos(P, R, cvs); atualizar_txt()
            if vai: ir("RES")
        except Exception as ex:
            limpar(pR)
            tk.Label(pR, text="ERRO: %s" % ex, bg=BG, fg=BAD, font=fnt(13), wraplength=wc, justify="left").pack(padx=pd, pady=pd)
            if vai: ir("RES")

    # VALIDACAO
    botao(pV, "EXECUTAR VALIDAÇÃO", lambda: validar_ui())
    cvl = tk.Frame(pV, bg=BG); cvl.pack(fill="x")
    def validar_ui():
        L = validar(); est["val"] = 1; est["v"] = texto_validacao(L); atualizar_txt(); limpar(cvl)
        n_ok = sum(1 for x in L if x[4])
        bn = tk.Frame(cvl, bg=OK if n_ok == len(L) else BAD, pady=pd); bn.pack(fill="x", padx=pd, pady=pd // 2)
        tk.Label(bn, text="%d / %d testes aprovados" % (n_ok, len(L)), bg=bn["bg"], fg=BG, font=fnt(18, True)).pack()
        for n, v, r, e, ok in L:
            rw = tk.Frame(cvl, bg=CARD, padx=int(8 * ESC), pady=int(5 * ESC)); rw.pack(fill="x", padx=pd, pady=1)
            tk.Label(rw, text="✔" if ok else "✘", bg=CARD, fg=OK if ok else BAD, font=fnt(16, True), width=2).pack(side="left")
            tk.Label(rw, text="%s\n%.4g  vs  %.4g   (erro %.2g%%)" % (n, v, r, 100 * e), bg=CARD, fg=TXT, font=fnt(10),
                     anchor="w", justify="left", wraplength=wg - int(40 * ESC)).pack(side="left", fill="x")

    calc(False); ir("ENT")
    root.mainloop()

if __name__ == "__main__":
    if "--validar" in sys.argv:
        print(texto_validacao(validar())); print(); print(relatorio(calcular(PAD)))
    else:
        interface()

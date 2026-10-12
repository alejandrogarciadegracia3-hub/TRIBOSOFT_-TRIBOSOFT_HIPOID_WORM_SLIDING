# -*- coding: utf-8 -*-
"""
MODELO TRIBOLÓGICO INTEGRADO - V7.81
Pydroid 3 / Android / Tkinter

V0.2:
- superfícies paramétricas gerais r1(u,v,theta), r2(xi,eta,theta)
- derivadas numéricas das superfícies
- normais
- velocidade relativa por evolução em theta
- teste de coincidência geométrica
- teste cinemático v_rel . n = 0
- identificação numérica de candidatos à linha de contato
- construção da representação do volume tribológico em torno dos candidatos

IMPORTANTE:
Esta versão NÃO inventa uma geometria física específica.
A geometria deve ser fornecida pelo usuário como seis expressões:
x1(u,v,theta), y1(...), z1(...)
x2(xi,eta,theta), y2(...), z2(...)

Convenção:
- u, v, xi, eta: parâmetros geométricos
- theta: posição cinemática
- theta_dot: velocidade angular de theta [rad/s]
- coordenadas internas: metros
"""

import math
import os
import threading
import time
import tkinter as tk
from tkinter import font as tkfont


# ============================================================
# V7.77 — INTERFACE MOBILE CONSOLIDADA
# ============================================================


# ============================================================
# ERROS
# ============================================================

class MissingInputError(Exception):
    pass


class ValidationError(Exception):
    pass


# ============================================================
# AMBIENTE MATEMÁTICO CONTROLADO
# ============================================================

SAFE_MATH = {
    "sin": math.sin,
    "cos": math.cos,
    "tan": math.tan,
    "asin": math.asin,
    "acos": math.acos,
    "atan": math.atan,
    "atan2": math.atan2,
    "sqrt": math.sqrt,
    "exp": math.exp,
    "log": math.log,
    "pi": math.pi,
    "abs": abs,
    "pow": pow,
}


# ============================================================
# VETORES
# ============================================================


# ============================================================
# SUPERFÍCIE PARAMÉTRICA
# ============================================================


# ============================================================
# PONTO AMOSTRADO
# ============================================================


# ============================================================
# CONTATO CINEMÁTICO
# ============================================================


# ============================================================
# AMOSTRAGEM
# ============================================================


# ============================================================
# BUSCA DE CANDIDATOS
# ============================================================


# ============================================================
# VOLUME TRIBOLÓGICO
# ============================================================


# ============================================================
# INTERFACE
# ============================================================
# CONTINUAÇÃO DISCRETA DA LINHA DE CONTATO
# ============================================================


# ============================================================
# PRESSÃO DE CONTATO — V0.5
# ============================================================


# ============================================================
# CAMPO DE TENSÃO SUBSUPERFICIAL — V0.6
# ============================================================


# ============================================================
# TENSOR DE TENSÕES — V0.7
# ============================================================


def von_mises_from_tensor(t):
    """
    Von Mises para tensor cartesiano:

        sigma_VM =
        sqrt(
          0.5[
            (sx-sy)^2 +
            (sy-sz)^2 +
            (sz-sx)^2
            + 6(txy²+txz²+tyz²)
          ]
        )
    """
    sx = t["sx"]
    sy = t["sy"]
    sz = t["sz"]

    txy = t["txy"]
    txz = t["txz"]
    tyz = t["tyz"]

    value = 0.5 * (
        (sx - sy)**2
        + (sy - sz)**2
        + (sz - sx)**2
        + 6.0 * (
            txy*txy
            + txz*txz
            + tyz*tyz
        )
    )

    return math.sqrt(
        max(0.0, value)
    )


# ============================================================
# TENSOR ESPACIAL COMPLETO — V0.9
# ============================================================

def boussinesq_point_tensor(
    x,
    y,
    z,
    P,
    nu
):
    """
    Tensor cartesiano de Boussinesq para uma carga pontual
    normal à superfície de um semi-espaço elástico isotrópico.

    Convenção:
      z > 0: profundidade
      P > 0: magnitude da carga normal
      tensões normais compressivas: sinal negativo

    Formulação cartesiana:

      R = sqrt(x²+y²+z²)
      r² = x²+y²

    Para r > 0:

      sx = P/(2π) [
             (1-2ν)/r² *
             ((1-z/R)(x²-y²)/r² + z y²/R³)
             - 3 z x²/R⁵
           ]

      sy = P/(2π) [
             (1-2ν)/r² *
             ((1-z/R)(y²-x²)/r² + z x²/R³)
             - 3 z y²/R⁵
           ]

      sz = -3P z³/(2πR⁵)

      txy = P/(2π) [
              (1-2ν)/r² *
              (2xy(1-z/R)/r² - xyz/R³)
              - 3xyz/R⁵
            ]

      txz = -3P x z²/(2πR⁵)

      tyz = -3P y z²/(2πR⁵)

    No eixo r=0 são usados os limites analíticos:
      sx = sy = P(1-2ν)/(4πz²)
      sz = -3P/(2πz²)
      txy = txz = tyz = 0.

    Referência da formulação:
    Boussinesq para carga normal em semi-espaço elástico.
    """

    if z <= 0.0:
        raise ValidationError(
            "A profundidade z deve ser > 0."
        )

    if P < 0.0:
        raise ValidationError(
            "P deve ser >= 0 nesta convenção."
        )

    if not (0.0 <= nu < 0.5):
        raise ValidationError(
            "Poisson nu deve satisfazer 0 <= nu < 0.5."
        )

    R2 = x*x + y*y + z*z
    R = math.sqrt(R2)
    r2 = x*x + y*y

    c = P / (2.0 * math.pi)

    if r2 < 1.0e-30:
        sx = P * (1.0 - 2.0*nu) / (
            4.0 * math.pi * z*z
        )
        sy = sx
        sz = -3.0 * P / (
            2.0 * math.pi * z*z
        )

        return {
            "sx": sx,
            "sy": sy,
            "sz": sz,
            "txy": 0.0,
            "txz": 0.0,
            "tyz": 0.0,
        }

    R3 = R2 * R
    R5 = R3 * R2

    one_minus = 1.0 - 2.0*nu

    sx = c * (
        one_minus / r2 * (
            (1.0 - z/R) * (x*x - y*y)/r2
            + z*y*y/R3
        )
        - 3.0*z*x*x/R5
    )

    sy = c * (
        one_minus / r2 * (
            (1.0 - z/R) * (y*y - x*x)/r2
            + z*x*x/R3
        )
        - 3.0*z*y*y/R5
    )

    sz = -3.0 * c * z**3 / R5

    txy = c * (
        one_minus / r2 * (
            2.0*x*y*(1.0 - z/R)/r2
            - x*y*z/R3
        )
        - 3.0*x*y*z/R5
    )

    txz = -3.0 * c * x*z*z / R5
    tyz = -3.0 * c * y*z*z / R5

    return {
        "sx": sx,
        "sy": sy,
        "sz": sz,
        "txy": txy,
        "txz": txz,
        "tyz": tyz,
    }


def boussinesq_distributed_tensor(
    x0, y0, z, a, b, p0, nu, nx=40
):
    """
    Superposição numérica de Boussinesq em coordenadas elípticas.

    r = sin(theta)
    x = a*r*cos(phi)
    y = b*r*sin(phi)

    dA = a*b*r*cos(theta) dtheta dphi
    p  = p0*cos(theta)

    Mesma soma de boussinesq_point_tensor (mesma ordem de operações),
    com os invariantes do anel e o laço interno inlinados: sem dict
    nem chamada de função por ponto.
    """
    if z <= 0.0:
        raise ValidationError("A profundidade z deve ser > 0.")

    if not (0.0 <= nu < 0.5):
        raise ValidationError(
            "Poisson nu deve satisfazer 0 <= nu < 0.5."
        )

    ntheta = max(8, int(nx))
    if ntheta % 2:
        ntheta += 1
    nphi = 2 * ntheta

    dtheta = (math.pi / 2.0) / ntheta
    dphi = 2.0 * math.pi / nphi

    cos_phi = [math.cos((ip + 0.5) * dphi) for ip in range(nphi)]
    sin_phi = [math.sin((ip + 0.5) * dphi) for ip in range(nphi)]

    one_minus = 1.0 - 2.0*nu
    z3 = z**3
    z_3 = 3.0*z

    sx_t = sy_t = sz_t = 0.0
    txy_t = txz_t = tyz_t = 0.0

    for it in range(ntheta):
        theta = (it + 0.5) * dtheta
        r = math.sin(theta)
        c = math.cos(theta)
        p = p0 * c
        dA = a * b * r * c * dtheta * dphi
        dP = p * dA

        if dP < 0.0:
            raise ValidationError(
                "P deve ser >= 0 nesta convenção."
            )

        cc = dP / (2.0 * math.pi)
        m3c = -3.0 * cc
        a_r = a * r
        b_r = b * r

        for ip in range(nphi):
            x = a_r * cos_phi[ip] - x0
            y = b_r * sin_phi[ip] - y0

            r2 = x*x + y*y

            if r2 < 1.0e-30:
                sx = dP * one_minus / (
                    4.0 * math.pi * z*z
                )
                sx_t += sx
                sy_t += sx
                sz_t += -3.0 * dP / (
                    2.0 * math.pi * z*z
                )
                continue

            R2 = r2 + z*z
            R = math.sqrt(R2)
            R3 = R2 * R
            R5 = R3 * R2

            sx_t += cc * (
                one_minus / r2 * (
                    (1.0 - z/R) * (x*x - y*y)/r2
                    + z*y*y/R3
                )
                - z_3*x*x/R5
            )

            sy_t += cc * (
                one_minus / r2 * (
                    (1.0 - z/R) * (y*y - x*x)/r2
                    + z*x*x/R3
                )
                - z_3*y*y/R5
            )

            sz_t += m3c * z3 / R5

            txy_t += cc * (
                one_minus / r2 * (
                    2.0*x*y*(1.0 - z/R)/r2
                    - x*y*z/R3
                )
                - 3.0*x*y*z/R5
            )

            txz_t += m3c * x*z*z / R5
            tyz_t += m3c * y*z*z / R5

    return {
        "sx": sx_t,
        "sy": sy_t,
        "sz": sz_t,
        "txy": txy_t,
        "txz": txz_t,
        "tyz": tyz_t,
    }


def spatial_full_tensor_field(
    a,
    b,
    p0,
    nu,
    zmax,
    nx_xy,
    nz,
    integration_grid,
    progress=None
):
    """
    Campo completo:
        sigma_x, sigma_y, sigma_z,
        tau_xy, tau_xz, tau_yz,
        sigma_VM

    em uma malha discreta (x,y,z).
    """

    nx_xy = max(3, int(nx_xy))
    nz = max(2, int(nz))

    field = []

    for ix in range(nx_xy):
        x0 = -a + 2.0*a*ix/(nx_xy - 1)

        for iy in range(nx_xy):
            y0 = -b + 2.0*b*iy/(nx_xy - 1)

            for iz in range(1, nz + 1):
                z = zmax * iz / nz

                tensor = boussinesq_distributed_tensor(
                    x0,
                    y0,
                    z,
                    a,
                    b,
                    p0,
                    nu,
                    integration_grid
                )

                vm = von_mises_from_tensor(
                    tensor
                )

                field.append({
                    "x": x0,
                    "y": y0,
                    "z": z,
                    "tensor": tensor,
                    "vm": vm,
                })

            if progress is not None:
                progress(ix * nx_xy + iy + 1, nx_xy * nx_xy)

    return field


_V782_FIELD_CACHE = {}


def spatial_full_tensor_field_cached(
    a,
    b,
    p0,
    nu,
    zmax,
    nx_xy,
    nz,
    integration_grid,
    progress=None
):
    """
    Mesmo campo de spatial_full_tensor_field.

    As tensões de Boussinesq são lineares em p0 e invariantes de escala
    em a: sigma/p0 em (x/a, y/a, z/a) depende só de b/a, nu e da malha.
    O campo normalizado é calculado uma vez e reaproveitado; a e p0 apenas
    reescalam. Mudar F, uₑ, v_slide, T0 etc. não recalcula o campo.
    """
    key = (
        b / a, float(nu), zmax / a,
        int(nx_xy), int(nz), int(integration_grid)
    )
    unit = _V782_FIELD_CACHE.get(key)

    if unit is None:
        unit = spatial_full_tensor_field(
            1.0, b / a, 1.0, nu, zmax / a,
            nx_xy, nz, integration_grid, progress
        )

        if len(_V782_FIELD_CACHE) >= 8:
            _V782_FIELD_CACHE.clear()

        _V782_FIELD_CACHE[key] = unit

    elif progress is not None:
        progress(1, 1)

    return [
        {
            "x": q["x"] * a,
            "y": q["y"] * a,
            "z": q["z"] * a,
            "tensor": {k: v * p0 for k, v in q["tensor"].items()},
            "vm": q["vm"] * p0,
        }
        for q in unit
    ]


def full_tensor_extrema(field):
    if not field:
        return None, None

    vm_item = max(
        field,
        key=lambda q: q["vm"]
    )

    comp_item = min(
        field,
        key=lambda q: q["tensor"]["sz"]
    )

    return vm_item, comp_item


# ============================================================
# CAMPO ESPACIAL SUBSUPERFICIAL — V0.8
# ============================================================


# ============================================================
# INTEGRIDADE — V1.0
# ============================================================


# ============================================================
# CRITÉRIO DE FADIGA — V1.1
# ============================================================


# ============================================================
# DANO ACUMULADO E VIDA — V1.2
# ============================================================

def sn_cycles_from_stress(
    sigma_a,
    sigma_ref,
    N_ref,
    exponent
):
    """
    Curva S-N de potência, explicitamente fornecida pelo usuário:

        sigma_a / sigma_ref =
        (N / N_ref)^(-exponent)

    Portanto:

        N =
        N_ref *
        (sigma_a / sigma_ref)^(-1/exponent)

    Esta função não inventa uma curva S-N.
    Todos os parâmetros são entradas do caso.
    """
    if sigma_a <= 0:
        return float("inf")

    if sigma_ref <= 0:
        raise ValidationError(
            "sigma_ref deve ser > 0."
        )

    if N_ref <= 0:
        raise ValidationError(
            "N_ref deve ser > 0."
        )

    if exponent <= 0:
        raise ValidationError(
            "O expoente S-N deve ser > 0."
        )

    ratio = sigma_a / sigma_ref

    return (
        N_ref
        * ratio ** (-1.0 / exponent)
    )


def life_from_field(
    field,
    sigma_ref,
    N_ref,
    exponent
):
    """
    Para cada ponto do campo:

        sigma_a := sigma_VM

        N_f(x) = curva S-N
        D = N / N_f

    O ponto crítico é aquele com menor N_f.
    """
    result = []

    for item in field:
        sigma_a = item["vm"]

        N_f = sn_cycles_from_stress(
            sigma_a,
            sigma_ref,
            N_ref,
            exponent
        )

        result.append({
            "x": item["x"],
            "y": item["y"],
            "z": item["z"],
            "sigma_a": sigma_a,
            "N_f": N_f,
        })

    return result


def life_extrema(result):
    if not result:
        return None

    return min(
        result,
        key=lambda q: q["N_f"]
    )


# ============================================================
# V7.55 — E1 → E2 → CURVATURA RELATIVA → EHL
# ============================================================


# ============================================================
# V7.54 — CURVATURA RELATIVA NO PONTO DE CONTATO
# ============================================================


# ============================================================
# V7.53 — CURVATURAS PRINCIPAIS A PARTIR DE E1
# ============================================================


# ============================================================
# V7.52 — RAIOS DE CURVATURA EQUIVALENTES EHL
# ============================================================


# ============================================================
# V7.51 — E' REDUZIDO / EQUIVALENTE
# ============================================================


# ============================================================
# EHL — V1.3
# ============================================================

def ehl_dimensionless_groups(
    eta0,
    ue,
    Eprime,
    Rx,
    force,
    alpha_p
):
    if eta0 <= 0:
        raise ValidationError("eta0 deve ser > 0.")
    if ue < 0:
        raise ValidationError("u_e deve ser >= 0.")
    if Eprime <= 0:
        raise ValidationError("E' deve ser > 0.")
    if Rx <= 0:
        raise ValidationError("Rx deve ser > 0.")
    if force <= 0:
        raise ValidationError("F deve ser > 0.")
    if alpha_p < 0:
        raise ValidationError("alpha_p deve ser >= 0.")

    U = eta0 * ue / (Eprime * Rx)
    G = alpha_p * Eprime
    W = force / (Eprime * Rx * Rx)

    return U, G, W


def ehl_film_thickness(
    eta0,
    ue,
    Eprime,
    Rx,
    Ry,
    force,
    alpha_p
):
    """
    Correlações EHL elípticas já definidas no modelo:

      U = eta0*ue/(E'*Rx)
      G = alpha_p*E'
      W = F/(E'*Rx²)
      k = 1.03*(Ry/Rx)^0.64

      hmin/Rx =
      3.63 U^0.68 G^0.49 W^-0.073
      (1-exp(-0.68k))

      hc/Rx =
      2.69 U^0.67 G^0.53 W^-0.067
      (1-0.61 exp(-0.73k))
    """
    if Ry <= 0:
        raise ValidationError("Ry deve ser > 0.")

    U, G, W = ehl_dimensionless_groups(
        eta0,
        ue,
        Eprime,
        Rx,
        force,
        alpha_p
    )

    # Ellipticity parameter from the Hamrock-Dowson point-contact correlation.
    k = 1.03 * (Ry / Rx) ** 0.64

    hmin = (
        Rx
        * 3.63
        * U**0.68
        * G**0.49
        * W**(-0.073)
        * (1.0 - math.exp(-0.68*k))
    )

    hc = (
        Rx
        * 2.69
        * U**0.67
        * G**0.53
        * W**(-0.067)
        * (1.0 - 0.61*math.exp(-0.73*k))
    )

    return {
        "U": U,
        "G": G,
        "W": W,
        "k": k,
        "hmin": hmin,
        "hc": hc,
    }


def film_parameter(
    hmin,
    Rq1,
    Rq2
):
    if hmin < 0:
        raise ValidationError(
            "hmin deve ser >= 0."
        )

    if Rq1 < 0 or Rq2 < 0:
        raise ValidationError(
            "Rq deve ser >= 0."
        )

    Rq = math.sqrt(
        Rq1*Rq1 + Rq2*Rq2
    )

    if Rq == 0:
        return float("inf"), Rq

    return hmin / Rq, Rq


# ============================================================
# ATRITO E TÉRMICA — V1.4
# ============================================================


def flash_temperature_rise(
    heat_flux,
    thermal_conductivity,
    thermal_diffusivity,
    contact_time
):
    """
    Modelo térmico local simplificado para elevação de
    temperatura por fonte superficial transitória:

        DeltaT_flash =
        q / k * sqrt(alpha * t / pi)

    Este módulo não substitui uma solução térmica completa.
    """
    if heat_flux < 0:
        raise ValidationError("q deve ser >= 0.")
    if thermal_conductivity <= 0:
        raise ValidationError("k térmica deve ser > 0.")
    if thermal_diffusivity <= 0:
        raise ValidationError("alpha térmica deve ser > 0.")
    if contact_time < 0:
        raise ValidationError("tempo de contato deve ser >= 0.")

    return (
        heat_flux
        / thermal_conductivity
        * math.sqrt(
            thermal_diffusivity
            * contact_time
            / math.pi
        )
    )


def bulk_temperature_rise(
    friction_power,
    thermal_conductance
):
    """
    Elevação de temperatura bulk em modelo concentrado:

        DeltaT_bulk = P_fric / H

    H é a condutância térmica equivalente do sistema.
    """
    if friction_power < 0:
        raise ValidationError("P_fric deve ser >= 0.")
    if thermal_conductance <= 0:
        raise ValidationError(
            "condutância térmica deve ser > 0."
        )

    return friction_power / thermal_conductance


def thermal_state(
    T0,
    heat_flux,
    thermal_conductivity,
    thermal_diffusivity,
    contact_time,
    friction_power,
    thermal_conductance
):
    dT_flash = flash_temperature_rise(
        heat_flux,
        thermal_conductivity,
        thermal_diffusivity,
        contact_time
    )

    dT_bulk = bulk_temperature_rise(
        friction_power,
        thermal_conductance
    )

    T_flash = T0 + dT_flash
    T_bulk = T0 + dT_bulk
    T_total = T0 + dT_bulk + dT_flash

    return {
        "dT_flash": dT_flash,
        "dT_bulk": dT_bulk,
        "T_flash": T_flash,
        "T_bulk": T_bulk,
        "T_total": T_total,
    }


# ============================================================
# ACOPLAMENTO TÉRMICO-EHL — V1.5
# ============================================================


# ============================================================
# DESGASTE — V7.78 / CARGA DE ASPEREZAS EM EHL MISTO
# ============================================================

def point_contact_asperity_load_ratio(
    dimensionless_load,
    dimensionless_speed,
    dimensionless_material,
    dimensionless_roughness,
    dimensionless_hardness,
):
    """
    Correlação de Masjedi & Khonsari para contato pontual/ellíptico:

        L_a = 10 W^(-0.083) U^(0.143) G^(0.314)
              ln[1 + σ̄^(4.689) V^(0.509)
                     W^(-0.501) U^(-2.90) G^(-2.870)]

    L_a é a razão de carga sustentada pelas asperezas, expressa em %.

    Referência:
        Masjedi & Khonsari, Tribology International 82 (2015) 228–244,
        Eq. (23), DOI 10.1016/j.triboint.2014.09.010.
    """
    vals = [
        dimensionless_load,
        dimensionless_speed,
        dimensionless_material,
        dimensionless_roughness,
        dimensionless_hardness,
    ]
    if any(v <= 0 for v in vals):
        raise ValidationError(
            "W, U, G, sigma_bar e V devem ser > 0 para calcular L_a."
        )

    term = (
        1.0
        + dimensionless_roughness ** 4.689
        * dimensionless_hardness ** 0.509
        * dimensionless_load ** (-0.501)
        * dimensionless_speed ** (-2.90)
        * dimensionless_material ** (-2.870)
    )

    La_percent = (
        10.0
        * dimensionless_load ** (-0.083)
        * dimensionless_speed ** 0.143
        * dimensionless_material ** 0.314
        * math.log(term)
    )

    if La_percent < 0:
        raise ValidationError("L_a calculado < 0%.")

    # A correlação fornece percentual; esta é a fração usada no balanço de carga.
    La_fraction = La_percent / 100.0

    return {
        "La_percent": La_percent,
        "La_fraction": La_fraction,
    }


def mixed_ehl_archard_wear_volume(
    wear_coefficient,
    normal_force,
    asperity_load_fraction,
    sliding_distance,
    hardness,
):
    """
    Archard modificado para o desgaste associado ao contato de asperezas:

        F_asp = L_a F_n
        V_w   = k_w F_asp s / H

    O coeficiente de desgaste K permanece o mesmo parâmetro de entrada
    do V7.77; somente a força efetivamente aplicada à lei de Archard
    passa a ser a parcela sustentada pelas asperezas.
    """
    if wear_coefficient < 0:
        raise ValidationError("k_w deve ser >= 0.")
    if normal_force < 0:
        raise ValidationError("F_n deve ser >= 0.")
    if not 0.0 <= asperity_load_fraction <= 1.0:
        raise ValidationError("A fração de carga de asperezas deve estar entre 0 e 1.")
    if sliding_distance < 0:
        raise ValidationError("s deve ser >= 0.")
    if hardness <= 0:
        raise ValidationError("H deve ser > 0.")

    asperity_force = asperity_load_fraction * normal_force
    volume = (
        wear_coefficient
        * asperity_force
        * sliding_distance
        / hardness
    )

    return {
        "asperity_force": asperity_force,
        "volume": volume,
    }


# ============================================================
# DESGASTE — V1.6
# ============================================================


def archard_wear_depth(
    wear_volume,
    contact_area
):
    """
    Profundidade média equivalente:

        d_w = V_w / A

    Não representa a distribuição espacial real do desgaste.
    """
    if wear_volume < 0:
        raise ValidationError(
            "V_w deve ser >= 0."
        )
    if contact_area <= 0:
        raise ValidationError(
            "A deve ser > 0."
        )

    return wear_volume / contact_area


# ============================================================
# INTEGRIDADE GLOBAL — V1.7
# ============================================================


# ============================================================
# INTEGRAÇÃO DOS RESULTADOS — V1.8
# ============================================================


# ============================================================
# PIPELINE AUTOMÁTICO — V1.9
# ============================================================


# ============================================================
# PIPELINE AUTOMÁTICO DOS MÓDULOS — V7.19
# ============================================================


# ============================================================
# ESTADO DE LIMITES + INTEGRIDADE AUTOMÁTICA — V7.19
# ============================================================


# ============================================================
# ESTADO FINAL NORMALIZADO — V7.19
# ============================================================


# ============================================================
# VIDA LIMITE INTEGRADA — V7.19
# ============================================================


# ============================================================
# CÁLCULO AUTOMÁTICO DAS VIDAS — V7.19
# ============================================================


# ============================================================
# CONVERSÃO TEMPO ↔ CICLOS — V7.19
# ============================================================

def time_to_cycles(time_value, cycle_frequency):
    """
    Converte tempo em ciclos somente quando a frequência de ciclo
    é explicitamente fornecida:

        N = f t
    """
    if time_value < 0:
        raise ValidationError("tempo deve ser >= 0.")
    if cycle_frequency <= 0:
        raise ValidationError(
            "frequência de ciclo deve ser > 0."
        )
    return cycle_frequency * time_value


def cycles_to_time(cycles, cycle_frequency):
    """
        t = N / f
    """
    if cycles < 0:
        raise ValidationError("ciclos devem ser >= 0.")
    if cycle_frequency <= 0:
        raise ValidationError(
            "frequência de ciclo deve ser > 0."
        )
    return cycles / cycle_frequency


# ============================================================
# ESTADO GLOBAL DO MODELO — V7.19
# ============================================================


# ============================================================
# EVOLUÇÃO DISCRETA DO ESTADO — V7.19
# ============================================================


# ============================================================
# AVALIAÇÃO FÍSICA POR ESTADO — V7.19
# ============================================================


# ============================================================
# INTEGRIDADE POR ESTADO E LIMITE DE VIDA — V7.19
# ============================================================


# ============================================================
# HISTÓRICO DE DANO E VIDA — V7.19
# ============================================================


# ============================================================
# ESTADO GLOBAL DE VIDA — V7.19
# ============================================================


# ============================================================
# RELATÓRIO DO ESTADO GLOBAL — V7.19
# ============================================================


# ============================================================
# INTEGRAÇÃO HERTZ → CAMPO 3D — V7.19
# ============================================================


# ============================================================
# HERTZ ELÍPTICO GERAL — V7.19
# ============================================================


# ============================================================
# SOLUÇÃO ELÁSTICA NUMÉRICA DO CONTATO — V7.19
# ============================================================


# ============================================================
# CONTATO ELÍPTICO EXPLÍCITO — V7.19
# ============================================================


# ============================================================
# PRESSÃO DE CONTATO DERIVADA — V7.19
# ============================================================


# ============================================================
# TENSÃO DERIVADA DO ESTADO DE CONTATO — V7.19
# ============================================================


# ============================================================
# CONTRATO DO ESTADO FÍSICO — V7.19
# ============================================================

PHYSICAL_STATE_REQUIRED = (
    "N", "t", "dt",
    "force", "Rx", "Ry", "ue", "v_slide", "Eprime", "spherical_contact",
    "eta0", "alpha_p", "Rq1", "Rq2",
    "T0", "beta", "kth", "alpha_th",
    "contact_time", "H_bulk",
    "wear_coefficient", "hardness", "contact_area",
    "wear_limit_depth",
)


# ============================================================
# VERIFICAÇÃO E RASTREABILIDADE — V7.19
# ============================================================


# ============================================================
# INTERFACE MULTIPÁGINA + ROLAGEM
# ============================================================

def bind_touch_scroll(canvas):
    """Rolagem por gesto vertical com o dedo, sem depender da scrollbar."""
    state = {"dragging": False, "last_y": 0}

    def inside_canvas(widget):
        while widget is not None:
            if widget is canvas:
                return True
            try:
                widget = widget.master
            except Exception:
                return False
        return False

    def press(event):
        if not inside_canvas(event.widget):
            return
        state["dragging"] = True
        state["last_y"] = event.y_root

    def motion(event):
        if not state["dragging"]:
            return
        dy = event.y_root - state["last_y"]
        state["last_y"] = event.y_root
        if dy:
            top, bottom = canvas.yview()
            span = max(1.0, bottom - top)
            canvas.yview_moveto(top - (dy / max(1, canvas.winfo_height())) * span)

    def release(event):
        state["dragging"] = False

    canvas.bind_all("<ButtonPress-1>", press, add="+")
    canvas.bind_all("<B1-Motion>", motion, add="+")
    canvas.bind_all("<ButtonRelease-1>", release, add="+")
    canvas.bind_all("<MouseWheel>", lambda e: canvas.yview_scroll(int(-e.delta / 120), "units"), add="+")


class App:
    BG = "#f2f2f2"

    def __init__(self, root):
        self.root = root
        root.title("Modelo Tribológico Integrado")

        # Pydroid 3 / Android: usar o tamanho real da tela e evitar uma
        # janela fixa maior que a área disponível. Não usar fullscreen, pois
        # em algumas versões do Tkinter no Android isso provoca tela preta.
        try:
            sw = int(root.winfo_screenwidth())
            sh = int(root.winfo_screenheight())
            w = max(360, min(sw, 900))
            h = max(520, min(sh, 1400))
            root.geometry(f"{w}x{h}+0+0")
            root.minsize(320, 480)
        except Exception:
            root.geometry("720x1200")
            root.minsize(320, 480)

        self.e = {}
        self.summary_vars = {}
        root.configure(bg=self.BG)

        # Página única rolável (mesmo padrão do final6): cabeçalho, entradas,
        # execução, resultados, gráficos e relatório em cartões.
        outer = tk.Frame(root, bg=self.BG)
        outer.pack(fill="both", expand=True)
        self.canvas = tk.Canvas(outer, bg=self.BG, highlightthickness=0, bd=0)
        scroll = tk.Scrollbar(outer, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(yscrollcommand=scroll.set)
        scroll.pack(side="right", fill="y")
        self.canvas.pack(side="left", fill="both", expand=True)

        self.frame = tk.Frame(self.canvas, bg=self.BG)
        self.window_id = self.canvas.create_window((0, 0), window=self.frame, anchor="nw")
        self.frame.bind("<Configure>", self._on_frame_configure)
        self.canvas.bind("<Configure>", self._on_canvas_configure)
        bind_touch_scroll(self.canvas)

        self._header()
        self.build_input_page()

    # ---------- layout ----------
    def _on_frame_configure(self, event=None):
        self.canvas.configure(scrollregion=self.canvas.bbox("all"))

    def _on_canvas_configure(self, event):
        self.canvas.itemconfigure(self.window_id, width=max(300, event.width))

    def _header(self):
        f = tk.Frame(self.frame, bg="#263238", padx=8, pady=6)
        f.pack(fill="x", padx=6, pady=(6, 5))
        tk.Button(f, text="VERBETES TÉCNICOS", command=self.show_glossary,
                  font=("Arial", 8, "bold"), bg="#455a64", fg="white",
                  relief="flat", padx=5, pady=4).pack(side="right", padx=(5, 0))

        title = tk.Frame(f, bg="#263238")
        title.pack(side="left", fill="x", expand=True)
        tk.Label(title, text="TRIBOLOGIA", fg="white", bg="#263238",
                 font=("Arial", 18, "bold"), anchor="w").pack(fill="x")
        tk.Label(title, text="Hertz  •  tensões  •  EHL  •  térmica  •  desgaste  •  vida",
                 fg="#e0e0e0", bg="#263238", anchor="w",
                 font=("Arial", 8)).pack(fill="x", pady=(1, 0))

    def _card(self, title):
        box = tk.Frame(self.frame, bg="#ffffff", bd=1, relief="solid")
        box.pack(fill="x", padx=6, pady=4)
        tk.Label(box, text=title, bg="#ffffff", anchor="w", justify="left", wraplength=330,
                 font=("Arial", 10, "bold")).pack(fill="x", padx=9, pady=(7, 4))
        return box

    def _note(self, parent, text):
        tk.Label(parent, text=text, bg="#ffffff", fg="#666666", anchor="w", justify="left",
                 wraplength=330, font=("Arial", 8)).pack(fill="x", padx=10, pady=(0, 7))

    def add_entry(self, parent, key, label, default):
        row = tk.Frame(parent, bg="#ffffff")
        row.pack(fill="x", padx=9, pady=2)
        tk.Label(row, text=label, width=28, anchor="w", bg="#ffffff",
                 font=("Arial", 9)).pack(side="left")
        ent = tk.Entry(row, width=11, justify="right", font=("Arial", 9))
        ent.pack(side="right")
        ent.insert(0, default)
        self.e[key] = ent

    # --------------------------------------------------------
    # PÁGINA ÚNICA
    # --------------------------------------------------------

    def build_input_page(self):
        # V7.77 — as entradas abaixo são a fonte do caso executado pela GUI.
        # Unidades de interface são convertidas para SI antes do cálculo.
        groups = [
            ("CONTATO E MATERIAIS",
             "R1 e R2 são os raios dos dois corpos, não os semi-eixos a e b do contato.",
             [
                 ("v777_force", "F normal [N]", "1000"),
                 ("v777_R1", "R1 [mm]", "1000000"),
                 ("v777_R2", "R2 [mm]", "72.5"),
                 ("v777_E1", "E1 [GPa]", "210"),
                 ("v777_E2", "E2 [GPa]", "210"),
                 ("v777_nu1", "ν1 [-]", "0.30"),
                 ("v777_nu2", "ν2 [-]", "0.30"),
             ]),
            ("LUBRIFICAÇÃO E CINEMÁTICA",
             "Valores estimados (V7.78): substitua pelos dados reais do seu caso.",
             [
                 ("v777_ue", "velocidade média uₑ [m/s]", "15.5"),
                 ("v777_vslide", "velocidade deslizamento [m/s]", "8"),
                 ("v777_eta0", "viscosidade η₀ [Pa.s]", "0.12"),
                 ("v777_alpha", "coef. pressão αp [Pa⁻¹]", "2e-8"),
                 ("v777_Rq1", "rugosidade Rq1 [µm]", "0.30"),
                 ("v777_Rq2", "rugosidade Rq2 [µm]", "0.30"),
             ]),
            ("TÉRMICA E ATRITO",
             "β é entrada do caso, mas não entra nas equações atuais.",
             [
                 ("v777_T0", "temperatura T₀ [°C]", "60"),
                 ("v777_mu", "coeficiente μ [-]", "0.05"),
                 ("v777_beta", "coef. térmico β [K⁻¹]", "0.03"),
                 ("v777_kcond", "condutividade [W/(m.K)]", "40"),
                 ("v777_diff", "difusividade [m²/s]", "1e-5"),
                 ("v777_tcontact", "tempo de contato [s]", "1e-3"),
                 ("v777_hcond", "condutância [W/K]", "100"),
             ]),
            ("DESGASTE", None,
             [
                 ("v777_kwear", "coeficiente de desgaste K [-]", "1e-8"),
                 ("v777_H", "dureza H [GPa]", "7"),
                 ("v777_tslide", "tempo de deslizamento [s]", "3600"),
             ]),
            ("LIMITES ADMISSÍVEIS", None,
             [
                 ("v777_siglim", "limite σ [GPa]", "1.80"),
                 ("v777_Tlim", "limite temperatura [°C]", "180"),
                 ("v777_hlim", "limite hmin [µm]", "0.50"),
                 ("v777_wlim", "limite desgaste [µm]", "100"),
             ]),
            ("CICLOS, FADIGA E PITTING",
             "N de avaliação é o número de ciclos em que o estado acumulado é mostrado. "
             "A vida N_lim é calculada pelo próprio modelo.",
             [
                 ("v777_freq", "frequência de ciclos [Hz]", "100"),
                 ("v777_N_eval", "ciclos para avaliação N [-]", "100000"),
                 ("v777_fat_sigref", "fadiga σ_ref [GPa]", "1.80"),
                 ("v777_fat_nref", "fadiga N_ref [ciclos]", "1000000"),
                 ("v777_fat_m", "expoente fadiga m [-]", "0.12"),
                 ("v779_pit_siglim", "σHlim pitting [GPa]", "1.80"),
                 ("v779_pit_shmin", "S_H mínimo [-]", "1.0"),
             ]),
        ]
        for title, note, fields in groups:
            box = self._card(title)
            for key, label, default in fields:
                self.add_entry(box, key, label, default)
            if note:
                self._note(box, note)
            else:
                tk.Frame(box, bg="#ffffff", height=4).pack()

        self._run_card()
        self._summary_card()
        self._charts_cards()
        self._report_card()

    def _run_card(self):
        box = self._card("EXECUÇÃO")
        # Barra própria em Canvas: mais visível e confiável no Tkinter/Android
        self.progress_canvas = tk.Canvas(
            box, height=34, bg="#e6e6e6", highlightthickness=1,
            highlightbackground="#bdbdbd", bd=0
        )
        self.progress_canvas.pack(fill="x", padx=10, pady=(3, 3))
        self.progress_fill = self.progress_canvas.create_rectangle(
            0, 0, 0, 34, fill="#1976d2", outline=""
        )
        self.progress_text = self.progress_canvas.create_text(
            0, 17, text="0.0%   |   pronto", fill="#222222",
            font=("Arial", 9, "bold"), anchor="center"
        )
        self.progress_canvas.bind("<Configure>", self._resize_progress)

        self.status_label = tk.Label(
            box, text="Pronto para executar", bg="#ffffff", fg="#333333",
            anchor="w", justify="left", font=("Arial", 8), wraplength=330
        )
        self.status_label.pack(fill="x", padx=10, pady=(0, 4))

        buttons = tk.Frame(box, bg="#ffffff")
        buttons.pack(fill="x", padx=8, pady=(0, 8))
        self.run_button = tk.Button(
            buttons, text="EXECUTAR", command=self.run_v777_gui,
            font=("Arial", 10, "bold"), height=1
        )
        self.run_button.pack(side="left", fill="x", expand=True, padx=(2, 3))
        tk.Button(
            buttons, text="GRAVAR TXT", command=self.save_results_txt,
            font=("Arial", 9), height=1
        ).pack(side="left", fill="x", expand=True, padx=(3, 2))

    def _resize_progress(self, event=None):
        w = max(1, self.progress_canvas.winfo_width())
        h = max(1, self.progress_canvas.winfo_height())
        self.progress_canvas.coords(self.progress_fill, 0, 0, 0, h)
        self.progress_canvas.coords(self.progress_text, w / 2, h / 2)

    def _set_progress(self, p, text):
        p = max(0.0, min(1.0, float(p)))
        w = max(1, self.progress_canvas.winfo_width())
        h = max(1, self.progress_canvas.winfo_height())
        self.progress_canvas.coords(self.progress_fill, 0, 0, w * p, h)
        self.progress_canvas.coords(self.progress_text, w / 2, h / 2)
        self.progress_canvas.itemconfigure(
            self.progress_text, text=f"{100*p:5.1f}%   |   {text}"
        )

    SUMMARY_DEFS = [
        ("sigma", "Tensão de von Mises máx."),
        ("p", "Pressão p0 / pm"),
        ("a", "Semi-eixo a = b"),
        ("h", "Filme h_min / h_c"),
        ("lam", "Λ = h_min / Rq"),
        ("T", "Temperatura total"),
        ("ug", "U_global"),
        ("life", "Vida limite"),
        ("cpu", "Tempo computacional"),
    ]

    def _summary_card(self):
        box = self._card("RESULTADOS PRINCIPAIS")
        for key, label in self.SUMMARY_DEFS:
            row = tk.Frame(box, bg="#ffffff")
            row.pack(fill="x", padx=9, pady=1)
            tk.Label(row, text=label, bg="#ffffff", fg="#555555",
                     anchor="w", font=("Arial", 8)).pack(side="left", fill="x", expand=True)
            v = tk.StringVar(value="—")
            self.summary_vars[key] = v
            tk.Label(row, textvariable=v, bg="#ffffff", fg="#111111",
                     anchor="e", justify="right", font=("Arial", 8, "bold"),
                     wraplength=190).pack(side="right", padx=(4, 0))
        tk.Frame(box, bg="#ffffff", height=6).pack()

    def _clear_summary(self):
        for v in self.summary_vars.values():
            v.set("—")

    def _fill_summary(self, result):
        h = result["hertz"]
        g = result["integrity"]
        lam = result["ehl"]["Lambda"]
        life = result["life"]
        self.summary_vars["sigma"].set(
            "%.1f MPa  (U_σ = %.3f)" % (result["stress"]["sigma_VM_max"] / 1e6, g["U_sigma"]))
        self.summary_vars["p"].set("%.1f / %.1f MPa" % (h["p0"] / 1e6, h["pm"] / 1e6))
        self.summary_vars["a"].set("%.4f mm" % (h["a"] * 1e3))
        self.summary_vars["h"].set("%.3f / %.3f µm" % (result["ehl"]["hmin"] * 1e6, result["ehl"]["hc"] * 1e6))
        self.summary_vars["lam"].set("%.3f" % lam if math.isfinite(lam) else "∞ (Rq = 0)")
        self.summary_vars["T"].set("%.1f °C" % (result["thermal"]["T_total"] - 273.15))
        self.summary_vars["ug"].set(
            "%.3f — %s" % (g["U_global"], "limite atingido" if g["U_global"] >= 1.0 else "admissível"))
        self.summary_vars["life"].set(
            "%.3e ciclos  (%s)" % (life["N_lim"], life["limiting_mechanism"]))

    def _report_card(self):
        box = self._card("RELATÓRIO COMPLETO")
        tk.Label(box, text="GEOMETRIA DE CONTATO", bg="#ffffff", fg="#607d8b", anchor="w",
                 font=("Arial", 7, "bold")).pack(fill="x", padx=10)
        self.v777_derived = tk.Text(box, height=7, width=40, wrap="word", font=("Courier", 8),
                                    bg="#fafafa", bd=0, relief="flat", padx=6, pady=5)
        self.v777_derived.pack(fill="x", padx=7, pady=(2, 6))
        tk.Label(box, text="RELATÓRIO TÉCNICO", bg="#ffffff", fg="#607d8b", anchor="w",
                 font=("Arial", 7, "bold")).pack(fill="x", padx=10)
        frame = tk.Frame(box, bg="#fafafa")
        frame.pack(fill="both", expand=True, padx=7, pady=(2, 4))
        self.v777_report = tk.Text(frame, height=24, width=40, wrap="word", font=("Courier", 8),
                                   bg="#fafafa", bd=0, relief="flat", padx=6, pady=5,
                                   yscrollcommand=lambda first, last: report_scroll.set(first, last))
        report_scroll = tk.Scrollbar(frame, orient="vertical", command=self.v777_report.yview)
        self.v777_report.pack(side="left", fill="both", expand=True)
        report_scroll.pack(side="right", fill="y")
        self._note(box, "A tela mostra o relatório técnico. O TXT inclui também as tabelas numéricas dos 9 gráficos.")

    CHART_DEFS = [
        ("p", "1 — Pressão de contato"),
        ("sz", "2 — Tensões × profundidade (eixo)"),
        ("h", "3 — Espessura de filme × velocidade"),
        ("lam", "4 — Λ × velocidade"),
        ("T", "5 — Temperatura × velocidade de deslizamento"),
        ("U", "6 — Integridade (U)"),
        ("sn", "7 — Curva S-N"),
        ("F", "8 — Margens × força"),
        ("N", "9 — Evolução por ciclos (N / N_lim)"),
    ]

    def _charts_cards(self):
        self.charts = {}
        self._chart_data = None
        self._busy = False
        self._chart_w = 0
        self._chart_job = None
        self._txt = ""
        self._fonts = {"t": tkfont.Font(family="Arial", size=7),
                       "l": tkfont.Font(family="Arial", size=8),
                       "b": tkfont.Font(family="Arial", size=8, weight="bold")}
        # todas as folgas escalam com a altura de linha da fonte (não são pixels fixos)
        self._ls = max(14, self._fonts["l"].metrics("linespace"))
        for key, title in self.CHART_DEFS:
            box = self._card(title)
            cv = tk.Canvas(box, height=260, bg="#fafafa", highlightthickness=1,
                           highlightbackground="#bdbdbd", bd=0)
            cv.pack(fill="x", padx=8, pady=(0, 8))
            cv.bind("<Configure>", self._on_chart_configure)
            self.charts[key] = cv

    def _on_chart_configure(self, event):
        """Redesenha só se a largura mudou, agrupando eventos (evita laço e travamento)."""
        if abs(event.width - self._chart_w) < 2:
            return
        self._chart_w = event.width
        if self._chart_job is not None:
            self.root.after_cancel(self._chart_job)
        self._chart_job = self.root.after(150, self._redraw_charts)

    def _redraw_charts(self):
        self._chart_job = None
        self._draw_charts()

    @staticmethod
    def _ticks(lo, hi, n=4):
        span = hi - lo
        if span <= 0:
            return [lo]
        step = 10 ** math.floor(math.log10(span / n))
        for m in (1, 2, 5, 10):
            if span / (step * m) <= n:
                step *= m
                break
        t0 = math.ceil(lo / step - 1e-9) * step
        out, v = [], t0
        while v <= hi + 1e-9 * span:
            out.append(0.0 if abs(v) < 1e-12 * max(1.0, abs(hi)) else v)
            v += step
        return out

    @staticmethod
    def _fmt_ticks(ticks):
        """Mesmo número de casas decimais em todas as marcas do eixo."""
        d = 0
        while d < 6 and any(abs(round(t, d) - t) > 1e-9 * max(1.0, abs(t)) for t in ticks):
            d += 1
        return [f"{t:.{d}f}" for t in ticks]

    @staticmethod
    def _ds(xs, ys, maxn=220):
        n = len(xs)
        if n <= maxn:
            return list(xs), list(ys)
        k = int(math.ceil(n / maxn))
        ix = list(range(0, n, k))
        if ix[-1] != n - 1:
            ix.append(n - 1)
        return [xs[i] for i in ix], [ys[i] for i in ix]

    def _tw(self, cv, text, font):
        """Largura REAL do texto, medida pelo próprio Canvas (bbox de um item de texto temporário)."""
        i = cv.create_text(-5000, -5000, text=text, anchor="nw", font=font)
        b = cv.bbox(i)
        cv.delete(i)
        return float(b[2] - b[0]) if b else float(font.measure(text))

    def _legend_pos(self, cv, items, w):
        """Posições dos itens da legenda; larguras medidas no Canvas, com folga de segurança."""
        ls = self._ls; f = self._fonts["l"]
        pad = max(10, 0.8 * ls)
        sq, tg = max(9, 0.55 * ls), max(8, 0.55 * ls)
        gap = max(30, 2.4 * ls)                       # espaço horizontal mínimo entre itens
        pitch = max(28, 2.1 * ls)                     # passo vertical entre linhas da legenda
        x, row, pos = pad, 0, []
        for col, lab in items:
            iw = sq + tg + self._tw(cv, lab, f) * 1.08 + 2
            if x > pad and x + iw > w - pad:
                x, row = pad, row + 1
            pos.append((x, row, iw))
            x += iw + gap
        return pos, (row + 1) * pitch + 0.7 * ls, (sq, tg, pitch)

    def _draw_legend(self, cv, items, pos, y0, dims):
        f = self._fonts["l"]
        sq, tg, pitch = dims
        for (col, lab), (x, row, iw) in zip(items, pos):
            yc = y0 + 0.35 * self._ls + row * pitch + pitch / 2
            cv.create_rectangle(x, yc - sq / 2, x + sq, yc + sq / 2, fill=col, outline=col)
            cv.create_text(x + sq + tg, yc, text=lab, anchor="w", fill="#333333", font=f)

    def _plot(self, cv, series, logy=False, xrange=None, markers=False,
              xlab="", yl=None, yr=None, empty="sem dados"):
        """series: (xs, ys, cor, rótulo, eixo 'l'|'r'[, dash]). yl/yr: (texto, cor) do rótulo de unidade."""
        ft, fl, fb = self._fonts["t"], self._fonts["l"], self._fonts["b"]
        ls = self._ls
        cv.delete("all")
        w = max(240, cv.winfo_width())
        plot_h = max(150, 8.5 * ls)
        series = [s for s in series if len(s[0]) >= 1]
        if not series:
            cv.configure(height=plot_h * 0.45 + 40)
            cv.create_text(w / 2, (plot_h * 0.45 + 40) / 2, text=empty, fill="#888888", font=fl)
            return
        has_r = any(s[4] == "r" for s in series)
        xs_all = [x for s in series for x in s[0]]
        x0, x1 = xrange if xrange else (min(xs_all), max(xs_all))
        if x1 - x0 < 1e-12:
            x1 = x0 + 1.0

        def rng(axis):
            vals = [y for s in series if s[4] == axis for y in s[1] if (y > 0 or not (logy and axis == "l"))]
            if not vals:
                return None
            if logy and axis == "l":
                vals = [math.log10(v) for v in vals]
            lo, hi = min(vals), max(vals)
            minspan = max(0.015 * max(abs(lo), abs(hi)), 1e-9)      # evita rótulos com casas decimais demais
            if hi - lo < minspan:
                mid = 0.5 * (hi + lo)
                lo, hi = mid - 0.5 * minspan, mid + 0.5 * minspan
            if hi - lo < 1e-9:
                hi, lo = hi + 0.5, lo - 0.5
            pad = 0.06 * (hi - lo)
            return lo - pad, hi + pad

        rl, rr = rng("l"), rng("r")
        tl = self._ticks(*rl, 5) if rl else []
        trr = self._ticks(*rr, 5) if rr else []
        lab_l = ([f"1e{t:.0f}" for t in tl] if logy else self._fmt_ticks(tl)) if tl else []
        lab_r = self._fmt_ticks(trr) if trr else []
        wl = max([self._tw(cv, s, ft) for s in lab_l] + [0])
        wr = max([self._tw(cv, s, ft) for s in lab_r] + [0])
        L = max(0.9 * ls + 22, wl + 0.8 * ls)
        R = (wr + 0.8 * ls) if has_r else 0.8 * ls
        items = [(s[2], s[3]) for s in series]
        lpos, lh, ldims = self._legend_pos(cv, items, w)
        unit_h = 1.7 * ls if (yl or yr) else 0
        T_ = 0.4 * ls + lh + unit_h + 0.5 * ls
        B = 0.5 * ls + 1.2 * ls + 0.5 * ls + (1.4 * ls if xlab else 0) + 0.4 * ls
        H = T_ + plot_h + B
        cv.configure(height=int(H))
        pw, ph_ = w - L - R, plot_h
        X = lambda v: L + (v - x0) / (x1 - x0) * pw

        def Y(v, r, lg):
            if lg:
                v = math.log10(max(v, 1e-30))
            return T_ + ph_ - (v - r[0]) / (r[1] - r[0]) * ph_

        # legenda e rótulos de unidade: faixas separadas, sem tocar nas curvas
        self._draw_legend(cv, items, lpos, 0.4 * ls, ldims)
        ytxt = 0.4 * ls + lh + 0.85 * ls
        if yl:
            cv.create_text(8, ytxt, text=yl[0], anchor="w", fill=yl[1], font=fb)
        if yr:
            cv.create_text(w - 8, ytxt, text=yr[0], anchor="e", fill=yr[1], font=fb)
        cv.create_rectangle(L, T_, w - R, T_ + ph_, outline="#bdbdbd")
        # eixo x
        nx = max(2, int(pw // (4.2 * ls)))
        tx_ticks = self._ticks(x0, x1, nx)
        for tx, lab in zip(tx_ticks, self._fmt_ticks(tx_ticks)):
            cv.create_line(X(tx), T_, X(tx), T_ + ph_, fill="#ececec")
            cv.create_line(X(tx), T_ + ph_, X(tx), T_ + ph_ + 4, fill="#757575")
            cv.create_text(X(tx), T_ + ph_ + 0.5 * ls, text=lab, anchor="n", fill="#444444", font=ft)
        if xlab:
            cv.create_text(L + pw / 2, H - 0.4 * ls, text=xlab, anchor="s", fill="#666666", font=ft)
        # eixos y
        lcols = [s[2] for s in series if s[4] == "l"]
        lcol = lcols[0] if len(lcols) == 1 else "#333333"
        for tv, lab in zip(tl, lab_l):
            yy = Y(10 ** tv if logy else tv, rl, logy)
            cv.create_line(L, yy, w - R, yy, fill="#ececec")
            cv.create_line(L - 4, yy, L, yy, fill="#757575")
            cv.create_text(L - 7, yy, text=lab, anchor="e", fill=lcol, font=ft)
        rcols = [s[2] for s in series if s[4] == "r"]
        rcol = rcols[0] if len(rcols) == 1 else "#333333"
        for tv, lab in zip(trr, lab_r):
            yy = Y(tv, rr, False)
            cv.create_line(w - R, yy, w - R + 4, yy, fill="#757575")
            cv.create_text(w - R + 7, yy, text=lab, anchor="w", fill=rcol, font=ft)
        # curvas
        for s in series:
            xs, ys, col, lab, ax = s[:5]
            dash = s[5] if len(s) > 5 else None
            r = rl if ax == "l" else rr
            lg = logy and ax == "l"
            if r is None:
                continue
            xd, yd = self._ds(xs, ys)
            pts = [q for a, b in zip(xd, yd) if (b > 0 or not lg) for q in (X(a), Y(b, r, lg))]
            if len(pts) >= 4:
                cv.create_line(*pts, fill=col, width=2, dash=dash)
            if markers or len(pts) == 2:
                for i in range(0, len(pts), 2):
                    cv.create_oval(pts[i] - 3.5, pts[i + 1] - 3.5, pts[i] + 3.5, pts[i + 1] + 3.5, fill=col, outline=col)

    def _draw_bars(self, cv, items, empty="sem dados"):
        """Barras horizontais de U (limite em 1). items: (nome, valor)."""
        cv.delete("all")
        ls = self._ls
        w = max(240, cv.winfo_width())
        ft, fl = self._fonts["t"], self._fonts["l"]
        vals = [(k, v) for k, v in items if v is not None and math.isfinite(v)]
        if not vals:
            cv.configure(height=int(max(150, 8.5 * ls) * 0.45 + 40))
            cv.create_text(w / 2, (max(150, 8.5 * ls) * 0.45 + 40) / 2, text=empty, fill="#888888", font=fl)
            return
        rowh = max(30, 2.2 * ls)
        bh = max(10, 0.7 * ls)
        top = 1.6 * ls
        H = top + len(vals) * rowh + 0.6 * ls
        cv.configure(height=int(H))
        x0 = 12 + max(self._tw(cv, k, ft) for k, _ in vals) * 1.08 + 12
        vw = self._tw(cv, "0.000", ft) * 1.08
        bw = max(40, w - x0 - vw - 24)
        vmax = max(1.2, max(v for _, v in vals) * 1.05)
        X = lambda v: x0 + bw * v / vmax
        cv.create_text(X(1.0), top - 0.6 * ls, text="limite = 1", anchor="s", fill="#c62828", font=ft)
        for i, (k, v) in enumerate(vals):
            y = top + i * rowh
            col = "#c62828" if v >= 1.0 else "#1976d2"
            cv.create_text(x0 - 10, y + rowh / 2, text=k, anchor="e", fill="#333333", font=ft)
            cv.create_rectangle(x0, y + (rowh - bh) / 2, max(x0 + 1, X(v)), y + (rowh + bh) / 2, fill=col, outline=col)
            cv.create_text(max(x0 + 1, X(v)) + 6, y + rowh / 2, text=f"{v:.3f}", anchor="w", fill="#555555", font=ft)
        cv.create_line(X(1.0), top - 0.3 * ls, X(1.0), top + len(vals) * rowh, fill="#c62828", dash=(4, 3))

    def _draw_charts(self, only=None):
        if not hasattr(self, "charts"):
            return
        d = self._chart_data
        wait = "calculando..." if self._busy else "sem dados"
        B_, O_, G_, R_, P_, K_, N_ = "#1976d2", "#ef6c00", "#2e7d32", "#c62828", "#6a1b9a", "#212121", "#8d6e63"
        DASH = (4, 3)
        ends = lambda xs, v: ([xs[0], xs[-1]], [v, v])

        def p_job(c, q):
            self._plot(c, [(q["x"], q["p"], B_, "p(r)", "l"),
                           (*ends(q["x"], q["pm"]), O_, "p_m", "l", DASH)],
                       xlab="r (mm)", yl=("MPa", B_))

        def sz_job(c, q):
            self._plot(c, [(q["z"], q["sz"], B_, "σz", "l"), (q["z"], q["sr"], G_, "σr", "l"),
                           (q["z"], q["vm"], R_, "σVM", "l"), (q["z"], q["tm"], P_, "τmax", "l"),
                           (*ends(q["z"], q["lim"]), K_, "σ lim", "l", DASH)],
                       xlab="profundidade z (mm)", yl=("MPa", "#333333"))

        def h_job(c, q):
            ser = [(q["u"], q["hmin"], B_, "h_min", "l"), (q["u"], q["hc"], O_, "h_c", "l")]
            if q["hmin_lim"] > 0:
                ser.append((*ends(q["u"], q["hmin_lim"]), R_, "h_min lim", "l", DASH))
            self._plot(c, ser, xlab="velocidade de arrastamento ue (m/s)", yl=("µm", "#333333"))

        def lam_job(c, q):
            if not q["lam"]:
                return self._plot(c, [], empty="Rq = 0: Λ indefinido")
            ser = [(q["u"], q["lam"], G_, "Λ", "l"), (*ends(q["u"], 1.0), R_, "Λ = 1", "l", DASH)]
            if q["lam_min"] > 0:
                ser.append((*ends(q["u"], q["lam_min"]), O_, "Λ mín", "l", DASH))
            self._plot(c, ser, xlab="velocidade de arrastamento ue (m/s)", yl=("Λ", G_))

        def t_job(c, q):
            self._plot(c, [(q["v"], q["tb"], B_, "T_bulk", "l"), (q["v"], q["tf"], O_, "T_flash", "l"),
                           (q["v"], q["tt"], R_, "T_total", "l"),
                           (*ends(q["v"], q["lim"]), K_, "T lim", "l", DASH)],
                       xlab="velocidade de deslizamento (m/s)", yl=("°C", "#333333"))

        def sn_job(c, q):
            ser = [(q["s"], q["n"], B_, "S-N", "l"), ([q["ref"][0]], [q["ref"][1]], G_, "ref", "l")]
            if q["op"]:
                ser.append(([q["op"][0]], [q["op"][1]], R_, "operação", "l"))
            self._plot(c, ser, logy=True, xlab="σ (MPa)", yl=("N (ciclos)", B_))

        def n_job(c, q):
            self._plot(c, [(q["x"], q["d"], B_, "D fadiga", "l"),
                           (q["x"], q["uw"], O_, "U_desgaste", "l"),
                           (q["x"], q["ug"], R_, "U_global", "l"),
                           (*ends(q["x"], 1.0), K_, "U = 1", "l", DASH)],
                       xlab="N / N_lim   (N_lim = %.3e ciclos)" % q["nlim"],
                       yl=("D, U", "#333333"))

        def f_job(c, q):
            ser = [(q["f"], q["usig"], B_, "U_σ", "l"), (q["f"], q["uglo"], R_, "U_global", "l"),
                   (*ends(q["f"], 1.0), K_, "U = 1", "l", DASH)]
            if q["sh"]:
                ser.append((q["f"], q["sh"], G_, "S_H", "r"))
            self._plot(c, ser, xlab="força F (N)", yl=("U", "#333333"),
                       yr=("S_H", G_) if q["sh"] else None)

        jobs = {"p": p_job, "sz": sz_job, "h": h_job, "lam": lam_job,
                "T": t_job, "sn": sn_job, "F": f_job, "N": n_job}
        for key, cv in self.charts.items():
            if only is not None and key not in only:
                continue
            try:
                src = "h" if key == "lam" else key
                q = d.get(src) if d else None
                if q is None:
                    if key == "U":
                        self._draw_bars(cv, [], empty=wait)
                    else:
                        self._plot(cv, [], empty=wait)
                elif "erro" in q:
                    self._plot(cv, [], empty="erro: " + q["erro"])
                elif key == "U":
                    self._draw_bars(cv, q["items"])
                else:
                    jobs[key](cv, q)
            except Exception as e:
                cv.delete("all")
                cv.create_text(6, 6, text=f"erro no gráfico: {e}", anchor="nw", fill="#c62828", font=self._fonts["l"])

    def _read_case(self):
        """Lê todas as entradas V7.78 da GUI e converte para SI."""
        def f(key):
            return float(self.e[key].get().strip().replace(",", "."))

        case = dict(V776_BASE_CASE)
        case.update({
            "force": f("v777_force"),
            "R1": f("v777_R1") * 1e-3,
            "R2": f("v777_R2") * 1e-3,
            "E1": f("v777_E1") * 1e9,
            "E2": f("v777_E2") * 1e9,
            "nu1": f("v777_nu1"),
            "nu2": f("v777_nu2"),
            "ue": f("v777_ue"),
            "v_slide": f("v777_vslide"),
            "eta0": f("v777_eta0"),
            "alpha_p": f("v777_alpha"),
            "Rq1": f("v777_Rq1") * 1e-6,
            "Rq2": f("v777_Rq2") * 1e-6,
            "T0": f("v777_T0") + 273.15,
            "mu": f("v777_mu"),
            "beta": f("v777_beta"),
            "thermal_conductivity": f("v777_kcond"),
            "thermal_diffusivity": f("v777_diff"),
            "contact_time": f("v777_tcontact"),
            "thermal_conductance": f("v777_hcond"),
            "wear_coefficient": f("v777_kwear"),
            "hardness": f("v777_H") * 1e9,
            "sliding_time": f("v777_tslide"),
            "sigma_limit": f("v777_siglim") * 1e9,
            "temperature_limit": f("v777_Tlim") + 273.15,
            "hmin_limit": f("v777_hlim") * 1e-6,
            "wear_limit_depth": f("v777_wlim") * 1e-6,
            "cycle_frequency": f("v777_freq"),
            "cycles_evaluation": f("v777_N_eval"),
            "fatigue_sigma_ref": f("v777_fat_sigref") * 1e9,
            "fatigue_N_ref": f("v777_fat_nref"),
            "fatigue_exponent_m": f("v777_fat_m"),
            "pitting_sigma_Hlim": f("v779_pit_siglim") * 1e9,
            "pitting_SH_min": f("v779_pit_shmin"),
        })
        return case

    def run_v777_gui(self):
        """Lê as entradas e executa o caso em segundo plano (a interface continua ativa)."""
        if self._busy:
            return
        try:
            case = self._read_case()
        except Exception as exc:
            self._v782_finish(None, None, exc)
            return
        self._busy = True
        self._chart_data = {}
        self._txt = ""
        self._sections_done = 0
        self.run_button.configure(state="disabled")
        self.v777_derived.delete("1.0", "end")
        self.v777_report.delete("1.0", "end")
        self._clear_summary()
        self._draw_charts()
        self._set_progress(0.0, "calculando o caso")
        self.status_label.configure(text="Calculando...")
        threading.Thread(target=self._v782_worker, args=(case,), daemon=True).start()

    def _v782_worker(self, case):
        try:
            t0 = time.perf_counter()
            result = execute_v776_reference_case(
                case,
                lambda done, total: self.root.after(0, self._v782_progress, done, total)
            )
            t_case = time.perf_counter() - t0
            self.root.after(0, self._v782_result, result, t_case)
            data = build_v782_chart_data(
                result,
                lambda key, q: self.root.after(0, self._v782_partial, key, q)
            )
            self.root.after(0, self._v782_finish, result, data, None,
                            t_case, time.perf_counter() - t0 - t_case)
        except Exception as exc:
            self.root.after(0, self._v782_finish, None, None, exc)

    def _v782_progress(self, done, total):
        if self._busy:
            self._set_progress(0.28 * done / max(1, total),
                               "campo de tensões %d/%d" % (done, total))

    def _v782_result(self, result, t_case):
        eq = result["equivalent"]
        self.v777_derived.delete("1.0", "end")
        self.v777_derived.insert("1.0",
            "E'  = %.8e Pa  (%.8f GPa)\n"
            "R'  = %.8e m   (%.8f mm)\n"
            "a=b = %.8e m   (%.8f mm)\n"
            "A   = %.8e m²\n"
            "pH  = %.8e Pa\n"
            "pm  = %.8e Pa" % (
                eq["Eprime"], eq["Eprime"] / 1e9,
                eq["Rprime"], eq["Rprime"] * 1e3,
                eq["a"], eq["a"] * 1e3,
                eq["contact_area"], eq["p0"], eq["pm"]
            )
        )
        self.v777_report.delete("1.0", "end")
        self.v777_report.insert("1.0", format_v776_reference_report(result))
        self._fill_summary(result)
        self.summary_vars["cpu"].set("caso %.2f s" % t_case)
        self._set_progress(0.3, "calculando gráficos")
        self.status_label.configure(text="Caso calculado; gerando gráficos...")

    def _v782_partial(self, key, q):
        self._sections_done += 1
        total = len(self.CHART_DEFS) - 1
        self._set_progress(0.3 + 0.7 * self._sections_done / total,
                           "bloco de gráficos %d/%d" % (self._sections_done, total))
        self._chart_data[key] = q
        self._draw_charts(only={key, "lam"} if key == "h" else {key})

    def _v782_finish(self, result, data, exc, t_case=0.0, t_charts=0.0):
        self._busy = False
        self.run_button.configure(state="normal")
        if exc is not None:
            self.v777_derived.delete("1.0", "end")
            self.v777_report.delete("1.0", "end")
            self.v777_report.insert("1.0", "ERRO V7.78\n\n" + str(exc))
            self._chart_data = None
            self._txt = ""
            self._clear_summary()
            self._draw_charts()
            self._set_progress(0.0, "erro")
            self.status_label.configure(text="Erro de entrada/cálculo (veja o relatório).")
            return
        self._chart_data = data
        timing = "caso %.2f s + gráficos %.2f s" % (t_case, t_charts)
        self.summary_vars["cpu"].set(timing)
        self._txt = format_v782_results_txt(result, data, timing)
        self._set_progress(1.0, "concluído")
        self.status_label.configure(text="Caso físico calculado.")

    def save_results_txt(self):
        text = getattr(self, "_txt", "")
        if not text:
            self.status_label.configure(text="Nenhum resultado disponível para gravação.")
            return
        name = "Tribologico_resultado_V781.txt"
        paths = [
            "/storage/emulated/0/Download/" + name,
            "/sdcard/Download/" + name,
            os.path.join(os.path.expanduser("~"), "Download", name),
        ]
        saved = None
        for path in paths:
            try:
                os.makedirs(os.path.dirname(path), exist_ok=True)
                with open(path, "w", encoding="utf-8") as fh:
                    fh.write(text)
                saved = path
                break
            except Exception:
                pass
        if saved:
            msg = "TXT gravado: " + saved
        else:
            msg = "TXT não pôde ser gravado no armazenamento."
        self.status_label.configure(text=msg)

    # ---------- verbetes técnicos (padrão do final6) ----------
    def show_glossary(self):
        # Página própria dentro da janela principal; não usa janela suspensa.
        if hasattr(self, "glossary_page") and self.glossary_page.winfo_exists():
            self.glossary_page.lift()
            return

        self.glossary_page = tk.Frame(self.root, bg="#f2f2f2")
        self.glossary_page.place(x=0, y=0, relwidth=1, relheight=1)

        outer = tk.Frame(self.glossary_page, bg="#f2f2f2")
        outer.pack(fill="both", expand=True)
        canvas = tk.Canvas(outer, bg="#f2f2f2", highlightthickness=0)
        canvas.pack(fill="both", expand=True)
        frame = tk.Frame(canvas, bg="#f2f2f2")
        wid = canvas.create_window((0, 0), window=frame, anchor="nw")
        frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))

        top = tk.Frame(frame, bg="#263238")
        top.pack(fill="x", padx=6, pady=(6, 4))
        tk.Button(top, text="← VOLTAR", command=self.close_glossary,
                  font=("Arial", 9, "bold"), bg="#455a64", fg="white",
                  relief="flat", padx=8, pady=5).pack(side="left", padx=6, pady=6)
        tk.Label(top, text="VERBETES TÉCNICOS", bg="#263238", fg="white",
                 font=("Arial", 14, "bold"), anchor="w"
                 ).pack(side="left", padx=(4, 6), pady=6)
        tk.Label(frame, text="Significado das grandezas e regras efetivamente presentes no código.",
                 bg="#f2f2f2", fg="#444444", justify="left", anchor="w",
                 font=("Arial", 8), wraplength=330).pack(fill="x", padx=10, pady=(0, 5))

        terms = GLOSSARY_TERMS

        # Filtro por seção para navegação dos verbetes.
        control = tk.Frame(frame, bg="#ffffff", bd=1, relief="solid")
        control.pack(fill="x", padx=6, pady=(0, 5))
        tk.Label(control, text="SEÇÃO", bg="#ffffff", fg="#555555",
                 font=("Arial", 8, "bold")).pack(anchor="w", padx=8, pady=(6, 1))
        section_names = ["TODAS AS SEÇÕES"] + sorted(dict.fromkeys(x[0] for x in terms), key=lambda x: int(x.split(" — ", 1)[0]))
        section_var = tk.StringVar(value=section_names[0])
        section_menu = tk.OptionMenu(control, section_var, *section_names)
        section_menu.configure(font=("Arial", 8), anchor="w", width=38)
        section_menu.pack(fill="x", padx=8, pady=(0, 4))
        count_var = tk.StringVar(value=f"{len(terms)} verbetes")
        tk.Label(control, textvariable=count_var, bg="#ffffff", fg="#546e7a",
                 font=("Arial", 8), anchor="w").pack(fill="x", padx=8, pady=(0, 6))

        results_frame = tk.Frame(frame, bg="#f2f2f2")
        results_frame.pack(fill="x", expand=True)

        def render_terms(*_):
            for child in results_frame.winfo_children():
                child.destroy()
            selected = section_var.get()
            shown = 0
            for sec, name, meaning in terms:
                if selected != "TODAS AS SEÇÕES" and sec != selected:
                    continue
                card = tk.Frame(results_frame, bg="#ffffff", bd=1, relief="solid")
                card.pack(fill="x", padx=6, pady=2)
                tk.Label(card, text=name, bg="#ffffff", fg="#111111", anchor="w",
                         font=("Arial", 9, "bold"), wraplength=320).pack(fill="x", padx=8, pady=(5, 1))
                tk.Label(card, text=sec, bg="#ffffff", fg="#607d8b", anchor="w",
                         font=("Arial", 7, "bold"), wraplength=320).pack(fill="x", padx=8, pady=(0, 1))
                tk.Label(card, text=meaning, bg="#ffffff", fg="#444444", anchor="w",
                         justify="left", wraplength=320, font=("Arial", 8)).pack(fill="x", padx=8, pady=(0, 5))
                shown += 1
            count_var.set(f"{shown} de {len(terms)} verbetes")
            canvas.configure(scrollregion=canvas.bbox("all"))
            canvas.yview_moveto(0)

        section_var.trace_add("write", render_terms)
        render_terms()

        # Texto responsivo: ocupa toda a largura útil em retrato e paisagem.
        def resize_glossary(event=None):
            w = max(260, canvas.winfo_width() - 32)
            for card in results_frame.winfo_children():
                card.configure(width=w)
                for sub in card.winfo_children():
                    sub.configure(wraplength=max(220, w - 24))
        canvas.bind("<Configure>", lambda e: (canvas.itemconfigure(wid, width=max(300, e.width)), resize_glossary(e)))
        bind_touch_scroll(canvas)
        canvas.after(50, resize_glossary)

    def close_glossary(self):
        if hasattr(self, "glossary_page") and self.glossary_page.winfo_exists():
            self.glossary_page.destroy()


    # --------------------------------------------------------
    # PÁGINA GEOMETRIA
    # --------------------------------------------------------


    # --------------------------------------------------------
    # PÁGINA CONTATO
    # --------------------------------------------------------


    # --------------------------------------------------------
    # PÁGINA PRESSÃO
    # --------------------------------------------------------


    # --------------------------------------------------------
    # PÁGINA TENSÃO
    # --------------------------------------------------------


    # --------------------------------------------------------
    # PÁGINA TENSOR
    # --------------------------------------------------------


    # --------------------------------------------------------
    # PÁGINA CAMPO ESPACIAL
    # --------------------------------------------------------


    # --------------------------------------------------------
    # PÁGINA INTEGRIDADE
    # --------------------------------------------------------


    # --------------------------------------------------------
    # PÁGINA FADIGA
    # --------------------------------------------------------


    # --------------------------------------------------------
    # PÁGINA VIDA
    # --------------------------------------------------------


    # --------------------------------------------------------
    # PÁGINA EHL
    # --------------------------------------------------------


    # --------------------------------------------------------
    # PÁGINA TÉRMICA
    # --------------------------------------------------------


    # --------------------------------------------------------
    # PÁGINA ACOPLAMENTO TÉRMICO-EHL
    # --------------------------------------------------------


    # --------------------------------------------------------
    # PÁGINA DESGASTE
    # --------------------------------------------------------


    # --------------------------------------------------------
    # PÁGINA INTEGRIDADE GLOBAL
    # --------------------------------------------------------


    # --------------------------------------------------------
    # PÁGINA INTEGRAÇÃO
    # --------------------------------------------------------


    # --------------------------------------------------------
    # PÁGINA PIPELINE
    # --------------------------------------------------------


    # --------------------------------------------------------
    # PÁGINA PIPELINE AUTOMÁTICO
    # --------------------------------------------------------


    # --------------------------------------------------------
    # PÁGINA VOLUME
    # --------------------------------------------------------


    # --------------------------------------------------------
    # EXECUÇÃO
    # --------------------------------------------------------


    # --------------------------------------------------------
    # SAÍDAS
    # --------------------------------------------------------


# ============================================================
# V7.62 — CONSTRUCAO DO CAMPO COMUM A PARTIR DE E2 + E4 + EHL
# ============================================================


if False and __name__ == "__main__":
    root = tk.Tk()
    App(root)
    root.mainloop()


# ========================= V7.66 — FADIGA MECÂNICA EFETIVA =========================


# ========================= V7.67 — INTEGRIDADE EVOLUTIVA UNIFICADA =========================


# ========================= V7.68 — EVOLUÇÃO COMPLETA TEMPO/CICLOS =========================


# ========================= V7.69 — FECHAMENTO DO PIPELINE E1→E5 + VIDA =========================


# ============================================================
# V7.70 — AUDITORIA FÍSICO-MATEMÁTICA DO CASO FECHADO
# ============================================================


# ============================================================
# V7.77 — FECHAMENTO DAS LACUNAS FÍSICO-MATEMÁTICAS
# ============================================================


# ============================================================
# V7.77 — CASO FÍSICO DE REFERÊNCIA
# ============================================================
#
# Esta extensão não substitui os módulos anteriores.
# Ela estabelece um caso físico central e deriva os resultados
# a partir de entradas físicas primárias.
#
# Entradas centrais:
#     F  = 1000 N
#     R1 = 25 mm
#     R2 = 25 mm
#     E1 = E2 = 210 GPa
#     nu1 = nu2 = 0.30
#
# O raio relativo é DERIVADO:
#     1/R' = 1/R1 + 1/R2
#
# Entradas tribológicas não determinadas por Hertz permanecem
# explicitamente classificadas como ESTIMATIVAS.
# ============================================================

V777_VERSION = "V7.78"
V779_VERSION = "V7.79"
V780_VERSION = "V7.81"

V776_BASE_CASE = {
    # -------------------------
    # Entradas físicas primárias
    # -------------------------
    "force": 1000.0,                 # N
    "R1": 1.0e3,                      # m (1e6 mm = infinito numerico)
    "R2": 72.5e-3,                    # m (72.5 mm)
    "E1": 210.0e9,                   # Pa
    "E2": 210.0e9,                   # Pa
    "nu1": 0.30,
    "nu2": 0.30,

    # -------------------------
    # Entradas tribológicas estimadas
    # -------------------------
    "ue": 15.5,                      # m/s
    "v_slide": 8.0,                   # m/s
    "eta0": 0.12,                    # Pa.s
    "alpha_p": 2.0e-8,               # Pa^-1
    "Rq1": 0.30e-6,                  # m
    "Rq2": 0.30e-6,                  # m
    "T0": 60.0 + 273.15,             # K
    "mu": 0.05,

    # -------------------------
    # Parâmetros térmicos estimados
    # -------------------------
    "beta": 0.03,                    # K^-1
    "thermal_conductivity": 40.0,    # W/(m.K)
    "thermal_diffusivity": 1.0e-5,   # m²/s
    "contact_time": 1.0e-3,          # s
    "thermal_conductance": 100.0,    # W/K

    # -------------------------
    # Desgaste estimado
    # -------------------------
    "wear_coefficient": 1.0e-8,
    "hardness": 7.0e9,               # Pa
    "sliding_time": 3600.0,          # s
    "cycle_frequency": 100.0,        # Hz
    "cycles_evaluation": 100000.0,   # ciclos

    # -------------------------
    # Limites de referência estimados
    # -------------------------
    "sigma_limit": 1.80e9,           # Pa
    "fatigue_sigma_ref": 1.80e9,     # Pa
    "fatigue_N_ref": 1.0e6,          # ciclos
    "fatigue_exponent_m": 0.12,

    # -------------------------
    # Pitting — V7.79
    # Criterio convencional de resistência superficial de worm:
    # S_H = sigma_HG / sigma_Hm ; S_H,min = 1.0.
    # sigma_HG deve ser fornecido para o material/condição reais.
    # -------------------------
    "pitting_sigma_Hlim": 1.80e9,   # Pa; referência de contato já usada no caso-base
    "pitting_SH_min": 1.0,           # fator mínimo de segurança
    "temperature_limit": 180.0 + 273.15,
    "hmin_limit": 0.50e-6,
    "wear_limit_depth": 100.0e-6,    # m
}


def derive_v776_equivalent_properties(case):
    """Deriva E' e R' a partir dos dois corpos do caso-base."""
    E1 = float(case["E1"])
    E2 = float(case["E2"])
    nu1 = float(case["nu1"])
    nu2 = float(case["nu2"])
    R1 = float(case["R1"])
    R2 = float(case["R2"])

    if E1 <= 0 or E2 <= 0:
        raise ValidationError("E1 e E2 devem ser > 0.")
    if not (0.0 <= nu1 < 0.5 and 0.0 <= nu2 < 0.5):
        raise ValidationError("nu1 e nu2 devem estar em 0 <= nu < 0.5.")
    if R1 <= 0 or R2 <= 0:
        raise ValidationError("R1 e R2 devem ser > 0.")

    Eprime = 1.0 / (
        (1.0 - nu1**2) / E1
        + (1.0 - nu2**2) / E2
    )
    Rprime = 1.0 / (1.0 / R1 + 1.0 / R2)

    return {
        "Eprime": Eprime,
        "Rprime": Rprime,
    }


def derive_v776_hertz(case):
    """Calcula contato Hertziano circular a partir de F e R'."""
    eq = derive_v776_equivalent_properties(case)
    F = float(case["force"])
    Eprime = eq["Eprime"]
    Rprime = eq["Rprime"]

    if F < 0:
        raise ValidationError("F deve ser >= 0.")

    a = ((3.0 * F * Rprime) / (4.0 * Eprime)) ** (1.0 / 3.0)
    b = a
    area = math.pi * a * b
    p0 = 3.0 * F / (2.0 * area)
    pm = F / area

    return {
        **eq,
        "a": a,
        "b": b,
        "contact_area": area,
        "p0": p0,
        "pm": pm,
    }


def derive_v781_pitting(c, hertz):
    """Pitting V7.79: S_H = sigma_HG / sigma_Hm, sem inventar sigma_HG."""
    pitting_sigma_Hm = float(hertz["pm"])
    pitting_sigma_Hlim = c.get("pitting_sigma_Hlim")
    pitting_SH_min = float(c.get("pitting_SH_min", 1.0))
    if pitting_sigma_Hlim is None:
        pitting = {
            "enabled": False,
            "criterion": "S_H = sigma_HG / sigma_Hm",
            "sigma_Hm": pitting_sigma_Hm,
            "sigma_Hlim": None,
            "S_H": None,
            "S_H_min": pitting_SH_min,
            "U_pitting": None,
            "admissible": None,
            "reason": "sigma_Hlim de pitting nao informado para o material real",
        }
    else:
        pitting_sigma_Hlim = float(pitting_sigma_Hlim)
        if pitting_sigma_Hlim <= 0.0:
            raise ValidationError("pitting_sigma_Hlim deve ser > 0 quando informado.")
        if pitting_SH_min <= 0.0:
            raise ValidationError("pitting_SH_min deve ser > 0.")
        pitting_SH = pitting_sigma_Hlim / pitting_sigma_Hm if pitting_sigma_Hm > 0.0 else float("inf")
        pitting_U = pitting_SH_min / pitting_SH if pitting_SH > 0.0 else float("inf")
        pitting = {
            "enabled": True,
            "criterion": "S_H = sigma_HG / sigma_Hm",
            "sigma_Hm": pitting_sigma_Hm,
            "sigma_Hlim": pitting_sigma_Hlim,
            "S_H": pitting_SH,
            "S_H_min": pitting_SH_min,
            "U_pitting": pitting_U,
            "admissible": pitting_SH >= pitting_SH_min,
            "contact_pressure_basis": "pm (Hertz mean contact stress)",
            "load_factor_to_limit": (pitting_SH / pitting_SH_min) ** 3,
        }
    return pitting


def derive_v781_surface_capacity(c, pitting, ehl, Lambda, Rq):
    """Capacidade de superfície V7.81: margens macro e micropitting."""
    macro_capacity = {
        "enabled": bool(pitting.get("enabled", False)),
        "sigma_Hm": pitting.get("sigma_Hm"),
        "sigma_Hlim": pitting.get("sigma_Hlim"),
        "S_H": pitting.get("S_H"),
        "S_H_min": pitting.get("S_H_min"),
        "load_factor_to_limit": None,
        "force_limit": None,
    }
    if macro_capacity["enabled"] and macro_capacity["S_H"] is not None:
        SH = float(macro_capacity["S_H"])
        SHmin = float(macro_capacity["S_H_min"])
        if SH > 0.0 and SHmin > 0.0:
            # Exato para o Hertz circular atualmente implementado:
            # p_m proporcional a F^(1/3).
            macro_capacity["load_factor_to_limit"] = (SH / SHmin) ** 3
            macro_capacity["force_limit"] = float(c["force"]) * macro_capacity["load_factor_to_limit"]

    lambda_current = float(Lambda)
    rq_current = float(Rq)
    hmin_limit = float(c["hmin_limit"])
    lambda_min = hmin_limit / rq_current if rq_current > 0.0 else float("inf")

    micro_capacity = {
        "criterion": "Lambda = h_min / sqrt(Rq1^2 + Rq2^2)",
        "h_min": float(ehl["hmin"]),
        "h_min_limit": hmin_limit,
        "Rq": rq_current,
        "Lambda": lambda_current,
        "Lambda_min": lambda_min,
        "S_lambda": (lambda_current / lambda_min) if lambda_min > 0.0 and math.isfinite(lambda_min) else float("inf"),
        "admissible": lambda_current >= lambda_min,
        "load_factor_to_limit": None,
        "force_limit": None,
        "life_predictable": False,
        "time_to_first_micropitting": None,
        "reason_no_time_prediction": "Lambda define margem de filme; nao fornece lei de nucleacao de micropitting em ciclos.",
    }
    if lambda_current > 0.0 and lambda_min > 0.0 and math.isfinite(lambda_min):
        # h_min = C * F^(-0.073), demais entradas mantidas constantes.
        micro_capacity["load_factor_to_limit"] = (lambda_current / lambda_min) ** (1.0 / 0.073)
        micro_capacity["force_limit"] = float(c["force"]) * micro_capacity["load_factor_to_limit"]
    return macro_capacity, micro_capacity


def derive_v781_friction_thermal(c, hertz, ehl):
    """Atrito derivado de F, mu e v_slide; térmica a partir do cisalhamento do filme."""
    friction_force = float(c["mu"]) * float(c["force"])
    friction_power = friction_force * abs(float(c["v_slide"]))
    tau_lub = (
        float(c["eta0"])
        * abs(float(c["v_slide"]))
        / ehl["hmin"]
    )
    q_fric = tau_lub * abs(float(c["v_slide"]))
    p_fric_ehl = q_fric * hertz["contact_area"]

    # Térmica: fluxo derivado do cisalhamento do filme; potência global
    # também é reportada separadamente pela aproximação de Coulomb.
    thermal = thermal_state(
        T0=float(c["T0"]),
        heat_flux=q_fric,
        thermal_conductivity=float(c["thermal_conductivity"]),
        thermal_diffusivity=float(c["thermal_diffusivity"]),
        contact_time=float(c["contact_time"]),
        friction_power=p_fric_ehl,
        thermal_conductance=float(c["thermal_conductance"]),
    )
    return friction_force, friction_power, tau_lub, q_fric, p_fric_ehl, thermal


def derive_v781_wear(c, hertz, ehl, Rq):
    """Desgaste V7.78: Archard sobre a parcela de carga das asperezas."""
    sigma_bar = Rq / hertz["Rprime"]
    V_dimless = float(c["hardness"]) / hertz["Eprime"]
    asperity = point_contact_asperity_load_ratio(
        ehl["W"],
        ehl["U"],
        ehl["G"],
        sigma_bar,
        V_dimless,
    )
    mixed_wear = mixed_ehl_archard_wear_volume(
        float(c["wear_coefficient"]),
        float(c["force"]),
        asperity["La_fraction"],
        abs(float(c["v_slide"])) * float(c["sliding_time"]),
        float(c["hardness"]),
    )
    sliding_distance = abs(float(c["v_slide"])) * float(c["sliding_time"])
    wear_volume = mixed_wear["volume"]
    wear_depth = archard_wear_depth(
        wear_volume,
        hertz["contact_area"],
    )
    return asperity, mixed_wear, sliding_distance, wear_volume, wear_depth


def derive_v781_life_limits(c, stress_field, wear_depth):
    """Vida por fadiga (S-N do caso) e por desgaste, em ciclos e tempo."""
    fatigue_field = life_from_field(
        stress_field,
        float(c["fatigue_sigma_ref"]),
        float(c["fatigue_N_ref"]),
        float(c["fatigue_exponent_m"]),
    )
    fatigue_limit = life_extrema(fatigue_field)
    N_fatigue = (
        float(fatigue_limit["N_f"])
        if fatigue_limit is not None else float("inf")
    )

    frequency = float(c.get("cycle_frequency", 0.0))
    if frequency <= 0:
        raise ValidationError("A frequência de ciclos deve ser > 0 para calcular vida em ciclos.")

    t_fatigue = cycles_to_time(N_fatigue, frequency) if math.isfinite(N_fatigue) else float("inf")

    # O desgaste calculado acima corresponde ao tempo de referência
    # sliding_time. Sob carga e velocidade constantes, a taxa média é constante.
    t_slide_ref = float(c["sliding_time"])
    wear_depth_rate = wear_depth / t_slide_ref if t_slide_ref > 0 else 0.0
    t_wear = (
        float(c["wear_limit_depth"]) / wear_depth_rate
        if wear_depth_rate > 0 else float("inf")
    )
    N_wear = time_to_cycles(t_wear, frequency) if math.isfinite(t_wear) else float("inf")
    return N_fatigue, frequency, t_fatigue, wear_depth_rate, t_wear, N_wear


def derive_v781_instant_integrity(c, vm_max, ehl, Lambda, thermal, wear_depth):
    """Integridade imediata no instante do caso-base."""
    sigma_vm = float(vm_max["vm"])
    G_sigma = float(c["sigma_limit"]) - sigma_vm
    U_sigma = sigma_vm / float(c["sigma_limit"])
    U_hmin = float(c["hmin_limit"]) / ehl["hmin"]
    U_lambda = 1.0 / Lambda if Lambda > 0 else float("inf")
    U_temperature = thermal["T_total"] / float(c["temperature_limit"])
    U_wear = wear_depth / float(c["wear_limit_depth"])
    U_global = max(U_sigma, U_lambda, U_temperature, U_wear)
    return sigma_vm, G_sigma, U_sigma, U_hmin, U_lambda, U_temperature, U_wear, U_global


def derive_v781_limiting_life(U_global, N_fatigue, t_fatigue, N_wear, t_wear):
    """Vida global: zero se U>=1; senão a primeira falha acumulada."""
    if U_global >= 1.0:
        N_limit = 0.0
        t_limit = 0.0
        limiting_mechanism = "integridade instantanea"
    elif N_fatigue <= N_wear:
        N_limit = N_fatigue
        t_limit = t_fatigue
        limiting_mechanism = "fadiga"
    else:
        N_limit = N_wear
        t_limit = t_wear
        limiting_mechanism = "desgaste"
    return N_limit, t_limit, limiting_mechanism


def derive_v781_evaluation_state(
    c, frequency, N_fatigue, wear_depth_rate,
    U_sigma, U_lambda, U_temperature
):
    """Estado acumulado no número de ciclos de avaliação do usuário."""
    N_eval = float(c.get("cycles_evaluation", 0.0))
    if N_eval < 0:
        raise ValidationError("Os ciclos de avaliação devem ser >= 0.")
    t_eval = cycles_to_time(N_eval, frequency)
    D_fatigue_eval = (N_eval / N_fatigue) if math.isfinite(N_fatigue) else 0.0
    wear_depth_eval = wear_depth_rate * t_eval
    U_wear_eval = wear_depth_eval / float(c["wear_limit_depth"])
    U_global_eval = max(
        U_sigma,
        U_lambda,
        U_temperature,
        U_wear_eval,
        D_fatigue_eval,
    )
    return N_eval, t_eval, D_fatigue_eval, U_global_eval


def execute_v776_reference_case(case=None, progress=None):
    """Executa o caso-base V7.77 sem substituir entradas ausentes por valores."""
    c = dict(V776_BASE_CASE if case is None else case)
    hertz = derive_v776_hertz(c)
    pitting = derive_v781_pitting(c, hertz)

    # Campo espacial de tensões: malha moderada para o relatório-base.
    stress_field = spatial_full_tensor_field_cached(
        a=hertz["a"],
        b=hertz["b"],
        p0=hertz["p0"],
        nu=float(c["nu1"]),
        zmax=3.0 * hertz["a"],
        nx_xy=9,
        nz=12,
        integration_grid=32,
        progress=progress,
    )
    vm_max, sz_min = full_tensor_extrema(stress_field)

    # EHL: Rx/Ry são os raios relativos derivados do mesmo contato.
    ehl = ehl_film_thickness(
        float(c["eta0"]),
        float(c["ue"]),
        hertz["Eprime"],
        hertz["Rprime"],
        hertz["Rprime"],
        float(c["force"]),
        float(c["alpha_p"]),
    )
    Lambda, Rq = film_parameter(
        ehl["hmin"],
        float(c["Rq1"]),
        float(c["Rq2"]),
    )

    macro_capacity, micro_capacity = derive_v781_surface_capacity(
        c, pitting, ehl, Lambda, Rq
    )
    (
        friction_force, friction_power, tau_lub, q_fric,
        p_fric_ehl, thermal
    ) = derive_v781_friction_thermal(c, hertz, ehl)
    (
        asperity, mixed_wear, sliding_distance,
        wear_volume, wear_depth
    ) = derive_v781_wear(c, hertz, ehl, Rq)
    (
        N_fatigue, frequency, t_fatigue,
        wear_depth_rate, t_wear, N_wear
    ) = derive_v781_life_limits(c, stress_field, wear_depth)
    (
        sigma_vm, G_sigma, U_sigma, U_hmin,
        U_lambda, U_temperature, U_wear, U_global
    ) = derive_v781_instant_integrity(
        c, vm_max, ehl, Lambda, thermal, wear_depth
    )
    N_limit, t_limit, limiting_mechanism = derive_v781_limiting_life(
        U_global, N_fatigue, t_fatigue, N_wear, t_wear
    )
    N_eval, t_eval, D_fatigue_eval, U_global_eval = derive_v781_evaluation_state(
        c, frequency, N_fatigue, wear_depth_rate,
        U_sigma, U_lambda, U_temperature
    )

    return {
        "version": V780_VERSION,
        "inputs": c,
        "equivalent": hertz,
        "hertz": hertz,
        "stress": {
            "field_points": len(stress_field),
            "sigma_VM_max": sigma_vm,
            "sigma_VM_point": {
                "x": vm_max["x"], "y": vm_max["y"], "z": vm_max["z"]
            },
            "sigma_z_min": sz_min["tensor"]["sz"],
            "sigma_z_point": {
                "x": sz_min["x"], "y": sz_min["y"], "z": sz_min["z"]
            },
        },
        "ehl": {
            **ehl,
            "Rq": Rq,
            "Lambda": Lambda,
        },
        "friction": {
            "friction_force": friction_force,
            "friction_power_coulomb": friction_power,
            "tau_lub": tau_lub,
            "q_fric": q_fric,
            "P_fric_ehl": p_fric_ehl,
        },
        "thermal": thermal,
        "pitting": pitting,
        "surface_capacity": {
            "macro": macro_capacity,
            "micro": micro_capacity,
        },
        "wear": {
            "sliding_distance": sliding_distance,
            "normal_force": float(c["force"]),
            "asperity_load_ratio_percent": asperity["La_percent"],
            "asperity_load_fraction": asperity["La_fraction"],
            "asperity_force": mixed_wear["asperity_force"],
            "volume": wear_volume,
            "depth": wear_depth,
        },
        "integrity": {
            "G_sigma": G_sigma,
            "U_sigma": U_sigma,
            "U_hmin": U_hmin,
            "U_lambda": U_lambda,
            "U_temperature": U_temperature,
            "U_wear": U_wear,
            "U_global": U_global,
            "U_global_eval": U_global_eval,
        },
        "fatigue": {
            "sigma_ref": float(c["fatigue_sigma_ref"]),
            "N_ref": float(c["fatigue_N_ref"]),
            "exponent_m": float(c["fatigue_exponent_m"]),
            "N_f": N_fatigue,
            "t_f": t_fatigue,
            "D_eval": D_fatigue_eval,
        },
        "life": {
            "frequency": frequency,
            "N_eval": N_eval,
            "t_eval": t_eval,
            "wear_depth_rate": wear_depth_rate,
            "N_wear": N_wear,
            "t_wear": t_wear,
            "N_fatigue": N_fatigue,
            "t_fatigue": t_fatigue,
            "N_lim": N_limit,
            "t_lim": t_limit,
            "limiting_mechanism": limiting_mechanism,
        },
        "status": "calculado",
    }


def _v782_table(title, headers, cols):
    out = ["", title, "  ".join("%-14s" % h for h in headers)]
    for row in zip(*cols):
        out.append("  ".join("%-14s" % ("%.6e" % v) for v in row))
    return out


def format_v782_results_txt(result, data, timing=""):
    """Relatório completo + tabelas numéricas dos gráficos, para gravação em TXT."""
    lines = [format_v776_reference_report(result)]
    if timing:
        lines += ["", "Tempo computacional: " + timing]
    lines += ["", "=== TABELAS DOS GRÁFICOS ==="]

    def put(key, fn):
        q = data.get(key)
        if q is None:
            return
        if "erro" in q:
            lines.append("")
            lines.append("[%s] erro: %s" % (key, q["erro"]))
            return
        lines.extend(fn(q))

    put("p", lambda q: _v782_table(
        "1 — Pressão de contato (p_m = %.6e MPa)" % q["pm"],
        ["r [mm]", "p [MPa]"], [q["x"], q["p"]]))
    put("sz", lambda q: _v782_table(
        "2 — Tensões no eixo × profundidade (σ_lim = %.6e MPa)" % q["lim"],
        ["z [mm]", "σz [MPa]", "σr [MPa]", "σVM [MPa]", "τmax [MPa]"],
        [q["z"], q["sz"], q["sr"], q["vm"], q["tm"]]))

    def film_rows(q):
        heads, cols = ["ue [m/s]", "h_min [µm]", "h_c [µm]"], [q["u"], q["hmin"], q["hc"]]
        if q["lam"]:
            heads.append("Λ [-]")
            cols.append(q["lam"])
        return _v782_table(
            "3/4 — Filme e Λ × velocidade (h_min lim = %.6e µm)" % q["hmin_lim"], heads, cols)

    put("h", film_rows)
    put("T", lambda q: _v782_table(
        "5 — Temperatura × velocidade de deslizamento (T lim = %.6e °C)" % q["lim"],
        ["v_slide [m/s]", "T_bulk [°C]", "T_flash [°C]", "T_total [°C]"],
        [q["v"], q["tb"], q["tf"], q["tt"]]))
    put("U", lambda q: ["", "6 — Integridade (limite = 1)"]
        + ["%-16s %.6e" % (k, v) for k, v in q["items"]])
    put("sn", lambda q: _v782_table(
        "7 — Curva S-N",
        ["σ [MPa]", "N [ciclos]"], [q["s"], q["n"]]))

    def force_rows(q):
        heads, cols = ["F [N]", "U_sigma", "U_global"], [q["f"], q["usig"], q["uglo"]]
        if q["sh"]:
            heads.append("S_H")
            cols.append(q["sh"])
        return _v782_table("8 — Margens × força", heads, cols)

    put("F", force_rows)
    put("N", lambda q: _v782_table(
        "9 — Evolução por ciclos (N_lim = %.6e ciclos, t_lim = %.6e s)" % (q["nlim"], q["tlim"]),
        ["N/N_lim", "D fadiga", "U_desgaste", "U_global"],
        [q["x"], q["d"], q["uw"], q["ug"]]))
    return "\n".join(lines) + "\n"


S01 = "01 — ARQUITETURA E FLUXO DO CÁLCULO"
S02 = "02 — ENTRADAS E UNIDADES"
S03 = "03 — CONTATO HERTZIANO"
S04 = "04 — CAMPO DE TENSÕES"
S05 = "05 — LUBRIFICAÇÃO EHL"
S06 = "06 — ATRITO E TÉRMICA"
S07 = "07 — DESGASTE"
S08 = "08 — FADIGA E PITTING"
S09 = "09 — INTEGRIDADE E VIDA"
S10 = "10 — GRÁFICOS"
S11 = "11 — PREMISSAS E LIMITAÇÕES"
S12 = "12 — REFERÊNCIAS BIBLIOGRÁFICAS"

GLOSSARY_TERMS = [
(S01, "Fluxo do cálculo", "O caso é calculado em cadeia: Hertz (E', R', a, p0, pm) → campo de tensões → EHL (h_min, h_c, Λ) → atrito e térmica → desgaste → vida por fadiga e por desgaste → integridade (U). Cada etapa usa apenas as entradas e os resultados das anteriores. Entradas ausentes não são substituídas por valores."),
(S01, "Sistema de unidades", "A interface usa unidades práticas (mm, GPa, µm, °C). Tudo é convertido para SI antes do cálculo; o relatório e as tabelas do TXT usam SI, exceto onde a unidade vem indicada."),
(S01, "Execução em segundo plano", "O cálculo roda numa thread separada. O relatório aparece primeiro e cada gráfico surge assim que fica pronto; a página continua rolando durante a simulação. O botão EXECUTAR fica desativado até o fim."),
(S01, "Cache do campo de tensões", "σ/p0 em função de (x/a, y/a, z/a) depende só de ν1, de b/a e da malha. O campo é calculado uma vez e reaproveitado quando só F, uₑ, v_slide, T₀ e demais entradas mudam; os valores são reescalados por a e p0. A primeira execução, ou uma mudança de ν1, recalcula o campo."),
(S01, "Tempo computacional", "Tempo medido no aparelho: o caso (campo de tensões e demais etapas) e os gráficos. Com o cache ativo, o caso fica bem mais curto na segunda execução."),
(S01, "Gravação em TXT", "Grava o relatório completo e as tabelas numéricas dos 9 gráficos em Download/Tribologico_resultado_V781.txt (Android). O arquivo é sobrescrito a cada gravação."),
(S01, "Caso de referência", "V776_BASE_CASE: conjunto de entradas padrão do modelo. A interface parte dele e substitui os campos preenchidos. Valores marcados [ESTIMADA] no relatório não vêm de medição."),
(S02, "F — força normal [N]", "Carga que comprime os dois corpos. Define a, p0 e pm (a ∝ F^(1/3)) e a parcela de carga que as asperezas suportam."),
(S02, "R1, R2 — raios dos corpos [mm]", "Raios de curvatura dos dois corpos. Não são os semi-eixos a e b do contato. Um raio muito grande (como 1 000 000 mm) representa uma superfície praticamente plana."),
(S02, "E1, E2 — módulos de elasticidade [GPa]", "Rigidez de cada corpo. Entram em E'."),
(S02, "ν1, ν2 — coeficientes de Poisson [-]", "Entram em E'. O campo de tensões usa ν1. Devem estar em 0 ≤ ν < 0,5."),
(S02, "uₑ — velocidade média de arrastamento [m/s]", "Velocidade com que o lubrificante é arrastado para o contato. Alimenta o parâmetro U e, portanto, h_min e h_c."),
(S02, "v_slide — velocidade de deslizamento [m/s]", "Velocidade relativa entre as superfícies. Alimenta o cisalhamento do filme, a potência de atrito e a distância de deslizamento do desgaste."),
(S02, "η₀ — viscosidade [Pa·s]", "Viscosidade do lubrificante à pressão ambiente. Entra em U e no cisalhamento do filme."),
(S02, "αp — coeficiente pressão-viscosidade [Pa⁻¹]", "Quanto a viscosidade cresce com a pressão. Entra no parâmetro G = αp·E'."),
(S02, "Rq1, Rq2 — rugosidades RMS [µm]", "Rugosidade quadrática média de cada superfície. Compõem Rq = √(Rq1² + Rq2²) para Λ e para o desgaste."),
(S02, "T₀ — temperatura de referência [°C]", "Temperatura de partida. As temperaturas do modelo são T₀ mais as elevações bulk e flash."),
(S02, "μ — coeficiente de atrito [-]", "Usado na potência de atrito de Coulomb (μ·F·v_slide), reportada separadamente. A térmica usa o cisalhamento do filme, não μ."),
(S02, "β — coeficiente térmico [K⁻¹]", "Entrada do caso, mas não entra em nenhuma equação do modelo atual."),
(S02, "k, α — condutividade e difusividade térmicas", "k [W/(m·K)] e α [m²/s]. Entram na elevação de temperatura flash."),
(S02, "Tempo de contato [s]", "Duração do aquecimento superficial usada na temperatura flash. É uma entrada fixa: não depende da velocidade."),
(S02, "Condutância térmica H [W/K]", "Condutância equivalente do sistema para dissipar a potência de atrito. Define a elevação bulk."),
(S02, "K — coeficiente de desgaste [-]", "Coeficiente adimensional da lei de Archard. Com K = 0 não há desgaste e U_desgaste fica em zero."),
(S02, "H — dureza [GPa]", "Dureza do material que desgasta. Entra na lei de Archard e em V = H/E' da correlação de asperezas."),
(S02, "Tempo de deslizamento [s]", "Tempo ao qual o desgaste calculado se refere. A taxa de desgaste é a profundidade dividida por esse tempo."),
(S02, "Limites: σ, temperatura, h_min, desgaste", "Valores admissíveis informados por você. Cada um define o denominador de um índice U. O modelo não os inventa."),
(S02, "Frequência de ciclos [Hz]", "Converte ciclos em tempo: t = N/f."),
(S02, "N de avaliação [ciclos]", "Número de ciclos em que o estado acumulado é mostrado (dano de fadiga e U_global_eval). Não é a vida."),
(S02, "σ_ref, N_ref, m — curva S-N", "Ponto de referência (σ_ref, N_ref) e expoente m da curva de fadiga. Todos são entradas do caso."),
(S02, "σHlim e S_H,min — pitting", "Tensão limite de contato (σHlim) e margem mínima admissível S_H,min do critério de pitting."),
(S03, "E' — módulo equivalente", "1/E' = (1 − ν1²)/E1 + (1 − ν2²)/E2. Representa a rigidez conjunta dos dois corpos."),
(S03, "R' — raio equivalente", "1/R' = 1/R1 + 1/R2. O contato circular equivale a uma esfera de raio R' sobre um plano."),
(S03, "a = b — semi-eixo de contato", "a = (3·F·R'/(4·E'))^(1/3). O contato é circular, então b = a."),
(S03, "A — área de contato", "A = π·a·b."),
(S03, "p0 — pressão máxima de Hertz", "p0 = 3F/(2A), no centro do contato. A distribuição é p(r) = p0·√(1 − (r/a)²)."),
(S03, "pm — pressão média", "pm = F/A = (2/3)·p0. É a tensão σHm usada no critério de pitting."),
(S04, "Campo de tensões (Boussinesq)", "O tensor de tensões é a soma numérica das soluções de Boussinesq de carga pontual, integrando a pressão de Hertz em coordenadas elípticas sobre a área de contato. Considera só a pressão normal; não inclui a tração superficial do atrito."),
(S04, "σz e σr no eixo", "Tensões normais ao longo do eixo do contato. σz é de compressão (negativa). Sobre o eixo σr = σx = σy."),
(S04, "σVM — tensão de von Mises", "Tensão equivalente do tensor completo. σVM_max é o maior valor do campo e alimenta U_sigma, G_sigma e a vida por fadiga."),
(S04, "τmax — cisalhamento máximo", "Mostrado no gráfico 2 como |σz − σr|/2 sobre o eixo, onde esses são os principais."),
(S04, "Malha do campo", "O relatório usa 9×9 pontos em x–y, 12 níveis em z (até 3a) e integração 32×32. O máximo é amostrado nessa malha, então sua posição em z tem passo de 0,25a."),
(S05, "U, G, W — grupos adimensionais", "U = η0·uₑ/(E'·R'), G = αp·E' e W = F/(E'·R'²). Entram nas correlações de filme e na razão de carga das asperezas."),
(S05, "h_min — espessura mínima do filme", "Correlação de Hamrock–Dowson: h_min/R = 3,63·U^0,68·G^0,49·W^−0,073·(1 − e^(−0,68k)), com k = 1,03·(Ry/Rx)^0,64. Cresce com a velocidade e quase não depende da carga (W^−0,073)."),
(S05, "h_c — espessura central do filme", "h_c/R = 2,69·U^0,67·G^0,53·W^−0,067·(1 − 0,61·e^(−0,73k)). Sempre maior que h_min."),
(S05, "Rq — rugosidade composta", "Rq = √(Rq1² + Rq2²)."),
(S05, "Λ — parâmetro de filme", "Λ = h_min/Rq. Quanto menor, mais as asperezas se tocam. Se Rq = 0, Λ é indefinido (infinito)."),
(S05, "Λ mín", "Λ mín = h_min,lim/Rq. É o valor de Λ em que o filme atinge o limite informado."),
(S06, "τ_lub — cisalhamento do filme", "τ = η0·v_slide/h_min (filme newtoniano com viscosidade η0)."),
(S06, "q — fluxo de calor por atrito", "q = τ·v_slide. Alimenta a temperatura flash."),
(S06, "P_fric — potência de atrito no filme", "P = q·A. Alimenta a temperatura bulk."),
(S06, "ΔT_flash", "ΔT = (q/k)·√(α·t/π), forma simplificada para fonte superficial transitória. Não substitui uma solução térmica completa."),
(S06, "ΔT_bulk", "ΔT = P_fric/H, modelo concentrado."),
(S06, "T_total", "T_total = T₀ + ΔT_bulk + ΔT_flash. O gráfico 5 mostra T_bulk = T₀ + ΔT_bulk, T_flash = T₀ + ΔT_flash e T_total, em °C."),
(S06, "Atrito de Coulomb", "Força μ·F e potência μ·F·v_slide. Reportado à parte; não alimenta a térmica."),
(S07, "L_a — razão de carga das asperezas", "Correlação de Masjedi & Khonsari (2015, Eq. 23): fração da carga suportada pelas asperezas no regime misto. Usa W, U, G, σ̄ = Rq/R' e V = H/E'."),
(S07, "Lei de Archard modificada", "F_asp = L_a·F e V_w = K·F_asp·s/H, com s = v_slide·t_deslizamento. A força aplicada à lei é só a parcela das asperezas."),
(S07, "Profundidade de desgaste", "d = V_w/A. A taxa de desgaste é d/t_deslizamento e é tratada como constante; t_desgaste = limite/taxa."),
(S08, "Curva S-N de potência", "N = N_ref·(σ/σ_ref)^(−1/m). Usa σ = σVM em cada ponto do campo."),
(S08, "N_f — vida por fadiga", "Menor N_f entre os pontos do campo (ponto crítico). Em geral coincide com o ponto de σVM máxima."),
(S08, "D — dano de fadiga (Miner)", "D = N/N_f, regra linear. D ≥ 1 indica falha por fadiga."),
(S08, "S_H — segurança contra pitting", "S_H = σHlim/σHm, com σHm = pm. Admissível se S_H ≥ S_H,min. U_pitting = S_H,min/S_H."),
(S08, "Fator de carga até o limite", "Pitting: (S_H/S_H,min)³, exato aqui porque pm ∝ F^(1/3). Micropitting: (Λ/Λ mín)^(1/0,073), porque h_min ∝ F^(−0,073)."),
(S09, "U — índice de utilização", "Razão entre a solicitação e o limite. U ≥ 1 significa limite atingido."),
(S09, "U_sigma e G_sigma", "U_sigma = σVM_max/σ_lim. G_sigma = σ_lim − σVM_max (margem em Pa)."),
(S09, "U_hmin e U_lambda", "U_hmin = h_min,lim/h_min. U_lambda = 1/Λ, ou seja, o limite é Λ = 1."),
(S09, "U_temperatura e U_desgaste", "U_temperatura = T_total/T_lim (em kelvin). U_desgaste = profundidade/limite de desgaste."),
(S09, "U_global", "Maior entre U_sigma, U_lambda, U_temperatura e U_desgaste. U_hmin e U_pitting são mostrados mas não entram no máximo."),
(S09, "U_global_eval", "Como U_global, mas com o desgaste e o dano de fadiga acumulados até N de avaliação."),
(S09, "N_lim e mecanismo limitante", "Se U_global ≥ 1: N_lim = 0 (integridade instantânea). Senão, o menor entre N de fadiga e N de desgaste. t_lim = N_lim/f."),
(S10, "1 — Pressão de contato", "p(r) de −a a +a, com a linha de pm."),
(S10, "2 — Tensões × profundidade", "σz, σr, σVM e τmax no eixo de z = 0,1a a 3a, com o limite σ_lim. Integração com malha mais fina que a do relatório."),
(S10, "3 e 4 — Filme e Λ × velocidade", "uₑ varia de 0,2 a 2 vezes o valor da entrada; o resto fica como na entrada. O gráfico 4 mostra a linha Λ = 1 e, se houver limite de filme, Λ mín."),
(S10, "5 — Temperatura × velocidade", "v_slide varia de 0,1 a 2 vezes a entrada, com o mesmo filme (h_min) da entrada e o mesmo tempo de contato."),
(S10, "6 — Integridade", "Barras dos índices U e do dano de fadiga no estado atual, com a linha de limite em 1."),
(S10, "7 — Curva S-N", "N × σ em escala logarítmica, com o ponto de operação (σVM_max, N_f) e o ponto de referência."),
(S10, "8 — Margens × força", "F varia de 0,2 a 3 vezes a entrada. σVM_max escala com p0 (exato para o mesmo ν); filme, térmica e desgaste são recalculados com as funções do modelo."),
(S10, "9 — Evolução por ciclos", "N de 0 a 1,2·N_lim, com carga e velocidade constantes. O dano D = N/N_f e o U_desgaste crescem linearmente; U_global é o maior entre os índices estáticos, o dano e o desgaste. O primeiro a chegar em 1 define o mecanismo limitante."),
(S11, "Correlações empíricas", "As expressões de h_min, h_c e L_a vêm de ajustes da literatura. O modelo não verifica se as entradas estão dentro das faixas em que foram obtidas."),
(S11, "Térmica simplificada", "Flash e bulk são fórmulas fechadas. Não há acoplamento térmico–EHL: a viscosidade não muda com a temperatura calculada."),
(S11, "Desgaste constante", "A taxa de desgaste é constante no tempo, sob carga e velocidade constantes. Não há evolução do perfil do contato."),
(S11, "Tensões só da pressão normal", "O campo ignora a tração superficial do atrito. Em contatos com μ elevado, as tensões reais próximas à superfície são maiores."),
(S11, "Sem interpolação", "O máximo de σVM vem de pontos de malha; o valor exato fica ligeiramente acima ou ao lado do reportado."),
(S12, "Hamrock & Dowson (1977)", "Isothermal elastohydrodynamic lubrication of point contacts. ASME Journal of Lubrication Technology, 99. Origem das correlações de h_min e h_c."),
(S12, "Masjedi & Khonsari (2015)", "Tribology International 82, 228–244. DOI 10.1016/j.triboint.2014.09.010. Origem da correlação de L_a."),
(S12, "Archard (1953)", "Contact and rubbing of flat surfaces. Journal of Applied Physics 24, 981–988."),
(S12, "Miner (1945)", "Cumulative damage in fatigue. Journal of Applied Mechanics 12, A159–A164."),
(S12, "Johnson (1985)", "Contact Mechanics. Cambridge University Press. Contato de Hertz e solução de Boussinesq."),
]


def _v782_linspace(lo, hi, n):
    return [lo + (hi - lo) * i / (n - 1) for i in range(n)]


def build_v782_chart_data(result, publish=None):
    """
    Séries dos gráficos a partir do resultado e das entradas do caso.
    Reutiliza as funções do modelo; cada gráfico é isolado: um erro em um
    não impede os demais.
    """
    c = result["inputs"]
    h = result["hertz"]
    a, b, p0, pm = h["a"], h["b"], h["p0"], h["pm"]
    nu = float(c["nu1"])
    F = float(c["force"])
    ue = float(c["ue"])
    v_slide = abs(float(c["v_slide"]))
    data = {}

    def section(key, fn):
        try:
            data[key] = fn()
        except Exception as exc:
            data[key] = {"erro": str(exc)}
        if publish is not None:
            publish(key, data[key])

    def pressure():
        xs = _v782_linspace(-a, a, 61)
        return {
            "x": [x * 1e3 for x in xs],
            "p": [p0 * math.sqrt(max(0.0, 1.0 - (x / a) ** 2)) / 1e6 for x in xs],
            "pm": pm / 1e6,
        }

    def axis_stress():
        zs = _v782_linspace(0.1 * a, 3.0 * a, 30)
        sz, sr, vm, tm = [], [], [], []
        for z in zs:
            t = boussinesq_distributed_tensor(0.0, 0.0, z, a, b, p0, nu, nx=48)
            sz.append(t["sz"] / 1e6)
            sr.append(t["sx"] / 1e6)
            vm.append(von_mises_from_tensor(t) / 1e6)
            tm.append(abs(t["sz"] - t["sx"]) / 2.0 / 1e6)
        return {
            "z": [z * 1e3 for z in zs],
            "sz": sz, "sr": sr, "vm": vm, "tm": tm,
            "lim": float(c["sigma_limit"]) / 1e6,
        }

    def film():
        us = _v782_linspace(0.2 * ue, 2.0 * ue, 25)
        hmin, hc, lam = [], [], []
        for u in us:
            e = ehl_film_thickness(
                float(c["eta0"]), u, h["Eprime"], h["Rprime"], h["Rprime"],
                F, float(c["alpha_p"])
            )
            hmin.append(e["hmin"] * 1e6)
            hc.append(e["hc"] * 1e6)
            L, Rq = film_parameter(e["hmin"], float(c["Rq1"]), float(c["Rq2"]))
            lam.append(L)
        finite = all(math.isfinite(L) for L in lam)
        hlim = float(c["hmin_limit"])
        Rq = result["ehl"]["Rq"]
        return {
            "u": us, "hmin": hmin, "hc": hc,
            "lam": lam if finite else [],
            "hmin_lim": hlim * 1e6,
            "lam_min": (hlim / Rq) if (Rq > 0 and hlim > 0) else 0.0,
        }

    def thermal():
        vs = _v782_linspace(0.1 * v_slide, 2.0 * v_slide, 25)
        tb, tf, tt = [], [], []
        for v in vs:
            th = derive_v781_friction_thermal(
                dict(c, v_slide=v), h, result["ehl"]
            )[5]
            tb.append(th["T_bulk"] - 273.15)
            tf.append(th["T_flash"] - 273.15)
            tt.append(th["T_total"] - 273.15)
        return {
            "v": vs, "tb": tb, "tf": tf, "tt": tt,
            "lim": float(c["temperature_limit"]) - 273.15,
        }

    def integrity():
        g = result["integrity"]
        items = [
            ("U_sigma", g["U_sigma"]),
            ("U_hmin", g["U_hmin"]),
            ("U_lambda", g["U_lambda"]),
            ("U_temperatura", g["U_temperature"]),
            ("U_desgaste", g["U_wear"]),
        ]
        if result["pitting"].get("enabled"):
            items.append(("U_pitting", result["pitting"]["U_pitting"]))
        items.append(("D_fadiga", result["fatigue"]["D_eval"]))
        items.append(("U_global", g["U_global"]))
        return {"items": items}

    def sn():
        f = result["fatigue"]
        sref, nref, m = f["sigma_ref"], f["N_ref"], f["exponent_m"]
        svm = result["stress"]["sigma_VM_max"]
        lo = 0.5 * svm
        hi = max(1.2 * sref, 1.5 * svm)
        sg = _v782_linspace(lo, hi, 30)
        return {
            "s": [x / 1e6 for x in sg],
            "n": [sn_cycles_from_stress(x, sref, nref, m) for x in sg],
            "op": (svm / 1e6, f["N_f"]) if math.isfinite(f["N_f"]) else None,
            "ref": (sref / 1e6, nref),
        }

    def force_sweep():
        if F <= 0 or p0 <= 0:
            raise ValidationError("Força nula: sem varredura.")
        c_vm = result["stress"]["sigma_VM_max"] / p0
        fs = _v782_linspace(0.2 * F, 3.0 * F, 16)
        usig, uglo, sh = [], [], []
        for Fi in fs:
            cF = dict(c, force=Fi)
            hF = derive_v776_hertz(cF)
            e = ehl_film_thickness(
                float(c["eta0"]), ue, hF["Eprime"], hF["Rprime"], hF["Rprime"],
                Fi, float(c["alpha_p"])
            )
            L, Rq = film_parameter(e["hmin"], float(c["Rq1"]), float(c["Rq2"]))
            th = derive_v781_friction_thermal(cF, hF, e)[5]
            wd = derive_v781_wear(cF, hF, e, Rq)[4]
            u = derive_v781_instant_integrity(
                cF, {"vm": c_vm * hF["p0"]}, e, L, th, wd
            )
            usig.append(u[2])
            uglo.append(u[7])
            pit = derive_v781_pitting(cF, hF)
            sh.append(pit["S_H"] if pit.get("enabled") else None)
        return {
            "f": fs, "usig": usig, "uglo": uglo,
            "sh": sh if all(x is not None and math.isfinite(x) for x in sh) else [],
        }

    def evolution():
        life = result["life"]
        N_lim = life["N_lim"]
        if not (math.isfinite(N_lim) and N_lim > 0):
            raise ValidationError(
                "Vida limite nula ou infinita: sem evolução por ciclos."
            )
        g = result["integrity"]
        xs = _v782_linspace(0.0, 1.2, 25)
        d, ug, uw = [], [], []
        for x in xs:
            _, t_e, D, Ug = derive_v781_evaluation_state(
                dict(c, cycles_evaluation=x * N_lim),
                life["frequency"], life["N_fatigue"], life["wear_depth_rate"],
                g["U_sigma"], g["U_lambda"], g["U_temperature"]
            )
            d.append(D)
            ug.append(Ug)
            uw.append(life["wear_depth_rate"] * t_e / float(c["wear_limit_depth"]))
        return {"x": xs, "d": d, "ug": ug, "uw": uw,
                "nlim": N_lim, "tlim": life["t_lim"]}

    section("p", pressure)
    section("sz", axis_stress)
    section("h", film)
    section("T", thermal)
    section("U", integrity)
    section("sn", sn)
    section("F", force_sweep)
    section("N", evolution)
    return data


def _v776_fmt(x):
    return f"{float(x):.8e}"


def format_v776_reference_report(result):
    """Relatório integral em texto puro, sem perda de dados numéricos."""
    c = result["inputs"]
    h = result["hertz"]
    s = result["stress"]
    e = result["ehl"]
    f = result["friction"]
    t = result["thermal"]
    w = result["wear"]
    g = result["integrity"]
    cap = result["surface_capacity"]
    macro = cap["macro"]
    micro = cap["micro"]

    lines = [
        "=== MODELO TRIBOLÓGICO INTEGRADO V7.81 ===",
        "=== CASO FÍSICO DE REFERÊNCIA ===",
        "",
        "ENTRADAS PRIMÁRIAS",
        f"F                  = {_v776_fmt(c['force'])} N",
        f"R1                 = {_v776_fmt(c['R1'])} m",
        f"R2                 = {_v776_fmt(c['R2'])} m",
        f"E1                 = {_v776_fmt(c['E1'])} Pa",
        f"E2                 = {_v776_fmt(c['E2'])} Pa",
        f"nu1                = {_v776_fmt(c['nu1'])}",
        f"nu2                = {_v776_fmt(c['nu2'])}",
        "",
        "PROPRIEDADES EQUIVALENTES",
        f"E'                 = {_v776_fmt(h['Eprime'])} Pa",
        f"R'                 = {_v776_fmt(h['Rprime'])} m",
        "",
        "CONTATO HERTZIANO",
        f"a                  = {_v776_fmt(h['a'])} m",
        f"b                  = {_v776_fmt(h['b'])} m",
        f"area               = {_v776_fmt(h['contact_area'])} m2",
        f"p_media            = {_v776_fmt(h['pm'])} Pa",
        f"p_H_max            = {_v776_fmt(h['p0'])} Pa",
        "",
        "CAMPO DE TENSOES",
        f"pontos             = {s['field_points']}",
        f"sigma_VM_max       = {_v776_fmt(s['sigma_VM_max'])} Pa",
        f"VM_x               = {_v776_fmt(s['sigma_VM_point']['x'])} m",
        f"VM_y               = {_v776_fmt(s['sigma_VM_point']['y'])} m",
        f"VM_z               = {_v776_fmt(s['sigma_VM_point']['z'])} m",
        f"sigma_z_min        = {_v776_fmt(s['sigma_z_min'])} Pa",
        "",
        "EHL",
        f"ue                 = {_v776_fmt(c['ue'])} m/s   [ESTIMADO]",
        f"eta0               = {_v776_fmt(c['eta0'])} Pa.s [ESTIMADO]",
        f"alpha_p            = {_v776_fmt(c['alpha_p'])} Pa^-1 [ESTIMADO]",
        f"U                  = {_v776_fmt(e['U'])}",
        f"G                  = {_v776_fmt(e['G'])}",
        f"W                  = {_v776_fmt(e['W'])}",
        f"k                  = {_v776_fmt(e['k'])}",
        f"h_min              = {_v776_fmt(e['hmin'])} m",
        f"h_c                = {_v776_fmt(e['hc'])} m",
        f"Rq                 = {_v776_fmt(e['Rq'])} m",
        f"Lambda             = {_v776_fmt(e['Lambda'])}",
        "",
        "ATRITO / CISALHAMENTO",
        f"mu                 = {_v776_fmt(c['mu'])} [ESTIMADO]",
        f"F_t                = {_v776_fmt(f['friction_force'])} N",
        f"P_fric_Coulomb     = {_v776_fmt(f['friction_power_coulomb'])} W",
        f"tau_lub            = {_v776_fmt(f['tau_lub'])} Pa",
        f"q_fric             = {_v776_fmt(f['q_fric'])} W/m2",
        f"P_fric_EHL         = {_v776_fmt(f['P_fric_ehl'])} W",
        "",
        "TERMICA",
        f"T0                 = {_v776_fmt(c['T0'])} K [ESTIMADO]",
        f"DeltaT_bulk        = {_v776_fmt(t['dT_bulk'])} K",
        f"DeltaT_flash       = {_v776_fmt(t['dT_flash'])} K",
        f"T_total            = {_v776_fmt(t['T_total'])} K",
        "",
        "DESGASTE — V7.78 / EHL MISTO",
        f"F_normal           = {_v776_fmt(w['normal_force'])} N",
        f"L_a                = {_v776_fmt(w['asperity_load_ratio_percent'])} %",
        f"L_a(fracao)        = {_v776_fmt(w['asperity_load_fraction'])}",
        f"F_asp              = {_v776_fmt(w['asperity_force'])} N",
        f"sliding_distance   = {_v776_fmt(w['sliding_distance'])} m",
        f"V_w                = {_v776_fmt(w['volume'])} m3",
        f"d_w                = {_v776_fmt(w['depth'])} m",
        "",
        "PITTING — V7.79",
        f"sigma_Hm          = {_v776_fmt(result['pitting']['sigma_Hm'])} Pa",
        f"sigma_Hlim        = {('NAO INFORMADO' if result['pitting']['sigma_Hlim'] is None else _v776_fmt(result['pitting']['sigma_Hlim']) + ' Pa')}",
        f"S_H               = {('NAO CALCULADO' if result['pitting']['S_H'] is None else _v776_fmt(result['pitting']['S_H']))}",
        f"S_H_min           = {_v776_fmt(result['pitting']['S_H_min'])}",
        f"U_pitting         = {('NAO CALCULADO' if result['pitting']['U_pitting'] is None else _v776_fmt(result['pitting']['U_pitting']))}",
        f"fator_carga_limite = {('NAO CALCULADO' if result['pitting']['S_H'] is None else _v776_fmt(result['pitting']['load_factor_to_limit']))}",
        f"pitting_admissible= {result['pitting']['admissible']}",
        "",
        "MACROPITTING — MARGEM DE CARGA",
        f"fator_carga_macro  = {('NAO CALCULADO' if macro['load_factor_to_limit'] is None else _v776_fmt(macro['load_factor_to_limit']))}",
        f"forca_limite_macro = {('NAO CALCULADO' if macro['force_limit'] is None else _v776_fmt(macro['force_limit']))} N",
        "",
        "MICROPITTING — MARGEM DE FILME",
        f"h_min              = {_v776_fmt(micro['h_min'])} m",
        f"h_min_limite       = {_v776_fmt(micro['h_min_limit'])} m",
        f"Rq_composto        = {_v776_fmt(micro['Rq'])} m",
        f"Lambda             = {_v776_fmt(micro['Lambda'])}",
        f"Lambda_min         = {_v776_fmt(micro['Lambda_min'])}",
        f"S_lambda           = {_v776_fmt(micro['S_lambda'])}",
        f"micro_admissible   = {micro['admissible']}",
        f"fator_carga_micro  = {('NAO CALCULADO' if micro['load_factor_to_limit'] is None else _v776_fmt(micro['load_factor_to_limit']))}",
        f"forca_limite_micro = {('NAO CALCULADO' if micro['force_limit'] is None else _v776_fmt(micro['force_limit']))} N",
        "vida_micro         = NAO PREVISTA PELO MODELO",
        "tempo_1o_micro      = NAO PREVISTO PELO MODELO",
        "",
        "INTEGRIDADE",
        f"U_sigma            = {_v776_fmt(g['U_sigma'])}",
        f"U_hmin             = {_v776_fmt(g['U_hmin'])}",
        f"U_Lambda           = {_v776_fmt(g['U_lambda'])}",
        f"U_temperature      = {_v776_fmt(g['U_temperature'])}",
        f"U_wear             = {_v776_fmt(g['U_wear'])}",
        f"U_global           = {_v776_fmt(g['U_global'])}",
        f"U_global_N_eval    = {_v776_fmt(g['U_global_eval'])}",
        "",
        "FADIGA / VIDA",
        f"frequencia         = {_v776_fmt(result['life']['frequency'])} Hz",
        f"N_avaliacao        = {_v776_fmt(result['life']['N_eval'])} ciclos",
        f"t_avaliacao        = {_v776_fmt(result['life']['t_eval'])} s",
        f"N_fadiga           = {_v776_fmt(result['life']['N_fatigue'])} ciclos",
        f"t_fadiga           = {_v776_fmt(result['life']['t_fatigue'])} s",
        f"D_fadiga(N_eval)   = {_v776_fmt(result['fatigue']['D_eval'])}",
        "",
        "VIDA POR DESGASTE",
        f"taxa desgaste      = {_v776_fmt(result['life']['wear_depth_rate'])} m/s",
        f"N_desgaste         = {_v776_fmt(result['life']['N_wear'])} ciclos",
        f"t_desgaste         = {_v776_fmt(result['life']['t_wear'])} s",
        "",
        "VIDA GLOBAL",
        f"N_lim              = {_v776_fmt(result['life']['N_lim'])} ciclos",
        f"t_lim              = {_v776_fmt(result['life']['t_lim'])} s",
        f"mecanismo limite   = {result['life']['limiting_mechanism']}",
        "",
        "CLASSIFICACAO DAS ENTRADAS",
        "PRIMARIAS: F, R1, R2, E1, E2, nu1, nu2",
        "DERIVADAS: E', R', a, b, area, p_media, p_H, campo de tensoes",
        "ESTIMADAS: ue, v_slide, eta0, alpha_p, Rq1, Rq2, T0, mu, beta, parametros termicos, desgaste e limites",
        "",
        "OBSERVACAO",
        "A vida por fadiga depende explicitamente de sigma_ref, N_ref e m fornecidos no caso. A vida por desgaste depende da lei de Archard, do limite de desgaste, da velocidade de deslizamento e da frequencia de ciclos. Esses parametros permanecem entradas de referencia e nao sao derivados do contato Hertziano.",
        "",
        "=== FIM V7.81 ===",
    ]
    return "\n".join(lines)


def run_v776_reference_case():
    result = execute_v776_reference_case()
    report = format_v776_reference_report(result)
    return result, report


if __name__ == "__main__":
    # Execução normal no Pydroid 3: interface Tkinter.
    # O caso de referência continua disponível por run_v776_reference_case().
    root = tk.Tk()
    try:
        root.option_add("*Font", ("TkDefaultFont", 10))
    except Exception:
        pass
    App(root)
    root.mainloop()

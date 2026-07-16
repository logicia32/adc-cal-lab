"""センサ -> 単電源アンプ -> ADC -> マイコン の鎖で、読み値をどこまで較正できるかを見る。

前段（アンプ）は ngspice で解く（spice/afe.cir）。オフセット・ゲイン・低電圧/高電圧の
振れ切りを持つ、単電源の非反転アンプ。そこに ADC の量子化を足し、
「精度の基準電圧を"正"として突き合わせる」較正を、はしごの各段で当てて残差を見る。

- 個体ばらつき: オペアンプ Vos と帰還抵抗を振って複数基板を作る
- 較正のはしご: (1)オフセット (2)+ゲイン (3)+多項式/折れ線
- 温度: 基準電圧の tempco（良い品 5ppm/℃ vs 安い品 50ppm/℃）と Vos ドリフトを 0〜55℃で振る

評価は使用レンジ USABLE 内で行う。単電源アンプが 0V 近傍で振れ切る最下部は、そもそも
較正で救えないので範囲外に置く（実務でも仕様範囲を出力が飽和する手前までにする）。

ngspice が要る（Debian/Ubuntu: sudo apt install ngspice）。NGSPICE 環境変数でパス上書き可。
DC 特性値・tempco は代表値。値そのものでなく傾向を見るためのもの。
"""

import os
import re
import shutil
import subprocess
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
AFE = os.path.join(HERE, "spice", "afe.cir")
NGSPICE = os.environ.get("NGSPICE", "ngspice")

VREF_NOM = 5.0      # ADC 基準電圧の公称 [V]（マイコンはこれで換算すると思っている）
N_BITS = 14         # ADC 分解能。1mV/5V=200ppm を刻むには 12bit では足りず 14bit 級が要る
G_NOM = 2.0         # 前段の公称ゲイン（Rf=Rg=10k -> 1+Rf/Rg=2）
CAL_T = 25.0        # 較正した温度 [℃]
USABLE = (0.10, 2.30)   # 評価する入力レンジ [V]（アンプが振れ切る最下部・最上部は外す）
CAL_PTS = [0.10, 0.4, 0.8, 1.3, 1.8, 2.30]   # 較正の基準点（精度電圧源を当てる点。端は使用レンジ端）
INL_LSB = 8             # ADC の積分非直線性(INL)の目安 [LSB]。全域に弓なりで散る分布誤差

FULL = (1 << N_BITS) - 1


def _ng_ok():
    return shutil.which(NGSPICE) is not None or os.path.exists(NGSPICE)


def afe_transfer(vos=2e-3, rf=10e3, rg=10e3):
    """spice/afe.cir を Vos・帰還抵抗を差し替えて実行し、(v_true, v_afe) を返す。"""
    src = open(AFE).read()
    src = re.sub(r"VOS=\S+", f"VOS={vos}", src, count=1)
    src = re.sub(r"RF=\S+", f"RF={rf}", src, count=1)
    src = re.sub(r"RG=\S+", f"RG={rg}", src, count=1)
    tmp = os.path.join(HERE, "_afe_run.cir")
    open(tmp, "w").write(src)
    subprocess.run([NGSPICE, "-b", tmp], cwd=HERE, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    d = np.loadtxt(os.path.join(HERE, "afe.data"))
    return d[:, 0], d[:, 1]


def adc_read(v_afe, vref_actual=VREF_NOM, inl_lsb=INL_LSB):
    """アンプ出力を ADC で量子化し、マイコンが公称 VREF で戻した電圧を返す。
    vref_actual を公称からずらすと基準電圧ドリフト（ゲイン誤差）を、
    inl_lsb で ADC の INL（全域に散る弓なりの非直線性）を表せる。"""
    code_i = v_afe / vref_actual * FULL
    inl = inl_lsb * np.sin(np.pi * np.clip(code_i, 0, FULL) / FULL)   # 中央で最大の弓なり
    code = np.clip(np.round(code_i + inl), 0, FULL)
    return code / FULL * VREF_NOM     # マイコンは公称 VREF_NOM で戻す


def _nearest(v_true, targets):
    return [int(np.argmin(np.abs(v_true - t))) for t in targets]


def calibrate(v_true, v_adc, level):
    """v_adc（ADC 換算電圧）から真の入力を推定する較正を当て、推定値を返す。
    level: 'raw' | 'offset' | 'gain' | 'poly2' | 'poly3' | 'piecewise'."""
    if level == "raw":
        return v_adc / G_NOM                       # 理想ゲインと決め打ち、オフセットも無視
    if level == "offset":
        i = _nearest(v_true, [CAL_PTS[0]])[0]      # 低い基準点で、入力換算のズレを引く
        off = v_adc[i] / G_NOM - v_true[i]
        return v_adc / G_NOM - off
    pts = _nearest(v_true, CAL_PTS)
    x, y = v_adc[pts], v_true[pts]
    if level == "gain":                            # 2 点で 1 次（オフセット＋ゲイン）
        return np.polyval(np.polyfit(x[[0, -1]], y[[0, -1]], 1), v_adc)
    if level == "poly2":
        return np.polyval(np.polyfit(x, y, 2), v_adc)
    if level == "poly3":
        return np.polyval(np.polyfit(x, y, 3), v_adc)
    if level == "piecewise":
        return np.interp(v_adc, x, y)
    raise ValueError(level)


def err_mV(v_true, v_est):
    return (v_est - v_true) * 1e3


def _usable(v_true):
    return (v_true >= USABLE[0]) & (v_true <= USABLE[1])


def _worst(v_true, v_est):
    m = _usable(v_true)
    return np.max(np.abs(err_mV(v_true, v_est)[m]))


# ---------- 図 ----------

def fig_transfer():
    """図1: 前段の伝達特性。全域と、低電圧が曲がる様子（出力が 0V に届かない）。"""
    v, o = afe_transfer()
    ideal = G_NOM * v
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.4))
    ax1.plot(v, o, lw=2, label="AFE output (SPICE)")
    ax1.plot(v, ideal, "--", color="gray", lw=1.3, label="ideal  Vout = 2*Vin")
    ax1.axvspan(0, USABLE[0], color="0.85"); ax1.axvspan(USABLE[1], v.max(), color="0.85")
    ax1.set_xlabel("true input Vin [V]"); ax1.set_ylabel("AFE output [V]")
    ax1.set_title("front-end transfer (single-supply amp)")
    ax1.legend(fontsize=8); ax1.grid(alpha=.3)

    lo = v <= 0.35
    ax2.plot(v[lo], o[lo] * 1e3, lw=2, label="AFE output (SPICE)")
    ax2.plot(v[lo], ideal[lo] * 1e3, "--", color="gray", lw=1.3, label="ideal line")
    ax2.axhline(50, color="#c0392b", ls=":", lw=1)
    ax2.text(0.12, 62, "output floors ~50 mV\n(can't reach 0 V)", fontsize=8, color="#a03020")
    ax2.axvspan(0, USABLE[0], color="0.85")
    ax2.set_xlabel("true input Vin [V]"); ax2.set_ylabel("AFE output [mV]")
    ax2.set_title("low-end zoom: the curve bends away from the line")
    ax2.legend(fontsize=8, loc="lower right"); ax2.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "transfer.png"), dpi=130)
    print("saved transfer.png")


def _boards(n=6, seed=1):
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        vos = rng.normal(0, 2e-3)                  # Vos ばらつき ~ N(0, 2mV)
        rf = 10e3 * (1 + rng.normal(0, 0.01))      # 帰還抵抗 ±1%
        v, o = afe_transfer(vos=vos, rf=rf)
        out.append((v, o))
    return out


def fig_per_unit():
    """図2: 個体差（生の誤差）と、オフセット＋ゲイン較正後の残差。"""
    boards = _boards()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.4), sharex=True)
    for k, (v, o) in enumerate(boards):
        m = _usable(v)
        v_adc = adc_read(o)
        ax1.plot(v[m], err_mV(v, calibrate(v, v_adc, "raw"))[m], lw=1.2, alpha=.8)
        ax2.plot(v[m], err_mV(v, calibrate(v, v_adc, "gain"))[m], lw=1.2, alpha=.8,
                 label=f"board {k+1}")
    for ax, t in [(ax1, "raw reading (no calibration)"),
                  (ax2, "after per-board offset + gain")]:
        ax.axhline(0, color="gray", lw=.8); ax.set_xlabel("true input Vin [V]")
        ax.set_title(t); ax.grid(alpha=.3)
    ax1.set_ylabel("error [mV]")
    ax2.legend(fontsize=7, ncol=2)
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "per_unit.png"), dpi=130)
    print("saved per_unit.png")


def fig_ladder():
    """図3: 較正のはしご。オフセット＋ゲイン＋曲がりを持つ 1 枚で、段ごとの残差を見る。"""
    # オフセット(Vos)・ゲイン誤差(帰還抵抗+基準電圧)・曲がり を併せ持つ代表機
    v, o = afe_transfer(vos=3e-3, rf=10.2e3)
    v_adc = adc_read(o, vref_actual=VREF_NOM * 1.003)   # 基準 +0.3% でゲイン誤差
    m = _usable(v)
    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    for level, name in [("raw", "raw (assume ideal)"),
                        ("offset", "+ offset (low ref pt)"),
                        ("gain", "+ gain (2-point linear)"),
                        ("poly2", "+ 2nd-order polynomial")]:
        ax.plot(v[m], err_mV(v, calibrate(v, v_adc, level))[m], lw=1.8, label=name)
    ax.axhline(1, color="#888", ls=":"); ax.axhline(-1, color="#888", ls=":")
    ax.text(0.12, 1.3, "+/-1 mV target", color="#555", fontsize=8)
    ax.set_xlabel("true input Vin [V]"); ax.set_ylabel("residual error [mV]")
    ax.set_title("calibration ladder: each rung buys accuracy")
    ax.legend(fontsize=9); ax.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "ladder.png"), dpi=130)
    for level in ["raw", "offset", "gain", "poly2", "poly3", "piecewise"]:
        print(f"  {level:9s}: worst|err| over range = {_worst(v, calibrate(v, v_adc, level)):6.2f} mV")
    print("saved ladder.png")


def fig_temperature():
    """図4: 0〜55℃で、良い基準 vs 安い基準、温度補正あり/なし。"""
    temps = np.linspace(0, 55, 12)
    v_nom, o25 = afe_transfer(vos=2e-3)
    v_adc25 = adc_read(o25)
    pts = _nearest(v_nom, CAL_PTS)
    coeff25 = np.polyfit(v_adc25[pts], v_nom[pts], 2)   # 25℃で 2 次較正
    m = _usable(v_nom)

    def curve(tc_ref_ppm, comp=False):
        errs = []
        for T in temps:
            dT = T - CAL_T
            vos_T = 2e-3 + 2e-6 * dT                        # Vos ドリフト ~2uV/℃
            # 基準電圧 tempco（1 次 + わずかな 2 次で、温度補正が完全には効かない現実味）
            vref_T = VREF_NOM * (1 + tc_ref_ppm * 1e-6 * dT + 0.3e-6 * (dT / 30) ** 2 * tc_ref_ppm)
            _, oT = afe_transfer(vos=vos_T)
            est = np.polyval(coeff25, adc_read(oT, vref_actual=vref_T))
            errs.append(est - v_nom)
        errs = np.array(errs)                               # [温度, 入力]
        if comp:                                            # 温度センサで T の 1 次で戻す
            fit = np.polyfit(temps, errs, 1)
            errs = errs - (np.outer(temps, fit[0]) + fit[1])
        return np.max(np.abs(errs[:, m]), axis=1) * 1e3

    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    ax.plot(temps, curve(5), lw=2, label="good ref 5ppm/C, 25C cal")
    ax.plot(temps, curve(50), lw=2, label="cheap ref 50ppm/C, 25C cal")
    ax.plot(temps, curve(50, comp=True), lw=2, ls="--", label="cheap ref + temp compensation")
    ax.axhline(1, color="#888", ls=":"); ax.text(1, 1.15, "1 mV target", color="#555", fontsize=8)
    ax.set_xlabel("board temperature [C]"); ax.set_ylabel("worst error over range [mV]")
    ax.set_title("drift over 0-55C: reference choice vs temp compensation")
    ax.legend(fontsize=9); ax.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "temperature.png"), dpi=130)
    for tc, comp in [(5, False), (50, False), (50, True)]:
        print(f"  ref {tc:2d}ppm/C comp={comp}: worst over 0-55C = {np.max(curve(tc, comp)):.2f} mV")
    print("saved temperature.png")


def main():
    if not _ng_ok():
        sys.exit(f"ngspice not found ('{NGSPICE}'). apt install ngspice or set NGSPICE=")
    fig_transfer()
    fig_per_unit()
    fig_ladder()
    fig_temperature()


if __name__ == "__main__":
    main()

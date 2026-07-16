# adc-cal-lab

センサ → 単電源アンプ → ADC → マイコン の鎖で、読み値の誤差をどこまで較正できるかを見る小さな実験。
前段のアンプは ngspice で解き（`spice/afe.cir`）、そこに ADC の量子化と INL を足して、
「精度の基準電圧を"正"として突き合わせる」較正を、はしごの各段で当てて残差を見る。

## 何を見るか

- **前段の伝達（`transfer.png`）** … 単電源アンプは出力が 0V まで振れず、低電圧で伝達が曲がる。
- **個体差（`per_unit.png`）** … オペアンプ Vos と帰還抵抗を振った複数基板。生の誤差は基板ごとに
  ばらつくが、オフセット＋ゲインを個体ごとに較正すると、残るのは全基板に共通の弓なり（ADC の INL）。
- **較正のはしご（`ladder.png`）** … raw → オフセット → ＋ゲイン → ＋2次多項式。
  この例では、オフセット＋ゲイン（1次）で約 **0.95mV**、2次多項式で約 **0.12mV**。
- **温度（`temperature.png`）** … 0〜55℃。基準電圧 5ppm/℃ なら室温1点較正のまま **0.35mV** に収まるが、
  50ppm/℃ だと **3.5mV** に広がり、温度センサでの補正が要る（**0.11mV** に回復）。

結論の目安: 制御盤 0〜55℃ で 1mV クラスなら、**個体ごとのオフセット＋ゲイン較正＋良い基準の選定**で届く。
温度センサによる温度補正は、それより厳しい仕様か、安い基準をソフトで救うときの一段上。

## 使い方

```bash
pip install -r requirements.txt
# ngspice が要る（Debian/Ubuntu: sudo apt install ngspice）
python adc_cal.py        # transfer.png / per_unit.png / ladder.png / temperature.png
```

`spice/afe.cir` 単体でも実行できる（`ngspice -b spice/afe.cir` → `afe.data`）。
ngspice のパスは `NGSPICE=/path/to/ngspice python adc_cal.py` で上書き可。

## モデルについて

前段は level-1 的な代役ではなく、帰還ループを解いた素朴なオペアンプのマクロモデル（有限ゲイン・
入力オフセット・出力の振れ切り）。ADC の INL は代表的な弓なりを与えたもの。DC 特性値・tempco・INL は
どれも代表値で、値そのものでなく傾向を見るためのもの。実機の実測ではなくモデル。実際に設計するときは、
使う型番のデータシートを引くこと。

## ライセンス

MIT（`LICENSE`）。すべて自作のコードとネットリスト。

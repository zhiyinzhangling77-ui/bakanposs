"""**旗156 の後知恵（post-hoc）**：印字済みの 2×2 表から**対応の強さ（φ 係数）**だけを出す。

**★これは事前登録 `PREREGISTRATION_step156.md` に無い量である**（旗149 の「印字済みの `p` に
多重性の補正を当てた（後知恵と明記）」と同じ扱い）。**主判定 D-1 の結論には使わない。**
**新しい合成も実データも触らない**——入力は `logs/step154_20260914_200009.txt` に印字済みの
`both / b / c / neither` の 4 つだけである。

**なぜ出すか**：事前登録 2 節の検出力見積もりは、潜在相関 ρ をコピュラで動かして
「n=600 で 0.79（ρ=0）〜1.00（ρ=0.95）」と書いた。**実際の対応がどちら寄りだったかは、
表が出てからでないと分からない。** 出た表で φ を測れば、**本周の McNemar がどのくらいの
検出力で回っていたのかを、上ぶれ・下ぶれの向きだけでも言える。**
"""
import math

TABLES = {
    "D-1 主判定：全 600 本": dict(both=11, b=33, c=47, neither=509),
    "D-3 併記：先頭 200 本": dict(both=6, b=7, c=16, neither=171),
}


def report(tag: str, both: int, b: int, c: int, neither: int) -> None:
    n = both + b + c + neither
    kh, kle = both + b, both + c
    ph, ple = kh / n, kle / n
    phi = (both * neither - b * c) / math.sqrt(kh * (n - kh) * kle * (n - kle))
    exp_both = n * ph * ple
    exp_disc = n * (ph * (1 - ple) + ple * (1 - ph))
    print(f"\n── {tag}（n={n}）")
    print(f"   θ→γH {kh}/{n} = {ph:.3f}／θ→γLE {kle}/{n} = {ple:.3f}"
          f"／差 = {ple - ph:+.3f}")
    print(f"   両方有意 {both} 本（独立なら期待 {exp_both:.1f} 本）")
    print(f"   discordant {b + c} 本（独立なら期待 {exp_disc:.1f} 本）")
    print(f"   **φ 係数 = {phi:+.3f}**（0 が独立・1 が完全一致）")


if __name__ == "__main__":
    print("**旗156 の後知恵：印字済みの 2×2 表から対応の強さだけを出す**")
    print("**事前登録に無い量である（後知恵と明記）。主判定 D-1 の結論には使わない。**")
    for tag, t in TABLES.items():
        report(tag, **t)

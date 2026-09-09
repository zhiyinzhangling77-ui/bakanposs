"""**出力を最初からファイルに残す**（旗110 の反省）。

## なぜ要るのか

**旗108 は 146 サイトを走らせ、出力がターミナルの履歴を超えて遡れなくなった。**
**必要だった数字（CN-Du2 の 4 群の日数）が読めず、走らせ直すことになった。**

**「長くなりそうなときは `tee` してください」と頼むのは設計の押し付けである。**
**道具の側が、最初からファイルに残す。**

## 何をするか

  ・`research/logs/<名前>_<日時>.txt` を開き、**標準出力と標準エラーの両方**を
    **画面とファイルの両方**へ流す（`tee` と同じ）。
  ・**標準エラーも入れる**——**pandas の警告や `RuntimeWarning` は stderr に出る**ので、
    **画面では出力に混ざるのに、`> file` では落ちる**。**混ざったまま残す方が再現に役立つ。**
  ・**最初と最後にパスを印字する**（**最後だけだと、途中で止めたときに分からない**）。
  ・**失敗しても走行を止めない**——**記録が取れないことは、解析を止める理由にならない。**

    from runlog import tee_stdout
    tee_stdout("step108")        # main() の parse_args の直後に置く
"""
from __future__ import annotations

import atexit
import sys
from datetime import datetime
from pathlib import Path

LOG_DIR = Path(__file__).resolve().parent / "logs"


class _Tee:
    """画面とファイルの両方へ書く。**片方が失敗しても、もう片方は生かす。**"""

    def __init__(self, stream, fh):
        self._stream, self._fh = stream, fh

    def write(self, s):
        try:
            self._stream.write(s)
        except Exception:
            pass
        try:
            self._fh.write(s)
        except Exception:
            pass
        return len(s)

    def flush(self):
        for x in (self._stream, self._fh):
            try:
                x.flush()
            except Exception:
                pass

    def isatty(self):
        return getattr(self._stream, "isatty", lambda: False)()


def tee_stdout(tag: str, quiet: bool = False) -> Path | None:
    """**標準出力と標準エラーを、画面とファイルの両方へ流す。**

    戻り値はログのパス（**開けなければ None**）。**開けなくても走行は止めない。**

    **道具の欠陥 #81（旗148 で踏んだ）**：名前は秒までしか無いので、
    **同じ秒に 2 つ走ると 2 つのプロセスが同じファイルを `"w"` で開き、互いの行を潰し合う。**
    旗148 では `--gatesB` と `--gatesD` の出力が 1 つのログに混ざり、
    **どちらのログも証拠として使えなくなった**（判定そのものは無事）。
    **`"x"`（排他作成）にして、埋まっていたら `_2`・`_3` と後ろに足す。**
    """
    try:
        LOG_DIR.mkdir(parents=True, exist_ok=True)
        stem = f"{tag}_{datetime.now():%Y%m%d_%H%M%S}"
        fh, path, i = None, None, 1
        while fh is None:
            path = LOG_DIR / (f"{stem}.txt" if i == 1 else f"{stem}_{i}.txt")
            try:                                              # 行バッファ＝途中で止めても残る
                fh = open(path, "x", encoding="utf-8", buffering=1)
            except FileExistsError:
                i += 1
                if i > 99:                                    # ここまで来たら名前の付け方が壊れている
                    raise
    except Exception as e:                                    # 書けない環境でも解析は続ける
        print(f"  （**記録を残せない**：{type(e).__name__}: {str(e)[:60]}）")
        return None

    sys.stdout = _Tee(sys.stdout, fh)
    sys.stderr = _Tee(sys.stderr, fh)                         # **警告も残す**

    if not quiet:
        print(f"  【記録】この実行の出力は **{path}** にも残ります"
              f"（**ターミナルを遡らなくてよい**）。")

    def _close():
        try:
            print(f"\n  【記録】出力を **{path}** に保存しました。")
        except Exception:
            pass
        try:
            fh.flush(); fh.close()
        except Exception:
            pass

    atexit.register(_close)
    return path

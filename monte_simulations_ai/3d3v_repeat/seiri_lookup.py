from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_TARGET_DIR = SCRIPT_DIR / "3d3v_lookup"

SCORING_RE = re.compile(
    r"^\s*(CL|TL|TLLookup):\s+"
    r"([+-]?\d+\.\d+E[+-]\d+)\s+"
    r"([+-]?\d+\.\d+E[+-]\d+)\s+"
    r"([+-]?\d+\.\d+E[+-]\d+)\s+"
    r"([+-]?\d+\.\d+E[+-]\d+)\s*$"
)
RUN_ID_RE = re.compile(r"run_(\d+)")
RESULT_SIZE_RE = re.compile(r"result_(\d+)")

OUTPUT_HEADER = ["Run_ID", "CL_CX", "CL_EL", "TL_CX", "TL_EL", "TLLu_CX", "TLLu_EL"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "3d3v_lookup の run.log から、既存 ensemble_results_*.csv と互換の"
            "CSV を作る。TLLu_CX と TLLu_EL を追加出力する。"
        )
    )
    parser.add_argument(
        "--target-dir",
        type=Path,
        default=DEFAULT_TARGET_DIR,
        help="result_* を含むディレクトリ。デフォルトは 3d3v_lookup。",
    )
    return parser.parse_args()


def get_run_id(path: Path) -> int:
    match = RUN_ID_RE.search(str(path))
    if match is None:
        raise ValueError(f"run ID を解釈できません: {path}")
    return int(match.group(1))


def get_result_size(path: Path) -> int:
    match = RESULT_SIZE_RE.search(str(path))
    if match is None:
        raise ValueError(f"result サイズを解釈できません: {path}")
    return int(match.group(1))


def parse_run_log(log_path: Path) -> dict[str, float | int | None]:
    row: dict[str, float | int | None] = {
        "Run_ID": get_run_id(log_path),
        "CL_CX": None,
        "CL_EL": None,
        "TL_CX": None,
        "TL_EL": None,
        "TLLu_CX": None,
        "TLLu_EL": None,
    }

    with log_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            match = SCORING_RE.match(line)
            if match is None:
                continue

            label = match.group(1)
            cx_value = float(match.group(3))
            el_value = float(match.group(4))

            if label == "CL":
                row["CL_CX"] = cx_value
                row["CL_EL"] = el_value
            elif label == "TL":
                row["TL_CX"] = cx_value
                row["TL_EL"] = el_value
            elif label == "TLLookup":
                row["TLLu_CX"] = cx_value
                row["TLLu_EL"] = el_value

    return row


def collect_result_dirs(target_dir: Path) -> list[Path]:
    return sorted(
        [path for path in target_dir.iterdir() if path.is_dir() and RESULT_SIZE_RE.fullmatch(path.name)],
        key=lambda path: get_result_size(path),
    )


def write_result_csv(result_dir: Path) -> tuple[Path, int]:
    log_files = sorted(result_dir.glob("run_*/run.log"), key=get_run_id)
    if not log_files:
        raise FileNotFoundError(f"run.log が見つかりません: {result_dir}")

    rows = [parse_run_log(log_path) for log_path in log_files]
    output_path = result_dir.parent / f"ensemble_results_{get_result_size(result_dir)}.csv"

    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_HEADER)
        writer.writeheader()
        writer.writerows(rows)

    return output_path, len(rows)


def main() -> None:
    args = parse_args()
    target_dir = args.target_dir.resolve()
    if not target_dir.exists():
        raise SystemExit(f"ディレクトリが存在しません: {target_dir}")

    result_dirs = collect_result_dirs(target_dir)
    if not result_dirs:
        raise SystemExit(f"result_* ディレクトリが見つかりません: {target_dir}")

    for result_dir in result_dirs:
        output_path, count = write_result_csv(result_dir)
        print(f"{output_path} <- {count} runs")


if __name__ == "__main__":
    main()

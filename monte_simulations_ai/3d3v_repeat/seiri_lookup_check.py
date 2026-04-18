from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_TARGET_DIR = SCRIPT_DIR / "3d3v_lookup_check"

SCIENTIFIC_RE = r"[+-]?\d+(?:\.\d+)?E[+-]?\d+"
SCORE_ROW_RE = re.compile(
    rf"^\s*(CL|TR):\s+({SCIENTIFIC_RE})\s+({SCIENTIFIC_RE})\s+({SCIENTIFIC_RE})\s+({SCIENTIFIC_RE})\s*$"
)
EL_ONLY_RE = re.compile(
    rf"^\s*(A\(EL\)|CLInner\(EL\)|TRInner\(EL\)|TRPretab\(EL\)):\s+({SCIENTIFIC_RE})\s*$"
)
RUN_ID_RE = re.compile(r"run_(\d+)")
RESULT_SIZE_RE = re.compile(r"result_fixedcdf_new_(\d+)")

OUTPUT_HEADER = [
    "Run_ID",
    "A_EL",
    "CL_CX",
    "CL_EL",
    "CLInner_EL",
    "TR_CX",
    "TR_EL",
    "TRInner_EL",
    "TRPretab_EL",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "3d3v_lookup_check の run.log から、ensemble_results_*.csv を作る。"
        )
    )
    parser.add_argument(
        "--target-dir",
        type=Path,
        default=DEFAULT_TARGET_DIR,
        help="result_* を含むディレクトリ。デフォルトは 3d3v_lookup_check。",
    )
    return parser.parse_args()


def parse_float(value: str) -> float:
    return float(value.replace("D", "E"))


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
    row: dict[str, float | int | None] = {column: None for column in OUTPUT_HEADER}
    row["Run_ID"] = get_run_id(log_path)

    with log_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            score_match = SCORE_ROW_RE.match(line)
            if score_match is not None:
                label = score_match.group(1)
                cx_value = parse_float(score_match.group(3))
                el_value = parse_float(score_match.group(4))

                if label == "CL":
                    row["CL_CX"] = cx_value
                    row["CL_EL"] = el_value
                elif label == "TR":
                    row["TR_CX"] = cx_value
                    row["TR_EL"] = el_value
                continue

            el_only_match = EL_ONLY_RE.match(line)
            if el_only_match is None:
                continue

            label = el_only_match.group(1)
            value = parse_float(el_only_match.group(2))

            if label == "A(EL)":
                row["A_EL"] = value
            elif label == "CLInner(EL)":
                row["CLInner_EL"] = value
            elif label == "TRInner(EL)":
                row["TRInner_EL"] = value
            elif label == "TRPretab(EL)":
                row["TRPretab_EL"] = value

    return row


def collect_result_dirs(target_dir: Path) -> list[Path]:
    return sorted(
        [path for path in target_dir.iterdir() if path.is_dir() and RESULT_SIZE_RE.fullmatch(path.name)],
        key=get_result_size,
    )


def write_result_csv(result_dir: Path) -> tuple[Path, int]:
    log_files = sorted(result_dir.glob("run_*/run.log"), key=get_run_id)
    if not log_files:
        raise FileNotFoundError(f"run.log が見つかりません: {result_dir}")

    rows = [parse_run_log(log_path) for log_path in log_files]
    output_path = result_dir.parent / f"ensemble_results_fixedcdf_new_{get_result_size(result_dir)}.csv"

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

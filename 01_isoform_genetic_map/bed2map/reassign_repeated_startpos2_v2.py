import pandas as pd
import argparse
import sys

MAX_PASSES = 20


def auto_rmrepeat_pos(df):
    """
    Ensure every isoform has a unique position in column [1] (TSS).

    Resolution priority for each group of colliding TSS positions:
      1. Replace with TTS (column [2]) if all TTS values in the group are unique.
      2. Otherwise apply the smallest integer offset that yields a globally
         unique position.

    Parameters
    ----------
    df : pd.DataFrame
        BED-like dataframe (no header). Column [1] = TSS, [2] = TTS, [3] = isoform ID.

    Returns
    -------
    pd.DataFrame
        Updated dataframe with collisions resolved where possible.
    """
    # Build a set of all current positions for O(1) membership checks.
    # This is updated as we assign new positions so later groups see
    # already-committed values.
    existing_positions = set(df[1].tolist())

    collision_groups = df[1].value_counts()

    for tss_value, count in collision_groups.items():
        if count < 2:
            continue  # no collision

        group = df[df[1] == tss_value]

        # ── Tier 1: swap to TTS if all TTS values are distinct ──────────────
        tts_values = group[2].tolist()
        if len(set(tts_values)) == len(tts_values):
            # Check none of the TTS values clash with positions outside the group
            outside = existing_positions - set(group[1].tolist())
            if not any(t in outside for t in tts_values):
                existing_positions -= set(group[1].tolist())
                df.loc[group.index, 1] = tts_values
                existing_positions.update(tts_values)
                continue

        # ── Tier 2: add smallest valid integer offset ────────────────────────
        # Remove this group's current positions so we don't collide with ourselves.
        existing_positions -= set(group[1].tolist())

        new_positions = []
        for rank, idx in enumerate(group.index):
            base = int(df.loc[idx, 1])
            offset = 0
            while True:
                # Spread candidates: rank shifts each isoform slightly apart
                candidate = base + rank + offset
                if candidate not in existing_positions and candidate not in new_positions:
                    new_positions.append(candidate)
                    break
                offset += 1

        df.loc[group.index, 1] = new_positions
        existing_positions.update(new_positions)

    return df


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Deduplicate isoform TSS positions in a BED file so every isoform "
            "has a unique coordinate, as required for genetic map construction."
        )
    )
    parser.add_argument("bedfile", help="Path to the input BED file (tab-separated, no header)")
    parser.add_argument(
        "--max-passes", type=int, default=MAX_PASSES,
        help=f"Maximum deduplication passes (default: {MAX_PASSES})"
    )
    args = parser.parse_args()

    # ── Load ────────────────────────────────────────────────────────────────
    df = pd.read_csv(args.bedfile, sep="\t", index_col=None, header=None)
    df.index = df[3]
    df[1] = df[1].astype(int)

    total_isoforms = df.shape[0]
    print(f"Loaded {total_isoforms} isoforms from '{args.bedfile}'.")

    # ── Deduplicate ──────────────────────────────────────────────────────────
    passes = 0
    while passes < args.max_passes:
        n_unique = len(set(df[1].tolist()))
        if n_unique == total_isoforms:
            break  # all positions unique — done
        duplicates = total_isoforms - n_unique
        print(f"  Pass {passes + 1}: {duplicates} duplicate position(s) remaining.")
        df = auto_rmrepeat_pos(df)
        passes += 1

    # ── Final check ──────────────────────────────────────────────────────────
    n_unique_final = len(set(df[1].tolist()))
    remaining = total_isoforms - n_unique_final

    if remaining == 0:
        print(f"All positions unique after {passes} pass(es). Writing output.")
    else:
        print(
            f"WARNING: {remaining} duplicate position(s) still remain after "
            f"{passes} pass(es). Please review the output file manually.",
            file=sys.stderr
        )

    # ── Sort and write ───────────────────────────────────────────────────────
    df = df.sort_values(by=[0, 1])
    out_path = "out_" + args.bedfile
    df.to_csv(out_path, header=None, sep="\t", index=None)
    print(f"Output written to '{out_path}'.")


if __name__ == "__main__":
    main()

import sys
import pandas as pd


def main():
    fn = sys.argv[1]
    df = pd.read_csv(fn)

    df["One"] = 1
    ndf = pd.pivot_table(
        df,
        values=["One"],
        index=["n", "graph", "d", "class"],
        aggfunc="count",
    ).reset_index()

    pivoted_df = ndf.pivot(index=["n", "graph", "d"], columns="class", values="One")
    pivoted_df.fillna(0, inplace=True)
    pivoted_df["Sum"] = pivoted_df[
        ["BAD", "BAD_CLIP", "FULL_CLIP", "GOOD", "HALF_BAD", "HALF_CLIP"]
    ].sum(axis=1)

    for c in ["BAD", "BAD_CLIP", "FULL_CLIP", "GOOD", "HALF_BAD", "HALF_CLIP"]:
        new_c = c + " (%)"
        pivoted_df[new_c] = pivoted_df[c] / pivoted_df["Sum"]
        pivoted_df[new_c] = pivoted_df[new_c].round(3)

    print(pivoted_df, file=sys.stderr)
    # pivoted_df.to_csv(sys.stdout)


if __name__ == "__main__":
    main()

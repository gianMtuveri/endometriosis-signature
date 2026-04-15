from pathlib import Path

from endometriosis_signature.dataset import (
    get_feature_columns,
    load_expression_matrix,
    make_clean_subset,
    merge_expression_and_metadata,
    parse_family_xml,
    save_dataset,
    save_feature_columns,
)


def main() -> None:
    raw_dir = Path("data/raw")
    processed_dir = Path("data/processed")
    processed_dir.mkdir(parents=True, exist_ok=True)

    expression_df = load_expression_matrix(raw_dir / "GSE51981_series_matrix.txt")
    metadata_df = parse_family_xml(raw_dir / "GSE51981_family.xml")
    data = merge_expression_and_metadata(expression_df, metadata_df)

    clean_data = make_clean_subset(data)
    feature_cols = get_feature_columns(data)

    save_dataset(data, processed_dir / "endometriosis_full.parquet")
    save_dataset(clean_data, processed_dir / "endometriosis_clean_filtered.parquet")
    save_feature_columns(feature_cols, processed_dir / "feature_columns.json")

    print("Full dataset:", data.shape)
    print(data["target"].value_counts())
    print("\nClean dataset:", clean_data.shape)
    print(clean_data["target"].value_counts())


if __name__ == "__main__":
    main()
#!/usr/bin/env python3
#
# VocProfile - 02_describe_accession_means.py
#
# Descriptive statistics of the accession means per compound, main class and subclass,
# for all accessions and by K group, on the data as they are and after the accession
# means without detection are filled (Table 1, Table S5).
#
# in   MEANS_FILE, DATASETS
# out  results/02_describe_accession_means/statistics_<dataset>/
#      results/02_describe_accession_means/statistics_by_K_<dataset>/
#      results/02_describe_accession_means/02_qc.txt
#
# Usage:
#   python scripts/02_describe_accession_means.py <config.txt>
#
# Notes: docs/NOTES.md, step 02
#
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gmean

from lib_config import load_config, input_path, require_sheets, output_dir
from lib_run import start_step, finish_step, save_figure, write_table

cfg, step, args = load_config("02")
start_step(cfg, 'means')
LEVELS = cfg["significance_levels"]
TREATMENT = (
    'first report: values without detection are left out of every statistic; nothing is filled',
    'second report: accession means without detection: set to NONDETECT_FRACTION of the smallest positive accession mean of the compound',
    'statistics of a class or subclass pool all accession x compound values of the class',
)
DPI = step["dpi"]


def compute_descriptive_stats(data: pd.DataFrame) -> pd.DataFrame:
    """Descriptive statistics of each numeric column."""
    numeric_data = data.select_dtypes(include=[np.number])
    if numeric_data.empty or numeric_data.size == 0:
        columns = ["Count", "Missing", "Mean", "SE_Mean", "Median", "StDev", "Variance",
                   "CoefVar", "Min", "Q1", "Q3", "Max", "Range", "IQR", "MAD",
                   "Skewness", "Kurtosis", "GeometricMean"]

        stats_df = pd.DataFrame([np.nan] * len(columns), index=columns).T
        stats_df['Mean ± SD'] = 'nan ± nan'
        return stats_df

    mean_val = numeric_data.mean(axis=0)
    std_val = numeric_data.std(axis=0)
    coef_var = (std_val / mean_val).replace([np.inf, -np.inf], np.nan)

    median_val = numeric_data.median(axis=0)
    mad_val = np.abs(numeric_data.subtract(median_val, axis=1)).median(axis=0)

    q1 = numeric_data.quantile(0.25, axis=0)
    q3 = numeric_data.quantile(0.75, axis=0)
    iqr = q3 - q1

    epsilon = 1e-9
    geometric_mean = gmean(numeric_data.abs().fillna(epsilon) + epsilon)

    stats_df = pd.DataFrame({
        "Count": numeric_data.notnull().sum(),
        "Missing": numeric_data.isnull().sum(),
        "Mean": mean_val,
        "SE_Mean": numeric_data.sem(axis=0),
        "Median": median_val,
        "StDev": std_val,
        "Variance": numeric_data.var(axis=0),
        "CoefVar": coef_var,
        "Min": numeric_data.min(axis=0),
        "Q1": q1,
        "Q3": q3,
        "Max": numeric_data.max(axis=0),
        "Range": numeric_data.max(axis=0) - numeric_data.min(axis=0),
        "IQR": iqr,
        "MAD": mad_val,
        "Skewness": numeric_data.skew(axis=0),
        "Kurtosis": numeric_data.kurt(axis=0),
        "GeometricMean": geometric_mean
    })

    stats_df['Mean ± SD'] = stats_df.apply(
        lambda row: f"{row['Mean']:.3f} ± {row['StDev']:.3f}", axis=1
    )

    return stats_df


def generate_unified_report(dataframe: pd.DataFrame) -> pd.DataFrame:
    """Statistics per compound, main class and subclass in one table."""
    print("Computing statistics for Compounds...")
    compound_stats = compute_descriptive_stats(dataframe)
    compound_stats['Level'] = 'Compound'
    compound_stats = compound_stats.reset_index(names=['Main_Class', 'Subclass', 'Name'])

    print("Computing statistics for Main Classes...")
    class_stats_list = []
    unique_classes = dataframe.columns.get_level_values('Main_Class').unique()
    for class_name in unique_classes:
        class_data = dataframe.columns[dataframe.columns.get_level_values('Main_Class') == class_name]
        if not class_data.empty:
            all_values = dataframe.loc[:, class_data].values.flatten()
            stats = compute_descriptive_stats(pd.DataFrame(all_values, columns=['Value']))

            class_row = stats.loc['Value'].to_frame().T
            class_row['Level'] = 'Main_Class'
            class_row['Main_Class'] = class_name
            class_row['Subclass'] = 'N/A'
            class_row['Name'] = class_name
            class_stats_list.append(class_row)

    class_stats = pd.concat(class_stats_list, ignore_index=True)

    print("Computing statistics for Subclasses...")
    subclass_stats_list = []
    unique_subclasses = dataframe.columns.get_level_values('Subclass').unique()
    for subclass_name in unique_subclasses:
        subclass_data = dataframe.columns[dataframe.columns.get_level_values('Subclass') == subclass_name]
        if not subclass_data.empty:
            all_values = dataframe.loc[:, subclass_data].values.flatten()
            stats = compute_descriptive_stats(pd.DataFrame(all_values, columns=['Value']))

            subclass_row = stats.loc['Value'].to_frame().T
            subclass_row['Level'] = 'Subclass'
            subclass_row['Subclass'] = subclass_name
            main_class_name = dataframe.columns.get_level_values('Main_Class')[dataframe.columns.get_level_values('Subclass') == subclass_name][0]
            subclass_row['Main_Class'] = main_class_name
            subclass_row['Name'] = subclass_name
            subclass_stats_list.append(subclass_row)

    subclass_stats = pd.concat(subclass_stats_list, ignore_index=True)

    final_columns = ['Level', 'Main_Class', 'Subclass', 'Name', 'Count', 'Missing', 'Mean', 'StDev', 'Mean ± SD', 'SE_Mean',
                     'Median', 'Variance', 'CoefVar', 'Min', 'Q1', 'Q3', 'Max',
                     'Range', 'IQR', 'MAD', 'Skewness', 'Kurtosis', 'GeometricMean']

    unified_report = pd.concat([compound_stats, class_stats, subclass_stats], ignore_index=True)
    unified_report = unified_report.reindex(columns=final_columns)

    unified_report.sort_values(by=['Level', 'Main_Class', 'Subclass', 'Name'], inplace=True)

    return unified_report


def process_data(dataframe: pd.DataFrame) -> (pd.DataFrame, pd.DataFrame):
    """Set non-detects to half the smallest positive value of the compound; also return the natural log."""
    print("Setting non-detects to half the smallest positive value...")

    df_filled = dataframe.copy()

    for col in df_filled.columns:
        if df_filled[col].dtype.kind not in 'biufc':
            continue

        min_val = df_filled.loc[df_filled[col] > 0, col].min(skipna=True)

        if pd.notna(min_val):
            fill_value = min_val * cfg["settings"].value("NONDETECT_FRACTION")
        else:
            fill_value = 1e-9

        df_filled[col] = df_filled[col].fillna(fill_value)

    df_log_transformed = df_filled.apply(np.log)

    df_binary_info = pd.DataFrame(columns=['Compound Name', 'Missing Percentage (%)'])

    return df_filled, df_log_transformed, df_binary_info


def plot_missing_percentage(dataframe: pd.DataFrame, output_folder: str, group_name: str = ""):
    """Bar plot of the percentage of accessions in which each compound was not detected."""
    print("Calculating and plotting missing data percentage...")
    total_samples = dataframe.shape[0]
    compounds = dataframe.columns.get_level_values('Acc')
    missing_percentage = (dataframe.isnull().sum() / total_samples) * 100

    plt.figure(figsize=(14, 6))
    missing_percentage.plot(kind="bar", color='skyblue')

    plt.xlabel("Compound", fontsize=12)
    plt.ylabel("Percentage of Missing Samples (%)", fontsize=12)
    plt.title(f"Missing Data Percentage per Compound ({group_name})", fontsize=14)
    plt.xticks(range(len(compounds)), compounds, rotation=90, fontsize=8)
    plt.tight_layout()

    plot_path = os.path.join(output_folder, f"missing_data_percentage_{group_name}.png")
    save_figure(plot_path, dpi=DPI)
    plt.close()
    print(f"Missing data percentage plot saved to: {plot_path}")


def generate_missing_data_report(dataframe: pd.DataFrame) -> pd.DataFrame:
    """Percentage of missing and present values per compound, sorted."""
    print("Generating missing data percentage report...")
    total_samples = dataframe.shape[0]
    missing_count = dataframe.isnull().sum()
    present_count = dataframe.notna().sum()

    report_df = pd.DataFrame({
        "Compound": missing_count.index.get_level_values('Acc'),
        "Missing Percentage (%)": (missing_count / total_samples) * 100,
        "Present Percentage (%)": (present_count / total_samples) * 100
    }).reset_index(drop=True)

    report_df = report_df.sort_values(by="Missing Percentage (%)", ascending=False)

    return report_df


if __name__ == "__main__":
    file_path = input_path(cfg, step["input"])
    sheets = step["sheets"]
    require_sheets(file_path, [cfg["sheet_of"][name] for name in sheets])
    out_root = output_dir(cfg, "02_describe_accession_means")
    folder_all = step["output_folder"]
    folder_k = step["output_folder"] + "_by_K"

    for sheet_name in sheets:
        print(f"\nProcessing dataset: {sheet_name}")

        try:
            df_full = pd.read_excel(file_path, sheet_name=cfg["sheet_of"][sheet_name], header=cfg["header_rows"], index_col=[0, 1])
            df_full.index.names = ['Acc', 'K']
            df_full.columns.names = ['Main_Class', 'Subclass', 'Acc']
            print(f"Sheet '{sheet_name}' loaded successfully.")
        except Exception as e:
            print(f"Error loading sheet '{sheet_name}': {e}")
            continue

        df_compounds = df_full.reset_index(level=[0, 1], drop=True)

        output_folder_all = os.path.join(out_root, f"{folder_all}_{sheet_name}")
        os.makedirs(output_folder_all, exist_ok=True)

        print(f"Analysis of the data as they are ({sheet_name})")
        unified_report_original = generate_unified_report(df_compounds)
        missing_report = generate_missing_data_report(df_compounds)

        output_path_xlsx_original = os.path.join(output_folder_all, f"unified_statistical_report_original_{sheet_name}.xlsx")
        with pd.ExcelWriter(output_path_xlsx_original) as writer:
            write_table(unified_report_original, writer, sheet_name='Original Data Analysis', index=False)
            write_table(missing_report, writer, sheet_name='Missing Data Report', index=False)
        print(f"Saved: {output_path_xlsx_original}")
        plot_missing_percentage(df_compounds, output_folder_all, group_name=f"{sheet_name}_all")

        df_filled, df_log_transformed, binary_info = process_data(df_compounds)

        output_path_filled = os.path.join(output_folder_all, f"processed_and_filled_data_{sheet_name}.xlsx")
        write_table(df_filled, output_path_filled, data_names=True)
        print(f"Saved: {output_path_filled}")

        print(f"Analysis after filling non-detects ({sheet_name})")
        unified_report_processed = generate_unified_report(df_filled)

        output_path_processed_xlsx = os.path.join(output_folder_all, f"unified_statistical_report_processed_{sheet_name}.xlsx")
        with pd.ExcelWriter(output_path_processed_xlsx) as writer:
            write_table(unified_report_processed, writer, sheet_name='Processed Data Analysis', index=False)
            if not binary_info.empty:
                write_table(binary_info, writer, sheet_name='Binary Converted Compounds', index=False)
        print(f"Saved: {output_path_processed_xlsx}")

        output_path_log = os.path.join(output_folder_all, f"log_transformed_data_{sheet_name}.xlsx")
        write_table(df_log_transformed, output_path_log, data_names=True)
        print(f"Saved: {output_path_log}")

        print(f"Analysis by K group ({sheet_name})")
        k_groups = df_full.groupby(level='K')
        all_k_reports = []

        for k_value, group_df in k_groups:
            print(f"  Analyzing K={k_value}...")
            group_df_compounds = group_df.reset_index(level=[0, 1], drop=True)
            unified_report_k = generate_unified_report(group_df_compounds)
            unified_report_k['K_Value'] = k_value
            all_k_reports.append(unified_report_k)

            output_folder_k = os.path.join(out_root, f"{folder_k}_{sheet_name}", str(k_value))
            os.makedirs(output_folder_k, exist_ok=True)
            missing_report_k = generate_missing_data_report(group_df_compounds)

            output_path_k_xlsx = os.path.join(output_folder_k, f"statistical_report_K_{k_value}_{sheet_name}.xlsx")
            with pd.ExcelWriter(output_path_k_xlsx) as writer:
                write_table(unified_report_k, writer, sheet_name=f'Analysis_K_{k_value}', index=False)
                write_table(missing_report_k, writer, sheet_name=f'Missing_Data_K_{k_value}', index=False)
            plot_missing_percentage(group_df_compounds, output_folder_k, group_name=f"K_{k_value}_{sheet_name}")

        final_k_report = pd.concat(all_k_reports, ignore_index=True)
        final_k_report.sort_values(by=['Main_Class', 'Subclass', 'Name', 'K_Value'], inplace=True)
        final_k_report.set_index(['Main_Class', 'Subclass', 'Name', 'K_Value'], inplace=True)

        output_path_final_k = os.path.join(out_root, f"{folder_k}_{sheet_name}", f"consolidated_k_analysis_report_{sheet_name}.xlsx")
        write_table(final_k_report, output_path_final_k, sheet_name='K-Group Analysis')
        print(f"Saved: {output_path_final_k}")

    print("\nAll datasets processed.")

finish_step(cfg, TREATMENT)

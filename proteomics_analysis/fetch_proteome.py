__author__ = 'Liam TrinhNguyen'
__email__ = 'LiamTrinhNguyen@gmail.com'
__version__ = 'ProteinPipeline_v1.4'

import os
import sys
import logging
import webbrowser
from datetime import datetime
from pathlib import Path
import polars as pl
from Bio import Entrez, SeqIO
from tqdm import tqdm

class ProteinPipeline:
    def __init__(self):
        self.run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.sample_query = "Homo sapiens[Organism] AND Insulin Receptor[Protein Name]"
        self.max_results = 5
        
        self.script_dir = Path(__file__).resolve().parent
        self.repo_root = self.script_dir.parent
        
        self._log_dir = self.repo_root / "log"
        self.data_dir = self.repo_root / "data" / f"run_{self.run_id}"
        self.results_dir = self.repo_root / "result" / f"run_{self.run_id}"
        
        self.logger = None
        self.auto_mode = False
        self.id_list = []
        self.fasta_path = None
        self.cleaned_df = None

    def _setup_logger(self) -> logging.Logger:
        os.makedirs(self._log_dir, exist_ok=True)
        logger_name = f"ProteinPipeline.{self.__class__.__name__}"
        logger = logging.getLogger(logger_name)
        
        if logger.hasHandlers():
            return logger
            
        logger.propagate = False
        logger.setLevel(logging.INFO)
        
        log_path = self._log_dir / f"pipeline_{self.run_id}.log"
        formatter = logging.Formatter('%(asctime)s - [%(name)s] - %(levelname)s: %(message)s')
        
        fh = logging.FileHandler(log_path, mode='a', encoding='utf-8', delay=False)
        fh.setFormatter(formatter)
        logger.addHandler(fh)
        
        ch = logging.StreamHandler()
        ch.setFormatter(formatter)
        logger.addHandler(ch)
        
        return logger

    def log_step_analysis(self, step_num: int, title: str, analysis: str):
        self.logger.info("=" * 80)
        self.logger.info(f"STEP {step_num}: {title.upper()} - DETAILED ANALYSIS")
        self.logger.info("=" * 80)
        self.logger.info(analysis.strip())
        self.logger.info("=" * 80 + "\n")

    def clear_screen(self):
        os.system('clear' if os.name == 'posix' else 'cls')

    def print_header(self, title: str):
        print("=" * 60)
        print(f" PROTEIN PIPELINE: {title.upper()} ")
        print("=" * 60)

    def prompt_step(self, step_num: int, title: str, description: str, command_desc: str, execution_callback):
        self.clear_screen()
        self.print_header(f"Step {step_num}: {title}")
        print(f"\n* WHAT IS HAPPENING:\n{description}\n")
        print(f"-> COMMAND:\n-> {command_desc}\n")
        
        if step_num == 1 or not self.auto_mode:
            input("-> Press ENTER to run this step... ")
        else:
            print("-> Running automatically...")
            
        self.logger.info(f"Starting Step {step_num}: {title}")
        print("\n-> Processing... Please wait.\n")
        
        try:
            execution_callback()
            self.logger.info(f"Step {step_num} completed successfully.")
            print("\n-> STEP COMPLETED SUCCESSFULLY!")
        except Exception as e:
            self.logger.error(f"Error in Step {step_num}: {e}")
            print(f"\n-> ERROR: {e}")
            sys.exit(1)
            
        if step_num == 1:
            self.auto_mode = True
            input("\n-> Press ENTER to continue with AUTOMATIC execution of remaining steps... ")
        elif not self.auto_mode:
            input("\n-> Press ENTER to continue... ")

    # ====================== HELPER METHODS ======================
    def _get_cleavage_rules(self):
        return {
            "trypsin": {"cut_after": ["K", "R"], "no_before": ["P"], "desc": "Cuts after Lysine (K) or Arginine (R), except before Proline (P)."},
            "chymotrypsin": {"cut_after": ["F", "W", "Y", "L"], "no_before": ["P"], "desc": "Cuts after aromatic/large hydrophobic residues (F, W, Y, L), except before Proline (P)."},
            "pepsin": {"cut_after": ["F", "L", "W", "Y"], "no_before": [], "desc": "Cuts after aromatic/leucine residues (F, L, W, Y) in acidic conditions."},
            "cnbr": {"cut_after": ["M"], "no_before": [], "desc": "Cyanogen Bromide chemical cleavage: Cuts strictly after Methionine (M)."}
        }

    def generate_mock_ms_intensities(self) -> Path:
        mock_file_path = self.results_dir / "mock_ms_output.tsv"
        self.logger.info(f"Manufacturing mock mass spec expression table -> {mock_file_path}")
        
        input_csv = self.results_dir / "proteomics_data.csv"
        if not input_csv.exists():
            raise FileNotFoundError("Run Steps 1-4 first to collect NCBI Accession references.")
            
        protein_df = pl.read_csv(input_csv)
        accessions = protein_df["Accession_ID"].to_list()
        n = len(accessions)
        
        # Create data with consistent lengths
        mock_data = {
            "Accession_ID": accessions,
            "Potential_contaminant": [None] * n,
            "Reverse": [None] * n,
            "Q_value": [0.001] * n,
            # Generate mock data lists of length n
            "Control_Rep1": [5000000.0] * n,
            "Control_Rep2": [5100000.0] * n,
            "Control_Rep3": [4900000.0] * n,
            "Disease_Rep1": [9500000.0] * n,
            "Disease_Rep2": [9900000.0] * n,
            "Disease_Rep3": [9200000.0] * n
        }
        
        mock_df = pl.DataFrame(mock_data)
        mock_df.write_csv(mock_file_path, separator="\t")
        self.logger.info("Mock mass spectrometry output table successfully prepared.")
        return mock_file_path

    # ====================== ADVANCED PROTEOMICS METHODS ======================
    def clean_proteomics_data(self, search_engine_output_path: Path) -> pl.DataFrame:
        self.logger.info("Cleaning raw proteomics search engine output...")
        df = pl.read_csv(search_engine_output_path, separator="\t")
        
        cleaned_df = df.filter(
            (pl.col("Potential_contaminant") != "+") & 
            (pl.col("Reverse") != "+") & 
            (pl.col("Q_value") <= 0.01)
        )
        self.logger.info(f"Filtered from {df.height} to {cleaned_df.height} high-confidence proteins.")
        return cleaned_df

    def normalize_intensities(self, df: pl.DataFrame, intensity_cols: list) -> pl.DataFrame:
        self.logger.info("Executing Log2 transformation and median normalization...")
        normalized_df = df.with_columns([
            pl.col(col).log(2).alias(col) for col in intensity_cols
        ])
        
        for col in intensity_cols:
            col_median = normalized_df[col].median()
            normalized_df = normalized_df.with_columns((pl.col(col) - col_median).alias(col))
        
        self.logger.info("Normalization complete.")
        return normalized_df

    def calculate_differential_expression(self, df: pl.DataFrame, control_cols: list, disease_cols: list) -> pl.DataFrame:
        self.logger.info("Computing differential expression...")
        expression_df = df.with_columns([
            pl.mean_horizontal(control_cols).alias("Mean_Control"),
            pl.mean_horizontal(disease_cols).alias("Mean_Disease")
        ])
        
        expression_df = expression_df.with_columns(
            (pl.col("Mean_Disease") - pl.col("Mean_Control")).alias("Log2_Fold_Change")
        )
        
        significant_markers = expression_df.filter(pl.col("Log2_Fold_Change").abs() >= 1.0)
        self.logger.info(f"Isolated {significant_markers.height} highly altered marker candidates.")
        return significant_markers

    def generate_pathway_enrichment_link(self, significant_df: pl.DataFrame) -> str:
        self.logger.info("Mapping protein network linkages...")
        protein_list = significant_df["Accession_ID"].to_list()
        query_string = "%0d".join(protein_list)
        string_db_url = f"https://string-db.org/cgi/network?identifiers={query_string}&species=9606"
        self.logger.info(f"Pathway Mapping URL: {string_db_url}")
        return string_db_url

    # ====================== PIPELINE STEPS ======================
    def run_step1(self):
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self.results_dir.mkdir(parents=True, exist_ok=True)
        Entrez.email = __email__
        self.logger.info(f"Run isolated working directories created for ID: {self.run_id}")
        analysis = f"""Step 1: Workspace & NCBI Setup\nDirectories created successfully."""
        self.log_step_analysis(1, "Setup Directories & NCBI", analysis)

    def run_step2(self):
        self.logger.info(f"Querying NCBI Protein database for target: '{self.sample_query}'")
        print("Connecting to NCBI E-Utilities...")
        handle = Entrez.esearch(db="protein", term=self.sample_query, retmax=self.max_results)
        record = Entrez.read(handle)
        handle.close()
        self.id_list = record.get("IdList", [])
        self.logger.info(f"Found {len(self.id_list)} unique matching protein Accession IDs.")
        analysis = f"""Step 2: NCBI Protein Search Results\nTotal Hits: {len(self.id_list)}"""
        self.log_step_analysis(2, "Search NCBI Protein Database", analysis)

    def run_step3(self):
        if not self.id_list:
            raise ValueError("No IDs found from search step.")
        output_path = self.data_dir / "downloaded_proteins.fasta"
        self.logger.info(f"Streaming {len(self.id_list)} entries to {output_path}")
        ids_csv = ",".join(self.id_list)
        handle = Entrez.efetch(db="protein", id=ids_csv, rettype="fasta", retmode="text")
        fasta_data = handle.read()
        handle.close()
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(fasta_data)
        self.fasta_path = str(output_path)
        self.logger.info("Raw FASTA stream complete.")
        analysis = f"""Step 3: Download Completed\nFile: {self.fasta_path}"""
        self.log_step_analysis(3, "Download FASTA Sequences", analysis)

    def run_step4(self):
        if not self.fasta_path or not Path(self.fasta_path).exists():
            raise FileNotFoundError("FASTA file not found.")
        output_csv = self.results_dir / "proteomics_data.csv"
        self.logger.info(f"Parsing FASTA into Polars DataFrame -> {output_csv}")
        protein_records = []
        raw_records = list(SeqIO.parse(self.fasta_path, "fasta"))
        for record in tqdm(raw_records, desc="Parsing FASTA Matrix", ncols=80, ascii=' -#'):
            protein_records.append({
                "Accession_ID": record.id,
                "Description": record.description,
                "Sequence_Length": len(record.seq),
                "Sequence": str(record.seq)
            })
        df = pl.DataFrame(protein_records)
        df.write_csv(output_csv)
        self.logger.info(f"Parsed {df.height} proteins.")
        analysis = f"""Step 4: Parse FASTA to Polars DataFrame\nRecords: {df.height}"""
        self.log_step_analysis(4, "Parse to Analysis-Ready Polars DataFrame", analysis)

    def run_step5(self):
        input_csv = self.results_dir / "proteomics_data.csv"
        output_csv = self.results_dir / f"digested_peptides.csv"
        
        rules_registry = self._get_cleavage_rules()
        print("\n--- Available Proteomics Digestion Enzymes ---")
        for key, info in rules_registry.items():
            print(f" • [{key.upper()}]: {info['desc']}")
            
        chosen = input(f"\nEnter enzyme choice (default 'trypsin'): ").strip().lower()
        if chosen not in rules_registry:
            chosen = "trypsin"
            
        rule = rules_registry[chosen]
        self.logger.info(f"Simulating in-silico {chosen.upper()} digestion -> {output_csv}")
        
        df = pl.read_csv(input_csv)
        peptide_records = []
        
        for row in tqdm(df.iter_rows(named=True), desc="Digesting Proteins", ncols=80, ascii=' -#'):
            seq = row["Sequence"]
            acc = row["Accession_ID"]
            current_peptide = []
            
            for i, amino_acid in enumerate(seq):
                current_peptide.append(amino_acid)
                if amino_acid in rule["cut_after"]:
                    if rule["no_before"] and (i + 1 < len(seq)) and (seq[i + 1] in rule["no_before"]):
                        continue
                    pep_str = "".join(current_peptide)
                    if len(pep_str) >= 6:
                        peptide_records.append({
                            "Parent_Accession": acc,
                            "Enzyme_Used": chosen.upper(),
                            "Peptide_Sequence": pep_str,
                            "Peptide_Length": len(pep_str)
                        })
                    current_peptide = []
            
            if current_peptide:
                pep_str = "".join(current_peptide)
                if len(pep_str) >= 6:
                    peptide_records.append({
                        "Parent_Accession": acc,
                        "Enzyme_Used": chosen.upper(),
                        "Peptide_Sequence": pep_str,
                        "Peptide_Length": len(pep_str)
                    })

        pep_df = pl.DataFrame(peptide_records)
        pep_df.write_csv(output_csv)
        self.logger.info(f"Generated {pep_df.height} unique {chosen.upper()} peptides.")
        
        analysis = f"""Step 5: Dynamic In-Silico Digestion Complete\nEnzyme: {chosen.upper()}\nFragments: {pep_df.height}"""
        self.log_step_analysis(5, "Simulate Enzymatic Digestion", analysis)
        print(f"\n{chosen.upper()} digestion completed: {pep_df.height} peptides generated.")

    def run_step6(self):
        print("\nDo you have a real search engine output file (e.g. proteinGroups.txt)?")
        use_real = input("Enter path to file or press ENTER to use mock data: ").strip()
        
        if use_real and Path(use_real).exists():
            input_file = Path(use_real)
            self.logger.info(f"Using user-provided file: {input_file}")
        else:
            input_file = self.generate_mock_ms_intensities()
        
        self.cleaned_df = self.clean_proteomics_data(input_file)
        output_path = self.results_dir / "cleaned_proteomics.tsv"
        self.cleaned_df.write_csv(output_path, separator="\t")
        print(f"Cleaned data saved to: {output_path}")

    def run_step7(self):
        if not hasattr(self, 'cleaned_df') or self.cleaned_df is None:
            print("Error: Run Step 6 first.")
            return
            
        intensity_cols = [col for col in self.cleaned_df.columns if any(x in col for x in ["Rep", "Control", "Disease"])]
        control_cols = [c for c in intensity_cols if "Control" in c]
        disease_cols = [c for c in intensity_cols if "Disease" in c]
        
        normalized_df = self.normalize_intensities(self.cleaned_df, intensity_cols)
        diff_df = self.calculate_differential_expression(normalized_df, control_cols, disease_cols)
        
        output_diff = self.results_dir / "differential_markers.tsv"
        diff_df.write_csv(output_diff, separator="\t")
        
        url = self.generate_pathway_enrichment_link(diff_df)
        print(f"\nDifferential markers saved: {output_diff}")
        print(f"Pathway network URL: {url}")
        
        open_browser = input("Open STRING-DB network in browser? (y/N): ").strip().lower()
        if open_browser == 'y':
            webbrowser.open(url)

    def run_pipeline(self):
        self.clear_screen()
        print("=" * 60)
        print(" INTERACTIVE NCBI PROTEIN FETCH & ANALYSIS PIPELINE ")
        print("=" * 60)
        print(f"\nVersion : {__version__} | Run ID: {self.run_id}\n")
        
        print(f"Default Query Target: {self.sample_query}")
        change_query = input("-> Change search query parameters? (y/N): ").strip().lower()
        
        if change_query == 'y':
            self.sample_query = input("-> Enter updated search term: ").strip()
            max_res = input(f"-> Max results constraint (default {self.max_results}): ").strip()
            if max_res.isdigit():
                self.max_results = int(max_res)
                
        self.logger = self._setup_logger()
        self.logger.info(f"Pipeline initialized | Target Query: {self.sample_query} | Limit: {self.max_results}")
        
        input("\n-> Configuration locked. Press ENTER to open pipeline loop... ")
        
        self.prompt_step(1, "Setup Workspace", "Generate isolated run storage frameworks.", "os.makedirs", self.run_step1)
        self.prompt_step(2, "Search NCBI Database", "Locate accession targets matching criteria.", "Entrez.esearch", self.run_step2)
        self.prompt_step(3, "Fetch Stream Vectors", "Download continuous FASTA structural files.", "Entrez.efetch", self.run_step3)
        self.prompt_step(4, "Structure Polars Matrix", "Parse sequences into organized data models via Polars.", "SeqIO.parse + pl.DataFrame", self.run_step4)
        self.prompt_step(5, "Simulate Enzymatic Digestion", "Perform dynamic in-silico digestion on protein strings.", "Custom Cleavage Registry Engine", self.run_step5)
        self.prompt_step(6, "Load & Clean MS Data", "Process real or mock mass spec output", "Data Cleaning & Filtering", self.run_step6)
        self.prompt_step(7, "Differential Expression & Pathway Analysis", "Normalization, biomarkers, and network mapping", "Advanced Analytics", self.run_step7)
        
        final_summary = f"""
============================================================
 PROTEIN PIPELINE: PIPELINE EXECUTION COMPLETED SUCCESSFULLY!
============================================================
Run Tracking Token : {self.run_id}
Search Constraints : {self.sample_query}
Raw FASTA Stream   : data/run_{self.run_id}/downloaded_proteins.fasta
Protein Matrix     : result/run_{self.run_id}/proteomics_data.csv
Digested Peptides  : result/run_{self.run_id}/digested_peptides.csv
Cleaned MS Data    : result/run_{self.run_id}/cleaned_proteomics.tsv
Differential Markers: result/run_{self.run_id}/differential_markers.tsv
System Log File    : log/pipeline_{self.run_id}.log
Pipeline finished. Ready for downstream proteomics / mass spectrometry analysis!
"""
        self.logger.info("=" * 60)
        self.logger.info("PIPELINE EXECUTION COMPLETED SUCCESSFULLY!")
        self.logger.info("=" * 60)
        self.logger.info(final_summary.strip())
        
        self.clear_screen()
        self.print_header("Pipeline Execution Completed Successfully!")
        print(final_summary)
        print("Pipeline finished. Ready for downstream proteomics / mass spectrometry analysis!\n")


if __name__ == "__main__":
    pipeline = ProteinPipeline()
    pipeline.run_pipeline()
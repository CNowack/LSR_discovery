import os
import subprocess
import argparse
import shutil
import zipfile

def main():
    parser = argparse.ArgumentParser(description="Download genomes using NCBI datasets CLI")
    parser.add_argument("--batch", required=True, help="Batch ID")
    parser.add_argument("--mapping", required=True, help="Comma-separated Label:Accession mapping")
    parser.add_argument("--outdir", required=True, help="Output directory")
    args = parser.parse_args()

    # Parse the mapping string back into a dictionary
    # Example: "E_coli_K12:GCF_000005845.2,P_aeruginosa:GCF_000006765.1"
    mapping = dict(item.split(":") for item in args.mapping.split(","))

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    for label, accession in mapping.items():
        print(f"--- Processing {label} ({accession}) ---")
        
        zip_path = os.path.join(args.outdir, f"{accession}.zip")
        final_fna = os.path.join(args.outdir, f"{label}.fna")

        # 1. Download the genome package
        # We only need the genomic sequence (fna), not the gff or protein files yet
        cmd = [
            "datasets", "download", "genome", "accession", accession,
            "--include", "genome",
            "--filename", zip_path
        ]
        
        try:
            subprocess.run(cmd, check=True)
            
            # 2. Extract the .fna file from the zip
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                # NCBI datasets structure is usually: 
                # ncbi_dataset/data/{accession}/{accession}_{name}_genomic.fna
                for file_info in zip_ref.infolist():
                    if file_info.filename.endswith(".fna"):
                        # Extract and rename to our workflow's {label}.fna
                        with zip_ref.open(file_info) as source, open(final_fna, "wb") as target:
                            shutil.copyfileobj(source, target)
                        break
            
            # 3. Cleanup the zip file to save space
            os.remove(zip_path)
            print(f"Successfully saved to {final_fna}")

        except subprocess.CalledProcessError as e:
            print(f"Error downloading {accession}: {e}")
        except Exception as e:
            print(f"Error processing {accession}: {e}")

    # Remove the empty ncbi_dataset folder if it exists
    dataset_folder = os.path.join(os.getcwd(), "ncbi_dataset")
    if os.path.exists(dataset_folder):
        shutil.rmtree(dataset_folder)

if __name__ == "__main__":
    main()
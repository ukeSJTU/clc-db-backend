from django.test import TestCase

# Create your tests here.
from django.core.management.base import BaseCommand
from core.models import Molecule, Category, Smile
import pandas as pd


class Command(BaseCommand):
    help = "Load molecule data into the database"

    fill_values = {
        "PubChem CID": 0,  # 对于 'PubChem CID' 列，用 0 填充
        "URL": "",  # 对于 'URL' 列，用空字符串填充
        "SMILES": "N/A",  # 对于 'SMILES' 列，用 'N/A' 填充
        "SMILES_IUPAC": "N/A",  # 对于 'SMILES_IUPAC' 列，用 'N/A' 填充
        "MF": "N/A",  # 对于 'MF' 列，用 'N/A' 填充
        "MW": 0.0,  # 对于 'MW' 列，用 0 填充
        "备注": "N/A",  # 对于 '备注' 列，用 'N/A' 填充
        "手性种类": "无",  # 对于 '手性种类' 列，用 '无' 填充
    }

    def add_arguments(self, parser):
        # Adding an optional argument to specify the path to the CSV file
        parser.add_argument(
            "-p",
            "--path",
            type=str,
            help="Path to the CSV file",
        )

    def handle(self, *args, **kwargs):
        csv_path = kwargs["path"]  # Load your dataframe as needed
        df = pd.read_csv(csv_path).fillna(self.fill_values)

        categories = {
            name: Category.objects.get_or_create(name=name)[0]
            for name in df["类别"].unique()
        }
        smiles_types = {
            smile: Smile.objects.get_or_create(smile=smile)[0]
            for smile in df["手性种类"].unique()
        }

        for idx, row in df.iterrows():
            molecule, created = Molecule.objects.update_or_create(
                cas_id=row["CAS号"],
                defaults={
                    "name": row.get("Name", ""),
                    "pubchem_cid": str(
                        int(row.get("PubChem CID", 0))
                    ),  # Convert to int, if not exist set to 0; then convert to string because PubChem CID is CharField
                    "url": row.get("URL", ""),
                    "pubchem_url": row.get("PubChemURL", ""),
                    "smiles": row.get("SMILES", ""),
                    "smiles_iupac": row.get("SMILES_IUPAC", ""),
                    "molecule_formula": row.get("MF", ""),
                    "molecular_weight": float(row.get("MW", 0)),
                    "remark": row.get("备注", ""),
                },
            )
            # Assign categories and smile types
            molecule.class_type.set([categories[row["类别"]]])
            molecule.smiles_type.set([smiles_types[row["手性种类"]]])

            molecule.save()
            self.stdout.write(
                self.style.SUCCESS(
                    f"Successfully added/updated molecule: {molecule.name}"
                )
            )

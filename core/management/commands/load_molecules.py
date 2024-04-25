from django.test import TestCase

# Create your tests here.
from django.core.management.base import BaseCommand
from core.models import Molecule, Category, Smile
import pandas as pd


class Command(BaseCommand):
    help = "Load molecule data into the database"

    def handle(self, *args, **kwargs):
        df = pd.read_csv(
            "/Users/uke/Downloads/all_data/merged.csv"
        )  # Load your dataframe as needed

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

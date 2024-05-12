from django.core.management.base import BaseCommand
from core.models import Molecule, Category, Smile
import pandas as pd


class Command(BaseCommand):
    help = "Load molecule data into the database"

    fill_values = {
        "PubChem CID": 0,
        "URL": "",
        "SMILES": "N/A",
        "SMILES_IUPAC": "N/A",
        "MF": "N/A",
        "MW": 0.0,
        "备注": "N/A",
        "Heavy Atom Count": 0,
        "Ring Count": 0,
        "Hydrogen Bond Acceptor Count": 0,
        "Hydrogen Bond Donor Count": 0,
        "Rotatable Bond Count": 0,
        "Description": "N/A",
        "types of chirality": "N/A",
    }

    def add_arguments(self, parser):
        parser.add_argument(
            "-p",
            "--path",
            type=str,
            help="Path to the CSV file",
        )

    def handle(self, *args, **kwargs):
        csv_path = kwargs["path"]
        df = pd.read_csv(csv_path).fillna(self.fill_values)

        # Create or get categories and smile types from the dataframe
        category_objs = {}
        for category_list in df["Category"].unique():
            categories = [cat.strip() for cat in category_list.split("、")]
            for category in categories:
                if category not in category_objs:
                    category_objs[category], _ = Category.objects.get_or_create(
                        name=category
                    )

        smile_type_objs = {}
        for smile_type in df["types of chirality"].unique():
            smile_type_objs[smile_type], _ = Smile.objects.get_or_create(
                smile=smile_type
            )

        for idx, row in df.iterrows():
            molecule, created = Molecule.objects.update_or_create(
                cas_id=row["CAS"],
                defaults={
                    "name": row["Name"],
                    "pubchem_cid": str(int(row.get("PubChem CID", 0))),
                    "pubchem_url": row.get("PubChemURL", ""),
                    "smiles": row.get("SMILES", ""),
                    "smiles_iupac": row.get("SMILES_IUPAC", ""),
                    "molecule_formula": row.get("MF", ""),
                    "molecular_weight": float(row.get("MW", 0)),
                    "description": row.get("Description", ""),
                    "heavy_atom_count": int(row.get("Heavy Atom Count", 0)),
                    "ring_count": int(row.get("Ring Count", 0)),
                    "hydrogen_bond_acceptor_count": int(
                        row.get("Hydrogen Bond Acceptor Count", 0)
                    ),
                    "hydrogen_bond_donor_count": int(
                        row.get("Hydrogen Bond Donor Count", 0)
                    ),
                    "rotatable_bond_count": int(row.get("Rotatable Bond Count", 0)),
                },
            )

            categories = [cat.strip() for cat in row["Category"].split("、")]
            molecule.class_type.set(
                [category_objs[cat] for cat in categories if cat in category_objs]
            )

            smile_type = row.get("types of chirality")
            if smile_type in smile_type_objs:
                molecule.chirality.set([smile_type_objs[smile_type]])

            molecule.save()

            self.stdout.write(
                self.style.SUCCESS(
                    f"Successfully added/updated molecule: {molecule.name}"
                )
            )

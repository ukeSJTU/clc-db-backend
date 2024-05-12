from django.db import models


# Define the Category model to store the category of the molecule
class Category(models.Model):
    name = models.CharField(max_length=100)

    def __str__(self):
        return self.name


# Define the Smile model to store the smile representation of the molecule
class Smile(models.Model):
    smile = models.CharField(max_length=100)

    def __str__(self):
        return self.smile


class Molecule(models.Model):
    IUPAC_NAME_MAX_LENGTH = 200
    SMILE_MAX_LENGTH = 200
    MOLECULE_FORMULA_MAX_LENGTH = 100

    name = models.CharField(
        max_length=IUPAC_NAME_MAX_LENGTH
    )  # explicitly set max_length to 200 in case of long names

    cas_id = models.CharField(max_length=20, unique=True)  # this should be unique

    class_type = models.ManyToManyField(
        Category, related_name="molecules"
    )  # class type of the molecule, metadata is in chinese format

    pubchem_url = models.URLField(
        blank=True, null=True
    )  # https://pubchem.ncbi.nlm.nih.gov/compound/46939810

    smiles = models.CharField(max_length=SMILE_MAX_LENGTH)

    chirality = models.ManyToManyField(
        Smile, related_name="molecules"
    )  # smile type of the molecule, metadata is in chinese format

    description = models.TextField(
        blank=True,
        null=True,
    )  # [optional] any additional information about the molecule

    pubchem_cid = models.CharField(
        max_length=20, blank=True, null=True
    )  # this should be unique

    smiles_iupac = models.CharField(
        max_length=IUPAC_NAME_MAX_LENGTH, blank=True, null=True
    )  # IUPAC name of the molecule

    molecule_formula = models.CharField(
        max_length=MOLECULE_FORMULA_MAX_LENGTH,
        blank=True,
        null=True,
    )  # C19H23O3P

    molecular_weight = models.FloatField(blank=True, null=True)  # 330.364

    # Calculated statistics
    heavy_atom_count = models.IntegerField(blank=True, null=True)  # 23
    ring_count = models.IntegerField(blank=True, null=True)  # 2
    hydrogen_bond_acceptor_count = models.IntegerField(blank=True, null=True)  # 3
    hydrogen_bond_donor_count = models.IntegerField(blank=True, null=True)  # 0
    rotatable_bond_count = models.IntegerField(blank=True, null=True)  # 2

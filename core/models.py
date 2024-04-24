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

    # name of the molecule: 3-(tert-Butyl)-4-(2,6-dimethoxyphenyl)-2,3-dihydrobenzo[d][1,3]oxaphosphole
    name = models.CharField(
        max_length=IUPAC_NAME_MAX_LENGTH
    )  # explicitly set max_length to 200 in case of long names

    # cas_id of the molecule: 1246888-90-3
    cas_id = models.CharField(max_length=20, unique=True)  # this should be unique

    # pubchem_cid of the molecule: 46939810
    pubchem_cid = models.CharField(
        max_length=20, blank=True, null=True
    )  # this should be unique

    # class type of the molecule: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    class_type = models.ManyToManyField(
        Category, related_name="molecules"
    )  # class type of the molecule, metadata is in chinese format

    # url to the molecule's external information page, could be empty, if not available default is empty string
    url = models.URLField(
        blank=True, null=True
    )  # https://www.bidepharm.com/products/1246888-90-3.html

    # url to the molecule's pubchem page, could be empty, if not available default is empty string
    pubchem_url = models.URLField(
        blank=True, null=True
    )  # https://pubchem.ncbi.nlm.nih.gov/compound/46939810

    # smiles representation of the molecule: CC(C)(C)P1COC2=CC=CC(=C21)C3=C(C=CC=C3OC)OC
    smiles = models.CharField(max_length=SMILE_MAX_LENGTH)

    # smile type of the smiles representation of the molecule: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    smiles_type = models.ManyToManyField(
        Smile, related_name="molecules"
    )  # smile type of the molecule, metadata is in chinese format

    # iupac name of the molecule: 3-tert-butyl-4-(2,6-dimethoxyphenyl)-2H-1,3-benzoxaphosphole
    smiles_iupac = models.CharField(
        max_length=IUPAC_NAME_MAX_LENGTH, blank=True, null=True
    )  # IUPAC name of the molecule

    # molecular formula of the molecule: C19H23O3P
    # TODO: remember to change blank=False later on
    molecule_formula = models.CharField(
        max_length=MOLECULE_FORMULA_MAX_LENGTH,
        blank=True,
        null=True,
    )  # C19H23O3P

    # molecular weight of the molecule: 330.364
    # TODO: remember to change blank=False later on
    molecular_weight = models.FloatField(blank=True, null=True)  # 330.364

    remark = models.TextField(
        blank=True,
        null=True,
    )  # [optional] any additional information about the molecule

    def __str__(self):
        return "Molecule: {} | CAS ID: {} | Class: {}".format(
            self.name, self.cas_id, self.class_type.all()
        )

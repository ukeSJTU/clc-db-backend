from django.db import models


class Category(models.Model):
    # Category model to store the category of the molecule
    name = models.CharField(max_length=100)

    def __str__(self):
        return self.name


class Smile(models.Model):
    # Smile model to store the smile representation of the molecule
    smile = models.CharField(max_length=100)

    def __str__(self):
        return self.smile


class Molecule(models.Model):
    name = models.CharField(max_length=200)  # observe extra long name in given dataset
    cas_id = models.CharField(max_length=20)  # format: ######-##-#
    class_type = models.ManyToManyField(
        Category, related_name="molecules"
    )  # class type of the molecule, metadata is in chinese format

    url = models.URLField()  # url to the molecule's external information page
    pubchem_url = models.URLField()  # url to the molecule's pubchem page
    smiles = models.CharField(max_length=200)  # smiles representation of the molecule

    smiles_type = models.ManyToManyField(
        Smile, related_name="molecules"
    )  # type of the smiles representation of the molecule

    remark = models.TextField(
        blank=True,
        null=True,
    )  # [optional] any additional information about the molecule

    def __str__(self):
        return "Molecule: {} | CAS ID: {} | Class: {}".format(
            self.name, self.cas_id, self.class_type.all()
        )

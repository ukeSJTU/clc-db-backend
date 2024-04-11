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

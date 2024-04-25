from core.models import Molecule
from core.serializers import CategorySerializer, SmilesSerializer
from rest_framework import serializers


class OverviewSerializer(serializers.ModelSerializer):
    class_type = CategorySerializer(many=True, read_only=True)

    class Meta:
        model = Molecule
        fields = ["name", "cas_id", "class_type"]


class DownloadSerializer(serializers.ModelSerializer):
    class_type = CategorySerializer(many=True, read_only=True)
    smiles_type = SmilesSerializer(many=True, read_only=True)

    class Meta:
        model = Molecule
        fields = [
            "name",
            "cas_id",
            "class_type",
            "smiles",
            "smiles_type",
            "remark",
        ]


class MoleculeSerializer(serializers.ModelSerializer):
    class_type = CategorySerializer(many=True, read_only=True)
    smiles_type = SmilesSerializer(many=True, read_only=True)

    class Meta:
        model = Molecule
        fields = [
            "name",
            "cas_id",
            "class_type",
            "url",
            "pubchem_url",
            "smiles",
            "smiles_type",
            "remark",
        ]


class CardOverviewMoleculeSerializer(serializers.ModelSerializer):
    class_type = CategorySerializer(many=True, read_only=True)

    class Meta:
        model = Molecule
        fields = [
            "name",
            "cas_id",
            "class_type",
            "molecule_formula",
            "molecular_weight",
        ]


class SheetOverviewMoleculeSerializer(serializers.ModelSerializer):
    class_type = CategorySerializer(many=True, read_only=True)

    class Meta:
        model = Molecule
        fields = ["name", "cas_id", "class_type"]


class CompleteMoleculeSerializer(serializers.ModelSerializer):
    class_type = CategorySerializer(many=True, read_only=True)
    smiles_type = SmilesSerializer(many=True, read_only=True)

    class Meta:
        model = Molecule
        fields = "__all__"

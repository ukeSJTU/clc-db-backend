import django_filters
from core.models import Molecule


class MoleculeFilter(django_filters.FilterSet):
    name = django_filters.CharFilter(field_name="name", lookup_expr="icontains")
    cas_id = django_filters.CharFilter(field_name="cas_id", lookup_expr="icontains")
    smiles = django_filters.CharFilter(field_name="smiles", lookup_expr="icontains")

    class Meta:
        model = Molecule
        fields = ["name", "cas_id", "smiles"]

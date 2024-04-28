import django_filters
from core.models import Molecule
from rdkit import Chem
from rdkit.Chem import rdchem, AllChem


# support for filtering with django_filters
class MoleculeFilter(django_filters.FilterSet):
    name = django_filters.CharFilter(field_name="name", lookup_expr="icontains")
    cas_id = django_filters.CharFilter(field_name="cas_id", lookup_expr="icontains")
    smiles = django_filters.CharFilter(field_name="smiles", method="smiles_search")

    class Meta:
        model = Molecule
        fields = ["name", "cas_id", "smiles"]

    def smiles_search(self, queryset, name, value):
        # Convert the input SMILES string to an RDKit molecule
        query_mol = Chem.MolFromSmiles(value)
        if query_mol is None:
            return queryset.none()  # If invalid SMILES, return no results

        # Filter molecules by SMILES using RDKit substructure search
        matching_ids = []
        for molecule in queryset:
            db_mol = Chem.MolFromSmiles(molecule.smiles)
            if db_mol is not None and db_mol.HasSubstructMatch(query_mol):
                matching_ids.append(molecule.id)

        # Return the queryset of molecules that match the substructure search
        return queryset.filter(id__in=matching_ids)

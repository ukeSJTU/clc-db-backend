from core.models import Category, Molecule
from core.views import OverviewPagination
from django.db.models import Prefetch
from rest_framework import filters, generics, pagination, viewsets
from rest_framework.decorators import action
from rest_framework.response import Response

from .serializers import OverviewSerializer, DownloadSerializer, MoleculeSerializer


# cards grid layout with pagination
class OverviewCardViewSet(viewsets.ReadOnlyModelViewSet):
    queryset = Molecule.objects.all()
    serializer_class = OverviewSerializer
    pagination_class = OverviewPagination


# table layout with pagination
class OverviewTableView(viewsets.ReadOnlyModelViewSet):
    queryset = Molecule.objects.all()
    serializer_class = OverviewSerializer
    pagination_class = OverviewPagination


# download api based on molecule's class_type
class DownloadClassesViewSet(viewsets.ModelViewSet):
    queryset = Molecule.objects.all()
    serializer_class = DownloadSerializer

    @action(detail=False, methods=["get"], url_path="all")
    def all_classes(self, request):
        categories = Category.objects.all()
        result = []

        for category in categories:
            molecules = Molecule.objects.filter(class_type=category)
            serializer = self.get_serializer(molecules, many=True)
            result.append({"class_type": category.name, "molecules": serializer.data})

        return Response(result)


# download api based on molecule's cas_id
class DownloadMoleculesViewSet(viewsets.ModelViewSet):
    queryset = Molecule.objects.all()
    serializer_class = DownloadSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["cas_id"]  # TODO: Add other fields here to search across keys

    @action(detail=False, methods=["get"], url_path="all")
    def all_molecules(self, request):
        queryset = self.get_queryset()  # You can filter queryset here if you need
        serializer = self.get_serializer(queryset, many=True)
        return Response(serializer.data)


# api for fetching (search api) single molecule data
class MoleculeViewSet(viewsets.ModelViewSet):
    queryset = Molecule.objects.all()
    serializer_class = MoleculeSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["cas_id"]  # Add other fields here to search across keys

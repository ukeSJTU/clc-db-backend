from core.models import Category, Molecule
from core.views import OverviewPagination, DownloadClassesPagination
from django.db.models import Prefetch
from rest_framework import filters, generics, pagination, viewsets
from rest_framework.decorators import action
from rest_framework.response import Response
from .filters import MoleculeFilter
from django_filters.rest_framework import DjangoFilterBackend

from .serializers import (
    OverviewSerializer,
    DownloadSerializer,
    MoleculeSerializer,
    CompleteMoleculeSerializer,
    CardOverviewMoleculeSerializer,
    CategorySerializer,
)


# cards grid layout with pagination
class OverviewCardViewSet(viewsets.ReadOnlyModelViewSet):
    queryset = Molecule.objects.all()
    serializer_class = CompleteMoleculeSerializer
    # serializer_class = CardOverviewMoleculeSerializer
    pagination_class = OverviewPagination


# table layout with pagination
class OverviewTableView(viewsets.ReadOnlyModelViewSet):
    queryset = Molecule.objects.all()
    serializer_class = CompleteMoleculeSerializer
    # serializer_class = OverviewSerializer
    pagination_class = OverviewPagination


# download api based on molecule's class_type
class DownloadClassesViewSet(viewsets.ModelViewSet):
    queryset = Molecule.objects.all()
    serializer_class = CompleteMoleculeSerializer
    # serializer_class = DownloadSerializer
    pagination_class = DownloadClassesPagination

    @action(detail=False, methods=["get"], url_path="all")
    def all_classes(self, request):
        queryset = Molecule.objects.all()

    serializer_class = CompleteMoleculeSerializer
    # serializer_class = DownloadSerializer
    pagination_class = DownloadClassesPagination

    @action(detail=False, methods=["get"], url_path="all")
    def all_classes(self, request):
        categories = Category.objects.all()
        paginated_categories = self.paginate_queryset(
            categories
        )  # Apply pagination to categories
        result = []

        if paginated_categories is not None:
            for category in paginated_categories:
                molecules = Molecule.objects.filter(class_type=category)
                serializer = self.get_serializer(molecules, many=True)
                result.append(
                    {"class_type": category.name, "molecules": serializer.data}
                )

            return self.get_paginated_response(result)  # Return paginated response
        else:
            # Handle the case if pagination is not applicable or error
            return Response({"error": "Pagination error or no categories found."})


# download api based on molecule's cas_id
class DownloadMoleculesViewSet(viewsets.ModelViewSet):
    queryset = Molecule.objects.all()
    serializer_class = CompleteMoleculeSerializer
    # serializer_class = DownloadSerializer
    pagination_class = DownloadClassesPagination  # Reuse the same pagination settings
    filter_backends = [filters.SearchFilter]
    search_fields = ["cas_id"]  # Extend search fields as needed

    @action(detail=False, methods=["get"], url_path="all")
    def all_molecules(self, request):
        queryset = self.get_queryset()  # You can apply filters here
        paginated_queryset = self.paginate_queryset(
            queryset
        )  # Apply pagination to queryset

        if paginated_queryset is not None:
            serializer = self.get_serializer(paginated_queryset, many=True)
            return self.get_paginated_response(
                serializer.data
            )  # Return paginated response
        else:
            # Handle the case if pagination is not applicable or error
            return Response({"error": "Pagination error or no molecules found."})


# api for fetching (search api) single molecule data
# class MoleculeViewSet(viewsets.ModelViewSet):
#     queryset = Molecule.objects.all()
#     serializer_class = MoleculeSerializer
#     filter_backends = [filters.SearchFilter]
#     search_fields = ["cas_id"]  # Add other fields here to search across keys


# api for searching with multiple query params single molecule data
class MoleculeViewSet(viewsets.ModelViewSet):
    queryset = Molecule.objects.all()
    serializer_class = CompleteMoleculeSerializer  # TODO: this serves search api, but should have its own serializer
    # serializer_class = CardOverviewMoleculeSerializer  # TODO: this serves search api, but should have its own serializer
    filter_backends = (DjangoFilterBackend,)
    filterset_class = MoleculeFilter

    pagination_class = OverviewPagination


class CategoryViewSet(viewsets.ViewSet):
    def list(self, request):
        categories = Category.objects.all()
        serializer = CategorySerializer(categories, many=True)
        return Response(serializer.data)

    def retrieve(self, request, pk=None):
        molecules = Molecule.objects.filter(class_type__id=pk)
        serialized_molecules = MoleculeSerializer(molecules, many=True)
        return Response(serialized_molecules.data)


# api for statistics
class StatisticsViewSet(viewsets.ViewSet):
    def list(self, request):
        total_molecules = Molecule.objects.count()
        total_categories = Category.objects.count()

        return Response(
            {
                "total_molecules": total_molecules,
                "total_categories": total_categories,
            }
        )


# api for weight distribution stats
class WeightDistributionViewSet(viewsets.ViewSet):
    def weight_distribution(self):
        weight_ranges = ["0-100", "100-200", "200-300", "300-400", "400-500", "500+"]

        weight_bins = {
            "0-100": 0,
            "100-200": 0,
            "200-300": 0,
            "300-400": 0,
            "400-500": 0,
            "500+": 0,
        }

        # Adjust weights and classify into bins
        for molecule in Molecule.objects.all():
            weight = molecule.molecular_weight
            if weight is not None:
                if 0 < weight < 100:
                    weight_bins["0-100"] += 1
                elif weight < 200:
                    weight_bins["100-200"] += 1
                elif weight < 300:
                    weight_bins["200-300"] += 1
                elif weight < 400:
                    weight_bins["300-400"] += 1
                elif weight < 500:
                    weight_bins["400-500"] += 1
                else:
                    weight_bins["500+"] += 1

        # Prepare data for the chart
        data = {
            "labels": weight_ranges,
            "values": [weight_bins[range_] for range_ in weight_ranges],
        }

        return data

    def list(self, request):
        return Response(self.weight_distribution())
